#!/usr/bin/env python3
"""
Step 05 -- Map reads.

Supports multiple mapping approaches (mapper backend + option set), for
example:
  - STAR default vs STAR strict_unique
  - STAR vs HISAT2

Produces mapper-tagged BAM files, BAM indices, and a mapping summary.
"""

from __future__ import annotations

import csv
import logging
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.metadata import read_samples_tsv
from src.utils import (
    ensure_dirs,
    get_enabled_methods,
    get_run_id,
    load_config,
    resolve_results_dir,
    resolve_work_dir,
    run_cmd,
    setup_logging,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# STAR mapping
# ---------------------------------------------------------------------------

def map_sample_star(
    sample_name: str,
    r1: Path,
    r2: Path | None,
    genome_dir: str,
    out_prefix: str,
    cfg: Dict[str, Any],
    option_set: Dict[str, Any],
) -> Path:
    """Run STAR for a single sample and return the BAM path."""
    exe = cfg["tools"].get("star", "STAR")
    threads = cfg["project"].get("threads", 4)
    star_p = cfg.get("star_params", {})

    cmd = [
        exe,
        "--genomeDir", genome_dir,
        "--readFilesIn", str(r1),
    ]

    if r2 is not None:
        cmd.append(str(r2))

    cmd.extend([
        "--readFilesCommand", star_p.get("readFilesCommand", "zcat"),
        "--outSAMtype", "BAM", "SortedByCoordinate",
        "--outFileNamePrefix", out_prefix,
        "--runThreadN", str(threads),
        "--outFilterMultimapNmax", str(
            option_set.get("outFilterMultimapNmax", star_p.get("outFilterMultimapNmax", 20))
        ),
    ])

    extra = option_set.get("extra_args", star_p.get("extra_args", ""))
    if extra:
        cmd.extend(extra.split())

    run_cmd(cmd, description=f"STAR {sample_name}")

    bam = Path(f"{out_prefix}Aligned.sortedByCoord.out.bam")
    return bam


def map_sample_hisat2(
    sample_name: str,
    r1: Path,
    r2: Path | None,
    index_prefix: str,
    out_dir: Path,
    cfg: Dict[str, Any],
    option_set: Dict[str, Any],
) -> Tuple[Path, Path]:
    """Run HISAT2 + samtools sort/index. Return (bam_path, log_path)."""
    hisat2 = cfg["tools"].get("hisat2", "hisat2")
    samtools = cfg["tools"].get("samtools", "samtools")
    threads = cfg["project"].get("threads", 4)
    sam_path = out_dir / f"{sample_name}.hisat2.sam"
    bam_path = out_dir / f"{sample_name}.hisat2.sorted.bam"
    log_path = out_dir / f"{sample_name}.hisat2.log"

    cmd = [
        hisat2,
        "-x", index_prefix,
        "-p", str(threads),
        "-S", str(sam_path),
    ]
    if r2 is None:
        cmd.extend(["-U", str(r1)])
    else:
        cmd.extend(["-1", str(r1), "-2", str(r2)])

    extra = option_set.get("extra_args", "")
    if extra:
        cmd.extend(extra.split())

    result = run_cmd(cmd, description=f"HISAT2 {sample_name}")
    log_path.write_text(result.stderr or "", encoding="utf-8")

    run_cmd(
        [
            samtools, "sort",
            "-@", str(threads),
            "-o", str(bam_path),
            str(sam_path),
        ],
        description=f"samtools sort {sample_name}",
    )
    if sam_path.exists():
        sam_path.unlink()

    return bam_path, log_path


def _bam_context(run_id: str, sample_name: str, method: str, mapper: str, mapper_opt: str) -> str:
    return (
        f"run_id={run_id}, sample={sample_name}, "
        f"branch={method}/{mapper}/{mapper_opt}"
    )


def validate_bam_integrity(
    bam: Path,
    cfg: Dict[str, Any],
    *,
    run_id: str,
    sample_name: str,
    method: str,
    mapper: str,
    mapper_opt: str,
) -> None:
    """Fail early if a BAM appears truncated/corrupt before indexing."""
    context = _bam_context(run_id, sample_name, method, mapper, mapper_opt)
    if not bam.exists():
        raise RuntimeError(f"BAM missing before validation ({context}): {bam}")

    samtools = cfg["tools"].get("samtools", "samtools")
    samtools_exe = shutil.which(samtools)
    if samtools_exe is not None:
        quick = run_cmd(
            [samtools_exe, "quickcheck", "-v", str(bam)],
            description=f"samtools quickcheck {bam.name}",
            check=False,
        )
        if quick.returncode != 0:
            details = (quick.stdout or "") + "\n" + (quick.stderr or "")
            raise RuntimeError(
                "BAM integrity check failed before indexing. "
                f"The BAM appears incomplete/corrupt ({context}). Path: {bam}\n"
                "Likely causes: interrupted write, storage/I/O issue, or filesystem problem.\n"
                "Please verify storage health and rerun this mapping branch.\n"
                f"quickcheck output:\n{details.strip() or '(empty)'}"
            )
        return

    # Fallback when samtools is unavailable: force a full read pass with pysam.
    try:
        import pysam
        with pysam.AlignmentFile(str(bam), "rb") as fh:
            for _ in fh.fetch(until_eof=True):
                pass
    except Exception as exc:
        raise RuntimeError(
            "BAM integrity check failed before indexing. "
            f"The BAM appears incomplete/corrupt ({context}). Path: {bam}\n"
            "Likely causes: interrupted write, storage/I/O issue, or filesystem problem.\n"
            "Please verify storage health and rerun this mapping branch."
        ) from exc


def index_bam(
    bam: Path,
    cfg: Dict[str, Any],
    *,
    run_id: str,
    sample_name: str,
    method: str,
    mapper: str,
    mapper_opt: str,
) -> None:
    """Index a BAM file with actionable diagnostics on failure."""
    context = _bam_context(run_id, sample_name, method, mapper, mapper_opt)
    try:
        import pysam
        logger.info("  Indexing %s with pysam...", bam.name)
        pysam.index(str(bam))
        return
    except ImportError:
        pass
    except Exception as exc:
        raise RuntimeError(
            "BAM indexing failed after integrity validation. "
            f"Context: {context}. BAM: {bam}\n"
            "Likely causes: transient I/O failure, partial write on filesystem, or index-write permission issue."
        ) from exc

    samtools = cfg["tools"].get("samtools", "samtools")
    threads = cfg["project"].get("threads", 4)
    try:
        run_cmd(
            [samtools, "index", "-@", str(threads), str(bam)],
            description=f"samtools index {bam.name}",
        )
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            "samtools index failed after integrity validation. "
            f"Context: {context}. BAM: {bam}\n"
            "Likely causes: transient I/O failure, partial write on filesystem, or filesystem corruption.\n"
            f"Original stderr:\n{exc.stderr or '(empty)'}"
        ) from exc


def filter_bam(
    bam: Path,
    cfg: Dict[str, Any],
    *,
    run_id: str,
    sample_name: str,
    method: str,
    mapper: str,
    mapper_opt: str,
) -> Path:
    """Create a quality-filtered BAM: MAPQ>=255, properly paired, no secondary/supplementary.

    Returns the path to the filtered BAM (indexed).
    """
    import pysam

    filtered = bam.parent / bam.name.replace(".bam", ".filtered.bam")
    if filtered.exists():
        logger.info("  Filtered BAM already exists: %s", filtered.name)
        return filtered

    logger.info("  Filtering %s -> %s", bam.name, filtered.name)
    # flag 2 = properly paired; exclude 2316 = secondary(256) + supplementary(2048) + not_passing_qc(512)
    with pysam.AlignmentFile(str(bam), "rb") as infile:
        with pysam.AlignmentFile(str(filtered), "wb", header=infile.header) as outfile:
            kept = 0
            total = 0
            for read in infile:
                total += 1
                if (read.mapping_quality >= 255
                        and not read.is_secondary
                        and not read.is_supplementary
                        and not read.is_qcfail):
                    outfile.write(read)
                    kept += 1

    validate_bam_integrity(
        filtered,
        cfg,
        run_id=run_id,
        sample_name=sample_name,
        method=method,
        mapper=mapper,
        mapper_opt=mapper_opt,
    )
    index_bam(
        filtered,
        cfg,
        run_id=run_id,
        sample_name=sample_name,
        method=method,
        mapper=mapper,
        mapper_opt=mapper_opt,
    )
    logger.info("    Kept %d / %d reads (%.1f%%)", kept, total,
                100.0 * kept / total if total else 0)
    return filtered


# ---------------------------------------------------------------------------
# Parse STAR Log.final.out
# ---------------------------------------------------------------------------

def parse_star_log(log_path: Path) -> Dict[str, str]:
    """Parse STAR Log.final.out into a dict of key -> value."""
    stats: Dict[str, str] = {}
    if not log_path.exists():
        logger.warning("STAR log not found: %s", log_path)
        return stats

    with open(log_path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if "|" in line:
                parts = line.split("|")
                if len(parts) == 2:
                    key = parts[0].strip()
                    val = parts[1].strip()
                    stats[key] = val
    return stats


def extract_mapping_summary(star_stats: Dict[str, str]) -> Dict[str, str]:
    """Extract the most important STAR mapping metrics."""
    keys_of_interest = {
        "Number of input reads": "input_reads",
        "Uniquely mapped reads number": "uniquely_mapped",
        "Uniquely mapped reads %": "uniquely_mapped_pct",
        "Number of reads mapped to multiple loci": "multi_mapped",
        "% of reads mapped to multiple loci": "multi_mapped_pct",
        "Number of reads unmapped: too many mismatches": "unmapped_mismatch",
        "Number of reads unmapped: too short": "unmapped_short",
        "Number of reads unmapped: other": "unmapped_other",
        "% of reads unmapped: too many mismatches": "unmapped_mismatch_pct",
        "% of reads unmapped: too short": "unmapped_short_pct",
    }
    summary: Dict[str, str] = {}
    for star_key, short_key in keys_of_interest.items():
        summary[short_key] = star_stats.get(star_key, "N/A")
    return summary


def parse_hisat2_log(log_path: Path) -> Dict[str, str]:
    """Parse HISAT2 stderr log into a STAR-like summary schema."""
    summary = {
        "input_reads": "N/A",
        "uniquely_mapped": "N/A",
        "uniquely_mapped_pct": "N/A",
        "multi_mapped": "N/A",
        "multi_mapped_pct": "N/A",
        "unmapped_mismatch": "N/A",
        "unmapped_short": "N/A",
        "unmapped_other": "N/A",
        "unmapped_mismatch_pct": "N/A",
        "unmapped_short_pct": "N/A",
    }
    if not log_path.exists():
        return summary

    text = log_path.read_text(encoding="utf-8", errors="ignore")
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if line.endswith("reads; of these:"):
            # Example: "12345 reads; of these:"
            summary["input_reads"] = line.split(" reads;")[0].replace(",", "").strip()
        elif "aligned exactly 1 time" in line:
            # Capture the first exact-1-time bucket as unique proxy.
            if summary["uniquely_mapped"] == "N/A":
                parts = line.split()
                if parts:
                    summary["uniquely_mapped"] = parts[0].replace(",", "")
                if "(" in line and "%" in line:
                    pct = line.split("(")[1].split(")")[0]
                    summary["uniquely_mapped_pct"] = pct
        elif "aligned >1 times" in line:
            parts = line.split()
            if parts:
                summary["multi_mapped"] = parts[0].replace(",", "")
            if "(" in line and "%" in line:
                pct = line.split("(")[1].split(")")[0]
                summary["multi_mapped_pct"] = pct
        elif "aligned 0 times" in line:
            # HISAT2 doesn't break out mismatch/short/other categories.
            parts = line.split()
            if parts:
                summary["unmapped_other"] = parts[0].replace(",", "")
            if "(" in line and "%" in line:
                pct = line.split("(")[1].split(")")[0]
                summary["unmapped_short_pct"] = pct
    return summary


def get_mapping_runs(cfg: Dict[str, Any]) -> List[Tuple[str, str, Dict[str, Any]]]:
    """Return enabled (mapper, option_set_name, option_set_params) tuples."""
    mapping_cfg = cfg.get("mapping", {}).get("backends", {})
    runs: List[Tuple[str, str, Dict[str, Any]]] = []

    if mapping_cfg:
        for mapper, mapper_cfg in mapping_cfg.items():
            if not isinstance(mapper_cfg, dict) or not mapper_cfg.get("enabled", False):
                continue
            option_sets = mapper_cfg.get("option_sets", {"default": {}})
            for opt_name, opt_params in option_sets.items():
                runs.append((mapper, opt_name, dict(opt_params or {})))

    if not runs:
        # Backward-compatible fallback: STAR single run from old config.
        runs.append(("star", "default", {}))

    return runs


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(cfg: Dict[str, Any], methods_override: List[str] | None = None) -> None:
    """Execute step 05."""
    logger.info("=" * 60)
    logger.info("STEP 05: Read mapping")
    logger.info("=" * 60)

    work_dir = resolve_work_dir(cfg)
    results_dir = Path(cfg["_results_dir"])
    run_id = str(cfg.get("_run_id", "unknown"))
    samples_tsv = Path(cfg.get("_samples_tsv", results_dir / "samples.tsv"))
    samples = read_samples_tsv(samples_tsv)
    methods = methods_override or get_enabled_methods(cfg)
    mapping_runs = get_mapping_runs(cfg)
    mapping_cfg = cfg.get("mapping", {}).get("backends", {})

    all_mapping_stats: List[Dict[str, str]] = []

    for method in methods:
        logger.info("--- Trimming method: %s ---", method)
        trim_dir = work_dir / "trimmed" / method
        for mapper, map_opt_name, map_opt_params in mapping_runs:
            logger.info("  Mapping approach: %s / %s", mapper, map_opt_name)
            map_out_dir = results_dir / method / "mapping" / mapper / map_opt_name
            ensure_dirs(map_out_dir)

            mapper_cfg = mapping_cfg.get(mapper, {})
            for s in samples:
                logger.info("    Mapping %s ...", s.sample_name)

                if s.layout == "paired":
                    r1 = trim_dir / f"{s.sample_name}_1.fastq.gz"
                    r2 = trim_dir / f"{s.sample_name}_2.fastq.gz"
                else:
                    r1 = trim_dir / f"{s.sample_name}.fastq.gz"
                    r2 = None

                if mapper == "hisat2":
                    index_prefix = str(mapper_cfg.get("index_prefix", "")).strip()
                    if not index_prefix:
                        logger.warning(
                            "HISAT2 enabled but index_prefix is empty; skipping %s/%s",
                            method, s.sample_name,
                        )
                        continue
                    bam, log_path = map_sample_hisat2(
                        s.sample_name, r1, r2, index_prefix, map_out_dir, cfg, map_opt_params,
                    )
                    if bam.exists():
                        validate_bam_integrity(
                            bam,
                            cfg,
                            run_id=run_id,
                            sample_name=s.sample_name,
                            method=method,
                            mapper=mapper,
                            mapper_opt=map_opt_name,
                        )
                        index_bam(
                            bam,
                            cfg,
                            run_id=run_id,
                            sample_name=s.sample_name,
                            method=method,
                            mapper=mapper,
                            mapper_opt=map_opt_name,
                        )
                    summary = parse_hisat2_log(log_path)
                else:
                    genome_dir = str(
                        mapper_cfg.get("genome_index", cfg["references"]["genome_index"])
                    )
                    prefix = str(map_out_dir / f"{s.sample_name}_")
                    bam = map_sample_star(
                        s.sample_name, r1, r2, genome_dir, prefix, cfg, map_opt_params,
                    )
                    if bam.exists():
                        validate_bam_integrity(
                            bam,
                            cfg,
                            run_id=run_id,
                            sample_name=s.sample_name,
                            method=method,
                            mapper=mapper,
                            mapper_opt=map_opt_name,
                        )
                        index_bam(
                            bam,
                            cfg,
                            run_id=run_id,
                            sample_name=s.sample_name,
                            method=method,
                            mapper=mapper,
                            mapper_opt=map_opt_name,
                        )
                    else:
                        logger.error("BAM not found after STAR: %s", bam)
                        continue
                    star_log = Path(f"{prefix}Log.final.out")
                    raw_stats = parse_star_log(star_log)
                    summary = extract_mapping_summary(raw_stats)

                filtered_bam = filter_bam(
                    bam,
                    cfg,
                    run_id=run_id,
                    sample_name=s.sample_name,
                    method=method,
                    mapper=mapper,
                    mapper_opt=map_opt_name,
                )

                summary["trim_method"] = method
                summary["mapper"] = mapper
                summary["mapper_option_set"] = map_opt_name
                summary["sample"] = s.sample_name
                summary["bam_path"] = str(bam)
                summary["filtered_bam_path"] = str(filtered_bam)
                all_mapping_stats.append(summary)

    # Write mapping summary table
    summary_path = results_dir / "mapping_summary.tsv"
    if all_mapping_stats:
        fields = [
            "trim_method", "mapper", "mapper_option_set", "sample", "bam_path", "filtered_bam_path",
            "input_reads", "uniquely_mapped", "uniquely_mapped_pct",
            "multi_mapped", "multi_mapped_pct",
            "unmapped_mismatch", "unmapped_short", "unmapped_other",
            "unmapped_mismatch_pct", "unmapped_short_pct",
        ]
        with open(summary_path, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                                    extrasaction="ignore")
            writer.writeheader()
            for row in all_mapping_stats:
                writer.writerow(row)
        logger.info("Mapping summary -> %s", summary_path)

    logger.info("STEP 05 complete.\n")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Step 05: STAR mapping")
    parser.add_argument("--config", required=True)
    parser.add_argument("--run-id", default=None)
    parser.add_argument("--methods", nargs="*", default=None)
    args = parser.parse_args()

    cfg = load_config(args.config)
    run_id = args.run_id or get_run_id(cfg)
    setup_logging(cfg["project"]["logs_dir"], run_id)
    cfg["_run_id"] = run_id
    cfg["_results_dir"] = str(resolve_results_dir(cfg, run_id))
    cfg["_work_dir"] = str(resolve_work_dir(cfg))

    main(cfg, methods_override=args.methods)
