#!/usr/bin/env python3
"""
Step 05 -- Map reads with STAR.

Maps reads for each trimming method separately.  Produces sorted BAM
files, BAM indices, and a mapping summary table parsed from STAR
``Log.final.out``.
"""

from __future__ import annotations

import csv
import logging
import re
import sys
from pathlib import Path
from typing import Any, Dict, List

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
        "--outFilterMultimapNmax", str(star_p.get("outFilterMultimapNmax", 20)),
    ])

    extra = star_p.get("extra_args", "")
    if extra:
        cmd.extend(extra.split())

    run_cmd(cmd, description=f"STAR {sample_name}")

    bam = Path(f"{out_prefix}Aligned.sortedByCoord.out.bam")
    return bam


def index_bam(bam: Path, cfg: Dict[str, Any]) -> None:
    """Index a BAM file with samtools."""
    samtools = cfg["tools"].get("samtools", "samtools")
    threads = cfg["project"].get("threads", 4)
    run_cmd(
        [samtools, "index", "-@", str(threads), str(bam)],
        description=f"samtools index {bam.name}",
    )


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


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(cfg: Dict[str, Any], methods_override: List[str] | None = None) -> None:
    """Execute step 05."""
    logger.info("=" * 60)
    logger.info("STEP 05: STAR mapping")
    logger.info("=" * 60)

    work_dir = resolve_work_dir(cfg)
    results_dir = Path(cfg["_results_dir"])
    genome_dir = cfg["references"]["genome_index"]

    samples_tsv = Path(cfg.get("_samples_tsv", results_dir / "samples.tsv"))
    samples = read_samples_tsv(samples_tsv)
    methods = methods_override or get_enabled_methods(cfg)

    all_mapping_stats: List[Dict[str, str]] = []

    for method in methods:
        logger.info("--- Mapping method: %s ---", method)
        trim_dir = work_dir / "trimmed" / method
        star_out_dir = results_dir / method / "star"
        ensure_dirs(star_out_dir)

        for s in samples:
            logger.info("  Mapping %s ...", s.sample_name)

            if s.layout == "paired":
                r1 = trim_dir / f"{s.sample_name}_1.fastq.gz"
                r2 = trim_dir / f"{s.sample_name}_2.fastq.gz"
            else:
                r1 = trim_dir / f"{s.sample_name}.fastq.gz"
                r2 = None

            prefix = str(star_out_dir / f"{s.sample_name}_")

            bam = map_sample_star(s.sample_name, r1, r2, genome_dir, prefix, cfg)

            # Index BAM
            if bam.exists():
                index_bam(bam, cfg)
            else:
                logger.error("BAM not found after STAR: %s", bam)
                continue

            # Parse STAR log
            star_log = Path(f"{prefix}Log.final.out")
            raw_stats = parse_star_log(star_log)
            summary = extract_mapping_summary(raw_stats)
            summary["method"] = method
            summary["sample"] = s.sample_name
            all_mapping_stats.append(summary)

    # Write mapping summary table
    summary_path = results_dir / "mapping_summary.tsv"
    if all_mapping_stats:
        fields = ["method", "sample"] + [
            k for k in all_mapping_stats[0] if k not in ("method", "sample")
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
