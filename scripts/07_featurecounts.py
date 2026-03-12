#!/usr/bin/env python3
"""
Step 07 -- Read count aggregation.

Supports two backends (set ``featurecounts.backend`` in config):

  - ``featurecounts``: calls the subread featureCounts binary (DEFAULT on
    servers where it is installed).
  - ``htseq``: uses HTSeq-count via Python -- ``pip install HTSeq``.
    No compiled binaries needed beyond pysam.

Both backends produce the same output: a clean gene x sample count matrix.
"""

from __future__ import annotations

import csv
import logging
import re
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


# =========================================================================
# Backend 1: featureCounts (subread binary)
# =========================================================================

def run_featurecounts(
    bam_files: List[Path],
    gtf: str,
    out_path: Path,
    cfg: Dict[str, Any],
    option_set: Dict[str, Any],
    paired_end: bool = False,
) -> Path:
    """Run featureCounts with a specific option set."""
    exe = cfg["tools"].get("featurecounts", "featureCounts")
    threads = cfg["project"].get("threads", 4)
    fc_cfg = cfg["featurecounts"]

    cmd = [
        exe,
        "-a", gtf,
        "-o", str(out_path),
        "-T", str(threads),
        "-t", fc_cfg.get("feature_type", "exon"),
        "-g", fc_cfg.get("attribute", "gene_id"),
        "-s", str(fc_cfg.get("strandedness", 2)),
    ]

    if paired_end:
        cmd.extend(["-p", "--countReadPairs"])

    if option_set.get("B", False):
        cmd.append("-B")
    if option_set.get("P", False):
        cmd.append("-P")
    if option_set.get("C", False):
        cmd.append("-C")
    if option_set.get("M", False):
        cmd.append("-M")
    if option_set.get("fraction", False):
        cmd.append("--fraction")
    q_val = option_set.get("Q", 0)
    if q_val and int(q_val) > 0:
        cmd.extend(["-Q", str(q_val)])

    extra = option_set.get("extra_args", "")
    if extra:
        cmd.extend(extra.split())

    cmd.extend(str(b) for b in bam_files)

    run_cmd(cmd, description=f"featureCounts -> {out_path.name}")
    return out_path


# =========================================================================
# Backend 2: HTSeq (pure Python)
# =========================================================================

def run_htseq_count(
    bam_files: List[Path],
    gtf: str,
    out_path: Path,
    cfg: Dict[str, Any],
    option_set: Dict[str, Any],
    sample_names: List[str],
) -> Path:
    """Count reads using HTSeq (pure Python, no compiled binary needed).

    Produces a count matrix in the same format as featureCounts clean output:
    Geneid<TAB>sample1<TAB>sample2<TAB>...

    Parameters
    ----------
    bam_files : list of Path
        Sorted, indexed BAM files (one per sample).
    gtf : str
        Path to GTF annotation.
    out_path : Path
        Output path for the combined count matrix.
    cfg : dict
        Pipeline config.
    option_set : dict
        Option set (Q for MAPQ filter, etc.).
    sample_names : list of str
        Sample names corresponding to bam_files.
    """
    import HTSeq

    fc_cfg = cfg["featurecounts"]
    strandedness = int(fc_cfg.get("strandedness", 2))
    feature_type = fc_cfg.get("feature_type", "exon")
    attribute = fc_cfg.get("attribute", "gene_id")
    minaqual = int(option_set.get("Q", 0))

    # Map strandedness int to HTSeq string
    strand_map = {0: "no", 1: "yes", 2: "reverse"}
    stranded = strand_map.get(strandedness, "reverse")

    logger.info("  HTSeq: loading GTF features (type=%s, attr=%s)...",
                feature_type, attribute)

    # Build feature array from GTF
    gtf_features = HTSeq.GenomicArrayOfSets("auto", stranded=(stranded != "no"))
    gene_ids = set()

    for feature in HTSeq.GFF_Reader(gtf):
        if feature.type == feature_type:
            gid = feature.attr.get(attribute, None)
            if gid:
                gene_ids.add(gid)
                if stranded != "no":
                    gtf_features[feature.iv] += gid
                else:
                    gtf_features[feature.iv] += gid

    gene_ids = sorted(gene_ids)
    logger.info("  HTSeq: %d genes from GTF", len(gene_ids))

    # Count each BAM file
    all_counts: Dict[str, Dict[str, int]] = {}

    for bam_path, sname in zip(bam_files, sample_names):
        logger.info("  HTSeq: counting %s ...", sname)
        counts: Dict[str, int] = {g: 0 for g in gene_ids}
        n_assigned = 0
        n_total = 0

        bam_reader = HTSeq.BAM_Reader(str(bam_path))

        for alnmt in bam_reader:
            n_total += 1

            if not alnmt.aligned:
                continue
            if alnmt.aQual < minaqual:
                continue

            gene_set = set()
            for iv, val in gtf_features[alnmt.iv].steps():
                gene_set |= val

            if len(gene_set) == 1:
                gene = list(gene_set)[0]
                counts[gene] += 1
                n_assigned += 1

        all_counts[sname] = counts
        logger.info("    %s: %d/%d reads assigned", sname, n_assigned, n_total)

    # Write combined count matrix (Geneid + sample columns)
    ensure_dirs(out_path.parent)
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write("Geneid\t" + "\t".join(sample_names) + "\n")
        for gene in gene_ids:
            row = [gene] + [str(all_counts[s].get(gene, 0)) for s in sample_names]
            fh.write("\t".join(row) + "\n")

    logger.info("  HTSeq count matrix -> %s", out_path)
    return out_path


# =========================================================================
# Clean count matrix (for featureCounts raw output)
# =========================================================================

def clean_count_matrix(raw_path: Path, clean_path: Path, samples: list) -> None:
    """Clean featureCounts output to a simple gene x sample matrix.

    Removes comment lines, strips annotation columns and .bam suffix.
    """
    with open(raw_path, "r", encoding="utf-8") as fin, \
         open(clean_path, "w", encoding="utf-8") as fout:

        for line in fin:
            if line.startswith("#"):
                continue

            parts = line.strip().split("\t")

            if parts[0] == "Geneid":
                header = ["Geneid"]
                for col in parts[6:]:
                    name = Path(col).name
                    name = re.sub(r"_?Aligned\.sortedByCoord\.out(?:\.filtered)?\.bam$", "", name)
                    name = re.sub(r"\.filtered\.bam$", "", name)
                    name = re.sub(r"\.bam$", "", name)
                    name = name.rstrip("_")
                    header.append(name)
                fout.write("\t".join(header) + "\n")
            else:
                gene_id = parts[0]
                counts = parts[6:]
                fout.write(gene_id + "\t" + "\t".join(counts) + "\n")

    logger.info("Cleaned count matrix -> %s", clean_path)


def parse_fc_summary(summary_path: Path) -> Dict[str, Dict[str, str]]:
    """Parse featureCounts .summary file."""
    stats: Dict[str, Dict[str, str]] = {}
    if not summary_path.exists():
        return stats

    with open(summary_path, encoding="utf-8") as fh:
        header = fh.readline().strip().split("\t")
        sample_names = header[1:]
        for line in fh:
            parts = line.strip().split("\t")
            status = parts[0]
            stats[status] = {}
            for i, sname in enumerate(sample_names):
                stats[status][sname] = parts[i + 1] if i + 1 < len(parts) else "0"

    return stats


def build_mapping_units(
    results_dir: Path,
    samples: list,
    methods: List[str],
    use_filtered_bam: bool = False,
) -> List[Dict[str, Any]]:
    """Collect mapping units from mapping_summary.tsv (or legacy STAR layout).

    If use_filtered_bam is True and mapping_summary has filtered_bam_path,
    use that instead of bam_path for counting (post-step-5 filtered BAMs).
    """
    mapping_summary = results_dir / "mapping_summary.tsv"
    units: Dict[Tuple[str, str, str], Dict[str, Any]] = {}

    if mapping_summary.exists():
        with open(mapping_summary, encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                method = row.get("method") or row.get("trim_method") or ""
                mapper = row.get("mapper", "star")
                mapper_opt = row.get("mapper_option_set", "default")
                sample = row.get("sample", "")
                bam_path = row.get("bam_path", "")
                if use_filtered_bam and row.get("filtered_bam_path"):
                    bam_path = row.get("filtered_bam_path", "")
                if method not in methods or not sample or not bam_path:
                    continue
                key = (method, mapper, mapper_opt)
                if key not in units:
                    units[key] = {
                        "method": method,
                        "mapper": mapper,
                        "mapper_option_set": mapper_opt,
                        "sample_to_bam": {},
                    }
                units[key]["sample_to_bam"][sample] = Path(bam_path)

    if not units:
        # Legacy fallback (single STAR run per trimming method)
        for method in methods:
            sample_to_bam: Dict[str, Path] = {}
            star_dir = results_dir / method / "star"
            for s in samples:
                bam = star_dir / f"{s.sample_name}_Aligned.sortedByCoord.out.bam"
                if bam.exists():
                    sample_to_bam[s.sample_name] = bam
            if sample_to_bam:
                units[(method, "star", "default")] = {
                    "method": method,
                    "mapper": "star",
                    "mapper_option_set": "default",
                    "sample_to_bam": sample_to_bam,
                }

    return [units[k] for k in sorted(units.keys())]


def _prepare_gtf(cfg: Dict[str, Any], results_dir: Path) -> str:
    """Return the GTF path to use, optionally pre-filtering for an attribute.

    When ``featurecounts.gtf_filter`` is set (e.g. ``"gene_name"``), rows
    lacking that attribute are stripped out.  The filtered file is cached in
    the work directory so it only runs once.
    """
    raw_gtf = cfg["references"]["gtf"]
    gtf_filter = cfg.get("featurecounts", {}).get("gtf_filter", "")
    if not gtf_filter:
        return raw_gtf

    work_dir = Path(cfg.get("_work_dir", "work"))
    ensure_dirs(work_dir)
    filtered = work_dir / f"filtered_{gtf_filter}_{Path(raw_gtf).stem.replace('.', '_')}.gtf"
    if filtered.exists():
        logger.info("Using cached filtered GTF: %s", filtered)
        return str(filtered)

    import gzip
    import shutil

    logger.info("Pre-filtering GTF for attribute '%s' ...", gtf_filter)
    opener = gzip.open if raw_gtf.endswith(".gz") else open
    n_in = 0
    n_out = 0
    with opener(raw_gtf, "rt", encoding="utf-8", errors="replace") as fin, \
         open(filtered, "w", encoding="utf-8") as fout:
        for line in fin:
            n_in += 1
            if line.startswith("#") or gtf_filter in line:
                fout.write(line)
                n_out += 1

    logger.info("  GTF filtered: %d -> %d lines (kept rows with '%s')",
                n_in, n_out, gtf_filter)
    return str(filtered)


# =========================================================================
# Main
# =========================================================================

def main(cfg: Dict[str, Any], methods_override: List[str] | None = None) -> None:
    """Execute step 07."""
    logger.info("=" * 60)
    logger.info("STEP 07: Read count aggregation")
    logger.info("=" * 60)

    results_dir = Path(cfg["_results_dir"])
    gtf = _prepare_gtf(cfg, results_dir)

    samples_tsv = Path(cfg.get("_samples_tsv", results_dir / "samples.tsv"))
    samples = read_samples_tsv(samples_tsv)
    methods = methods_override or get_enabled_methods(cfg)
    option_sets = cfg["featurecounts"].get("option_sets", {"default": {}})

    # Determine counting backend
    backend = cfg["featurecounts"].get("backend", "featurecounts")
    logger.info("Counting backend: %s", backend)

    all_summary_rows: List[Dict[str, str]] = []
    use_filtered_bam = cfg.get("featurecounts", {}).get("use_filtered_bam", False)
    if use_filtered_bam:
        logger.info("Using filtered BAMs for counting (featurecounts.use_filtered_bam=true)")
    mapping_units = build_mapping_units(
        results_dir, samples, methods, use_filtered_bam=use_filtered_bam
    )

    for unit in mapping_units:
        method = unit["method"]
        mapper = unit["mapper"]
        mapper_opt = unit["mapper_option_set"]
        logger.info(
            "--- Counting for trim=%s mapper=%s mapper_option=%s ---",
            method, mapper, mapper_opt,
        )
        fc_dir = results_dir / method / "featurecounts" / mapper / mapper_opt
        ensure_dirs(fc_dir)

        # Collect BAM files in deterministic order
        bam_files: List[Path] = []
        sample_names: List[str] = []
        for s in samples:
            bam = unit["sample_to_bam"].get(s.sample_name)
            if bam is not None and bam.exists():
                bam_files.append(bam)
                sample_names.append(s.sample_name)
            else:
                logger.warning(
                    "BAM not found for trim=%s mapper=%s opt=%s sample=%s",
                    method, mapper, mapper_opt, s.sample_name,
                )

        if not bam_files:
            logger.error("No BAM files for method '%s', skipping", method)
            continue

        paired_end = any(getattr(s, "layout", "paired") == "paired" for s in samples)

        for opt_name, opt_params in option_sets.items():
            logger.info("  Option set: %s", opt_name)

            if backend == "htseq":
                # HTSeq produces a clean matrix directly
                clean_counts = fc_dir / f"count_matrix_{opt_name}.tsv"
                run_htseq_count(
                    bam_files, gtf, clean_counts, cfg, opt_params, sample_names,
                )
                # No separate summary file for HTSeq

            else:
                # featureCounts binary
                raw_counts = fc_dir / f"counts_{opt_name}.tsv"
                run_featurecounts(
                    bam_files,
                    gtf,
                    raw_counts,
                    cfg,
                    opt_params,
                    paired_end=paired_end,
                )

                # Clean matrix
                clean_counts = fc_dir / f"count_matrix_{opt_name}.tsv"
                clean_count_matrix(raw_counts, clean_counts, samples)

                # Parse summary — capture ALL status categories
                summary_file = Path(str(raw_counts) + ".summary")
                summary = parse_fc_summary(summary_file)

                if summary:
                    first_status = next(iter(summary.values()))
                    raw_sample_names = list(first_status.keys())
                else:
                    raw_sample_names = []

                for raw_sname in raw_sample_names:
                    clean_sname = Path(raw_sname).name
                    clean_sname = re.sub(
                        r"_?Aligned\.sortedByCoord\.out(?:\.filtered)?\.bam$", "", clean_sname
                    )
                    clean_sname = clean_sname.rstrip("_")
                    row_data: Dict[str, str] = {
                        "trim_method": method,
                        "mapper": mapper,
                        "mapper_option_set": mapper_opt,
                        "option_set": opt_name,
                        "sample": clean_sname,
                    }
                    for status, sname_vals in summary.items():
                        col = status.replace(" ", "_")
                        row_data[col] = sname_vals.get(raw_sname, "0")
                    all_summary_rows.append(row_data)

    # Write featureCounts summary with all status categories
    fc_summary_path = results_dir / "featurecounts_summary.tsv"
    if all_summary_rows:
        key_cols = [
            "trim_method", "mapper", "mapper_option_set",
            "option_set", "sample",
        ]
        stat_cols = sorted(
            {k for row in all_summary_rows for k in row if k not in key_cols}
        )
        # Put Assigned first among stat columns
        if "Assigned" in stat_cols:
            stat_cols.remove("Assigned")
            stat_cols = ["Assigned"] + stat_cols
        fieldnames = key_cols + stat_cols

        with open(fc_summary_path, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(
                fh, fieldnames=fieldnames, delimiter="\t",
                extrasaction="ignore",
            )
            writer.writeheader()
            for row in all_summary_rows:
                writer.writerow(row)
        logger.info("featureCounts summary -> %s", fc_summary_path)

    logger.info("STEP 07 complete.\n")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Step 07: Read counting")
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
