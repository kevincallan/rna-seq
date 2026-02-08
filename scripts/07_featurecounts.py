#!/usr/bin/env python3
"""
Step 07 -- featureCounts read aggregation.

Runs featureCounts for each trimming method and each option set defined
in config.  Produces raw count matrices and parsed summary statistics.
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
# featureCounts runner
# ---------------------------------------------------------------------------

def run_featurecounts(
    bam_files: List[Path],
    gtf: str,
    out_path: Path,
    cfg: Dict[str, Any],
    option_set: Dict[str, Any],
) -> Path:
    """Run featureCounts with a specific option set.

    Returns the path to the count matrix output.
    """
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

    # Optional flags from option set
    if option_set.get("B", False):
        cmd.append("-B")
    if option_set.get("P", False):
        cmd.append("-P")
    if option_set.get("C", False):
        cmd.append("-C")
    q_val = option_set.get("Q", 0)
    if q_val and int(q_val) > 0:
        cmd.extend(["-Q", str(q_val)])

    extra = option_set.get("extra_args", "")
    if extra:
        cmd.extend(extra.split())

    # Add BAM files
    cmd.extend(str(b) for b in bam_files)

    run_cmd(cmd, description=f"featureCounts -> {out_path.name}")
    return out_path


# ---------------------------------------------------------------------------
# Clean count matrix (remove .bam suffix, skip comment lines)
# ---------------------------------------------------------------------------

def clean_count_matrix(raw_path: Path, clean_path: Path, samples: list) -> None:
    """Clean the featureCounts output to a simple gene x sample matrix.

    - Removes comment lines (starting with #)
    - Keeps only Geneid and count columns (skips Chr, Start, End, Strand, Length)
    - Strips ``.bam`` and STAR prefix from column headers
    """
    with open(raw_path, "r", encoding="utf-8") as fin, \
         open(clean_path, "w", encoding="utf-8") as fout:

        for line in fin:
            if line.startswith("#"):
                continue

            parts = line.strip().split("\t")

            if parts[0] == "Geneid":
                # Header line -- columns: Geneid, Chr, Start, End, Strand, Length, sample1.bam, ...
                header = ["Geneid"]
                for col in parts[6:]:
                    # Clean: remove path and _Aligned.sortedByCoord.out.bam
                    name = Path(col).name
                    name = re.sub(r"_?Aligned\.sortedByCoord\.out\.bam$", "", name)
                    name = re.sub(r"\.bam$", "", name)
                    # Remove trailing underscore from STAR prefix
                    name = name.rstrip("_")
                    header.append(name)
                fout.write("\t".join(header) + "\n")
            else:
                # Data line: Geneid + counts (skip annotation columns 1-5)
                gene_id = parts[0]
                counts = parts[6:]
                fout.write(gene_id + "\t" + "\t".join(counts) + "\n")

    logger.info("Cleaned count matrix -> %s", clean_path)


# ---------------------------------------------------------------------------
# Parse featureCounts summary
# ---------------------------------------------------------------------------

def parse_fc_summary(summary_path: Path) -> Dict[str, Dict[str, str]]:
    """Parse featureCounts .summary file.

    Returns dict of {category: {sample: count}}.
    """
    stats: Dict[str, Dict[str, str]] = {}
    if not summary_path.exists():
        return stats

    with open(summary_path, encoding="utf-8") as fh:
        header = fh.readline().strip().split("\t")
        sample_names = header[1:]  # first col is "Status"
        for line in fh:
            parts = line.strip().split("\t")
            status = parts[0]
            stats[status] = {}
            for i, sname in enumerate(sample_names):
                stats[status][sname] = parts[i + 1] if i + 1 < len(parts) else "0"

    return stats


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(cfg: Dict[str, Any], methods_override: List[str] | None = None) -> None:
    """Execute step 07."""
    logger.info("=" * 60)
    logger.info("STEP 07: featureCounts")
    logger.info("=" * 60)

    results_dir = Path(cfg["_results_dir"])
    gtf = cfg["references"]["gtf"]

    samples_tsv = Path(cfg.get("_samples_tsv", results_dir / "samples.tsv"))
    samples = read_samples_tsv(samples_tsv)
    methods = methods_override or get_enabled_methods(cfg)
    option_sets = cfg["featurecounts"].get("option_sets", {"default": {}})

    all_summary_rows: List[Dict[str, str]] = []

    for method in methods:
        logger.info("--- featureCounts for method: %s ---", method)
        star_dir = results_dir / method / "star"
        fc_dir = results_dir / method / "featurecounts"
        ensure_dirs(fc_dir)

        # Collect BAM files in deterministic order
        bam_files: List[Path] = []
        for s in samples:
            bam = star_dir / f"{s.sample_name}_Aligned.sortedByCoord.out.bam"
            if bam.exists():
                bam_files.append(bam)
            else:
                logger.warning("BAM not found: %s", bam)

        if not bam_files:
            logger.error("No BAM files for method '%s', skipping", method)
            continue

        for opt_name, opt_params in option_sets.items():
            logger.info("  Option set: %s", opt_name)

            raw_counts = fc_dir / f"counts_{opt_name}.tsv"
            run_featurecounts(bam_files, gtf, raw_counts, cfg, opt_params)

            # Clean matrix
            clean_counts = fc_dir / f"count_matrix_{opt_name}.tsv"
            clean_count_matrix(raw_counts, clean_counts, samples)

            # Parse summary
            summary_file = Path(str(raw_counts) + ".summary")
            summary = parse_fc_summary(summary_file)

            # Record assigned reads for comparison
            assigned = summary.get("Assigned", {})
            for sname, count in assigned.items():
                clean_sname = Path(sname).name
                clean_sname = re.sub(r"_?Aligned\.sortedByCoord\.out\.bam$", "", clean_sname)
                clean_sname = clean_sname.rstrip("_")
                all_summary_rows.append({
                    "method": method,
                    "option_set": opt_name,
                    "sample": clean_sname,
                    "assigned_reads": count,
                })

    # Write featureCounts summary
    fc_summary_path = results_dir / "featurecounts_summary.tsv"
    if all_summary_rows:
        with open(fc_summary_path, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=["method", "option_set", "sample", "assigned_reads"],
                delimiter="\t",
            )
            writer.writeheader()
            for row in all_summary_rows:
                writer.writerow(row)
        logger.info("featureCounts summary -> %s", fc_summary_path)

    logger.info("STEP 07 complete.\n")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Step 07: featureCounts")
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
