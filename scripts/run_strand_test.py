#!/usr/bin/env python3
"""
Strandedness test helper.

Runs featureCounts on a single mapping unit with -s 0, -s 1, -s 2 and
compares assigned-read percentages to recommend the correct strandedness
setting for the dataset.

Usage
-----
# Generic (uses main config):
./py scripts/run_strand_test.py --config config/config.yaml \\
    --dataset GSE48519 --species mouse

# Explicit run-id (reuses existing BAMs from a prior pipeline run):
./py scripts/run_strand_test.py --config config/config.yaml \\
    --dataset GSE48519 --species mouse --run-id my_run_id

Prerequisites
-------------
The mapping step (step 5) must have completed so that BAMs exist in the
results directory.  This script runs only the featureCounts binary (not
the full pipeline) on a single mapping unit to keep it fast.
"""

from __future__ import annotations

import argparse
import csv
import logging
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
from src.analysis_unit import build_mapping_units_with_bams

logger = logging.getLogger(__name__)


def _run_fc_strand(
    bam_files: List[Path],
    gtf: str,
    out_path: Path,
    cfg: Dict[str, Any],
    strandedness: int,
    paired_end: bool,
) -> Path:
    """Run featureCounts with a specific strandedness setting."""
    exe = cfg["tools"].get("featurecounts", "featureCounts")
    threads = cfg["project"].get("threads", 4)
    fc_cfg = cfg.get("featurecounts", {})

    cmd = [
        exe,
        "-a", gtf,
        "-o", str(out_path),
        "-T", str(threads),
        "-t", fc_cfg.get("feature_type", "exon"),
        "-g", fc_cfg.get("attribute", "gene_id"),
        "-s", str(strandedness),
    ]
    if paired_end:
        cmd.extend(["-p", "--countReadPairs"])
    cmd.extend(str(b) for b in bam_files)

    run_cmd(cmd, description=f"featureCounts -s {strandedness}")
    return out_path


def _parse_summary(summary_path: Path) -> Dict[str, int]:
    """Parse a featureCounts .summary file into status -> total count."""
    totals: Dict[str, int] = {}
    if not summary_path.exists():
        return totals
    with open(summary_path, encoding="utf-8") as fh:
        header = fh.readline().strip().split("\t")
        for line in fh:
            parts = line.strip().split("\t")
            status = parts[0]
            total = sum(int(v) for v in parts[1:])
            totals[status] = total
    return totals


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Strandedness test: compare featureCounts -s 0/1/2",
    )
    parser.add_argument("--config", "-c", required=True)
    parser.add_argument("--dataset", default=None)
    parser.add_argument("--species", default=None)
    parser.add_argument("--data-root", default="/data")
    parser.add_argument("--run-id", default=None)
    parser.add_argument("--methods", nargs="*", default=None,
                        help="Trimming methods to test (default: first enabled)")
    args = parser.parse_args()

    cfg = load_config(args.config)

    # Apply dataset/species overrides (reuse logic from run_pipeline)
    if args.dataset:
        dataset_dir = Path(args.data_root) / args.dataset
        cfg["data"]["fastq_dir"] = str(dataset_dir)
        cfg["data"]["metadata_csv"] = str(dataset_dir / "metadata.csv")
    if args.species:
        from scripts.run_pipeline import _apply_data_pointer, SPECIES_MAP
        class _NS:
            pass
        ns = _NS()
        ns.dataset = args.dataset
        ns.species = args.species
        ns.data_root = args.data_root
        ns.metadata = None
        ns.outdir = None
        _apply_data_pointer(ns, cfg)

    run_id = args.run_id or get_run_id(cfg)
    setup_logging(cfg["project"]["logs_dir"], run_id)
    cfg["_run_id"] = run_id
    cfg["_results_dir"] = str(resolve_results_dir(cfg, run_id))
    cfg["_work_dir"] = str(resolve_work_dir(cfg))

    results_dir = Path(cfg["_results_dir"])
    gtf = cfg["references"]["gtf"]
    samples_tsv = results_dir / "samples.tsv"
    if not samples_tsv.exists():
        sys.exit(f"samples.tsv not found at {samples_tsv}. Run steps 0-1 first.")
    samples = read_samples_tsv(samples_tsv)

    methods = args.methods or get_enabled_methods(cfg)[:1]
    mapping_units = build_mapping_units_with_bams(results_dir, samples, methods)
    if not mapping_units:
        sys.exit("No mapping units found. Run step 5 (mapping) first.")

    unit = mapping_units[0]
    bam_files = [unit["sample_to_bam"][s.sample_name]
                 for s in samples
                 if s.sample_name in unit["sample_to_bam"]]
    if not bam_files:
        sys.exit("No BAM files found for the selected mapping unit.")

    paired_end = any(getattr(s, "layout", "paired") == "paired" for s in samples)

    work_dir = Path(cfg["_work_dir"]) / "strand_test"
    ensure_dirs(work_dir)

    print("=" * 60)
    print("Strandedness Test")
    print(f"Unit: {unit['method']}/{unit['mapper']}/{unit['mapper_option_set']}")
    print(f"Samples: {len(bam_files)}")
    print("=" * 60)

    results_table: List[Dict[str, Any]] = []
    for strand_val in [0, 1, 2]:
        out = work_dir / f"strand_{strand_val}.tsv"
        _run_fc_strand(bam_files, gtf, out, cfg, strand_val, paired_end)

        summary = _parse_summary(Path(f"{out}.summary"))
        assigned = summary.get("Assigned", 0)
        total = sum(summary.values()) if summary else 0
        pct = 100.0 * assigned / total if total > 0 else 0.0

        results_table.append({
            "strandedness": strand_val,
            "assigned": assigned,
            "total": total,
            "assigned_pct": pct,
        })
        print(f"  -s {strand_val}: {assigned:,} / {total:,} assigned ({pct:.1f}%)")

    print()
    best = max(results_table, key=lambda r: r["assigned_pct"])
    print(f"Recommendation: use strandedness = {best['strandedness']} "
          f"({best['assigned_pct']:.1f}% assigned)")
    print(f"  -> Set featurecounts.strandedness: {best['strandedness']} in your config YAML")
    print()

    strand_map = {0: "unstranded", 1: "stranded (sense)", 2: "reverse-stranded (antisense)"}
    print(f"Interpretation: {strand_map.get(best['strandedness'], 'unknown')}")
    print("=" * 60)


if __name__ == "__main__":
    main()
