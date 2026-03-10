#!/usr/bin/env python3
"""
Step 06 -- Generate BigWig coverage tracks.

Produces per-sample BigWig files for each trimming method:
  - CPM-normalised with MAPQ filter (unique mappers only)
  - Optionally DESeq2 size-factor scaled (if size factors file exists)
"""

from __future__ import annotations

import logging
import csv
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

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


def build_mapping_units(
    results_dir: Path,
    samples: list,
    methods: List[str],
    prefer_filtered: bool = False,
) -> List[Dict[str, Any]]:
    """Collect mapping units and per-sample BAMs.

    When *prefer_filtered* is True, use the pre-filtered BAM
    (``filtered_bam_path`` in mapping_summary.tsv) if it exists.
    """
    mapping_summary = results_dir / "mapping_summary.tsv"
    units: Dict[tuple[str, str, str], Dict[str, Any]] = {}
    if mapping_summary.exists():
        with open(mapping_summary, encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                method = row.get("trim_method") or row.get("method") or ""
                mapper = row.get("mapper", "star")
                mapper_opt = row.get("mapper_option_set", "default")
                sample = row.get("sample", "")
                bam_path = row.get("bam_path", "")
                if method not in methods or not sample or not bam_path:
                    continue

                if prefer_filtered:
                    filtered = row.get("filtered_bam_path", "")
                    if filtered and Path(filtered).exists():
                        bam_path = filtered

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
        for method in methods:
            star_dir = results_dir / method / "star"
            sample_to_bam: Dict[str, Path] = {}
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


def make_bigwig(
    bam: Path,
    out_bw: Path,
    cfg: Dict[str, Any],
    scale_factor: Optional[float] = None,
) -> None:
    """Run bamCoverage to create a BigWig file.

    Parameters
    ----------
    bam : Path
        Input BAM file (must be sorted and indexed).
    out_bw : Path
        Output BigWig path.
    cfg : dict
        Pipeline config.
    scale_factor : float, optional
        If provided, use ``--scaleFactor`` instead of ``--normalizeUsing``.
    """
    exe = cfg["tools"].get("bamcoverage", "bamCoverage")
    bw_cfg = cfg.get("bigwig", {})
    threads = cfg["project"].get("threads", 4)

    cmd = [
        exe,
        "-b", str(bam),
        "-o", str(out_bw),
        "--minMappingQuality", str(bw_cfg.get("mapq_filter", 255)),
        "--binSize", str(bw_cfg.get("bin_size", 50)),
        "-p", str(threads),
    ]

    if scale_factor is not None:
        cmd.extend(["--scaleFactor", str(scale_factor)])
    else:
        cmd.extend(["--normalizeUsing", bw_cfg.get("normalization", "CPM")])

    run_cmd(cmd, description=f"bamCoverage -> {out_bw.name}")


def load_size_factors(path: Path) -> Dict[str, float]:
    """Load DESeq2 size factors from a TSV (sample<TAB>factor)."""
    factors: Dict[str, float] = {}
    if not path.exists():
        return factors
    with open(path, encoding="utf-8") as fh:
        header = fh.readline()  # skip header
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                factors[parts[0]] = float(parts[1])
    return factors


def main(cfg: Dict[str, Any], methods_override: List[str] | None = None) -> None:
    """Execute step 06."""
    logger.info("=" * 60)
    logger.info("STEP 06: BigWig generation")
    logger.info("=" * 60)

    results_dir = Path(cfg["_results_dir"])
    bw_cfg = cfg.get("bigwig", {})

    samples_tsv = Path(cfg.get("_samples_tsv", results_dir / "samples.tsv"))
    samples = read_samples_tsv(samples_tsv)
    methods = methods_override or get_enabled_methods(cfg)

    mapping_units = build_mapping_units(results_dir, samples, methods,
                                        prefer_filtered=True)
    for unit in mapping_units:
        method = unit["method"]
        mapper = unit["mapper"]
        mapper_opt = unit["mapper_option_set"]
        logger.info(
            "--- BigWigs for trim=%s mapper=%s mapper_option=%s ---",
            method, mapper, mapper_opt,
        )
        bw_dir = results_dir / method / "bigwig" / mapper / mapper_opt
        ensure_dirs(bw_dir)

        # Try to load DESeq2 size factors for this mapping unit.
        sf_path = results_dir / method / "deseq2" / mapper / mapper_opt / "size_factors.tsv"
        size_factors = load_size_factors(sf_path)
        use_sf = bw_cfg.get("use_deseq2_sizefactors", True) and bool(size_factors)

        for s in samples:
            bam = unit["sample_to_bam"].get(s.sample_name)
            if bam is None or not bam.exists():
                logger.warning("BAM not found for %s (%s/%s/%s)",
                               s.sample_name, method, mapper, mapper_opt)
                continue

            # CPM normalised BigWig
            bw_cpm = bw_dir / f"{s.sample_name}.CPM.MAPQ255.bw"
            logger.info("  %s -> CPM BigWig", s.sample_name)
            make_bigwig(bam, bw_cpm, cfg)

            # DESeq2 size-factor scaled BigWig
            if use_sf and s.sample_name in size_factors:
                sf = size_factors[s.sample_name]
                # bamCoverage scaleFactor = 1/sf (inverse of DESeq2 factor)
                inv_sf = 1.0 / sf if sf != 0 else 1.0
                bw_deseq = bw_dir / f"{s.sample_name}.DESeq2scaled.MAPQ255.bw"
                logger.info("  %s -> DESeq2-scaled BigWig (sf=%.4f, 1/sf=%.4f)",
                            s.sample_name, sf, inv_sf)
                make_bigwig(bam, bw_deseq, cfg, scale_factor=inv_sf)

    logger.info("STEP 06 complete.\n")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Step 06: BigWig generation")
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
