#!/usr/bin/env python3
"""
Step 06 -- Generate BigWig coverage tracks.

Produces per-sample BigWig files for each trimming method:
  - CPM-normalised with MAPQ filter (unique mappers only)
  - Optionally DESeq2 size-factor scaled (if size factors file exists)
"""

from __future__ import annotations

import logging
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


from src.analysis_unit import (
    build_mapping_units_with_bams,
    read_selected_analysis,
    read_selected_visualisation,
)


def make_bigwig(
    bam: Path,
    out_bw: Path,
    cfg: Dict[str, Any],
    scale_factor: Optional[float] = None,
) -> None:
    """Run bamCoverage to create a BigWig file.

    Skips generation if ``out_bw`` already exists unless ``cfg["_force"]``
    is set.

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
    if out_bw.exists() and not cfg.get("_force", False):
        logger.info("  Skipping %s (already exists; use --force to regenerate)", out_bw.name)
        return

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


def _find_size_factors(
    results_dir: Path,
    method: str,
    mapper: str,
    mapper_opt: str,
    count_opt: str,
) -> Dict[str, float]:
    """Search for size_factors.tsv in count-option-aware then legacy paths."""
    from src.analysis_unit import resolve_de_dir

    new_path = resolve_de_dir(results_dir, method, mapper, mapper_opt, count_opt) / "size_factors.tsv"
    if new_path.exists():
        return load_size_factors(new_path)

    legacy = results_dir / method / "deseq2" / mapper / mapper_opt / "size_factors.tsv"
    if legacy.exists():
        return load_size_factors(legacy)

    return {}


def main(cfg: Dict[str, Any], methods_override: List[str] | None = None) -> None:
    """Execute BigWig generation."""
    logger.info("=" * 60)
    logger.info("STEP: BigWig generation")
    logger.info("=" * 60)

    results_dir = Path(cfg["_results_dir"])
    bw_cfg = cfg.get("bigwig", {})

    if not bw_cfg.get("enabled", True):
        logger.info("BigWig generation disabled in config (bigwig.enabled=false). Skipping.")
        return

    mode = bw_cfg.get("mode", "all_units")
    samples_tsv = Path(cfg.get("_samples_tsv", results_dir / "samples.tsv"))
    samples = read_samples_tsv(samples_tsv)
    methods = methods_override or get_enabled_methods(cfg)

    mapping_units = build_mapping_units_with_bams(results_dir, samples, methods,
                                                  prefer_filtered=True)

    selected = read_selected_visualisation(results_dir) if mode == "selected_only" else None
    if mode == "selected_only" and selected is None:
        # Backward compatibility for runs that only produced selected_analysis.tsv.
        selected = read_selected_analysis(results_dir)
    if mode == "selected_only" and selected is None:
        logger.warning(
            "BigWig mode is 'selected_only' but no selected_visualisation.tsv/selected_analysis.tsv found. "
            "Falling back to all_units mode."
        )
        mode = "all_units"

    # Determine which count_opt to use for size factors
    option_sets = cfg.get("featurecounts", {}).get("option_sets", {"default": {}})
    sf_count_opt = list(option_sets.keys())[0] if option_sets else "default"
    if selected:
        sf_count_opt = selected.count_option_set

    for unit in mapping_units:
        method = unit["method"]
        mapper = unit["mapper"]
        mapper_opt = unit["mapper_option_set"]

        if mode == "selected_only" and selected is not None:
            if (method, mapper, mapper_opt) != selected.mapping_key:
                logger.info(
                    "Skipping BigWig for %s/%s/%s (not the selected branch)",
                    method, mapper, mapper_opt,
                )
                continue

        logger.info(
            "--- BigWigs for trim=%s mapper=%s mapper_option=%s ---",
            method, mapper, mapper_opt,
        )
        bw_dir = results_dir / method / "bigwig" / mapper / mapper_opt
        ensure_dirs(bw_dir)

        size_factors: Dict[str, float] = {}
        if bw_cfg.get("use_deseq2_sizefactors", True):
            size_factors = _find_size_factors(
                results_dir, method, mapper, mapper_opt, sf_count_opt
            )
        use_sf = bool(size_factors)

        for s in samples:
            bam = unit["sample_to_bam"].get(s.sample_name)
            if bam is None or not bam.exists():
                logger.warning("BAM not found for %s (%s/%s/%s)",
                               s.sample_name, method, mapper, mapper_opt)
                continue

            bw_cpm = bw_dir / f"{s.sample_name}.CPM.MAPQ255.bw"
            logger.info("  %s -> CPM BigWig", s.sample_name)
            make_bigwig(bam, bw_cpm, cfg)

            if use_sf and s.sample_name in size_factors:
                sf = size_factors[s.sample_name]
                inv_sf = 1.0 / sf if sf != 0 else 1.0
                bw_deseq = bw_dir / f"{s.sample_name}.DESeq2scaled.MAPQ255.bw"
                logger.info("  %s -> DESeq2-scaled BigWig (sf=%.4f, 1/sf=%.4f)",
                            s.sample_name, sf, inv_sf)
                make_bigwig(bam, bw_deseq, cfg, scale_factor=inv_sf)

    logger.info("BigWig generation complete.\n")


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
