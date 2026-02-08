#!/usr/bin/env python3
"""
RNA-seq Pipeline Orchestrator
==============================

Runs the full pipeline or individual steps.  All configuration comes
from ``config/config.yaml`` (override with ``--config``).

Usage examples
--------------
# Full pipeline:
python scripts/run_pipeline.py --config config/config.yaml run

# Specific steps only:
python scripts/run_pipeline.py --config config/config.yaml run --steps 0 1 2 5

# Override trimming methods:
python scripts/run_pipeline.py --config config/config.yaml run --methods none cutadapt fastp

# Use a specific subset filter:
python scripts/run_pipeline.py --config config/config.yaml run --subset day3_wt_vs_tet1

# Single step:
python scripts/run_pipeline.py --config config/config.yaml run --steps 5
"""

from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

# Allow running from repository root
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.utils import (
    get_enabled_methods,
    get_run_id,
    load_config,
    resolve_results_dir,
    resolve_work_dir,
    setup_logging,
)

logger = logging.getLogger("pipeline")


# ---------------------------------------------------------------------------
# Step registry
# ---------------------------------------------------------------------------

STEPS = {
    0:  ("00_validate_env",     "Validate environment"),
    1:  ("01_prepare_samples",  "Prepare samples"),
    2:  ("02_trim_reads",       "Trim reads"),
    3:  ("03_qc_fastqc",       "FastQC"),
    4:  ("04_multiqc",         "MultiQC"),
    5:  ("05_map_star",        "STAR mapping"),
    6:  ("06_bigwig",          "BigWig generation"),
    7:  ("07_featurecounts",   "featureCounts"),
    8:  ("08_filter_matrix",   "Filter count matrix"),
    9:  ("09_deseq2",         "DESeq2"),
    10: ("10_compare_methods", "Compare methods"),
    11: ("11_make_report",     "Generate report"),
}

ALL_STEPS = sorted(STEPS.keys())


def import_step(step_num: int):
    """Dynamically import a step module."""
    module_name, _ = STEPS[step_num]
    # Import from scripts package
    import importlib
    mod = importlib.import_module(f"scripts.{module_name}")
    return mod


# ---------------------------------------------------------------------------
# Pipeline runner
# ---------------------------------------------------------------------------

def run_pipeline(
    cfg: Dict[str, Any],
    run_id: str,
    steps: Optional[List[int]] = None,
    methods_override: Optional[List[str]] = None,
    subset_name: Optional[str] = None,
) -> None:
    """Run the pipeline steps in order.

    Parameters
    ----------
    cfg : dict
        Loaded config.
    run_id : str
        Unique run identifier.
    steps : list of int, optional
        Which steps to run.  Default: all.
    methods_override : list of str, optional
        Override trimming methods from config.
    subset_name : str, optional
        Override subset filter from config.
    """
    steps = steps if steps is not None else ALL_STEPS

    logger.info("=" * 70)
    logger.info("RNA-seq Pipeline   |   Run ID: %s", run_id)
    logger.info("Steps to run: %s", steps)
    logger.info("=" * 70)

    t0 = time.time()

    for step_num in steps:
        if step_num not in STEPS:
            logger.error("Unknown step: %d (valid: %s)", step_num, ALL_STEPS)
            continue

        _, description = STEPS[step_num]
        logger.info("\n>>> STEP %d: %s", step_num, description)
        step_start = time.time()

        mod = import_step(step_num)

        try:
            if step_num == 0:
                cfg = mod.main(cfg, run_id)
            elif step_num == 1:
                mod.main(cfg, subset_name=subset_name)
            elif step_num in (2, 3, 4, 5, 6, 7, 8):
                mod.main(cfg, methods_override=methods_override)
            elif step_num == 9:
                mod.main(cfg, methods_override=methods_override)
            elif step_num == 10:
                mod.main(cfg, methods_override=methods_override)
            elif step_num == 11:
                mod.main(cfg)
            else:
                mod.main(cfg)

        except Exception as exc:
            logger.error(
                "STEP %d FAILED: %s: %s", step_num, type(exc).__name__, exc,
                exc_info=True,
            )
            raise

        elapsed = time.time() - step_start
        logger.info("<<< STEP %d completed in %.1f s", step_num, elapsed)

    total = time.time() - t0
    logger.info("\n" + "=" * 70)
    logger.info("Pipeline complete.  Total time: %.1f s", total)
    logger.info("Results directory: %s", cfg.get("_results_dir", "N/A"))
    logger.info("=" * 70)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser."""
    parser = argparse.ArgumentParser(
        description="RNA-seq Analysis Pipeline Orchestrator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    parser.add_argument(
        "--config", "-c",
        default="config/config.yaml",
        help="Path to config.yaml (default: config/config.yaml)",
    )
    parser.add_argument(
        "--run-id",
        default=None,
        help="Override run ID (default: auto-generated from project name + timestamp)",
    )

    sub = parser.add_subparsers(dest="command", help="Pipeline commands")

    # --- run subcommand -----------------------------------------------------
    run_p = sub.add_parser("run", help="Run the pipeline (all steps or selected)")
    run_p.add_argument(
        "--steps",
        nargs="*",
        type=int,
        default=None,
        help=f"Step numbers to run (default: all). Valid: {ALL_STEPS}",
    )
    run_p.add_argument(
        "--methods",
        nargs="*",
        default=None,
        help="Override trimming methods (e.g., none cutadapt fastp)",
    )
    run_p.add_argument(
        "--subset",
        default=None,
        help="Subset filter name from config (e.g., day3_wt_vs_tet1)",
    )

    # --- list subcommand ----------------------------------------------------
    sub.add_parser("list", help="List available pipeline steps")

    return parser


def main() -> None:
    """CLI entry point."""
    parser = build_parser()
    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    if args.command == "list":
        print("\nAvailable pipeline steps:")
        print("-" * 50)
        for num in ALL_STEPS:
            module, desc = STEPS[num]
            print(f"  {num:2d}  {desc:30s}  ({module}.py)")
        print()
        sys.exit(0)

    # --- Load config --------------------------------------------------------
    cfg = load_config(args.config)
    run_id = args.run_id or get_run_id(cfg)

    # Set up logging
    setup_logging(cfg["project"]["logs_dir"], run_id)

    # Pre-enrich config with runtime paths (step 00 will enrich further)
    cfg["_run_id"] = run_id
    cfg["_results_dir"] = str(resolve_results_dir(cfg, run_id))
    cfg["_work_dir"] = str(resolve_work_dir(cfg))

    # --- Run ----------------------------------------------------------------
    if args.command == "run":
        run_pipeline(
            cfg,
            run_id,
            steps=args.steps,
            methods_override=args.methods,
            subset_name=args.subset,
        )


if __name__ == "__main__":
    main()
