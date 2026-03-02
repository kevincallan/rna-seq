#!/usr/bin/env python3
"""
Step 01 -- Prepare samples.

Parse metadata CSV, apply subset filters, build deterministic design
table, create FASTQ symlinks, and write samples.tsv.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Any, Dict

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.metadata import (
    apply_subset_filters,
    build_design_table,
    create_symlinks,
    parse_metadata,
    write_sample_description,
    write_samples_tsv,
)
from src.utils import (
    get_run_id,
    load_config,
    resolve_results_dir,
    resolve_work_dir,
    setup_logging,
)

logger = logging.getLogger(__name__)


def main(cfg: Dict[str, Any], subset_name: str | None = None) -> None:
    """Execute step 01.

    Parameters
    ----------
    cfg : dict
        Loaded (and enriched) config.
    subset_name : str, optional
        Override the subset filter key.
    """
    logger.info("=" * 60)
    logger.info("STEP 01: Prepare samples")
    logger.info("=" * 60)

    run_id = cfg["_run_id"]
    results_dir = Path(cfg["_results_dir"])
    work_dir = resolve_work_dir(cfg)

    # 1. Parse metadata
    rows = parse_metadata(cfg)

    # 2. Apply subset filters
    filtered = apply_subset_filters(rows, cfg, subset_name=subset_name)

    # 3. Build design table with deterministic replicate numbering
    samples = build_design_table(filtered, cfg)

    # 4. Write samples.tsv
    samples_tsv = results_dir / "samples.tsv"
    write_samples_tsv(samples, samples_tsv)

    # 5. Create FASTQ symlinks
    create_symlinks(samples, work_dir)

    # 6. Write sample_description.txt for DESeq2
    desc_path = results_dir / "sample_description.txt"
    write_sample_description(samples, desc_path)

    # Store sample list in config for downstream steps
    cfg["_samples"] = samples
    cfg["_samples_tsv"] = str(samples_tsv)

    logger.info("STEP 01 complete: %d samples prepared.\n", len(samples))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Step 01: Prepare samples")
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    parser.add_argument("--run-id", default=None, help="Override run ID")
    parser.add_argument("--subset", default=None, help="Subset filter name")
    args = parser.parse_args()

    cfg = load_config(args.config)
    run_id = args.run_id or get_run_id(cfg)
    setup_logging(cfg["project"]["logs_dir"], run_id)

    # Minimal enrichment for standalone run
    cfg["_run_id"] = run_id
    cfg["_results_dir"] = str(resolve_results_dir(cfg, run_id))
    cfg["_work_dir"] = str(resolve_work_dir(cfg))

    main(cfg, subset_name=args.subset)
