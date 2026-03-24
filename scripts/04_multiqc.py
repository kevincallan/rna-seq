#!/usr/bin/env python3
"""
Step 04 -- Run MultiQC.

Generates a MultiQC report per trimming method (combining FastQC, STAR,
featureCounts, etc.) plus a combined index linking to each report.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Any, Dict, List

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.utils import (
    ensure_dirs,
    get_effective_trim_methods,
    get_trim_config_summary,
    get_run_id,
    load_config,
    resolve_results_dir,
    resolve_work_dir,
    run_cmd,
    setup_logging,
)

logger = logging.getLogger(__name__)


def run_multiqc(
    search_dir: Path,
    out_dir: Path,
    title: str,
    exe: str = "multiqc",
) -> Path:
    """Run MultiQC over a directory and return the report path."""
    ensure_dirs(out_dir)

    cmd = [
        exe,
        str(search_dir),
        "-o", str(out_dir),
        "--title", title,
        "--force",
        "--no-data-dir",
    ]
    run_cmd(cmd, description=f"MultiQC: {title}")

    report = out_dir / "multiqc_report.html"
    return report


def main(cfg: Dict[str, Any], methods_override: List[str] | None = None) -> None:
    """Execute step 04."""
    logger.info("=" * 60)
    logger.info("STEP 04: MultiQC reports")
    logger.info("=" * 60)

    results_dir = Path(cfg["_results_dir"])
    exe = cfg["tools"].get("multiqc", "multiqc")
    methods = get_effective_trim_methods(cfg, methods_override)
    trim_cfg = get_trim_config_summary(cfg)
    logger.info(
        "Trim config: primary=%s compare_methods=%s effective=%s",
        trim_cfg["primary_method"],
        trim_cfg["compare_methods"],
        ",".join(methods),
    )
    project_name = cfg["project"]["name"]

    report_index: List[Dict[str, str]] = []

    # --- Raw QC report -------------------------------------------------------
    raw_qc_dir = results_dir / "raw_qc"
    if raw_qc_dir.exists():
        logger.info("MultiQC for raw reads...")
        rpt = run_multiqc(
            raw_qc_dir,
            results_dir / "reports" / "multiqc_raw",
            f"{project_name} - Raw Reads QC",
            exe,
        )
        report_index.append({"method": "raw", "report": str(rpt)})

    # --- Per-method reports --------------------------------------------------
    for method in methods:
        method_dir = results_dir / method
        if not method_dir.exists():
            logger.warning("No results directory for method '%s', skipping", method)
            continue

        logger.info("MultiQC for method '%s'...", method)
        rpt = run_multiqc(
            method_dir,
            results_dir / "reports" / f"multiqc_{method}",
            f"{project_name} - {method}",
            exe,
        )
        report_index.append({"method": method, "report": str(rpt)})

    # --- Write index table ---------------------------------------------------
    index_path = results_dir / "reports" / "multiqc_index.tsv"
    with open(index_path, "w", encoding="utf-8") as fh:
        fh.write("method\treport_path\n")
        for entry in report_index:
            fh.write(f"{entry['method']}\t{entry['report']}\n")

    logger.info("MultiQC index -> %s", index_path)
    logger.info("STEP 04 complete.\n")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Step 04: MultiQC")
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
