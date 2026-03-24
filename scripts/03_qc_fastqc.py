#!/usr/bin/env python3
"""
Step 03 -- Run FastQC.

Runs FastQC on raw and trimmed FASTQs for each enabled trimming method.
Output goes to ``results/<run_id>/<method>/qc/fastqc/``.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Any, Dict, List

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.metadata import read_samples_tsv
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


def run_fastqc(
    fastq_files: List[Path],
    out_dir: Path,
    threads: int,
    exe: str = "fastqc",
) -> None:
    """Run FastQC on a list of FASTQ files."""
    if not fastq_files:
        logger.warning("No FASTQ files provided for FastQC")
        return

    ensure_dirs(out_dir)

    cmd = [
        exe,
        "-t", str(threads),
        "-o", str(out_dir),
        "--noextract",
    ] + [str(f) for f in fastq_files]

    run_cmd(cmd, description=f"FastQC -> {out_dir.name}")


def main(cfg: Dict[str, Any], methods_override: List[str] | None = None) -> None:
    """Execute step 03."""
    logger.info("=" * 60)
    logger.info("STEP 03: FastQC quality control")
    logger.info("=" * 60)

    work_dir = resolve_work_dir(cfg)
    results_dir = Path(cfg["_results_dir"])
    threads = cfg["project"].get("threads", 4)
    exe = cfg["tools"].get("fastqc", "fastqc")

    samples_tsv = Path(cfg.get("_samples_tsv", results_dir / "samples.tsv"))
    samples = read_samples_tsv(samples_tsv)

    methods = get_effective_trim_methods(cfg, methods_override)
    trim_cfg = get_trim_config_summary(cfg)
    logger.info(
        "Trim config: primary=%s compare_methods=%s effective=%s",
        trim_cfg["primary_method"],
        trim_cfg["compare_methods"],
        ",".join(methods),
    )

    # --- FastQC on raw (linked) FASTQs --------------------------------------
    logger.info("Running FastQC on raw reads...")
    raw_dir = results_dir / "raw_qc" / "fastqc"
    raw_fastqs: List[Path] = []
    link_dir = work_dir / "fastq_links"
    for s in samples:
        if s.layout == "paired":
            raw_fastqs.append(link_dir / f"{s.sample_name}_1.fastq.gz")
            raw_fastqs.append(link_dir / f"{s.sample_name}_2.fastq.gz")
        else:
            raw_fastqs.append(link_dir / f"{s.sample_name}.fastq.gz")

    existing = [f for f in raw_fastqs if f.exists()]
    if existing:
        run_fastqc(existing, raw_dir, threads, exe)
    else:
        logger.warning("No raw FASTQ symlinks found in %s", link_dir)

    # --- FastQC on trimmed FASTQs per method --------------------------------
    for method in methods:
        logger.info("Running FastQC on %s-trimmed reads...", method)
        trim_dir = work_dir / "trimmed" / method
        qc_dir = results_dir / method / "qc" / "fastqc"

        trimmed_fastqs: List[Path] = []
        for s in samples:
            if s.layout == "paired":
                trimmed_fastqs.append(trim_dir / f"{s.sample_name}_1.fastq.gz")
                trimmed_fastqs.append(trim_dir / f"{s.sample_name}_2.fastq.gz")
            else:
                trimmed_fastqs.append(trim_dir / f"{s.sample_name}.fastq.gz")

        existing_trimmed = [f for f in trimmed_fastqs if f.exists()]
        if existing_trimmed:
            run_fastqc(existing_trimmed, qc_dir, threads, exe)
        else:
            logger.warning("No trimmed FASTQs found for method '%s' in %s",
                           method, trim_dir)

    logger.info("STEP 03 complete.\n")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Step 03: FastQC")
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
