#!/usr/bin/env python3
"""
Optional helper -- Download FASTQs from SRA.

Downloads FASTQ files for a list of SRA Run IDs using ``fasterq-dump``
(from the SRA Toolkit).  This is optional and may not be available on
exam machines without internet access.

Usage:
    ./py scripts/download_sra.py --ids SRR925874 SRR925875 --outdir /data/fastqs
    ./py scripts/download_sra.py --config config/config.yaml
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from pathlib import Path
from typing import List

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.utils import ensure_dirs, load_config, run_cmd, setup_logging

logger = logging.getLogger(__name__)


def download_run(
    run_id: str,
    out_dir: Path,
    threads: int = 4,
) -> None:
    """Download a single SRA run using fasterq-dump."""
    ensure_dirs(out_dir)

    cmd = [
        "fasterq-dump",
        run_id,
        "-O", str(out_dir),
        "-e", str(threads),
        "--split-files",
        "--progress",
    ]

    run_cmd(cmd, description=f"fasterq-dump {run_id}")

    # Compress with pigz/gzip
    for fq in out_dir.glob(f"{run_id}*.fastq"):
        logger.info("Compressing %s ...", fq.name)
        try:
            run_cmd(["pigz", "-p", str(threads), str(fq)],
                    description=f"pigz {fq.name}")
        except FileNotFoundError:
            run_cmd(["gzip", str(fq)], description=f"gzip {fq.name}")


def get_run_ids_from_config(config_path: str) -> List[str]:
    """Extract Run IDs from the metadata CSV referenced in config."""
    cfg = load_config(config_path)
    csv_path = cfg["data"]["metadata_csv"]
    run_col = cfg["column_mapping"]["run_id_col"]

    ids: List[str] = []
    with open(csv_path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            run_id = row.get(run_col, "").strip()
            if run_id:
                ids.append(run_id)

    return ids


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Download FASTQs from SRA (optional helper)"
    )
    parser.add_argument("--ids", nargs="*", help="SRA Run IDs to download")
    parser.add_argument("--config", help="Config file to extract Run IDs from metadata")
    parser.add_argument("--outdir", default="fastqs", help="Output directory")
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    setup_logging("logs", "download_sra")

    run_ids: List[str] = []
    if args.ids:
        run_ids = args.ids
    elif args.config:
        run_ids = get_run_ids_from_config(args.config)
    else:
        parser.error("Provide --ids or --config")

    out_dir = Path(args.outdir)
    logger.info("Downloading %d runs to %s", len(run_ids), out_dir)

    for rid in run_ids:
        logger.info("--- Downloading %s ---", rid)
        try:
            download_run(rid, out_dir, args.threads)
        except Exception as exc:
            logger.error("Failed to download %s: %s", rid, exc)
            continue

    logger.info("Download complete.")


if __name__ == "__main__":
    main()
