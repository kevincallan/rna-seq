#!/usr/bin/env python3
"""
Colab + Google Drive setup: create output dirs, optionally clone repo, and point
config results_dir/work_dir to Drive so pipeline outputs persist across sessions.

Run from a Colab cell after mounting Drive, e.g.:

  from google.colab import drive
  drive.mount("/content/drive")

  %cd /content/drive/MyDrive
  !git clone https://github.com/kevincallan/rna-seq.git
  %cd rna-seq
  !bash scripts/setup_colab.sh

  !./py scripts/colab_drive_setup.py

Then run the pipeline as usual; results and work will be on Drive.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

# Default: Colab mounts Drive at /content/drive/MyDrive (no space in "MyDrive")
DRIVE_BASE = os.environ.get("COLAB_DRIVE_BASE", "/content/drive/MyDrive")
RESULTS_DIR = Path(DRIVE_BASE) / "rna_seq_results"
WORK_DIR = Path(DRIVE_BASE) / "rna_seq_work"


def main() -> None:
    # Repo root = parent of scripts/
    repo_root = Path(__file__).resolve().parent.parent
    config_path = repo_root / "config" / "config_colab.yaml"

    if not config_path.exists():
        print("Config not found at", config_path, file=sys.stderr)
        print("Run setup_colab.sh first (from the repo root).", file=sys.stderr)
        sys.exit(1)

    # Create Drive dirs
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {RESULTS_DIR}")
    print(f"Work directory:   {WORK_DIR}")

    # Patch config: set results_dir, work_dir, logs_dir to Drive paths (line-by-line to keep rest of file)
    logs_dir = Path(DRIVE_BASE) / "rna_seq_logs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    with open(config_path) as f:
        lines = f.readlines()
    out = []
    for line in lines:
        if line.strip().startswith("results_dir:"):
            out.append(f'  results_dir: "{RESULTS_DIR}"\n')
        elif line.strip().startswith("work_dir:"):
            out.append(f'  work_dir: "{WORK_DIR}"\n')
        elif line.strip().startswith("logs_dir:"):
            out.append(f'  logs_dir: "{logs_dir}"\n')
        else:
            out.append(line)
    with open(config_path, "w") as f:
        f.writelines(out)

    print("Config updated: results_dir, work_dir, logs_dir now point to Drive.")
    print("Run: ./py scripts/run_pipeline.py --config config/config_colab.yaml run")


if __name__ == "__main__":
    main()
