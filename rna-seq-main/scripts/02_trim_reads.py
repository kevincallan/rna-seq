#!/usr/bin/env python3
"""
Step 02 -- Trim reads.

For each enabled trimming method, produce trimmed (or untrimmed) FASTQs
into ``work/trimmed/<method>/``.  Records per-sample trimming metrics.
"""

from __future__ import annotations

import csv
import logging
import os
import sys
from pathlib import Path
from typing import Any, Dict, List

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.metadata import Sample, read_samples_tsv
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
# Per-method trimming functions
# ---------------------------------------------------------------------------

def _trim_none(sample: Sample, work_dir: Path) -> Dict[str, Any]:
    """No trimming -- symlink raw FASTQs into the 'none' directory."""
    out_dir = work_dir / "trimmed" / "none"
    ensure_dirs(out_dir)
    link_dir = work_dir / "fastq_links"

    metrics: Dict[str, Any] = {"method": "none", "sample": sample.sample_name}

    if sample.layout == "paired":
        src1 = link_dir / f"{sample.sample_name}_1.fastq.gz"
        src2 = link_dir / f"{sample.sample_name}_2.fastq.gz"
        dst1 = out_dir / f"{sample.sample_name}_1.fastq.gz"
        dst2 = out_dir / f"{sample.sample_name}_2.fastq.gz"
        for src, dst in [(src1, dst1), (src2, dst2)]:
            if dst.exists() or dst.is_symlink():
                dst.unlink()
            os.symlink(src.resolve(), dst)
    else:
        src = link_dir / f"{sample.sample_name}.fastq.gz"
        dst = out_dir / f"{sample.sample_name}.fastq.gz"
        if dst.exists() or dst.is_symlink():
            dst.unlink()
        os.symlink(src.resolve(), dst)

    metrics["reads_in"] = "N/A"
    metrics["reads_out"] = "N/A"
    metrics["bases_trimmed"] = "N/A"
    return metrics


def _trim_cutadapt(
    sample: Sample,
    work_dir: Path,
    cfg: Dict[str, Any],
    results_dir: Path,
) -> Dict[str, Any]:
    """Trim with cutadapt."""
    params = cfg["trimming"]["cutadapt"]
    out_dir = work_dir / "trimmed" / "cutadapt"
    log_dir = results_dir / "cutadapt" / "qc"
    ensure_dirs(out_dir, log_dir)

    link_dir = work_dir / "fastq_links"
    metrics: Dict[str, Any] = {"method": "cutadapt", "sample": sample.sample_name}

    exe = cfg["tools"].get("cutadapt", "cutadapt")

    if sample.layout == "paired":
        in1 = link_dir / f"{sample.sample_name}_1.fastq.gz"
        in2 = link_dir / f"{sample.sample_name}_2.fastq.gz"
        out1 = out_dir / f"{sample.sample_name}_1.fastq.gz"
        out2 = out_dir / f"{sample.sample_name}_2.fastq.gz"

        cmd = [
            exe,
            "-a", params.get("adapter_fwd", "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"),
            "-A", params.get("adapter_rev", "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"),
            "-q", str(params.get("quality", 20)),
            "-m", str(params.get("min_length", 25)),
            "-o", str(out1),
            "-p", str(out2),
            str(in1), str(in2),
        ]
    else:
        in1 = link_dir / f"{sample.sample_name}.fastq.gz"
        out1 = out_dir / f"{sample.sample_name}.fastq.gz"

        cmd = [
            exe,
            "-a", params.get("adapter_fwd", "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"),
            "-q", str(params.get("quality", 20)),
            "-m", str(params.get("min_length", 25)),
            "-o", str(out1),
            str(in1),
        ]

    # Add extra args
    extra = params.get("extra_args", "")
    if extra:
        cmd.extend(extra.split())

    result = run_cmd(cmd, description=f"cutadapt {sample.sample_name}")

    # Parse cutadapt log from stderr
    metrics["reads_in"] = _parse_cutadapt_metric(result.stderr, "Total read pairs processed")
    metrics["reads_out"] = _parse_cutadapt_metric(result.stderr, "Pairs written")
    metrics["bases_trimmed"] = "see_log"

    # Save log
    log_path = log_dir / f"{sample.sample_name}.cutadapt.log"
    log_path.write_text(result.stderr or "", encoding="utf-8")

    return metrics


def _parse_cutadapt_metric(text: str, label: str) -> str:
    """Extract a value from cutadapt log output."""
    if not text:
        return "N/A"
    for line in text.splitlines():
        if label in line:
            # e.g. "Total read pairs processed:          1,234"
            parts = line.split(":")
            if len(parts) >= 2:
                return parts[-1].strip().replace(",", "")
    return "N/A"


def _trim_fastp(
    sample: Sample,
    work_dir: Path,
    cfg: Dict[str, Any],
    results_dir: Path,
) -> Dict[str, Any]:
    """Trim with fastp."""
    params = cfg["trimming"]["fastp"]
    out_dir = work_dir / "trimmed" / "fastp"
    log_dir = results_dir / "fastp" / "qc"
    ensure_dirs(out_dir, log_dir)

    link_dir = work_dir / "fastq_links"
    metrics: Dict[str, Any] = {"method": "fastp", "sample": sample.sample_name}

    exe = cfg["tools"].get("fastp", "fastp")

    json_report = log_dir / f"{sample.sample_name}.fastp.json"
    html_report = log_dir / f"{sample.sample_name}.fastp.html"

    if sample.layout == "paired":
        in1 = link_dir / f"{sample.sample_name}_1.fastq.gz"
        in2 = link_dir / f"{sample.sample_name}_2.fastq.gz"
        out1 = out_dir / f"{sample.sample_name}_1.fastq.gz"
        out2 = out_dir / f"{sample.sample_name}_2.fastq.gz"

        cmd = [
            exe,
            "-i", str(in1),
            "-I", str(in2),
            "-o", str(out1),
            "-O", str(out2),
            "-q", str(params.get("quality", 20)),
            "-l", str(params.get("min_length", 25)),
            "-j", str(json_report),
            "-h", str(html_report),
            "-w", str(cfg["project"].get("threads", 4)),
        ]
    else:
        in1 = link_dir / f"{sample.sample_name}.fastq.gz"
        out1 = out_dir / f"{sample.sample_name}.fastq.gz"

        cmd = [
            exe,
            "-i", str(in1),
            "-o", str(out1),
            "-q", str(params.get("quality", 20)),
            "-l", str(params.get("min_length", 25)),
            "-j", str(json_report),
            "-h", str(html_report),
            "-w", str(cfg["project"].get("threads", 4)),
        ]

    if not params.get("detect_adapter", True):
        cmd.append("--disable_adapter_trimming")

    extra = params.get("extra_args", "")
    if extra:
        cmd.extend(extra.split())

    result = run_cmd(cmd, description=f"fastp {sample.sample_name}")

    # Parse fastp JSON for metrics
    import json as json_mod
    if json_report.exists():
        with open(json_report) as fh:
            rpt = json_mod.load(fh)
        metrics["reads_in"] = rpt.get("summary", {}).get("before_filtering", {}).get("total_reads", "N/A")
        metrics["reads_out"] = rpt.get("summary", {}).get("after_filtering", {}).get("total_reads", "N/A")
        metrics["bases_trimmed"] = "see_json"
    else:
        metrics["reads_in"] = "N/A"
        metrics["reads_out"] = "N/A"
        metrics["bases_trimmed"] = "N/A"

    return metrics


def _trim_trimmomatic(
    sample: Sample,
    work_dir: Path,
    cfg: Dict[str, Any],
    results_dir: Path,
) -> Dict[str, Any]:
    """Trim with Trimmomatic."""
    params = cfg["trimming"]["trimmomatic"]
    out_dir = work_dir / "trimmed" / "trimmomatic"
    log_dir = results_dir / "trimmomatic" / "qc"
    ensure_dirs(out_dir, log_dir)

    link_dir = work_dir / "fastq_links"
    metrics: Dict[str, Any] = {"method": "trimmomatic", "sample": sample.sample_name}

    exe = cfg["tools"].get("trimmomatic", "trimmomatic")
    threads = cfg["project"].get("threads", 4)

    if sample.layout == "paired":
        in1 = link_dir / f"{sample.sample_name}_1.fastq.gz"
        in2 = link_dir / f"{sample.sample_name}_2.fastq.gz"
        out1_p = out_dir / f"{sample.sample_name}_1.fastq.gz"
        out1_u = out_dir / f"{sample.sample_name}_1.unpaired.fastq.gz"
        out2_p = out_dir / f"{sample.sample_name}_2.fastq.gz"
        out2_u = out_dir / f"{sample.sample_name}_2.unpaired.fastq.gz"
        log_file = log_dir / f"{sample.sample_name}.trimmomatic.log"

        cmd = [
            exe, "PE",
            "-threads", str(threads),
            "-trimlog", str(log_file),
            str(in1), str(in2),
            str(out1_p), str(out1_u),
            str(out2_p), str(out2_u),
            f"ILLUMINACLIP:{params.get('adapter_file', 'TruSeq3-PE-2.fa')}:2:30:10",
            f"LEADING:{params.get('leading', 3)}",
            f"TRAILING:{params.get('trailing', 3)}",
            f"SLIDINGWINDOW:{params.get('slidingwindow', '4:20')}",
            f"MINLEN:{params.get('min_length', 25)}",
        ]
    else:
        in1 = link_dir / f"{sample.sample_name}.fastq.gz"
        out1 = out_dir / f"{sample.sample_name}.fastq.gz"
        log_file = log_dir / f"{sample.sample_name}.trimmomatic.log"

        cmd = [
            exe, "SE",
            "-threads", str(threads),
            "-trimlog", str(log_file),
            str(in1), str(out1),
            f"ILLUMINACLIP:{params.get('adapter_file', 'TruSeq3-SE.fa')}:2:30:10",
            f"LEADING:{params.get('leading', 3)}",
            f"TRAILING:{params.get('trailing', 3)}",
            f"SLIDINGWINDOW:{params.get('slidingwindow', '4:20')}",
            f"MINLEN:{params.get('min_length', 25)}",
        ]

    extra = params.get("extra_args", "")
    if extra:
        cmd.extend(extra.split())

    result = run_cmd(cmd, description=f"trimmomatic {sample.sample_name}")

    metrics["reads_in"] = "see_log"
    metrics["reads_out"] = "see_log"
    metrics["bases_trimmed"] = "see_log"
    return metrics


# ---------------------------------------------------------------------------
# Dispatcher
# ---------------------------------------------------------------------------

_TRIMMERS = {
    "none": _trim_none,
    "cutadapt": _trim_cutadapt,
    "fastp": _trim_fastp,
    "trimmomatic": _trim_trimmomatic,
}


def main(cfg: Dict[str, Any], methods_override: List[str] | None = None) -> None:
    """Execute step 02 for all enabled trimming methods."""
    logger.info("=" * 60)
    logger.info("STEP 02: Trim reads")
    logger.info("=" * 60)

    work_dir = resolve_work_dir(cfg)
    results_dir = Path(cfg["_results_dir"])

    # Load samples
    samples_tsv = Path(cfg.get("_samples_tsv", results_dir / "samples.tsv"))
    samples = read_samples_tsv(samples_tsv)
    # Restore FASTQ paths from links
    link_dir = work_dir / "fastq_links"
    for s in samples:
        if s.layout == "paired":
            s.r1 = str(link_dir / f"{s.sample_name}_1.fastq.gz")
            s.r2 = str(link_dir / f"{s.sample_name}_2.fastq.gz")
        else:
            s.r1 = str(link_dir / f"{s.sample_name}.fastq.gz")

    methods = methods_override or get_enabled_methods(cfg)
    logger.info("Trimming methods: %s", methods)

    all_metrics: List[Dict[str, Any]] = []

    for method in methods:
        logger.info("--- Trimming method: %s ---", method)
        trim_fn = _TRIMMERS.get(method)
        if trim_fn is None:
            logger.error("Unknown trimming method: %s", method)
            continue

        for sample in samples:
            logger.info("  Processing %s ...", sample.sample_name)
            if method == "none":
                m = trim_fn(sample, work_dir)
            else:
                m = trim_fn(sample, work_dir, cfg, results_dir)
            all_metrics.append(m)

    # Write trimming summary
    summary_path = results_dir / "trimming_summary.tsv"
    _write_trim_summary(all_metrics, summary_path)

    logger.info("STEP 02 complete.\n")


def _write_trim_summary(metrics: List[Dict[str, Any]], path: Path) -> None:
    """Write trimming metrics to TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["method", "sample", "reads_in", "reads_out", "bases_trimmed"]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        for m in metrics:
            writer.writerow(m)
    logger.info("Trimming summary -> %s", path)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Step 02: Trim reads")
    parser.add_argument("--config", required=True)
    parser.add_argument("--run-id", default=None)
    parser.add_argument("--methods", nargs="*", default=None,
                        help="Override trimming methods")
    args = parser.parse_args()

    cfg = load_config(args.config)
    run_id = args.run_id or get_run_id(cfg)
    setup_logging(cfg["project"]["logs_dir"], run_id)
    cfg["_run_id"] = run_id
    cfg["_results_dir"] = str(resolve_results_dir(cfg, run_id))
    cfg["_work_dir"] = str(resolve_work_dir(cfg))

    main(cfg, methods_override=args.methods)
