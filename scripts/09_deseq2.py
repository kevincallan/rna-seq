#!/usr/bin/env python3
"""
Step 09 -- Differential expression with DESeq2.

Runs DESeq2 for each trimming method and each contrast defined in
config.  Supports two backends:

  - ``rscript``: calls ``deseq2_run.R`` via Rscript (default)
  - ``wrapper``: calls ``DESeq2_wrapper`` (available on course server)

Produces DE results tables, PCA/MA plots, size factors, and DE summaries.
"""

from __future__ import annotations

import csv
import logging
import sys
from pathlib import Path
from typing import Any, Dict, List

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.metadata import read_samples_tsv, write_sample_description
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

# Path to the R script shipped with this pipeline
DESEQ2_R_SCRIPT = Path(__file__).resolve().parent / "deseq2_run.R"


# ---------------------------------------------------------------------------
# DESeq2 via custom R script
# ---------------------------------------------------------------------------

def run_deseq2_rscript(
    count_matrix: Path,
    sample_desc: Path,
    contrast_name: str,
    numerator: str,
    denominator: str,
    out_dir: Path,
    cfg: Dict[str, Any],
) -> None:
    """Run DESeq2 analysis using deseq2_run.R."""
    rscript = cfg["tools"].get("rscript", "Rscript")
    deseq_cfg = cfg.get("deseq2", {})
    fdr = deseq_cfg.get("fdr_threshold", 0.05)
    lfc = deseq_cfg.get("lfc_threshold", 0.0)
    ref_level = deseq_cfg.get("reference_level", "")

    ensure_dirs(out_dir)

    cmd = [
        rscript, str(DESEQ2_R_SCRIPT),
        "--counts", str(count_matrix),
        "--samples", str(sample_desc),
        "--contrast_name", contrast_name,
        "--numerator", numerator,
        "--denominator", denominator,
        "--fdr", str(fdr),
        "--lfc", str(lfc),
        "--outdir", str(out_dir),
    ]

    if ref_level:
        cmd.extend(["--reference_level", ref_level])

    run_cmd(cmd, description=f"DESeq2 {contrast_name}")


# ---------------------------------------------------------------------------
# DESeq2 via DESeq2_wrapper (course server)
# ---------------------------------------------------------------------------

def run_deseq2_wrapper(
    count_matrix: Path,
    sample_desc: Path,
    out_dir: Path,
    cfg: Dict[str, Any],
) -> None:
    """Run DESeq2 using the DESeq2_wrapper tool."""
    ensure_dirs(out_dir)

    cmd = ["DESeq2_wrapper", str(count_matrix), str(sample_desc)]
    run_cmd(cmd, description="DESeq2_wrapper", cwd=out_dir)


# ---------------------------------------------------------------------------
# Summary helpers
# ---------------------------------------------------------------------------

def count_degs(results_path: Path, fdr: float = 0.05) -> Dict[str, int]:
    """Count DEGs from a DESeq2 results TSV."""
    total = 0
    sig_up = 0
    sig_down = 0

    if not results_path.exists():
        return {"total_tested": 0, "sig_up": 0, "sig_down": 0, "sig_total": 0}

    with open(results_path, encoding="utf-8") as fh:
        header = fh.readline().strip().split("\t")
        # Find padj and log2FoldChange column indices
        try:
            padj_idx = header.index("padj")
        except ValueError:
            # Try without header row name
            padj_idx = -1
            for i, h in enumerate(header):
                if "padj" in h.lower():
                    padj_idx = i
                    break

        try:
            lfc_idx = header.index("log2FoldChange")
        except ValueError:
            lfc_idx = -1
            for i, h in enumerate(header):
                if "log2foldchange" in h.lower():
                    lfc_idx = i
                    break

        if padj_idx < 0:
            return {"total_tested": 0, "sig_up": 0, "sig_down": 0, "sig_total": 0}

        for line in fh:
            parts = line.strip().split("\t")
            total += 1
            try:
                padj = float(parts[padj_idx])
            except (ValueError, IndexError):
                continue

            if padj < fdr:
                if lfc_idx >= 0:
                    try:
                        lfc = float(parts[lfc_idx])
                        if lfc > 0:
                            sig_up += 1
                        else:
                            sig_down += 1
                    except (ValueError, IndexError):
                        sig_up += 1  # count as significant without direction
                else:
                    sig_up += 1

    return {
        "total_tested": total,
        "sig_up": sig_up,
        "sig_down": sig_down,
        "sig_total": sig_up + sig_down,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(cfg: Dict[str, Any], methods_override: List[str] | None = None) -> None:
    """Execute step 09."""
    logger.info("=" * 60)
    logger.info("STEP 09: DESeq2 differential expression")
    logger.info("=" * 60)

    results_dir = Path(cfg["_results_dir"])
    deseq_cfg = cfg.get("deseq2", {})
    method_backend = deseq_cfg.get("method", "rscript")
    fdr = deseq_cfg.get("fdr_threshold", 0.05)
    comparisons = cfg.get("comparisons", [])

    samples_tsv = Path(cfg.get("_samples_tsv", results_dir / "samples.tsv"))
    samples = read_samples_tsv(samples_tsv)
    methods = methods_override or get_enabled_methods(cfg)
    option_sets = cfg["featurecounts"].get("option_sets", {"default": {}})

    # Use the first option set for DE by default (typically "default")
    primary_opt = list(option_sets.keys())[0]

    all_de_stats: List[Dict[str, Any]] = []

    for method in methods:
        logger.info("--- DESeq2 for method: %s ---", method)
        fc_dir = results_dir / method / "featurecounts"
        de_dir = results_dir / method / "deseq2"
        ensure_dirs(de_dir)

        # Count matrix
        count_matrix = fc_dir / f"clean_matrix_{primary_opt}.tsv"
        if not count_matrix.exists():
            logger.error("Count matrix not found: %s", count_matrix)
            continue

        # Write sample description for this method
        sample_desc = de_dir / "sample_description.txt"
        write_sample_description(samples, sample_desc)

        if method_backend == "wrapper":
            # DESeq2_wrapper runs a single default contrast
            logger.info("  Running DESeq2_wrapper...")
            run_deseq2_wrapper(count_matrix, sample_desc, de_dir, cfg)

            # Try to parse results
            de_all = de_dir / "DESeq2.de_all.tsv"
            if de_all.exists():
                stats = count_degs(de_all, fdr)
                stats["method"] = method
                stats["contrast"] = "auto"
                all_de_stats.append(stats)

        else:
            # Custom R script -- run each contrast
            for contrast in comparisons:
                cname = contrast["name"]
                num = contrast["numerator"]
                den = contrast["denominator"]
                contrast_dir = de_dir / cname

                logger.info("  Contrast: %s (%s vs %s)", cname, num, den)
                run_deseq2_rscript(
                    count_matrix, sample_desc,
                    cname, num, den,
                    contrast_dir, cfg,
                )

                # Parse DE results
                de_all = contrast_dir / "de_all.tsv"
                stats = count_degs(de_all, fdr)
                stats["method"] = method
                stats["contrast"] = cname
                all_de_stats.append(stats)

    # Write DE summary
    de_summary_path = results_dir / "de_summary.tsv"
    if all_de_stats:
        with open(de_summary_path, "w", newline="", encoding="utf-8") as fh:
            fields = ["method", "contrast", "total_tested", "sig_up", "sig_down", "sig_total"]
            writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                                    extrasaction="ignore")
            writer.writeheader()
            for row in all_de_stats:
                writer.writerow(row)
        logger.info("DE summary -> %s", de_summary_path)

    logger.info("STEP 09 complete.\n")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Step 09: DESeq2")
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
