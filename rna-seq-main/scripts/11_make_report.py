#!/usr/bin/env python3
"""
Step 11 -- Generate final analysis report.

Produces ``results/<run_id>/reports/report.md`` (and optionally HTML)
combining all analysis summaries into a single document suitable for
submission or reference.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Any, Dict, List

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.metadata import read_samples_tsv
from src.reporting import MarkdownReport, render_html
from src.utils import (
    get_enabled_methods,
    get_run_id,
    load_config,
    resolve_results_dir,
    resolve_work_dir,
    setup_logging,
)

logger = logging.getLogger(__name__)


def main(cfg: Dict[str, Any]) -> None:
    """Execute step 11."""
    logger.info("=" * 60)
    logger.info("STEP 11: Generate report")
    logger.info("=" * 60)

    results_dir = Path(cfg["_results_dir"])
    methods = get_enabled_methods(cfg)
    comparisons = cfg.get("comparisons", [])
    project_name = cfg["project"]["name"]

    rpt = MarkdownReport(f"RNA-seq Analysis Report: {project_name}")

    # === Dataset description ================================================
    rpt.h2("Dataset Description")
    rpt.paragraph(f"**Project:** {project_name}")
    rpt.paragraph(f"**Metadata CSV:** `{cfg['data']['metadata_csv']}`")
    rpt.paragraph(f"**FASTQ directory:** `{cfg['data']['fastq_dir']}`")
    rpt.paragraph(f"**Genome index:** `{cfg['references']['genome_index']}`")
    rpt.paragraph(f"**GTF:** `{cfg['references']['gtf']}`")
    rpt.paragraph(f"**Layout:** {cfg['data'].get('layout', 'paired')}")
    rpt.paragraph(f"**Active subset:** {cfg.get('active_subset', 'default')}")

    # === Sample table =======================================================
    rpt.h2("Sample Table")
    samples_tsv = results_dir / "samples.tsv"
    if samples_tsv.exists():
        rpt.table_from_tsv(samples_tsv)
    else:
        rpt.paragraph("*(samples.tsv not found)*")

    # === QC summary =========================================================
    rpt.h2("Quality Control")
    rpt.paragraph("FastQC was run on raw and trimmed reads. MultiQC reports:")

    mqc_index = results_dir / "reports" / "multiqc_index.tsv"
    if mqc_index.exists():
        rpt.table_from_tsv(mqc_index)
    else:
        rpt.paragraph("*(MultiQC index not available)*")

    # === Trimming summary ===================================================
    rpt.h2("Trimming Summary")
    rpt.paragraph(f"Methods compared: {', '.join(methods)}")
    trim_summary = results_dir / "trimming_summary.tsv"
    if trim_summary.exists():
        rpt.table_from_tsv(trim_summary)
    else:
        rpt.paragraph("*(Trimming summary not available)*")

    # === Mapping summary ====================================================
    rpt.h2("Mapping Summary (STAR)")
    mapping_summary = results_dir / "mapping_summary.tsv"
    if mapping_summary.exists():
        rpt.table_from_tsv(mapping_summary)
    else:
        rpt.paragraph("*(Mapping summary not available)*")

    # === featureCounts summary ==============================================
    rpt.h2("featureCounts Summary")
    fc_summary = results_dir / "featurecounts_summary.tsv"
    if fc_summary.exists():
        rpt.table_from_tsv(fc_summary)
    else:
        rpt.paragraph("*(featureCounts summary not available)*")

    # === Filtering summary ==================================================
    rpt.h2("Count Matrix Filtering")
    filter_summary = results_dir / "filtering_summary.tsv"
    if filter_summary.exists():
        rpt.table_from_tsv(filter_summary)
    else:
        rpt.paragraph("*(Filtering summary not available)*")

    # === DE results per method / contrast ===================================
    rpt.h2("Differential Expression Results")
    de_summary = results_dir / "de_summary.tsv"
    if de_summary.exists():
        rpt.table_from_tsv(de_summary)
    else:
        rpt.paragraph("*(DE summary not available)*")

    # Per-method detail
    for method in methods:
        rpt.h3(f"Method: {method}")

        for contrast in comparisons:
            cname = contrast["name"]
            rpt.h4(f"Contrast: {cname}")

            de_all = results_dir / method / "deseq2" / cname / "de_all.tsv"
            if de_all.exists():
                rpt.paragraph(f"Full results: `{de_all}`")

                # Show top 10 genes
                rpt.paragraph("**Top 10 genes by adjusted p-value:**")
                lines: List[str] = []
                with open(de_all, encoding="utf-8") as fh:
                    header_line = fh.readline().strip()
                    headers = header_line.split("\t")
                    for i, line in enumerate(fh):
                        if i >= 10:
                            break
                        lines.append(line.strip().split("\t"))

                if lines:
                    rpt.table(headers, lines)
            else:
                rpt.paragraph(f"*(Results not found for {method}/{cname})*")

    # === Trimming comparison ================================================
    rpt.h2("Trimming Method Comparison")
    comparison_md = results_dir / "reports" / "method_comparison.md"
    if comparison_md.exists():
        rpt.paragraph(
            f"See detailed comparison: [{comparison_md.name}]({comparison_md.name})"
        )
        # Inline key sections
        content = comparison_md.read_text(encoding="utf-8")
        # Include everything after the title
        idx = content.find("\n## ")
        if idx >= 0:
            rpt._lines.append(content[idx:])
    else:
        rpt.paragraph("*(Method comparison not available -- run step 10 first)*")

    # === Provenance =========================================================
    rpt.add_config_provenance(cfg)

    # === Write report =======================================================
    report_path = results_dir / "reports" / "report.md"
    rpt.write(report_path)

    # Optional HTML
    html_path = render_html(report_path)
    logger.info("Report written: %s", report_path)
    if html_path != report_path:
        logger.info("HTML report: %s", html_path)

    logger.info("STEP 11 complete.\n")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Step 11: Generate report")
    parser.add_argument("--config", required=True)
    parser.add_argument("--run-id", default=None)
    args = parser.parse_args()

    cfg = load_config(args.config)
    run_id = args.run_id or get_run_id(cfg)
    setup_logging(cfg["project"]["logs_dir"], run_id)
    cfg["_run_id"] = run_id
    cfg["_results_dir"] = str(resolve_results_dir(cfg, run_id))
    cfg["_work_dir"] = str(resolve_work_dir(cfg))

    main(cfg)
