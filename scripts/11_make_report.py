#!/usr/bin/env python3
"""
Step 11 -- Generate final analysis report.

Produces ``results/<run_id>/reports/report.md`` (and optionally HTML)
combining all analysis summaries into a single document suitable for
submission or reference.
"""

from __future__ import annotations

import csv
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


def _infer_methods_from_outputs(results_dir: Path, cfg: Dict[str, Any]) -> List[str]:
    """Infer methods that actually have outputs for this run."""
    methods: set[str] = set()

    # Prefer methods observed in summaries
    for summary in ["mapping_summary.tsv", "featurecounts_summary.tsv", "de_summary.tsv"]:
        p = results_dir / summary
        if not p.exists():
            continue
        with open(p, encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                m = row.get("method") or row.get("trim_method")
                if m:
                    methods.add(m)

    # Fallback: detect method dirs in results tree
    if not methods:
        for p in results_dir.iterdir():
            if not p.is_dir():
                continue
            if p.name in {"reports"}:
                continue
            if (p / "star").exists() or (p / "deseq2").exists() or (p / "qc").exists():
                methods.add(p.name)

    # Final fallback: config-enabled methods
    if not methods:
        methods.update(get_enabled_methods(cfg))

    return sorted(methods)


def _infer_analysis_units(results_dir: Path, methods: List[str]) -> List[Dict[str, str]]:
    """Infer trim/mapping analysis units from DE or mapping summaries."""
    units: set[tuple[str, str, str]] = set()
    de_summary = results_dir / "de_summary.tsv"
    if de_summary.exists():
        with open(de_summary, encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                method = row.get("method") or row.get("trim_method") or ""
                if method not in methods:
                    continue
                mapper = row.get("mapper", "star")
                mapper_opt = row.get("mapper_option_set", "default")
                units.add((method, mapper, mapper_opt))

    if not units:
        mapping_summary = results_dir / "mapping_summary.tsv"
        if mapping_summary.exists():
            with open(mapping_summary, encoding="utf-8") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    method = row.get("method") or row.get("trim_method") or ""
                    if method not in methods:
                        continue
                    mapper = row.get("mapper", "star")
                    mapper_opt = row.get("mapper_option_set", "default")
                    units.add((method, mapper, mapper_opt))

    if not units:
        for method in methods:
            units.add((method, "star", "default"))

    return [
        {"method": m, "mapper": p, "mapper_option_set": o}
        for (m, p, o) in sorted(units)
    ]


def _de_base_for_unit(results_dir: Path, unit: Dict[str, str]) -> Path:
    method = unit["method"]
    mapper = unit.get("mapper", "star")
    mapper_opt = unit.get("mapper_option_set", "default")
    nested = results_dir / method / "deseq2" / mapper / mapper_opt
    if nested.exists():
        return nested
    return results_dir / method / "deseq2"


def main(cfg: Dict[str, Any], methods_override: List[str] | None = None) -> None:
    """Execute step 11."""
    logger.info("=" * 60)
    logger.info("STEP 11: Generate report")
    logger.info("=" * 60)

    results_dir = Path(cfg["_results_dir"])
    methods = methods_override or _infer_methods_from_outputs(results_dir, cfg)
    analysis_units = _infer_analysis_units(results_dir, methods)
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
    if methods == ["none"]:
        rpt.paragraph(
            "Only the `none` baseline was run: reads were passed through "
            "without adapter/quality trimming."
        )
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

    # Per-analysis-unit detail
    for unit in analysis_units:
        method = unit["method"]
        mapper = unit.get("mapper", "star")
        mapper_opt = unit.get("mapper_option_set", "default")
        rpt.h3(f"Analysis Unit: {method} / {mapper} / {mapper_opt}")

        de_base = _de_base_for_unit(results_dir, unit)
        for contrast in comparisons:
            cname = contrast["name"]
            rpt.h4(f"Contrast: {cname}")

            de_all = de_base / cname / "de_all.tsv"
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
                rpt.paragraph(
                    f"*(Results not found for {method}/{mapper}/{mapper_opt}/{cname})*"
                )

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

    # === IGV Loading Guide ===================================================
    rpt.h2("IGV Coverage Tracks")
    rpt.paragraph(
        "BigWig files for viewing read coverage in IGV are located in the "
        "results tree under `<trim_method>/bigwig/<mapper>/<mapper_option>/`. "
        "Load these into IGV alongside the reference genome to inspect "
        "coverage at top DE gene loci."
    )
    rpt.paragraph("**How to use:**")
    rpt.bullet("Open IGV and load the correct genome (e.g. hg38 or mm39).")
    rpt.bullet("File > Load from File > select `.bw` files for each sample.")
    rpt.bullet(
        "Navigate to a top DE gene (from `top20_de_genes.tsv`) to see "
        "coverage differences between conditions."
    )
    rpt.bullet("Take screenshots of 1-2 interesting genes for the report.")
    rpt.paragraph("")

    # List BigWig files found in this run
    bw_files = sorted(results_dir.rglob("*.bw"))
    if bw_files:
        rpt.paragraph(f"**BigWig files produced ({len(bw_files)}):**")
        for bw in bw_files[:20]:
            rpt.bullet(f"`{bw.relative_to(results_dir)}`")
        if len(bw_files) > 20:
            rpt.paragraph(f"*(... and {len(bw_files) - 20} more)*")
    else:
        rpt.paragraph("*(No BigWig files found -- run step 06)*")

    # === Top DE genes for IGV ================================================
    rpt.h3("Suggested Genes to Inspect in IGV")
    rpt.paragraph(
        "The top differentially expressed genes from each contrast are listed "
        "in `top20_de_genes.tsv` in each DE output directory. Navigate to "
        "these genes in IGV to visualise the coverage difference."
    )

    # === Functional Enrichment Checklist =====================================
    rpt.h2("Functional Enrichment Analysis Checklist")
    rpt.paragraph(
        "The pipeline exports gene lists ready for external enrichment tools. "
        "Complete this checklist to earn advanced-task marks."
    )
    rpt.paragraph("**Enrichment gene list files (per contrast):**")
    rpt.bullet("`top_genes_for_enrichment.txt` -- all significant DEGs")
    rpt.bullet("`top_up_genes_for_enrichment.txt` -- upregulated only")
    rpt.bullet("`top_down_genes_for_enrichment.txt` -- downregulated only")
    rpt.bullet("`top200_*.txt` -- capped versions for tools with upload limits")
    rpt.paragraph("")

    rpt.paragraph("**Step-by-step workflow:**")
    rpt.paragraph(
        "1. Locate the enrichment files in "
        "`results/<run_id>/<trim_method>/deseq2/<mapper>/<mapper_option>/<contrast>/`"
    )
    rpt.paragraph(
        "2. **g:Profiler** (https://biit.cs.ut.ee/gprofiler/gost): "
        "Paste gene list, select organism, run. Screenshot the top enriched "
        "GO terms and pathways."
    )
    rpt.paragraph(
        "3. **Enrichr** (https://maayanlab.cloud/Enrichr/): "
        "Paste gene list. Check GO Biological Process, KEGG, and "
        "WikiPathways tabs. Screenshot combined-score bar charts."
    )
    rpt.paragraph(
        "4. **iDEP** (https://bioinformatics.sdstate.edu/idep/): "
        "Upload `normalized_counts.tsv` from the same DE directory. "
        "Provides PCA, clustering, enrichment, and pathway analysis."
    )
    rpt.paragraph(
        "5. Include 1-2 key enrichment findings in the report with "
        "adjusted p-values and GO/pathway IDs."
    )

    # List enrichment files found
    enrich_files = sorted(results_dir.rglob("*_for_enrichment.txt"))
    if enrich_files:
        rpt.paragraph(f"\n**Enrichment files produced ({len(enrich_files)}):**")
        for ef in enrich_files[:30]:
            rpt.bullet(f"`{ef.relative_to(results_dir)}`")
        if len(enrich_files) > 30:
            rpt.paragraph(f"*(... and {len(enrich_files) - 30} more)*")

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
    parser.add_argument("--methods", nargs="*", default=None)
    args = parser.parse_args()

    cfg = load_config(args.config)
    run_id = args.run_id or get_run_id(cfg)
    setup_logging(cfg["project"]["logs_dir"], run_id)
    cfg["_run_id"] = run_id
    cfg["_results_dir"] = str(resolve_results_dir(cfg, run_id))
    cfg["_work_dir"] = str(resolve_work_dir(cfg))

    main(cfg, methods_override=args.methods)
