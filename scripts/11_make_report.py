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
    get_effective_trim_methods,
    get_trim_config_summary,
    get_run_id,
    load_config,
    resolve_results_dir,
    resolve_work_dir,
    setup_logging,
)
from src.analysis_unit import (
    AnalysisUnit,
    infer_analysis_units_from_de_summary,
    read_selected_analysis,
    read_selected_count_comparison,
    read_selected_visualisation,
    resolve_de_dir,
    resolve_de_base_legacy,
)

logger = logging.getLogger(__name__)


def _infer_methods_from_outputs(results_dir: Path, cfg: Dict[str, Any]) -> List[str]:
    """Infer methods that actually have outputs for this run."""
    methods: set[str] = set()

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

    if not methods:
        for p in results_dir.iterdir():
            if not p.is_dir():
                continue
            if p.name in {"reports"}:
                continue
            if (p / "star").exists() or (p / "deseq2").exists() or (p / "qc").exists():
                methods.add(p.name)

    if not methods:
        methods.update(get_effective_trim_methods(cfg))

    return sorted(methods)


def _resolve_de_all_for_unit(results_dir: Path, au: AnalysisUnit, contrast: str) -> Path:
    """Find de_all.tsv for an analysis unit, checking new then legacy paths."""
    de_dir = resolve_de_dir(
        results_dir, au.method, au.mapper,
        au.mapper_option_set, au.count_option_set,
    )
    de_path = de_dir / contrast / "de_all.tsv"
    if de_path.exists():
        return de_path

    legacy_base = resolve_de_base_legacy(
        results_dir, au.method, au.mapper, au.mapper_option_set,
    )
    legacy_path = legacy_base / contrast / "de_all.tsv"
    if legacy_path.exists():
        return legacy_path

    return de_path  # return expected (new) path even if missing


def _branch4_from_row(row: Dict[str, str], count_key: str = "count_option_set") -> tuple[str, str, str, str]:
    return (
        row.get("trim_method", row.get("method", "")),
        row.get("mapper", "star"),
        row.get("mapper_option_set", "default"),
        row.get(count_key, "default"),
    )


def _read_tsv(path: Path) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    if not path.exists():
        return rows
    with open(path, encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append(dict(row))
    return rows


def main(cfg: Dict[str, Any], methods_override: List[str] | None = None) -> None:
    """Execute step 11."""
    logger.info("=" * 60)
    logger.info("STEP 11: Generate report")
    logger.info("=" * 60)

    results_dir = Path(cfg["_results_dir"])
    methods = get_effective_trim_methods(cfg, methods_override)
    analysis_units = infer_analysis_units_from_de_summary(results_dir, methods)
    comparisons = cfg.get("comparisons", [])
    project_name = cfg["project"]["name"]
    trim_cfg = get_trim_config_summary(cfg)

    selected_count = read_selected_count_comparison(results_dir)
    selected_vis = read_selected_visualisation(results_dir)
    selected_legacy = read_selected_analysis(results_dir)
    selected_primary = selected_legacy
    selected = selected_primary or selected_count

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
    rpt.paragraph(f"**Primary trim method:** `{trim_cfg['primary_method']}`")
    rpt.paragraph(f"**Trim-method comparison enabled:** `{trim_cfg['compare_methods']}`")
    rpt.paragraph(f"**Effective trim methods this run:** `{', '.join(trim_cfg['effective_methods'])}`")
    if trim_cfg["non_selected_methods"]:
        rpt.paragraph(
            "**Intentionally not run:** "
            f"`{', '.join(trim_cfg['non_selected_methods'])}` "
            "(trim comparison disabled or excluded)."
        )

    cutadapt_cfg = cfg.get("trimming", {}).get("cutadapt", {})
    adapter_fwd = str(cutadapt_cfg.get("adapter_fwd", "")).strip()
    adapter_rev = str(cutadapt_cfg.get("adapter_rev", "")).strip()
    quality = int(cutadapt_cfg.get("quality", 20))
    min_length = int(cutadapt_cfg.get("min_length", 25))
    has_adapter = bool(adapter_fwd or adapter_rev)
    quality_active = quality > 0 or min_length > 0
    if has_adapter and quality_active:
        cutadapt_mode = "adapter+quality"
    elif has_adapter:
        cutadapt_mode = "adapter-only"
    else:
        cutadapt_mode = "quality-only"
    rpt.paragraph(f"**Cutadapt mode:** `{cutadapt_mode}`")

    rpt.paragraph(f"Methods run: {', '.join(methods)}")
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
        with open(fc_summary, encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            modes = sorted({row.get("counting_mode", "").strip() for row in reader if row.get("counting_mode", "").strip()})
        if modes:
            rpt.paragraph(f"**Counting mode:** {', '.join(modes)}")
        else:
            layout = str(cfg.get("data", {}).get("layout", "paired")).lower()
            inferred_mode = "read_pairs" if layout == "paired" else "reads"
            rpt.paragraph(f"**Counting mode (inferred from layout):** {inferred_mode}")
    else:
        rpt.paragraph("*(featureCounts summary not available)*")

    # === Filtering summary ==================================================
    rpt.h2("Count Matrix Filtering")
    filter_summary = results_dir / "filtering_summary.tsv"
    if filter_summary.exists():
        rpt.table_from_tsv(filter_summary)
    else:
        rpt.paragraph("*(Filtering summary not available)*")

    # === Selected analysis branches =========================================
    if selected_primary or selected_count or selected_vis:
        rpt.h2("Selected Analysis Branches")
        if selected_primary:
            rpt.paragraph(
                f"**Primary analysis branch:** `{selected_primary.label}`"
            )
            sel_tsv = results_dir / "selected_analysis.tsv"
            if sel_tsv.exists():
                rpt.table_from_tsv(sel_tsv)
        if selected_count:
            rpt.paragraph(
                f"**Count-comparison branch:** `{selected_count.label}`"
            )
            cc_tsv = results_dir / "selected_count_comparison.tsv"
            if cc_tsv.exists():
                rpt.table_from_tsv(cc_tsv)
        if selected_vis:
            rpt.paragraph(
                f"**Visualisation branch (BigWig/IGV):** `{selected_vis.label}`"
            )
            vis_tsv = results_dir / "selected_visualisation.tsv"
            if vis_tsv.exists():
                rpt.table_from_tsv(vis_tsv)
        else:
            rpt.paragraph("*(selected_visualisation.tsv not found)*")

    # === DE results per analysis unit / contrast ============================
    rpt.h2("Differential Expression Results")
    de_summary = results_dir / "de_summary.tsv"
    if de_summary.exists():
        rpt.h3("DE Summary Table")
        rpt.table_from_tsv(de_summary)
    else:
        rpt.paragraph("*(DE summary not available)*")

    for au in analysis_units:
        markers: List[str] = []
        if selected_primary is not None and au == selected_primary:
            markers.append("PRIMARY")
        if selected_count is not None and au == selected_count:
            markers.append("COUNT_COMPARISON")
        if selected_vis is not None and au == selected_vis:
            markers.append("VISUALISATION")
        marker = f" [{' & '.join(markers)}]" if markers else ""
        rpt.h3(f"Analysis Unit: {au.label}{marker}")

        for contrast in comparisons:
            cname = contrast["name"]
            rpt.h4(f"Contrast: {cname}")

            de_all = _resolve_de_all_for_unit(results_dir, au, cname)

            if not de_all.exists():
                rpt.paragraph(
                    f"DE analysis did not complete for "
                    f"`{au.method}/{au.mapper}/{au.mapper_option_set}"
                    f"/{au.count_option_set}/{cname}`."
                )
                rpt.paragraph(f"Expected path: `{de_all}`")
                continue

            rpt.paragraph(f"Full results: `{de_all}`")

            lines: List[List[str]] = []
            headers: List[str] = []
            sig_count = 0
            with open(de_all, encoding="utf-8") as fh:
                header_line = fh.readline().strip()
                headers = header_line.split("\t")
                padj_idx = -1
                for idx, h in enumerate(headers):
                    if "padj" in h.lower():
                        padj_idx = idx
                        break

                for i, line in enumerate(fh):
                    parts = line.strip().split("\t")
                    if padj_idx >= 0:
                        try:
                            pval = float(parts[padj_idx])
                            if pval < 0.05:
                                sig_count += 1
                        except (ValueError, IndexError):
                            pass
                    if i < 10:
                        lines.append(parts)

            if sig_count == 0:
                rpt.paragraph(
                    f"DE completed; **0 genes** reached significance "
                    f"at FDR < 0.05."
                )
            else:
                rpt.paragraph(
                    f"**{sig_count} significant genes** (FDR < 0.05)."
                )

            if lines:
                rpt.paragraph("**Top 10 genes by adjusted p-value:**")
                rpt.table(headers, lines)

    # === Trimming comparison ================================================
    rpt.h2("Method and Parameter Comparison")
    comparison_md = results_dir / "reports" / "method_comparison.md"
    if comparison_md.exists():
        rpt.paragraph(
            f"See detailed comparison: [{comparison_md.name}]({comparison_md.name})"
        )
        content = comparison_md.read_text(encoding="utf-8")
        idx = content.find("\n## ")
        if idx >= 0:
            rpt._lines.append(content[idx:])
    else:
        rpt.paragraph("*(Method comparison not available -- run step 10 first)*")

    # === Pipeline Health ====================================================
    rpt.h2("Pipeline Health")
    health_checks: List[Dict[str, str]] = []

    for name, path in [
        ("mapping_summary.tsv", results_dir / "mapping_summary.tsv"),
        ("featurecounts_summary.tsv", results_dir / "featurecounts_summary.tsv"),
        ("filtering_summary.tsv", results_dir / "filtering_summary.tsv"),
        ("de_summary.tsv", results_dir / "de_summary.tsv"),
        ("redundancy_summary.tsv", results_dir / "redundancy_summary.tsv"),
        ("selected_count_comparison.tsv", results_dir / "selected_count_comparison.tsv"),
        ("selected_visualisation.tsv", results_dir / "selected_visualisation.tsv"),
        ("selected_analysis.tsv", results_dir / "selected_analysis.tsv"),
    ]:
        status = "PASS" if path.exists() else "MISSING"
        health_checks.append({"check": name, "status": status, "path": str(path)})

    fc_rows = _read_tsv(results_dir / "featurecounts_summary.tsv")
    filt_rows = _read_tsv(results_dir / "filtering_summary.tsv")
    de_rows = _read_tsv(results_dir / "de_summary.tsv")
    branch4_fc = {
        (
            r.get("trim_method", r.get("method", "")),
            r.get("mapper", "star"),
            r.get("mapper_option_set", "default"),
            r.get("option_set", "default"),
        )
        for r in fc_rows
    }
    branch4_filt = {
        (
            r.get("trim_method", r.get("method", "")),
            r.get("mapper", "star"),
            r.get("mapper_option_set", "default"),
            r.get("option_set", "default"),
        )
        for r in filt_rows
    }
    branch4_de_status: Dict[tuple[str, str, str, str], List[str]] = {}
    for row in de_rows:
        key = _branch4_from_row(row)
        branch4_de_status.setdefault(key, []).append(str(row.get("status", "completed")))

    for key in sorted(branch4_fc | branch4_filt):
        statuses = branch4_de_status.get(key, [])
        key_label = "/".join(key)
        if not statuses:
            health_checks.append({
                "check": f"Branch summary sync: {key_label}",
                "status": "MISSING",
                "path": "de_summary.tsv",
            })
            continue
        unique = {s.lower() for s in statuses}
        if "malformed" in unique:
            status = "MALFORMED"
        elif unique <= {"redundant_skipped"}:
            status = "REDUNDANT_SKIPPED"
        elif unique <= {"completed"}:
            sig_sum = 0
            for row in de_rows:
                if _branch4_from_row(row) != key:
                    continue
                try:
                    sig_sum += int(float(row.get("sig_total", "0") or 0))
                except ValueError:
                    pass
            status = "COMPLETED_ZERO_DEG" if sig_sum == 0 else "PASS"
        elif "completed" in unique:
            status = "PASS"
        elif "de_failed" in unique or "missing_output" in unique:
            status = "MISSING"
        else:
            status = "MALFORMED"
        health_checks.append({
            "check": f"Branch summary sync: {key_label}",
            "status": status,
            "path": "de_summary.tsv",
        })

    selected_refs = [
        ("selected_primary_analysis", selected_primary),
        ("selected_count_comparison", selected_count),
        ("selected_visualisation", selected_vis),
    ]
    real_branches = set(branch4_de_status.keys())
    for label, selected_unit in selected_refs:
        if selected_unit is None:
            continue
        key = (
            selected_unit.method,
            selected_unit.mapper,
            selected_unit.mapper_option_set,
            selected_unit.count_option_set,
        )
        health_checks.append({
            "check": f"{label} exists in de_summary",
            "status": "PASS" if key in real_branches else "MISSING",
            "path": "de_summary.tsv",
        })

    bw_cfg = cfg.get("bigwig", {})
    if bw_cfg.get("enabled", False) and bw_cfg.get("mode", "all_units") == "selected_only" and selected_vis:
        bw_dir = results_dir / selected_vis.method / "bigwig" / selected_vis.mapper / selected_vis.mapper_option_set
        has_bw = bw_dir.exists() and any(bw_dir.glob("*.bw"))
        health_checks.append({
            "check": "Selected visualisation BigWigs present",
            "status": "PASS" if has_bw else "MISSING",
            "path": str(bw_dir),
        })

    for au in analysis_units:
        for contrast in comparisons:
            cname = contrast["name"]
            de_all = _resolve_de_all_for_unit(results_dir, au, cname)
            if not de_all.exists():
                health_checks.append({
                    "check": f"DE: {au.label}/{cname}",
                    "status": "MISSING",
                    "path": str(de_all),
                })
            else:
                health_checks.append({
                    "check": f"DE: {au.label}/{cname}",
                    "status": "PASS",
                    "path": str(de_all),
                })

    rpt.table(
        ["Check", "Status", "Path"],
        [[c["check"], c["status"], f"`{c['path']}`"] for c in health_checks],
    )

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

    bw_files = sorted(results_dir.rglob("*.bw"))
    if bw_files:
        rpt.paragraph(f"**BigWig files produced ({len(bw_files)}):**")
        for bw in bw_files[:20]:
            rpt.bullet(f"`{bw.relative_to(results_dir)}`")
        if len(bw_files) > 20:
            rpt.paragraph(f"*(... and {len(bw_files) - 20} more)*")
    else:
        rpt.paragraph("*(No BigWig files found -- run BigWig step if needed)*")

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
        "The pipeline exports gene lists ready for external enrichment tools."
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
        "`results/<run_id>/<trim_method>/deseq2/<mapper>/<mapper_option>/<count_option>/<contrast>/`"
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
