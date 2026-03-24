#!/usr/bin/env python3
"""
Step 09 -- Differential expression analysis.

Backends (set ``deseq2.method`` in config):

  - ``pydeseq2``: Pure Python -- no R needed.  DEFAULT and recommended.
  - ``wrapper``: calls ``DESeq2_wrapper`` (available on course server).
"""

from __future__ import annotations

import csv
import logging
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.metadata import Sample, read_samples_tsv, write_sample_description
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


# =========================================================================
# Backend 1: PyDESeq2 (pure Python -- DEFAULT)
# =========================================================================

def run_pydeseq2(
    count_matrix_path: Path,
    samples: List[Sample],
    contrast_name: str,
    numerator: str,
    denominator: str,
    out_dir: Path,
    cfg: Dict[str, Any],
) -> None:
    """Run differential expression using PyDESeq2 (pure Python).

    Produces the same output files as the R backend:
      - de_all.tsv, de_significant.tsv
      - normalized_counts.tsv, size_factors.tsv
      - pca.pdf, ma_plot.pdf
      - session_info.txt
    """
    import numpy as np
    import pandas as pd

    ensure_dirs(out_dir)
    deseq_cfg = cfg.get("deseq2", {})
    fdr = deseq_cfg.get("fdr_threshold", 0.05)
    lfc_threshold = deseq_cfg.get("lfc_threshold", 0.0)
    ref_level = deseq_cfg.get("reference_level", denominator)

    logger.info("  PyDESeq2: loading count matrix from %s", count_matrix_path)

    # --- Load count matrix (genes x samples) -> transpose to (samples x genes)
    raw = pd.read_csv(count_matrix_path, sep="\t", index_col=0)
    # featureCounts produces genes-as-rows; PyDESeq2 wants samples-as-rows
    counts_df = raw.T.copy()
    # Ensure integer counts
    counts_df = counts_df.round().astype(int)

    logger.info("  Count matrix: %d samples x %d genes", *counts_df.shape)

    # --- Build metadata DataFrame (index = sample names matching count columns)
    meta_records = []
    for s in samples:
        meta_records.append({"sample": s.sample_name, "condition": s.condition})
    metadata = pd.DataFrame(meta_records).set_index("sample")

    # Keep only samples present in both count matrix and metadata
    common = sorted(set(counts_df.index) & set(metadata.index))
    if not common:
        raise ValueError(
            "No matching sample names between count matrix columns and "
            "sample metadata. Check that featureCounts column headers match "
            "sample names in samples.tsv."
        )
    counts_df = counts_df.loc[common]
    metadata = metadata.loc[common]

    # Filter to only conditions in the contrast
    mask = metadata["condition"].isin([numerator, denominator])
    counts_df = counts_df.loc[mask]
    metadata = metadata.loc[mask]

    logger.info("  Samples in contrast (%s vs %s): %d",
                numerator, denominator, len(counts_df))

    if len(counts_df) < 2:
        raise ValueError(
            f"Not enough samples for contrast {numerator} vs {denominator}"
        )

    # --- Import PyDESeq2 ------------------------------------------------
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats

    n_cpus = cfg["project"].get("threads", 4)
    inference = DefaultInference(n_cpus=n_cpus)

    # --- Fit model -------------------------------------------------------
    logger.info("  Fitting DESeq2 model...")
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=metadata,
        design="~condition",
        refit_cooks=True,
        inference=inference,
        ref_level=["condition", ref_level],
    )
    dds.deseq2()

    # --- Extract size factors --------------------------------------------
    size_factors = dds.obs["size_factors"]
    sf_df = pd.DataFrame({
        "sample": size_factors.index,
        "size_factor": size_factors.values,
    })
    sf_df.to_csv(out_dir / "size_factors.tsv", sep="\t", index=False)
    # Also write to parent deseq2 dir for BigWig step
    parent_sf = out_dir.parent / "size_factors.tsv"
    sf_df.to_csv(parent_sf, sep="\t", index=False)

    # --- Normalized counts -----------------------------------------------
    norm_counts = pd.DataFrame(
        dds.layers["normed_counts"],
        index=dds.obs_names,
        columns=dds.var_names,
    ).T  # Back to genes-as-rows for compatibility
    norm_counts.to_csv(out_dir / "normalized_counts.tsv", sep="\t")

    # --- Statistical test ------------------------------------------------
    logger.info("  Running Wald test: %s vs %s", numerator, denominator)
    ds = DeseqStats(
        dds,
        contrast=["condition", numerator, denominator],
        alpha=fdr,
        inference=inference,
    )
    ds.summary()

    results = ds.results_df.copy()
    results = results.sort_values("padj", na_position="last")

    # --- Write all results -----------------------------------------------
    results.to_csv(out_dir / "de_all.tsv", sep="\t")

    # Significant results
    sig_mask = results["padj"].notna() & (results["padj"] < fdr)
    if lfc_threshold > 0:
        sig_mask = sig_mask & (results["log2FoldChange"].abs() > lfc_threshold)
    sig = results.loc[sig_mask]
    sig.to_csv(out_dir / "de_significant.tsv", sep="\t")

    n_up = int((sig["log2FoldChange"] > 0).sum())
    n_down = int((sig["log2FoldChange"] < 0).sum())
    logger.info("  Significant DEGs (FDR < %s): %d (up=%d, down=%d)",
                fdr, len(sig), n_up, n_down)

    # --- Gene lists for functional enrichment --------------------------------
    sig_genes = sig.index.tolist()
    up_genes = sig.loc[sig["log2FoldChange"] > 0].index.tolist()
    down_genes = sig.loc[sig["log2FoldChange"] < 0].index.tolist()

    enrichment_path = out_dir / "top_genes_for_enrichment.txt"
    enrichment_up_path = out_dir / "top_up_genes_for_enrichment.txt"
    enrichment_down_path = out_dir / "top_down_genes_for_enrichment.txt"
    enrichment_top200_path = out_dir / "top200_genes_for_enrichment.txt"
    enrichment_top200_up_path = out_dir / "top200_up_genes_for_enrichment.txt"
    enrichment_top200_down_path = out_dir / "top200_down_genes_for_enrichment.txt"

    enrichment_path.write_text("\n".join(sig_genes), encoding="utf-8")
    enrichment_up_path.write_text("\n".join(up_genes), encoding="utf-8")
    enrichment_down_path.write_text("\n".join(down_genes), encoding="utf-8")
    enrichment_top200_path.write_text("\n".join(sig_genes[:200]), encoding="utf-8")
    enrichment_top200_up_path.write_text("\n".join(up_genes[:200]), encoding="utf-8")
    enrichment_top200_down_path.write_text("\n".join(down_genes[:200]), encoding="utf-8")

    logger.info(
        "  Enrichment gene lists written (all=%d, up=%d, down=%d) -> %s",
        len(sig_genes), len(up_genes), len(down_genes), out_dir,
    )

    # --- Top-20 DE gene table ---------------------------------------------------
    top20 = results.head(20)[["baseMean", "log2FoldChange", "padj"]].copy()
    top20.to_csv(out_dir / "top20_de_genes.tsv", sep="\t")
    logger.info("  Top-20 DE gene table -> %s", out_dir / "top20_de_genes.tsv")

    # --- PCA plot --------------------------------------------------------
    _pydeseq2_pca_plot(dds, metadata, contrast_name, out_dir)

    # --- MA plot ---------------------------------------------------------
    _pydeseq2_ma_plot(results, contrast_name, fdr, out_dir)

    # --- Volcano plot ----------------------------------------------------
    _pydeseq2_volcano_plot(results, contrast_name, fdr, lfc_threshold, out_dir)

    # --- Session info ----------------------------------------------------
    _pydeseq2_session_info(out_dir, contrast_name)

    logger.info("  PyDESeq2 complete -> %s", out_dir)


def _pydeseq2_pca_plot(dds, metadata, contrast_name: str, out_dir: Path) -> None:
    """Generate a PCA plot from the variance-stabilised counts."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
        import pandas as pd
        from sklearn.decomposition import PCA

        # Use log2(normalized + 1) as a simple VST-like transform
        norm = dds.layers["normed_counts"]
        log_norm = np.log2(norm + 1)

        pca = PCA(n_components=2)
        pcs = pca.fit_transform(log_norm)
        pct_var = pca.explained_variance_ratio_ * 100

        pca_df = pd.DataFrame({
            "PC1": pcs[:, 0],
            "PC2": pcs[:, 1],
            "condition": metadata["condition"].values,
            "sample": metadata.index.values,
        })

        fig, ax = plt.subplots(figsize=(8, 6))
        for cond, grp in pca_df.groupby("condition"):
            ax.scatter(grp["PC1"], grp["PC2"], label=cond, s=60)
            for _, row in grp.iterrows():
                ax.annotate(row["sample"], (row["PC1"], row["PC2"]),
                            fontsize=8, ha="left", va="bottom")

        ax.set_xlabel(f"PC1: {pct_var[0]:.1f}% variance")
        ax.set_ylabel(f"PC2: {pct_var[1]:.1f}% variance")
        ax.set_title(f"PCA - {contrast_name}")
        ax.legend()
        fig.tight_layout()
        fig.savefig(out_dir / "pca.pdf")
        plt.close(fig)
        logger.info("  PCA plot saved.")
    except Exception as exc:
        logger.warning("  PCA plot failed: %s", exc)


def _pydeseq2_ma_plot(results, contrast_name: str, fdr: float, out_dir: Path) -> None:
    """Generate an MA plot."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np

        fig, ax = plt.subplots(figsize=(8, 6))

        mask_sig = results["padj"].notna() & (results["padj"] < fdr)
        mask_ns = ~mask_sig

        ax.scatter(
            np.log10(results.loc[mask_ns, "baseMean"] + 1),
            results.loc[mask_ns, "log2FoldChange"],
            s=4, alpha=0.4, color="grey", label="NS",
        )
        ax.scatter(
            np.log10(results.loc[mask_sig, "baseMean"] + 1),
            results.loc[mask_sig, "log2FoldChange"],
            s=6, alpha=0.6, color="red", label=f"FDR < {fdr}",
        )
        ax.axhline(0, color="black", linewidth=0.5)
        ax.set_xlabel("log10(baseMean + 1)")
        ax.set_ylabel("log2 Fold Change")
        ax.set_title(f"MA Plot - {contrast_name}")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(out_dir / "ma_plot.pdf")
        plt.close(fig)
        logger.info("  MA plot saved.")
    except Exception as exc:
        logger.warning("  MA plot failed: %s", exc)


def _pydeseq2_volcano_plot(
    results, contrast_name: str, fdr: float, lfc_threshold: float, out_dir: Path,
) -> None:
    """Generate a volcano plot: log2FC (x) vs -log10(padj) (y)."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np

        fig, ax = plt.subplots(figsize=(8, 6))

        padj = results["padj"].values.astype(float)
        lfc = results["log2FoldChange"].values.astype(float)
        neg_log_p = -np.log10(np.where(padj > 0, padj, 1e-300))

        sig = np.isfinite(padj) & (padj < fdr)
        if lfc_threshold > 0:
            sig = sig & (np.abs(lfc) > lfc_threshold)

        ax.scatter(lfc[~sig], neg_log_p[~sig], s=4, alpha=0.4,
                   color="grey", label="NS")
        ax.scatter(lfc[sig], neg_log_p[sig], s=6, alpha=0.6,
                   color="red", label=f"FDR < {fdr}")
        ax.axhline(-np.log10(fdr), color="blue", linestyle="--", linewidth=0.5,
                   label=f"p_adj = {fdr}")
        if lfc_threshold > 0:
            ax.axvline(lfc_threshold, color="blue", linestyle="--", linewidth=0.5)
            ax.axvline(-lfc_threshold, color="blue", linestyle="--", linewidth=0.5)

        ax.set_xlabel("log2 Fold Change")
        ax.set_ylabel("-log10(adjusted p-value)")
        ax.set_title(f"Volcano Plot - {contrast_name}")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(out_dir / "volcano.pdf")
        plt.close(fig)
        logger.info("  Volcano plot saved.")
    except Exception as exc:
        logger.warning("  Volcano plot failed: %s", exc)


def _pydeseq2_session_info(out_dir: Path, contrast_name: str) -> None:
    """Write a Python equivalent of R's sessionInfo()."""
    import platform
    from datetime import datetime, timezone

    lines = [
        f"PyDESeq2 analysis: {contrast_name}",
        f"Date: {datetime.now(timezone.utc).isoformat()}",
        f"Python: {sys.version}",
        f"Platform: {platform.platform()}",
        "",
    ]

    for pkg in ["pydeseq2", "pandas", "numpy", "scipy", "sklearn",
                "matplotlib", "anndata"]:
        try:
            mod = __import__(pkg)
            lines.append(f"  {pkg}: {getattr(mod, '__version__', 'unknown')}")
        except ImportError:
            lines.append(f"  {pkg}: NOT INSTALLED")

    (out_dir / "session_info.txt").write_text("\n".join(lines), encoding="utf-8")


# =========================================================================
# Backend 2: DESeq2_wrapper (course server)
# =========================================================================

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


# =========================================================================
# Summary helpers
# =========================================================================

def count_degs(results_path: Path, fdr: float = 0.05) -> Dict[str, int]:
    """Count DEGs from a DESeq2/PyDESeq2 results TSV."""
    total = 0
    sig_up = 0
    sig_down = 0

    if not results_path.exists():
        return {"total_tested": 0, "sig_up": 0, "sig_down": 0, "sig_total": 0}

    with open(results_path, encoding="utf-8") as fh:
        header = fh.readline().strip().split("\t")
        # Find padj and log2FoldChange column indices
        padj_idx = -1
        for i, h in enumerate(header):
            if h.strip().lower() == "padj" or "padj" in h.lower():
                padj_idx = i
                break

        lfc_idx = -1
        for i, h in enumerate(header):
            if "log2foldchange" in h.strip().lower():
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
                        sig_up += 1
                else:
                    sig_up += 1

    return {
        "total_tested": total,
        "sig_up": sig_up,
        "sig_down": sig_down,
        "sig_total": sig_up + sig_down,
    }


from src.analysis_unit import build_mapping_units, resolve_count_matrix, resolve_de_dir


# =========================================================================
# Main
# =========================================================================


def _empty_de_stats() -> Dict[str, int]:
    return {"total_tested": 0, "sig_up": 0, "sig_down": 0, "sig_total": 0}


def _read_redundancy_summary(
    results_dir: Path,
) -> Dict[Tuple[str, str, str, str], Dict[str, str]]:
    """Load redundant-branch mapping from redundancy_summary.tsv if present."""
    path = results_dir / "redundancy_summary.tsv"
    out: Dict[Tuple[str, str, str, str], Dict[str, str]] = {}
    if not path.exists():
        return out

    with open(path, encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            is_redundant = str(row.get("is_redundant", "")).strip().lower()
            if is_redundant not in {"1", "true", "yes"}:
                continue
            key = (
                row.get("trim_method", ""),
                row.get("mapper", "star"),
                row.get("mapper_option_set", "default"),
                row.get("count_option_set", "default"),
            )
            out[key] = {
                "canonical_trim_method": row.get("canonical_trim_method", ""),
                "canonical_mapper": row.get("canonical_mapper", ""),
                "canonical_mapper_option_set": row.get("canonical_mapper_option_set", ""),
                "canonical_count_option_set": row.get("canonical_count_option_set", ""),
                "fingerprint": row.get("fingerprint", ""),
                "redundancy_reason": row.get("redundancy_reason", "identical_clean_matrix"),
            }
    if out:
        logger.info("Loaded %d redundant branch mappings from %s", len(out), path)
    return out

def main(cfg: Dict[str, Any], methods_override: List[str] | None = None) -> None:
    """Execute step 09."""
    logger.info("=" * 60)
    logger.info("STEP 09: Differential expression analysis")
    logger.info("=" * 60)

    results_dir = Path(cfg["_results_dir"])
    deseq_cfg = cfg.get("deseq2", {})
    method_backend = deseq_cfg.get("method", "pydeseq2")
    fdr = deseq_cfg.get("fdr_threshold", 0.05)
    comparisons = cfg.get("comparisons", [])
    strict_mode = bool(cfg.get("_strict", False) or deseq_cfg.get("require_complete_summary", False))

    logger.info("Backend: %s", method_backend)

    samples_tsv = Path(cfg.get("_samples_tsv", results_dir / "samples.tsv"))
    samples = read_samples_tsv(samples_tsv)
    methods = methods_override or get_enabled_methods(cfg)
    option_sets = cfg["featurecounts"].get("option_sets", {"default": {}})
    contrast_names = ["auto"] if method_backend == "wrapper" else [c["name"] for c in comparisons]
    if not contrast_names:
        raise ValueError("No DE contrasts configured. Please set config.comparisons.")

    redundant_map = _read_redundancy_summary(results_dir)
    all_de_stats: Dict[Tuple[str, str, str, str, str], Dict[str, Any]] = {}
    integrity_issues: List[str] = []

    mapping_units = build_mapping_units(results_dir, methods)
    for unit in mapping_units:
        method = unit["method"]
        mapper = unit["mapper"]
        mapper_opt = unit["mapper_option_set"]

        for opt_name in option_sets:
            branch_key_4 = (method, mapper, mapper_opt, opt_name)
            redundancy = redundant_map.get(branch_key_4)
            logger.info(
                "--- DE for trim=%s mapper=%s mapper_option=%s count_option=%s ---",
                method, mapper, mapper_opt, opt_name,
            )
            de_dir = resolve_de_dir(results_dir, method, mapper, mapper_opt, opt_name)
            ensure_dirs(de_dir)

            if redundancy is not None:
                logger.info(
                    "Skipping redundant DE branch %s/%s/%s/%s",
                    method, mapper, mapper_opt, opt_name,
                )
                for cname in contrast_names:
                    key = (method, mapper, mapper_opt, opt_name, cname)
                    row = {
                        "trim_method": method,
                        "mapper": mapper,
                        "mapper_option_set": mapper_opt,
                        "count_option_set": opt_name,
                        "contrast": cname,
                        "status": "redundant_skipped",
                        "message": redundancy.get("redundancy_reason", "identical_clean_matrix"),
                        "canonical_trim_method": redundancy.get("canonical_trim_method", ""),
                        "canonical_mapper": redundancy.get("canonical_mapper", ""),
                        "canonical_mapper_option_set": redundancy.get("canonical_mapper_option_set", ""),
                        "canonical_count_option_set": redundancy.get("canonical_count_option_set", ""),
                        "fingerprint": redundancy.get("fingerprint", ""),
                    }
                    row.update(_empty_de_stats())
                    all_de_stats[key] = row
                continue

            count_matrix_path = resolve_count_matrix(
                results_dir, method, mapper, mapper_opt, opt_name
            )
            if count_matrix_path is None:
                logger.error(
                    "Count matrix not found for %s/%s/%s/%s -- skipping",
                    method, mapper, mapper_opt, opt_name,
                )
                for cname in contrast_names:
                    key = (method, mapper, mapper_opt, opt_name, cname)
                    row = {
                        "trim_method": method,
                        "mapper": mapper,
                        "mapper_option_set": mapper_opt,
                        "count_option_set": opt_name,
                        "contrast": cname,
                        "status": "count_matrix_missing",
                        "message": "clean count matrix not found",
                        "canonical_trim_method": "",
                        "canonical_mapper": "",
                        "canonical_mapper_option_set": "",
                        "canonical_count_option_set": "",
                        "fingerprint": "",
                    }
                    row.update(_empty_de_stats())
                    all_de_stats[key] = row
                continue

            sample_desc = de_dir / "sample_description.txt"
            write_sample_description(samples, sample_desc)

            if method_backend == "wrapper":
                logger.info("  Running DESeq2_wrapper...")
                key = (method, mapper, mapper_opt, opt_name, "auto")
                try:
                    run_deseq2_wrapper(count_matrix_path, sample_desc, de_dir, cfg)
                    de_all = de_dir / "DESeq2.de_all.tsv"
                    if de_all.exists():
                        stats = count_degs(de_all, fdr)
                        row = {
                            "trim_method": method,
                            "mapper": mapper,
                            "mapper_option_set": mapper_opt,
                            "count_option_set": opt_name,
                            "contrast": "auto",
                            "status": "completed",
                            "message": "",
                            "canonical_trim_method": "",
                            "canonical_mapper": "",
                            "canonical_mapper_option_set": "",
                            "canonical_count_option_set": "",
                            "fingerprint": "",
                        }
                        row.update(stats)
                        all_de_stats[key] = row
                    else:
                        msg = f"DE output missing after wrapper run: {de_all}"
                        integrity_issues.append(msg)
                        row = {
                            "trim_method": method,
                            "mapper": mapper,
                            "mapper_option_set": mapper_opt,
                            "count_option_set": opt_name,
                            "contrast": "auto",
                            "status": "missing_output",
                            "message": msg,
                            "canonical_trim_method": "",
                            "canonical_mapper": "",
                            "canonical_mapper_option_set": "",
                            "canonical_count_option_set": "",
                            "fingerprint": "",
                        }
                        row.update(_empty_de_stats())
                        all_de_stats[key] = row
                except Exception as exc:
                    msg = f"DE wrapper failed for {method}/{mapper}/{mapper_opt}/{opt_name}/auto: {exc}"
                    integrity_issues.append(msg)
                    row = {
                        "trim_method": method,
                        "mapper": mapper,
                        "mapper_option_set": mapper_opt,
                        "count_option_set": opt_name,
                        "contrast": "auto",
                        "status": "de_failed",
                        "message": str(exc),
                        "canonical_trim_method": "",
                        "canonical_mapper": "",
                        "canonical_mapper_option_set": "",
                        "canonical_count_option_set": "",
                        "fingerprint": "",
                    }
                    row.update(_empty_de_stats())
                    all_de_stats[key] = row

            elif method_backend == "pydeseq2":
                for contrast in comparisons:
                    cname = contrast["name"]
                    num = contrast["numerator"]
                    den = contrast["denominator"]
                    contrast_dir = de_dir / cname
                    logger.info("  Contrast: %s (%s vs %s)", cname, num, den)
                    key = (method, mapper, mapper_opt, opt_name, cname)
                    try:
                        run_pydeseq2(
                            count_matrix_path, samples,
                            cname, num, den, contrast_dir, cfg,
                        )
                        de_all = contrast_dir / "de_all.tsv"
                        if de_all.exists():
                            stats = count_degs(de_all, fdr)
                            row = {
                                "trim_method": method,
                                "mapper": mapper,
                                "mapper_option_set": mapper_opt,
                                "count_option_set": opt_name,
                                "contrast": cname,
                                "status": "completed",
                                "message": "",
                                "canonical_trim_method": "",
                                "canonical_mapper": "",
                                "canonical_mapper_option_set": "",
                                "canonical_count_option_set": "",
                                "fingerprint": "",
                            }
                            row.update(stats)
                            all_de_stats[key] = row
                        else:
                            msg = f"DE output missing after pydeseq2 run: {de_all}"
                            integrity_issues.append(msg)
                            row = {
                                "trim_method": method,
                                "mapper": mapper,
                                "mapper_option_set": mapper_opt,
                                "count_option_set": opt_name,
                                "contrast": cname,
                                "status": "missing_output",
                                "message": msg,
                                "canonical_trim_method": "",
                                "canonical_mapper": "",
                                "canonical_mapper_option_set": "",
                                "canonical_count_option_set": "",
                                "fingerprint": "",
                            }
                            row.update(_empty_de_stats())
                            all_de_stats[key] = row
                    except Exception as exc:
                        msg = f"DE failed for {method}/{mapper}/{mapper_opt}/{opt_name}/{cname}: {exc}"
                        integrity_issues.append(msg)
                        row = {
                            "trim_method": method,
                            "mapper": mapper,
                            "mapper_option_set": mapper_opt,
                            "count_option_set": opt_name,
                            "contrast": cname,
                            "status": "de_failed",
                            "message": str(exc),
                            "canonical_trim_method": "",
                            "canonical_mapper": "",
                            "canonical_mapper_option_set": "",
                            "canonical_count_option_set": "",
                            "fingerprint": "",
                        }
                        row.update(_empty_de_stats())
                        all_de_stats[key] = row
            else:
                raise ValueError(
                    f"Unknown deseq2.method: '{method_backend}'. "
                    "Use 'pydeseq2' or 'wrapper'."
                )

    # Ensure every expected branch/contrast is represented exactly once.
    expected_keys: List[Tuple[str, str, str, str, str]] = []
    for unit in mapping_units:
        for opt_name in option_sets:
            for cname in contrast_names:
                expected_keys.append(
                    (unit["method"], unit["mapper"], unit["mapper_option_set"], opt_name, cname)
                )

    for key in expected_keys:
        if key in all_de_stats:
            continue
        method, mapper, mapper_opt, opt_name, cname = key
        msg = "DE summary slot missing from run accounting"
        integrity_issues.append(f"{msg}: {method}/{mapper}/{mapper_opt}/{opt_name}/{cname}")
        row = {
            "trim_method": method,
            "mapper": mapper,
            "mapper_option_set": mapper_opt,
            "count_option_set": opt_name,
            "contrast": cname,
            "status": "malformed",
            "message": msg,
            "canonical_trim_method": "",
            "canonical_mapper": "",
            "canonical_mapper_option_set": "",
            "canonical_count_option_set": "",
            "fingerprint": "",
        }
        row.update(_empty_de_stats())
        all_de_stats[key] = row

    de_summary_path = results_dir / "de_summary.tsv"
    fields = [
        "trim_method",
        "mapper",
        "mapper_option_set",
        "count_option_set",
        "contrast",
        "status",
        "message",
        "total_tested",
        "sig_up",
        "sig_down",
        "sig_total",
        "canonical_trim_method",
        "canonical_mapper",
        "canonical_mapper_option_set",
        "canonical_count_option_set",
        "fingerprint",
    ]
    with open(de_summary_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for key in sorted(all_de_stats.keys()):
            writer.writerow(all_de_stats[key])
    logger.info("DE summary -> %s", de_summary_path)

    status_counts: Dict[str, int] = {}
    for row in all_de_stats.values():
        status = str(row.get("status", "unknown"))
        status_counts[status] = status_counts.get(status, 0) + 1
    logger.info(
        "DE summary statuses: %s",
        ", ".join(f"{k}={v}" for k, v in sorted(status_counts.items())),
    )

    if integrity_issues:
        for msg in integrity_issues:
            logger.warning("DE summary integrity: %s", msg)
        if strict_mode:
            raise RuntimeError(
                "DE summary integrity checks failed in strict mode. "
                "Inspect de_summary.tsv status/message columns."
            )

    logger.info("STEP 09 complete.\n")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Step 09: Differential expression")
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
