#!/usr/bin/env python3
"""
Step 00 -- Validate environment.

Checks config validity, verifies tools are installed, confirms reference
files exist, creates output directories, and writes a run manifest.
"""

from __future__ import annotations

import json
import logging
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List

# Allow running as ``./py scripts/00_validate_env.py``
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.utils import (
    check_tool,
    ensure_dirs,
    get_enabled_methods,
    get_run_id,
    get_tool_version,
    hash_config,
    load_config,
    resolve_results_dir,
    resolve_work_dir,
    setup_logging,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Validation helpers
# ---------------------------------------------------------------------------

def validate_references(cfg: Dict[str, Any]) -> List[str]:
    """Check that reference files/directories exist.  Return list of errors."""
    errors: List[str] = []
    refs = cfg.get("references", {})

    genome_idx = refs.get("genome_index", "")
    if genome_idx and not Path(genome_idx).exists():
        errors.append(f"STAR genome index not found: {genome_idx}")

    gtf = refs.get("gtf", "")
    if gtf and not Path(gtf).exists():
        errors.append(f"GTF annotation not found: {gtf}")

    fasta = refs.get("genome_fasta", "")
    if fasta and not Path(fasta).exists():
        errors.append(f"Genome FASTA not found: {fasta}")

    return errors


def validate_data(cfg: Dict[str, Any]) -> List[str]:
    """Check that input data paths exist."""
    errors: List[str] = []
    data = cfg.get("data", {})

    csv_path = data.get("metadata_csv", "")
    if csv_path and not Path(csv_path).exists():
        errors.append(f"Metadata CSV not found: {csv_path}")

    fq_dir = data.get("fastq_dir", "")
    if fq_dir and not Path(fq_dir).is_dir():
        errors.append(f"FASTQ directory not found: {fq_dir}")

    return errors


def validate_tools(cfg: Dict[str, Any]) -> Dict[str, str]:
    """Check tool availability and collect versions.

    Returns
    -------
    dict
        Mapping of tool name -> version string.  Tools not found raise errors.
    """
    tools_cfg = cfg.get("tools", {})
    methods = get_enabled_methods(cfg)
    versions: Dict[str, str] = {}

    # Always-required compiled tools -- just the mapper
    required = ["star"]

    # These are Python packages but validate they're on PATH if used as CLI
    # (multiqc, cutadapt, bamcoverage are pip-installable)
    required.extend(["fastqc", "multiqc", "bamcoverage"])

    # samtools only needed if pysam is not available
    try:
        import pysam
        versions["pysam"] = pysam.__version__
        logger.info("  %-16s OK  (v%s -- replaces samtools)", "pysam", pysam.__version__)
    except ImportError:
        required.append("samtools")

    # featureCounts only needed if not using htseq backend
    counting_backend = cfg.get("featurecounts", {}).get("backend", "featurecounts")
    if counting_backend == "htseq":
        try:
            import HTSeq
            versions["HTSeq"] = HTSeq.__version__
            logger.info("  %-16s OK  (v%s -- replaces featureCounts)", "HTSeq", HTSeq.__version__)
        except ImportError:
            logger.error("  HTSeq not installed. Run: pip install HTSeq")
            raise FileNotFoundError(
                "Python package 'HTSeq' not found. "
                "Install with: pip install HTSeq"
            )
    else:
        required.append("featurecounts")

    # Add trimming tools for enabled methods
    if "cutadapt" in methods:
        required.append("cutadapt")
    if "fastp" in methods:
        required.append("fastp")
    if "trimmomatic" in methods:
        required.append("trimmomatic")

    # R is only needed for R-based DESeq2 backends
    de_method = cfg.get("deseq2", {}).get("method", "pydeseq2")
    if de_method == "rscript":
        required.append("rscript")
    elif de_method == "pydeseq2":
        # Verify pydeseq2 is importable (pure Python, no system tool)
        try:
            import pydeseq2
            versions["pydeseq2"] = getattr(pydeseq2, "__version__", "installed")
            logger.info("  %-16s OK  (v%s)", "pydeseq2", versions["pydeseq2"])
        except ImportError:
            logger.error("  pydeseq2 not installed. Run: pip install pydeseq2")
            raise FileNotFoundError(
                "Python package 'pydeseq2' not found. "
                "Install with: pip install pydeseq2"
            )

    for tool_key in required:
        exe = tools_cfg.get(tool_key, tool_key)
        try:
            resolved = check_tool(tool_key, exe)
            versions[tool_key] = get_tool_version(resolved)
            logger.info("  %-16s OK  (%s)", tool_key, versions[tool_key][:60])
        except FileNotFoundError as exc:
            logger.error("  %-16s MISSING  %s", tool_key, exc)
            raise

    return versions


# ---------------------------------------------------------------------------
# Run manifest
# ---------------------------------------------------------------------------

def write_run_manifest(
    cfg: Dict[str, Any],
    run_id: str,
    tool_versions: Dict[str, str],
    results_dir: Path,
) -> Path:
    """Write a JSON manifest capturing run provenance."""
    manifest = {
        "run_id": run_id,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "config_hash": hash_config(cfg),
        "config_path": cfg.get("_config_path", ""),
        "python_version": sys.version,
        "tool_versions": tool_versions,
        "enabled_trim_methods": get_enabled_methods(cfg),
    }

    manifest_path = results_dir / "run_manifest.json"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with open(manifest_path, "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2, default=str)

    logger.info("Run manifest -> %s", manifest_path)
    return manifest_path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(cfg: Dict[str, Any], run_id: str) -> Dict[str, Any]:
    """Execute step 00.

    Parameters
    ----------
    cfg : dict
        Loaded config.
    run_id : str
        Unique run identifier.

    Returns
    -------
    dict
        Enriched config with ``_run_id``, ``_results_dir``, ``_work_dir``,
        ``_tool_versions``.
    """
    logger.info("=" * 60)
    logger.info("STEP 00: Validate environment")
    logger.info("=" * 60)

    # --- Validate data & references -----------------------------------------
    errors: List[str] = []
    errors.extend(validate_data(cfg))

    ref_errors = validate_references(cfg)
    if ref_errors:
        for e in ref_errors:
            logger.warning("Reference warning: %s", e)
        # Warnings only -- references may be on a different machine

    if errors:
        for e in errors:
            logger.error("Validation error: %s", e)
        raise RuntimeError(
            f"Environment validation failed with {len(errors)} error(s). "
            "See log for details."
        )

    # --- Check tools --------------------------------------------------------
    logger.info("Checking tool availability...")
    tool_versions = validate_tools(cfg)

    # --- Create directories -------------------------------------------------
    results_dir = resolve_results_dir(cfg, run_id)
    work_dir = resolve_work_dir(cfg)
    logs_dir = Path(cfg["project"]["logs_dir"])

    ensure_dirs(
        results_dir,
        work_dir,
        work_dir / "fastq_links",
        work_dir / "trimmed",
        logs_dir,
    )

    for method in get_enabled_methods(cfg):
        ensure_dirs(
            work_dir / "trimmed" / method,
            results_dir / method / "star",
            results_dir / method / "featurecounts",
            results_dir / method / "deseq2",
            results_dir / method / "bigwig",
            results_dir / method / "qc",
        )

    ensure_dirs(results_dir / "reports")
    logger.info("Directory structure created under %s", results_dir)

    # --- Write manifest -----------------------------------------------------
    write_run_manifest(cfg, run_id, tool_versions, results_dir)

    # --- Enrich config with runtime metadata --------------------------------
    cfg["_run_id"] = run_id
    cfg["_results_dir"] = str(results_dir)
    cfg["_work_dir"] = str(work_dir)
    cfg["_tool_versions"] = tool_versions

    logger.info("STEP 00 complete.\n")
    return cfg


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Step 00: Validate environment")
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    parser.add_argument("--run-id", default=None, help="Override run ID")
    args = parser.parse_args()

    cfg = load_config(args.config)
    run_id = args.run_id or get_run_id(cfg)
    setup_logging(cfg["project"]["logs_dir"], run_id)
    main(cfg, run_id)
