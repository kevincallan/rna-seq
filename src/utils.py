"""
Shared utilities for the RNA-seq pipeline.

Provides logging setup, subprocess execution, hashing, path helpers,
tool availability checks, and YAML config loading/validation.
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Union

import yaml

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

def setup_logging(
    log_dir: Union[str, Path],
    run_id: str,
    level: int = logging.INFO,
) -> logging.Logger:
    """Configure root logger with file + console handlers.

    Parameters
    ----------
    log_dir : path-like
        Directory to write the log file into.
    run_id : str
        Used in the log filename: ``<run_id>.log``.
    level : int
        Logging level (default INFO).

    Returns
    -------
    logging.Logger
        The configured root logger.
    """
    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"{run_id}.log"

    fmt = logging.Formatter(
        "%(asctime)s [%(levelname)-7s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    root = logging.getLogger()
    root.setLevel(level)
    # Clear pre-existing handlers to avoid duplicates on re-import
    root.handlers.clear()

    fh = logging.FileHandler(log_file, mode="a", encoding="utf-8")
    fh.setLevel(level)
    fh.setFormatter(fmt)
    root.addHandler(fh)

    ch = logging.StreamHandler(sys.stderr)
    ch.setLevel(level)
    ch.setFormatter(fmt)
    root.addHandler(ch)

    root.info("Logging initialised -> %s", log_file)
    return root


# ---------------------------------------------------------------------------
# Subprocess execution
# ---------------------------------------------------------------------------

def run_cmd(
    cmd: Sequence[str],
    *,
    description: str = "",
    cwd: Optional[Union[str, Path]] = None,
    env: Optional[Dict[str, str]] = None,
    capture: bool = True,
    check: bool = True,
) -> subprocess.CompletedProcess:
    """Run an external command safely (no shell=True).

    Parameters
    ----------
    cmd : sequence of str
        Command and arguments, e.g. ``["samtools", "index", "file.bam"]``.
    description : str
        Human-readable label for log messages.
    cwd : path-like, optional
        Working directory for the command.
    env : dict, optional
        Additional environment variables (merged with ``os.environ``).
    capture : bool
        If True, capture stdout and stderr.
    check : bool
        If True, raise ``subprocess.CalledProcessError`` on non-zero exit.

    Returns
    -------
    subprocess.CompletedProcess
    """
    logger = logging.getLogger("run_cmd")
    cmd_str = " ".join(str(c) for c in cmd)
    label = f"[{description}] " if description else ""
    logger.info("%sRunning: %s", label, cmd_str)

    merged_env = None
    if env:
        merged_env = {**os.environ, **env}

    result = subprocess.run(
        [str(c) for c in cmd],
        cwd=cwd,
        env=merged_env,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.PIPE if capture else None,
        encoding="utf-8",
        check=False,  # we handle check manually for better logging
    )

    if result.stdout:
        for line in result.stdout.strip().splitlines():
            logger.debug("%sSTDOUT: %s", label, line)
    if result.stderr:
        for line in result.stderr.strip().splitlines():
            logger.debug("%sSTDERR: %s", label, line)

    if check and result.returncode != 0:
        logger.error(
            "%sCommand failed (exit %d): %s\nSTDERR:\n%s",
            label,
            result.returncode,
            cmd_str,
            result.stderr or "(empty)",
        )
        raise subprocess.CalledProcessError(
            result.returncode, cmd_str, result.stdout, result.stderr
        )

    logger.info("%sCompleted (exit %d)", label, result.returncode)
    return result


def run_cmd_shell(
    cmd: str,
    *,
    description: str = "",
    cwd: Optional[Union[str, Path]] = None,
    check: bool = True,
) -> subprocess.CompletedProcess:
    """Run a shell command string (use sparingly -- pipes, redirects).

    Only use when piping between commands is unavoidable.
    """
    logger = logging.getLogger("run_cmd_shell")
    label = f"[{description}] " if description else ""
    logger.info("%sRunning (shell): %s", label, cmd)

    result = subprocess.run(
        cmd,
        shell=True,
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        encoding="utf-8",
        check=False,
    )

    if check and result.returncode != 0:
        logger.error(
            "%sShell command failed (exit %d): %s\nSTDERR:\n%s",
            label,
            result.returncode,
            cmd,
            result.stderr or "(empty)",
        )
        raise subprocess.CalledProcessError(
            result.returncode, cmd, result.stdout, result.stderr
        )

    logger.info("%sCompleted (exit %d)", label, result.returncode)
    return result


# ---------------------------------------------------------------------------
# Hashing helpers
# ---------------------------------------------------------------------------

def hash_file(path: Union[str, Path], algorithm: str = "sha256") -> str:
    """Return hex digest of a file."""
    h = hashlib.new(algorithm)
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def hash_string(text: str, algorithm: str = "sha256") -> str:
    """Return hex digest of a string."""
    return hashlib.new(algorithm, text.encode("utf-8")).hexdigest()


def hash_config(cfg: Dict[str, Any]) -> str:
    """Deterministic SHA-256 of a config dictionary."""
    canon = json.dumps(cfg, sort_keys=True, default=str)
    return hash_string(canon)


# ---------------------------------------------------------------------------
# Tool availability
# ---------------------------------------------------------------------------

def check_tool(name: str, path: str = "") -> str:
    """Verify that an external tool is available.

    Parameters
    ----------
    name : str
        Human-readable tool name (for error messages).
    path : str
        Executable name or absolute path.  Falls back to *name* if empty.

    Returns
    -------
    str
        The resolved path (from ``shutil.which``).

    Raises
    ------
    FileNotFoundError
        If the tool cannot be found on ``$PATH``.
    """
    exe = path or name
    resolved = shutil.which(exe)
    if resolved is None:
        raise FileNotFoundError(
            f"Tool '{name}' not found (tried '{exe}'). "
            "Install it or set the correct path in config.yaml."
        )
    return resolved


def get_tool_version(exe: str) -> str:
    """Try to retrieve a version string from a tool."""
    for flag in ("--version", "-version", "-v"):
        try:
            r = subprocess.run(
                [exe, flag],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                encoding="utf-8",
                check=False,
            )
            out = (r.stdout + r.stderr).strip()
            if out:
                return out.splitlines()[0]
        except Exception:
            continue
    return "unknown"


# ---------------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------------

def ensure_dirs(*paths: Union[str, Path]) -> None:
    """Create directories if they do not exist."""
    for p in paths:
        Path(p).mkdir(parents=True, exist_ok=True)


def safe_name(text: str) -> str:
    """Sanitise a string for use in file/directory names."""
    import re
    text = text.strip().lower()
    text = re.sub(r"[/\\:\s]+", "_", text)
    text = re.sub(r"[^a-z0-9_\-.]", "", text)
    text = re.sub(r"_+", "_", text)
    return text.strip("_")


# ---------------------------------------------------------------------------
# Config loading & validation
# ---------------------------------------------------------------------------

_REQUIRED_SECTIONS = [
    "project",
    "data",
    "column_mapping",
    "references",
    "tools",
    "trimming",
    "featurecounts",
    "deseq2",
    "comparisons",
]


def load_config(path: Union[str, Path]) -> Dict[str, Any]:
    """Load and validate a YAML config file.

    Returns
    -------
    dict
        Parsed config with ``_config_path`` injected.

    Raises
    ------
    FileNotFoundError
        If the YAML file does not exist.
    ValueError
        If required sections are missing.
    """
    path = Path(path).resolve()
    if not path.is_file():
        raise FileNotFoundError(f"Config file not found: {path}")

    with open(path, "r", encoding="utf-8") as fh:
        cfg: Dict[str, Any] = yaml.safe_load(fh)

    missing = [s for s in _REQUIRED_SECTIONS if s not in cfg]
    if missing:
        raise ValueError(
            f"Config is missing required sections: {', '.join(missing)}"
        )

    # Inject metadata
    cfg["_config_path"] = str(path)
    cfg["_config_hash"] = hash_config(cfg)

    return cfg


def get_enabled_methods(cfg: Dict[str, Any]) -> List[str]:
    """Return list of trimming method names that are enabled in config."""
    methods = []
    control_keys = {"primary_method", "compare_methods", "comparison_methods"}
    for name, params in cfg.get("trimming", {}).items():
        if name in control_keys:
            continue
        if isinstance(params, dict) and params.get("enabled", True):
            methods.append(name)
        elif params is None:
            # 'none' method may have no params dict
            methods.append(name)
    return sorted(methods)


def _configured_trim_methods(cfg: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    trimming = cfg.get("trimming", {})
    control_keys = {"primary_method", "compare_methods", "comparison_methods"}
    out: Dict[str, Dict[str, Any]] = {}
    for name, params in trimming.items():
        if name in control_keys:
            continue
        if params is None:
            out[name] = {}
        elif isinstance(params, dict):
            out[name] = dict(params)
    return out


def _dedupe_preserve_order(items: List[str]) -> List[str]:
    seen: Set[str] = set()
    out: List[str] = []
    for x in items:
        if x not in seen:
            out.append(x)
            seen.add(x)
    return out


def get_trim_config_summary(cfg: Dict[str, Any]) -> Dict[str, Any]:
    """Validate and summarize trim-method config for the current run."""
    trimming = cfg.get("trimming", {})
    configured = _configured_trim_methods(cfg)

    primary = str(trimming.get("primary_method", "")).strip()
    if not primary:
        raise ValueError(
            "Invalid trimming config: trimming.primary_method must be set."
        )
    if primary not in configured:
        raise ValueError(
            f"Invalid trimming config: primary_method '{primary}' is not configured under trimming."
        )

    compare = bool(trimming.get("compare_methods", False))
    raw_cmp = trimming.get("comparison_methods", [])
    if raw_cmp is None:
        raw_cmp = []
    if not isinstance(raw_cmp, list):
        raise ValueError(
            "Invalid trimming config: trimming.comparison_methods must be a list."
        )
    cmp_methods = [str(m).strip() for m in raw_cmp if str(m).strip()]
    for m in cmp_methods:
        if m not in configured:
            raise ValueError(
                f"Invalid trimming config: comparison method '{m}' is not configured under trimming."
            )

    effective = [primary]
    if compare:
        effective.extend(cmp_methods)
    effective = _dedupe_preserve_order(effective)
    if not effective:
        raise ValueError(
            "Invalid trimming config: effective trim method list is empty."
        )

    non_selected = [m for m in configured.keys() if m not in effective]
    return {
        "primary_method": primary,
        "compare_methods": compare,
        "comparison_methods": cmp_methods,
        "effective_methods": effective,
        "configured_methods": list(configured.keys()),
        "non_selected_methods": non_selected,
    }


def get_effective_trim_methods(
    cfg: Dict[str, Any],
    methods_override: Optional[List[str]] = None,
) -> List[str]:
    """Return effective trim methods honoring config model and CLI override."""
    configured = _configured_trim_methods(cfg)
    if methods_override:
        cleaned = [str(m).strip() for m in methods_override if str(m).strip()]
        cleaned = _dedupe_preserve_order(cleaned)
        unknown = [m for m in cleaned if m not in configured]
        if unknown:
            raise ValueError(
                "Invalid --methods override: unknown/unconfigured method(s): "
                + ", ".join(unknown)
            )
        if not cleaned:
            raise ValueError("Invalid --methods override: no methods provided.")
        return cleaned
    return get_trim_config_summary(cfg)["effective_methods"]


def get_run_id(cfg: Dict[str, Any]) -> str:
    """Generate a deterministic run identifier from project name + timestamp."""
    ts = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    project = safe_name(cfg["project"]["name"])
    return f"{project}_{ts}"


def resolve_results_dir(cfg: Dict[str, Any], run_id: str) -> Path:
    """Return the results directory for a given run."""
    return Path(cfg["project"]["results_dir"]) / run_id


def resolve_work_dir(cfg: Dict[str, Any], run_id: Optional[str] = None) -> Path:
    """Return the run-scoped work directory.

    The base work root remains ``project.work_dir``. Per-run intermediates are
    isolated under ``<work_root>/<run_id>`` to avoid cross-run interference.
    """
    base = Path(cfg["project"]["work_dir"])
    rid = run_id or str(cfg.get("_run_id", "")).strip()
    return base / rid if rid else base


def resolve_cache_dir(cfg: Dict[str, Any]) -> Path:
    """Return shared cache directory for reusable artifacts."""
    project = cfg.get("project", {})
    return Path(project.get("cache_dir", "cache"))
