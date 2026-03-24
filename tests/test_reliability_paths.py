"""Reliability regressions: run-scoped work, cache, and BAM diagnostics."""

from __future__ import annotations

import importlib
import sys
import types
from pathlib import Path

import pytest

from src.utils import resolve_cache_dir, resolve_work_dir


def test_resolve_work_dir_is_run_scoped():
    cfg = {"project": {"work_dir": "work"}, "_run_id": "run_abc"}
    assert resolve_work_dir(cfg) == Path("work") / "run_abc"
    assert resolve_work_dir(cfg, "override") == Path("work") / "override"


def test_resolve_cache_dir_default_and_override():
    cfg_default = {"project": {}}
    assert resolve_cache_dir(cfg_default) == Path("cache")
    cfg_custom = {"project": {"cache_dir": "/tmp/my_cache"}}
    assert resolve_cache_dir(cfg_custom) == Path("/tmp/my_cache")


def test_prepare_gtf_uses_cache_and_reuses_existing(tmp_path):
    mod = importlib.import_module("scripts.07_featurecounts")
    raw_gtf = tmp_path / "raw.gtf"
    raw_gtf.write_text(
        "#header\nchr1\tsrc\texon\t1\t10\t.\t+\t.\tgene_id \"g1\"; gene_name \"A\";\n",
        encoding="utf-8",
    )
    cache_dir = tmp_path / "cache"
    cfg = {
        "references": {"gtf": str(raw_gtf)},
        "featurecounts": {"gtf_filter": "gene_name"},
        "_cache_dir": str(cache_dir),
    }

    out1 = Path(mod._prepare_gtf(cfg, tmp_path))
    assert out1.exists()
    assert cache_dir / "filtered_gtf" in out1.parents
    content_1 = out1.read_text(encoding="utf-8")

    # Second call should reuse exactly the same cached path/content.
    out2 = Path(mod._prepare_gtf(cfg, tmp_path))
    assert out2 == out1
    assert out2.read_text(encoding="utf-8") == content_1
    # Lock file should not be left behind after successful build.
    assert not out2.with_suffix(out2.suffix + ".lock").exists()


def test_prepare_gtf_reuses_existing_when_lock_present(tmp_path):
    mod = importlib.import_module("scripts.07_featurecounts")
    raw_gtf = tmp_path / "raw.gtf"
    raw_gtf.write_text("chr1\tsrc\texon\t1\t10\t.\t+\t.\tgene_name \"A\";\n", encoding="utf-8")
    cfg = {
        "references": {"gtf": str(raw_gtf)},
        "featurecounts": {"gtf_filter": "gene_name"},
        "_cache_dir": str(tmp_path / "cache"),
    }
    out = Path(mod._prepare_gtf(cfg, tmp_path))
    lock = out.with_suffix(out.suffix + ".lock")
    lock.write_text("other-run", encoding="utf-8")
    reused = Path(mod._prepare_gtf(cfg, tmp_path))
    assert reused == out


def test_validate_bam_integrity_reports_corruption(monkeypatch, tmp_path):
    mod = importlib.import_module("scripts.05_map_star")
    bam = tmp_path / "bad.bam"
    bam.write_text("not-a-bam", encoding="utf-8")

    monkeypatch.setattr(mod.shutil, "which", lambda *_: None)

    class BadAlignmentFile:
        def __init__(self, *args, **kwargs):
            raise OSError("BGZF read error")

    fake_pysam = types.SimpleNamespace(AlignmentFile=BadAlignmentFile)
    monkeypatch.setitem(sys.modules, "pysam", fake_pysam)

    with pytest.raises(RuntimeError, match="appears incomplete/corrupt"):
        mod.validate_bam_integrity(
            bam,
            {"tools": {}, "project": {}},
            run_id="rid",
            sample_name="s1",
            method="none",
            mapper="star",
            mapper_opt="default",
        )


def test_index_bam_wraps_pysam_index_failures(monkeypatch, tmp_path):
    mod = importlib.import_module("scripts.05_map_star")
    bam = tmp_path / "ok.bam"
    bam.write_bytes(b"bam")

    def boom(*args, **kwargs):
        raise RuntimeError("index failed")

    fake_pysam = types.SimpleNamespace(index=boom)
    monkeypatch.setitem(sys.modules, "pysam", fake_pysam)

    with pytest.raises(RuntimeError, match="BAM indexing failed after integrity validation"):
        mod.index_bam(
            bam,
            {"tools": {}, "project": {}},
            run_id="rid",
            sample_name="s1",
            method="none",
            mapper="star",
            mapper_opt="default",
        )
