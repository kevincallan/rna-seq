"""Tests for primary trim method control and optional comparison behavior."""

from __future__ import annotations

from pathlib import Path
import types
import importlib

import pytest

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.metadata import Sample
from src.utils import get_effective_trim_methods, get_trim_config_summary


def _base_cfg() -> dict:
    return {
        "project": {"name": "x", "results_dir": "results", "work_dir": "work", "logs_dir": "logs"},
        "trimming": {
            "primary_method": "cutadapt",
            "compare_methods": False,
            "comparison_methods": [],
            "none": {"enabled": True},
            "cutadapt": {"enabled": True, "adapter_fwd": "", "adapter_rev": "", "quality": 20, "min_length": 25, "extra_args": ""},
        },
    }


def test_primary_only_resolution():
    cfg = _base_cfg()
    assert get_effective_trim_methods(cfg) == ["cutadapt"]
    summary = get_trim_config_summary(cfg)
    assert summary["primary_method"] == "cutadapt"
    assert summary["compare_methods"] is False
    assert summary["effective_methods"] == ["cutadapt"]


def test_compare_enabled_resolution_and_dedupe_order():
    cfg = _base_cfg()
    cfg["trimming"]["compare_methods"] = True
    cfg["trimming"]["comparison_methods"] = ["none", "cutadapt", "none"]
    assert get_effective_trim_methods(cfg) == ["cutadapt", "none"]


def test_invalid_primary_method_fails_fast():
    cfg = _base_cfg()
    cfg["trimming"]["primary_method"] = "missing"
    with pytest.raises(ValueError, match="primary_method"):
        get_effective_trim_methods(cfg)


def test_invalid_comparison_method_fails_fast():
    cfg = _base_cfg()
    cfg["trimming"]["compare_methods"] = True
    cfg["trimming"]["comparison_methods"] = ["none", "bogus"]
    with pytest.raises(ValueError, match="comparison method"):
        get_effective_trim_methods(cfg)


def test_methods_override_precedence_and_validation():
    cfg = _base_cfg()
    assert get_effective_trim_methods(cfg, ["none", "cutadapt", "none"]) == ["none", "cutadapt"]
    with pytest.raises(ValueError, match="unknown/unconfigured"):
        get_effective_trim_methods(cfg, ["none", "bad_method"])


def test_step00_creates_only_effective_trim_dirs(monkeypatch, tmp_path):
    step00 = importlib.import_module("scripts.00_validate_env")

    cfg = _base_cfg()
    cfg["project"]["results_dir"] = str(tmp_path / "results")
    cfg["project"]["work_dir"] = str(tmp_path / "work")
    cfg["project"]["cache_dir"] = str(tmp_path / "cache")
    cfg["project"]["logs_dir"] = str(tmp_path / "logs")
    cfg["data"] = {"metadata_csv": str(tmp_path / "m.csv"), "fastq_dir": str(tmp_path)}
    cfg["references"] = {"genome_index": str(tmp_path), "gtf": str(tmp_path / "a.gtf"), "genome_fasta": str(tmp_path / "a.fa")}
    cfg["tools"] = {}
    cfg["mapping"] = {"backends": {"star": {"enabled": True}}}
    cfg["featurecounts"] = {"backend": "featurecounts"}
    cfg["deseq2"] = {}
    cfg["comparisons"] = []

    monkeypatch.setattr(step00, "validate_data", lambda _cfg: [])
    monkeypatch.setattr(step00, "validate_references", lambda _cfg: [])
    monkeypatch.setattr(step00, "validate_tools", lambda _cfg: {"dummy": "1"})
    monkeypatch.setattr(step00, "write_run_manifest", lambda *args, **kwargs: tmp_path / "manifest.json")

    run_id = "RID"
    out = step00.main(cfg, run_id)
    work_dir = Path(out["_work_dir"])
    assert (work_dir / "trimmed" / "cutadapt").exists()
    assert not (work_dir / "trimmed" / "none").exists()


def test_cutadapt_blank_adapters_omit_flags(monkeypatch, tmp_path):
    step02 = importlib.import_module("scripts.02_trim_reads")

    sample = Sample("SRR1", "ctrl", 1, layout="paired")
    work_dir = tmp_path / "work"
    link_dir = work_dir / "fastq_links"
    out_dir = tmp_path / "results"
    link_dir.mkdir(parents=True)
    (link_dir / f"{sample.sample_name}_1.fastq.gz").write_text("", encoding="utf-8")
    (link_dir / f"{sample.sample_name}_2.fastq.gz").write_text("", encoding="utf-8")

    cfg = _base_cfg()
    cfg["tools"] = {"cutadapt": "cutadapt"}
    cfg["trimming"]["cutadapt"]["adapter_fwd"] = ""
    cfg["trimming"]["cutadapt"]["adapter_rev"] = ""

    captured = {}

    def fake_run(cmd, **kwargs):
        captured["cmd"] = list(cmd)
        return types.SimpleNamespace(stdout="", stderr="")

    monkeypatch.setattr(step02, "run_cmd", fake_run)
    metrics = step02._trim_cutadapt(sample, work_dir, cfg, out_dir)
    assert "-a" not in captured["cmd"]
    assert "-A" not in captured["cmd"]
    assert metrics["trim_mode"] == "quality-only"


def test_cutadapt_with_adapters_includes_flags(monkeypatch, tmp_path):
    step02 = importlib.import_module("scripts.02_trim_reads")

    sample = Sample("SRR1", "ctrl", 1, layout="paired")
    work_dir = tmp_path / "work"
    link_dir = work_dir / "fastq_links"
    out_dir = tmp_path / "results"
    link_dir.mkdir(parents=True)
    (link_dir / f"{sample.sample_name}_1.fastq.gz").write_text("", encoding="utf-8")
    (link_dir / f"{sample.sample_name}_2.fastq.gz").write_text("", encoding="utf-8")

    cfg = _base_cfg()
    cfg["tools"] = {"cutadapt": "cutadapt"}
    cfg["trimming"]["cutadapt"]["adapter_fwd"] = "AAA"
    cfg["trimming"]["cutadapt"]["adapter_rev"] = "TTT"

    captured = {}

    def fake_run(cmd, **kwargs):
        captured["cmd"] = list(cmd)
        return types.SimpleNamespace(stdout="", stderr="")

    monkeypatch.setattr(step02, "run_cmd", fake_run)
    metrics = step02._trim_cutadapt(sample, work_dir, cfg, out_dir)
    assert "-a" in captured["cmd"]
    assert "-A" in captured["cmd"]
    assert metrics["trim_mode"] == "adapter+quality"


def test_fastqc_raw_always_runs_even_primary_only(monkeypatch, tmp_path):
    step03 = importlib.import_module("scripts.03_qc_fastqc")

    results_dir = tmp_path / "results"
    work_dir = tmp_path / "work" / "RID"
    link_dir = work_dir / "fastq_links"
    link_dir.mkdir(parents=True)

    sample = Sample("SRR1", "ctrl", 1, layout="paired")
    samples_tsv = results_dir / "samples.tsv"
    samples_tsv.parent.mkdir(parents=True)
    samples_tsv.write_text(
        "run_id\tcondition\treplicate\tsample_name\tlayout\tr1\tr2\n"
        f"{sample.run_id}\t{sample.condition}\t{sample.replicate}\t{sample.sample_name}\tpaired\t\t\n",
        encoding="utf-8",
    )
    (link_dir / f"{sample.sample_name}_1.fastq.gz").write_text("", encoding="utf-8")
    (link_dir / f"{sample.sample_name}_2.fastq.gz").write_text("", encoding="utf-8")

    cfg = _base_cfg()
    cfg["project"]["work_dir"] = str(tmp_path / "work")
    cfg["project"]["threads"] = 1
    cfg["tools"] = {"fastqc": "fastqc"}
    cfg["_results_dir"] = str(results_dir)
    cfg["_work_dir"] = str(work_dir)
    cfg["_samples_tsv"] = str(samples_tsv)
    cfg["_run_id"] = "RID"

    calls = []

    def fake_run_fastqc(fastq_files, out_dir, threads, exe="fastqc"):
        calls.append((list(fastq_files), Path(out_dir)))

    monkeypatch.setattr(step03, "run_fastqc", fake_run_fastqc)
    step03.main(cfg, methods_override=None)

    assert calls, "Expected raw FastQC call"
    # First call is raw by implementation order.
    assert "raw_qc" in str(calls[0][1])
