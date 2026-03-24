"""Tests for src.analysis_unit helpers."""

from __future__ import annotations

import csv
import tempfile
from pathlib import Path

import pytest

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.analysis_unit import (
    AnalysisUnit,
    build_full_analysis_units,
    build_mapping_units,
    infer_analysis_units_from_de_summary,
    read_selected_analysis,
    resolve_count_matrix,
    resolve_de_dir,
    resolve_fc_dir,
    unit_label,
    write_selected_analysis,
)


class TestAnalysisUnit:
    def test_label(self):
        au = AnalysisUnit("cutadapt", "star", "default", "strict")
        assert au.label == "cutadapt/star/default/strict"

    def test_mapping_key(self):
        au = AnalysisUnit("cutadapt", "star", "default", "strict")
        assert au.mapping_key == ("cutadapt", "star", "default")

    def test_ordering(self):
        a = AnalysisUnit("cutadapt", "star", "default", "default")
        b = AnalysisUnit("none", "star", "default", "default")
        assert a < b

    def test_as_dict(self):
        au = AnalysisUnit("none", "star", "default", "strict")
        d = au.as_dict()
        assert d == {
            "method": "none",
            "mapper": "star",
            "mapper_option_set": "default",
            "count_option_set": "strict",
        }


class TestResolveCountMatrix:
    def test_primary_path_exists(self, tmp_path):
        fc_dir = tmp_path / "cutadapt" / "featurecounts" / "star" / "default"
        fc_dir.mkdir(parents=True)
        matrix = fc_dir / "clean_matrix_strict.tsv"
        matrix.write_text("Geneid\tS1\n")

        result = resolve_count_matrix(tmp_path, "cutadapt", "star", "default", "strict")
        assert result == matrix

    def test_legacy_fallback(self, tmp_path):
        legacy_dir = tmp_path / "cutadapt" / "featurecounts"
        legacy_dir.mkdir(parents=True)
        legacy = legacy_dir / "clean_matrix_default.tsv"
        legacy.write_text("Geneid\tS1\n")

        result = resolve_count_matrix(tmp_path, "cutadapt", "star", "default", "default")
        assert result == legacy

    def test_returns_none_when_missing(self, tmp_path):
        result = resolve_count_matrix(tmp_path, "cutadapt", "star", "default", "missing")
        assert result is None


class TestResolveDeDirs:
    def test_resolve_de_dir(self, tmp_path):
        result = resolve_de_dir(tmp_path, "cutadapt", "star", "default", "strict")
        expected = tmp_path / "cutadapt" / "deseq2" / "star" / "default" / "strict"
        assert result == expected

    def test_resolve_fc_dir(self, tmp_path):
        result = resolve_fc_dir(tmp_path, "cutadapt", "star", "default")
        expected = tmp_path / "cutadapt" / "featurecounts" / "star" / "default"
        assert result == expected


class TestBuildMappingUnits:
    def test_from_mapping_summary(self, tmp_path):
        summary = tmp_path / "mapping_summary.tsv"
        with open(summary, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["trim_method", "mapper", "mapper_option_set", "sample"],
                               delimiter="\t")
            w.writeheader()
            w.writerow({"trim_method": "cutadapt", "mapper": "star",
                        "mapper_option_set": "default", "sample": "S1"})

        units = build_mapping_units(tmp_path, ["cutadapt"])
        assert len(units) == 1
        assert units[0]["method"] == "cutadapt"

    def test_fallback_when_no_summary(self, tmp_path):
        units = build_mapping_units(tmp_path, ["none"])
        assert len(units) == 1
        assert units[0] == {"method": "none", "mapper": "star", "mapper_option_set": "default"}


class TestBuildFullAnalysisUnits:
    def test_cross_product(self, tmp_path):
        summary = tmp_path / "mapping_summary.tsv"
        with open(summary, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["trim_method", "mapper", "mapper_option_set", "sample"],
                               delimiter="\t")
            w.writeheader()
            w.writerow({"trim_method": "cutadapt", "mapper": "star",
                        "mapper_option_set": "default", "sample": "S1"})

        option_sets = {"default": {}, "strict": {}}
        units = build_full_analysis_units(tmp_path, ["cutadapt"], option_sets)
        assert len(units) == 2
        assert units[0].count_option_set == "default"
        assert units[1].count_option_set == "strict"


class TestInferFromDeSummary:
    def test_reads_count_option_set(self, tmp_path):
        de_summary = tmp_path / "de_summary.tsv"
        with open(de_summary, "w", newline="") as fh:
            w = csv.DictWriter(
                fh,
                fieldnames=["trim_method", "mapper", "mapper_option_set",
                            "count_option_set", "contrast"],
                delimiter="\t",
            )
            w.writeheader()
            w.writerow({"trim_method": "cutadapt", "mapper": "star",
                        "mapper_option_set": "default", "count_option_set": "strict",
                        "contrast": "a_vs_b"})
            w.writerow({"trim_method": "cutadapt", "mapper": "star",
                        "mapper_option_set": "default", "count_option_set": "default",
                        "contrast": "a_vs_b"})

        units = infer_analysis_units_from_de_summary(tmp_path, ["cutadapt"])
        assert len(units) == 2
        opts = {u.count_option_set for u in units}
        assert opts == {"default", "strict"}


class TestSelectedAnalysis:
    def test_write_and_read(self, tmp_path):
        au = AnalysisUnit("cutadapt", "star", "default", "strict")
        write_selected_analysis(tmp_path, au, "test_reason")

        result = read_selected_analysis(tmp_path)
        assert result == au

    def test_read_missing(self, tmp_path):
        result = read_selected_analysis(tmp_path)
        assert result is None


class TestUnitLabel:
    def test_with_count_option(self):
        d = {"method": "cutadapt", "mapper": "star",
             "mapper_option_set": "default", "count_option_set": "strict"}
        assert unit_label(d) == "cutadapt/star/default/strict"

    def test_without_count_option(self):
        d = {"method": "cutadapt", "mapper": "star", "mapper_option_set": "default"}
        assert unit_label(d) == "cutadapt/star/default"
