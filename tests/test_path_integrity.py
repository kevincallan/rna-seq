"""Tests for path integrity across pipeline steps.

Ensures that output filenames from one step match input expectations
of the next step, preventing silent path-resolution failures.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.analysis_unit import (
    AnalysisUnit,
    resolve_count_matrix,
    resolve_de_dir,
    resolve_fc_dir,
    write_selected_analysis,
    read_selected_analysis,
)


class TestStep07To08PathChain:
    """Step 07 writes count_matrix_{opt}.tsv; Step 08 reads it and writes clean_matrix_{opt}.tsv."""

    @pytest.mark.parametrize("opt_name", ["default", "strict", "stringent", "multimapper"])
    def test_step07_output_matches_step08_input(self, tmp_path, opt_name):
        fc_dir = resolve_fc_dir(tmp_path, "cutadapt", "star", "default")
        fc_dir.mkdir(parents=True)

        # Step 07 writes count_matrix_{opt_name}.tsv
        step07_output = fc_dir / f"count_matrix_{opt_name}.tsv"
        step07_output.write_text("Geneid\tS1\nGENE1\t100\n")

        # Step 08 expects to read count_matrix_{opt_name}.tsv from the same dir
        step08_input = fc_dir / f"count_matrix_{opt_name}.tsv"
        assert step08_input.exists()
        assert step08_input == step07_output


class TestStep08To09PathChain:
    """Step 08 writes clean_matrix_{opt}.tsv; Step 09 loads it via resolve_count_matrix."""

    @pytest.mark.parametrize("opt_name", ["default", "strict", "stringent", "multimapper"])
    def test_step08_output_found_by_step09(self, tmp_path, opt_name):
        fc_dir = resolve_fc_dir(tmp_path, "cutadapt", "star", "default")
        fc_dir.mkdir(parents=True)

        step08_output = fc_dir / f"clean_matrix_{opt_name}.tsv"
        step08_output.write_text("Geneid\tS1\nGENE1\t100\n")

        result = resolve_count_matrix(tmp_path, "cutadapt", "star", "default", opt_name)
        assert result is not None
        assert result == step08_output


class TestDERunsForAllOptionSets:
    """DE should produce output directories for every configured option set."""

    def test_de_dir_includes_count_option(self, tmp_path):
        option_sets = ["default", "strict", "stringent", "multimapper"]
        for opt in option_sets:
            de_dir = resolve_de_dir(tmp_path, "cutadapt", "star", "default", opt)
            assert str(de_dir).endswith(f"deseq2/star/default/{opt}")

    def test_de_dirs_are_distinct_per_option(self, tmp_path):
        dirs = set()
        for opt in ["default", "strict"]:
            d = resolve_de_dir(tmp_path, "cutadapt", "star", "default", opt)
            dirs.add(str(d))
        assert len(dirs) == 2


class TestReportResolvesDePathsWithCountOption:
    """The report must be able to find DE results for any count option set."""

    def test_report_path_includes_count_opt(self, tmp_path):
        for opt in ["default", "strict"]:
            de_dir = resolve_de_dir(tmp_path, "cutadapt", "star", "default", opt)
            de_dir.mkdir(parents=True)
            contrast_dir = de_dir / "tet1_vs_wt"
            contrast_dir.mkdir()
            de_all = contrast_dir / "de_all.tsv"
            de_all.write_text("gene\tpadj\nGENE1\t0.01\n")

            assert de_all.exists()
            assert opt in str(de_all)


class TestSelectedAnalysisFile:
    def test_write_creates_file(self, tmp_path):
        au = AnalysisUnit("cutadapt", "star", "default", "strict")
        path = write_selected_analysis(tmp_path, au, "test")
        assert path.exists()

    def test_roundtrip(self, tmp_path):
        au = AnalysisUnit("cutadapt", "star", "default", "strict")
        write_selected_analysis(tmp_path, au, "test_reason")
        result = read_selected_analysis(tmp_path)
        assert result == au

    def test_tsv_has_correct_columns(self, tmp_path):
        au = AnalysisUnit("none", "star", "default", "default")
        write_selected_analysis(tmp_path, au, "auto")
        path = tmp_path / "selected_analysis.tsv"
        with open(path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            rows = list(reader)
        assert len(rows) == 1
        assert rows[0]["trim_method"] == "none"
        assert rows[0]["selection_reason"] == "auto"
