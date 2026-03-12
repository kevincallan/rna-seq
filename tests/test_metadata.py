"""
Tests for metadata parsing, condition mapping, and replicate numbering.

Run with:  pytest tests/ -v
"""

from __future__ import annotations

import csv
import os
import tempfile
from pathlib import Path

import pytest

# Ensure the project root is on the path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.metadata import (
    Sample,
    apply_subset_filters,
    build_condition,
    build_design_table,
    create_symlinks,
    validate_metadata_columns,
    write_samples_tsv,
    read_samples_tsv,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def sample_metadata_rows():
    """Simulated metadata rows from an SRA CSV."""
    return [
        {"Run": "SRR925874", "Assay Type": "RNA-Seq", "Genotype": "wild-type", "differentiation": ""},
        {"Run": "SRR925875", "Assay Type": "RNA-Seq", "Genotype": "wild-type", "differentiation": ""},
        {"Run": "SRR925876", "Assay Type": "RNA-Seq", "Genotype": "tet1-/-", "differentiation": ""},
        {"Run": "SRR925877", "Assay Type": "RNA-Seq", "Genotype": "tet1-/-", "differentiation": ""},
        {"Run": "SRR925878", "Assay Type": "RNA-Seq", "Genotype": "tet2-/-", "differentiation": ""},
        {"Run": "SRR925879", "Assay Type": "RNA-Seq", "Genotype": "tet2-/-", "differentiation": ""},
        {"Run": "SRR931799", "Assay Type": "RNA-Seq", "Genotype": "wild-type", "differentiation": "day 3"},
        {"Run": "SRR931800", "Assay Type": "ChIP-Seq", "Genotype": "wild-type", "differentiation": ""},
    ]


@pytest.fixture
def base_config():
    """Minimal config for testing."""
    return {
        "data": {
            "metadata_csv": "/tmp/test.csv",
            "fastq_dir": "/tmp/fastqs",
            "layout": "paired",
        },
        "column_mapping": {
            "run_id_col": "Run",
            "condition_cols": ["Genotype"],
            "condition_map": {
                "wild-type": "wt",
                "tet1-/-": "tet1",
                "tet2-/-": "tet2",
            },
            "replicate_strategy": "sorted_run_id",
        },
        "subset_filters": {
            "default": {
                "Assay Type": "RNA-Seq",
                "differentiation": "",
            },
            "day3": {
                "Assay Type": "RNA-Seq",
                "differentiation": "day 3",
            },
        },
        "active_subset": "default",
    }


# ---------------------------------------------------------------------------
# Tests: subset filtering
# ---------------------------------------------------------------------------

class TestSubsetFiltering:
    def test_default_filter(self, sample_metadata_rows, base_config):
        """Default filter keeps only RNA-Seq undifferentiated."""
        result = apply_subset_filters(sample_metadata_rows, base_config)
        assert len(result) == 6  # 6 RNA-Seq, undifferentiated
        assert all(r["differentiation"] == "" for r in result)
        assert all(r["Assay Type"] == "RNA-Seq" for r in result)

    def test_day3_filter(self, sample_metadata_rows, base_config):
        """Day 3 filter keeps only day 3 RNA-Seq."""
        result = apply_subset_filters(
            sample_metadata_rows, base_config, subset_name="day3"
        )
        assert len(result) == 1
        assert result[0]["Run"] == "SRR931799"

    def test_empty_filter_returns_all(self, sample_metadata_rows, base_config):
        """No matching filter key uses all rows."""
        result = apply_subset_filters(
            sample_metadata_rows, base_config, subset_name="nonexistent"
        )
        assert len(result) == len(sample_metadata_rows)

    def test_no_matches_raises(self, sample_metadata_rows, base_config):
        """Filter matching zero rows raises ValueError."""
        base_config["subset_filters"]["impossible"] = {
            "Assay Type": "ATAC-Seq",
        }
        with pytest.raises(ValueError, match="matched 0 rows"):
            apply_subset_filters(
                sample_metadata_rows, base_config, subset_name="impossible"
            )


# ---------------------------------------------------------------------------
# Tests: condition building
# ---------------------------------------------------------------------------

class TestConditionBuilding:
    def test_single_condition_col(self, base_config):
        row = {"Genotype": "wild-type"}
        cond = build_condition(row, base_config)
        assert cond == "wt"

    def test_unmapped_condition(self, base_config):
        """Values not in condition_map are used as-is (sanitised)."""
        row = {"Genotype": "triple-ko"}
        cond = build_condition(row, base_config)
        assert cond == "triple-ko"

    def test_multi_condition_cols(self, base_config):
        """Multiple condition columns are joined with underscore."""
        base_config["column_mapping"]["condition_cols"] = ["Genotype", "differentiation"]
        row = {"Genotype": "wild-type", "differentiation": "day 3"}
        cond = build_condition(row, base_config)
        assert cond == "wt_day_3"


# ---------------------------------------------------------------------------
# Tests: design table
# ---------------------------------------------------------------------------

class TestDesignTable:
    def test_deterministic_replicates(self, sample_metadata_rows, base_config):
        """Replicate numbers are deterministic based on sorted Run IDs."""
        filtered = apply_subset_filters(sample_metadata_rows, base_config)
        samples = build_design_table(filtered, base_config)

        # 6 samples: 2 wt, 2 tet1, 2 tet2
        assert len(samples) == 6

        # Check that replicates are 1, 2 within each condition
        by_cond = {}
        for s in samples:
            by_cond.setdefault(s.condition, []).append(s)

        for cond, group in by_cond.items():
            reps = [s.replicate for s in group]
            assert reps == [1, 2], f"Condition {cond}: expected [1,2], got {reps}"

    def test_deterministic_across_runs(self, sample_metadata_rows, base_config):
        """Running build_design_table twice gives identical results."""
        filtered = apply_subset_filters(sample_metadata_rows, base_config)
        samples_a = build_design_table(filtered, base_config)
        samples_b = build_design_table(filtered, base_config)

        for a, b in zip(samples_a, samples_b):
            assert a.sample_name == b.sample_name
            assert a.run_id == b.run_id
            assert a.replicate == b.replicate

    def test_sample_names(self, sample_metadata_rows, base_config):
        """Sample names follow condition_rep pattern."""
        filtered = apply_subset_filters(sample_metadata_rows, base_config)
        samples = build_design_table(filtered, base_config)

        names = {s.sample_name for s in samples}
        assert "wt_1" in names
        assert "wt_2" in names
        assert "tet1_1" in names
        assert "tet2_1" in names


# ---------------------------------------------------------------------------
# Tests: samples.tsv round-trip
# ---------------------------------------------------------------------------

class TestSamplesTSV:
    def test_write_and_read(self, sample_metadata_rows, base_config, tmp_path):
        """Writing and reading samples.tsv preserves data."""
        filtered = apply_subset_filters(sample_metadata_rows, base_config)
        samples = build_design_table(filtered, base_config)

        tsv_path = tmp_path / "samples.tsv"
        write_samples_tsv(samples, tsv_path)

        loaded = read_samples_tsv(tsv_path)
        assert len(loaded) == len(samples)

        for orig, back in zip(samples, loaded):
            assert orig.sample_name == back.sample_name
            assert orig.condition == back.condition
            assert orig.run_id == back.run_id
            assert orig.replicate == back.replicate


# ---------------------------------------------------------------------------
# Tests: validation and edge cases (exam-day hardening)
# ---------------------------------------------------------------------------

class TestColumnValidation:
    def test_missing_run_id_col(self, sample_metadata_rows, base_config):
        """Config references a run_id_col that doesn't exist in the CSV."""
        base_config["column_mapping"]["run_id_col"] = "NonExistentCol"
        with pytest.raises(ValueError, match="run_id_col 'NonExistentCol' not found"):
            validate_metadata_columns(sample_metadata_rows, base_config)

    def test_missing_condition_col(self, sample_metadata_rows, base_config):
        """Config references a condition column that doesn't exist."""
        base_config["column_mapping"]["condition_cols"] = ["MissingCol"]
        with pytest.raises(ValueError, match="condition_col 'MissingCol' not found"):
            validate_metadata_columns(sample_metadata_rows, base_config)


class TestEmptyCondition:
    def test_empty_condition_raises(self, base_config):
        """Condition mapping that produces an empty name should raise."""
        base_config["column_mapping"]["condition_cols"] = ["Genotype"]
        base_config["column_mapping"]["condition_map"] = {"wild-type": ""}
        row = {"Run": "SRR000001", "Genotype": "wild-type"}
        with pytest.raises(ValueError, match="Empty condition"):
            build_condition(row, base_config)


class TestUnequalReplicates:
    def test_three_plus_two(self, base_config):
        """Unequal replicate counts produce correct numbering."""
        rows = [
            {"Run": "SRR1", "Assay Type": "RNA-Seq", "Genotype": "wild-type", "differentiation": ""},
            {"Run": "SRR2", "Assay Type": "RNA-Seq", "Genotype": "wild-type", "differentiation": ""},
            {"Run": "SRR3", "Assay Type": "RNA-Seq", "Genotype": "wild-type", "differentiation": ""},
            {"Run": "SRR4", "Assay Type": "RNA-Seq", "Genotype": "tet1-/-", "differentiation": ""},
            {"Run": "SRR5", "Assay Type": "RNA-Seq", "Genotype": "tet1-/-", "differentiation": ""},
        ]
        samples = build_design_table(rows, base_config)
        assert len(samples) == 5

        by_cond = {}
        for s in samples:
            by_cond.setdefault(s.condition, []).append(s)

        assert [s.replicate for s in by_cond["wt"]] == [1, 2, 3]
        assert [s.replicate for s in by_cond["tet1"]] == [1, 2]

        names = [s.sample_name for s in samples]
        assert len(names) == len(set(names)), "All sample names must be unique"


class TestDuplicateRunId:
    def test_duplicate_raises(self, base_config):
        """Duplicate run IDs in filtered metadata should raise."""
        rows = [
            {"Run": "SRR1", "Genotype": "wild-type"},
            {"Run": "SRR1", "Genotype": "wild-type"},
        ]
        with pytest.raises(ValueError, match="Duplicate run IDs"):
            build_design_table(rows, base_config)


class TestGSE104853Integration:
    def test_mirna125a_design(self):
        """End-to-end: GSE104853-style metadata produces expected design."""
        rows = [
            {"Run": "SRR6162944", "Assay Type": "RNA-Seq",
             "treatment": "transfected with negative control microRNA mimic"},
            {"Run": "SRR6162945", "Assay Type": "RNA-Seq",
             "treatment": "transfected with negative control microRNA mimic"},
            {"Run": "SRR6162946", "Assay Type": "RNA-Seq",
             "treatment": "transfected with negative control microRNA mimic"},
            {"Run": "SRR6162947", "Assay Type": "RNA-Seq",
             "treatment": "transfected with miR-125a-5p mimic"},
            {"Run": "SRR6162948", "Assay Type": "RNA-Seq",
             "treatment": "transfected with miR-125a-5p mimic"},
            {"Run": "SRR6162949", "Assay Type": "RNA-Seq",
             "treatment": "transfected with miR-125a-5p mimic"},
        ]
        cfg = {
            "data": {
                "metadata_csv": "/data/GSE104853/metadata.csv",
                "fastq_dir": "/data/GSE104853",
                "layout": "paired",
            },
            "column_mapping": {
                "run_id_col": "Run",
                "condition_cols": ["treatment"],
                "condition_map": {
                    "transfected with negative control microRNA mimic": "Ctrl",
                    "transfected with miR-125a-5p mimic": "miRNA125a",
                },
                "replicate_strategy": "sorted_run_id",
            },
            "subset_filters": {"default": {"Assay Type": "RNA-Seq"}},
            "active_subset": "default",
        }

        filtered = apply_subset_filters(rows, cfg)
        assert len(filtered) == 6

        samples = build_design_table(filtered, cfg)
        assert len(samples) == 6

        names = sorted(s.sample_name for s in samples)
        assert names == ["ctrl_1", "ctrl_2", "ctrl_3",
                         "mirna125a_1", "mirna125a_2", "mirna125a_3"]

        conditions = set(s.condition for s in samples)
        assert conditions == {"ctrl", "mirna125a"}
