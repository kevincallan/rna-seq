"""
Tests for naming conventions and symlink generation.

Run with:  pytest tests/ -v
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.metadata import Sample, create_symlinks
from src.utils import safe_name


# ---------------------------------------------------------------------------
# Tests: safe_name
# ---------------------------------------------------------------------------

class TestSafeName:
    def test_basic(self):
        assert safe_name("wild-type") == "wild-type"

    def test_spaces(self):
        assert safe_name("Some Condition") == "some_condition"

    def test_slashes(self):
        """Slashes are converted to underscores by safe_name."""
        result = safe_name("tet1-/-")
        assert "/" not in result
        assert result == "tet1-_-"

    def test_special_chars(self):
        result = safe_name("condition (A) + B")
        assert "(" not in result
        assert ")" not in result
        assert "+" not in result

    def test_multiple_spaces(self):
        result = safe_name("  spaced   out  ")
        assert "  " not in result
        assert result == "spaced_out"

    def test_empty(self):
        assert safe_name("") == ""

    def test_numeric(self):
        assert safe_name("day3") == "day3"

    def test_case_lowered(self):
        assert safe_name("WT_Day3") == "wt_day3"


# ---------------------------------------------------------------------------
# Tests: symlink naming
# ---------------------------------------------------------------------------

class TestSymlinkNaming:
    def test_paired_symlink_names(self, tmp_path):
        """Paired-end symlinks follow condition_rep_mate.fastq.gz pattern."""
        # Create fake source files
        src_dir = tmp_path / "source"
        src_dir.mkdir()
        (src_dir / "SRR001_1.fastq.gz").touch()
        (src_dir / "SRR001_2.fastq.gz").touch()

        sample = Sample(
            run_id="SRR001",
            condition="wt",
            replicate=1,
            layout="paired",
            r1=str(src_dir / "SRR001_1.fastq.gz"),
            r2=str(src_dir / "SRR001_2.fastq.gz"),
        )

        work_dir = tmp_path / "work"
        create_symlinks([sample], work_dir)

        link_dir = work_dir / "fastq_links"
        assert (link_dir / "wt_1_1.fastq.gz").is_symlink()
        assert (link_dir / "wt_1_2.fastq.gz").is_symlink()

    def test_single_symlink_names(self, tmp_path):
        """Single-end symlinks follow condition_rep.fastq.gz pattern."""
        src_dir = tmp_path / "source"
        src_dir.mkdir()
        (src_dir / "SRR002.fastq.gz").touch()

        sample = Sample(
            run_id="SRR002",
            condition="ko",
            replicate=2,
            layout="single",
            r1=str(src_dir / "SRR002.fastq.gz"),
        )

        work_dir = tmp_path / "work"
        create_symlinks([sample], work_dir)

        link_dir = work_dir / "fastq_links"
        assert (link_dir / "ko_2.fastq.gz").is_symlink()
        # Should NOT have mate 2
        assert not (link_dir / "ko_2_2.fastq.gz").exists()

    def test_symlinks_are_idempotent(self, tmp_path):
        """Running create_symlinks twice doesn't error or duplicate."""
        src_dir = tmp_path / "source"
        src_dir.mkdir()
        (src_dir / "SRR003_1.fastq.gz").touch()
        (src_dir / "SRR003_2.fastq.gz").touch()

        sample = Sample(
            run_id="SRR003",
            condition="treated",
            replicate=1,
            layout="paired",
            r1=str(src_dir / "SRR003_1.fastq.gz"),
            r2=str(src_dir / "SRR003_2.fastq.gz"),
        )

        work_dir = tmp_path / "work"
        create_symlinks([sample], work_dir)
        create_symlinks([sample], work_dir)  # should not raise

        link_dir = work_dir / "fastq_links"
        assert (link_dir / "treated_1_1.fastq.gz").is_symlink()

    def test_multiple_samples_no_collision(self, tmp_path):
        """Multiple samples in same condition get distinct symlinks."""
        src_dir = tmp_path / "source"
        src_dir.mkdir()
        for rid in ["SRR010", "SRR011"]:
            for mate in ["1", "2"]:
                (src_dir / f"{rid}_{mate}.fastq.gz").touch()

        samples = [
            Sample("SRR010", "wt", 1, "paired",
                   str(src_dir / "SRR010_1.fastq.gz"),
                   str(src_dir / "SRR010_2.fastq.gz")),
            Sample("SRR011", "wt", 2, "paired",
                   str(src_dir / "SRR011_1.fastq.gz"),
                   str(src_dir / "SRR011_2.fastq.gz")),
        ]

        work_dir = tmp_path / "work"
        create_symlinks(samples, work_dir)

        link_dir = work_dir / "fastq_links"
        assert (link_dir / "wt_1_1.fastq.gz").is_symlink()
        assert (link_dir / "wt_2_1.fastq.gz").is_symlink()
        # Ensure they point to different targets
        t1 = os.readlink(link_dir / "wt_1_1.fastq.gz")
        t2 = os.readlink(link_dir / "wt_2_1.fastq.gz")
        assert t1 != t2


# ---------------------------------------------------------------------------
# Tests: Sample class
# ---------------------------------------------------------------------------

class TestSampleClass:
    def test_sample_name_format(self):
        s = Sample("SRR123", "ctrl", 3)
        assert s.sample_name == "ctrl_3"

    def test_repr(self):
        s = Sample("SRR123", "ctrl", 3)
        assert "ctrl_3" in repr(s)
        assert "SRR123" in repr(s)
