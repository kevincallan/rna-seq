"""Tests for pipeline execution profiles."""

from __future__ import annotations

from pathlib import Path

import pytest

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.run_pipeline import PROFILES, STEPS, ALL_STEPS


class TestProfileDefinitions:
    def test_primary_excludes_bigwig(self):
        bigwig_step = None
        for num, (module, desc) in STEPS.items():
            if "bigwig" in module.lower():
                bigwig_step = num
                break
        assert bigwig_step is not None, "BigWig step should exist in STEPS"
        assert bigwig_step not in PROFILES["primary"], \
            "primary profile should not include BigWig step"

    def test_primary_bw_includes_bigwig(self):
        bigwig_step = None
        for num, (module, desc) in STEPS.items():
            if "bigwig" in module.lower():
                bigwig_step = num
                break
        assert bigwig_step in PROFILES["primary_bw"]

    def test_full_includes_all_steps(self):
        assert set(PROFILES["full"]) == set(ALL_STEPS)

    def test_primary_includes_essential_steps(self):
        primary = set(PROFILES["primary"])
        for num, (module, desc) in STEPS.items():
            if "validate" in desc.lower() or "prepare" in desc.lower():
                assert num in primary, f"Step {num} ({desc}) should be in primary"
            if "featurecounts" in module.lower():
                assert num in primary, f"featureCounts step should be in primary"
            if "deseq2" in module.lower():
                assert num in primary, f"DESeq2 step should be in primary"
            if "report" in module.lower():
                assert num in primary, f"Report step should be in primary"
            if "compare" in module.lower():
                assert num in primary, f"Compare step should be in primary"


class TestStepOrdering:
    def test_bigwig_after_deseq2(self):
        deseq2_step = None
        bigwig_step = None
        for num, (module, desc) in STEPS.items():
            if "deseq2" in module.lower():
                deseq2_step = num
            if "bigwig" in module.lower():
                bigwig_step = num
        assert deseq2_step is not None and bigwig_step is not None
        assert bigwig_step > deseq2_step, \
            "BigWig step must run after DESeq2 so size factors are available"

    def test_featurecounts_before_filter(self):
        fc_step = None
        filter_step = None
        for num, (module, desc) in STEPS.items():
            if "featurecounts" in module.lower():
                fc_step = num
            if "filter" in module.lower():
                filter_step = num
        assert fc_step is not None and filter_step is not None
        assert fc_step < filter_step

    def test_filter_before_deseq2(self):
        filter_step = None
        deseq2_step = None
        for num, (module, desc) in STEPS.items():
            if "filter" in module.lower():
                filter_step = num
            if "deseq2" in module.lower():
                deseq2_step = num
        assert filter_step is not None and deseq2_step is not None
        assert filter_step < deseq2_step

    def test_compare_before_report(self):
        compare_step = None
        report_step = None
        for num, (module, desc) in STEPS.items():
            if "compare" in module.lower():
                compare_step = num
            if "report" in module.lower():
                report_step = num
        assert compare_step is not None and report_step is not None
        assert compare_step < report_step


class TestProfileStepsAreValid:
    @pytest.mark.parametrize("profile_name", PROFILES.keys())
    def test_all_profile_steps_exist(self, profile_name):
        for step_num in PROFILES[profile_name]:
            assert step_num in STEPS, \
                f"Profile '{profile_name}' references non-existent step {step_num}"
