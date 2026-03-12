# GSE104853 0-DEG fix: strandedness and optional filtered BAM

## What was changed

1. **Strandedness (tested, reverted to 0)**  
   - **Config:** `config/config_GSE104853.yaml`: `featurecounts.strandedness` kept at `0` (unstranded).  
   - **Reason:** Pipeline reported 0 significant DEGs; wrong strandedness was a candidate. After running `run_strand_test.py`, **strand 0 gave ~42% assigned vs ~22% for strand 1/2**, supporting unstranded libraries. Strand 1 and 2 did not yield any DEGs either. Config was reverted to `strandedness: 0` with a comment referencing the strand test.

2. **Strand-test configs and script**  
   - **Added:** `config/config_GSE104853_strand1.yaml` and `config/config_GSE104853_strand2.yaml` (same as GSE104853 but `strandedness: 1` and `2`, main branch only).  
   - **Added:** `scripts/run_strand_test.py`: copies main-branch BAMs from `gse104853_test` into `gse104853_strand1` and `gse104853_strand2`, runs steps 7–9 with the strand1/strand2 configs, and collects the comparison table.

3. **Optional: count from filtered BAMs (second fix)**  
   - **Config:** `featurecounts.use_filtered_bam` (default `false`). When `true`, step 7 uses `filtered_bam_path` from `mapping_summary.tsv` instead of `bam_path` (post–step‑5 filtered BAMs).  
   - **Code:** `scripts/07_featurecounts.py`: `build_mapping_units()` accepts `use_filtered_bam`; when true and `filtered_bam_path` is present, that path is used for counting.  
   - **When to use:** Only if correcting strandedness does not restore a plausible DE signal (e.g. still 0 DEGs or wrong direction for IL10RA/CD163). Then set `use_filtered_bam: true` in `config_GSE104853.yaml`, re-run steps 7–9 for the main branch, and compare.

4. **Config clean-up**  
   - **Config:** Removed duplicate `extra_args` under `featurecounts.option_sets.multimapper` in `config_GSE104853.yaml`.

## Strand test outcome (2026-03-12)

| setting tested   | assigned % | genes after filter | DEGs padj<0.05 | IL10RA log2FC/padj | CD163 log2FC/padj | comments |
|------------------|------------|---------------------|----------------|---------------------|-------------------|----------|
| strand 0 (baseline) | 42.2%   | 6694                | 0              | -0.798 / 1.0        | -1.188 / 1.0      | baseline |
| strand 1         | 22.3%      | 4660                | 0              | -0.601 / 1.0        | (filtered out)    |          |
| strand 2         | 21.8%      | 4693                | 0              | -1.057 / 1.0        | (filtered out)   |          |

**Interpretation:** Assignment rate is highest with strand 0 (unstranded), so the library is treated as unstranded; config stays at `strandedness: 0`. Strandedness was **not** the cause of 0 DEGs. IL10RA and CD163 are in the expected direction (negative log2FC in mirna125a vs ctrl) but none reach FDR &lt; 0.05. Next options: try `use_filtered_bam: true` and re-run steps 7–9, or attribute 0 DEGs to limited power (n=3 per group) and use effect sizes/direction for interpretation.

## How to run the fixed analysis

1. **Strand test (optional but recommended)**  
   From the pipeline directory, with `results/gse104853_test` already present:
   ```bash
   ./py scripts/run_strand_test.py
   ```
   Update `featurecounts.strandedness` in `config/config_GSE104853.yaml` to the recommended value (0, 1, or 2) if it differs from the current setting.

2. **Fixed run (new run-id)**  
   ```bash
   ./py scripts/run_pipeline.py \
     --config config/config_GSE104853.yaml \
     --dataset GSE104853 \
     --species human \
     --run-id gse104853_fix_strand \
     run --methods cutadapt --steps 7 8 9 10 11
   ```
   This uses the existing BAMs from a previous full run (e.g. `gse104853_test`); ensure that the results directory for that run still exists, or that you point to it (e.g. by copying/symlinking the main branch into `results/gse104853_fix_strand` and running steps 7–9 only). For a clean fixed run from scratch, run all steps with `--run-id gse104853_fix_strand` and `--methods cutadapt` (and no `--steps` or steps 0–11).

3. **If strandedness alone is not enough**  
   Set in `config_GSE104853.yaml`:
   ```yaml
   featurecounts:
     use_filtered_bam: true
   ```
   Re-run steps 7–9 (and 10–11 if desired) for the main branch and compare DEG counts and IL10RA/CD163 again.

## Answers to the validation questions

- **Was strandedness the actual cause?**  
  No. The strand test showed strand 0 (unstranded) has the highest assignment (~42% vs ~22% for 1/2). Strand 1 and 2 did not yield any DEGs. Config was reverted to `strandedness: 0`.

- **Did correcting it restore a plausible transcriptomic signal?**  
  Strand change did not change significance (0 DEGs in all). Direction is plausible: IL10RA and CD163 are negative (down in mirna125a vs ctrl), consistent with the paper; padj remains high (limited power with n=3).

- **If not, is counting from filtered BAMs the next justified fix?**  
  Yes. Set `use_filtered_bam: true` in `config_GSE104853.yaml`, re-run steps 7–9 for the main branch, and compare. Otherwise treat 0 DEGs as a power/sample-size limitation and use effect direction (e.g. IL10RA, CD163) for interpretation.

- **Does the pipeline still work operationally after the fix?**  
  Yes. Strand test ran successfully; only config and optional step‑7 behaviour were added. Pipeline is operational.
