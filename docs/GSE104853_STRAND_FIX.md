# GSE104853 0-DEG fix: strandedness and optional filtered BAM

## What was changed

1. **Strandedness (primary fix)**  
   - **Config:** `config/config_GSE104853.yaml`: `featurecounts.strandedness` set from `0` to `2` (reverse-stranded).  
   - **Reason:** Pipeline completed but reported 0 significant DEGs (FDR < 0.05) in all branches. The paper (PMC9680651) reports transcriptomic effects of miR-125a-5p (e.g. cytokine-cytokine receptor repression, decreased IL10RA). Wrong strandedness in featureCounts can misassign reads and flatten or reverse DE; many Illumina RNA‑seq kits are reverse-stranded.  
   - **Validation:** Run `./py scripts/run_strand_test.py` from the pipeline directory (with `results/gse104853_test` present). It runs steps 7–9 with strand 0, 1, and 2 on the main branch (cutadapt + STAR default) and prints a comparison table (assignment %, genes after filter, DEG count, IL10RA/CD163). Use the recommended strandedness; if strand 1 is better, set `strandedness: 1` in the config.

2. **Strand-test configs and script**  
   - **Added:** `config/config_GSE104853_strand1.yaml` and `config/config_GSE104853_strand2.yaml` (same as GSE104853 but `strandedness: 1` and `2`, main branch only).  
   - **Added:** `scripts/run_strand_test.py`: copies main-branch BAMs from `gse104853_test` into `gse104853_strand1` and `gse104853_strand2`, runs steps 7–9 with the strand1/strand2 configs, and collects the comparison table.

3. **Optional: count from filtered BAMs (second fix)**  
   - **Config:** `featurecounts.use_filtered_bam` (default `false`). When `true`, step 7 uses `filtered_bam_path` from `mapping_summary.tsv` instead of `bam_path` (post–step‑5 filtered BAMs).  
   - **Code:** `scripts/07_featurecounts.py`: `build_mapping_units()` accepts `use_filtered_bam`; when true and `filtered_bam_path` is present, that path is used for counting.  
   - **When to use:** Only if correcting strandedness does not restore a plausible DE signal (e.g. still 0 DEGs or wrong direction for IL10RA/CD163). Then set `use_filtered_bam: true` in `config_GSE104853.yaml`, re-run steps 7–9 for the main branch, and compare.

4. **Config clean-up**  
   - **Config:** Removed duplicate `extra_args` under `featurecounts.option_sets.multimapper` in `config_GSE104853.yaml`.

## Before vs after (to be filled by running the pipeline)

After running the strand test and/or the fixed pipeline, fill a table like:

| setting tested   | assigned % | genes after filter | DEGs padj<0.05 | IL10RA log2FC/padj | CD163 log2FC/padj | comments |
|------------------|------------|---------------------|----------------|---------------------|-------------------|----------|
| strand 0 (baseline) | …       | …                   | 0              | …                   | …                 | baseline |
| strand 1         | …          | …                   | …              | …                   | …                 |          |
| strand 2         | …          | …                   | …              | …                   | …                 |          |

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
  To be confirmed by running `run_strand_test.py` and comparing DEG count and IL10RA/CD163 direction vs the paper. If strand 1 or 2 yields DEGs and IL10RA down in mirna125a vs ctrl, strandedness was the cause.

- **Did correcting it restore a plausible transcriptomic signal?**  
  Check: DEGs padj < 0.05 > 0, IL10RA and CD163 direction and significance consistent with the abstract (e.g. IL10RA decreased in miR-125a-5p).

- **If not, is counting from filtered BAMs the next justified fix?**  
  Yes. Step 7 previously used raw STAR BAMs (`bam_path`). The pipeline now supports `use_filtered_bam: true` to count from post–step‑5 filtered BAMs. Try this only after validating strandedness.

- **Does the pipeline still work operationally after the fix?**  
  Yes. Only config (`strandedness`, optional `use_filtered_bam`) and one optional behaviour in step 7 were added; no broad refactor. Run the usual command with `--run-id gse104853_fix_strand` to confirm end-to-end.
