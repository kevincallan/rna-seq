#!/bin/bash
# =============================================================================
# Colab Setup Script
# =============================================================================
# Run this ONCE after cloning the repo on Google Colab.
#
# Option A - Ephemeral (lost on disconnect):
#   !git clone https://github.com/kevincallan/rna-seq.git
#   %cd rna-seq
#   !bash scripts/setup_colab.sh
#
# Option B - Persist repo and outputs on Google Drive (survives disconnect):
#   from google.colab import drive
#   drive.mount('/content/drive')
#   %cd /content/drive/MyDrive
#   !git clone https://github.com/kevincallan/rna-seq.git
#   %cd rna-seq
#   !bash scripts/setup_colab.sh
#   # References and data will be under /content/... (re-run setup after reconnect).
#   # To save results to Drive after a run: cp -r results /content/drive/MyDrive/rna-seq_results
#
# Total time: ~5-10 min.  Disk: ~5 GB for genome index + test FASTQs.
# =============================================================================

set -e
echo "========================================="
echo "  RNA-seq Pipeline -- Colab Setup"
echo "========================================="

# -----------------------------------------------------------------
# 1. Install Python packages
# -----------------------------------------------------------------
echo ""
echo ">>> Installing Python packages..."
pip install -q pyyaml pandas numpy scipy matplotlib scikit-learn \
    pydeseq2 pysam HTSeq cutadapt multiqc deeptools markdown anndata

echo "    Python packages OK."

# -----------------------------------------------------------------
# 2. Install STAR (the only required compiled binary)
# -----------------------------------------------------------------
echo ""
echo ">>> Installing STAR aligner..."
# conda is available on Colab via condacolab
if ! command -v STAR &> /dev/null; then
    # Install via apt (fastest on Colab)
    apt-get update -qq && apt-get install -y -qq rna-star 2>/dev/null || {
        # Fallback: download pre-built binary
        echo "    apt failed, downloading STAR binary..."
        STAR_VERSION="2.7.11b"
        wget -q "https://github.com/alexdobin/STAR/releases/download/${STAR_VERSION}/STAR_${STAR_VERSION}.zip" \
            -O /tmp/star.zip
        cd /tmp && unzip -q star.zip && cd -
        cp /tmp/STAR_${STAR_VERSION}/Linux_x86_64_static/STAR /usr/local/bin/
        chmod +x /usr/local/bin/STAR
    }
fi
echo "    STAR: $(STAR --version 2>&1 | head -1 || echo 'installed')"

# Also install samtools + FastQC as backup (small, fast)
apt-get install -y -qq samtools fastqc 2>/dev/null || true

# -----------------------------------------------------------------
# 3. Download genome index + GTF (mouse mm39, chr19 subset for speed)
# -----------------------------------------------------------------
echo ""
echo ">>> Setting up reference genome..."
REFDIR="/content/references/mm39"
mkdir -p "$REFDIR"

# We'll generate a small STAR index from chr19 for testing.
# For full analysis, you'd use a complete genome index.
if [ ! -d "$REFDIR/STAR" ]; then
    echo "    Downloading mm39 chr19 FASTA + GTF..."
    # Download chr19 FASTA from UCSC
    wget -q "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/chromosomes/chr19.fa.gz" \
        -O "$REFDIR/chr19.fa.gz"
    gunzip -f "$REFDIR/chr19.fa.gz"

    # Download GTF from Ensembl (GRCm39)
    wget -q "https://ftp.ensembl.org/pub/release-113/gtf/mus_musculus/Mus_musculus.GRCm39.113.chr.gtf.gz" \
        -O "$REFDIR/mm39.gtf.gz"
    gunzip -f "$REFDIR/mm39.gtf.gz"

    # Filter GTF to chr19 only (for the small index)
    grep "^#\|^19" "$REFDIR/mm39.gtf" > "$REFDIR/chr19.gtf"
    # Also fix chromosome naming (Ensembl uses "19", UCSC uses "chr19")
    sed -i 's/^19/chr19/' "$REFDIR/chr19.gtf"

    echo "    Building STAR index for chr19 (takes ~2 min)..."
    mkdir -p "$REFDIR/STAR"
    STAR --runMode genomeGenerate \
         --genomeDir "$REFDIR/STAR" \
         --genomeFastaFiles "$REFDIR/chr19.fa" \
         --sjdbGTFfile "$REFDIR/chr19.gtf" \
         --genomeSAindexNbases 11 \
         --runThreadN 2 \
         --outTmpDir /tmp/star_tmp \
         2>&1 | tail -3

    echo "    STAR index built."
else
    echo "    STAR index already exists."
fi

# -----------------------------------------------------------------
# 4. Download test FASTQs from SRA (small subset)
# -----------------------------------------------------------------
echo ""
echo ">>> Setting up test FASTQ data..."
FQDIR="/content/fastqs"
mkdir -p "$FQDIR"

# Download the 5% subset used in the course practicals
# These are ~30-70 MB each -- manageable on Colab
# We'll use the undifferentiated subset (6 samples, paired-end)
RUNS=("SRR925874" "SRR925875" "SRR925876" "SRR925877" "SRR925878" "SRR925879")

# Check if already downloaded
NEED_DOWNLOAD=false
for run in "${RUNS[@]}"; do
    if [ ! -f "$FQDIR/${run}_1.fastq.gz" ]; then
        NEED_DOWNLOAD=true
        break
    fi
done

if [ "$NEED_DOWNLOAD" = true ]; then
    echo "    Installing SRA toolkit..."
    pip install -q ffq 2>/dev/null || true

    # Try fasterq-dump first, fall back to wget from ENA
    for run in "${RUNS[@]}"; do
        if [ -f "$FQDIR/${run}_1.fastq.gz" ]; then
            echo "    $run already exists, skipping."
            continue
        fi
        echo "    Downloading $run from ENA..."
        # ENA direct download (faster than SRA toolkit on Colab)
        wget -q "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${run:0:6}/$run/${run}_1.fastq.gz" \
            -O "$FQDIR/${run}_1.fastq.gz" 2>/dev/null || \
        wget -q "https://ftp.sra.ebi.ac.uk/vol1/fastq/${run:0:6}/$run/${run}_1.fastq.gz" \
            -O "$FQDIR/${run}_1.fastq.gz" 2>/dev/null || \
            echo "    WARNING: Could not download $run mate 1"

        wget -q "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${run:0:6}/$run/${run}_2.fastq.gz" \
            -O "$FQDIR/${run}_2.fastq.gz" 2>/dev/null || \
        wget -q "https://ftp.sra.ebi.ac.uk/vol1/fastq/${run:0:6}/$run/${run}_2.fastq.gz" \
            -O "$FQDIR/${run}_2.fastq.gz" 2>/dev/null || \
            echo "    WARNING: Could not download $run mate 2"
    done
else
    echo "    FASTQs already present."
fi

echo "    FASTQs in $FQDIR:"
ls -lh "$FQDIR"/*.fastq.gz 2>/dev/null || echo "    (no FASTQs found)"

# -----------------------------------------------------------------
# 5. Create Colab metadata CSV
# -----------------------------------------------------------------
echo ""
echo ">>> Creating metadata CSV..."
METADIR="/content/fastqs"
cat > "$METADIR/metadata.csv" << 'METADATA'
Run,Assay Type,Genotype,differentiation
SRR925874,RNA-Seq,wild-type,
SRR925875,RNA-Seq,wild-type,
SRR925876,RNA-Seq,tet1-/-,
SRR925877,RNA-Seq,tet1-/-,
SRR925878,RNA-Seq,tet2-/-,
SRR925879,RNA-Seq,tet2-/-,
METADATA
echo "    Metadata CSV created at $METADIR/metadata.csv"

# -----------------------------------------------------------------
# 6. Write Colab config (in repo so it works when repo is on Drive)
# -----------------------------------------------------------------
REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
echo ""
echo ">>> Writing Colab config..."
cat > "$REPO_ROOT/config/config_colab.yaml" << 'COLABCFG'
# Config for Google Colab -- paths adjusted for /content/
project:
  name: "GSE48519_colab"
  results_dir: "results"
  work_dir: "work"
  logs_dir: "logs"
  threads: 2

data:
  metadata_csv: "/content/fastqs/metadata.csv"
  fastq_dir: "/content/fastqs"
  layout: "paired"

column_mapping:
  run_id_col: "Run"
  condition_cols:
    - "Genotype"
  condition_map:
    "wild-type": "wt"
    "tet1-/-": "tet1"
    "tet2-/-": "tet2"
  replicate_strategy: "sorted_run_id"

subset_filters:
  default:
    "Assay Type": "RNA-Seq"
    "differentiation": ""

active_subset: "default"

references:
  genome_index: "/content/references/mm39/STAR"
  gtf: "/content/references/mm39/chr19.gtf"
  genome_fasta: "/content/references/mm39/chr19.fa"

tools:
  fastqc: "fastqc"
  multiqc: "multiqc"
  cutadapt: "cutadapt"
  fastp: "fastp"
  trimmomatic: "trimmomatic"
  star: "STAR"
  samtools: "samtools"
  featurecounts: "featureCounts"
  bamcoverage: "bamCoverage"

star_params:
  outSAMtype: "BAM SortedByCoordinate"
  readFilesCommand: "zcat"
  outFilterMultimapNmax: 20
  extra_args: ""

trimming:
  none:
    enabled: true
  cutadapt:
    enabled: true
    adapter_fwd: "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
    adapter_rev: "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    quality: 20
    min_length: 25
    extra_args: ""
  fastp:
    enabled: false
  trimmomatic:
    enabled: false

featurecounts:
  backend: "htseq"               # Pure Python -- no subread binary needed
  strandedness: 2
  feature_type: "exon"
  attribute: "gene_id"
  option_sets:
    default:
      label: "Default"
      B: false
      P: false
      C: false
      Q: 0
      extra_args: ""

filtering:
  min_max_count: "auto"

deseq2:
  fdr_threshold: 0.05
  lfc_threshold: 0.0
  method: "pydeseq2"             # Pure Python -- no R needed
  reference_level: "wt"

comparisons:
  - name: "tet1_vs_wt"
    numerator: "tet1"
    denominator: "wt"
  - name: "tet2_vs_wt"
    numerator: "tet2"
    denominator: "wt"

bigwig:
  normalization: "CPM"
  mapq_filter: 255
  bin_size: 50
  use_deseq2_sizefactors: true
COLABCFG

echo "    Colab config at $REPO_ROOT/config/config_colab.yaml"

# -----------------------------------------------------------------
# 7. Verify setup
# -----------------------------------------------------------------
echo ""
echo "========================================="
echo "  Setup complete!"
echo "========================================="
echo ""
echo "  Python packages:  OK"
echo "  STAR:             $(which STAR 2>/dev/null || echo 'CHECK MANUALLY')"
echo "  Reference:        $REFDIR/STAR/"
echo "  FASTQs:           $FQDIR/"
echo "  Metadata:         $FQDIR/metadata.csv"
echo "  Config:           $REPO_ROOT/config/config_colab.yaml"
echo ""
echo "  To run the pipeline:"
echo "    cd $REPO_ROOT"
echo "    ./py scripts/run_pipeline.py --config config/config_colab.yaml run"
echo ""
echo "  To run just QC + mapping (faster test):"
echo "    ./py scripts/run_pipeline.py --config config/config_colab.yaml run --steps 0 1 2 5"
echo ""
echo "  To run with only 'none' trimming (fastest):"
echo "    ./py scripts/run_pipeline.py --config config/config_colab.yaml run --methods none"
echo ""
