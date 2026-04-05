# RNA-seq Nextflow Pipeline

[![CI](https://github.com/Ekin-Kahraman/rnaseq-nextflow-pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/Ekin-Kahraman/rnaseq-nextflow-pipeline/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Nextflow](https://img.shields.io/badge/Nextflow-%E2%89%A524.0-brightgreen)](https://www.nextflow.io/)

End-to-end bulk RNA-seq pipeline in Nextflow DSL2: raw FASTQ reads through quality control, adapter trimming, genome alignment, gene quantification, differential expression, and aggregated QC reporting. Every step runs in its own Docker container.

Demonstrated on the [Himes et al. (2014)](https://doi.org/10.1371/journal.pone.0099625) airway smooth muscle dataset — dexamethasone-treated vs untreated human airway cells. For covariate-adjusted analysis on a larger COVID-19 cohort, see [bulk-rnaseq-differential-expression](https://github.com/Ekin-Kahraman/bulk-rnaseq-differential-expression).

## Workflow

```
FASTQ (paired-end)
    │
    ▼
 FastQC ──────── Raw read quality assessment
    │
    ▼
 fastp ────────── Adapter trimming, quality filtering
    │
    ▼
 HISAT2 ──────── Align to GRCh38 reference genome
    │
    ▼
 samtools ─────── Sort and index BAM
    │
    ▼
 featureCounts ── Gene-level quantification (Gencode v38)
    │
    ▼
 DESeq2 ──────── Differential expression + PCA + volcano plot
    │
    ▼
 MultiQC ─────── Aggregate QC report across all samples
```

## Processes

| Process | Tool | Container |
|---------|------|-----------|
| FASTQC_RAW | FastQC 0.12.1 | `quay.io/biocontainers/fastqc` |
| FASTP | fastp 0.23.4 | `quay.io/biocontainers/fastp` |
| HISAT2_ALIGN | HISAT2 2.2.1 | `quay.io/biocontainers/hisat2` |
| SAMTOOLS_SORT | samtools 1.21 | `quay.io/biocontainers/samtools` |
| FEATURECOUNTS | Subread 2.0.6 | `quay.io/biocontainers/subread` |
| DESEQ2 | DESeq2 1.42 + ggplot2 | `quay.io/biocontainers/bioconductor-deseq2` |
| MULTIQC | MultiQC 1.27 | `quay.io/biocontainers/multiqc` |

All containers sourced from [BioContainers](https://biocontainers.pro/).

## Dataset

**Himes et al. (2014)** — RNA-seq of human airway smooth muscle cells treated with dexamethasone (a glucocorticoid anti-inflammatory). GEO accession [GSE52778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778). This dataset is used in the [DESeq2 vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) and the [Bioconductor RNA-seq workflow](https://www.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html).

| Sample | Accession | Condition | Donor |
|--------|-----------|-----------|-------|
| N61311_untreated | SRR1039508 | untreated | N61311 |
| N61311_Dex | SRR1039509 | dexamethasone | N61311 |
| N052611_untreated | SRR1039512 | untreated | N052611 |
| N052611_Dex | SRR1039513 | dexamethasone | N052611 |

## Quick Start

**Prerequisites:** [Nextflow](https://www.nextflow.io/) (>=24.0), [Docker](https://www.docker.com/), Java (>=11)

### Test data (synthetic, ~2 minutes)

```bash
git clone https://github.com/Ekin-Kahraman/rnaseq-nextflow-pipeline.git
cd rnaseq-nextflow-pipeline
python test/create_test_data.py
nextflow run main.nf -profile test,docker \
    --genome_index "$(pwd)/test/genome" \
    --gtf "$(pwd)/test/genes.gtf"
```

### Real data (airway dataset)

```bash
# 1. Download HISAT2 GRCh38 index (~4GB)
mkdir -p genome && cd genome
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
tar xzf grch38_genome.tar.gz

# 2. Download Gencode v38 GTF
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
gunzip gencode.v38.annotation.gtf.gz
cd ..

# 3. Download FASTQ files from ENA (see assets/samplesheet.csv for accessions)
mkdir -p data
# Example: wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/008/SRR1039508/SRR1039508_1.fastq.gz -O data/SRR1039508_1.fastq.gz

# 4. Run
nextflow run main.nf -profile docker \
    --genome_index genome/grch38/genome \
    --gtf genome/gencode.v38.annotation.gtf
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--samplesheet` | `assets/samplesheet.csv` | CSV: sample_id, fastq_1, fastq_2, condition |
| `--genome_index` | required | HISAT2 index prefix |
| `--gtf` | required | Gene annotation GTF |
| `--outdir` | `results` | Output directory |
| `--strandedness` | `2` (reverse) | featureCounts strandedness (0/1/2) |
| `--ref_condition` | `untreated` | DESeq2 reference level |

## Output

```
results/
├── fastqc_raw/       Raw read QC reports
├── fastp/            Trimming reports (JSON)
├── hisat2/           Alignment logs
├── bam/              Sorted BAM files
├── counts/           Gene count matrix
├── deseq2/           DE results, volcano plot, PCA plot
└── multiqc/          Aggregated QC report
```

## Design Decisions

- **HISAT2 over STAR** — runs on 8GB RAM. STAR requires 32GB for the human genome. Accessible on any machine.
- **BioContainers** — published, maintained Docker containers. No custom builds.
- **Configurable reference level** — `--ref_condition` sets the DESeq2 baseline. Works with any experimental design.
- **Adaptive gene filter** — automatically adjusts minimum count threshold based on library size (stringent for real data, permissive for test data).
- **Test profile** — synthetic 50-gene genome with genome-sampled reads. Verifies the full pipeline in ~2 minutes.

## Licence

MIT
