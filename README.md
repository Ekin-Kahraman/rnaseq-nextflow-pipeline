# RNA-seq Nextflow Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Nextflow](https://img.shields.io/badge/Nextflow-%E2%89%A524.0-brightgreen)](https://www.nextflow.io/)

End-to-end bulk RNA-seq pipeline in Nextflow DSL2: raw FASTQ reads through quality control, alignment, gene quantification, and differential expression. Every step runs in its own Docker container.

Designed for the GSE152075 SARS-CoV-2 nasopharyngeal dataset (6 samples: 3 COVID-positive, 3 negative). For the full covariate-adjusted statistical analysis on 484 samples, see [bulk-rnaseq-differential-expression](https://github.com/Ekin-Kahraman/bulk-rnaseq-differential-expression).

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
 featureCounts ── Gene-level quantification
    │
    ▼
 DESeq2 ──────── Differential expression (positive vs negative)
    │
    ▼
 MultiQC ─────── Aggregate QC report
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

All containers sourced from [BioContainers](https://biocontainers.pro/). No custom Dockerfiles.

## Quick Start

**Prerequisites:** [Nextflow](https://www.nextflow.io/) (>=24.0), [Docker](https://www.docker.com/), Java (>=11)

### Test data (synthetic, runs in ~2 minutes)

```bash
git clone https://github.com/Ekin-Kahraman/rnaseq-nextflow-pipeline.git
cd rnaseq-nextflow-pipeline
python test/create_test_data.py
nextflow run main.nf -profile test,docker \
    --genome_index "$(pwd)/test/genome" \
    --gtf "$(pwd)/test/genes.gtf"
```

### Real data (GSE152075)

```bash
# Download HISAT2 GRCh38 index (~4GB)
# Download FASTQ files from SRA (see assets/samplesheet.csv for accessions)

nextflow run main.nf -profile docker \
    --genome_index /path/to/grch38/genome \
    --gtf /path/to/gencode.v38.annotation.gtf \
    --samplesheet assets/samplesheet.csv
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--samplesheet` | `assets/samplesheet.csv` | CSV with columns: sample_id, fastq_1, fastq_2, condition |
| `--genome_index` | required | HISAT2 index prefix |
| `--gtf` | required | Gene annotation GTF |
| `--outdir` | `results` | Output directory |
| `--strandedness` | `2` (reverse) | featureCounts strandedness (0/1/2) |

## Output

```
results/
├── fastqc_raw/       Raw read QC reports (HTML + ZIP)
├── fastp/            Trimming reports (JSON)
├── hisat2/           Alignment logs
├── bam/              Sorted BAM files
├── counts/           Gene count matrix (featureCounts)
├── deseq2/           DE results CSV, volcano plot, PCA plot
└── multiqc/          Aggregated QC report
```

## Project Structure

```
rnaseq-nextflow-pipeline/
├── main.nf              Pipeline (7 processes, Nextflow DSL2)
├── nextflow.config      Parameters, containers, profiles
├── assets/
│   └── samplesheet.csv  Sample metadata (SRA accessions)
├── test/
│   ├── create_test_data.py   Generate synthetic test data
│   ├── samplesheet.csv       Test sample metadata
│   ├── genome.fa             Synthetic genome (50 genes)
│   └── genes.gtf             Synthetic annotation
├── LICENSE
└── README.md
```

## Design Decisions

- **HISAT2 over STAR** — runs on 8GB RAM. STAR requires 32GB for the human genome index. Anyone can clone and run this pipeline.
- **BioContainers, not custom Dockerfiles** — industry standard, maintained by the community, reproducible without building.
- **Separate samtools process** — HISAT2 and samtools in their own containers. Clean separation of concerns.
- **Test profile** — synthetic 50-gene genome with reads sampled from the reference. Runs in ~2 minutes. Verifies the pipeline without downloading 30GB of real data.
- **DESeq2 dispersion fallback** — handles small test datasets where standard dispersion fitting fails. Uses gene-wise estimates when the mean-dispersion trend cannot be fitted.
- **Configurable strandedness** — `--strandedness` parameter for featureCounts. Default reverse-stranded (standard for Illumina dUTP protocols), unstranded for test data.

## Licence

MIT
