<h1 align="center">RNA-seq Nextflow pipeline</h1>

<p align="center">
  <strong>Take paired-end FASTQ reads through gene counting, differential expression and quality reporting in one containerised workflow.</strong>
</p>

<p align="center">
  Designed around the airway smooth-muscle data from
  <a href="https://doi.org/10.1371/journal.pone.0099625">Himes et al. (2014)</a>
  (<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778">GSE52778</a>).
</p>

<p align="center">
  <a href="https://github.com/Ekin-Kahraman/rnaseq-nextflow-pipeline/actions/workflows/ci.yml"><img alt="CI" src="https://github.com/Ekin-Kahraman/rnaseq-nextflow-pipeline/actions/workflows/ci.yml/badge.svg"></a>
  <a href="LICENSE"><img alt="MIT licence" src="https://img.shields.io/badge/License-MIT-blue.svg"></a>
  <a href="https://www.nextflow.io/"><img alt="Nextflow 24.0 or later" src="https://img.shields.io/badge/Nextflow-%E2%89%A524.0-brightgreen"></a>
  <a href="docs/cloud.md"><img alt="AWS Batch profile" src="https://img.shields.io/badge/AWS%20Batch-profile-orange"></a>
</p>

<p align="center">
  <a href="#workflow">Workflow</a> · <a href="#output">Output</a> · <a href="#quick-start">Quick start</a> · <a href="#reproducibility">Reproducibility</a> · <a href="#details">Details</a> · <a href="#limitations">Limitations</a>
</p>

## What this project adds

<table>
  <tr>
    <td valign="top"><strong>Input → output</strong></td>
    <td>Paired FASTQ reads, a sample sheet, HISAT2 index and GTF annotation become QC reports, trimmed reads, sorted BAMs, gene counts, DESeq2 results and plots.</td>
  </tr>
  <tr>
    <td valign="top"><strong>Preflight checks</strong></td>
    <td>Checks required columns, unique clean sample IDs, exactly two conditions with the reference present, paired reads and a complete eight-file HISAT2 index before analysis starts.</td>
  </tr>
  <tr>
    <td valign="top"><strong>Failure and restart handling</strong></td>
    <td>Normal DESeq2 errors stop and report the original error message; report files are overwritten safely on reruns, and Nextflow <code>-resume</code> can reuse cached tasks.</td>
  </tr>
  <tr>
    <td valign="top"><strong>Reproducibility</strong></td>
    <td>Versioned BioContainers and saved trace, timeline, report and DAG; CI runs the complete synthetic workflow and validates output sample and gene IDs.</td>
  </tr>
  <tr>
    <td valign="top"><strong>Optional report portal</strong></td>
    <td>A separate FastAPI service records run metadata and serves links to reports in S3; it does not run the analysis. No real AWS Batch run is claimed.</td>
  </tr>
</table>

## Workflow

`FASTQ → FastQC → fastp → HISAT2 → samtools → featureCounts → DESeq2 → MultiQC`

## Output

| Directory | Useful outputs |
| --- | --- |
| `results/fastqc_raw/`, `fastp/`, `multiqc/` | Read-level QC and aggregate report |
| `results/bam/`, `counts/` | Sorted alignments and gene count matrix |
| `results/deseq2/` | Differential-expression CSV, volcano plot and PCA plot |
| `results/pipeline_info/` | Nextflow report, timeline, trace and workflow DAG |

## Quick start

**Prerequisites:** [Nextflow](https://www.nextflow.io/) 24.0 or later, [Docker](https://www.docker.com/) and Java 21 (the CI runtime).

Generate synthetic reads and run the full test profile:

```bash
git clone https://github.com/Ekin-Kahraman/rnaseq-nextflow-pipeline.git
cd rnaseq-nextflow-pipeline
python test/create_test_data.py
nextflow run main.nf -profile test,docker \
    --genome_index "$(pwd)/test/genome" \
    --gtf "$(pwd)/test/genes.gtf"
python scripts/validate_outputs.py results --samplesheet test/samplesheet.csv
```

The synthetic sample sheet pairs two donors across two conditions and exercises the donor-aware design. The test profile permits a logged gene-wise dispersion fallback only if standard DESeq2 fitting fails on these synthetic reads; normal runs do not use that fallback.

Add `-resume` to reuse completed tasks after a change, while retaining the same Nextflow work directory.

## Reproducibility

All analysis steps use version-pinned BioContainers. Each run writes Nextflow execution records to `results/pipeline_info/`. CI tests input validation, the full synthetic workflow and output checks. The output validator checks sample alignment, duplicate IDs, count values and DESeq2 result fields.

The optional report portal has its own tests. The portal and AWS Batch profile are not evidence of a real cloud analysis; the synthetic run is the verified execution path.

## Details

<details>
<summary>Real data: Himes airway dataset</summary>

The workflow is designed around the [Himes et al. (2014)](https://doi.org/10.1371/journal.pone.0099625) airway smooth-muscle dataset (dexamethasone vs untreated, GEO [GSE52778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778)). This dataset is used in the [DESeq2 vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) and the [Bioconductor RNA-seq workflow](https://www.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html). For a covariate-adjusted analysis of a COVID-19 cohort, see [bulk-rnaseq-differential-expression](https://github.com/Ekin-Kahraman/bulk-rnaseq-differential-expression).


Himes et al. report using the [TruSeq RNA Sample Prep Kit v2](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0099625), which Illumina describes as non-stranded ([kit description](https://www.illumina.com/products/by-type/sequencing-kits/library-prep-kits/truseq-rna-v2.html), [technical note](https://knowledge.illumina.com/library-preparation/rna-library-prep/library-preparation-rna-library-prep-reference_material-list/000002238)). The real-data command therefore sets `--strandedness 0` explicitly; the general pipeline default remains `2` for reverse-stranded libraries.

Download the reference files and FASTQs, then run the workflow:

```bash
# HISAT2 GRCh38 index (~4GB)
mkdir -p genome && cd genome
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
tar xzf grch38_genome.tar.gz

# Gencode v38 GTF
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
gunzip gencode.v38.annotation.gtf.gz
cd ..

# Download the paired FASTQs listed in assets/samplesheet.csv
mkdir -p data
wget -P data https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/008/SRR1039508/SRR1039508_1.fastq.gz
wget -P data https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/008/SRR1039508/SRR1039508_2.fastq.gz
wget -P data https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/009/SRR1039509/SRR1039509_1.fastq.gz
wget -P data https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/009/SRR1039509/SRR1039509_2.fastq.gz
wget -P data https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/002/SRR1039512/SRR1039512_1.fastq.gz
wget -P data https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/002/SRR1039512/SRR1039512_2.fastq.gz
wget -P data https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/003/SRR1039513/SRR1039513_1.fastq.gz
wget -P data https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/003/SRR1039513/SRR1039513_2.fastq.gz

nextflow run main.nf -profile docker \
    --samplesheet assets/samplesheet.csv \
    --genome_index genome/grch38/genome \
    --gtf genome/gencode.v38.annotation.gtf \
    --strandedness 0 \
    --design donor_condition
```

The bundled sample sheet contains two donors, each with one untreated and one dexamethasone sample. Its paired design is explicit in the command. The original Himes dataset includes more samples; this four-sample subset is for illustrating the pipeline, not a powered biological comparison.

</details>

<details>
<summary>Samples and tool versions</summary>

The bundled sample sheet uses two paired donors from the Himes airway dataset:

| Sample | SRA | Condition | Donor |
| --- | --- | --- | --- |
| N61311_untreated | SRR1039508 | untreated | N61311 |
| N61311_Dex | SRR1039509 | dexamethasone | N61311 |
| N052611_untreated | SRR1039512 | untreated | N052611 |
| N052611_Dex | SRR1039513 | dexamethasone | N052611 |

Each pipeline stage runs in a pinned BioContainer from [BioContainers](https://biocontainers.pro/):

| Process | Tool and version | Container |
| --- | --- | --- |
| `FASTQC_RAW` | FastQC 0.12.1 | `quay.io/biocontainers/fastqc` |
| `FASTP` | fastp 0.23.4 | `quay.io/biocontainers/fastp` |
| `HISAT2_ALIGN` | HISAT2 2.2.1 | `quay.io/biocontainers/hisat2` |
| `SAMTOOLS_SORT` | samtools 1.21 | `quay.io/biocontainers/samtools` |
| `FEATURECOUNTS` | Subread 2.0.6 | `quay.io/biocontainers/subread` |
| `DESEQ2` | DESeq2 1.42 and ggplot2 | `quay.io/biocontainers/bioconductor-deseq2` |
| `MULTIQC` | MultiQC 1.27 | `quay.io/biocontainers/multiqc` |

</details>

<details>
<summary>Sample sheet, parameters and design</summary>

The CSV must contain `sample_id`, `fastq_1`, `fastq_2` and `condition`; `donor` is needed for `--design donor_condition`. IDs and condition labels must be free of surrounding whitespace. Sample IDs may contain letters, numbers, dots, underscores and hyphens. The workflow requires exactly two condition labels and a matching `--ref_condition`.

| Parameter | Default | Description |
| --- | --- | --- |
| `--samplesheet` | `assets/samplesheet.csv` | CSV sample manifest |
| `--genome_index` | required | Prefix for all eight `.ht2` or `.ht2l` index files |
| `--gtf` | required | Gene annotation GTF |
| `--outdir` | `results` | Output directory |
| `--strandedness` | `2` (reverse) | featureCounts strandedness: 0 unstranded, 1 forward, 2 reverse |
| `--ref_condition` | `untreated` | DESeq2 reference level for the contrast |
| `--design` | `condition` | `condition` or `donor_condition` (requires paired donor labels) |
| `--aws_queue` | none | AWS Batch queue for `-profile awsbatch` |
| `--aws_region` | `eu-west-2` | AWS region for `-profile awsbatch` |
| `--aws_workdir` | none | S3 work directory for `-profile awsbatch` |

`condition` uses a `~ condition` DESeq2 model and requires at least one residual degree of freedom. If either group has fewer than two samples, validation warns that dispersion estimates may be unreliable. `donor_condition` uses `~ donor + condition`; validation requires at least two donors, both conditions represented for every donor, and at least one residual degree of freedom. The reference level is controlled by `--ref_condition`. Choose the design explicitly for the study: the condition-only default does not adjust for donor.

Other implementation choices:

- **HISAT2** is the splice-aware aligner used by this CPU-based workflow; this repository has no controlled accuracy or memory comparison with STAR.
- **featureCounts** uses `gene_name`, paired-read counting and configurable strandedness. Counts depend on the annotation and settings; equivalence with other tools is not assumed.
- **fastp** performs adapter detection and quality filtering with Phred 20 and minimum read length 36. These are this pipeline's settings, not a claim that it reproduces Himes et al.'s trimming procedure.
- **Docker** is configured for local runs; **Singularity** is available for HPC environments where Docker is typically unavailable.
- Nextflow execution reports are overwritten on reruns so a previous run does not block new run records.

</details>

<details>
<summary>AWS Batch and Seqera Platform</summary>

See [docs/cloud.md](docs/cloud.md) for AWS Batch and Seqera Platform launch notes. The AWS Batch profile configures managed compute with S3 work and output paths; no real AWS run record is committed.

```bash
nextflow run Ekin-Kahraman/rnaseq-nextflow-pipeline \
    -profile awsbatch \
    --aws_queue rnaseq-job-queue \
    --aws_region eu-west-2 \
    --aws_workdir s3://my-rnaseq-bucket/work \
    --samplesheet s3://my-rnaseq-bucket/inputs/samplesheet.csv \
    --genome_index s3://my-rnaseq-bucket/reference/grch38/genome \
    --gtf s3://my-rnaseq-bucket/reference/gencode.v38.annotation.gtf \
    --outdir s3://my-rnaseq-bucket/results/airway
```

</details>

<details>
<summary>Optional report portal</summary>

The [report portal](cloud/report-portal/) records run details and gives collaborators temporary links to reports stored in S3. It uses FastAPI, with Postgres in the deployment configuration and SQLite for local tests. The dashboard and API do not run the RNA-seq analysis.

```bash
cd cloud/report-portal
pip install -r requirements.txt
uvicorn app.main:app --reload --port 8000
```

Run the local Postgres stack with `docker compose up --build` from `cloud/report-portal/`. `render.yaml` defines the Docker FastAPI service, managed Postgres, demo seed run and AWS secret placeholders.

Demo endpoints (availability may vary; the health check timed out during the 7 September 2026 documentation review):

- Dashboard: <https://rnaseq-report-portal.onrender.com/>
- Health: <https://rnaseq-report-portal.onrender.com/health>
- Seeded artefact metadata: <https://rnaseq-report-portal.onrender.com/runs/synthetic-ci-001/artifacts/report>

</details>

## Limitations

- The bundled airway sample sheet contains two donors with one sample per condition, giving one residual degree of freedom in the donor-aware model. It illustrates the workflow and does not support a well-powered biological conclusion.
- The default design is `~ condition`; use `--design donor_condition` for a paired study with at least two donors and both conditions for each donor. The condition-only model ignores donor blocking.
- ENA lists an additional unpaired FASTQ file for SRR1039513; this paired-end workflow consumes only its `_1`/`_2` files. The omitted singleton reads are not part of the analysis.
- CI uses deterministic synthetic reads and a synthetic 50-gene genome. Synthetic DESeq2 fitting may use the explicitly enabled gene-wise dispersion fallback after a logged failure; normal runs stop with the original fitting error.
- Real Himes/GSE52778 runs require external FASTQs, a GRCh38 HISAT2 index and Gencode annotation.
- No STAR option is implemented, and no controlled aligner or quantification benchmark is included.
- The AWS Batch profile and report portal are implemented, but no real AWS Batch analysis is claimed. A seeded portal demo is not evidence of a completed cloud run.

## Licence

MIT
