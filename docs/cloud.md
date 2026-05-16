# Cloud Execution

This pipeline can run locally with Docker or on cloud compute through Nextflow's AWS Batch executor. Keep FASTQs, references, work directories and results in the same region to avoid slow cross-region transfer and avoidable egress cost.

## AWS Batch

Prerequisites:

- An AWS Batch compute environment and job queue.
- An S3 bucket for the Nextflow work directory and published results.
- IAM permissions for Batch, ECS, ECR image pulls, CloudWatch logs and S3 read/write.
- The AWS CLI configured locally or a Seqera Platform workspace connected to AWS.

Example command:

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

The `awsbatch` profile sets `process.executor = 'awsbatch'`, enables Docker containers and writes the Nextflow work directory to `--aws_workdir`.

## Seqera Platform

Seqera Platform can launch this repository directly from GitHub. Use the `awsbatch` profile for AWS execution, or use the `docker` profile for a local/VM compute environment. The `nextflow_schema.json` file exposes the main parameters in the launch form.

Suggested launch fields:

| Field | Value |
| --- | --- |
| Repository | `https://github.com/Ekin-Kahraman/rnaseq-nextflow-pipeline` |
| Revision | `main` or a release tag |
| Profile | `awsbatch` |
| Work directory | S3 URI matching `--aws_workdir` |
| Output directory | S3 URI passed as `--outdir` |

## Run Evidence

Every run writes reproducibility artefacts under `results/pipeline_info/`:

- `report.html` - Nextflow execution report.
- `timeline.html` - task-level runtime timeline.
- `trace.txt` - machine-readable task trace.
- `dag.dot` - workflow graph.

The CI workflow also runs `scripts/validate_outputs.py` against the synthetic test run so broken or incomplete published artefacts fail the build.

## Report Portal

The repository includes a minimal FastAPI report portal in `cloud/report-portal/`. It stores run metadata in Postgres through `DATABASE_URL` and returns presigned S3 URLs for the main Nextflow artefacts:

- `pipeline_info/report.html`
- `pipeline_info/timeline.html`
- `pipeline_info/trace.txt`
- `pipeline_info/dag.dot`
- `multiqc/multiqc_report.html`

The root route renders a browser dashboard for registered runs, while `/docs` exposes the OpenAPI UI. This makes the cloud proof reviewable without needing to read raw JSON.

Local smoke run:

```bash
cd cloud/report-portal
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
uvicorn app.main:app --reload --port 8000
```

Local Postgres stack:

```bash
cd cloud/report-portal
docker compose up --build
```

Production shape:

```bash
docker build -t rnaseq-report-portal cloud/report-portal
docker run --rm -p 8000:8000 \
  -e DATABASE_URL=postgresql+psycopg://rnaseq:change_me@postgres:5432/rnaseq \
  -e AWS_REGION=eu-west-2 \
  rnaseq-report-portal
```

For ECS/Fargate, give the task role read-only access to the S3 result prefix and keep write access limited to the Postgres database. The service does not need permission to launch Batch jobs.

For Render, the repository root `render.yaml` provisions the Docker web service and a managed Postgres database. After pushing `render.yaml`, open:

```text
https://dashboard.render.com/blueprint/new?repo=https://github.com/Ekin-Kahraman/rnaseq-nextflow-pipeline
```

Fill the AWS secret fields only when presigned S3 URLs need to work against a real result bucket. The seeded demo run is enough to show the portal and metadata layer before wiring real S3 credentials.

Current live smoke deployment:

- Dashboard: <https://rnaseq-report-portal.onrender.com/>
- Health: <https://rnaseq-report-portal.onrender.com/health>
- Seeded artefact metadata: <https://rnaseq-report-portal.onrender.com/runs/synthetic-ci-001/artifacts/report>

This free Render service was created directly from the Docker portal root for fast public proof. The validated `render.yaml` remains the reproducible Postgres-backed deployment path.
