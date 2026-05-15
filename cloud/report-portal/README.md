# RNA-seq Report Portal

Small FastAPI service for registering cloud RNA-seq runs and serving signed links to reports stored in S3.

This is intentionally separate from the Nextflow pipeline. The pipeline remains responsible for compute and published artefacts; the portal gives reviewers and collaborators a minimal cloud-facing surface for run status and report access.

## Architecture

```text
Nextflow AWS Batch run
  -> s3://bucket/results/<run>/
      -> pipeline_info/report.html
      -> pipeline_info/timeline.html
      -> pipeline_info/trace.txt
      -> pipeline_info/dag.dot
      -> multiqc/multiqc_report.html
  -> POST /runs into this service
  -> collaborators request signed report URLs from /runs/{id}/artifacts/{artifact}/presign
```

## Configuration

| Variable | Example | Purpose |
| --- | --- | --- |
| `DATABASE_URL` | `postgresql+psycopg://rnaseq:change_me@db:5432/rnaseq` | Metadata database. Defaults to local SQLite for development. |
| `AWS_REGION` | `eu-west-2` | Region used by the AWS SDK. |
| AWS credentials | IAM role, env vars, or workload identity | Required only for signed S3 URLs. |

## Local Smoke Run

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
uvicorn app.main:app --reload --port 8000
```

Register a completed run:

```bash
curl -X POST http://localhost:8000/runs \
  -H 'content-type: application/json' \
  -d '{
    "run_id": "airway-test-001",
    "name": "Synthetic CI airway test",
    "status": "succeeded",
    "s3_prefix": "s3://my-rnaseq-bucket/results/airway-test-001"
  }'
```

Get a signed report URL:

```bash
curl http://localhost:8000/runs/airway-test-001/artifacts/report/presign
```

## Container

```bash
docker build -t rnaseq-report-portal .
docker run --rm -p 8000:8000 \
  -e DATABASE_URL=postgresql+psycopg://rnaseq:change_me@host.docker.internal:5432/rnaseq \
  rnaseq-report-portal
```

## Tests

```bash
pip install -r requirements.txt
pytest tests
```
