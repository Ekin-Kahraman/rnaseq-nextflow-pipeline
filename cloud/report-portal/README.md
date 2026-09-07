# RNA-seq Report Portal

Small FastAPI service for registering cloud RNA-seq runs and serving signed links to reports stored in S3.

The Nextflow pipeline runs the analysis. This service stores run details and
provides a dashboard and report access. The seeded demonstration is not a
record of a completed AWS Batch analysis.

The root route renders a small dashboard for browser review. The API remains available through `/docs`, `/runs` and `/runs/{id}/artifacts`.

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
| `DEMO_RUN_ID` | `synthetic-ci-001` | Optional seed run for live demos. |
| `DEMO_S3_PREFIX` | `s3://bucket/results/synthetic-ci-001` | Optional S3 prefix for the seeded demo run. |

## Local Smoke Run

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
uvicorn app.main:app --reload --port 8000
```

Open <http://localhost:8000> for the dashboard or <http://localhost:8000/docs> for the OpenAPI UI.

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

Run the full local stack with Postgres:

```bash
docker compose up --build
```

The compose stack seeds a synthetic run so the dashboard has a visible record immediately.

## Render Blueprint

The repository root contains `render.yaml` for a reproducible Render deployment:

- Docker web service built from `cloud/report-portal/Dockerfile`.
- Managed Postgres database connected through `DATABASE_URL`.
- Optional demo run seeded through `DEMO_RUN_ID` and `DEMO_S3_PREFIX`.
- AWS credentials stored as Dashboard secrets for S3 presigned URLs.

Open the Blueprint after the file is pushed:

```text
https://dashboard.render.com/blueprint/new?repo=https://github.com/Ekin-Kahraman/rnaseq-nextflow-pipeline
```

Demo endpoints (availability may vary; the health check timed out during the
7 September 2026 documentation review):

- Dashboard: <https://rnaseq-report-portal.onrender.com/>
- Health: <https://rnaseq-report-portal.onrender.com/health>
- Seeded artefact metadata: <https://rnaseq-report-portal.onrender.com/runs/synthetic-ci-001/artifacts/report>

## Tests

```bash
pip install -r requirements.txt
pytest tests
```
