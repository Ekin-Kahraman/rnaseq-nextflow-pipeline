from __future__ import annotations

import os
import sys
from pathlib import Path

from fastapi.testclient import TestClient

PORTAL_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PORTAL_ROOT))
os.environ.setdefault("DATABASE_URL", "sqlite:///:memory:")

from app.main import create_app, get_s3_client  # noqa: E402


def test_register_list_and_update_run(tmp_path):
    app = create_app(f"sqlite:///{tmp_path / 'runs.db'}")
    client = TestClient(app)

    assert client.get("/health").json() == {"status": "ok"}

    create_resp = client.post(
        "/runs",
        json={
            "run_id": "ci-run-001",
            "name": "Synthetic CI test",
            "status": "running",
            "s3_prefix": "s3://rnaseq-results/ci-run-001",
        },
    )
    assert create_resp.status_code == 201
    assert create_resp.json()["id"] == "ci-run-001"

    duplicate_resp = client.post(
        "/runs",
        json={"run_id": "ci-run-001", "name": "Duplicate"},
    )
    assert duplicate_resp.status_code == 409

    update_resp = client.patch("/runs/ci-run-001", json={"status": "succeeded"})
    assert update_resp.status_code == 200
    assert update_resp.json()["status"] == "succeeded"

    list_resp = client.get("/runs?status=succeeded")
    assert list_resp.status_code == 200
    assert [run["id"] for run in list_resp.json()] == ["ci-run-001"]


def test_presign_report_uses_expected_s3_key(tmp_path):
    app = create_app(f"sqlite:///{tmp_path / 'runs.db'}")

    class FakeS3:
        def generate_presigned_url(self, operation, Params, ExpiresIn):
            assert operation == "get_object"
            assert Params == {
                "Bucket": "rnaseq-results",
                "Key": "airway-001/pipeline_info/report.html",
            }
            assert ExpiresIn == 900
            return "https://signed.example/airway-001/report.html"

    app.dependency_overrides[get_s3_client] = lambda: FakeS3()
    client = TestClient(app)

    client.post(
        "/runs",
        json={
            "run_id": "airway-001",
            "name": "Airway AWS Batch run",
            "status": "succeeded",
            "s3_prefix": "s3://rnaseq-results/airway-001",
        },
    )

    resp = client.get("/runs/airway-001/artifacts/report/presign?expires=900")
    assert resp.status_code == 200
    body = resp.json()
    assert body["s3_uri"] == "s3://rnaseq-results/airway-001/pipeline_info/report.html"
    assert body["url"] == "https://signed.example/airway-001/report.html"
