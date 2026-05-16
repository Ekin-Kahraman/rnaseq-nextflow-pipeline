from __future__ import annotations

import os
from datetime import datetime, timezone
from html import escape
from typing import Literal
from urllib.parse import urlparse
from uuid import uuid4

from fastapi import Depends, FastAPI, HTTPException, Query, status
from fastapi.responses import HTMLResponse
from pydantic import BaseModel, Field
from sqlalchemy import Column, DateTime, String, create_engine, select, text
from sqlalchemy.orm import Session, declarative_base, sessionmaker


DEFAULT_DATABASE_URL = "sqlite:///./runs.db"
RunStatus = Literal["submitted", "running", "succeeded", "failed", "unknown"]
ArtifactName = Literal["report", "timeline", "trace", "dag", "multiqc"]
ARTIFACTS: tuple[ArtifactName, ...] = ("report", "timeline", "trace", "dag", "multiqc")
VALID_STATUSES = {"submitted", "running", "succeeded", "failed", "unknown"}

Base = declarative_base()


def utc_now() -> datetime:
    return datetime.now(timezone.utc)


class RunModel(Base):
    __tablename__ = "runs"

    id = Column(String, primary_key=True)
    name = Column(String, nullable=False)
    status = Column(String, nullable=False, default="submitted")
    s3_prefix = Column(String, nullable=True)
    nextflow_report = Column(String, nullable=False, default="pipeline_info/report.html")
    timeline = Column(String, nullable=False, default="pipeline_info/timeline.html")
    trace = Column(String, nullable=False, default="pipeline_info/trace.txt")
    dag = Column(String, nullable=False, default="pipeline_info/dag.dot")
    multiqc_report = Column(String, nullable=False, default="multiqc/multiqc_report.html")
    created_at = Column(DateTime(timezone=True), nullable=False, default=utc_now)
    updated_at = Column(DateTime(timezone=True), nullable=False, default=utc_now, onupdate=utc_now)


class RunCreate(BaseModel):
    run_id: str | None = Field(default=None, description="Stable external run id. Defaults to a UUID.")
    name: str = Field(min_length=1)
    status: RunStatus = "submitted"
    s3_prefix: str | None = Field(default=None, description="S3 prefix containing published Nextflow results.")
    nextflow_report: str = "pipeline_info/report.html"
    timeline: str = "pipeline_info/timeline.html"
    trace: str = "pipeline_info/trace.txt"
    dag: str = "pipeline_info/dag.dot"
    multiqc_report: str = "multiqc/multiqc_report.html"


class RunUpdate(BaseModel):
    name: str | None = Field(default=None, min_length=1)
    status: RunStatus | None = None
    s3_prefix: str | None = None
    nextflow_report: str | None = None
    timeline: str | None = None
    trace: str | None = None
    dag: str | None = None
    multiqc_report: str | None = None


class RunOut(BaseModel):
    id: str
    name: str
    status: str
    s3_prefix: str | None
    nextflow_report: str
    timeline: str
    trace: str
    dag: str
    multiqc_report: str
    created_at: datetime
    updated_at: datetime

    model_config = {"from_attributes": True}


class ArtifactOut(BaseModel):
    artifact: ArtifactName
    s3_uri: str | None
    presign_path: str


def model_data(model: BaseModel, *, exclude_unset: bool = False) -> dict:
    if hasattr(model, "model_dump"):
        return model.model_dump(exclude_unset=exclude_unset)
    return model.dict(exclude_unset=exclude_unset)


def normalise_database_url(database_url: str) -> str:
    if database_url.startswith("postgres://"):
        return "postgresql+psycopg://" + database_url.removeprefix("postgres://")
    if database_url.startswith("postgresql://"):
        return "postgresql+psycopg://" + database_url.removeprefix("postgresql://")
    return database_url


def engine_options(database_url: str) -> dict:
    if database_url.startswith("sqlite"):
        return {"connect_args": {"check_same_thread": False}}
    return {"pool_pre_ping": True}


def parse_s3_uri(uri: str) -> tuple[str, str]:
    parsed = urlparse(uri)
    if parsed.scheme != "s3" or not parsed.netloc or not parsed.path.strip("/"):
        raise HTTPException(status_code=400, detail=f"Invalid S3 URI: {uri}")
    return parsed.netloc, parsed.path.lstrip("/")


def artifact_path(run: RunModel, artifact: ArtifactName) -> str:
    attr = "nextflow_report" if artifact == "report" else "multiqc_report" if artifact == "multiqc" else artifact
    return getattr(run, attr)


def artifact_s3_uri(run: RunModel, artifact: ArtifactName) -> str | None:
    value = artifact_path(run, artifact)
    if value.startswith("s3://"):
        return value
    if not run.s3_prefix:
        return None
    return f"{run.s3_prefix.rstrip('/')}/{value.lstrip('/')}"


def resolve_artifact_uri(run: RunModel, artifact: ArtifactName) -> str:
    uri = artifact_s3_uri(run, artifact)
    if uri is None:
        raise HTTPException(status_code=404, detail="Run has no s3_prefix for relative artefact paths")
    return uri


def render_dashboard(runs: list[RunModel]) -> str:
    cards = []
    for run in runs:
        artifact_links = "".join(
            f'<li><a href="/runs/{escape(run.id)}/artifacts/{artifact}">{artifact}</a></li>'
            for artifact in ARTIFACTS
        )
        cards.append(
            f"""
            <article class="run">
              <div>
                <h2>{escape(run.name)}</h2>
                <p class="meta">{escape(run.id)} | {escape(run.status)}</p>
                <p>{escape(run.s3_prefix or "No S3 prefix configured")}</p>
              </div>
              <ul>{artifact_links}</ul>
            </article>
            """
        )
    body = "\n".join(cards) if cards else '<p class="empty">No runs registered yet.</p>'
    return f"""
    <!doctype html>
    <html lang="en">
      <head>
        <meta charset="utf-8">
        <meta name="viewport" content="width=device-width, initial-scale=1">
        <title>RNA-seq Report Portal</title>
        <style>
          :root {{
            color-scheme: light;
            --bg: #f7f8fb;
            --panel: #ffffff;
            --text: #111827;
            --muted: #5b6472;
            --line: #d8dee8;
            --accent: #0f766e;
          }}
          * {{ box-sizing: border-box; }}
          body {{
            margin: 0;
            background: var(--bg);
            color: var(--text);
            font-family: ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
          }}
          main {{ max-width: 1040px; margin: 0 auto; padding: 40px 20px; }}
          header {{ display: flex; justify-content: space-between; gap: 20px; align-items: flex-start; margin-bottom: 28px; }}
          h1 {{ margin: 0 0 8px; font-size: clamp(2rem, 5vw, 3.5rem); line-height: 1; }}
          h2 {{ margin: 0 0 6px; font-size: 1.05rem; }}
          p {{ margin: 0; color: var(--muted); overflow-wrap: anywhere; }}
          a {{ color: var(--accent); font-weight: 700; text-decoration: none; }}
          a:hover {{ text-decoration: underline; }}
          .actions {{ display: flex; gap: 12px; flex-wrap: wrap; justify-content: flex-end; }}
          .button {{ border: 1px solid var(--line); border-radius: 6px; padding: 9px 12px; background: var(--panel); color: var(--text); }}
          .run {{
            display: grid;
            grid-template-columns: minmax(0, 1fr) auto;
            gap: 24px;
            align-items: center;
            padding: 18px;
            background: var(--panel);
            border: 1px solid var(--line);
            border-radius: 8px;
            margin-bottom: 12px;
          }}
          .meta {{ font-size: 0.9rem; margin-bottom: 8px; }}
          ul {{ display: flex; gap: 10px; flex-wrap: wrap; justify-content: flex-end; padding: 0; margin: 0; list-style: none; }}
          li a {{ display: block; border: 1px solid var(--line); border-radius: 6px; padding: 7px 10px; background: #f9fafb; }}
          .empty {{ padding: 18px; background: var(--panel); border: 1px solid var(--line); border-radius: 8px; }}
          @media (max-width: 760px) {{
            header, .run {{ display: block; }}
            .actions, ul {{ justify-content: flex-start; margin-top: 14px; }}
          }}
        </style>
      </head>
      <body>
        <main>
          <header>
            <div>
              <h1>RNA-seq Report Portal</h1>
              <p>Registered Nextflow runs with S3-backed report artefacts.</p>
            </div>
            <nav class="actions">
              <a class="button" href="/docs">API docs</a>
              <a class="button" href="/health">Health</a>
            </nav>
          </header>
          <section>{body}</section>
        </main>
      </body>
    </html>
    """


def seed_demo_run(session_local: sessionmaker[Session]) -> None:
    run_id = os.getenv("DEMO_RUN_ID")
    if not run_id:
        return
    status_value = os.getenv("DEMO_RUN_STATUS", "succeeded")
    run_status = status_value if status_value in VALID_STATUSES else "unknown"
    with session_local() as db:
        if db.get(RunModel, run_id):
            return
        db.add(
            RunModel(
                id=run_id,
                name=os.getenv("DEMO_RUN_NAME", "Synthetic CI airway test"),
                status=run_status,
                s3_prefix=os.getenv("DEMO_S3_PREFIX"),
            )
        )
        db.commit()


def get_s3_client():
    import boto3

    return boto3.client("s3", region_name=os.getenv("AWS_REGION"))


def create_app(database_url: str | None = None) -> FastAPI:
    db_url = normalise_database_url(database_url or os.getenv("DATABASE_URL", DEFAULT_DATABASE_URL))
    engine = create_engine(db_url, **engine_options(db_url))
    session_local = sessionmaker(bind=engine, autoflush=False, autocommit=False)
    Base.metadata.create_all(bind=engine)
    seed_demo_run(session_local)

    app = FastAPI(
        title="RNA-seq Report Portal",
        version="0.1.0",
        description="Registers Nextflow RNA-seq runs and signs S3 report artefact URLs.",
    )

    def get_db():
        db = session_local()
        try:
            yield db
        finally:
            db.close()

    @app.get("/", response_class=HTMLResponse, include_in_schema=False)
    def dashboard(db: Session = Depends(get_db)) -> HTMLResponse:
        stmt = select(RunModel).order_by(RunModel.created_at.desc()).limit(50)
        return HTMLResponse(render_dashboard(list(db.scalars(stmt))))

    @app.get("/health")
    def health(db: Session = Depends(get_db)) -> dict:
        db.execute(text("SELECT 1"))
        return {"status": "ok"}

    @app.post("/runs", response_model=RunOut, status_code=status.HTTP_201_CREATED)
    def create_run(payload: RunCreate, db: Session = Depends(get_db)) -> RunModel:
        run_id = payload.run_id or str(uuid4())
        if db.get(RunModel, run_id):
            raise HTTPException(status_code=409, detail=f"Run already exists: {run_id}")
        data = model_data(payload)
        data.pop("run_id", None)
        run = RunModel(id=run_id, **data)
        db.add(run)
        db.commit()
        db.refresh(run)
        return run

    @app.get("/runs", response_model=list[RunOut])
    def list_runs(status_filter: RunStatus | None = Query(default=None, alias="status"), db: Session = Depends(get_db)):
        stmt = select(RunModel).order_by(RunModel.created_at.desc())
        if status_filter:
            stmt = stmt.where(RunModel.status == status_filter)
        return list(db.scalars(stmt))

    @app.get("/runs/{run_id}", response_model=RunOut)
    def get_run(run_id: str, db: Session = Depends(get_db)) -> RunModel:
        run = db.get(RunModel, run_id)
        if not run:
            raise HTTPException(status_code=404, detail=f"Unknown run: {run_id}")
        return run

    @app.get("/runs/{run_id}/artifacts", response_model=list[ArtifactOut])
    def list_artifacts(run_id: str, db: Session = Depends(get_db)) -> list[dict]:
        run = db.get(RunModel, run_id)
        if not run:
            raise HTTPException(status_code=404, detail=f"Unknown run: {run_id}")
        return [
            {
                "artifact": artifact,
                "s3_uri": artifact_s3_uri(run, artifact),
                "presign_path": f"/runs/{run.id}/artifacts/{artifact}/presign",
            }
            for artifact in ARTIFACTS
        ]

    @app.get("/runs/{run_id}/artifacts/{artifact}", response_model=ArtifactOut)
    def get_artifact(run_id: str, artifact: ArtifactName, db: Session = Depends(get_db)) -> dict:
        run = db.get(RunModel, run_id)
        if not run:
            raise HTTPException(status_code=404, detail=f"Unknown run: {run_id}")
        return {
            "artifact": artifact,
            "s3_uri": artifact_s3_uri(run, artifact),
            "presign_path": f"/runs/{run.id}/artifacts/{artifact}/presign",
        }

    @app.patch("/runs/{run_id}", response_model=RunOut)
    def update_run(run_id: str, payload: RunUpdate, db: Session = Depends(get_db)) -> RunModel:
        run = db.get(RunModel, run_id)
        if not run:
            raise HTTPException(status_code=404, detail=f"Unknown run: {run_id}")
        for key, value in model_data(payload, exclude_unset=True).items():
            setattr(run, key, value)
        run.updated_at = utc_now()
        db.add(run)
        db.commit()
        db.refresh(run)
        return run

    @app.get("/runs/{run_id}/artifacts/{artifact}/presign")
    def presign_artifact(
        run_id: str,
        artifact: ArtifactName,
        expires: int = Query(default=3600, ge=60, le=604800),
        db: Session = Depends(get_db),
        s3_client=Depends(get_s3_client),
    ) -> dict:
        run = db.get(RunModel, run_id)
        if not run:
            raise HTTPException(status_code=404, detail=f"Unknown run: {run_id}")
        s3_uri = resolve_artifact_uri(run, artifact)
        bucket, key = parse_s3_uri(s3_uri)
        try:
            url = s3_client.generate_presigned_url(
                "get_object",
                Params={"Bucket": bucket, "Key": key},
                ExpiresIn=expires,
            )
        except Exception as exc:
            from botocore.exceptions import BotoCoreError, ClientError, NoCredentialsError

            if isinstance(exc, (BotoCoreError, ClientError, NoCredentialsError)):
                raise HTTPException(status_code=503, detail=f"S3 presigning unavailable: {exc.__class__.__name__}") from exc
            raise
        return {"artifact": artifact, "s3_uri": s3_uri, "url": url, "expires_in": expires}

    return app


app = create_app()
