from __future__ import annotations

import os
from datetime import datetime, timezone
from typing import Literal
from urllib.parse import urlparse
from uuid import uuid4

from fastapi import Depends, FastAPI, HTTPException, Query, status
from pydantic import BaseModel, Field
from sqlalchemy import Column, DateTime, String, create_engine, select, text
from sqlalchemy.orm import Session, declarative_base, sessionmaker


DEFAULT_DATABASE_URL = "sqlite:///./runs.db"
RunStatus = Literal["submitted", "running", "succeeded", "failed", "unknown"]
ArtifactName = Literal["report", "timeline", "trace", "dag", "multiqc"]

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


def model_data(model: BaseModel, *, exclude_unset: bool = False) -> dict:
    if hasattr(model, "model_dump"):
        return model.model_dump(exclude_unset=exclude_unset)
    return model.dict(exclude_unset=exclude_unset)


def engine_options(database_url: str) -> dict:
    if database_url.startswith("sqlite"):
        return {"connect_args": {"check_same_thread": False}}
    return {"pool_pre_ping": True}


def parse_s3_uri(uri: str) -> tuple[str, str]:
    parsed = urlparse(uri)
    if parsed.scheme != "s3" or not parsed.netloc or not parsed.path.strip("/"):
        raise HTTPException(status_code=400, detail=f"Invalid S3 URI: {uri}")
    return parsed.netloc, parsed.path.lstrip("/")


def resolve_artifact_uri(run: RunModel, artifact: ArtifactName) -> str:
    attr = "nextflow_report" if artifact == "report" else "multiqc_report" if artifact == "multiqc" else artifact
    value = getattr(run, attr)
    if value.startswith("s3://"):
        return value
    if not run.s3_prefix:
        raise HTTPException(status_code=404, detail="Run has no s3_prefix for relative artefact paths")
    return f"{run.s3_prefix.rstrip('/')}/{value.lstrip('/')}"


def get_s3_client():
    import boto3

    return boto3.client("s3", region_name=os.getenv("AWS_REGION"))


def create_app(database_url: str | None = None) -> FastAPI:
    db_url = database_url or os.getenv("DATABASE_URL", DEFAULT_DATABASE_URL)
    engine = create_engine(db_url, **engine_options(db_url))
    session_local = sessionmaker(bind=engine, autoflush=False, autocommit=False)
    Base.metadata.create_all(bind=engine)

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
        url = s3_client.generate_presigned_url(
            "get_object",
            Params={"Bucket": bucket, "Key": key},
            ExpiresIn=expires,
        )
        return {"artifact": artifact, "s3_uri": s3_uri, "url": url, "expires_in": expires}

    return app


app = create_app()
