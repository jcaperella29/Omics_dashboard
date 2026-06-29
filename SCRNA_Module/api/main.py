from __future__ import annotations

import csv
import json
import os
import re
import shutil
import subprocess
import tempfile
import time
import uuid
import zipfile
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any, Literal

from fastapi import BackgroundTasks, FastAPI, File, Form, HTTPException, UploadFile, status
from fastapi.responses import FileResponse
from pydantic import BaseModel, Field

APP_VERSION = "0.2.0"


class JobState(str, Enum):
    queued = "queued"
    running = "running"
    complete = "complete"
    failed = "failed"
    deleted = "deleted"


class EnrichBackend(str, Enum):
    gprof = "gprof"
    enrichr = "enrichr"


class Species(str, Enum):
    hsapiens = "hsapiens"
    mmusculus = "mmusculus"


class IdType(str, Enum):
    Symbol = "Symbol"
    Ensembl = "Ensembl"
    Entrez = "Entrez"


@dataclass(frozen=True)
class Settings:
    base_dir: Path = Path(__file__).resolve().parents[1]
    max_upload_mb: int = int(os.environ.get("SCRNA_MAX_UPLOAD_MB", "500"))
    cli_timeout_seconds: int = int(os.environ.get("SCRNA_CLI_TIMEOUT_SECONDS", "7200"))
    jobs_dir: Path = Path(os.environ.get("SCRNA_JOBS_DIR", Path(__file__).resolve().parents[1] / "jobs"))
    cli_script: Path = Path(os.environ.get("SCRNA_CLI_SCRIPT", Path(__file__).resolve().parents[1] / "cli" / "run_scrna.R"))


settings = Settings()
settings.jobs_dir.mkdir(parents=True, exist_ok=True)

app = FastAPI(
    title="JCAP scRNA-seq API",
    version=APP_VERSION,
    description="Production-oriented FastAPI wrapper for the scRNA-seq Seurat CLI module.",
)


class HealthResponse(BaseModel):
    status: Literal["ok", "degraded"]
    version: str
    cli_script: str
    cli_exists: bool
    jobs_dir: str
    jobs_dir_writable: bool


class SubmitResponse(BaseModel):
    job_id: str
    status: JobState
    status_url: str
    download_url: str | None = None
    tables_url: str | None = None


class JobStatusResponse(BaseModel):
    job_id: str
    status: JobState
    created_at: float | None = None
    started_at: float | None = None
    finished_at: float | None = None
    returncode: int | None = None
    error: str | None = None
    summary: dict[str, Any] = Field(default_factory=dict)
    run_status: dict[str, Any] = Field(default_factory=dict)
    has_bundle: bool = False
    download_url: str | None = None
    tables_url: str | None = None


class ValidationResponse(BaseModel):
    job_id: str
    counts_file: str
    meta_file: str
    counts_columns_preview: list[str]
    metadata_columns_preview: list[str]
    warnings: list[str] = Field(default_factory=list)
    notes: list[str]


SAFE_FILENAME_RE = re.compile(r"[^A-Za-z0-9._-]+")
ALLOWED_TABLE_SUFFIXES = {".csv"}
ALLOWED_UPLOAD_SUFFIXES = {".csv", ".tsv", ".txt", ".rds", ".h5", ".h5ad", ".mtx", ".gz"}


def _now() -> float:
    return time.time()


def _job_id() -> str:
    return str(uuid.uuid4())


def _parse_job_id(job_id: str) -> str:
    try:
        return str(uuid.UUID(job_id))
    except ValueError as exc:
        raise HTTPException(status_code=400, detail="Invalid job_id") from exc


def _job_dir(job_id: str) -> Path:
    return settings.jobs_dir / _parse_job_id(job_id)


def _safe_name(name: str, fallback: str = "upload.dat") -> str:
    raw = Path(name or fallback).name.strip() or fallback
    cleaned = SAFE_FILENAME_RE.sub("_", raw)
    return cleaned[:180]


def _ensure_allowed_upload_name(filename: str) -> None:
    suffixes = [s.lower() for s in Path(filename).suffixes]
    if not suffixes or not any(s in ALLOWED_UPLOAD_SUFFIXES for s in suffixes):
        raise HTTPException(status_code=415, detail=f"Unsupported upload type for {filename}")


def _atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", dir=path.parent, delete=False, encoding="utf-8") as tmp:
        json.dump(payload, tmp, indent=2, sort_keys=True)
        tmp.write("\n")
        tmp_path = Path(tmp.name)
    tmp_path.replace(path)


def _read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def _write_state(job_dir: Path, **updates: Any) -> dict[str, Any]:
    state_path = job_dir / "job_state.json"
    state = _read_json(state_path)
    state.update(updates)
    _atomic_write_json(state_path, state)
    return state


def _save_upload(upload: UploadFile, dest: Path) -> int:
    filename = _safe_name(upload.filename or dest.name)
    _ensure_allowed_upload_name(filename)
    max_bytes = settings.max_upload_mb * 1024 * 1024
    size = 0
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        with dest.open("wb") as fh:
            while chunk := upload.file.read(1024 * 1024):
                size += len(chunk)
                if size > max_bytes:
                    raise HTTPException(status_code=413, detail=f"File exceeds {settings.max_upload_mb} MB limit")
                fh.write(chunk)
    finally:
        upload.file.close()
    if size == 0:
        raise HTTPException(status_code=400, detail=f"Empty upload: {filename}")
    return size


def _csv_header(path: Path, delimiter: str = ",") -> list[str]:
    try:
        with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as fh:
            sample = fh.readline()
            if not sample:
                return []
            return [col.strip() for col in next(csv.reader([sample], delimiter=delimiter))]
    except Exception:
        return []


def _zip_dir(src_dir: Path, zip_path: Path) -> None:
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for path in src_dir.rglob("*"):
            if path.is_file() and path != zip_path:
                zf.write(path, path.relative_to(src_dir))


def _build_cmd(
    counts_path: Path,
    meta_path: Path,
    out_dir: Path,
    species: Species,
    id_type: IdType,
    enrich_backend: EnrichBackend,
    gprof_sources: str,
    enrichr_db: str,
    top_n_features: int,
    max_pcs: int,
    skip_pca_umap: bool,
    skip_enrichment: bool,
    skip_classifier: bool,
    skip_power: bool,
) -> list[str]:
    cmd = [
        "Rscript", str(settings.cli_script),
        "--counts", str(counts_path),
        "--metadata", str(meta_path),
        "--outdir", str(out_dir),
        "--species", species.value,
        "--id-type", id_type.value,
        "--enrich-backend", enrich_backend.value,
        "--gprof-sources", gprof_sources,
        "--enrichr-db", enrichr_db,
        "--top-n-features", str(top_n_features),
        "--max-pcs", str(max_pcs),
    ]
    if skip_pca_umap:
        cmd.append("--skip-pca-umap")
    if skip_enrichment:
        cmd.append("--skip-enrichment")
    if skip_classifier:
        cmd.append("--skip-classifier")
    if skip_power:
        cmd.append("--skip-power")
    return cmd


def _run_job(job_id: str, cmd: list[str]) -> None:
    job_dir = _job_dir(job_id)
    out_dir = job_dir / "results"
    run_log = job_dir / "run.log"
    _write_state(job_dir, status=JobState.running.value, started_at=_now())

    try:
        result = subprocess.run(
            cmd,
            cwd=str(job_dir),
            capture_output=True,
            text=True,
            timeout=settings.cli_timeout_seconds,
            check=False,
        )
        run_log.write_text(
            "COMMAND:\n" + " ".join(cmd) +
            "\n\nSTDOUT:\n" + result.stdout +
            "\n\nSTDERR:\n" + result.stderr,
            encoding="utf-8",
        )
        if result.returncode != 0:
            _write_state(
                job_dir,
                status=JobState.failed.value,
                finished_at=_now(),
                returncode=result.returncode,
                error=result.stderr[-4000:] or "CLI failed without stderr",
            )
            return

        bundle_path = job_dir / f"scrna_results_{job_id}.zip"
        _zip_dir(out_dir, bundle_path)
        _write_state(
            job_dir,
            status=JobState.complete.value,
            finished_at=_now(),
            returncode=0,
            summary=_read_json(out_dir / "run_summary.json"),
            run_status=_read_json(out_dir / "run_status.json"),
            download_url=f"/jobs/{job_id}/download",
            tables_url=f"/jobs/{job_id}/tables",
        )
    except subprocess.TimeoutExpired:
        _write_state(
            job_dir,
            status=JobState.failed.value,
            finished_at=_now(),
            error=f"CLI timed out after {settings.cli_timeout_seconds} seconds",
        )
    except Exception as exc:
        _write_state(job_dir, status=JobState.failed.value, finished_at=_now(), error=str(exc))


@app.get("/", response_model=HealthResponse)
def root() -> HealthResponse:
    return health()


@app.get("/health", response_model=HealthResponse)
def health() -> HealthResponse:
    cli_exists = settings.cli_script.exists()
    jobs_dir_writable = os.access(settings.jobs_dir, os.W_OK)
    overall = "ok" if cli_exists and jobs_dir_writable else "degraded"
    return HealthResponse(
        status=overall,
        version=APP_VERSION,
        cli_script=str(settings.cli_script),
        cli_exists=cli_exists,
        jobs_dir=str(settings.jobs_dir),
        jobs_dir_writable=jobs_dir_writable,
    )


@app.post("/validate", response_model=ValidationResponse)
def validate_inputs(
    counts_file: UploadFile = File(...),
    meta_file: UploadFile = File(...),
) -> ValidationResponse:
    job_id = _job_id()
    job_dir = _job_dir(job_id)
    input_dir = job_dir / "inputs"
    input_dir.mkdir(parents=True, exist_ok=True)

    counts_name = "counts_" + _safe_name(counts_file.filename or "counts.csv")
    meta_name = "metadata_" + _safe_name(meta_file.filename or "metadata.csv")
    counts_path = input_dir / counts_name
    meta_path = input_dir / meta_name
    _save_upload(counts_file, counts_path)
    _save_upload(meta_file, meta_path)

    counts_header = _csv_header(counts_path)
    meta_header = _csv_header(meta_path)
    warnings: list[str] = []
    if not counts_header:
        warnings.append("Could not read a CSV-style header from counts_file.")
    if not meta_header:
        warnings.append("Could not read a CSV-style header from meta_file.")
    for required in ("stim", "cell_type"):
        if meta_header and required not in meta_header:
            warnings.append(f"Metadata header does not contain required column: {required}")

    _write_state(job_dir, status=JobState.complete.value, created_at=_now(), validation_only=True)
    return ValidationResponse(
        job_id=job_id,
        counts_file=counts_path.name,
        meta_file=meta_path.name,
        counts_columns_preview=counts_header[:10],
        metadata_columns_preview=meta_header[:20],
        warnings=warnings,
        notes=[
            "Counts CSV should have genes as rows and cells as columns.",
            "Metadata CSV should have cells as rows; stim and cell_type are required for the full pipeline.",
        ],
    )


@app.post("/jobs", response_model=SubmitResponse, status_code=status.HTTP_202_ACCEPTED)
def submit_scrna_job(
    background_tasks: BackgroundTasks,
    counts_file: UploadFile = File(...),
    meta_file: UploadFile = File(...),
    species: Species = Form(Species.hsapiens),
    id_type: IdType = Form(IdType.Symbol),
    enrich_backend: EnrichBackend = Form(EnrichBackend.gprof),
    gprof_sources: str = Form("GO:BP,GO:MF,GO:CC,KEGG,REAC"),
    enrichr_db: str = Form("GO_Biological_Process_2023"),
    top_n_features: int = Form(10, ge=1, le=100),
    max_pcs: int = Form(30, ge=2, le=100),
    skip_pca_umap: bool = Form(False),
    skip_enrichment: bool = Form(False),
    skip_classifier: bool = Form(False),
    skip_power: bool = Form(False),
) -> SubmitResponse:
    if not settings.cli_script.exists():
        raise HTTPException(status_code=503, detail="scRNA CLI script is not available")

    job_id = _job_id()
    job_dir = _job_dir(job_id)
    input_dir = job_dir / "inputs"
    out_dir = job_dir / "results"
    input_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    counts_path = input_dir / ("counts_" + _safe_name(counts_file.filename or "counts.csv"))
    meta_path = input_dir / ("metadata_" + _safe_name(meta_file.filename or "metadata.csv"))
    counts_size = _save_upload(counts_file, counts_path)
    meta_size = _save_upload(meta_file, meta_path)

    cmd = _build_cmd(
        counts_path=counts_path,
        meta_path=meta_path,
        out_dir=out_dir,
        species=species,
        id_type=id_type,
        enrich_backend=enrich_backend,
        gprof_sources=gprof_sources,
        enrichr_db=enrichr_db,
        top_n_features=top_n_features,
        max_pcs=max_pcs,
        skip_pca_umap=skip_pca_umap,
        skip_enrichment=skip_enrichment,
        skip_classifier=skip_classifier,
        skip_power=skip_power,
    )
    _write_state(
        job_dir,
        job_id=job_id,
        status=JobState.queued.value,
        created_at=_now(),
        counts_file=counts_path.name,
        meta_file=meta_path.name,
        counts_size_bytes=counts_size,
        meta_size_bytes=meta_size,
        params={
            "species": species.value,
            "id_type": id_type.value,
            "enrich_backend": enrich_backend.value,
            "gprof_sources": gprof_sources,
            "enrichr_db": enrichr_db,
            "top_n_features": top_n_features,
            "max_pcs": max_pcs,
            "skip_pca_umap": skip_pca_umap,
            "skip_enrichment": skip_enrichment,
            "skip_classifier": skip_classifier,
            "skip_power": skip_power,
        },
    )
    background_tasks.add_task(_run_job, job_id, cmd)
    return SubmitResponse(job_id=job_id, status=JobState.queued, status_url=f"/jobs/{job_id}")


@app.post("/run", response_model=SubmitResponse, status_code=status.HTTP_202_ACCEPTED)
def backward_compatible_run(
    background_tasks: BackgroundTasks,
    counts_file: UploadFile = File(...),
    meta_file: UploadFile = File(...),
    species: Species = Form(Species.hsapiens),
    id_type: IdType = Form(IdType.Symbol),
    enrich_backend: EnrichBackend = Form(EnrichBackend.gprof),
    gprof_sources: str = Form("GO:BP,GO:MF,GO:CC,KEGG,REAC"),
    enrichr_db: str = Form("GO_Biological_Process_2023"),
    top_n_features: int = Form(10, ge=1, le=100),
    max_pcs: int = Form(30, ge=2, le=100),
    skip_pca_umap: bool = Form(False),
    skip_enrichment: bool = Form(False),
    skip_classifier: bool = Form(False),
    skip_power: bool = Form(False),
) -> SubmitResponse:
    return submit_scrna_job(
        background_tasks=background_tasks,
        counts_file=counts_file,
        meta_file=meta_file,
        species=species,
        id_type=id_type,
        enrich_backend=enrich_backend,
        gprof_sources=gprof_sources,
        enrichr_db=enrichr_db,
        top_n_features=top_n_features,
        max_pcs=max_pcs,
        skip_pca_umap=skip_pca_umap,
        skip_enrichment=skip_enrichment,
        skip_classifier=skip_classifier,
        skip_power=skip_power,
    )


@app.get("/jobs/{job_id}", response_model=JobStatusResponse)
def job_status(job_id: str) -> JobStatusResponse:
    job_dir = _job_dir(job_id)
    if not job_dir.exists():
        raise HTTPException(status_code=404, detail="Job not found")
    state = _read_json(job_dir / "job_state.json")
    out_dir = job_dir / "results"
    has_bundle = any(job_dir.glob("scrna_results_*.zip"))
    return JobStatusResponse(
        job_id=_parse_job_id(job_id),
        status=JobState(state.get("status", JobState.failed.value)),
        created_at=state.get("created_at"),
        started_at=state.get("started_at"),
        finished_at=state.get("finished_at"),
        returncode=state.get("returncode"),
        error=state.get("error"),
        summary=state.get("summary") or _read_json(out_dir / "run_summary.json"),
        run_status=state.get("run_status") or _read_json(out_dir / "run_status.json"),
        has_bundle=has_bundle,
        download_url=f"/jobs/{job_id}/download" if has_bundle else None,
        tables_url=f"/jobs/{job_id}/tables" if (out_dir / "tables").exists() else None,
    )


@app.get("/jobs/{job_id}/tables")
def list_tables(job_id: str) -> dict[str, Any]:
    tables_dir = _job_dir(job_id) / "results" / "tables"
    if not tables_dir.exists():
        raise HTTPException(status_code=404, detail="Tables not found")
    return {"job_id": _parse_job_id(job_id), "tables": sorted(p.name for p in tables_dir.glob("*.csv"))}


@app.get("/jobs/{job_id}/tables/{filename}")
def download_table(job_id: str, filename: str) -> FileResponse:
    clean = _safe_name(filename)
    path = _job_dir(job_id) / "results" / "tables" / clean
    if path.suffix.lower() not in ALLOWED_TABLE_SUFFIXES or not path.exists() or not path.is_file():
        raise HTTPException(status_code=404, detail="Table not found")
    return FileResponse(path, media_type="text/csv", filename=path.name)


@app.get("/jobs/{job_id}/download")
def download_bundle(job_id: str) -> FileResponse:
    job_dir = _job_dir(job_id)
    bundles = sorted(job_dir.glob("scrna_results_*.zip"))
    if not bundles:
        raise HTTPException(status_code=404, detail="Bundle not found")
    return FileResponse(bundles[0], media_type="application/zip", filename=bundles[0].name)


@app.delete("/jobs/{job_id}")
def delete_job(job_id: str) -> dict[str, Any]:
    job_dir = _job_dir(job_id)
    if not job_dir.exists():
        raise HTTPException(status_code=404, detail="Job not found")
    shutil.rmtree(job_dir)
    return {"job_id": _parse_job_id(job_id), "deleted": True}
