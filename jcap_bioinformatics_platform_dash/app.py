from __future__ import annotations

import base64
import io
import json
import os
import shutil
import tempfile
import uuid
import zipfile
from pathlib import Path
from typing import Any
from urllib.parse import urljoin

import pandas as pd
import plotly.graph_objects as go
import requests
from dash import Dash, Input, Output, State, dcc, dash_table, html, no_update

# -----------------------------------------------------------------------------
# Config
# -----------------------------------------------------------------------------
RNA_API_BASE = os.environ.get("RNA_API_BASE", "http://127.0.0.1:8000").rstrip("/")
NETWORK_API_BASE = os.environ.get("NETWORK_API_BASE", "http://127.0.0.1:8001").rstrip("/")
TRIAGE_API_BASE = os.environ.get("TRIAGE_API_BASE", "http://127.0.0.1:8002").rstrip("/")
WORK_ROOT = Path(os.environ.get("JCAP_PLATFORM_WORKDIR", tempfile.gettempdir())) / "jcap_bioinformatics_platform"
WORK_ROOT.mkdir(parents=True, exist_ok=True)

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

def make_run_dir(prefix: str) -> Path:
    d = WORK_ROOT / f"{prefix}_{uuid.uuid4().hex[:10]}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def decode_upload(contents: str, filename: str, outdir: Path) -> Path:
    if not contents:
        raise ValueError("No uploaded file contents.")
    _content_type, content_string = contents.split(",", 1)
    data = base64.b64decode(content_string)
    safe = Path(filename or f"upload_{uuid.uuid4().hex}.csv").name
    path = outdir / safe
    path.write_bytes(data)
    return path


def read_csv_preview(path: Path, n: int = 200) -> pd.DataFrame:
    try:
        return pd.read_csv(path).head(n)
    except Exception:
        return pd.DataFrame({"message": [f"Could not preview {path.name}"]})


def datatable_from_df(df: pd.DataFrame, page_size: int = 12):
    if df is None or df.empty:
        return html.Div("No rows to display.", className="help")
    return dash_table.DataTable(
        data=df.to_dict("records"),
        columns=[{"name": c, "id": c} for c in df.columns],
        page_size=page_size,
        sort_action="native",
        filter_action="native",
        style_table={"overflowX": "auto", "maxHeight": "65vh", "overflowY": "auto"},
        style_cell={"fontFamily": "Inter, Segoe UI, Arial", "fontSize": "0.85rem", "padding": "7px", "textAlign": "left", "maxWidth": 360, "whiteSpace": "normal"},
        style_header={"fontWeight": "800", "background": "#eef2ff"},
    )


def status(text: str, kind: str = "ok"):
    cls = {"ok": "status-ok", "warn": "status-warn", "bad": "status-bad"}.get(kind, "status-ok")
    return html.Div(text, className=cls)


def api_get_json(url: str) -> dict[str, Any]:
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    return r.json()


def api_post_json(url: str, payload: dict[str, Any]) -> dict[str, Any]:
    r = requests.post(url, json=payload, timeout=600)
    r.raise_for_status()
    return r.json()


def api_post_file(url: str, file_path: Path, *, file_field: str = "file", data: dict[str, Any] | None = None) -> dict[str, Any]:
    with file_path.open("rb") as fh:
        r = requests.post(url, files={file_field: (file_path.name, fh, "text/csv")}, data=data or {}, timeout=1200)
    r.raise_for_status()
    return r.json()


def api_post_file_binary(url: str, file_path: Path, *, file_field: str = "file", data: dict[str, Any] | None = None) -> bytes:
    with file_path.open("rb") as fh:
        r = requests.post(url, files={file_field: (file_path.name, fh, "text/csv")}, data=data or {}, timeout=1800)
    r.raise_for_status()
    return r.content


def write_json(path: Path, obj: Any) -> None:
    path.write_text(json.dumps(obj, indent=2, default=str), encoding="utf-8")


def figure_to_svg_bytes(fig_dict: dict[str, Any]) -> bytes:
    fig = go.Figure(fig_dict)
    return fig.to_image(format="svg", width=1400, height=1000, scale=1)


def download_file(path_str: str):
    path = Path(path_str)
    if not path.exists():
        return no_update
    return dcc.send_file(str(path))


def json_pretty(obj: Any) -> str:
    return json.dumps(obj, indent=2, default=str)[:20000]

def load_store_path(store: dict[str, Any] | None, key: str) -> Path | None:
    if not store or not store.get(key):
        return None
    p = Path(store[key])
    return p if p.exists() else None


def rna_output_dirs(store: dict[str, Any] | None) -> list[Path]:
    dirs: list[Path] = []
    if not store:
        return dirs
    for key in ("bundle_dir", "unzip_dir", "run_dir"):
        p = load_store_path(store, key)
        if p and p.is_dir() and p not in dirs:
            dirs.append(p)
    return dirs


def find_first_file(dirs: list[Path], patterns: list[str]) -> Path | None:
    for d in dirs:
        for pat in patterns:
            hits = sorted(d.rglob(pat))
            if hits:
                return hits[0]
    return None


def find_files(dirs: list[Path], patterns: list[str]) -> list[Path]:
    out: list[Path] = []
    seen = set()
    for d in dirs:
        for pat in patterns:
            for hit in sorted(d.rglob(pat)):
                if hit.is_file() and hit not in seen:
                    out.append(hit)
                    seen.add(hit)
    return out


def render_csv_file(path: Path | None, *, title: str | None = None):
    if not path or not path.exists():
        return html.Div("No matching table found yet. Fetch/extract the RNA-seq bundle or wait for the job to finish.", className="help")
    try:
        df = pd.read_csv(path)
    except Exception as exc:
        return html.Div([html.H4(path.name), html.Pre(str(exc), className="markdown-box")])
    return html.Div([
        html.H4(title or path.name),
        html.Div(str(path), className="help"),
        datatable_from_df(df, page_size=15),
    ])


def render_multiple_csvs(paths: list[Path]):
    if not paths:
        return html.Div("No matching CSV files found yet.", className="help")
    children = []
    for path in paths:
        children.append(render_csv_file(path))
    return html.Div(children)


def render_html_plot(path: Path | None):
    if not path or not path.exists():
        return html.Div("No matching plot HTML found yet. Fetch/extract the RNA-seq bundle or wait for the job to finish.", className="help")
    text = path.read_text(encoding="utf-8", errors="replace")
    return html.Div([
        html.H4(path.name),
        html.Div(str(path), className="help"),
        html.Iframe(srcDoc=text, style={"width": "100%", "height": "78vh", "border": "1px solid #cbd5e1", "borderRadius": "10px", "background": "white"}),
    ])


def save_and_extract_zip(content: bytes, run_dir: Path, zip_name: str = "rna_results_bundle.zip") -> tuple[Path, Path]:
    bundle_path = run_dir / zip_name
    bundle_path.write_bytes(content)
    bundle_dir = run_dir / "rna_bundle_unzipped"
    if bundle_dir.exists():
        shutil.rmtree(bundle_dir)
    bundle_dir.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(bundle_path) as z:
        z.extractall(bundle_dir)
    return bundle_path, bundle_dir


def candidate_rna_bundle_urls(job: dict[str, Any] | None, status_obj: dict[str, Any] | None = None) -> list[str]:
    job = job or {}
    status_obj = status_obj or {}
    candidates: list[str] = []
    for obj in (status_obj, job):
        for key in ("download_url", "bundle_url", "result_url", "results_url", "zip_url"):
            val = obj.get(key)
            if val:
                candidates.append(urljoin(RNA_API_BASE + "/", str(val).lstrip("/")))
    job_id = job.get("job_id") or status_obj.get("job_id")
    if job_id:
        candidates += [
            f"{RNA_API_BASE}/rnaseq/jobs/{job_id}/download",
            f"{RNA_API_BASE}/rnaseq/jobs/{job_id}/bundle",
            f"{RNA_API_BASE}/rnaseq/jobs/{job_id}/results.zip",
            f"{RNA_API_BASE}/rnaseq/jobs/{job_id}/results",
        ]
    # de-duplicate while preserving order
    seen = set()
    out = []
    for url in candidates:
        if url not in seen:
            out.append(url)
            seen.add(url)
    return out


def try_fetch_rna_bundle(job: dict[str, Any] | None, status_obj: dict[str, Any] | None, run_dir: Path) -> tuple[Path | None, Path | None, str | None]:
    errors = []
    for url in candidate_rna_bundle_urls(job, status_obj):
        try:
            r = requests.get(url, timeout=1200)
            r.raise_for_status()
            content_type = r.headers.get("content-type", "")
            # Prefer zip responses. Some APIs may return JSON status from /results; skip those.
            if r.content[:2] != b"PK" and "zip" not in content_type.lower():
                errors.append(f"{url}: response was not a zip")
                continue
            bundle_path, bundle_dir = save_and_extract_zip(r.content, run_dir)
            return bundle_path, bundle_dir, url
        except Exception as exc:
            errors.append(f"{url}: {exc}")
    return None, None, " | ".join(errors[:3]) if errors else "No RNA bundle URL candidates were available."

# -----------------------------------------------------------------------------
# Layout fragments
# -----------------------------------------------------------------------------

def upload_box(id_: str, label: str, multiple: bool = False):
    return dcc.Upload(
        id=id_,
        multiple=multiple,
        className="upload-box",
        children=html.Div([html.Strong(label), html.Div("Drag/drop or click to select file", className="help")]),
    )


def uploaded_file_badges(*names):
    clean = [n for n in names if n]
    if not clean:
        return html.Div("No files selected yet.", className="help")
    return html.Div([html.Span(Path(n).name, className="file-badge") for n in clean], className="uploaded-list")


def omics_layout():
    return html.Div([
        dcc.Tabs(id="omics-subtabs", value="rna-seq", children=[dcc.Tab(label="RNA-seq", value="rna-seq")]),
        html.Div(className="card", children=[
            html.H3("RNA-seq module"),
            html.Div("Upload counts and phenotype metadata, submit the job to the RNA-seq API, then inspect and download the returned outputs.", className="help"),
            html.Div(className="grid-2", children=[
                html.Div([html.Label("Counts matrix", className="label"), upload_box("rna-counts-upload", "Upload counts_data.csv")]),
                html.Div([html.Label("Phenotype metadata", className="label"), upload_box("rna-pheno-upload", "Upload phenotype_data.csv")]),
            ]),
            html.Div(id="rna-uploaded-files", style={"marginTop": "8px"}),
            html.Div(className="grid-3", style={"marginTop": "12px"}, children=[
                html.Div([html.Label("Phenotype column", className="label"), dcc.Input(id="rna-phenotype-column", value="dex", type="text", style={"width": "100%"})]),
                html.Div([html.Label("Species", className="label"), dcc.Dropdown(id="rna-species", value="hsapiens", clearable=False, options=[{"label": x, "value": x} for x in ["hsapiens", "mmusculus", "drerio", "dmelanogaster"]])]),
                html.Div([html.Label("CLI script", className="label"), dcc.Input(id="rna-cli-script", value="rnaseq_cli.R", type="text", style={"width": "100%"})]),
            ]),
            html.Div(style={"marginTop": "12px"}, children=[
                dcc.Checklist(id="rna-run-options", value=[], options=[
                    {"label": " Skip enrichment", "value": "skip_enrichment"},
                    {"label": " Skip classifier", "value": "skip_classifier"},
                    {"label": " Skip plots", "value": "skip_plots"},
                ], inline=True),
                html.Button("Run RNA-seq API", id="rna-run-btn", n_clicks=0, className="btn"),
                html.Button("Poll job", id="rna-poll-btn", n_clicks=0, className="btn btn-secondary"),
                html.Button("Fetch/extract RNA bundle", id="rna-fetch-bundle-btn", n_clicks=0, className="btn btn-secondary"),
                html.Button("Download saved bundle", id="rna-download-bundle-btn", n_clicks=0, className="btn btn-secondary"),
                html.Button("Download current table CSV", id="rna-download-current-csv-btn", n_clicks=0, className="btn btn-secondary"),
                html.Button("Download current plot HTML", id="rna-download-current-plot-btn", n_clicks=0, className="btn btn-secondary"),
                dcc.Download(id="rna-download"),
            ]),
            html.Div(id="rna-status", style={"marginTop": "10px"}),
            dcc.Store(id="rna-store"),
        ]),
        html.Div(className="card", children=[
            html.H3("RNA-seq outputs"),
            dcc.Tabs(id="rna-output-tabs", value="rna-de-results", children=[
                dcc.Tab(label="DE results", value="rna-de-results"),
                dcc.Tab(label="Full limma table", value="rna-limma"),
                dcc.Tab(label="Enrichment", value="rna-enrichment"),
                dcc.Tab(label="Classifier", value="rna-classifier"),
                dcc.Tab(label="Power analysis", value="rna-power"),
                dcc.Tab(label="Plots", value="rna-plots"),
                dcc.Tab(label="Output files", value="rna-output-files"),
                dcc.Tab(label="Job JSON", value="rna-job-json"),
                dcc.Tab(label="Status JSON", value="rna-status-json"),
            ]),
            html.Div(id="rna-output-panel", style={"marginTop": "12px"}),
        ]),
    ])


def network_layout():
    return html.Div([
        html.Div(className="card", children=[
            html.H3("Network visualization module"),
            html.Div("Upload an enrichment edge file and run graph build, filtering, diffusion, projection, and consensus candidate scoring through the Network API.", className="help"),
            html.Div(className="grid-3", children=[
                html.Div([html.Label("Enrichment edge file", className="label"), upload_box("network-edge-upload", "Upload *_network_edges.csv")]),
                html.Div([html.Label("Max groups in plot", className="label"), dcc.Input(id="network-max-groups", value=50, type="number", style={"width": "100%"})]),
                html.Div([html.Label("Layout", className="label"), dcc.Dropdown(id="network-layout-mode", value="bipartite", clearable=False, options=[{"label": "Bipartite", "value": "bipartite"}, {"label": "Force-directed", "value": "force"}])]),
            ]),
            html.Div(id="network-uploaded-files", style={"marginTop": "8px"}),
            html.Div(style={"marginTop": "12px"}, children=[
                dcc.Checklist(id="network-options", value=["thickness"], inline=True, options=[
                    {"label": " Edge thickness by weight", "value": "thickness"},
                    {"label": " Show labels", "value": "labels"},
                    {"label": " Largest component only", "value": "largest"},
                ]),
                html.Button("Run Network API", id="network-run-btn", n_clicks=0, className="btn"),
                html.Button("Download graph SVG", id="network-download-main-svg-btn", n_clicks=0, className="btn btn-secondary"),
                html.Button("Download projection SVG", id="network-download-proj-svg-btn", n_clicks=0, className="btn btn-secondary"),
                html.Button("Download current network table CSV", id="network-download-current-csv-btn", n_clicks=0, className="btn btn-secondary"),
                html.Button("Download all network tables ZIP", id="network-download-tables-zip-btn", n_clicks=0, className="btn btn-secondary"),
                dcc.Download(id="network-download"),
            ]),
            html.Div(id="network-status", style={"marginTop": "10px"}),
            dcc.Store(id="network-store"),
        ]),
        html.Div(className="card", children=[
            html.H3("Network outputs"),
            dcc.Tabs(id="network-output-tabs", value="network-main-graph", children=[
                dcc.Tab(label="Main graph", value="network-main-graph"),
                dcc.Tab(label="Projection graph", value="network-projection-graph"),
                dcc.Tab(label="Bipartite candidates", value="network-bip-candidates"),
                dcc.Tab(label="Projection candidates", value="network-proj-candidates"),
                dcc.Tab(label="Consensus", value="network-consensus"),
                dcc.Tab(label="Nodes/edges", value="network-nodes-edges"),
                dcc.Tab(label="Raw JSON", value="network-json"),
            ]),
            html.Div(id="network-output-panel", style={"marginTop": "12px"}),
        ]),
    ])


def triage_layout():
    return html.Div([
        html.Div(className="card", children=[
            html.H3("LLM triage module"),
            html.Div("Upload enrichment data, provide experiment context, and run deterministic pretriage or full PubMed-aware LLM interpretation.", className="help"),
            html.Div([html.Label("Enrichment / network edge file", className="label"), upload_box("triage-upload", "Upload enrichment CSV or *_network_edges.csv")]),
            html.Div(id="triage-uploaded-files", style={"marginTop": "8px"}),
            html.Div(className="grid-3", style={"marginTop": "12px"}, children=[
                html.Div([html.Label("Phenotype", className="label"), dcc.Input(id="triage-phenotype", value="dexamethasone treatment response vs untreated control", type="text", style={"width": "100%"})]),
                html.Div([html.Label("Organism", className="label"), dcc.Input(id="triage-organism", value="hsapiens", type="text", style={"width": "100%"})]),
                html.Div([html.Label("Assay", className="label"), dcc.Input(id="triage-assay", value="bulk_rnaseq", type="text", style={"width": "100%"})]),
                html.Div([html.Label("Tissue", className="label"), dcc.Input(id="triage-tissue", value="airway smooth muscle", type="text", style={"width": "100%"})]),
                html.Div([html.Label("Cell type", className="label"), dcc.Input(id="triage-cell-type", value="human airway smooth muscle cell lines", type="text", style={"width": "100%"})]),
                html.Div([html.Label("Perturbation", className="label"), dcc.Input(id="triage-perturbation", value="dexamethasone", type="text", style={"width": "100%"})]),
                html.Div([html.Label("Timepoint", className="label"), dcc.Input(id="triage-timepoint", value="18 hours", type="text", style={"width": "100%"})]),
                html.Div([html.Label("Input format", className="label"), dcc.Dropdown(id="triage-input-format", value="long_edges", clearable=False, options=[{"label": "Auto", "value": "auto"}, {"label": "Long edges", "value": "long_edges"}, {"label": "Enrichment table", "value": "enrichment"}])]),
                html.Div([html.Label("Mode", className="label"), dcc.Dropdown(id="triage-mode", value="pretriage", clearable=False, options=[{"label": "Pretriage/no LLM", "value": "pretriage"}, {"label": "Full PubMed + LLM", "value": "full"}])]),
            ]),
            html.Div(style={"marginTop": "12px"}, children=[
                html.Button("Run Triage API", id="triage-run-btn", n_clicks=0, className="btn"),
                html.Button("Download bundle", id="triage-download-bundle-btn", n_clicks=0, className="btn btn-secondary"),
                html.Button("Download markdown", id="triage-download-md-btn", n_clicks=0, className="btn btn-secondary"),
                html.Button("Download PDF", id="triage-download-pdf-btn", n_clicks=0, className="btn btn-secondary"),
                dcc.Download(id="triage-download"),
            ]),
            html.Div(id="triage-status", style={"marginTop": "10px"}),
            dcc.Store(id="triage-store"),
        ]),
        html.Div(className="card", children=[
            html.H3("Triage outputs"),
            dcc.Tabs(id="triage-output-tabs", value="triage-report-md", children=[
                dcc.Tab(label="Markdown report", value="triage-report-md"),
                dcc.Tab(label="Programs", value="triage-programs"),
                dcc.Tab(label="Triage rows", value="triage-rows"),
                dcc.Tab(label="Claims", value="triage-claims"),
                dcc.Tab(label="Files", value="triage-files"),
                dcc.Tab(label="Raw JSON", value="triage-json"),
            ]),
            html.Div(id="triage-output-panel", style={"marginTop": "12px"}),
        ]),
    ])

# -----------------------------------------------------------------------------
# App
# -----------------------------------------------------------------------------
app = Dash(__name__, suppress_callback_exceptions=True, title="JCAP Bioinformatics Platform")
server = app.server

app.layout = html.Div(className="platform-shell", children=[
    html.Div(className="hero", children=[
        html.H1("JCAP Bioinformatics Platform"),
        html.P("A modular multi-omics dashboard for running reproducible analysis engines, network visualization, and evidence-aware LLM triage. Upload data once, send it through specialized APIs, inspect tables and plots, and download complete result bundles."),
    ]),
    dcc.Tabs(id="main-tabs", value="omics", children=[
        dcc.Tab(label="Omics", value="omics"),
        dcc.Tab(label="Network Visualization", value="network"),
        dcc.Tab(label="LLM Triage", value="triage"),
    ]),
    html.Div(id="main-content"),
])

@app.callback(Output("main-content", "children"), Input("main-tabs", "value"))
def render_main(tab):
    if tab == "network":
        return network_layout()
    if tab == "triage":
        return triage_layout()
    return omics_layout()


# -----------------------------------------------------------------------------
# Uploaded filename display callbacks
# -----------------------------------------------------------------------------
@app.callback(
    Output("rna-uploaded-files", "children"),
    Input("rna-counts-upload", "filename"),
    Input("rna-pheno-upload", "filename"),
)
def show_rna_uploaded_files(counts_name, pheno_name):
    return html.Div([html.Strong("Selected files: "), uploaded_file_badges(counts_name, pheno_name)])


@app.callback(Output("network-uploaded-files", "children"), Input("network-edge-upload", "filename"))
def show_network_uploaded_files(edge_name):
    return html.Div([html.Strong("Selected file: "), uploaded_file_badges(edge_name)])


@app.callback(Output("triage-uploaded-files", "children"), Input("triage-upload", "filename"))
def show_triage_uploaded_files(name):
    return html.Div([html.Strong("Selected file: "), uploaded_file_badges(name)])

# -----------------------------------------------------------------------------
# RNA callbacks
# -----------------------------------------------------------------------------
@app.callback(
    Output("rna-status", "children"),
    Output("rna-store", "data"),
    Input("rna-run-btn", "n_clicks"),
    State("rna-counts-upload", "contents"), State("rna-counts-upload", "filename"),
    State("rna-pheno-upload", "contents"), State("rna-pheno-upload", "filename"),
    State("rna-phenotype-column", "value"), State("rna-species", "value"), State("rna-cli-script", "value"), State("rna-run-options", "value"),
    prevent_initial_call=True,
)
def run_rna(n, counts_contents, counts_name, pheno_contents, pheno_name, phenotype_col, species, cli_script, options):
    if not counts_contents or not pheno_contents:
        return status("Upload both counts and phenotype files.", "bad"), no_update
    run_dir = make_run_dir("rna")
    try:
        counts_path = decode_upload(counts_contents, counts_name, run_dir)
        pheno_path = decode_upload(pheno_contents, pheno_name, run_dir)
        form = {
            "phenotype_column": phenotype_col or "dex",
            "species": species or "hsapiens",
            "cli_script": cli_script or "rnaseq_cli.R",
            "skip_enrichment": "true" if "skip_enrichment" in (options or []) else "false",
            "skip_classifier": "true" if "skip_classifier" in (options or []) else "false",
            "skip_plots": "true" if "skip_plots" in (options or []) else "false",
        }
        with counts_path.open("rb") as cfh, pheno_path.open("rb") as pfh:
            r = requests.post(
                f"{RNA_API_BASE}/rnaseq/jobs",
                files={"counts": (counts_path.name, cfh, "text/csv"), "phenotype": (pheno_path.name, pfh, "text/csv")},
                data=form,
                timeout=600,
            )
        r.raise_for_status()
        job = r.json()
        write_json(run_dir / "rna_job_response.json", job)
        store = {"run_dir": str(run_dir), "job": job, "status": None, "bundle_path": None, "api_base": RNA_API_BASE}
        return status(f"RNA-seq job submitted. job_id={job.get('job_id', 'unknown')}", "ok"), store
    except Exception as e:
        return status(f"RNA-seq API error: {e}", "bad"), {"run_dir": str(run_dir), "error": str(e)}

@app.callback(
    Output("rna-status", "children", allow_duplicate=True),
    Output("rna-store", "data", allow_duplicate=True),
    Input("rna-poll-btn", "n_clicks"),
    State("rna-store", "data"),
    prevent_initial_call=True,
)
def poll_rna(n, store):
    if not store or not store.get("job"):
        return status("No RNA-seq job has been submitted yet.", "warn"), no_update
    job = store.get("job", {})
    job_id = job.get("job_id")
    candidates = []
    if job.get("status_url"):
        candidates.append(urljoin(RNA_API_BASE + "/", str(job["status_url"]).lstrip("/")))
    if job_id:
        candidates += [f"{RNA_API_BASE}/rnaseq/jobs/{job_id}", f"{RNA_API_BASE}/rnaseq/jobs/{job_id}/status"]
    errors = []
    for url in candidates:
        try:
            st = api_get_json(url)
            store["status"] = st
            write_json(Path(store["run_dir"]) / "rna_status_response.json", st)
            bundle_path, bundle_dir, bundle_msg = try_fetch_rna_bundle(store.get("job"), st, Path(store["run_dir"]))
            if bundle_path and bundle_dir:
                store["bundle_path"] = str(bundle_path)
                store["bundle_dir"] = str(bundle_dir)
                return status(f"RNA-seq status updated and bundle extracted from {bundle_msg}", "ok"), store
            return status(f"RNA-seq status updated from {url}. Bundle not available yet: {bundle_msg}", "warn"), store
        except Exception as e:
            errors.append(f"{url}: {e}")
    return status("Could not poll RNA-seq job. Configure status_url or endpoint. " + " | ".join(errors[:2]), "bad"), store

@app.callback(
    Output("rna-status", "children", allow_duplicate=True),
    Output("rna-store", "data", allow_duplicate=True),
    Input("rna-fetch-bundle-btn", "n_clicks"),
    State("rna-store", "data"),
    prevent_initial_call=True,
)
def fetch_rna_bundle(n, store):
    if not store or not store.get("job"):
        return status("No RNA-seq job has been submitted yet.", "warn"), no_update
    run_dir = Path(store.get("run_dir", ""))
    if not run_dir.exists():
        return status("Local RNA run directory no longer exists.", "bad"), no_update
    bundle_path, bundle_dir, msg = try_fetch_rna_bundle(store.get("job"), store.get("status"), run_dir)
    if bundle_path and bundle_dir:
        store["bundle_path"] = str(bundle_path)
        store["bundle_dir"] = str(bundle_dir)
        return status(f"RNA-seq bundle extracted from {msg}", "ok"), store
    return status(f"Could not fetch RNA-seq bundle yet. {msg}", "bad"), store

@app.callback(Output("rna-output-panel", "children"), Input("rna-output-tabs", "value"), State("rna-store", "data"))
def render_rna_output(tab, store):
    if not store:
        return html.Div("Run an RNA-seq job first.", className="help")
    if store.get("error"):
        return html.Pre(store["error"], className="markdown-box")
    dirs = rna_output_dirs(store)

    if tab == "rna-de-results":
        return render_csv_file(find_first_file(dirs, ["DE_results_*.csv", "DE_results.csv"]), title="Differential expression final results")
    if tab == "rna-limma":
        return render_csv_file(find_first_file(dirs, ["limma_all_results_*.csv", "limma_all_results.csv"]), title="Full limma/voom results")
    if tab == "rna-enrichment":
        paths = find_files(dirs, ["enrichment*.csv", "*network_edges.csv"])
        return render_multiple_csvs(paths)
    if tab == "rna-classifier":
        paths = find_files(dirs, ["rf_predictions_*.csv", "rf_metrics_*.csv", "rf*.csv"])
        return render_multiple_csvs(paths)
    if tab == "rna-power":
        paths = find_files(dirs, ["power_summary_*.csv", "power*.csv"])
        return render_multiple_csvs(paths)
    if tab == "rna-plots":
        plot_files = find_files(dirs, ["pca_plot.html", "umap_plot.html", "volcano_plot.html", "heatmap_plot.html", "power_curve_plot.html", "rf_roc_plot.html", "*plot.html"])
        if not plot_files:
            return html.Div("No plot HTML files found yet.", className="help")
        return html.Div([
            html.Div("Showing the first plot found. Use Output files to see all available plot files.", className="help"),
            render_html_plot(plot_files[0]),
        ])
    if tab == "rna-status-json":
        return html.Pre(json_pretty(store.get("status") or {}), className="markdown-box")
    if tab == "rna-output-files":
        rows = []
        for d in dirs:
            for p in sorted(d.rglob("*")):
                if p.is_file():
                    rows.append({"file": p.name, "relative_path": str(p.relative_to(d)), "size_kb": round(p.stat().st_size / 1024, 2), "root": str(d)})
        return datatable_from_df(pd.DataFrame(rows), page_size=20)
    return html.Pre(json_pretty(store.get("job") or store), className="markdown-box")

@app.callback(Output("rna-download", "data"), Input("rna-download-bundle-btn", "n_clicks"), State("rna-store", "data"), prevent_initial_call=True)
def download_rna_bundle(n, store):
    if not store:
        return no_update
    existing = load_store_path(store, "bundle_path")
    if existing and existing.exists():
        return dcc.send_file(str(existing))
    # Last chance: try to fetch and download immediately. This will not update UI state,
    # but it lets the user download when the API endpoint is ready.
    run_dir = Path(store.get("run_dir", tempfile.gettempdir()))
    bundle_path, _bundle_dir, _msg = try_fetch_rna_bundle(store.get("job"), store.get("status"), run_dir)
    if bundle_path:
        return dcc.send_file(str(bundle_path))
    return no_update


@app.callback(Output("rna-download", "data", allow_duplicate=True), Input("rna-download-current-csv-btn", "n_clicks"), State("rna-output-tabs", "value"), State("rna-store", "data"), prevent_initial_call=True)
def download_current_rna_csv(n, tab, store):
    dirs = rna_output_dirs(store)
    path = None
    if tab == "rna-de-results":
        path = find_first_file(dirs, ["DE_results_*.csv", "DE_results.csv"])
    elif tab == "rna-limma":
        path = find_first_file(dirs, ["limma_all_results_*.csv", "limma_all_results.csv"])
    elif tab == "rna-power":
        path = find_first_file(dirs, ["power_summary_*.csv", "power*.csv"])
    elif tab == "rna-classifier":
        path = find_first_file(dirs, ["rf_metrics_*.csv", "rf_predictions_*.csv", "rf*.csv"])
    elif tab == "rna-enrichment":
        path = find_first_file(dirs, ["enrichment*.csv", "*network_edges.csv"])
    if path and path.exists():
        return dcc.send_file(str(path))
    return no_update


@app.callback(Output("rna-download", "data", allow_duplicate=True), Input("rna-download-current-plot-btn", "n_clicks"), State("rna-store", "data"), prevent_initial_call=True)
def download_current_rna_plot(n, store):
    dirs = rna_output_dirs(store)
    path = find_first_file(dirs, ["pca_plot.html", "umap_plot.html", "volcano_plot.html", "heatmap_plot.html", "power_curve_plot.html", "rf_roc_plot.html", "*plot.html"])
    if path and path.exists():
        return dcc.send_file(str(path))
    return no_update

# -----------------------------------------------------------------------------
# Network callbacks
# -----------------------------------------------------------------------------
@app.callback(
    Output("network-status", "children"), Output("network-store", "data"),
    Input("network-run-btn", "n_clicks"),
    State("network-edge-upload", "contents"), State("network-edge-upload", "filename"),
    State("network-max-groups", "value"), State("network-layout-mode", "value"), State("network-options", "value"),
    prevent_initial_call=True,
)
def run_network(n, contents, filename, max_groups, layout_mode, opts):
    if not contents:
        return status("Upload an enrichment edge file first.", "bad"), no_update
    run_dir = make_run_dir("network")
    opts = opts or []
    try:
        edge_path = decode_upload(contents, filename, run_dir)
        build = api_post_file(
            f"{NETWORK_API_BASE}/network/build",
            edge_path,
            data={"options_json": json.dumps({"apply_preset": False, "item_col": "gene", "group_col": "term", "weight_col": "adjusted_pvalue"})},
        )
        write_json(run_dir / "01_build.json", build)
        graph = build["graph"]
        filtered = api_post_json(f"{NETWORK_API_BASE}/network/filter", {
            "graph": graph,
            "options": {
                "search": "", "min_degree": 1, "min_weight": 0.0,
                "max_groups": int(max_groups or 50), "largest_component_only": "largest" in opts,
                "layout_mode": layout_mode or "bipartite", "show_labels": "labels" in opts,
                "thickness_by_weight": "thickness" in opts, "return_figure": True,
            },
        })
        write_json(run_dir / "02_filter.json", filtered)
        bip = api_post_json(f"{NETWORK_API_BASE}/diffusion/bipartite", {"graph": filtered["graph"], "top_n": 50, "candidate_top_n": 30})
        write_json(run_dir / "03_bipartite.json", bip)
        proj = api_post_json(f"{NETWORK_API_BASE}/projection/build", {"graph": filtered["graph"], "method": "jaccard", "return_figure": True, "show_labels": True})
        write_json(run_dir / "04_projection.json", proj)
        pdiff = api_post_json(f"{NETWORK_API_BASE}/diffusion/projection", {"projection_graph": proj["projection_graph"], "top_n": 50, "candidate_top_n": 30})
        write_json(run_dir / "05_projection_diffusion.json", pdiff)
        consensus = api_post_json(f"{NETWORK_API_BASE}/consensus", {"bipartite_results": bip, "projection_results": pdiff, "top_n": 30})
        write_json(run_dir / "06_consensus.json", consensus)
        # Save CSV outputs
        (run_dir / "bipartite_candidates.csv").write_text(bip.get("candidate_csv", ""), encoding="utf-8")
        (run_dir / "projection_candidates.csv").write_text(pdiff.get("candidate_csv", ""), encoding="utf-8")
        (run_dir / "consensus_candidates.csv").write_text(consensus.get("csv", ""), encoding="utf-8")
        pd.DataFrame(filtered.get("graph", {}).get("nodes", [])).to_csv(run_dir / "filtered_nodes.csv", index=False)
        pd.DataFrame(filtered.get("graph", {}).get("edges", [])).to_csv(run_dir / "filtered_edges.csv", index=False)
        pd.DataFrame(proj.get("projection_graph", {}).get("nodes", [])).to_csv(run_dir / "projection_nodes.csv", index=False)
        pd.DataFrame(proj.get("projection_graph", {}).get("edges", [])).to_csv(run_dir / "projection_edges.csv", index=False)
        pd.DataFrame(bip.get("rows", [])).to_csv(run_dir / "bipartite_diffusion_rows.csv", index=False)
        pd.DataFrame(pdiff.get("rows", [])).to_csv(run_dir / "projection_diffusion_rows.csv", index=False)
        store = {"run_dir": str(run_dir), "build": build, "filtered": filtered, "bip": bip, "projection": proj, "pdiff": pdiff, "consensus": consensus}
        return status("Network API run complete.", "ok"), store
    except Exception as e:
        return status(f"Network API error: {e}", "bad"), {"run_dir": str(run_dir), "error": str(e)}

@app.callback(Output("network-output-panel", "children"), Input("network-output-tabs", "value"), State("network-store", "data"))
def render_network_output(tab, store):
    if not store:
        return html.Div("Run the Network API first.", className="help")
    if store.get("error"):
        return html.Pre(store["error"], className="markdown-box")
    if tab == "network-main-graph":
        return dcc.Graph(figure=store.get("filtered", {}).get("figure", {}), style={"height": "75vh"})
    if tab == "network-projection-graph":
        return dcc.Graph(figure=store.get("projection", {}).get("figure", {}), style={"height": "75vh"})
    if tab == "network-bip-candidates":
        return datatable_from_df(pd.DataFrame(store.get("bip", {}).get("candidate_rows", [])))
    if tab == "network-proj-candidates":
        return datatable_from_df(pd.DataFrame(store.get("pdiff", {}).get("candidate_rows", [])))
    if tab == "network-consensus":
        return datatable_from_df(pd.DataFrame(store.get("consensus", {}).get("candidate_rows", [])))
    if tab == "network-nodes-edges":
        nodes = pd.DataFrame(store.get("filtered", {}).get("graph", {}).get("nodes", []))
        edges = pd.DataFrame(store.get("filtered", {}).get("graph", {}).get("edges", []))
        return html.Div([html.H4("Nodes"), datatable_from_df(nodes), html.H4("Edges"), datatable_from_df(edges)])
    return html.Pre(json_pretty(store), className="markdown-box")

@app.callback(Output("network-download", "data"), Input("network-download-main-svg-btn", "n_clicks"), State("network-store", "data"), prevent_initial_call=True)
def download_main_svg(n, store):
    if not store or not store.get("filtered", {}).get("figure"):
        return no_update
    path = Path(store["run_dir"]) / "main_network.svg"
    path.write_bytes(figure_to_svg_bytes(store["filtered"]["figure"]))
    return dcc.send_file(str(path))

@app.callback(Output("network-download", "data", allow_duplicate=True), Input("network-download-proj-svg-btn", "n_clicks"), State("network-store", "data"), prevent_initial_call=True)
def download_proj_svg(n, store):
    if not store or not store.get("projection", {}).get("figure"):
        return no_update
    path = Path(store["run_dir"]) / "projection_network.svg"
    path.write_bytes(figure_to_svg_bytes(store["projection"]["figure"]))
    return dcc.send_file(str(path))


@app.callback(
    Output("network-download", "data", allow_duplicate=True),
    Input("network-download-current-csv-btn", "n_clicks"),
    State("network-output-tabs", "value"),
    State("network-store", "data"),
    prevent_initial_call=True,
)
def download_current_network_csv(n, tab, store):
    if not store or not store.get("run_dir"):
        return no_update
    run_dir = Path(store["run_dir"])
    mapping = {
        "network-bip-candidates": "bipartite_candidates.csv",
        "network-proj-candidates": "projection_candidates.csv",
        "network-consensus": "consensus_candidates.csv",
        "network-nodes-edges": "filtered_nodes.csv",
    }
    path = run_dir / mapping.get(tab, "consensus_candidates.csv")
    if path.exists():
        return dcc.send_file(str(path))
    return no_update


@app.callback(
    Output("network-download", "data", allow_duplicate=True),
    Input("network-download-tables-zip-btn", "n_clicks"),
    State("network-store", "data"),
    prevent_initial_call=True,
)
def download_all_network_tables_zip(n, store):
    if not store or not store.get("run_dir"):
        return no_update
    run_dir = Path(store["run_dir"])
    names = [
        "bipartite_candidates.csv",
        "projection_candidates.csv",
        "consensus_candidates.csv",
        "filtered_nodes.csv",
        "filtered_edges.csv",
        "projection_nodes.csv",
        "projection_edges.csv",
        "bipartite_diffusion_rows.csv",
        "projection_diffusion_rows.csv",
    ]
    zip_path = run_dir / "network_tables.zip"
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as z:
        for name in names:
            table_path = run_dir / name
            if table_path.exists():
                z.write(table_path, arcname=name)
    return dcc.send_file(str(zip_path))

# -----------------------------------------------------------------------------
# Triage callbacks
# -----------------------------------------------------------------------------
@app.callback(
    Output("triage-status", "children"), Output("triage-store", "data"),
    Input("triage-run-btn", "n_clicks"),
    State("triage-upload", "contents"), State("triage-upload", "filename"),
    State("triage-phenotype", "value"), State("triage-organism", "value"), State("triage-assay", "value"),
    State("triage-tissue", "value"), State("triage-cell-type", "value"), State("triage-perturbation", "value"), State("triage-timepoint", "value"),
    State("triage-input-format", "value"), State("triage-mode", "value"),
    prevent_initial_call=True,
)
def run_triage(n, contents, filename, phenotype, organism, assay, tissue, cell_type, perturbation, timepoint, input_format, mode):
    if not contents:
        return status("Upload an enrichment file first.", "bad"), no_update
    run_dir = make_run_dir("triage")
    try:
        input_path = decode_upload(contents, filename, run_dir)
        form = {
            "input_format": input_format or "auto", "mode": mode or "pretriage", "make_pdf": "true",
            "phenotype": phenotype or "", "organism": organism or "", "assay": assay or "",
            "tissue": tissue or "", "cell_type": cell_type or "", "perturbation": perturbation or "", "timepoint": timepoint or "",
        }
        bundle = api_post_file_binary(f"{TRIAGE_API_BASE}/analyze-bundle", input_path, data=form)
        bundle_path = run_dir / "triage_bundle.zip"
        bundle_path.write_bytes(bundle)
        unzip_dir = run_dir / "bundle"
        unzip_dir.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(bundle_path) as z:
            z.extractall(unzip_dir)
        store = {"run_dir": str(run_dir), "bundle_path": str(bundle_path), "unzip_dir": str(unzip_dir), "mode": mode}
        return status(f"Triage API run complete. mode={mode}", "ok"), store
    except Exception as e:
        return status(f"Triage API error: {e}", "bad"), {"run_dir": str(run_dir), "error": str(e)}

@app.callback(Output("triage-output-panel", "children"), Input("triage-output-tabs", "value"), State("triage-store", "data"))
def render_triage_output(tab, store):
    if not store:
        return html.Div("Run the Triage API first.", className="help")
    if store.get("error"):
        return html.Pre(store["error"], className="markdown-box")
    d = Path(store.get("unzip_dir", ""))
    if tab == "triage-report-md":
        p = d / "report.md"
        text = p.read_text(encoding="utf-8") if p.exists() else "No report.md found."
        return html.Pre(text, className="markdown-box")
    if tab == "triage-programs":
        p = d / "programs.csv"
        return datatable_from_df(pd.read_csv(p) if p.exists() else pd.DataFrame())
    if tab == "triage-rows":
        p = d / "triage_rows.csv"
        return datatable_from_df(pd.read_csv(p) if p.exists() else pd.DataFrame())
    if tab == "triage-claims":
        p = d / "claims.csv"
        return datatable_from_df(pd.read_csv(p) if p.exists() else pd.DataFrame())
    if tab == "triage-files":
        files = sorted([p.name for p in d.glob("*")]) if d.exists() else []
        return html.Div([html.Div(f, className="file-list") for f in files])
    p = d / "analysis_result.json"
    return html.Pre(p.read_text(encoding="utf-8")[:25000] if p.exists() else "No analysis_result.json found.", className="markdown-box")

@app.callback(Output("triage-download", "data"), Input("triage-download-bundle-btn", "n_clicks"), State("triage-store", "data"), prevent_initial_call=True)
def download_triage_bundle(n, store):
    if not store:
        return no_update
    return download_file(store.get("bundle_path", ""))

@app.callback(Output("triage-download", "data", allow_duplicate=True), Input("triage-download-md-btn", "n_clicks"), State("triage-store", "data"), prevent_initial_call=True)
def download_triage_md(n, store):
    if not store:
        return no_update
    return download_file(str(Path(store.get("unzip_dir", "")) / "report.md"))

@app.callback(Output("triage-download", "data", allow_duplicate=True), Input("triage-download-pdf-btn", "n_clicks"), State("triage-store", "data"), prevent_initial_call=True)
def download_triage_pdf(n, store):
    if not store:
        return no_update
    return download_file(str(Path(store.get("unzip_dir", "")) / "enrichment_llm_report.pdf"))

if __name__ == "__main__":
    app.run(host="127.0.0.1", port=int(os.environ.get("DASH_PORT", "8050")), debug=True)
