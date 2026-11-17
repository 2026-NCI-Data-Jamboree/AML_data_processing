# /// script
# requires-python = ">=3.10"
# dependencies = ["pandas>=2"]
# ///
"""Prepare cell type / cell state proportions of bulk RNA samples for the portal.

    results/ + data/  ->  deconvolution layers added to results/patients/<DATASET_ID>/ and results/cell_lines/
        values/deconv_<reference>_<method>_<level>.tsv.gz   rows = cell types or states, columns = sample key
        features/deconv_*.tsv, rows in datasets.tsv (modality deconvolution, unit proportion), manifest.tsv
        deconvolution_runs.tsv       one row per layer
        deconvolution_warnings.tsv   flags on runs, samples and labels (never hides a result)
      results/single_cell/<REFERENCE_ID>/{reference,markers,provenance}.tsv

The deconvolution tool is in deconvolution/ (R, own renv) and is launched as a separate Rscript process.
Layout and decisions: docs/results_contract.md (Deconvolution).

RUN LAST: prepare_cell_line_data.R and prepare_patient_data.R wipe their results/ subfolders.

    Rscript prepare_cell_line_data.R
    Rscript prepare_patient_data.R
    uv run prepare_deconvolution.py [--only DATASET_ID [DATASET_ID ...]]

Run from the repo root. Inputs: data/reference/aml-scatlas.rds, data/reference/hgnc_complete_set.txt and the
raw files named in SOURCES. The tool's renv must be restored first (deconvolution/README.md).

Sources done so far: WANG_2022_PED. To do: BEAT_AML, ZENODO_14160480, cell_lines (log inversion and
symbol -> Ensembl mapping with HGNC), MALANI_2022 (once harmonized).
"""

import argparse
import csv
import json
import os
import shutil
import subprocess
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Callable

import pandas as pd

# ---- Reference configuration: everything specific to the reference lives here ----------------------
# Swapping the reference means rebuilding the .rds and changing this block, then rerunning.
REFERENCE = {
    "id": "aml_scatlas",  # used in layer keys and in results/single_cell/<id>/
    "dataset_id": "AML_SCATLAS",  # row in dataset_sources.tsv
    "rds": Path("data/reference/aml-scatlas.rds"),
    "malignant": "Malignant",  # malignant cell type label in the reference
    "name": "AML scAtlas",
    "source": "AML scAtlas AnnData, https://doi.org/10.48420/27269946; "
    "built by 2026-nci-data-jamboree/analyses/aml-scatlas-reference",
    "citation": "Whittle et al., eLife 2026, https://doi.org/10.7554/eLife.104978",
    "limitations": (
        "Bone marrow, diagnosis samples only; Smart-seq2 libraries (Velten, Stetson) dropped because their "
        "counts are reads, not UMIs. The HSPC label is set only for cells whose author label is HSPC, so blasts "
        "in monocytic AML are labelled CD14+ Mono or ProMono: the normal myeloid profiles carry leukaemic signal "
        "and the malignant profile under-represents monocytic and more committed AML."
    ),
}

# ---- Settings ---------------------------------------------------------------------------------------
N_CORES = 2
N_MARKERS = 20  # markers per label in markers.tsv
MIN_BULK_SAMPLES = 10  # below this InstaPrism's filter lets one sample decide which genes are dropped
LOW_MALIGNANT = 0.1  # initial malignant proportion below this: the updated fit is unreliable
HEMOGLOBIN_GENES = ["HBB", "HBA1", "HBA2", "HBD"]
HEMOGLOBIN_SHARE = 0.10  # these genes above this share of a sample's signal
MIN_REFERENCE_GENES_PRESENT = 0.80  # share of the reference's genes present in the bulk matrix

HGNC_PATH = Path("data/reference/hgnc_complete_set.txt")
TOOL_DIR = Path("deconvolution")
RUNS_DIR = TOOL_DIR / "runs"  # bulk matrices, logs and tool outputs per run (gitignored)

# One run per cohort (and one for the cell lines). Each run: InstaPrism (initial + updated fit) and NNLS.
TOOL_METHODS = {
    "instaprism": ["--method", "instaprism", "--malignant-cell-type", REFERENCE["malignant"]],
    "nnls": ["--method", "nnls"],
}
METHOD_LABEL = {
    "instaprism": "InstaPrism initial fit",
    "instaprism_update_all": "InstaPrism updated fit",
    "instaprism_update_malignant": "InstaPrism updated fit (malignant profile only)",
    "instaprism_update_pooled": "InstaPrism updated fit (pooled)",
    "nnls": "NNLS",
}
UPDATE_METHODS = {"instaprism_update_all", "instaprism_update_malignant", "instaprism_update_pooled"}

WARNING_COLUMNS = ["key", "scope", "sample_id", "label", "severity", "code", "message", "source"]


# ---- Sources ----------------------------------------------------------------------------------------
@dataclass
class Bulk:
    """Linear bulk expression for the tool: genes x samples, columns = the folder's sample key."""

    matrix: pd.DataFrame  # rescaled to 1e6 per sample
    source: str  # input file, for deconvolution_runs.tsv
    conversion: str  # what was done to it, for deconvolution_runs.tsv
    tool_input: str  # the tool's --bulk-input: "counts" or "normalized"


def rescale(matrix: pd.DataFrame) -> pd.DataFrame:
    return matrix / matrix.sum() * 1e6


def read_wang(sample_ids: list[str]) -> Bulk:
    path = Path("data/patients/Blood_Cancer_Discovery_2022/GSE192638_Gene_count_matrix_updated.txt")
    table = pd.read_csv(path, sep="\t").set_index("gene_id")
    matrix = table[[s for s in sample_ids if s in table.columns]].astype(float)
    return Bulk(
        rescale(matrix),
        str(path),
        "none (linear TPM-like values rescaled to 1e6 per sample)",
        "normalized",
    )


@dataclass
class Source:
    dataset_id: str
    folder: Path
    key: str  # sample key column: sample_id (patients) or cell_line
    count_column: str  # datasets.tsv column holding the number of samples
    reader: Callable[[list[str]], Bulk]


SOURCES = [
    Source("WANG_2022_PED", Path("results/patients/WANG_2022_PED"), "sample_id", "n_samples", read_wang),
]


# ---- Helpers ----------------------------------------------------------------------------------------
def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, sep="\t", index=False)


def git_commit() -> str:
    try:
        return subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True, check=True
        ).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return ""


def run_tool(script: str, args: list, log: Path) -> None:
    """Run a tool script with Rscript from deconvolution/, so the tool activates its own renv."""
    env = {k: v for k, v in os.environ.items() if k not in ("RENV_PROJECT", "R_LIBS_USER")}
    log.parent.mkdir(parents=True, exist_ok=True)
    with open(log, "w") as f:
        result = subprocess.run(
            ["Rscript", script, *map(str, args)],
            cwd=TOOL_DIR, env=env, stdout=f, stderr=subprocess.STDOUT,
        )
    if result.returncode != 0:
        raise RuntimeError(f"{script} failed; see {log}")


def read_hgnc() -> tuple[pd.DataFrame, str]:
    """Ensembl ID <-> approved symbol, and the newest record date of the snapshot."""
    hgnc = pd.read_csv(
        HGNC_PATH, sep="\t", dtype=str, keep_default_na=False, quoting=csv.QUOTE_NONE,
        usecols=["symbol", "ensembl_gene_id", "date_modified"],
    )
    snapshot = hgnc["date_modified"].max()
    return hgnc[hgnc["ensembl_gene_id"] != ""], snapshot


def warning_rows(severity: str, code: str, message, source: str, key="", scope="run", sample_id="") -> pd.DataFrame:
    """Rows for deconvolution_warnings.tsv; key, sample_id and message may be scalars or lists."""
    return pd.DataFrame(
        {"key": key, "scope": scope, "sample_id": sample_id, "label": "", "severity": severity,
         "code": code, "message": message, "source": source},
        columns=WARNING_COLUMNS,
    )


def classify_tool_warning(message: str) -> tuple[str, str]:
    """A stable code and a severity for a warning raised by the tool."""
    if message.startswith("Initial '") and "proportion is below" in message:
        return "low_malignant_initial", "caution"
    if "not all integers" in message:
        return "non_integer_counts", "info"
    if "appear more than once" in message:
        return "duplicate_gene_ids", "info"
    if "overlap" in message.lower():
        return "low_gene_overlap", "warning"
    if "min-cells" in message or "dropp" in message.lower():
        return "labels_dropped", "info"
    return "tool_warning", "warning"


# ---- One source -------------------------------------------------------------------------------------
def run_methods(source: Source, bulk: Bulk) -> tuple[pd.DataFrame, dict[str, dict], dict[str, str]]:
    """Run the tool once per method. Returns the long proportions, run_info per tool run, and which
    tool run produced each method (InstaPrism's run also yields the updated fit)."""
    run_dir = RUNS_DIR / source.dataset_id
    run_dir.mkdir(parents=True, exist_ok=True)
    bulk_path = (run_dir / "bulk.tsv").resolve()
    bulk.matrix.rename_axis("gene_id").to_csv(bulk_path, sep="\t")

    proportions, run_infos, run_of_method = [], {}, {}
    for run, method_args in TOOL_METHODS.items():
        print(f"running {run}")
        prefix = f"{source.dataset_id}_{run}"
        run_tool(
            "deconvolve.R",
            ["--reference", REFERENCE["rds"].resolve(), "--bulk", bulk_path, "--bulk-input", bulk.tool_input,
             *method_args, "--output-dir", run_dir.resolve(), "--output-prefix", prefix, "--n-cores", N_CORES],
            run_dir.resolve() / f"{run}.log",
        )
        long = pd.read_csv(run_dir / f"{prefix}_proportions.tsv", sep="\t")
        proportions.append(long)
        run_infos[run] = json.loads((run_dir / f"{prefix}_run_info.json").read_text())
        run_of_method.update({method: run for method in long["method"].unique()})
    return pd.concat(proportions), run_infos, run_of_method


def write_layers(source: Source, props: pd.DataFrame, run_infos, run_of_method, bulk: Bulk,
                 sample_ids: list[str], commit: str, hgnc_snapshot: str) -> dict[str, tuple[str, str]]:
    """Write the layer files and the runs table. Returns {layer key: (method, level)}."""
    folder = source.folder
    for old in [*(folder / "values").glob("deconv_*"), *(folder / "features").glob("deconv_*")]:
        old.unlink()

    datasets = read_tsv(folder / "datasets.tsv")
    datasets = datasets[~datasets["key"].str.startswith("deconv_")]
    manifest = read_tsv(folder / "manifest.tsv")
    manifest = manifest[~manifest["path"].str.contains("deconv")]

    layers, dataset_rows, manifest_rows, run_rows = {}, [], [], []
    for (method, level, quantity), layer in props.groupby(["method", "level", "quantity"], sort=False):
        key = f"deconv_{REFERENCE['id']}_{method}_{level}"
        wide = layer.pivot_table(index="label", columns="sample_id", values="proportion", aggfunc="sum")
        wide = wide.sort_index()[[s for s in sample_ids if s in wide.columns]]
        wide.round(6).rename_axis("feature").to_csv(folder / "values" / f"{key}.tsv.gz", sep="\t", compression="gzip")
        features = pd.DataFrame({"feature": wide.index, "label": wide.index, "gene_symbol": "", "entrez_id": ""})
        write_tsv(features, folder / "features" / f"{key}.tsv")

        level_name = "cell types" if level == "cell_type" else "cell states"
        label = f"{METHOD_LABEL[method]}, {level_name} ({REFERENCE['name']})"
        dataset_rows.append({
            **{column: "" for column in datasets.columns},
            "key": key, "label": label, "modality": "deconvolution", "unit": "proportion",
            "source": f"{REFERENCE['name']} reference, {METHOD_LABEL[method].split(' (')[0]}",
            "n_features": len(wide), source.count_column: wide.shape[1],
        })
        manifest_rows += [
            {"path": f"features/{key}.tsv", "kind": "features", "description": f"{label} - feature annotation"},
            {"path": f"values/{key}.tsv.gz", "kind": "values", "description": f"{label} - values per sample"},
        ]
        info = run_infos[run_of_method[method]]
        is_instaprism = method.startswith("instaprism")
        run_rows.append({
            "key": key, "reference": REFERENCE["id"], "method": method, "level": level, "quantity": quantity,
            "instaprism_filter": str(info["method"]["options"]["instaprism_filter"]).lower() if is_instaprism else "",
            "update_mode": info["method"]["options"]["update_mode"] if method in UPDATE_METHODS else "",
            "min_cells": info["parameters"]["min_cells"],
            "malignant_cell_type": REFERENCE["malignant"] if is_instaprism else "",
            "n_samples": wide.shape[1], "input_source": bulk.source, "input_conversion": bulk.conversion,
            "tool_version": commit, "run_date": date.today().isoformat(), "hgnc_snapshot": hgnc_snapshot,
        })
        layers[key] = (method, level)

    manifest_rows += [
        {"path": "deconvolution_runs.tsv", "kind": "deconvolution",
         "description": "Deconvolution runs, one row per layer"},
        {"path": "deconvolution_warnings.tsv", "kind": "deconvolution",
         "description": "Deconvolution warnings on runs, samples and labels"},
    ]
    write_tsv(pd.concat([datasets, pd.DataFrame(dataset_rows)]), folder / "datasets.tsv")
    write_tsv(pd.concat([manifest, pd.DataFrame(manifest_rows)]), folder / "manifest.tsv")
    write_tsv(pd.DataFrame(run_rows), folder / "deconvolution_runs.tsv")
    return layers


def collect_warnings(bulk: Bulk, props: pd.DataFrame, run_infos, run_of_method, layers, hgnc) -> pd.DataFrame:
    keys_of_run = {run: [k for k, (m, _) in layers.items() if run_of_method[m] == run] for run in run_infos}
    update_keys = [k for k, (m, _) in layers.items() if m in UPDATE_METHODS]
    rows = []

    # warnings raised by the tool, flagged on the layers of the run that raised them
    for run, info in run_infos.items():
        for entry in info["warnings"]:
            code, severity = classify_tool_warning(entry["message"])
            # the low-malignant warning is about the reference update, so it flags the updated layers only
            keys = [k for k in keys_of_run[run] if k in update_keys] if code == "low_malignant_initial" else keys_of_run[run]
            if keys:
                rows.append(warning_rows(severity, code, entry["message"], "tool", key=keys))

    if update_keys:
        rows.append(warning_rows(
            "caution", "updated_fit_batch_dependent", "The updated fit re-estimates the normal cell type profiles "
            "pooled over all samples in this run, so a sample's estimate depends on the other samples; it also uses "
            f"the {REFERENCE['name']} malignant profile, which under-represents monocytic and more committed AML. "
            "Prefer the initial fit unless the update is needed.", "processing", key=update_keys))
        initial = props[(props["method"] == "instaprism") & (props["level"] == "cell_type")
                        & (props["label"] == REFERENCE["malignant"])]
        low = initial[initial["proportion"] < LOW_MALIGNANT]
        for key in update_keys:
            rows.append(warning_rows(
                "caution", "low_malignant_initial",
                [f"Initial {REFERENCE['malignant']} proportion {p:.3f} is below {LOW_MALIGNANT:.2f}: the updated fit "
                 "may let the malignant profile absorb normal cell types" for p in low["proportion"]],
                "processing", key=key, scope="sample", sample_id=list(low["sample_id"])))

    n_samples = bulk.matrix.shape[1]
    if n_samples < MIN_BULK_SAMPLES:
        rows.append(warning_rows(
            "warning", "small_bulk_set",
            f"Only {n_samples} samples in this run (below {MIN_BULK_SAMPLES}): InstaPrism's gene filter then "
            "lets a single sample decide which genes are dropped", "processing"))

    hemoglobin_ids = hgnc.loc[hgnc["symbol"].isin(HEMOGLOBIN_GENES), "ensembl_gene_id"]
    share = bulk.matrix[bulk.matrix.index.isin(hemoglobin_ids)].sum() / bulk.matrix.sum()
    share = share[share > HEMOGLOBIN_SHARE]
    if len(share):
        rows.append(warning_rows(
            "caution", "hemoglobin_dominated",
            [f"{', '.join(HEMOGLOBIN_GENES)} make up {100 * s:.0f}% of this sample's signal (above "
             f"{100 * HEMOGLOBIN_SHARE:.0f}%); red cell contamination can distort the estimates" for s in share],
            "processing", scope="sample", sample_id=list(share.index)))

    funnel = next(iter(run_infos.values()))["genes"]["funnel"]
    present = funnel["shared"] / funnel["reference"]
    if present < MIN_REFERENCE_GENES_PRESENT:
        rows.append(warning_rows(
            "warning", "low_reference_gene_coverage",
            f"Only {100 * present:.0f}% of the reference's genes ({funnel['shared']} of {funnel['reference']}) are "
            f"present in the bulk matrix (below {100 * MIN_REFERENCE_GENES_PRESENT:.0f}%)", "processing"))

    if not rows:
        return pd.DataFrame(columns=WARNING_COLUMNS)
    return pd.concat(rows).sort_values(["key", "scope", "sample_id"], kind="stable").drop_duplicates()


def process_source(source: Source, hgnc: pd.DataFrame, hgnc_snapshot: str, commit: str) -> None:
    print(f"== {source.dataset_id}")
    samples = read_tsv(source.folder / "samples.tsv")
    sample_ids = list(samples[source.key])
    bulk = source.reader(sample_ids)
    print(f"{bulk.matrix.shape[1]} of {len(sample_ids)} samples have RNA; {bulk.matrix.shape[0]} genes")

    props, run_infos, run_of_method = run_methods(source, bulk)
    layers = write_layers(source, props, run_infos, run_of_method, bulk, sample_ids, commit, hgnc_snapshot)
    warnings = collect_warnings(bulk, props, run_infos, run_of_method, layers, hgnc)
    write_tsv(warnings, source.folder / "deconvolution_warnings.tsv")
    print(f"{len(layers)} layers; {len(warnings)} warning rows")


# ---- Reference summary ------------------------------------------------------------------------------
def write_reference_files(hgnc: pd.DataFrame, hgnc_snapshot: str, commit: str) -> None:
    """results/single_cell/<REFERENCE_ID>/: made by the tool's summarize_reference.R plus HGNC symbols."""
    print(f"== reference {REFERENCE['id']}")
    summary = (RUNS_DIR / "reference_summary").resolve()
    run_tool(
        "summarize_reference.R",
        ["--reference", REFERENCE["rds"].resolve(), "--malignant-cell-type", REFERENCE["malignant"],
         "--n-markers", N_MARKERS, "--output-dir", summary],
        RUNS_DIR.resolve() / "reference_summary.log",
    )
    out = Path("results/single_cell") / REFERENCE["id"]
    out.mkdir(parents=True, exist_ok=True)
    shutil.copy(summary / "reference.tsv", out / "reference.tsv")

    markers = pd.read_csv(summary / "markers.tsv", sep="\t", dtype={"gene_id": str})
    symbol_of = hgnc.drop_duplicates("ensembl_gene_id").set_index("ensembl_gene_id")["symbol"]  # first record wins
    markers["gene_symbol"] = markers["gene_id"].map(symbol_of).fillna("")
    write_tsv(markers[["level", "label", "rank", "gene_id", "gene_symbol", "score"]], out / "markers.tsv")

    build = read_tsv(summary / "build.tsv").set_index("key")["value"]
    provenance = {
        "reference_id": REFERENCE["id"], "reference_name": REFERENCE["name"], "source": REFERENCE["source"],
        "citation": REFERENCE["citation"], "build_date": build["build_date"], "processing_commit": commit,
        "hgnc_snapshot": hgnc_snapshot, "cell_type_column": build["cell_type_column"],
        "cell_state_column": build["cell_state_column"], "sample_column": build["sample_column"],
        "unlabelled": build["unlabelled"], "excluded_labels": build["exclude_labels"],
        "n_input_files": build["n_input_files"], "min_cells": build["min_cells"],
        "malignant_cell_type": REFERENCE["malignant"], "limitations": REFERENCE["limitations"],
    }
    write_tsv(pd.DataFrame({"key": provenance.keys(), "value": provenance.values()}), out / "provenance.tsv")
    print(f"wrote {out}")


# ---- Main -------------------------------------------------------------------------------------------
def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--only", nargs="+", metavar="DATASET_ID", help="rerun only these sources")
    args = parser.parse_args()
    known = {s.dataset_id for s in SOURCES}
    if args.only and set(args.only) - known:
        parser.error(f"unknown or not yet supported source(s): {', '.join(sorted(set(args.only) - known))}")

    if not Path("results").is_dir():
        raise SystemExit("results/ is missing; run prepare_cell_line_data.R and prepare_patient_data.R first")
    for path, hint in [(REFERENCE["rds"], "data/README.md"), (HGNC_PATH, "data/README.md"),
                       (TOOL_DIR / "renv" / "library", f"restore the tool's renv: {TOOL_DIR}/README.md")]:
        if not path.exists():
            raise SystemExit(f"{path} not found ({hint})")

    hgnc, hgnc_snapshot = read_hgnc()
    commit = git_commit()
    for source in SOURCES:
        if not args.only or source.dataset_id in args.only:
            process_source(source, hgnc, hgnc_snapshot, commit)
    write_reference_files(hgnc, hgnc_snapshot, commit)
    shutil.copy("dataset_sources.tsv", "results")
    print("Done.")


if __name__ == "__main__":
    main()
