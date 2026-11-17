# /// script
# requires-python = ">=3.10"
# dependencies = ["pandas>=2"]
# ///
"""Prepare cell type / cell state proportions of bulk RNA samples for the portal.

    results/ + data/  ->  deconvolution layers added to results/patients/<DATASET_ID>/
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

INPUT: raw or estimated counts, or linear CPM (library-size normalized only). Transcript-length-normalized data
(TPM, RPKM, FPKM) and log-scale values must not be deconvolved: the reference is built from counts, and dividing by
transcript length changes the relative weight of genes. A source is judged by the unit recorded for its data
(datasets.tsv `unit`, or the file's documentation); read_log_expression() refuses TPM, RPKM and FPKM and any unit
that is not in CPM_UNITS. Zenodo 14160480 (log2(TPM+1)) is therefore not a source.

Sources:
  raw counts   BEAT_AML, TCGA_LAML, TARGET_AML, LEUCEGENE, LUMC_AML (AMLmap's uncorrected_counts.csv; the batch-
               corrected counts only feed the RNA layer) and MALANI_2022 (File_8 raw read counts)
  log2 CPM     WANG_2022_PED (Table S13, CPM per the paper; log2(CPM + 0.001)) and AML_STJUDE_2024 (voom log2(CPM),
               then ComBat-corrected: approximate, flagged `caution`), inverted from the harmonized layer
Not sources: ZENODO_14160480 (TPM), LEE_2018 (microarray).
Cell lines are not deconvolved: a line is one cell type, so proportions say nothing.
"""

import argparse
import csv
import json
import re
import os
import shutil
import subprocess
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path
from typing import Callable

import numpy as np
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

# One run per cohort. Each run: InstaPrism (initial + updated fit) and NNLS.
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
    """Bulk RNA for the tool: genes (Ensembl IDs) x samples, columns = the folder's sample key.

    Raw or estimated counts, or linear CPM-like values (library-size normalized only). Transcript-length-normalized
    data (TPM, RPKM, FPKM) and log-scale values must not be deconvolved, so a cohort whose only RNA is such data
    is not deconvolved: there is no reader for it.
    """

    matrix: pd.DataFrame
    source: str  # input file, for deconvolution_runs.tsv
    conversion: str  # what was done to it, for deconvolution_runs.tsv
    tool_input: str  # the tool's --bulk-input: "counts" or "normalized"
    notes: list[tuple[str, str, str]] = field(default_factory=list)  # (severity, code, message) for the warnings table


def rescale(matrix: pd.DataFrame) -> pd.DataFrame:
    return matrix / matrix.sum() * 1e6


AMLMAP_COUNTS = Path("data/patients/AMLmap/counts/uncorrected_counts.csv")
AMLMAP_META = Path("data/patients/AMLmap/meta.csv")


def read_amlmap(study: str, amlmap_id: Callable[[str], str] = lambda sample_id: sample_id) -> Callable:
    """Reader for one AMLmap study: raw integer counts of the pipeline that re-processed all five studies
    (`uncorrected_counts.csv`, not the batch-corrected files). `amlmap_id` maps a folder's sample ID to AMLmap's."""

    def read(sample_ids: list[str], hgnc: pd.DataFrame) -> Bulk:
        meta = pd.read_csv(AMLMAP_META, usecols=["sample_id", "study"])
        in_study = set(meta.loc[meta["study"] == study, "sample_id"])
        wanted = {sid: amlmap_id(sid) for sid in sample_ids if amlmap_id(sid) in in_study}
        header = pd.read_csv(AMLMAP_COUNTS, nrows=0).columns
        columns = [header[0], *[c for c in header[1:] if c in set(wanted.values())]]
        counts = pd.read_csv(AMLMAP_COUNTS, index_col=0, usecols=columns).astype(float)
        counts = counts.rename(columns={v: k for k, v in wanted.items()})
        counts.index = counts.index.str.replace(r"\.\d+$", "", regex=True)  # drop Ensembl version suffixes
        counts = counts.groupby(level=0).sum()
        return Bulk(counts, str(AMLMAP_COUNTS), f"none (raw integer counts, AMLmap {study})", "counts")

    return read


MALANI_COUNTS = Path("data/patients/Functional_Precision_Medicine_Tumor_Board_AML/File_8_RNA_seq_Raw_Reads_163S_4Healthy.csv")


def read_malani(sample_ids: list[str], hgnc: pd.DataFrame) -> Bulk:
    """Malani 2022 raw read counts (163 AML samples, 4 healthy CD34+ controls that are not in the cohort)."""
    header = pd.read_csv(MALANI_COUNTS, nrows=0).columns  # first column: Ensembl ID, then one column per sample
    wanted = [s for s in sample_ids if s in header]
    counts = pd.read_csv(MALANI_COUNTS, index_col=0, usecols=[header[0], *wanted]).astype(float)
    return Bulk(counts, str(MALANI_COUNTS), "none (raw integer counts)", "counts")


# The harmonized units that can be deconvolved after inverting to linear CPM: library-size-normalized only. The
# pseudocount is not in the unit string, so each cohort states its own (see the sources below).
CPM_UNITS = {"log2 CPM"}
LENGTH_NORMALIZED = re.compile(r"\b(TPM|RPKM|FPKM)\b", re.IGNORECASE)
GENE_LOSS_WARNING = 0.05  # warning (not info) if more than this share of genes cannot be mapped


def read_log_expression(folder: Path, pseudocount: float, caveat: str = "") -> Callable:
    """Reader for a folder's harmonized expression layer when its recorded unit is a log2 CPM: invert it with
    `2^x - pseudocount`, map the symbol-keyed genes to Ensembl IDs, sum genes that share an ID and rescale to
    1e6 per sample. Refuses TPM, RPKM and FPKM, and any unit that is not in CPM_UNITS. `caveat` is added to the
    `caution` warning when the layer was processed further after the CPM step."""

    def read(sample_ids: list[str], hgnc: pd.DataFrame) -> Bulk:
        unit = read_tsv(folder / "datasets.tsv").set_index("key").loc["expression", "unit"]
        if LENGTH_NORMALIZED.search(unit):
            raise ValueError(f"{folder}: the expression unit '{unit}' is transcript-length-normalized and must not be "
                             "deconvolved; use counts or CPM")
        if unit not in CPM_UNITS:
            raise ValueError(f"{folder}: no conversion to linear CPM for the expression unit '{unit}'")
        features = read_tsv(folder / "features" / "expression.tsv").set_index("feature")
        values = pd.read_csv(folder / "values" / "expression.tsv.gz", sep="\t", index_col=0)
        values = values[[s for s in sample_ids if s in values.columns]]
        linear = (np.exp2(values.astype(float)) - pseudocount).clip(lower=0)

        mapped = symbols_to_ensembl(features.loc[linear.index, "gene_symbol"], hgnc)
        resolved = mapped["ensembl"] != ""
        by_ensembl = linear[resolved].groupby(mapped.loc[resolved, "ensembl"]).sum()  # shared IDs are summed
        steps = mapped.loc[resolved, "step"].value_counts()
        lost = (~resolved).sum() / len(resolved)
        notes = [
            ("caution", "log_scale_input",
             f"The input is on the {unit} scale and was converted to linear CPM before deconvolution; the "
             "conversion is only as exact as the recorded unit and pseudocount, and the layer's gene filtering "
             "(protein-coding genes) means the rescaled values are not the original CPM. " + caveat),
            ("warning" if lost > GENE_LOSS_WARNING else "info", "gene_mapping_loss",
             f"{(~resolved).sum()} of {len(resolved)} genes ({100 * lost:.1f}%) could not be mapped to a single Ensembl "
             f"ID with HGNC and were dropped (resolved by: {', '.join(f'{n} {s}' for s, n in steps.items())}); "
             f"{resolved.sum() - len(by_ensembl)} genes shared an ID and were summed."),
        ]
        return Bulk(
            rescale(by_ensembl),
            f"{folder}/values/expression.tsv.gz",
            f"inverted from {unit} with 2^x - {pseudocount:g}; symbols mapped to Ensembl IDs with HGNC; rescaled to 1e6 per sample",
            "normalized",
            notes,
        )

    return read


@dataclass
class Source:
    dataset_id: str
    folder: Path
    key: str  # sample key column: sample_id (patients) or cell_line
    count_column: str  # datasets.tsv column holding the number of samples
    reader: Callable[[list[str], pd.DataFrame], Bulk]  # (sample ids, HGNC table) -> linear bulk matrix


SOURCES = [
    # Wang's harmonized layer is log2 CPM; its raw GEO matrix is TPM-like and is not used
    Source("WANG_2022_PED", Path("results/patients/WANG_2022_PED"), "sample_id", "n_samples",
           read_log_expression(Path("results/patients/WANG_2022_PED"), pseudocount=0.001)),
    # St. Jude: voom log2(CPM) (Umeda et al. 2024: "The count data were transformed to log2(cpm) using Voom"), then
    # corrected with ComBat for library strand and for the TARGET batch. Voom's own prior is ~0.01 CPM, so none is
    # subtracted; the ComBat step cannot be undone without its parameters, hence the caveat.
    Source("AML_STJUDE_2024", Path("results/patients/AML_STJUDE_2024"), "sample_id", "n_samples",
           read_log_expression(
               Path("results/patients/AML_STJUDE_2024"), pseudocount=0.0,
               caveat="St. Jude's values are voom log2(CPM) after ComBat correction for library strand and the TARGET "
                      "batch (Umeda et al. 2024), so they are only approximately CPM (values reach -13.6, below voom's "
                      "floor); the batch correction cannot be undone.")),
    Source("MALANI_2022", Path("results/patients/MALANI_2022"), "sample_id", "n_samples", read_malani),
    # Beat AML RNA samples end in "R" in the folder and have no suffix in AMLmap
    Source("BEAT_AML", Path("results/patients/BEAT_AML"), "sample_id", "n_samples",
           read_amlmap("BEAT_AML", lambda sid: sid.removesuffix("R") if sid.endswith("R") else "")),
    # the four other AMLmap studies are cohorts of prepare_patient_data.R (dataset ID <- AMLmap study name)
    Source("TCGA_LAML", Path("results/patients/TCGA_LAML"), "sample_id", "n_samples", read_amlmap("TCGA_AML")),
    Source("TARGET_AML", Path("results/patients/TARGET_AML"), "sample_id", "n_samples", read_amlmap("TARGET_AML")),
    Source("LEUCEGENE", Path("results/patients/LEUCEGENE"), "sample_id", "n_samples", read_amlmap("LEUCEGENE_AML")),
    Source("LUMC_AML", Path("results/patients/LUMC_AML"), "sample_id", "n_samples", read_amlmap("100LUMC_AML")),
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
    """The approved genes with an Ensembl ID, and the newest record date of the snapshot."""
    hgnc = pd.read_csv(
        HGNC_PATH, sep="\t", dtype=str, keep_default_na=False, quoting=csv.QUOTE_NONE,
        usecols=["symbol", "alias_symbol", "prev_symbol", "ensembl_gene_id", "date_modified"],
    )
    snapshot = hgnc["date_modified"].max()
    for column in ["alias_symbol", "prev_symbol"]:
        hgnc[column] = hgnc[column].str.strip('"')
    return hgnc[hgnc["ensembl_gene_id"] != ""], snapshot


def symbols_to_ensembl(symbols: pd.Series, hgnc: pd.DataFrame) -> pd.DataFrame:
    """Map gene symbols to Ensembl IDs, stopping at the first step that resolves a gene:
    approved symbol; then a previous symbol or an alias that points to one gene only.
    Returns a frame indexed like `symbols` with `ensembl` (empty if unresolved or ambiguous) and `step`."""

    def unique_targets(column: str) -> pd.Series:
        pairs = hgnc.assign(term=hgnc[column].str.split("|")).explode("term")
        pairs = pairs[pairs["term"] != ""].drop_duplicates(["term", "ensembl_gene_id"])
        return pairs.groupby("term")["ensembl_gene_id"].agg(lambda ids: ids.iloc[0] if len(ids) == 1 else "")

    lookups = [("approved symbol", hgnc.drop_duplicates("symbol").set_index("symbol")["ensembl_gene_id"]),
               ("previous symbol", unique_targets("prev_symbol")), ("alias", unique_targets("alias_symbol"))]
    result = pd.DataFrame({"ensembl": "", "step": ""}, index=symbols.index)
    for step, table in lookups:
        found = symbols.map(table).fillna("")
        new = (result["ensembl"] == "") & (found != "")
        result.loc[new, "ensembl"] = found[new]
        result.loc[new, "step"] = step
    return result


def warning_rows(severity: str, code: str, message, source: str, key="", scope="run", sample_id="") -> pd.DataFrame:
    """Rows for deconvolution_warnings.tsv; key, sample_id and message may be scalars or lists."""
    fields = {"key": key, "scope": scope, "sample_id": sample_id, "label": "", "severity": severity,
              "code": code, "message": message, "source": source}
    n_rows = max((len(v) for v in fields.values() if isinstance(v, list)), default=1)
    return pd.DataFrame(fields, index=range(n_rows), columns=WARNING_COLUMNS)


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

    for severity, code, message in bulk.notes:
        rows.append(warning_rows(severity, code, message, "processing"))

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
    bulk = source.reader(sample_ids, hgnc)
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
