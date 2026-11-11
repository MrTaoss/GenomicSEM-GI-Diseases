#!/usr/bin/env python3
"""Core spatial genetic enrichment workflow for gastrointestinal Genomic SEM factors.

Companion to: Resolving the latent genetic architecture and spatial enrichment
patterns of common gastrointestinal diseases.

Spatial transcriptomics: Harnik et al., Nature (2024), eight adult proximal
small-intestine Visium sections (P1, P5, P8, P9, P10, P11, P12, P13).
Anatomical labeling follows the authors' MATLAB implementation:
https://github.com/yotamharnik/ST_Human_SI
Reference commit: cbccfecca66b7896e0595905ad3c01f4e6484ecf
Spatial enrichment: gsMap 1.73.8 (Song et al., Nature, 2025).

Set PROJECT_ROOT and the three factor GWAS files before running. Input source
files use the original Harnik repository naming. The GSMap stage consumes the
three final, Q_SNP/local-LD-filtered factor GWAS summary-statistics files.
Run stages in order: --stage prepare, --stage gsmap, --stage anatomy.
The gsmap stage uses official gsMap commands, not reimplemented algorithms.
"""

from __future__ import annotations

import argparse
import gc
import json
import math
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import anndata as ad
import numpy as np
import pandas as pd
from PIL import Image
from scipy import sparse
from scipy.spatial.distance import cdist


# Paths for prepared input data and gsMap reference resources.
PROJECT_ROOT = Path("path/to/gsmap_project")
RAW_DIR = PROJECT_ROOT / "Data" / "RAW"
H5AD_DIR = PROJECT_ROOT / "gsMap_example_data" / "ST"
WORKDIR = PROJECT_ROOT / "workdir"
ANATOMY_WORKDIR = PROJECT_ROOT / "workdir_anatomy"
ANNOTATION_DIR = PROJECT_ROOT / "anatomy_annotation"
RESOURCE_DIR = PROJECT_ROOT / "gsMap_resource"
GWAS_DIR = Path("path/to/final_factor_GWAS_sumstats")
SCRIPT_VERSION = "1.0.0"

SAMPLES = ("P1", "P5", "P8", "P9", "P10", "P11", "P12", "P13")
TRAITS = ("FUPPER", "FBOWEL", "FIBD")
PAPER_IDS = {
    "P1": "P1", "P5": "P2", "P8": "P3", "P9": "P4",
    "P10": "P5", "P11": "P6", "P12": "P7", "P13": "P8",
}
EXPECTED_TOTAL_RETAINED = 20_971
EXPECTED_TOTAL_LABELED = 17_960
EXPECTED_RETAINED_BY_SAMPLE = {
    "P1": 2281, "P5": 3018, "P8": 3105, "P9": 3110,
    "P10": 3741, "P11": 2577, "P12": 1581, "P13": 1558,
}
EXPECTED_LABELED_BY_SAMPLE = {
    "P1": 1950, "P5": 2677, "P8": 2788, "P9": 2662,
    "P10": 2707, "P11": 2395, "P12": 1332, "P13": 1449,
}
LOW_UMI_Z_THRESHOLD = -4.0
MITO_Z_THRESHOLD = 4.0
PIXEL_COLUMNS = ("Apical", "Basel", "Musmucosa", "Lumen", "Stroma", "Submucosa")
POSITION_COLUMNS = (
    "barcode", "in_tissue", "array_row", "array_col",
    "pxl_row_in_fullres", "pxl_col_in_fullres",
)

N_BINS = 7
VILLI_LEN_THRESH = 3
WITHIN_VILLUS_NEIGHBOUR_UM = 110.0
P10_MUSCULARIS_THRESHOLD_HIRES_PIXELS = 50.0
VILLUS_LENGTH_PASSES = ((0.0, 400.0), (0.0, math.inf))
REGION_BY_CODE = {
    -2: "SubMucosa", -1: "MusMucosa", 0: "Crypt",
    1: "V1", 2: "V2", 3: "V3", 4: "V4", 5: "V5", 6: "V6",
}
REGION_ORDER = [REGION_BY_CODE[code] for code in sorted(REGION_BY_CODE)]


class BuildError(RuntimeError):
    pass


class AnatomyError(RuntimeError):
    pass


def read_count_header(path: Path) -> tuple[str, list[str]]:
    columns = pd.read_csv(path, sep=";", nrows=0).columns.astype(str).tolist()
    if len(columns) < 2:
        raise BuildError(f"Count matrix has fewer than two columns: {path}")
    return columns[0], columns[1:]


def read_counts(path: Path, selected_barcodes: list[str] | None = None) -> pd.DataFrame:
    gene_column, all_barcodes = read_count_header(path)
    if selected_barcodes is None:
        selected_barcodes = all_barcodes
    missing = sorted(set(selected_barcodes).difference(all_barcodes))
    if missing:
        raise BuildError(f"Requested count columns are absent from {path}: {missing[:5]}")

    usecols = [gene_column, *selected_barcodes]
    dtype_map = {barcode: np.int32 for barcode in selected_barcodes}
    frame = pd.read_csv(
        path,
        sep=";",
        usecols=usecols,
        dtype=dtype_map,
        index_col=gene_column,
        engine="c",
        memory_map=True,
    )
    frame.index = frame.index.astype(str).str.strip()
    frame.columns = frame.columns.astype(str)
    if frame.index.hasnans or (frame.index == "").any():
        raise BuildError(f"Missing gene symbols in {path}")
    if frame.columns.has_duplicates:
        raise BuildError(f"Duplicated spot barcodes in {path}")
    values = frame.to_numpy(copy=False)
    if not np.issubdtype(values.dtype, np.integer) or np.any(values < 0):
        raise BuildError(f"Counts must be non-negative integers: {path}")
    return frame


def read_positions(path: Path) -> pd.DataFrame:
    positions = pd.read_csv(path, header=None, names=POSITION_COLUMNS, dtype={"barcode": str})
    positions["barcode"] = positions["barcode"].str.strip()
    if positions["barcode"].duplicated().any():
        raise BuildError(f"Duplicated barcodes in positions file: {path}")
    positions = positions.set_index("barcode", drop=False)
    return positions


def read_manual_exclusions(path: Path) -> set[str]:
    frame = pd.read_csv(path, dtype=str)
    if "Barcode" not in frame.columns:
        raise BuildError(f"Manual exclusion file lacks 'Barcode': {path}")
    barcodes = frame["Barcode"].dropna().astype(str).str.strip()
    return set(barcodes[barcodes != ""])


def read_pixel_metadata(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    required = [*PIXEL_COLUMNS, "barcodes"]
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise BuildError(f"Pixel-classifier file lacks columns: {missing}")
    if frame["barcodes"].duplicated().any():
        raise BuildError("Pixel-classifier file contains duplicated patient-prefixed barcodes")
    frame = frame.set_index("barcodes")
    for column in PIXEL_COLUMNS:
        frame[column] = pd.to_numeric(frame[column], errors="raise").astype(np.float32)
        if ((frame[column] < 0) | (frame[column] > 1)).any():
            raise BuildError(f"Pixel fraction outside [0,1] in column {column}")
    return frame


def read_metadata(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, encoding="utf-8-sig")
    required = {"Patient.Code.orig", "Patient.Code.paper"}
    if not required.issubset(frame.columns):
        raise BuildError(f"Metadata lacks required columns: {sorted(required)}")
    frame["Patient.Code.orig"] = frame["Patient.Code.orig"].astype(str)
    return frame.set_index("Patient.Code.orig", drop=False)


def zscore_matlab(values: np.ndarray, label: str) -> np.ndarray:
    """MATLAB-compatible sample standard deviation (normalization by N-1)."""
    values = np.asarray(values, dtype=np.float64)
    standard_deviation = values.std(ddof=1)
    if not np.isfinite(standard_deviation) or standard_deviation == 0:
        raise BuildError(f"Cannot compute a finite, non-zero standard deviation for {label}")
    return (values - values.mean()) / standard_deviation


def calculate_background_corrected_qc(
    utt_counts: np.ndarray,
    background_mean: np.ndarray,
    genes: np.ndarray,
    chunk_rows: int = 512,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Reproduce author QC without materializing a full corrected matrix."""
    n_genes, n_spots = utt_counts.shape
    if background_mean.shape != (n_genes,):
        raise BuildError("Background vector is not aligned with the UTT gene rows")

    total = np.zeros(n_spots, dtype=np.float64)
    mito_total = np.zeros(n_spots, dtype=np.float64)
    mito_mask = np.fromiter(("mt-" in gene.lower() for gene in genes), dtype=bool, count=n_genes)
    if not mito_mask.any():
        raise BuildError("No mitochondrial genes matched the authors' case-insensitive 'mt-' rule")

    for start in range(0, n_genes, chunk_rows):
        end = min(start + chunk_rows, n_genes)
        block = utt_counts[start:end, :].astype(np.float64, copy=True)
        block -= background_mean[start:end, None]
        np.maximum(block, 0.0, out=block)
        total += block.sum(axis=0)
        local_mito = mito_mask[start:end]
        if local_mito.any():
            mito_total += block[local_mito, :].sum(axis=0)

    if np.any(total <= 0) or not np.isfinite(total).all():
        raise BuildError("Background-corrected total UMI is non-positive or non-finite")
    log_total = np.log10(total)
    mito_fraction = mito_total / total
    z_log_total = zscore_matlab(log_total, "log10 corrected total UMI")
    z_mito_fraction = zscore_matlab(mito_fraction, "corrected mitochondrial fraction")
    return total, mito_fraction, z_log_total, z_mito_fraction


def duplicate_plan(genes: np.ndarray) -> tuple[np.ndarray, dict[int, int], pd.DataFrame]:
    first_position: dict[str, int] = {}
    keep: list[int] = []
    duplicate_to_first: dict[int, int] = {}
    records: list[dict[str, object]] = []
    for row, gene in enumerate(genes):
        if gene not in first_position:
            first_position[gene] = row
            keep.append(row)
        else:
            first = first_position[gene]
            duplicate_to_first[row] = first
            records.append(
                {
                    "gene_symbol": gene,
                    "first_row_zero_based": first,
                    "duplicate_row_zero_based": row,
                    "action": "sum_raw_UMI_into_first_row",
                }
            )
    return np.asarray(keep, dtype=np.int64), duplicate_to_first, pd.DataFrame(records)


def merge_duplicate_gene_rows(
    raw_counts: np.ndarray,
    background_mean: np.ndarray,
    keep_rows: np.ndarray,
    duplicate_to_first: dict[int, int],
) -> tuple[np.ndarray, np.ndarray]:
    original_sum = int(raw_counts.sum(dtype=np.int64))
    for duplicate, first in duplicate_to_first.items():
        merged = raw_counts[first, :].astype(np.int64) + raw_counts[duplicate, :].astype(np.int64)
        if merged.max(initial=0) > np.iinfo(np.int32).max:
            raise BuildError("Merged duplicate-gene UMI count exceeds int32 range")
        raw_counts[first, :] = merged.astype(np.int32)
        background_mean[first] += background_mean[duplicate]
    raw_unique = raw_counts[keep_rows, :]
    background_unique = background_mean[keep_rows]
    if int(raw_unique.sum(dtype=np.int64)) != original_sum:
        raise BuildError("Raw UMI total changed while merging duplicate gene symbols")
    return raw_unique, background_unique


def clinical_metadata_dict(metadata: pd.DataFrame, sample: str) -> dict[str, str]:
    if sample not in metadata.index:
        raise BuildError(f"Sample {sample} is absent from metadata.csv")
    row = metadata.loc[sample]
    if isinstance(row, pd.DataFrame):
        raise BuildError(f"Sample {sample} occurs more than once in metadata.csv")
    return {str(key): "" if pd.isna(value) else str(value) for key, value in row.items()}


def construct_obs(
    sample: str,
    utt_barcodes: list[str],
    retain: np.ndarray,
    positions: pd.DataFrame,
    pixel: pd.DataFrame,
    corrected_total: np.ndarray,
    mito_fraction: np.ndarray,
    z_log_total: np.ndarray,
    z_mito: np.ndarray,
) -> tuple[pd.DataFrame, np.ndarray, pd.DataFrame]:
    missing_positions = sorted(set(utt_barcodes).difference(positions.index))
    if missing_positions:
        raise BuildError(f"{sample}: UTT barcodes missing from positions: {missing_positions[:5]}")
    ordered_positions = positions.loc[utt_barcodes].copy()
    if not (ordered_positions["in_tissue"].astype(int) == 1).all():
        raise BuildError(f"{sample}: one or more UTT count columns are not marked in-tissue")

    all_obs = pd.DataFrame(index=pd.Index(utt_barcodes, name="barcode"))
    for column in POSITION_COLUMNS[1:]:
        all_obs[column] = ordered_positions[column].to_numpy()
    all_obs["qc_corrected_total_umi"] = corrected_total
    all_obs["qc_log10_corrected_total_umi"] = np.log10(corrected_total)
    all_obs["qc_z_log10_corrected_total_umi"] = z_log_total
    all_obs["qc_corrected_mito_fraction"] = mito_fraction
    all_obs["qc_z_corrected_mito_fraction"] = z_mito

    pixel_keys = pd.Index([f"{sample}_{barcode}" for barcode in utt_barcodes])
    pixel_subset = pixel.reindex(pixel_keys)
    for column in PIXEL_COLUMNS:
        all_obs[column] = pixel_subset[column].to_numpy(dtype=np.float32)

    obs = all_obs.loc[np.asarray(utt_barcodes)[retain]].copy()
    obs["patient_original"] = sample
    obs["patient_publication"] = PAPER_IDS[sample]
    obs["qc_pass"] = True
    obs["cauchy_all"] = pd.Categorical(["All_spots"] * len(obs), categories=["All_spots"])

    spatial = obs[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy(dtype=np.float64)
    if spatial.shape != (len(obs), 2) or not np.isfinite(spatial).all():
        raise BuildError(f"{sample}: invalid full-resolution spatial coordinates")
    return obs, spatial, all_obs


def build_anndata(
    sample: str,
    raw_unique: np.ndarray,
    unique_genes: np.ndarray,
    background_unique: np.ndarray,
    obs: pd.DataFrame,
    spatial: np.ndarray,
    raw_dir: Path,
    clinical: dict[str, str],
) -> ad.AnnData:
    matrix = sparse.csr_matrix(raw_unique.T, dtype=np.int32)
    matrix.eliminate_zeros()
    if np.any(matrix.data < 0) or not np.issubdtype(matrix.dtype, np.integer):
        raise BuildError(f"{sample}: final sparse count matrix is not non-negative integer data")

    var = pd.DataFrame(index=pd.Index(unique_genes, name="gene_symbol"))
    var["gene_symbol"] = unique_genes
    var["off_tissue_mean_umi"] = background_unique.astype(np.float64)
    if var.index.has_duplicates:
        raise BuildError(f"{sample}: duplicate var_names remain after merge")

    adata = ad.AnnData(X=matrix.copy(), obs=obs, var=var)
    adata.layers["count"] = matrix.copy()
    adata.obsm["spatial"] = spatial

    image_path = raw_dir / f"{sample}_tissue_hires_image.png"
    with Image.open(image_path) as image_handle:
        hires_image = np.asarray(image_handle.convert("RGB"))
    with (raw_dir / f"{sample}_scalefactors_json.json").open("r", encoding="utf-8") as handle:
        scalefactors = json.load(handle)
    scalefactors = {str(key): float(value) for key, value in scalefactors.items()}

    adata.uns["spatial"] = {
        sample: {
            "images": {"hires": hires_image},
            "scalefactors": scalefactors,
            "metadata": {"source": "Harnik_2024_RAW", "sample": sample},
        }
    }
    adata.uns["clinical_metadata"] = clinical
    adata.uns["build_provenance"] = {
        "script": "gsMap_GI_spatial_analysis.py",
        "script_version": SCRIPT_VERSION,
        "counts_in_X_and_count_layer": "original_raw_integer_UTT_UMI_after_spot_QC",
        "background_correction_use": "QC_only_not_model_input",
        "model_annotation": "none",
        "cauchy_column": "cauchy_all",
        "cauchy_value": "All_spots",
        "duplicate_gene_rule": "sum_raw_UMI_by_exact_gene_symbol",
        "spatial_coordinate_order": "x=pxl_col_in_fullres,y=pxl_row_in_fullres",
    }
    return adata


# Background correction selects spots only; native gsMap receives raw integer UMI counts.
def prepare_visium_inputs() -> None:
    H5AD_DIR.mkdir(parents=True, exist_ok=True)
    pixel = read_pixel_metadata(RAW_DIR / "Pixel_classifier_output.txt")
    metadata = read_metadata(RAW_DIR / "metadata.csv")
    n_total = 0

    for sample in SAMPLES:
        all_path = RAW_DIR / f"{sample}_counts_ALL.csv"
        utt_path = RAW_DIR / f"{sample}_counts_UTT.csv"
        _, all_barcodes = read_count_header(all_path)
        _, utt_barcodes = read_count_header(utt_path)
        off_tissue = [b for b in all_barcodes if b not in set(utt_barcodes)]
        off_counts = read_counts(all_path, off_tissue)
        genes = off_counts.index.to_numpy(dtype=str)
        background = off_counts.to_numpy(copy=False).mean(axis=1, dtype=np.float64)
        del off_counts

        utt = read_counts(utt_path)
        assert np.array_equal(utt.index.to_numpy(dtype=str), genes)
        umi = utt.to_numpy(copy=False)
        corrected_total, mito_fraction, z_total, z_mito = calculate_background_corrected_qc(
            umi, background, genes
        )
        excluded = read_manual_exclusions(
            RAW_DIR / "manual_annotation" / f"spots2exclude_{sample.lower()}.csv"
        )
        keep = (z_total > LOW_UMI_Z_THRESHOLD) & (z_mito < MITO_Z_THRESHOLD)
        keep &= np.fromiter((b not in excluded for b in utt_barcodes), dtype=bool)
        assert int(keep.sum()) == EXPECTED_RETAINED_BY_SAMPLE[sample]
        n_total += int(keep.sum())

        positions = read_positions(RAW_DIR / f"{sample}_tissue_positions_list.csv")
        obs, spatial, _ = construct_obs(
            sample, utt_barcodes, keep, positions, pixel,
            corrected_total, mito_fraction, z_total, z_mito,
        )
        kept_rows, duplicate_to_first, _ = duplicate_plan(genes)
        raw_unique, background_unique = merge_duplicate_gene_rows(
            umi[:, keep].copy(), background.copy(), kept_rows, duplicate_to_first
        )
        adata = build_anndata(
            sample, raw_unique, genes[kept_rows], background_unique,
            obs, spatial, RAW_DIR, clinical_metadata_dict(metadata, sample)
        )
        adata.write_h5ad(H5AD_DIR / f"{sample}.h5ad", compression="gzip")
        del adata, utt, raw_unique, umi
        gc.collect()

    assert n_total == EXPECTED_TOTAL_RETAINED, n_total


# Native gsMap 1.73.8 pipeline; parameters match the original project commands.
def gsmap(*args: str) -> None:
    subprocess.run([sys.executable, "-m", "gsMap", *map(str, args)], check=True)


def native_gsmap_analysis() -> None:
    WORKDIR.mkdir(parents=True, exist_ok=True)
    slice_mean = WORKDIR / "Harnik_8samples_slice_mean.parquet"
    gsmap(
        "create_slice_mean", "--sample_name_list", *SAMPLES,
        "--h5ad_list", *(H5AD_DIR / f"{s}.h5ad" for s in SAMPLES),
        "--slice_mean_output_file", slice_mean, "--data_layer", "count"
    )

    # Annotation-free training: do not pass cauchy_all as a modeling annotation.
    for sample in SAMPLES:
        gsmap(
            "run_find_latent_representations", "--workdir", WORKDIR,
            "--sample_name", sample, "--input_hdf5_path", H5AD_DIR / f"{sample}.h5ad",
            "--data_layer", "count", "--epochs", "300",
            "--feat_hidden1", "256", "--feat_hidden2", "128",
            "--gat_hidden1", "64", "--gat_hidden2", "30",
            "--p_drop", "0.1", "--gat_lr", "0.001",
            "--n_neighbors", "11", "--n_comps", "300",
            "--convergence_threshold", "0.0001",
        )

    for sample in SAMPLES:
        gsmap(
            "run_latent_to_gene", "--workdir", WORKDIR, "--sample_name", sample,
            "--latent_representation", "latent_GVAE", "--num_neighbour", "21",
            "--num_neighbour_spatial", "101", "--gM_slices", slice_mean,
        )

    bfile = RESOURCE_DIR / "LD_Reference_Panel" / "1000G_EUR_Phase3_plink" / "1000G.EUR.QC"
    hm3 = RESOURCE_DIR / "LDSC_resource" / "hapmap3_snps" / "hm"
    gtf = RESOURCE_DIR / "genome_annotation" / "gtf" / "gencode.v46lift37.basic.annotation.gtf"
    for sample in SAMPLES:
        for chromosome in range(1, 23):
            gsmap(
                "run_generate_ldscore", "--workdir", WORKDIR, "--sample_name", sample,
                "--chrom", chromosome, "--bfile_root", bfile, "--keep_snp_root", hm3,
                "--gtf_annotation_file", gtf, "--gene_window_size", "50000",
                "--spots_per_chunk", "250", "--ld_wind", "1", "--ld_unit", "CM",
            )

    weights = RESOURCE_DIR / "LDSC_resource" / "weights_hm3_no_hla" / "weights."
    for sample in SAMPLES:
        for trait in TRAITS:
            gsmap(
                "run_spatial_ldsc", "--workdir", WORKDIR, "--sample_name", sample,
                "--sumstats_file", GWAS_DIR / f"{trait}.sumstats.gz",
                "--w_file", weights, "--trait_name", trait,
                "--n_blocks", "200", "--num_processes", "1",
            )

    # First-round overall (three tests) and section-level (24 tests) ACAT.
    output = WORKDIR / "cauchy_combination"
    output.mkdir(exist_ok=True)
    for trait in TRAITS:
        gsmap(
            "run_cauchy_combination", "--workdir", WORKDIR,
            "--sample_name_list", *SAMPLES, "--trait_name", trait,
            "--annotation", "cauchy_all", "--output_file",
            output / f"Harnik_8samples_{trait}_cauchy_combination.csv.gz",
        )
        for sample in SAMPLES:
            gsmap(
                "run_cauchy_combination", "--workdir", WORKDIR,
                "--sample_name", sample, "--trait_name", trait,
                "--annotation", "cauchy_all",
            )
    # Bonferroni thresholds: 0.05/3 (overall); 0.05/24 (section-level).


@dataclass
class AnnotationResult:
    sample: str
    table: pd.DataFrame
    percentile_edges_by_pass: list[dict[str, Any]]
    nearest_spot_distance_hires_pixels: float
    input_anchor_count: int
    retained_anchor_count: int
    retained_mucosal_count: int
    retained_muscularis_count: int
    retained_submucosa_count: int


def require_columns(frame: pd.DataFrame, columns: Iterable[str], label: str) -> None:
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise AnatomyError(f"{label} lacks required columns: {missing}")


def matlab_percentile(values: np.ndarray, percentiles: np.ndarray) -> np.ndarray:
    """MATLAB R2023a prctile default: midpoint/exact = Hyndman-Fan type 5."""
    values = np.asarray(values, dtype=np.float64)
    if values.ndim != 1 or values.size == 0 or not np.isfinite(values).all():
        raise AnatomyError("Percentile input must be a non-empty finite vector.")
    return np.percentile(values, percentiles, method="hazen")


def pairwise_and_spacing(coordinates: np.ndarray) -> tuple[np.ndarray, float]:
    pairwise = cdist(coordinates, coordinates, metric="euclidean")
    np.fill_diagonal(pairwise, np.nan)
    nearest = np.nanmin(pairwise, axis=1)
    positive = nearest[nearest > 0]
    if positive.size == 0:
        raise AnatomyError("Cannot determine a positive nearest-spot distance.")
    minimum = float(np.min(positive))
    return pairwise * (100.0 / minimum), minimum


def apply_author_percentile_passes(
    spot_indices: np.ndarray,
    distance_to_base_um: np.ndarray,
    villus_length_um: np.ndarray,
) -> tuple[np.ndarray, list[dict[str, Any]]]:
    if not (
        len(spot_indices) == len(distance_to_base_um) == len(villus_length_um)
    ):
        raise AnatomyError("Mucosal annotation vectors have unequal lengths.")
    bin_ids = np.zeros(len(spot_indices), dtype=np.int8)
    percentiles = np.linspace(0.0, 100.0, N_BINS + 1)
    pass_records: list[dict[str, Any]] = []

    # Preserve both overlapping source passes.  The second [0,Inf] pass writes
    # over all final bin IDs, exactly as the MATLAB function does.
    for lower, upper in VILLUS_LENGTH_PASSES:
        in_pass = np.flatnonzero(
            (villus_length_um >= lower) & (villus_length_um <= upper)
        )
        if in_pass.size == 0:
            raise AnatomyError(f"Author villus-length pass [{lower}, {upper}] is empty.")
        edges = matlab_percentile(distance_to_base_um[in_pass], percentiles)
        for j in range(N_BINS):
            selected = in_pass[
                (distance_to_base_um[in_pass] >= edges[j])
                & (distance_to_base_um[in_pass] <= edges[j + 1])
            ]
            # Inclusive boundaries plus ascending assignment reproduce MATLAB:
            # shared-boundary spots are overwritten into the later bin.
            bin_ids[selected] = j + 1
        if np.any(bin_ids[in_pass] == 0):
            raise AnatomyError("One or more mucosal spots were not percentile-binned.")
        pass_records.append(
            {
                "lower_villus_length_um": lower,
                "upper_villus_length_um": "Inf" if math.isinf(upper) else upper,
                "spots": int(in_pass.size),
                "percentile_edges_um": [float(value) for value in edges],
            }
        )
    return bin_ids, pass_records


def empty_annotation_table(sample: str, barcodes: pd.Index) -> pd.DataFrame:
    table = pd.DataFrame(index=pd.Index(barcodes.astype(str), name="barcode"))
    table.insert(0, "patient_original", sample)
    table["villus_bin"] = np.nan
    table["anatomical_region"] = pd.NA
    table["anatomy_labeled"] = False
    table["anatomy_component"] = "Unassigned"
    table["distance_to_base_um"] = np.nan
    table["villus_length_um"] = np.nan
    table["closest_base_barcode"] = pd.NA
    table["author_percentile_bin_1_to_7"] = np.nan
    return table


def assign_labels(
    table: pd.DataFrame,
    mucosal_indices: np.ndarray,
    mucosal_bins: np.ndarray,
    distance_to_base_um: np.ndarray,
    villus_length_um: np.ndarray,
    closest_base_indices: np.ndarray,
    muscularis_indices: np.ndarray,
    submucosa_indices: np.ndarray,
) -> None:
    barcodes = table.index.to_numpy(dtype=str)
    if len(np.intersect1d(muscularis_indices, submucosa_indices)):
        raise AnatomyError("Muscularis and submucosa index sets overlap.")
    if len(np.intersect1d(mucosal_indices, muscularis_indices)) or len(
        np.intersect1d(mucosal_indices, submucosa_indices)
    ):
        raise AnatomyError("Mucosal and submucosal label index sets overlap.")

    mucosal_codes = mucosal_bins.astype(np.int16) - 1
    table.iloc[mucosal_indices, table.columns.get_loc("villus_bin")] = mucosal_codes
    table.iloc[mucosal_indices, table.columns.get_loc("anatomy_component")] = "Mucosa"
    table.iloc[mucosal_indices, table.columns.get_loc("distance_to_base_um")] = distance_to_base_um
    table.iloc[mucosal_indices, table.columns.get_loc("villus_length_um")] = villus_length_um
    table.iloc[mucosal_indices, table.columns.get_loc("closest_base_barcode")] = barcodes[
        closest_base_indices
    ]
    table.iloc[
        mucosal_indices, table.columns.get_loc("author_percentile_bin_1_to_7")
    ] = mucosal_bins

    table.iloc[muscularis_indices, table.columns.get_loc("villus_bin")] = -1
    table.iloc[muscularis_indices, table.columns.get_loc("anatomy_component")] = "Muscularis_mucosa"
    table.iloc[submucosa_indices, table.columns.get_loc("villus_bin")] = -2
    table.iloc[submucosa_indices, table.columns.get_loc("anatomy_component")] = "Submucosa"

    labeled = table["villus_bin"].notna()
    codes = table.loc[labeled, "villus_bin"].astype(int)
    if not set(codes).issubset(REGION_BY_CODE):
        raise AnatomyError("Unexpected anatomical code was generated.")
    table.loc[labeled, "anatomical_region"] = codes.map(REGION_BY_CODE).to_numpy()
    table.loc[labeled, "anatomy_labeled"] = True


def read_first_column_barcodes(path: Path, expected_header: str | None = None) -> list[str]:
    frame = pd.read_csv(path, dtype=str)
    if frame.shape[1] < 1:
        raise AnatomyError(f"Annotation file has no columns: {path}")
    if expected_header is not None and expected_header not in frame.columns:
        raise AnatomyError(f"{path} lacks expected column {expected_header!r}.")
    values = frame.iloc[:, 0].dropna().astype(str).str.strip()
    if (values == "").any() or values.duplicated().any():
        raise AnatomyError(f"Blank or duplicate barcodes in {path}.")
    return values.tolist()


def annotate_standard_sample(
    paths: dict[str, Path], sample: str, obs: pd.DataFrame, coordinates: np.ndarray
) -> AnnotationResult:
    table = empty_annotation_table(sample, obs.index)
    anchor_path = paths["raw"] / "manual_annotation" / f"crypt_base_{sample}.csv"
    input_anchors = read_first_column_barcodes(anchor_path, "Barcode")
    anchor_set = set(input_anchors)
    is_anchor = np.fromiter((barcode in anchor_set for barcode in obs.index), dtype=bool)
    submucosal_indices = np.flatnonzero(is_anchor)
    epithelial_indices = np.flatnonzero(~is_anchor)
    if submucosal_indices.size == 0 or epithelial_indices.size == 0:
        raise AnatomyError(f"{sample} has an empty anchor or epithelial set.")

    distance_all_to_anchor = cdist(coordinates, coordinates[submucosal_indices])
    distance_closest = np.min(distance_all_to_anchor, axis=1)
    closest_anchor = submucosal_indices[np.argmin(distance_all_to_anchor, axis=1)]

    distance_anchor_to_epithelium = cdist(
        coordinates[submucosal_indices], coordinates[epithelial_indices]
    )
    distance_to_epithelium = np.min(distance_anchor_to_epithelium, axis=1)
    median_distance = float(np.median(distance_to_epithelium))
    muscularis_indices = submucosal_indices[distance_to_epithelium < median_distance]
    submucosa_indices = submucosal_indices[distance_to_epithelium >= median_distance]

    pairwise_um, minimum_spacing = pairwise_and_spacing(coordinates)
    distance_closest_um = distance_closest * (100.0 / minimum_spacing)

    included_groups: list[tuple[int, np.ndarray]] = []
    for anchor_index in np.unique(closest_anchor[epithelial_indices]):
        group = np.intersect1d(
            np.flatnonzero(closest_anchor == anchor_index), epithelial_indices
        )
        if group.size == 1:
            nearest_within = np.array([np.nan])
        else:
            group_distances = pairwise_um[np.ix_(group, group)]
            nearest_within = np.nanmin(group_distances, axis=1)
        # Exact source order: connectivity check first, then length >= 3.
        if np.all(nearest_within <= WITHIN_VILLUS_NEIGHBOUR_UM) and (
            group.size >= VILLI_LEN_THRESH
        ):
            included_groups.append((int(anchor_index), group))

    if not included_groups:
        raise AnatomyError(f"{sample} has no villus groups after author filters.")

    mucosal_indices = np.concatenate([group for _, group in included_groups])
    closest_for_mucosa = np.concatenate(
        [np.full(group.size, anchor, dtype=np.int64) for anchor, group in included_groups]
    )
    mucosal_distance = distance_closest_um[mucosal_indices]
    mucosal_villus_length = np.concatenate(
        [
            np.full(group.size, np.max(distance_closest_um[group]), dtype=np.float64)
            for _, group in included_groups
        ]
    )
    mucosal_bins, pass_records = apply_author_percentile_passes(
        mucosal_indices, mucosal_distance, mucosal_villus_length
    )
    assign_labels(
        table,
        mucosal_indices,
        mucosal_bins,
        mucosal_distance,
        mucosal_villus_length,
        closest_for_mucosa,
        muscularis_indices,
        submucosa_indices,
    )
    del pairwise_um, distance_all_to_anchor, distance_anchor_to_epithelium
    gc.collect()

    return AnnotationResult(
        sample=sample,
        table=table,
        percentile_edges_by_pass=pass_records,
        nearest_spot_distance_hires_pixels=minimum_spacing,
        input_anchor_count=len(input_anchors),
        retained_anchor_count=int(submucosal_indices.size),
        retained_mucosal_count=int(mucosal_indices.size),
        retained_muscularis_count=int(muscularis_indices.size),
        retained_submucosa_count=int(submucosa_indices.size),
    )


def annotate_p10(
    paths: dict[str, Path], obs: pd.DataFrame, coordinates: np.ndarray
) -> AnnotationResult:
    sample = "P10"
    table = empty_annotation_table(sample, obs.index)
    dist_path = paths["raw"] / "manual_annotation" / "P10_dist2submuc.csv"
    graph_path = paths["raw"] / "manual_annotation" / "P10_Graph-Based.csv"

    distance_table = pd.read_csv(dist_path, dtype=str)
    require_columns(distance_table, ["Barcode", "dist2submuc"], str(dist_path))
    if distance_table["Barcode"].duplicated().any():
        raise AnatomyError("P10_dist2submuc has duplicate barcodes.")
    distance_map = dict(
        zip(
            distance_table["Barcode"].astype(str).str.strip(),
            distance_table["dist2submuc"].astype(str).str.strip(),
        )
    )
    # Exact source default for a retained spot absent from the table.
    distance_annotation = np.array(
        [distance_map.get(str(barcode), "non_epi_out") for barcode in obs.index],
        dtype=object,
    )

    graph_table = pd.read_csv(graph_path, dtype=str)
    require_columns(graph_table, ["Barcode"], str(graph_path))
    graph_barcodes = graph_table["Barcode"].dropna().astype(str).str.strip()
    if graph_barcodes.duplicated().any():
        raise AnatomyError("P10_Graph-Based has duplicate barcodes.")
    serosa_set = set(graph_barcodes)
    is_serosa = np.fromiter((barcode in serosa_set for barcode in obs.index), dtype=bool)

    # Source-compatibility note: the public MATLAB function accesses ``Var1``
    # although the public CSV's first-column header is ``Barcode``.  Reading
    # the first/barcode column preserves the code's evident operation: every
    # listed graph-based spot is removed from candidate submucosal anchors.
    pairwise_um, minimum_spacing = pairwise_and_spacing(coordinates)
    del pairwise_um  # P10 uses the global spacing only, not connectivity filtering.

    epithelial_all: list[np.ndarray] = []
    closest_all: list[np.ndarray] = []
    distance_all: list[np.ndarray] = []
    muscularis_all: list[np.ndarray] = []
    submucosa_all: list[np.ndarray] = []
    retained_anchor_total = 0

    for branch in ("inn", "out"):
        branch_labels = sorted(
            {str(value) for value in distance_annotation if branch in str(value)}
        )
        expected = [f"epi_{branch}", f"non_epi_{branch}"]
        if branch_labels != expected:
            raise AnatomyError(
                f"P10 {branch} labels are {branch_labels}; expected {expected}."
            )
        epithelial_indices = np.flatnonzero(distance_annotation == branch_labels[0])
        submucosal_indices = np.flatnonzero(distance_annotation == branch_labels[1])
        submucosal_indices = np.setdiff1d(
            submucosal_indices, np.flatnonzero(is_serosa)
        )
        if epithelial_indices.size == 0 or submucosal_indices.size == 0:
            raise AnatomyError(f"P10 {branch} has an empty epithelial or anchor set.")
        retained_anchor_total += int(submucosal_indices.size)

        distances = cdist(coordinates[epithelial_indices], coordinates[submucosal_indices])
        distance_closest_um = np.min(distances, axis=1) * (100.0 / minimum_spacing)
        closest_anchor = submucosal_indices[np.argmin(distances, axis=1)]

        anchor_to_epithelium = cdist(
            coordinates[submucosal_indices], coordinates[epithelial_indices]
        )
        distance_to_epithelium = np.min(anchor_to_epithelium, axis=1)
        muscularis_all.append(
            submucosal_indices[
                distance_to_epithelium < P10_MUSCULARIS_THRESHOLD_HIRES_PIXELS
            ]
        )
        submucosa_all.append(
            submucosal_indices[
                distance_to_epithelium >= P10_MUSCULARIS_THRESHOLD_HIRES_PIXELS
            ]
        )
        epithelial_all.append(epithelial_indices)
        closest_all.append(closest_anchor)
        distance_all.append(distance_closest_um)

    epithelial_merged = np.concatenate(epithelial_all)
    closest_merged = np.concatenate(closest_all)
    distance_merged = np.concatenate(distance_all)
    muscularis_indices = np.concatenate(muscularis_all)
    submucosa_indices = np.concatenate(submucosa_all)

    included_positions: list[np.ndarray] = []
    for anchor_index in np.unique(closest_merged):
        positions = np.flatnonzero(closest_merged == anchor_index)
        if positions.size >= VILLI_LEN_THRESH:
            included_positions.append(positions)
    if not included_positions:
        raise AnatomyError("P10 has no villus groups after the author length filter.")

    mucosal_positions = np.concatenate(included_positions)
    mucosal_indices = epithelial_merged[mucosal_positions]
    mucosal_distance = distance_merged[mucosal_positions]
    closest_for_mucosa = closest_merged[mucosal_positions]
    mucosal_villus_length = np.concatenate(
        [
            np.full(
                positions.size,
                np.max(distance_merged[positions]),
                dtype=np.float64,
            )
            for positions in included_positions
        ]
    )
    mucosal_bins, pass_records = apply_author_percentile_passes(
        mucosal_indices, mucosal_distance, mucosal_villus_length
    )
    assign_labels(
        table,
        mucosal_indices,
        mucosal_bins,
        mucosal_distance,
        mucosal_villus_length,
        closest_for_mucosa,
        muscularis_indices,
        submucosa_indices,
    )
    gc.collect()

    return AnnotationResult(
        sample=sample,
        table=table,
        percentile_edges_by_pass=pass_records,
        nearest_spot_distance_hires_pixels=minimum_spacing,
        input_anchor_count=int((distance_table["dist2submuc"].str.contains("non_epi")).sum()),
        retained_anchor_count=retained_anchor_total,
        retained_mucosal_count=int(mucosal_indices.size),
        retained_muscularis_count=int(muscularis_indices.size),
        retained_submucosa_count=int(submucosa_indices.size),
    )


# Regional analysis reuses Step 4 outputs without recalculating spatial LD scores.
def regional_anatomy_analysis() -> None:
    ANNOTATION_DIR.mkdir(parents=True, exist_ok=True)
    paths = {"raw": RAW_DIR}
    n_spots = 0
    n_labeled = 0

    for sample in SAMPLES:
        source = H5AD_DIR / f"{sample}.h5ad"
        a = ad.read_h5ad(source, backed="r")
        try:
            obs = pd.DataFrame(index=pd.Index(a.obs_names.astype(str), name="barcode"))
            xy = np.asarray(a.obsm["spatial"], dtype=np.float64)
        finally:
            a.file.close()
        with (RAW_DIR / f"{sample}_scalefactors_json.json").open(encoding="utf-8") as f:
            hires_scale = float(json.load(f)["tissue_hires_scalef"])
        result = (annotate_p10(paths, obs, xy * hires_scale) if sample == "P10"
                  else annotate_standard_sample(paths, sample, obs, xy * hires_scale))
        result.table.reset_index().to_csv(
            ANNOTATION_DIR / f"{sample}_harnik_anatomy.csv.gz", index=False
        )
        n_spots += len(result.table)
        sample_labeled = int(result.table["anatomy_labeled"].sum())
        assert sample_labeled == EXPECTED_LABELED_BY_SAMPLE[sample]
        n_labeled += sample_labeled

        latent_path = (WORKDIR / sample / "find_latent_representations" /
                       f"{sample}_add_latent.h5ad")
        a = ad.read_h5ad(latent_path)
        labels = result.table.reindex(a.obs_names.astype(str))
        subset = a[labels["anatomy_labeled"].to_numpy(dtype=bool)].copy()
        selected = labels.loc[subset.obs_names.astype(str)]
        subset.obs["villus_bin"] = selected["villus_bin"].astype(np.int8).to_numpy()
        subset.obs["anatomical_region"] = pd.Categorical(
            selected["anatomical_region"].astype(str), categories=REGION_ORDER, ordered=True
        )
        subset.obs["anatomy_component"] = pd.Categorical(
            selected["anatomy_component"].astype(str)
        )
        target = (ANATOMY_WORKDIR / sample / "find_latent_representations" /
                  f"{sample}_add_latent.h5ad")
        target.parent.mkdir(parents=True, exist_ok=True)
        subset.write_h5ad(target, compression="gzip")
        del a, subset

        # The full spot-level spatial-LDSC result remains unchanged.
        for trait in TRAITS:
            original = WORKDIR / sample / "spatial_ldsc" / f"{sample}_{trait}.csv.gz"
            destination = (ANATOMY_WORKDIR / sample / "spatial_ldsc" / original.name)
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(original, destination)

    assert n_spots == EXPECTED_TOTAL_RETAINED, n_spots
    assert n_labeled == EXPECTED_TOTAL_LABELED, n_labeled
    output = ANATOMY_WORKDIR / "cauchy_combination"
    output.mkdir(parents=True, exist_ok=True)
    for trait in TRAITS:
        gsmap(
            "run_cauchy_combination", "--workdir", ANATOMY_WORKDIR,
            "--sample_name_list", *SAMPLES, "--trait_name", trait,
            "--annotation", "anatomical_region", "--output_file",
            output / f"Harnik_8samples_anatomy_{trait}_cauchy_combination.csv.gz",
        )
    # Nine regions x three factors: Bonferroni threshold 0.05/27.
    # Unlabeled spots remain in the first-round analysis but not in this branch.


def main() -> None:
    parser = argparse.ArgumentParser(description="Genomic SEM factor gsMap workflow")
    parser.add_argument("--stage", required=True, choices=("prepare", "gsmap", "anatomy"))
    stage = parser.parse_args().stage
    if stage == "prepare":
        prepare_visium_inputs()
    elif stage == "gsmap":
        native_gsmap_analysis()
    else:
        regional_anatomy_analysis()


if __name__ == "__main__":
    main()
