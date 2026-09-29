"""BetaReg-specific Visium HD input, validation, and export helpers.

This module contains the compatibility code needed for these Space Ranger
outputs.  It is deliberately separate from the small SPARKLE runner so none of
these functions can be mistaken for functions supplied by ``stambient``.
"""

import csv
import gc
import gzip
import json
import logging
import shutil
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np
from scipy import sparse
from scipy.io import mmread, mmwrite


@dataclass
class BetaRegSparkleInput:
    """Validated inputs and identifiers for one BetaReg SPARKLE run."""

    sparkle_data: dict
    features: list
    barcodes: list
    feature_row_index: np.ndarray
    original_counts: sparse.csr_matrix
    pixel_size_um: float


def _read_gzip_tsv(path):
    with gzip.open(path, "rt") as handle:
        return list(csv.reader(handle, delimiter="\t"))


def _write_gzip_tsv(path, rows):
    with gzip.open(path, "wt") as handle:
        csv.writer(handle, delimiter="\t", lineterminator="\n").writerows(rows)


def make_json_serializable(value):
    """Convert NumPy values in SPARKLE diagnostics into JSON values."""
    if isinstance(value, dict):
        return {str(key): make_json_serializable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [make_json_serializable(item) for item in value]
    if isinstance(value, np.generic):
        return make_json_serializable(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def _load_betareg_feature_slice(feature_slice, pixel_size_um):
    """Load the capture units required for the BetaReg Space Ranger outputs.

    ``stambient.io_utils.load_visiumhd()`` version 0.1.3 allocates its lookup
    from the filtered-tissue bounds.  In these files, feature coordinates and
    segmented cells extend beyond those bounds, causing an IndexError and
    omitting boundary-cell counts.  This compatibility loader instead uses the
    full slide dimensions, includes all segmented units, and restricts
    background units to the filtered-tissue mask.
    """
    with h5py.File(feature_slice, "r") as handle:
        metadata = json.loads(handle.attrs["metadata_json"])
        tissue_rows = handle["masks/filtered/row"][:]
        tissue_cols = handle["masks/filtered/col"][:]
        segmentation = handle["segmentations/cell_segmentation_mask"]

        selected = np.zeros((metadata["nrows"], metadata["ncols"]), dtype=bool)
        selected[tissue_rows, tissue_cols] = True
        selected[segmentation["row"][:], segmentation["col"][:]] = True
        rows, cols = np.nonzero(selected)
        del selected

        coordinate_lookup = np.full(
            (metadata["nrows"], metadata["ncols"]), -1, dtype=np.int32
        )
        coordinate_lookup[rows, cols] = np.arange(len(rows), dtype=np.int32)
        coordinates = (np.column_stack((cols, rows)).astype(float) + 0.5) * pixel_size_um

        feature_keys = sorted(handle["feature_slices"], key=int)
        nonzero_count = 0
        for key in feature_keys:
            feature = handle["feature_slices"][key]
            indices = coordinate_lookup[feature["row"][:], feature["col"][:]]
            nonzero_count += np.count_nonzero(indices >= 0)

        matrix_rows = np.empty(nonzero_count, dtype=np.int32)
        matrix_cols = np.empty(nonzero_count, dtype=np.int32)
        matrix_values = np.empty(nonzero_count, dtype=np.float64)
        offset = 0
        for feature_index, key in enumerate(feature_keys):
            feature = handle["feature_slices"][key]
            indices = coordinate_lookup[feature["row"][:], feature["col"][:]]
            keep = indices >= 0
            end = offset + np.count_nonzero(keep)
            matrix_rows[offset:end] = feature_index
            matrix_cols[offset:end] = indices[keep]
            matrix_values[offset:end] = feature["data"][:][keep]
            offset = end

        expression = sparse.csr_matrix(
            (matrix_values, (matrix_rows, matrix_cols)),
            shape=(len(feature_keys), len(rows)),
        )

        segmentation_indices = coordinate_lookup[
            segmentation["row"][:], segmentation["col"][:]
        ]
        keep = segmentation_indices >= 0
        cell_ids, cell_index = np.unique(
            segmentation["data"][:][keep], return_inverse=True
        )
        if len(np.unique(segmentation_indices[keep])) != np.count_nonzero(keep):
            raise ValueError("Multiple segmentation assignments for a capture unit.")
        labels = np.full(len(rows), -1, dtype=np.int32)
        labels[segmentation_indices[keep]] = cell_index

    logging.info(
        "Loaded %s genes, %s capture units, %s cells, %s background units.",
        len(feature_keys), len(rows), len(cell_ids), np.count_nonzero(labels < 0),
    )
    return {
        "spot_expr": expression,
        "spot_coords": coordinates,
        "spot_labels": labels,
        "cell_ids": cell_ids,
    }


def prepare_betareg_sparkle_input(space_ranger_outs):
    """Load one sample and verify reconstructed counts against Space Ranger."""
    space_ranger_outs = Path(space_ranger_outs)
    feature_slice = space_ranger_outs / "feature_slice.h5"
    matrix_dir = space_ranger_outs / "segmented_outputs/raw_feature_cell_matrix"
    features = _read_gzip_tsv(matrix_dir / "features.tsv.gz")
    barcodes = [row[0] for row in _read_gzip_tsv(matrix_dir / "barcodes.tsv.gz")]
    feature_ids = [row[0] for row in features]
    if len(set(feature_ids)) != len(feature_ids) or len(set(barcodes)) != len(barcodes):
        raise ValueError("Space Ranger gene IDs and barcodes must be unique.")

    with h5py.File(feature_slice, "r") as handle:
        pixel_size_um = float(json.loads(handle.attrs["metadata_json"])["spot_pitch"])
        if not np.isfinite(pixel_size_um) or pixel_size_um <= 0:
            raise ValueError("Invalid spot pitch.")
        feature_indices = sorted(int(key) for key in handle["feature_slices"])
        h5_feature_ids = handle["features/id"].asstr()[:]
        loaded_feature_ids = [h5_feature_ids[index] for index in feature_indices]
        if np.any(handle["segmentations/cell_segmentation_mask/data"][:] <= 0):
            raise ValueError("Non-positive cell labels; inspect the segmentation mask.")

    sparkle_data = _load_betareg_feature_slice(feature_slice, pixel_size_um)
    feature_lookup = {feature_id: index for index, feature_id in enumerate(feature_ids)}
    feature_row_index = np.array(
        [feature_lookup[feature_id] for feature_id in loaded_feature_ids]
    )
    if len(np.unique(feature_row_index)) != len(feature_row_index):
        raise ValueError("Duplicate loaded gene IDs.")

    cell_barcodes = [
        f"cellid_{int(cell_id):09d}-1" for cell_id in sparkle_data["cell_ids"]
    ]
    barcode_lookup = {barcode: index for index, barcode in enumerate(barcodes)}
    barcode_index = np.array([barcode_lookup[barcode] for barcode in cell_barcodes])

    labels = sparkle_data["spot_labels"]
    intracellular = np.flatnonzero(labels >= 0)
    assignment = sparse.csr_matrix(
        (np.ones(len(intracellular)), (intracellular, labels[intracellular])),
        shape=(len(labels), len(cell_barcodes)),
    )
    reconstructed = (sparkle_data["spot_expr"] @ assignment).tocsr()
    with gzip.open(matrix_dir / "matrix.mtx.gz", "rb") as handle:
        space_ranger_counts = mmread(handle, spmatrix=True).tocsr()
    if space_ranger_counts.shape != (len(features), len(barcodes)):
        raise ValueError("Space Ranger matrix dimensions disagree with identifiers.")

    original_counts = space_ranger_counts[:, barcode_index].tocsr()
    difference = reconstructed - original_counts[feature_row_index, :]
    difference.eliminate_zeros()
    absent_features = np.setdiff1d(np.arange(len(features)), feature_row_index)
    if difference.nnz or original_counts[absent_features, :].count_nonzero():
        raise ValueError(
            f"Mask counts differ from Space Ranger ({difference.nnz} entries in "
            "loaded genes). Resolve this before correction."
        )

    logging.info(
        "Exact count agreement: %s loaded genes, %s cells, %s absent zero genes.",
        len(feature_row_index), len(cell_barcodes), len(absent_features),
    )
    logging.info(
        "One dense working cell matrix: %.2f GiB (not peak RAM).",
        8 * len(feature_row_index) * len(cell_barcodes) / 1024**3,
    )
    return BetaRegSparkleInput(
        sparkle_data=sparkle_data,
        features=features,
        barcodes=cell_barcodes,
        feature_row_index=feature_row_index,
        original_counts=original_counts,
        pixel_size_um=pixel_size_um,
    )


def validate_betareg_sparkle_counts(corrected_counts, run_input):
    """Check shape, finiteness, non-negativity, and monotonic subtraction."""
    expected_shape = (len(run_input.feature_row_index), len(run_input.barcodes))
    if corrected_counts.shape != expected_shape:
        raise ValueError("Unexpected corrected matrix shape.")
    for start in range(0, len(run_input.feature_row_index), 128):
        corrected_block = corrected_counts[start:start + 128]
        original_block = run_input.original_counts[
            run_input.feature_row_index[start:start + 128], :
        ].toarray()
        if (
            not np.isfinite(corrected_block).all()
            or (corrected_block < 0).any()
            or (corrected_block > original_block + 1e-8).any()
        ):
            raise ValueError("Invalid corrected counts or counts increased.")


def write_betareg_sparkle_export(
    output_dir, space_ranger_outs, corrected_counts, run_input, diagnostics,
    r2_threshold,
):
    """Write corrected counts and aligned metadata for the R importer."""
    output_dir = Path(output_dir)
    space_ranger_outs = Path(space_ranger_outs)

    compact = sparse.coo_matrix(corrected_counts)
    corrected_sparse = sparse.csc_matrix(
        (
            compact.data,
            (run_input.feature_row_index[compact.row], compact.col),
        ),
        shape=run_input.original_counts.shape,
    )
    del compact, corrected_counts
    gc.collect()

    for name, matrix in (
        ("corrected_counts", corrected_sparse),
        ("original_counts", run_input.original_counts),
    ):
        logging.info("Writing %s", name)
        with gzip.open(output_dir / f"{name}.mtx.gz", "wb") as handle:
            mmwrite(handle, matrix, field="real", precision=17, symmetry="general")

    _write_gzip_tsv(output_dir / "features.tsv.gz", run_input.features)
    _write_gzip_tsv(
        output_dir / "barcodes.tsv.gz", ((barcode,) for barcode in run_input.barcodes)
    )
    shutil.copyfile(
        space_ranger_outs
        / "segmented_outputs/filtered_feature_cell_matrix/features.tsv.gz",
        output_dir / "reference_features.tsv.gz",
    )

    fit_by_loaded_row = {
        int(index): (alpha, r2)
        for index, alpha, r2 in zip(
            diagnostics["gene_indices"], diagnostics["alphas"], diagnostics["r2_scores"]
        )
    }
    loaded_row_by_export_row = {
        int(export_row): loaded_row
        for loaded_row, export_row in enumerate(run_input.feature_row_index)
    }
    gene_rows = [["gene_id", "alpha", "r2", "correction_eligible"]]
    for export_row, feature in enumerate(run_input.features):
        alpha, r2 = fit_by_loaded_row.get(
            loaded_row_by_export_row.get(export_row), (np.nan, np.nan)
        )
        gene_rows.append(
            [
                feature[0],
                alpha if np.isfinite(alpha) else "NA",
                r2 if np.isfinite(r2) else "NA",
                int(np.isfinite(r2) and r2 >= r2_threshold),
            ]
        )
    _write_gzip_tsv(output_dir / "gene_diagnostics.tsv.gz", gene_rows)
