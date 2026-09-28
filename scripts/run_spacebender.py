#!/usr/bin/env python3
"""CPU pilot for Visium HD 8-um counts.

Usage: python scripts/run_spacebender.py INPUT.h5ad OUTPUT.h5ad
Keep the cloned SpaceBender directory in the project root (beside scripts/),
or run from its parent directory. Activate the spacebender environment first.
Parameters below are starting values, not validated HD-specific settings.
"""

import argparse
import gc
import importlib.metadata
import logging
import os
from pathlib import Path
import subprocess
import sys
import time


SETTINGS = {
    # On a complete square grid, eight neighbours are the surrounding 3x3 ring.
    # Filtering/crop edges can make these neighbours more distant.
    "k_neighbors": 8,
    # Coordinate units may be pixels; do not assume a physical radius.
    "ambient_radius": None,
    "n_epochs": 150,            # Upstream default; inspect loss in the run log.
    "batch_size": 64,           # Smaller than default to limit temporary memory.
    "learning_rate": 1e-3,      # Upstream default optimiser step size.
    "z_dim": 64,                # Default latent expression dimensions.
    "hidden_dims": [512, 512],  # Default neural-network layer widths.
    "n_posterior_samples": 25,  # Default Monte Carlo samples for denoising.
    "device": "cpu",
    "train_in_tissue_only": False,  # Include background bins when provided.
    "random_seed": 1234,
    "lambda_multiplier": 1.0,  # No extra amplification of estimated noise.
    "use_map": True,           # Most probable noise count, rather than a draw.
}


def validate_input(path):
    import anndata as ad
    import numpy as np
    from scipy import sparse

    # Read counts in blocks so validation does not retain a second full matrix.
    a = ad.read_h5ad(path, backed="r")
    try:
        if a.n_obs <= SETTINGS["k_neighbors"] or a.n_vars < 2:
            raise ValueError("Need at least 9 bins and 2 genes.")
        if not a.obs_names.is_unique or not a.var_names.is_unique:
            raise ValueError("Bin and gene identifiers must be unique.")
        if "in_tissue" not in a.obs:
            raise ValueError("Missing obs['in_tissue'].")
        tissue = a.obs["in_tissue"]
        if tissue.isna().any() or not tissue.isin([True, False, 0, 1]).all():
            raise ValueError("in_tissue must contain boolean or numeric 0/1 values.")
        tissue = tissue.to_numpy().astype(bool)
        if not tissue.any():
            raise ValueError("No in-tissue bins were supplied.")
        if "spatial" not in a.obsm:
            raise ValueError("Missing obsm['spatial'].")
        coords = np.asarray(a.obsm["spatial"])
        if coords.shape != (a.n_obs, 2) or not np.isfinite(coords).all():
            raise ValueError("Spatial coordinates must be a finite bins-by-2 array.")
        gene_totals = np.zeros(a.n_vars, dtype=np.float64)
        positive_tissue = False
        zero_bins = 0
        for start in range(0, a.n_obs, 1024):
            block = a.X[start:start + 1024, :]
            values = block.data if sparse.issparse(block) else np.asarray(block)
            if (not np.isfinite(values).all() or (values < 0).any()
                    or not np.equal(values, np.floor(values)).all()):
                raise ValueError("X must contain finite, non-negative integer counts.")
            totals = np.asarray(block.sum(axis=1)).ravel()
            zero_bins += int((totals == 0).sum())
            positive_tissue |= bool((totals[tissue[start:start + 1024]] > 0).any())
            gene_totals += np.asarray(block.sum(axis=0)).ravel()
        expressed = int((gene_totals > 0).sum())
        if not positive_tissue or expressed < 2:
            raise ValueError("Need positive in-tissue counts and at least 2 expressed genes.")
        logging.info("Input: %s bins, %s genes (%s expressed), %s in-tissue",
                     a.n_obs, a.n_vars, expressed, int(tissue.sum()))
        logging.info("One dense float32 ambient matrix: %.2f GiB; preprocessing "
                     "can hold about three such matrices, plus other allocations.",
                     4 * a.n_obs * expressed / 1024**3)
        if tissue.all():
            logging.warning("No out-of-tissue bins: SpaceBender uses its fallback "
                            "background prior. Treat this crop as a pilot.")
        if zero_bins:
            logging.warning("%s bins have zero total counts. Upstream posterior "
                            "code may fail on an entirely zero-count batch.", zero_bins)
    finally:
        a.file.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="H5AD with raw counts and spatial metadata")
    parser.add_argument("output", type=Path, help="New denoised H5AD file")
    args = parser.parse_args()
    src, dst = args.input.resolve(), args.output.resolve()
    if not src.is_file():
        parser.error(f"Input does not exist: {src}")
    if src == dst or dst.exists():
        parser.error("Choose a new output file; existing files will not be overwritten.")
    if not dst.parent.is_dir():
        parser.error(f"Create the output directory first: {dst.parent}")

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.StreamHandler(), logging.FileHandler(str(dst) + ".log", mode="a")],
    )
    # Python normally searches scripts/, not the project root, for this import.
    for parent in (Path(__file__).resolve().parent.parent, Path.cwd()):
        if (parent / "SpaceBender" / "spacebender.py").is_file():
            sys.path.insert(0, str(parent))
            break

    import torch
    from SpaceBender import spacebender

    # Honour a per-task Slurm allocation, otherwise use at most four CPUs.
    available = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)
    threads = max(1, min(available, int(os.environ.get("SLURM_CPUS_PER_TASK", "4"))))
    torch.set_num_threads(threads)
    torch.set_num_interop_threads(1)
    logging.info("Input: %s; output: %s; CPU threads: %s", src, dst, threads)
    logging.info("Settings: %s", SETTINGS)
    logging.info("Python: %s; SpaceBender source: %s", sys.version, spacebender.__file__)
    for package in ("torch", "pyro-ppl", "anndata", "numpy", "scipy"):
        logging.info("%s: %s", package, importlib.metadata.version(package))
    try:
        revision = subprocess.run(
            ["git", "-C", str(Path(spacebender.__file__).parent), "rev-parse", "HEAD"],
            capture_output=True, text=True, timeout=10, check=True,
        ).stdout.strip()
        logging.info("SpaceBender commit: %s", revision)
    except (OSError, subprocess.SubprocessError):
        logging.info("SpaceBender git revision unavailable")

    validate_input(src)
    gc.collect()
    started = time.monotonic()
    spacebender.run_spacebender(input_file=str(src), output_file=str(dst), **SETTINGS)
    if not dst.is_file():
        raise RuntimeError("SpaceBender returned without creating the output file.")
    logging.info("Finished in %.1f minutes: %s", (time.monotonic() - started) / 60, dst)


if __name__ == "__main__":
    try:
        main()
    except Exception:
        logging.exception("SpaceBender run failed")
        sys.exit(1)
