#!/usr/bin/env python3
"""Correct one Visium HD sample with all genes using SPARKLE.

Usage: python -u scripts/run_sparkle.py SPACE_RANGER_OUTS NEW_OUTPUT_DIR
Use --validate-only first: reconstruct cell counts and check Space Ranger equality.
Only directories containing diagnostics.json are completed, importable results.
Each output directory must be new. Validation outputs cannot be reused for a run.
The input includes all segmented capture units and in-tissue background units.
All available genes are fitted (n_high_genes=None); absent zero genes are retained
in the export. Outputs include original/corrected Matrix Market counts, gene and
cell identifiers, gene fit diagnostics, package versions, timings and peak RAM.
"""
import argparse
import importlib.metadata
import json
import logging
import platform
import resource
import sys
import time
from pathlib import Path

from stambient import SPARKLE

from betareg_sparkle_io import (
    make_json_serializable,
    prepare_betareg_sparkle_input,
    validate_betareg_sparkle_counts,
    write_betareg_sparkle_export,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('outs', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--validate-only', action='store_true')
    parser.add_argument('--bin-size', type=float, default=50.0, help='Background bin side, micrometres')
    parser.add_argument('--max-radius', type=float, default=200.0, help='Neighbour radius, micrometres')
    parser.add_argument('--lambda-grid', type=float, nargs='+', default=[10, 20, 30, 50, 70, 100, 150, 200, 300])
    parser.add_argument('--r2-threshold', type=float, default=0.01)
    parser.add_argument('--n-lambda-genes', type=int, default=50)
    parser.add_argument('--gpu', action='store_true', help='Request NVIDIA CUDA; default CPU')
    args = parser.parse_args()
    params = dict(bin_size=args.bin_size, max_radius=args.max_radius,
                  lambda_grid=args.lambda_grid, r2_threshold=args.r2_threshold,
                  n_lambda_genes=args.n_lambda_genes, n_high_genes=None,
                  use_gpu=args.gpu, self_confidence_penalty=True,
                  use_expr_weight=False, verbose=True)
    model = SPARKLE(**params)
    args.output.mkdir(parents=True, exist_ok=False)
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s',
                        handlers=[logging.StreamHandler(), logging.FileHandler(args.output / 'run.log')])
    started = time.perf_counter()
    try:
        run_input = prepare_betareg_sparkle_input(args.outs)
        load_seconds = time.perf_counter() - started
        if args.validate_only:
            logging.info('Validation passed in %.1f seconds. No correction run.', load_seconds)
            return
        corrected, diagnostics = model.fit_transform(
            run_input.sparkle_data['spot_expr'],
            run_input.sparkle_data['spot_coords'],
            run_input.sparkle_data['spot_labels'],
        )
        validate_betareg_sparkle_counts(corrected, run_input)
        # Fine-resolution input is no longer needed while the dense corrected
        # matrix is converted and written.
        run_input.sparkle_data = None
        write_betareg_sparkle_export(
            args.output, args.outs, corrected, run_input, diagnostics,
            args.r2_threshold,
        )
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        diagnostics.update(
            schema_version=1, status='complete', input_outs=str(args.outs.resolve()),
            parameters=params, pixel_size_um=run_input.pixel_size_um,
            original_counts_verified=True,
            export_shape=list(run_input.original_counts.shape),
            loading_validation_seconds=load_seconds,
            elapsed_seconds=time.perf_counter() - started,
            peak_rss_gib=peak / (1024**3 if sys.platform == 'darwin' else 1024**2),
            versions={p: importlib.metadata.version(p) for p in ['stambient', 'numpy', 'scipy', 'h5py']},
            python=platform.python_version())
        # Completion marker written last; R rejects unfinished directories.
        with open(args.output / 'diagnostics.json', 'w') as handle:
            json.dump(make_json_serializable(diagnostics), handle, indent=2, allow_nan=False)
        logging.info('Completed in %.1f minutes; peak RSS %.2f GiB.',
                     diagnostics['elapsed_seconds'] / 60, diagnostics['peak_rss_gib'])
    except Exception:
        logging.exception('Run failed; this is not a completed export.')
        raise


if __name__ == '__main__':
    main()
