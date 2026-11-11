"""Query AlphaGenome 2 conditions or run one public-reference smoke prediction.

Inputs: ALPHAGENOME_API_KEY environment variable; mode and output directory.
Outputs: conditions.tsv or RNA-seq values.npy and metadata.tsv; receipt.json.
Example: python scripts/alphagenome2.py conditions --output results
No custom sequences or cell data are accepted by this bounded helper.
"""
from __future__ import annotations

import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import sys


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Validate CLI before importing the SDK or accessing credentials."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['conditions', 'smoke'])
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--organism', choices=['human', 'mouse'], default='human')
    parser.add_argument('--timeout', type=float, default=60)
    parser.add_argument('--confirm-noncommercial', action='store_true',
                        help='Confirm eligibility and current API terms before network access.')
    args = parser.parse_args(argv)
    if not 0 < args.timeout <= 300:
        parser.error('--timeout must be >0 and <=300 seconds')
    if args.mode == 'smoke' and args.organism != 'human':
        parser.error('The fixed public smoke example is human only.')
    return args


def run(args: argparse.Namespace) -> None:
    """Make one explicitly requested RPC and write outputs only on success."""
    if not args.confirm_noncommercial:
        sys.exit('Confirm noncommercial eligibility with --confirm-noncommercial; review API terms.')
    key = os.environ.get('ALPHAGENOME_API_KEY', '').strip()
    if not key:
        sys.exit('Set ALPHAGENOME_API_KEY using your own key; do not pass keys as CLI arguments.')
    if args.output.exists() and any(args.output.iterdir()):
        sys.exit('Output directory must be new or empty; refusing to overwrite results.')
    try:
        from alphagenome.models.v2 import dna_client
        from alphagenome.data import genome
        import numpy as np
        import grpc
    except ImportError:
        sys.exit('Install with: python -m pip install -r requirements.txt (AlphaGenome 2 SDK required)')
    receipt = {'sdk': importlib.metadata.version('alphagenome'),
               'python': platform.python_version(), 'api_namespace': 'alphagenome.models.v2',
               'mode': args.mode, 'organism': args.organism,
               'model_selection': 'SDK server default; no immutable model ID asserted'}
    model = dna_client.create(key, timeout=args.timeout)
    organism = (dna_client.Organism.HOMO_SAPIENS if args.organism == 'human'
                else dna_client.Organism.MUS_MUSCULUS)
    try:
        if args.mode == 'conditions':
            conditions = model.supported_conditions(organism)
            args.output.mkdir(parents=True, exist_ok=True)
            conditions.to_csv(args.output / 'conditions.tsv', sep='\t', index=False)
            receipt.update(metadata_access=True, prediction_executed=False, rows=len(conditions))
        else:
            # Public tutorial locus GLB1L, resized to the smallest supported input.
            variant = genome.Variant(chromosome='chr2', position=219243701,
                                     reference_bases='C', alternate_bases='T')
            interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_16KB)
            result = model.predict_interval(interval, organism=organism,
                      requested_outputs=[dna_client.OutputType.RNA_SEQ],
                      condition='UBERON:0000178_CL:0001054')
            tracks = result.rna_seq
            if tracks is None:
                sys.exit('RPC returned no RNA-seq tracks; prediction not verified.')
            values = tracks.values
            if values.size == 0 or not np.isfinite(values).all():
                sys.exit('Empty/nonfinite tracks; prediction not verified.')
            args.output.mkdir(parents=True, exist_ok=True)
            np.save(args.output / 'rna_seq.npy', values, allow_pickle=False)
            tracks.metadata.to_csv(args.output / 'metadata.tsv', sep='\t', index=False)
            receipt.update(prediction_executed=True, shape=list(values.shape),
                           interval=str(interval), assembly='GRCh38',
                           condition='UBERON:0000178_CL:0001054',
                           interpretation='Technical smoke only; no variant effect computed')
    except grpc.RpcError as exc:
        # Server error text can echo request details: disclose code only, never credentials.
        code = exc.code().name
        sys.exit(f'AlphaGenome RPC failed: {code}. No automatic retry; check terms/auth/quota. '
                 'For RESOURCE_EXHAUSTED wait for provider Retry-After/retry guidance before a manual retry.')
    (args.output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(f'✓ AlphaGenome 2 {args.mode} completed: {args.output}')


def main() -> None:
    """CLI entry point; avoid printing potentially sensitive exception payloads."""
    args = parse_args()
    try:
        run(args)
    except (OSError, ValueError, TypeError, AttributeError) as exc:
        sys.exit(f'Local/SDK failure ({type(exc).__name__}); inspect inputs, permissions and SDK version.')


if __name__ == '__main__':
    main()
