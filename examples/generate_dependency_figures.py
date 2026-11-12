"""Render reproducible molecule and dependency-graph SVGs.

Run from the project root:
    python examples/generate_dependency_figures.py --output examples/figures

By default the arrows represent possible *positional* influences.
For exact configuration-dependent influence, add --exact (more costly).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from stereoapps.core import write_detection_examples


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output", type=Path, default=Path("examples/figures")
    )
    parser.add_argument("--exact", action="store_true")
    parser.add_argument(
        "--all-candidate-sites",
        action="store_true",
        help="Use all existing SMG tetrahedral descriptors instead of "
             "the RDKit-matched subset",
    )
    args = parser.parse_args()
    reports = write_detection_examples(
        args.output,
        exact=args.exact,
        matched_rdkit_sites=not args.all_candidate_sites,
    )
    print("EXAMPLE_REPORT=" + json.dumps(reports, sort_keys=True))


if __name__ == "__main__":
    main()
