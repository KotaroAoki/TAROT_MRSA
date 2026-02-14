#!/usr/bin/env python3
"""
Calculate kappa (transition/transversion ratio) from PhyML stats output.

Reads substitution rates from PhyML stats file and computes kappa as:
  kappa = mean(transitions) / mean(transversions)
       = ((AG + CT) / 2) / ((AC + AT + CG + GT) / 4)
"""

import argparse
import re
import sys


def extract_rate(filepath, pattern):
    """Extract a substitution rate value from PhyML stats file."""
    with open(filepath, "r") as f:
        for line in f:
            if pattern in line:
                # Extract the numeric value (last whitespace-separated field)
                parts = line.strip().split()
                for part in reversed(parts):
                    try:
                        return float(part)
                    except ValueError:
                        continue
    print(f"Error: Pattern '{pattern}' not found in {filepath}", file=sys.stderr)
    sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Calculate kappa from PhyML stats output"
    )
    parser.add_argument(
        "-i", "--input", required=True, help="PhyML stats file (*_phyml_stats.txt)"
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Output file to write kappa value"
    )
    args = parser.parse_args()

    # Extract substitution rates
    ag = extract_rate(args.input, "A <-> G")
    ct = extract_rate(args.input, "C <-> T")
    ac = extract_rate(args.input, "A <-> C")
    at = extract_rate(args.input, "A <-> T")
    cg = extract_rate(args.input, "C <-> G")
    gt = extract_rate(args.input, "G <-> T")

    # Calculate kappa
    transitions = (ag + ct) / 2.0
    transversions = (ac + at + cg + gt) / 4.0

    if transversions == 0:
        print("Error: Transversion rate is zero, cannot compute kappa.", file=sys.stderr)
        sys.exit(1)

    kappa = transitions / transversions

    # Write kappa to output file
    with open(args.output, "w") as f:
        f.write(f"{kappa:.5f}\n")

    print(f"Kappa = {kappa:.5f}", file=sys.stderr)


if __name__ == "__main__":
    main()
