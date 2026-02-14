#!/usr/bin/env python3
"""
Extract SNP positions from the consensus matrix.

Replaces the eval-based dynamic awk processing and diff-based SNP extraction
from the original bash script (lines 134-142).

Identifies positions where not all samples have the same base (i.e., SNP sites),
and generates a multifasta with only those positions.
"""

import argparse
import sys


def main():
    parser = argparse.ArgumentParser(
        description="Extract SNP positions from consensus matrix"
    )
    parser.add_argument(
        "-i", "--input", required=True, help="Original matrix file (positions x samples, tab-separated)"
    )
    parser.add_argument(
        "-n", "--name-list", required=True, help="Name list file"
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Output FASTA file with SNP positions only"
    )
    args = parser.parse_args()

    # Read name list
    with open(args.name_list, "r") as f:
        names = [line.strip() for line in f if line.strip()]
    num_samples = len(names)

    if num_samples == 0:
        print("Error: No sample names found.", file=sys.stderr)
        sys.exit(1)

    # Read original matrix and find SNP positions
    snp_positions = []  # list of lists of bases at SNP sites
    total_positions = 0

    with open(args.input, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            bases = line.split("\t")
            total_positions += 1
            # A position is a SNP if not all bases are the same
            if len(set(bases)) > 1:
                snp_positions.append(bases)

    if not snp_positions:
        print("Warning: No SNP positions found.", file=sys.stderr)
        # Write empty FASTA
        with open(args.output, "w") as f:
            for name in names:
                f.write(name + "\n")
                f.write("\n")
        sys.exit(0)

    # Transpose: rows become samples
    transposed = list(map(list, zip(*snp_positions)))

    # Write multifasta with only SNP positions
    with open(args.output, "w") as f:
        for name, seq_bases in zip(names, transposed):
            f.write(name + "\n")
            f.write("".join(seq_bases) + "\n")

    print(
        f"Found {len(snp_positions)} SNP positions out of {total_positions} total positions.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
