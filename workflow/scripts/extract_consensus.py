#!/usr/bin/env python3
"""
Extract core genome consensus sequences from VarScan mpileup2cns output.

Replaces the complex awk/perl pipeline in the original bash script (lines 130, 133).
Parses VarScan output, filters positions, extracts consensus bases per sample,
and produces transposed text and original matrix files.
"""

import argparse
import re
import sys


def is_standard_base(base):
    """Check if a base is a standard nucleotide (A, T, G, C)."""
    return base.upper() in {"A", "T", "G", "C"}


def parse_varscan_output(input_file, num_samples):
    """
    Parse VarScan mpileup2cns output and extract consensus bases.

    Filters:
    - Lines containing "N:" are excluded
    - Lines containing "/" are excluded
    - Only lines containing "Pass" are kept
    - Only standard bases (A, T, G, C) are retained

    Returns:
        consensus_matrix: list of lists, each inner list is one position's bases
                          across all samples
    """
    consensus_matrix = []

    with open(input_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # Apply filters
            if "N:" in line:
                continue
            if "/" in line:
                continue
            if "Pass" not in line:
                continue

            fields = line.split("\t")
            # Sample data starts at field index 10 (0-based)
            # Each sample has 6 fields in VarScan output:
            # Cons:Cov:Reads1:Reads2:Freq:P-value
            sample_fields = fields[10:]

            bases = []
            valid = True
            for i in range(num_samples):
                idx = i * 6  # Each sample occupies 6 fields
                if idx >= len(sample_fields):
                    valid = False
                    break
                sample_data = sample_fields[idx]
                # The consensus base is before the first ":"
                base = sample_data.split(":")[0] if ":" in sample_data else sample_data
                if not is_standard_base(base):
                    valid = False
                    break
                bases.append(base.upper())

            if valid and len(bases) == num_samples:
                consensus_matrix.append(bases)

    return consensus_matrix


def transpose_matrix(matrix):
    """Transpose a 2D matrix (list of lists)."""
    if not matrix:
        return []
    return list(map(list, zip(*matrix)))


def main():
    parser = argparse.ArgumentParser(
        description="Extract core genome consensus from VarScan output"
    )
    parser.add_argument(
        "-i", "--input", required=True, help="VarScan mpileup2cns output file"
    )
    parser.add_argument(
        "-n", "--name-list", required=True, help="Name list file (one >name per line)"
    )
    parser.add_argument(
        "-o",
        "--output-prefix",
        required=True,
        help="Output prefix for generated files",
    )
    args = parser.parse_args()

    # Read name list to determine number of samples
    with open(args.name_list, "r") as f:
        names = [line.strip() for line in f if line.strip()]
    num_samples = len(names)

    if num_samples == 0:
        print("Error: No sample names found in name list.", file=sys.stderr)
        sys.exit(1)

    # Parse VarScan output
    consensus_matrix = parse_varscan_output(args.input, num_samples)

    if not consensus_matrix:
        print("Warning: No valid positions found.", file=sys.stderr)
        # Create empty output files
        open(args.output_prefix + "_tran.fasta", "w").close()
        open(args.output_prefix + "_original.txt", "w").close()
        sys.exit(0)

    # Transpose: rows become samples, columns become positions
    transposed = transpose_matrix(consensus_matrix)

    # Write transposed FASTA (core genome multifasta)
    with open(args.output_prefix + "_tran.fasta", "w") as f:
        for name, seq_bases in zip(names, transposed):
            f.write(name + "\n")
            f.write("".join(seq_bases) + "\n")

    # Write original matrix (positions x samples, tab-separated)
    with open(args.output_prefix + "_original.txt", "w") as f:
        for position_bases in consensus_matrix:
            f.write("\t".join(position_bases) + "\n")

    print(
        f"Processed {len(consensus_matrix)} positions for {num_samples} samples.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
