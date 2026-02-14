#!/usr/bin/env python3
"""
Post-ClonalFrameML processing: extract SNPs from recombination-removed FASTA.

Replaces the complex perl/awk pipeline after ClonalFrameML (lines 167-175).
Reads the recRemoved.fasta, transposes, identifies SNP positions,
and generates a final SNP-only FASTA.
"""

import argparse
import sys


def read_fasta(filepath):
    """Read a FASTA file and return list of (name, sequence) tuples."""
    sequences = []
    current_name = None
    current_seq = []

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    sequences.append((current_name, "".join(current_seq)))
                current_name = line
                current_seq = []
            else:
                current_seq.append(line)
        if current_name is not None:
            sequences.append((current_name, "".join(current_seq)))

    return sequences


def main():
    parser = argparse.ArgumentParser(
        description="Extract SNPs from recombination-removed FASTA"
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Recombination-removed FASTA file (*.recRemoved.fasta)"
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output SNP-only FASTA file"
    )
    args = parser.parse_args()

    # Read input FASTA
    sequences = read_fasta(args.input)

    if not sequences:
        print("Error: No sequences found in input file.", file=sys.stderr)
        sys.exit(1)

    names = [s[0] for s in sequences]
    seqs = [s[1] for s in sequences]

    # Verify all sequences have the same length
    seq_len = len(seqs[0])
    for i, seq in enumerate(seqs):
        if len(seq) != seq_len:
            print(
                f"Error: Sequence {names[i]} has length {len(seq)}, "
                f"expected {seq_len}.",
                file=sys.stderr,
            )
            sys.exit(1)

    # Find SNP positions (where not all bases are the same)
    snp_positions = []
    for pos in range(seq_len):
        bases_at_pos = set(seq[pos].upper() for seq in seqs)
        if len(bases_at_pos) > 1:
            snp_positions.append(pos)

    # Extract SNP sequences
    with open(args.output, "w") as f:
        for name, seq in zip(names, seqs):
            f.write(name + "\n")
            snp_seq = "".join(seq[pos] for pos in snp_positions)
            f.write(snp_seq + "\n")

    print(
        f"Found {len(snp_positions)} SNP positions out of {seq_len} total positions "
        f"across {len(sequences)} sequences.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
