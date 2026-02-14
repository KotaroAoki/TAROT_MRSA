#!/usr/bin/env python3
"""
Convert FASTA to Relaxed Phylip format.
Uses Biopython to handle conversion.
"""

import argparse
from Bio import SeqIO
import sys


def main():
    parser = argparse.ArgumentParser(description="Convert FASTA to Phylip")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output Phylip file")
    args = parser.parse_args()

    count = SeqIO.convert(args.input, "fasta", args.output, "phylip-relaxed")
    print(f"Converted {count} sequences to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
