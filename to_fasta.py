import sys
import os
import pandas as pd
import numpy as np
import json

# Convert file format from fast to FASTA format


def stockholm_to_fasta(source_stockholm_alignment_path, target_fasta_alignment_path):
    source_stockholm_alignment_file = open(
        source_stockholm_alignment_path, 'r')
    target_fasta_alignment_file = open(target_fasta_alignment_path, 'w')

    stockholm_lines = source_stockholm_alignment_file.readlines()
    for i, line in enumerate(stockholm_lines):
        if i > 0:
            target_fasta_alignment_file.write('\n')
        intro = line.split()[0].strip()
        if len(intro) == 0:
            continue
        intro = '>' + intro
        seq = line.split()[1].strip()
        target_fasta_alignment_file.write(intro + '\n' + seq)
    source_stockholm_alignment_file.close()
    target_fasta_alignment_file.close()


def main(source_stockholm_alignment_path, target_fasta_alignment_path):
    stockholm_to_fasta(source_stockholm_alignment_path,
                       target_fasta_alignment_path)

    return 1


if __name__ == "__main__":
    main()
