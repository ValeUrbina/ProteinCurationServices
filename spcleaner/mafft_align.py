
import sys
import os

# Use mafft to align a fasta alignment
# https://www.ebi.ac.uk/Tools/msa/mafft/
# https://mafft.cbrc.jp/alignment/software/


def mafft_align(fasta_path, aligned_fasta_path):
    myCmd = 'mafft --auto ' + fasta_path + ' > '+aligned_fasta_path
    os.system(myCmd)


def main():
    fasta_path = sys.argv[1]
    aligned_fasta_path = sys.argv[2]

    mafft_align(fasta_path, aligned_fasta_path)
    return

# Execution: mafft_align.py input_file_path output_file_path
# Example: mafft_align.py  PF000/alignment PF000/trimmed_alignment


if __name__ == "__main__":
    main()
