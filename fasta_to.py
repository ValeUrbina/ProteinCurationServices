import sys
import os

# Convert file format from fast to PFAM CURATIONS format


def fasta_to_stockholm(source_fasta_alignment_path, target_stockholm_alignment_path):
    source_fasta_alignment_file = open(source_fasta_alignment_path, 'r')
    target_stockholm_alignment_file = open(
        target_stockholm_alignment_path, 'w')

    len_intro_stockholm = 22
    fasta_lines = source_fasta_alignment_file.readlines()
    seq = ''
    for i, line in enumerate(fasta_lines):
        if len(line) == 0:
            break
        if line[0] == '>':
            if i != 0:
                target_stockholm_alignment_file.write(seq+'\n')
            seq = line[1:-1] + (' ') * (len_intro_stockholm - len(line[1:-1]))
        else:
            seq += line[:-1]
    source_fasta_alignment_file.close()
    target_stockholm_alignment_file.close()


def main(source_fasta_alignment_path, target_stockholm_alignment_path):
    fasta_to_stockholm(source_fasta_alignment_path,
                       target_stockholm_alignment_path)
    return


# Execution: fasta_to_stockholm.py input_file_path(fasta format) output_file_path
# Example: fasta_to_stockholm.py  PF000/fast_alignment PF000/stockholm_alignment
if __name__ == "__main__":
    main()
