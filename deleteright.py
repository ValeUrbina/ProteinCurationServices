import sys
import os
from typing import Sequence
import to_fasta

# delete multiple sequences
# input: column, output: new align
# this file should be located in the directory specified by pfam_curation_tools docker
# file_name = ALIGN
# outputfile_name = delseqrightALIGN


def main(PATH, directory, pfam_code, column, file_name, outputfile_name):
    # Nos dirigimos a la carpeta donde se encuentra el archivo
    input_path = os.path.join(
        PATH, 'pfam_data', directory, pfam_code, file_name)
    output_path = os.path.join(
        PATH, 'pfam_data', directory, pfam_code, outputfile_name)
    input_file = open(input_path, 'r')
    output_file = open(output_path, 'w')
    sequences = input_file.readlines()
    for i, sequence in enumerate(sequences):
        simple_seq = sequence.split()
        index = sequence.find(simple_seq[1])
        output_file.write(sequence[:index])
        output_file.write(sequence[index+column:])

    input_file.close()
    output_file.close()

    fasta_path = os.path.join(
        PATH, 'pfam_data', directory, pfam_code, outputfile_name + '.fasta')
    to_fasta.main(output_path, fasta_path)

    return fasta_path


if __name__ == "__main__":
    main()
