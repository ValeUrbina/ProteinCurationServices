import sys
import os
import to_fasta

# delete fragments
# output: new align
# this file should be located in the directory specified by pfam_curation_tools docker

# file_name = ALIGN
# outputfile_name = delfragmentsALIGN


def is_fragment(seq):
    aux = seq.split()[1].strip('.')
    return aux[0] == '-' or aux[-1] == '-'


def main(PATH, directory, pfam_code, file_name, outputfile_name):
    # Nos dirigimos a la carpeta donde se encuentra el archivo
    input_path = os.path.join(
        PATH, 'pfam_data', directory, pfam_code, file_name)
    output_path = os.path.join(
        PATH, 'pfam_data', directory, pfam_code, outputfile_name)
    input_file = open(input_path, 'r')
    output_file = open(output_path, 'w')
    lines = input_file.readlines()
    for i, line in enumerate(lines):
        if (is_fragment(line)):
            continue
        else:
            output_file.write(line)

    input_file.close()
    output_file.close()

    fasta_path = os.path.join(
        PATH, 'pfam_data', directory, pfam_code, outputfile_name + '.fasta')
    to_fasta.main(output_path, fasta_path)

    return fasta_path


if __name__ == "__main__":
    main()
