import sys
import os
import to_fasta

# delete multiple sequences
# input: list of acc numbers, output: new align
# this file should be located in the directory specified by pfam_curation_tools docker
# file_name = ALIGN
# outputfile_name = delseqALIGN


def main(PATH, directory, pfam_code, accNumbers, file_name, outputfile_name):
    # Nos dirigimos a la carpeta donde se encuentra el archivo
    input_path = os.path.join(
        PATH, 'pfam_data', directory, pfam_code, file_name)
    output_path = os.path.join(
        PATH, 'pfam_data', directory, pfam_code, outputfile_name)
    input_file = open(input_path, 'r')
    output_file = open(output_path, 'w')
    lines = input_file.readlines()
    for i, line in enumerate(lines):
        if (i in accNumbers):
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
