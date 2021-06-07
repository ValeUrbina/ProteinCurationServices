import sys
import os

# return pfamout file in specific segments
# this file should be located in the directory specified by pfam_curation_tools docker

PATH = '/home/valeria/Documentos/Tesis_2/Docker/pfam_curation/'


def main(directory, pfam_code, evalue):
    # Nos dirigimos a la carpeta donde se encuentra el archivo
    output_path = os.path.join(
        PATH, 'pfam_data', directory, pfam_code, 'PFAMOUT')
    output_file = open(output_path, 'r')
    lines = output_file.readlines()
    sequences = lines[7]
    cutoffvalues = lines[7]
    i_domain = None
    i_cutoff = None
    for i, line in enumerate(lines[9:]):
        if line[0] == "#":
            i_domain = i + 9
            break
        else:
            sequences = sequences + line
            read_evalue = float(line.split()[-4])
            if i_cutoff == None and evalue < read_evalue:
                i_cutoff = i + 9

    cutoffvalues = cutoffvalues + lines[i_cutoff - 2] + lines[i_cutoff -
                                                              1] + lines[i_cutoff] + lines[i_cutoff + 1]

    domains = lines[i_domain + 4]
    for i, line in enumerate(lines[i_domain + 7:]):
        domains = domains + line

    output_file.close()
    return {"cutoffvalues": cutoffvalues, "sequences": sequences, "domains": domains}


if __name__ == "__main__":
    main()
