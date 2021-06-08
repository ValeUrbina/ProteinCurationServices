import sys
import os

# return the units, sequences and pdbs
# this file should be located in the directory specified by pfam_curation_tools docker


def filter_units(units_array):
    units_array = units_array.split()
    if units_array[0] == 'UNIT':
        return True
    else:
        return False


def main(PATH, directory, pfam_code, pdb_id):
    # Nos dirigimos a la carpeta donde se encuentra el archivo
    path_pdb = 'repeatsdb/PF03377-DBFILES/' + pdb_id + '.db'
    path_units = 'repeatsdb/PF03377-PDBS.txt'
    path_sequences = os.path.join(
        PATH, 'pfam_data', directory, pfam_code, 'ALIGN')
    output_file = open(path_pdb, 'r')
    lines = output_file.readlines()
    unit = list(filter(filter_units, lines))
    unit = list(map(lambda x: {'unit': x}, unit))
    output_file.close()

    output_file = open(path_units, 'r')
    lines = output_file.readlines()
    pdbs = list(map(lambda x: {'name': x}, lines))
    output_file.close()

    output_file = open(path_sequences, 'r')
    lines = output_file.readlines()
    sequences = list(
        map(lambda x: {'name': x.split()[0], 'seq': x.split()[1]}, lines))

    result = {'unit': unit, 'pdbs': pdbs, 'sequences': sequences}

    return result


if __name__ == "__main__":
    main()
