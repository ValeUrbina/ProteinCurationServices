from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import AlignInfo
import numpy as np
import pandas as pd
import sys
import time
import os
from tqdm import tqdm

HEADER_PATH = 'spcleaner/header.txt'
PATH = '/home/valeria/Documentos/Tesis_2/Docker/pfam_curation/'

sorted_amino = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
                'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z']

BLOSUM62 = [[4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0, -4],
            [-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -
                1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1, -4],
            [-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -
                2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1, -4],
            [-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -
                3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1, -4],
            [0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -
                1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4],
            [-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,
                0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1, -4],
            [-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -
                2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4],
            [0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -
                3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -4],
            [-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -
                2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1, -4],
            [-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,
                1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1, -4],
            [-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,
                2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1, -4],
            [-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -
                1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1, -4],
            [-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,
                5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1, -4],
            [-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,
                0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1, -4],
            [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -
                2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2, -4],
            [1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -
                1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0, -4],
            [0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -
                1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0, -4],
            [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -
                1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2, -4],
            [-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -
                1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1, -4],
            [0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,
                1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1, -4],
            [-2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -
                3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1, -4],
            [-1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -
                1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4],
            [0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -
                1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1, -4],
            [-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -1]]


def insert_sort(lista, n):
    index = 0
    if len(lista) == 0:
        lista.append(n)
        return lista
    # Searching for the position
    for i in range(len(lista)):
        if lista[i][0] < n[0]:
            index = i
            break
        index = i + 1
    # Inserting n in the list
    lista = lista[:index] + [n] + lista[index:]
    return lista


def blossum_results(aminoacids_count, n_seq):
    max_level = -100.0
    amino_list = []
    for j, aminoacid_j in enumerate(sorted_amino):
        if aminoacid_j not in aminoacids_count.keys():
            aminoacids_count[aminoacid_j] = 0

    for j, aminoacid_j in enumerate(sorted_amino):
        # Convert counts to similarity counts
        simCount = 0
        for k, aminoacid_k in enumerate(sorted_amino):
            if (aminoacid_j == aminoacid_k):
                simCount += (aminoacids_count[aminoacid_j] - 1) * \
                    aminoacids_count[aminoacid_k] * BLOSUM62[j][k]
            else:
                simCount += aminoacids_count[aminoacid_j] * \
                    aminoacids_count[aminoacid_k] * BLOSUM62[j][k]

        if n_seq < 2:
            level = 0.0
        else:
            level = simCount / (n_seq * (n_seq - 1))

        aux_list = []
        if level > max_level:
            # print(aminoacid_j)
            # print("simCount",simCount,level,n_seq)
            for l, aminoacid_l in enumerate(sorted_amino):
                if (BLOSUM62[j][l] > 0):
                    aux_list.append(aminoacid_l)
            max_level = level
            amino_list = aux_list

    return max_level, amino_list


def add_to_delete_rows(delete_rows, seq_col):
    for i, amino in enumerate(seq_col):
        if amino != '-':
            delete_rows.add(i)
    return

# Bio needs a header on files to consider them stockholm format


def create_file_with_header(input_file_path, file_with_header_path):
    with open(file_with_header_path, 'w') as new:
        with open(HEADER_PATH) as prefix:
            new.write(prefix.read())
        with open(input_file_path) as old:  # alineado_stockholm
            new.write(old.read())


def process_curation(n_seq, len_seq, align, uniprot_codes, empty_percentage, blossum_level):
    important_columns = []
    delete_rows = set()
    delete_cols = []
    summary_align = AlignInfo.SummaryInfo(align)
    for col in tqdm(range(len_seq)):
        aminoacids_count = summary_align.pos_specific_score_matrix()[col]
        if '-' not in aminoacids_count:
            aminoacids_count['-'] = 0
        if aminoacids_count['-'] / n_seq > empty_percentage:
            delete_cols.append(col)
            add_to_delete_rows(delete_rows, summary_align.get_column(col))
        result, important_amino = blossum_results(aminoacids_count, n_seq)
        if result > blossum_level:
            # print(col,result,important_amino)
            important_columns = insert_sort(
                important_columns, (result, col, important_amino))
    return important_columns, delete_rows, delete_cols


def curate_alignment(pfam_code, uniprot_codes, file_directory,
                     file_input_name, file_output_name, vertical_threshold,
                     empty_percentage, blossum_level):
    input_file_path = f"{file_directory}/{file_input_name}"
    file_with_header_path = f"{file_directory}/{file_input_name}_with_header"
    output_file_path = f"{file_directory}/{file_output_name}"

    create_file_with_header(input_file_path, file_with_header_path)

    align = AlignIO.read(file_with_header_path, "stockholm")
    len_seq = len(align._records[0])
    n_seq = len(align._records)

    important_columns, delete_rows, delete_cols = process_curation(n_seq, len_seq, align, uniprot_codes,
                                                                   empty_percentage, blossum_level)
    if len(important_columns) == 0:
        return 0
    new_file = open(output_file_path, 'w')
    curr_file = open(input_file_path, 'r')

    curr_seq = curr_file.readlines()
    new_n_seq = 0

    for i in range(n_seq):
        if i in delete_rows:
            #print("no considerada",i)
            continue
        count = 0
        seq = align._records[i]
        for (result, col, important_amino) in important_columns:
            count += int(seq[col] in important_amino)
        count /= len(important_columns)
        if count > vertical_threshold or (align._records[i].name.split('.')[0] in uniprot_codes):
            new_file.write(curr_seq[i])
            new_n_seq += 1

    print("Initial number of sequences:", n_seq)
    print("New number of sequences:", new_n_seq)
    new_file.close()
    curr_file.close()
    return new_n_seq

# Curate a pfam alignment


def get_uniprot_codes_from_pfam_code(pfam_code):
    # TODO(valeria) Implement get_uniprot_codes_from_pfam_code
    # this was my code but it should be dynamic
    # or you could update this excel daily and keep it that way
    # return ['Q72L02', 'Q2W8Q0', 'P38825']
    df_codes = pd.read_csv(
        'spcleaner/pdb_unprot_pfam_reviewed_alfasolenoids.csv')
    uniprot_codes = list(
        df_codes[df_codes['PFAM_ID'] == pfam_code]['SP_PRIMARY'].unique())
    return uniprot_codes

# Curate a pfam alignment


def main(pfam_code, file_directory, file_input_name, file_output_name):
    #pfam_code = sys.argv[1]
    #file_directory = sys.argv[2]
    #file_input_name = sys.argv[3]
    #file_output_name = sys.argv[4]
    # vertical_threshold = sys.argv[5]
    # empty_percentage = sys.argv[6]
    # blossum_level = sys.argv[7]
    uniprot_codes = get_uniprot_codes_from_pfam_code(pfam_code)
    vertical_threshold = 0.65
    empty_percentage = 0.70
    blossum_level = 0.50
    curate_alignment(pfam_code, uniprot_codes, file_directory, file_input_name, file_output_name,
                     vertical_threshold, empty_percentage, blossum_level)
    output_path = os.path.join(file_directory, file_output_name)
    return output_path


# Execution: alignment_curator.py pfam_code pfam_file_directory file_input_name(stockholm format) file_output_name(sotckholm format)
# Example: alignments_curator.py PF14559 ../pfam/PF14559/ PF14559_stockholm_alignment PF14559_stockholm_alignment_curated
if __name__ == "__main__":
    main()
