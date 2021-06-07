
import sys
import os
import pandas as pd
import numpy as np
import json
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import AlignInfo
from tqdm import tqdm

HEADER_PATH = '/home/valeria/Documentos/Tesis_2/spcleaner/curation/header.txt'


def cut_sequences(lines, final_idx_pfam_code):
    intro_len = 24
    new_lines = ''

    max_len = 0
    len_list = []

    for line in lines[:-1]:
        seq = line[24:]
        for beg, end in final_idx_pfam_code:
            short_seq = seq[beg:end+1]
            short_seq = short_seq.replace(
                '-', '').replace('.', '').replace('\n', '')
            len_list.append(len(short_seq))

    len_list = np.array(len_list)
    min_range = int(np.percentile(len_list, 20))
    max_range = max_len = int(np.percentile(len_list, 80))

    j = -1
    for i, line in enumerate(lines[:-1]):
        intro = line[:24]
        uniprot_code = intro.split(' ')[0].split('/')[0]
        idxs = [int(a) for a in intro.split(' ')[0].split('/')[1].split('-')]
        first_idx = idxs[0]
        curr_idx = first_idx
        seq = line[24:]

        for beg, end in final_idx_pfam_code:
            j += 1
            if len_list[j] < min_range or len_list[j] > max_range:
                continue
            short_seq = seq[beg:end+1]
            short_seq = short_seq.replace(
                '-', '').replace('.', '').replace('\n', '')
            if len(short_seq) < 5:
                continue
            new_idxs = str(curr_idx) + '-' + str(curr_idx + len(short_seq)-1)
            curr_idx = curr_idx + len(short_seq)
            short_seq += ('.')*(max_len-len(short_seq))

            new_intro = uniprot_code + '/' + new_idxs
            new_intro += (' ')*(intro_len-len(new_intro))

            new_lines += new_intro + short_seq + '\n'
    return new_lines


def create_file_with_header(input_file_path, file_with_header_path):
    with open(file_with_header_path, 'w') as new:
        with open(HEADER_PATH) as prefix:
            new.write(prefix.read())
        with open(input_file_path) as old:  # alineado_stockholm
            new.write(old.read())


def get_uniprot_codes_from_pfam_code(pfam_code):
    # TODO(valeria) Implement get_uniprot_codes_from_pfam_code
    # this was my code but it should be dynamic
    # or you could update this excel daily and keep it that way
    df_codes = pd.read_csv("./pdb_unprot_pfam_reviewed_alfasolenoids.csv")
    uniprot_codes = list(
        df_codes[df_codes['PFAM_ID'] == pfam_code]['SP_PRIMARY'].unique())
    return uniprot_codes
    # return ['Q72L02', 'Q2W8Q0', 'P38825']


def get_pdb_chains(uniprot_found, pfam_code):
    # TODO(valeria) Implement get_pdb_chains
    # Same as get_uniprot_codes_from_pfam_code
    # this was my code but it should be dynamic
    # or you could update this excel daily and keep it that way
    df_codes = pd.read_csv("./pdb_unprot_pfam_reviewed_alfasolenoids.csv")
    df_codes['PDB_chain'] = df_codes['PDB']+df_codes['CHAIN']
    df = df_codes[(df_codes['SP_PRIMARY'] == uniprot_found) &
                  (df_codes['PFAM_ID'] == pfam_code)].copy()
    pdb_chains = list(df['PDB_chain'].unique())
    return pdb_chains


def get_pdb_to_uniprot_map_df(uniprot_code):
    # TODO(valeria) get_pdb_to_uniprot_map_df get_pdb_chains
    # Same as get_uniprot_codes_from_pfam_code
    # this was my code but it should be dynamic
    # or you could update this excel daily and keep it that way
    pdb_to_uniprot = pd.read_csv('./pdb2uniprot.csv')
    pdb_to_uniprot['pdb position'] = pd.to_numeric(
        pdb_to_uniprot['pdb position'])
    pdb_to_uniprot = pdb_to_uniprot[pdb_to_uniprot['uniprot id']
                                    == uniprot_code]
    return pdb_to_uniprot


def get_units_info(pfam_code, input_file_path, pdb_directory_path):

    result = None
    file_with_header_path = f'{input_file_path}_with_header'
    create_file_with_header(input_file_path, file_with_header_path)

    align = AlignIO.read(file_with_header_path, "stockholm")
    uniprot_codes = get_uniprot_codes_from_pfam_code(
        pfam_code)  # TODO(valeria)

    print('uniprot_codes: ', uniprot_codes)
    # Find uniprot
    idx_found_list = []
    uniprot_found_list = []
    full_seq_list = []
    for i, al in enumerate(align):
        uniprot_in_align = al.id.split('.')[0]
        if uniprot_in_align in uniprot_codes:
            uniprot_found_list.append(uniprot_in_align)
            idx_found_list.append(i)
            full_seq_list.append(al)

    print('uniprot_found_list: ', uniprot_found_list)
    print('idx_found_list: ', idx_found_list)

    # Check if exists
    if idx_found_list == []:
        raise Exception('Pfam code without uniprots in alignmnet')

    units_uniprot = []
    uniprot_found = None
    # for each uniprot found in pfam alignment
    for i, idx_found in enumerate(idx_found_list):
        uniprot_found = uniprot_found_list[i]
        full_seq = full_seq_list[i]

        pdb_chains = get_pdb_chains(uniprot_found, pfam_code)

        # indices de inicio y fin en align
        seq_beg_idx, seq_end_idx = [
            int(a) for a in full_seq.id.split('/')[1].split('-')]
        pdb_to_uniprot_map_df = get_pdb_to_uniprot_map_df(uniprot_found)
        for chosen_pdb_chain in pdb_chains:
            # repeat_units index
            f = open(pdb_directory_path + chosen_pdb_chain + '.db')
            lines = f.readlines()  # each line represents a repeat unit
            units_uniprot = []

            pdb_to_uniprot_index = pdb_to_uniprot_map_df[(pdb_to_uniprot_map_df['pdb id'] == chosen_pdb_chain[:-1])
                                                         & (pdb_to_uniprot_map_df['pdb chain'] == chosen_pdb_chain[-1])]

            # getting the possition of each repeat unit in uniprot format
            for line in lines:
                if line[:4] != 'UNIT':
                    continue
                aux_db = [int(a) for a in line[5:-1].split(' ')]
                try:
                    # this is a tuple [being, end] of repeat unit
                    aux_uniprot = [pdb_to_uniprot_index[pdb_to_uniprot_index['pdb position']
                                                        == a].iloc[0]['uniprot position'] for a in aux_db]
                except:
                    #print("Error al matchear indices pdb a uniprot")
                    continue

                # if its inside our alignment add it
                if aux_uniprot[0] >= seq_beg_idx and aux_uniprot[1] <= seq_end_idx:
                    units_uniprot.append(aux_uniprot)

            if len(units_uniprot) > 0:
                break
        if len(units_uniprot) > 0:
            break

    print('units_uniprot: ', units_uniprot)
    print('uniprot_found: ', uniprot_found)

    if len(units_uniprot) > 0:
        result = {'uniprot': uniprot_found,
                  'units': units_uniprot,
                  'idx': idx_found,
                  'full_seq': full_seq}

    return result

# Our alignments have '-' and '.' and our repeat unit position idxs should consider them:


def get_final_idxs(units_info):
    #final_idx = {}
    result = units_info
    # for pfam_code in tqdm(result):
    #print('pfam_code', pfam_code, result[pfam_code])
    print('result', result['uniprot'])
    uniprot_code = result['uniprot']
    rec = result['full_seq']
    seq = rec.seq
    seq_beg_idx, seq_end_idx = [int(a)
                                for a in rec.id.split('/')[1].split('-')]
    uniprot_units = result['units']
    seq_unit = []

    idx = seq_beg_idx - 1
    pos_idx = 0
    x = 0
    y = 0
    final_units = []
    for i, s in enumerate(seq):
        if s in ['-', '.']:
            continue
        idx += 1
        if idx == uniprot_units[len(final_units)][pos_idx]:
            #print("encontro un match ",idx,pos_idx, i, s)
            if pos_idx == 0:
                x = i
                pos_idx = 1
            else:
                y = i
                pos_idx = 0
                final_units.append([x, y])
        if len(final_units) == len(uniprot_units):
            break
    #final_idx[pfam_code] = final_units

    return final_units


def create_new_alignments(final_idx, pfam_code, pfam_directory_path):
    # for pfam_code in tqdm(final_idx):
    file = open(pfam_directory_path + pfam_code+'/ALIGN', 'r')
    new_file = open(pfam_directory_path + pfam_code+'/ALIGN_repeat', 'w')
    lines = file.readlines()
    new_lines = cut_sequences(lines, final_idx)
    new_lines = ('').join(new_lines)
    new_file.write(new_lines)
    file.close()
    new_file.close()

# Divides alignment according to repeat units on pdb


def main():
    pfam_code = sys.argv[1]  # PF14559
    # alignment file (ALIGN) /home/valeria/Documentos/Tesis_2/Docker/pfam_curation/pfam_data/PF14559/ALIGN
    input_file_path = sys.argv[2]
    # /home/valeria/Documentos/Tesis_2/spcleaner/PDB
    pdb_directory_path = sys.argv[3]
    pfam_directory_path = sys.argv[4]
    # TODO(valeria) throw exception if units_info is None
    print('pfamcode', pfam_code)
    units_info = get_units_info(pfam_code, input_file_path, pdb_directory_path)
    # TODO(valeria) throw exception if units_info is None
    final_idx = get_final_idxs(units_info)
    create_new_alignments(final_idx, pfam_code, pfam_directory_path)


# Execution: divide_alignmnets.py pfam_code input_file_path pdb_directory_path
# Example: divide_alignmnets.py PF14559 .PF14559/ALIGN ./DbfilesAll_2017_Octubre_update/
if __name__ == "__main__":
    main()
