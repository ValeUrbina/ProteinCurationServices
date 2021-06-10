from os.path import join
import sys
import os
import to_fasta

# Create a hmm and search the similar domains using hmmer
# this file should be located in the directory specified by pfam_curation_tools docker

SCRIPT = "./exec_hmmersearch.sh {} {} {} {}"

# file_name = file_name
# outputfile_name = hmmALIGN, domtblout.txt


def main(PATH, pfamseq_path, directory, pfam_code, file_name):
    # Nos dirigimos a la carpeta donde se levantará el docker
    os.chdir(PATH + '/Scripts')
    # Se ejecuta hmmer
    project_path = os.path.join(
        PATH, 'pfam_data', directory, pfam_code)
    my_cmd = SCRIPT.format(project_path, 'hmmALIGN', file_name, pfamseq_path)
    os.system(my_cmd)
    dom_path = os.path.join(
        PATH, 'pfam_data', directory, pfam_code, 'domtblout.txt')
    # Se retorna la ruta del archivo
    if os.path.exists(dom_path):
        return dom_path
    return {"error": "There's no such domtblout.txt file"}


# Execution: python3 pfam_download.py pfam_code
# Example: python3 pfam_download.py PF00023
if __name__ == "__main__":
    main()
