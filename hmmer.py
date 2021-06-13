from os.path import join
import sys
import os
import to_fasta

# Create a hmm and search the similar domains using hmmer
# this file should be located in the directory specified by pfam_curation_tools docker

SCRIPT1 = "hmmbuild {} {}"
SCRIPT2 = "hmmsearch --tblout tblout.txt --domtblout domtblout.txt --pfamtblout pfamtblout.txt hmmALIGN {}"

# file_name = file_name
# outputfile_name = hmmALIGN, domtblout.txt


def main(PATH, pfamseq_path, directory, pfam_code, file_name):
    # Nos dirigimos a la carpeta donde se va a ejecutar hmmer
    project_path = os.path.join(
        PATH, 'pfam_data', directory, pfam_code)
    os.chdir(project_path)
    # Se ejecuta hmmer
    target_path = os.path.join(
        PATH, 'pfam_data', directory, pfam_code, 'hmmALIGN')
    source_path = os.path.join(
        PATH, 'pfam_data', directory, pfam_code, file_name)
    my_cmd_1 = SCRIPT1.format(target_path, source_path)
    my_cmd_2 = SCRIPT2.format(pfamseq_path)
    os.system(my_cmd_1)
    os.system(my_cmd_2)
    dom_path = os.path.join(
        PATH, 'pfam_data', directory, pfam_code, 'domtblout.txt')
    # Se retorna la ruta del archivo
    if os.path.exists(dom_path):
        return dom_path
    return "There's no such domtblout.txt file"


# Execution: python3 pfam_download.py pfam_code
# Example: python3 pfam_download.py PF00023
if __name__ == "__main__":
    main()
