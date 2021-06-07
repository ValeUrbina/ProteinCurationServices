from os.path import join
import sys
import os
import to_fasta

# Create a hmm and search the similar domains using hmmer
# this file should be located in the directory specified by pfam_curation_tools docker

SCRIPT = "./exec_hmmersearch.sh {} {} {} {}"
PATH = '/home/valeria/Documentos/Tesis_2/Docker/pfam_curation/'
pfamseq_path = '/home/valeria/Documentos/Tesis_2/Docker/pfam_curation/seqlib/pfamseq'


def main(directory, pfam_code):
    # Nos dirigimos a la carpeta donde se levantará el docker
    os.chdir(PATH + '/Scripts')
    # Se ejecuta hmmer
    project_path = os.path.join(
        PATH, 'pfam_data', directory, pfam_code)
    # align_path = os.path.join(PATH, 'pfam_data', directory, pfam_code, 'ALIGN')
    #hmmprofile_path = os.path.join(PATH, 'pfam_data', directory, pfam_code, 'hmmALIGN')
    my_cmd = SCRIPT.format(project_path, 'hmmALIGN', 'ALIGN', pfamseq_path)
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
