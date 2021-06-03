from os.path import join
import sys
import os
import to_fasta

# Download pfam file with pfco using pfam_curation_tools docker
# this file should be located in the directory specified by pfam_curation_tools docker

SCRIPT = "docker run --rm -it -v $(pwd)/Scripts:/Scripts -v $(pwd)/pfam_data:/home/pfam/pfam_data -v $(pwd)/pfam.conf:/home/pfam/pfam.conf -v $(pwd)/seqlib:/data/seqlib -v $(pwd)/Dictionary/dictionary:/home/pfam/Dictionary/dictionary -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix dockerhub.ebi.ac.uk/pfam/pfam-curation bash -c '/Scripts/exc_pfco.sh {} {}'"
PATH = '/home/valeria/Documentos/Tesis_2/Docker/pfam_curation/'


def main(directory, pfam_code):
    # se crea la carpeta si no existe
    try:
        os.mkdir(os.path.join(PATH, 'pfam_data', directory))
    except:
        pass

    # Nos dirigimos a la carpeta donde se levantará el docker
    os.chdir(PATH)
    # Se ejecuta pfco
    my_cmd = SCRIPT.format(directory, pfam_code)
    os.system(my_cmd)
    # Se cambia al formato FASTA
    align_path = os.path.join(
        PATH, 'pfam_data', directory, pfam_code, 'ALIGN')
    fasta_path = os.path.join(
        PATH, 'pfam_data', directory, pfam_code, 'ALIGN.fasta')
    to_fasta.main(align_path, fasta_path)
    # Se retorna la ruta del archivo
    if os.path.exists(fasta_path):
        return fasta_path
    return {"error": "There's no such ALIGN.fasta file"}


# Execution: python3 pfam_download.py pfam_code
# Example: python3 pfam_download.py PF00023
if __name__ == "__main__":
    main()
