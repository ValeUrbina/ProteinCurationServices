import sys
import os
import to_fasta

# Save a new a non existing family
# this file should be located in the directory specified by pfam_curation_tools docker

SCRIPT = "docker run --rm -it -v $(pwd)/Scripts:/Scripts -v $(pwd)/pfam_data:/home/pfam/pfam_data -v $(pwd)/pfam.conf:/home/pfam/pfam.conf -v $(pwd)/seqlib:/data/seqlib -v $(pwd)/Dictionary/dictionary:/home/pfam/Dictionary/dictionary -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix dockerhub.ebi.ac.uk/pfam/pfam-curation bash -c '/Scripts/exc_pfnew.sh {} {}'"
PATH = '/home/valeria/Documentos/Tesis_2/Docker/pfam_curation/'


def main(directory, accnumdir):
    # Nos dirigimos a la carpeta donde se levantará el docker
    os.chdir(PATH)
    my_cmd = SCRIPT.format(directory, accnumdir)
    # Se ejecuta pfnew
    try:
        os.system(my_cmd)
        # Se cambia al formato FASTA
        pfnew_path = os.path.join(
            PATH, 'pfam_data', directory, 'pfnew.txt')
        # Se retorna la ruta del archivo
        return pfnew_path
    except:
        return "pfnew error"


if __name__ == "__main__":
    main()
