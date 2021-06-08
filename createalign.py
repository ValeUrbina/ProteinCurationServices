import sys
import os
import to_fasta

# Align an alignment using createaling function
# this file should be located in the directory specified by pfam_curation_tools docker

SCRIPT = "docker run --rm -it -v $(pwd)/Scripts:/Scripts -v $(pwd)/pfam_data:/home/pfam/pfam_data -v $(pwd)/pfam.conf:/home/pfam/pfam.conf -v $(pwd)/seqlib:/data/seqlib -v $(pwd)/Dictionary/dictionary:/home/pfam/Dictionary/dictionary -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix dockerhub.ebi.ac.uk/pfam/pfam-curation bash -c '/Scripts/exc_createalign.sh {} {} {}'"


def main(PATH, directory, pfam_code, seed, method):
    # Nos dirigimos a la carpeta donde se levantará el docker
    os.chdir(PATH)
    # Se ejecuta wholeseq
    my_cmd = SCRIPT.format(directory + '/' + pfam_code, seed, method)
    os.system(my_cmd)
    # Se cambia al formato FASTA
    align_path = os.path.join(
        PATH, 'pfam_data', directory, pfam_code, 'newALIGN')
    fasta_path = os.path.join(
        PATH, 'pfam_data', directory, pfam_code, 'newALIGN.fasta')
    to_fasta.main(align_path, fasta_path)
    # Se retorna la ruta del archivo
    return fasta_path


if __name__ == "__main__":
    main()
