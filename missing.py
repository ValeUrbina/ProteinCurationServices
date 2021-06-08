import sys
import os
import codecs

# find the missing sequences
# this file should be located in the directory specified by pfam_curation_tools docker

SCRIPT = "docker run --rm -it -v $(pwd)/Scripts:/Scripts -v $(pwd)/pfam_data:/home/pfam/pfam_data -v $(pwd)/pfam.conf:/home/pfam/pfam.conf -v $(pwd)/seqlib:/data/seqlib -v $(pwd)/Dictionary/dictionary:/home/pfam/Dictionary/dictionary -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix dockerhub.ebi.ac.uk/pfam/pfam-curation bash -c '/Scripts/exc_missing.sh {} {}'"
PATH = '/home/valeria/Documentos/Tesis_2/Docker/pfam_curation/'


def main(directory, pfam_code):
    # Nos dirigimos a la carpeta donde se levantará el docker
    os.chdir(PATH)
    # Se ejecuta missing
    my_cmd = SCRIPT.format(directory, pfam_code)
    result = os.popen(my_cmd).read()
    # Se guarda la ruta del archivo de salida
    missing_path = os.path.join(
        PATH, 'pfam_data', directory, pfam_code, 'missing')
    found_path = os.path.join(
        PATH, 'pfam_data', directory, pfam_code, 'found')

    missing_file = "There's no such missing file"
    found_file = "There's no such found file"

    if os.path.exists(missing_path):
        # Get the data from the file
        with open(missing_path, 'rb') as fp:
            missing_file = fp.read()

        # The data is base64 encoded. Let's decode it.
        missing_file = codecs.decode(missing_file, 'base64')

    if os.path.exists(found_path):
        # Get the data from the file
        with open(found_path, 'rb') as fp:
            found_file = fp.read()

        # The data is base64 encoded. Let's decode it.
        found_file = codecs.decode(found_file, 'base64')

    # Se retorna la ruta del archivo
    return {"result": result, "missing_file": missing_file, "found_file": found_file}


if __name__ == "__main__":
    main()
