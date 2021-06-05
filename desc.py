import sys
import os

# return desc file
# this file should be located in the directory specified by pfam_curation_tools docker

PATH = '/home/valeria/Documentos/Tesis_2/Docker/pfam_curation/'


def main(directory, pfam_code):
    # Nos dirigimos a la carpeta donde se levantará el docker
    desc_path = os.path.join(
        PATH, 'pfam_data', directory, pfam_code, 'DESC')
    # Se retorna la ruta del archivo
    return desc_path


if __name__ == "__main__":
    main()
