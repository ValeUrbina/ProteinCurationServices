import sys
import os

# create a HMM over the existing alignment
# this file should be located in the directory specified by pfam_curation_tools docker

SCRIPT = "mkdir /home/valeria/Documentos/Tesis_2/Docker/pfam_curation/pfam_data/{}/"


def main(directory):
    my_cmd = SCRIPT.format(directory)
    os.system(my_cmd)
    return 1


if __name__ == "__main__":
    main()
