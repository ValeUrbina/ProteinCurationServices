import sys
import os

# makenonreduntand - aun falta
# this file should be located in the directory specified by pfam_curation_tools docker

SCRIPT = "/home/valeria/Documentos/Tesis_2/Docker/pfam_curation/Scripts/exc_makenonredundant.sh {} {} {}"


def main(directory, cutoff, seed):
    my_cmd = SCRIPT.format(directory, cutoff, seed)
    os.system(my_cmd)
    return 1


if __name__ == "__main__":
    main()
