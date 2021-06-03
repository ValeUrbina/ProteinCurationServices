import sys
import os

# Delete fragments using belvu

SCRIPT = "/home/valeria/Documentos/Tesis_2/Docker/pfam_curation/Scripts/exc_deletefragments.sh {} {}"


def main(directory, seed):
    my_cmd = SCRIPT.format(directory, seed)
    os.system(my_cmd)
    return 1


if __name__ == "__main__":
    main()
