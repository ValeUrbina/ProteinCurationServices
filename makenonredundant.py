import sys
import os

# makenonreduntand
# this file should be located in the directory specified by pfam_curation_tools docker
# USAR LA FUNCION DE skipredundant

# SCRIPT = "/home/valeria/Documentos/Tesis_2/Docker/pfam_curation/Scripts/exc_makenonredundant.sh {} {} {}"
SCRIPT = "/home/ubuntu/tesis/pfam_curation/Scripts/exc_makenonredundant.sh {} {} {} {}"


def main(directory, cutoff, file_name, outputfile_name):
    my_cmd = SCRIPT.format(directory, cutoff, file_name, outputfile_name)
    os.system(my_cmd)
    return 1


if __name__ == "__main__":
    main()
