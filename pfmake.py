import sys
import os

# create a HMM over the existing alignment
# this file should be located in the directory specified by pfam_curation_tools docker

SCRIPT = "docker run --rm -it -v $(pwd)/Scripts:/Scripts -v $(pwd)/pfam_data:/home/pfam/pfam_data -v $(pwd)/pfam.conf:/home/pfam/pfam.conf -v $(pwd)/seqlib:/data/seqlib -v $(pwd)/Dictionary/dictionary:/home/pfam/Dictionary/dictionary -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix dockerhub.ebi.ac.uk/pfam/pfam-curation bash -c '/Scripts/exc_pfmake.sh {} {} {}'"


def main(directory, tmayus, tminus):
    my_cmd = SCRIPT.format(directory, tmayus, tminus)
    os.system(my_cmd)
    return 1


if __name__ == "__main__":
    main()
