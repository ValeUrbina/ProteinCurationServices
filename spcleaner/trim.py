import sys

# http://trimal.cgenomics.org/downloads
# Delete gaps on alignments
# -A-D             AD
# -B-E >>becomes>> BE
# -C-F             CF


def trimmer(source_alignment_path, target_alignment_path):
    myCmd = 'trimal -in ' + source_alignment_path + \
        ' -out ' + target_alignment_path + ' -noallgaps'


def main():
    source_alignment_path = sys.argv[1]
    target_alignment_path = sys.argv[2]
    trimmer(source_alignment_path, target_alignment_path)
    return


# Execution: trim.py input_file_path output_file_path
# Example: trim.py  PF000/alignment PF000/trimmed_alignment
if __name__ == "__main__":
    main()
