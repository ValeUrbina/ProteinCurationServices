dir=$1
seed=$2
method=$3
outputfile_name=$4

cd /home/pfam/pfam_data/$dir
create_alignment -fasta $seed -$method > $outputfile_name