dir=$1
seed=$2

cd /home/pfam/pfam_data/$dir
create_alignment -fasta $seed -mu > newALIGN