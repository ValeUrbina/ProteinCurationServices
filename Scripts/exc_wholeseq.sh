dir=$1
seed=$2

cd /home/pfam/pfam_data/$dir
wholeseq -align $seed -mu > wholeSEED