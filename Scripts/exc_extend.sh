dir=$1
seed=$2
n=$3
c=$4

cd /home/pfam/pfam_data/$dir
extend -align $seed -n $n -c $c -mu > extendedSEED