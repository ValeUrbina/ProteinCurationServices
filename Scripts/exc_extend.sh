dir=$1
seed=$2
n=$3
c=$4
method=$5

cd /home/pfam/pfam_data/$dir
extend -align $seed -n $n -c $c -$method > extendALIGN