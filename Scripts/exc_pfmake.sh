dir=$1
Tmayus=$2
tminus=$3
e=$4

cd /home/pfam/pfam_data/$dir
pfmake -T $Tmayus -t $tminus -e $e