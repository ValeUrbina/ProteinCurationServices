dir=$1
pfam=$2
descripcion=$3

cd /home/pfam/pfam_data/$dir
pfci -onlydesc $pfam -i -m $descripcion