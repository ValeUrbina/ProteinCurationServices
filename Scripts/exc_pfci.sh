dir=$1
options=$2
pfam=$3
descripcion=$4

cd /home/pfam/pfam_data/$dir
pfci $options $pfam -i -m $descripcion