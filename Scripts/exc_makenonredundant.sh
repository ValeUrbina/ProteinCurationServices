dir=$1
cutoff=$2
seed=$3

cd /home/valeria/Documentos/Tesis_2/Docker/pfam_curation/pfam_data/$dir
belvu -n $cutoff $seed