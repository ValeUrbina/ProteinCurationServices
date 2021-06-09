dir=$1
seed=$2
outputfile_name=$3

cd /home/pfam/pfam_data/$dir
wholeseq -align $seed -mu > $outputfile_name