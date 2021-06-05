dir=$1
accnumber=$2

cd /home/pfam/pfam_data/$dir
efetch $accnumber -mu > efetch-$accnumber.txt