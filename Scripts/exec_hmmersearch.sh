path=$1
target_path_hmm_profile=$2
source_path_alignment=$3
pfamseq_path=$4

cd $path
hmmbuild $target_path_hmm_profile $source_path_alignment
hmmsearch --tblout tblout.txt --domtblout domtblout.txt --pfamtblout pfamtblout.txt hmmALIGN $pfamseq_path