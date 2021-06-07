module load python/2.7.11
module load tmalign
module load dssp
module load mustang
i=$SGE_TASK_ID
mkdir /home/layhir/benchmark/pid_$i
python  /projects/Reupred/runPredForTask.py targetListPred$i /home/layhir/benchmark/pid_$i/

