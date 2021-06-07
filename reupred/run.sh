module load python/2.7.13

i=$1
d=$2
e=$5
q=$4
c=$3
v=$6
#python /home/layla/PycharmProjects/Server/reupred.py -i $i -d $d -c $c -q $q -e $e -s /home/layla/PycharmProjects/Server/ -v $v
python /home/layla/PycharmProjects/RepeatsDB-lite0/reupred.py -i $i -d $d -c $c -q $q -e $e -s /home/layla/PycharmProjects/RepeatsDB-lite0/ -v $v
#python /home/lisanna/Projects/repeatsdb-lite-app/reupred.py -i $i -d $d -c $c -q $q -e $e -s /home/lisanna/Projects/repeatsdb-lite-app/ -v $v
#/home/lisanna/Projects/repeatsdb-lite-app/run.sh 3vsf /home/lisanna/Projects/repeatsdb-lite-ws/output/repeatsdb-lite-output/c3f5df29-6c48-4e18-bd34-65707ebe8289/ - all true false
