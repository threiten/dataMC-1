#!/bin/bash

# qsub all in steps

#k=0.025
#for i in `seq 1 39`;
#do
#    # q=`echo "$k*$i"|bc`
#    qsub jobba.sh "mc"   $q 0 1500000 -q short.q
#    qsub jobba.sh "data" $q 0 -1 -q short.q 
#done


# re-run only a few points
for i in 0.150 0.375;
do
    echo $i
    qsub jobba.sh "mc"   $i 0 1500000 -q short.q
done    




