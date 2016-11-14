#!/bin/bash

k=0.025

#for i in `seq 1 39`;
for i in 0.0150 0.275 0.350 0.375 0.5000 0.625;
do
    q=`echo "$k*$i"|bc`
    # echo $q
    qsub jobba.sh "mc"   $q 0 1500000 -q short.q
    #qsub jobba.sh "data" $q 0 -1 -q short.q 
done    




