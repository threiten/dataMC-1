#!/bin/bash

# qsub all in steps

#  "R9", "S4", "SigmaIeIe", "EtaWidth", "PhiWidth", "CovarianceIphiIphi", "SigmaRR"]

k=0.025

for i in `seq 1 1`;
do
    q=`echo "$k*$i"|bc`

    qsub jobba.sh "mc"   R9 $q 0 1500000 -q short.q
    qsub jobba.sh "data" R9 $q 0 -1 -q short.q 

    qsub jobba.sh "mc"   SigmaIeIe $q 0 1500000 -q short.q
    qsub jobba.sh "data" SigmaIeIe $q 0 -1 -q short.q 
done


# re-run only a few points
#for i in 0.98;
#do
#    echo $i
#
#done    
#



