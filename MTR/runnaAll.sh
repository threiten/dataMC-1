#!/bin/bash

k=0.1


for var in  "R9" # "S4" "SigmaIeIe" "EtaWidth" "PhiWidth" "CovarianceIphiIphi" "SigmaRR" "CovarianceIetaIphi" "PhoIso03" "ChIso03" "ChIso03worst";
do
    echo Jobs for $var
    for i in `seq 1 9`;
    do
	q=`echo "$k*$i"|bc`
	qsub jobbaAll.sh "mc"   $var $q 0 2000000 3 9  -q short.q
      	qsub jobbaAll.sh "data" $var $q 0 -1      3 9  -q short.q 
	done
done
