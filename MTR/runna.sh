#!/bin/bash

k=0.1

for var in  "R9" "S4" "SigmaIeIe" "EtaWidth" "PhiWidth" "CovarianceIphiIphi" "SigmaRR";
do
    echo Jobs for $var
    # for i in `seq 1 9`;
    for i in `seq 1 1`;
    do
	q=`echo "$k*$i"|bc`
	
	qsub jobba.sh "mc" $var $q 0 1500000 3 9 -q short.q
	qsub jobba.sh "data" $var $q 0 -1 3 9 -q short.q 

#	qsub jobba.sh "mc" $var $q 0 1500000  10 9 -q short.q
#	qsub jobba.sh "data" $var $q 0 -1  10 9 -q short.q 
#
#	qsub jobba.sh "mc" $var $q 0 1500000  3 50 -q short.q
#	qsub jobba.sh "data" $var $q 0 -1  3 50 -q short.q 
#
#	qsub jobba.sh "mc" $var $q 0 1500000  10 50 -q short.q
#	qsub jobba.sh "data" $var $q 0 -1  10 50 -q short.q 
	
	done
done
