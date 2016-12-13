#!/bin/bash

# k=0.1
# 
# 
# for var in  "R9" "S4" "SigmaIeIe" "EtaWidth" "PhiWidth" "CovarianceIphiIphi" "SigmaRR" "CovarianceIetaIphi" "PhoIso03" "ChIso03" "ChIso03worst";
# do
#     echo Jobs for $var
#     for i in `seq 1 9`;
#     do
# 	q=`echo "$k*$i"|bc`
# 	
# 	qsub jobba.sh "mc"   $var $q 0 2000000 3 9 "EB" -q short.q
# 	qsub jobba.sh "mc"   $var $q 0 2000000 3 9 "EE" -q short.q
# 	qsub jobba.sh "data" $var $q 0 -1      3 9 "EB" -q short.q 
# 	qsub jobba.sh "data" $var $q 0 -1      3 9 "EE" -q short.q 
# 
# 	done
# done

qsub jobba.sh "data" ChIso03 0.1 0 -1      3 9 "EB" -q short.q 
qsub jobba.sh "data" ChIso03 0.2 0 -1      3 9 "EB" -q short.q 
qsub jobba.sh "data" ChIso03 0.3 0 -1      3 9 "EB" -q short.q 
qsub jobba.sh "data" ChIso03 0.4 0 -1      3 9 "EB" -q short.q 
qsub jobba.sh "data" ChIso03 0.7 0 -1      3 9 "EB" -q short.q 
qsub jobba.sh "data" ChIso03 0.8 0 -1      3 9 "EB" -q short.q 
qsub jobba.sh "data" ChIso03 0.9 0 -1      3 9 "EB" -q short.q 

qsub jobba.sh "data" ChIso03worst 0.1 0 -1      3 9 "EB" -q short.q 
qsub jobba.sh "data" ChIso03worst 0.2 0 -1      3 9 "EB" -q short.q 
qsub jobba.sh "data" ChIso03worst 0.3 0 -1      3 9 "EB" -q short.q 
qsub jobba.sh "data" ChIso03worst 0.4 0 -1      3 9 "EB" -q short.q 
qsub jobba.sh "data" ChIso03worst 0.6 0 -1      3 9 "EB" -q short.q 
qsub jobba.sh "data" ChIso03worst 0.7 0 -1      3 9 "EB" -q short.q 
qsub jobba.sh "data" ChIso03worst 0.8 0 -1      3 9 "EB" -q short.q 
qsub jobba.sh "data" ChIso03worst 0.9 0 -1      3 9 "EB" -q short.q 

qsub jobba.sh "data" PhoIso03 0.1 0 -1      3 9 "EB" -q short.q 
qsub jobba.sh "data" PhoIso03 0.2 0 -1      3 9 "EB" -q short.q 
qsub jobba.sh "data" PhoIso03 0.3 0 -1      3 9 "EB" -q short.q 
qsub jobba.sh "data" PhoIso03 0.4 0 -1      3 9 "EB" -q short.q 
qsub jobba.sh "data" PhoIso03 0.5 0 -1      3 9 "EB" -q short.q 
qsub jobba.sh "data" PhoIso03 0.6 0 -1      3 9 "EB" -q short.q 
qsub jobba.sh "data" PhoIso03 0.7 0 -1      3 9 "EB" -q short.q 
qsub jobba.sh "data" PhoIso03 0.9 0 -1      3 9 "EB" -q short.q 
