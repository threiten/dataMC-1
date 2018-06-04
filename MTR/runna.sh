#!/bin/bash

k=0.1


for var in "probeSigmaIeIe" "probeR9" "probeS4" "probeEtaWidth" "probePhiWidth" "probeCovarianceIetaIphi" "probeCovarianceIphiIphi";
do
    echo Jobs for $var
    for q in 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95;
    do
	#q=`echo "$k*$i"|bc`
	
	qsub -q all.q -l h_vmem=6G jobba.sh "mc"   $var $q 0 -1 3 9 "EB" 
	qsub -q all.q -l h_vmem=6G jobba.sh "mc"   $var $q 0 -1 3 9 "EE"
      	qsub -q all.q -l h_vmem=6G jobba.sh "data" $var $q 0 -1 3 9 "EB"
	qsub -q all.q -l h_vmem=6G jobba.sh "data" $var $q 0 -1 3 9 "EE"
	#qsub -q all.q -l h_vmem=6G jobba.sh "data" $var $q 0 -1 3 9 "EBEE"
	#qsub -q all.q -l h_vmem=6G jobba.sh "mc"   $var $q 0 -1 3 9 "EBEE"
	#qsub -q all.q -l h_vmem=6G jobba_woRF.sh "mc"   $var $q 0 -1 3 9 "EB" 
	#qsub -q all.q -l h_vmem=6G jobba_woRF.sh "mc"   $var $q 0 -1 3 9 "EE"
      	#qsub -q all.q -l h_vmem=6G jobba_woRF.sh "data" $var $q 0 -1 3 9 "EB"
	#qsub -q all.q -l h_vmem=6G jobba_woRF.sh "data" $var $q 0 -1 3 9 "EE"
	#qsub -q all.q -l h_vmem=6G jobba_woRF.sh "data" $var $q 0 -1 3 9 "EBEE"
	#qsub -q all.q -l h_vmem=6G jobba_woRF.sh "mc"   $var $q 0 -1 3 9 "EBEE"
	#qsub -q all.q -l h_vmem=6G jobba_RunC.sh "mc"   $var $q 0 -1 3 9 "EB" 
	#qsub -q all.q -l h_vmem=6G jobba_RunC.sh "mc"   $var $q 0 -1 3 9 "EE"
      	#qsub -q all.q -l h_vmem=6G jobba_RunC.sh "data" $var $q 0 -1 3 9 "EB"
	#qsub -q all.q -l h_vmem=6G jobba_RunC.sh "data" $var $q 0 -1 3 9 "EE"
	#qsub -q all.q -l h_vmem=6G jobba_RunC.sh "data" $var $q 0 -1 3 9 "EBEE"
	#qsub -q all.q -l h_vmem=6G jobba_RunC.sh "mc"   $var $q 0 -1 3 9 "EBEE"
	done
done
