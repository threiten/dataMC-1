#!/bin/bash
#
echo SWITCH OFF DISPLAY
export DISPLAY=
echo DISPLAY = $DISPLAY

source /mnt/t3nfs01/data01/shome/giulioisac/setupython.sh

echo JOBBA: Training on $1 for Y = $2 , quantile = $3,  startEvt = $4, stopEvt = $5, maxDepth = $6, minLeaf = $7

python /mnt/t3nfs01/data01/shome/giulioisac/dataMC/MTR/train_quantileRegression_BatchAll.py $1 $2 $3 $4 $5 $6 $7

echo JOBBA DONE: Training on $1 for Y = $2 , quantile = $3,  startEvt = $4, stopEvt = $5, maxDepth = $6, minLeaf = $7
