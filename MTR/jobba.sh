#!/bin/bash
#
echo SWITCH OFF DISPLAY
export DISPLAY=
echo DISPLAY = $DISPLAY

source /mnt/t3nfs01/data01/shome/mdonega/bootJupyter.sh

echo JOBBA: Training for :  $1 , quantile = $2,  startEvt = $3, stopEvt = $4
python /mnt/t3nfs01/data01/shome/mdonega/dataMC/MTR/train_quantileRegression_Batch.py $1 $2 $3 $4
echo JOBBA DONE: Training for :  $1 , quantile = $2,  startEvt = $3, stopEvt = $4