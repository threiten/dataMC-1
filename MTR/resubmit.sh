#!/bin/bash

rm resubmit.txt*

for file in `ls ~/jobba.sh.o*`:
do
    if grep 'DONE' $file;
    then
	echo DONE
    else
	s=`grep 'JOBBA' $file`
	echo LINE = $s
	a=( $s )
	echo $s > bla.tmp
	y=`awk '{print $6}' bla.tmp | awk 'BEGIN{FS="="}{print $2}'`
	rm bla.tmp
	words=\"${a[3]}\"" "$y" "${a[9]}" "${a[12]}" "${a[15]}" "${a[18]}" "${a[21]} 
	echo $words > bla.tmp
	echo $words
	args=`sed 's/\,//g' bla.tmp`
	echo NO COMMAS $args
	echo qsub jobba.sh $args -q short.q >> resubmit.txt
	rm bla.tmp
    fi

done

	

# JOBBA: Training on mc for Y=R9 , quantile = .2, startEvt = 0, stopEvt = 1500000, maxDepth = 3, minLeaf = 9
