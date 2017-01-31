#!/bin/bash

obsid=$1

for s in ${obsid}*.stats;do
	gps=${s:11:10}
	freq=${s:24:6}
	gain=`sed -n '/#gain/{n;p;}' $s`

	echo $gps $gain >> ${1}_${freq}.gains
done 
