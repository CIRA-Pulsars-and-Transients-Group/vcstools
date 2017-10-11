#!/bin/bash

obsid=$1

for s in ${obsid}*.stats;do
    gps=${s:11:10}
    freq=${s%%MHz*}
    freq=${freq##*_}
    gain=`sed -n '/#gain/{n;p;}' $s`

    echo $gps $gain >> ${1}_${freq}MHz.gains
done 
