#!/bin/bash
group=`lfs quota -u $USER /group | tail -1 | awk '{print $2/1e9}'`
scratch2=`lfs quota -u $USER /scratch2 | tail -1 | awk '{print $2/1e9}'`
echo "You're using ${group} TB on /group"
echo "You're using ${scratch2} TB on /scratch2"
