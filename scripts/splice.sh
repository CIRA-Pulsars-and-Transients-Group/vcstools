#!/bin/bash -l

read -p "Project ID [G0024]: " project
read -p "Observation ID: " obsid
read -p "Pointing [None]: " pointing
read -p "Number of 200-second chunks [2]: " nsets
read -p "Lowest coarse channel number [133]: " lochan
read -p "Highest coarse channel number [156]: " hichan

PROJECT=${project:-G0024}
OBSID=${obsid}
POINTING=${pointing:-None}
NSETS=${nsets:-2}
LO=${lochan:-133}
HI=${hichan:-156}

## If no observation ID is given, abort the process.
if [ -z $OBSID ]; then
    echo "No observation ID provided. Aborting!"
    exit 1
fi

## Since we append to files in the next step, check to
## make sure we don't have any vestigial metafiles left
## over. If we do, delete them.
for f in ${PWD}/splice_*; do

    ## Check if the glob gets expanded to existing files.
    ## If not, f here will be exactly the pattern above
    ## and the exists test will evaluate to false.
    [ -e "$f" ] && rm splice_*

    ## This is all we needed to know, so we can break 
    ## after the first iteration.
    break
done

## For each 200-second chunk, loop over all channels and
## write a temporary metafile containing all of the files
## to be combined.
for j in $(seq -f "%04g" 1 1 $NSETS);do

    sname=splice_${j}

    for i in $(seq -f "%03g" $LO 1 $HI);do 
        if [ $POINTING == None ]; then
            echo ${PROJECT}_${OBSID}_ch${i}_${j}.fits >> $sname
        else
            echo ${PROJECT}_${OBSID}_${POINTING}_ch${i}_${j}.fits >> $sname
        fi
    done

    ## Combine the channels for each individual 200-second 
    ## chunk, then rename it to something more sensible.
    splice_psrfits $(cat $sname) ${OBSID}_tmp
    mv ${OBSID}_tmp_0001.fits ${OBSID}_ch${LO}-${HI}_${j}.fits
done
