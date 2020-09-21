#!/bin/bash

#Run this script from a directory with an RTS solution and it will create the analyasis plots for you

ATTEMPT_NUM=$(find . -mindepth 1 -maxdepth 1 -type d | wc -l)
ATTEMPT_NUM=$((ATTEMPT_NUM+1))
mkdir "attempt_number_${ATTEMPT_NUM}"

cp BandpassCalibration_node0* ./attempt_number_"${ATTEMPT_NUM}"
cp DI_JonesMatrices_node0*.dat  ./attempt_number_"${ATTEMPT_NUM}"
cp flagged_tiles.txt ./attempt_number_"${ATTEMPT_NUM}"
cp flagged_channels.txt ./attempt_number_"${ATTEMPT_NUM}"

cd ./attempt_number_"${ATTEMPT_NUM}"

for i in $(seq -w 1 24); do
    echo "Plotting channel $i"
    plot_BPcal_128T.py --file "BandpassCalibration_node0${i}.dat" --outname "plot_channel_${i}.png" > "chan_${i}_output.txt"
    plot_BPcal_128T.py --file "BandpassCalibration_node0${i}.dat" -p --outname "plot_phase_${i}.png" > "phase_${i}_output.txt"
done



cd ../
