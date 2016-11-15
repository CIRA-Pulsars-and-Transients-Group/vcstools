#!/bin/bash

if [ -f ch_100.dat ]; then
    echo 'ch_100.dat present - can run tests'
else
    echo 'ch_100.dat not present ... need to pull from repo - this may take a couple of minutes - its a ~3GB file'
    if ! [ -x "$(command -v git)" ]; then
        echo 'git is not installed! buildit' >&2
        exit
    fi
    rm -rf remote-data
    mkdir remote-data
    
    cd remote-data
    git init
    git remote add origin -t feature/vcs https://ord006@bitbucket.csiro.au/scm/~ord006/data.git
    git pull
    cd -
    ln -s remote-data/vcs/ch_100.dat 
    
    if [ -f ch_100.dat ]; then
        echo 'ch_100.dat present - can run tests'
    else
        echo 'failed retrieval ....'
        exit
    fi
fi

if ! [ -x "$(command -v make_beam)" ]; then
  echo 'make_beam is not installed! buildit' >&2
  exit
fi
if ! [ -x "$(command -v dspsr)" ]; then
  echo 'dspsr is not installed! buildit' >&2
  exit
fi


rm -rf 01
mkdir ./01
cd ./01
size=32768
runner=mpirun

dd if=../ch_100.dat of=xx_0_ch100.dat bs=$size count=10000
dd if=../ch_100.dat of=xx_1_ch100.dat bs=$size skip=10000 count=10000
dd if=../ch_100.dat of=xx_2_ch100.dat bs=$size skip=20000 count=10000
dd if=../ch_100.dat of=xx_3_ch100.dat bs=$size skip=30000 count=10000
dd if=../ch_100.dat of=xx_4_ch100.dat bs=$size skip=40000 count=10000
dd if=../ch_100.dat of=xx_5_ch100.dat bs=$size skip=50000 count=10000
dd if=../ch_100.dat of=xx_6_ch100.dat bs=$size skip=60000 count=10000

cp ../psrfits_header.txt ./
python ../make_auxfiles.py
dir=`pwd`

rm *.ar
rm *.fits
rm *.vdif
rm *.hdr
rm *.Fp

$runner make_beam -e dat -i -f psrfits_header.txt -d $dir -n 128 -a 128 -r 10000 -o xx -w flags.txt -c phases.txt -t 1 -D $dir/../
dspsr  -c 0.0064 -L 0.1 -A G0024_1118168248_01_0001.fits 
mv G0024_1118168248_01_0001.fits incoherent.fits
mv 2015-06-12-18\:17\:11.ar incoherent.ar
$runner make_beam -e dat -f psrfits_header.txt -d $dir -n 128 -a 128 -r 10000 -o xx -w flags.txt -c phases.txt -t 1 -D $dir/../
dspsr -c 0.0064 -L 0.1 -A G0024_1118168248_01_0001.fits
mv G0024_1118168248_01_0001.fits coherent.fits
mv 2015-06-12-18\:17\:11.ar coherent.ar
$runner make_beam -e dat -v psrfits_header.txt -d $dir -n 128 -a 128 -r 10000 -o xx  -w flags.txt -c phases.txt -t 1 -D $dir/../
dspsr -c 0.0064 -L 0.1 -A -D 0.0 G0024_1118168248_01.hdr
pam -e F -F coherent.ar
pam --site 7 -m coherent.F
pat -s ../coherent.std coherent.F -f tempo2 > coherent.tim
tempo2 -f ../coherent.par -gr plk coherent.tim
cd ../
