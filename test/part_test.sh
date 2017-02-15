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

if [ -f 1118168248.mf ]; then
    echo 'metafits file present can run tests'
else
    wget  http://mwa-metadata01.pawsey.org.au/metadata/fits?obs_id=1118168248
    mv fits?obs_id=1118168248 1118168248.mf
    
fi

rm -rf 01
mkdir ./01
cd ./01
cp ../1118168248.mf ./1118168248.mf
size=32768
runner=mpirun

dd if=../ch_100.dat of=xx_0_ch100.dat bs=$size count=10000
dd if=../ch_100.dat of=xx_1_ch100.dat bs=$size skip=10000 count=10000
dd if=../ch_100.dat of=xx_2_ch100.dat bs=$size skip=20000 count=10000
dd if=../ch_100.dat of=xx_3_ch100.dat bs=$size skip=30000 count=10000
dd if=../ch_100.dat of=xx_4_ch100.dat bs=$size skip=40000 count=10000
dd if=../ch_100.dat of=xx_5_ch100.dat bs=$size skip=50000 count=10000
dd if=../ch_100.dat of=xx_6_ch100.dat bs=$size skip=60000 count=10000


dir=`pwd`

rm *.ar
rm *.fits
rm *.vdif
rm *.hdr
rm *.Fp

get_delays -o 1118168248 -m 1118168248.mf -w 10000 -n 128 -p -b 1 -n 128 -z 2015-06-12T18:17:12 -f 128000000 -r 04:37:00 -d -47:15:00 -a .

python ../make_auxfiles.py
echo 100 > channel

$runner make_beam -N 0 -e dat -i -f psrfits_header.txt -d $dir -n 128 -a 128 -r 10000 -o xx -w flags.txt -c phases.txt -t 1 -D $dir/../
dspsr  -c 0.0064 -L 0.1 -A G0024_1118168248_01_0001.fits 
mv G0024_1118168248_01_0001.fits incoherent.fits
mv 2*.ar  incoherent.ar
$runner make_beam -N 0 -e dat -f psrfits_header.txt -d $dir -n 128 -a 128 -r 10000 -o xx -w flags.txt -c phases.txt -t 1 -D $dir/../
dspsr -c 0.0064 -L 0.1 -A G0024_1118168248_01_0001.fits
mv G0024_1118168248_01_0001.fits coherent.fits
mv 2*.ar coherent.ar
$runner make_beam -e dat -v psrfits_header.txt -d $dir -n 128 -a 128 -r 10000 -o xx  -w flags.txt -c phases.txt -t 1 -D $dir/../
dspsr -c 0.0064 -L 0.1 -A -D 0.0 G0024_1118168248_02.hdr

#pam -e F -F coherent.ar
#pam --site 7 -m coherent.F
#pat -s ../coherent.std coherent.F -f tempo2 > coherent.tim
#tempo2 -f ../coherent.par -gr plk coherent.tim
cd ../

rm -rf 02
mkdir ./02
cd ./02
cp ../1118168248.mf ./1118168248.mf
size=32768
runner=mpirun

dd if=../ch_100.dat of=xx_0_ch101.dat bs=$size count=10000
dd if=../ch_100.dat of=xx_1_ch101.dat bs=$size skip=10000 count=10000
dd if=../ch_100.dat of=xx_2_ch101.dat bs=$size skip=20000 count=10000
dd if=../ch_100.dat of=xx_3_ch101.dat bs=$size skip=30000 count=10000
dd if=../ch_100.dat of=xx_4_ch101.dat bs=$size skip=40000 count=10000
dd if=../ch_100.dat of=xx_5_ch101.dat bs=$size skip=50000 count=10000
dd if=../ch_100.dat of=xx_6_ch101.dat bs=$size skip=60000 count=10000


rm *.ar
rm *.fits
rm *.vdif
rm *.hdr
rm *.Fp

dir=`pwd`

get_delays -o 1118168248 -m 1118168248.mf -w 10000 -n 128 -p -b 1 -n 128 -z 2015-06-12T18:17:12 -f 129280000 -r 04:37:00 -d -47:15:00 -a .

python ../make_auxfiles.py
echo 101 > channel

$runner make_beam -N 1 -e dat -i -f psrfits_header.txt -d $dir -n 128 -a 128 -r 10000 -o xx -w flags.txt -c phases.txt -t 1 -D $dir/../
dspsr  -c 0.0064 -L 0.1 -A G0024_1118168248_02_0001.fits
mv G0024_1118168248_02_0001.fits incoherent.fits
mv 2*.ar incoherent.ar
$runner make_beam -N 1 -e dat -f psrfits_header.txt -d $dir -n 128 -a 128 -r 10000 -o xx -w flags.txt -c phases.txt -t 1 -D $dir/../
dspsr -c 0.0064 -L 0.1 -A G0024_1118168248_02_0001.fits
mv G0024_1118168248_02_0001.fits coherent.fits
mv 2*.ar coherent.ar

cd ../

rm splice_coherent.fits
rm splice_incoherent.fits
rm *.ar

splice_psrfits 01/coherent.fits 02/coherent.fits splice_coherent
splice_psrfits 01/incoherent.fits 02/incoherent.fits splice_incoherent

dspsr -c 0.0064 -L 0.1 -A splice_incoherent_0001.fits
mv 2*.ar incoherent_splice.ar

dspsr -c 0.0064 -L 0.1 -A splice_coherent_0001.fits
mv 2*.ar coherent_splice.ar
