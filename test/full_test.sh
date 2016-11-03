aprun test_mkfile -c 128 -i 256 -t 100000 -o ch_100.dat
rm -rf ./01
mkdir ./01
cd ./01
mv ../ch_100.dat ./
cp ../psrfits_header.txt ./
python ../make_auxfiles.py
dir=`pwd`

aprun make_beam -e dat -i -f psrfits_header.txt -d $dir -n 128 -a 128 -r 10000 -o xx -w flags.txt -c phases.txt -t 1 -D $dir/../
dspsr -c 0.0064 G0024_1118168248_01_0001.fits 
mv 2015-06-12-18\:17\:11.ar incoherent.ar
rm G0024_1118168248_01_0001.fits
aprun make_beam -e dat -f psrfits_header.txt -d $dir -n 128 -a 128 -r 10000 -o xx -w flags.txt -c phases.txt -t 1 -D $dir/../
dspsr -c 0.0064 G0024_1118168248_01_0001.fits
mv 2015-06-12-18\:17\:11.ar coherent.ar
rm G0024_1118168248_01_0001.fits
aprun make_beam -e dat -v psrfits_header.txt -d $dir -n 128 -a 128 -r 10000 -o xx -w flags.txt -c phases.txt -t 1 -D $dir/../
dspsr -c 0.0064 G0024_1118168248_01.hdr
