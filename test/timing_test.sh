cd ./01

pam -e Fp -Fp coherent.ar
pam --site 7 -m coherent.Fp
pat -s ../coherent.std coherent.Fp -f \"tempo2\" > coherent.tim
tempo2 -f ../coherent.par -gr plk coherent.tim
cd ..
