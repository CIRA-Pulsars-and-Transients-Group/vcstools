cd ./01
pam -e Fp -Fp coherent.ar
pam --site 1 -m coherent.Fp
pat -s coherent.std coherent.Fp > coherent.tim
tempo2 -f coherent.par -gr plk coherent.tim
