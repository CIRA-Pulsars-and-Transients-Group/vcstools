ninp = 256
nchan = 128
nstep = 1000

f = open('flags.txt','w')

for input in range(ninp):
    f.write('1.0\n')
f.close()

f = open ('phases.txt','w')
for step in range (nstep):
    for input in range(ninp):
        for chan in range(nchan):
            f.write('0.0')
f.close()

f = open ('channel','w')
f.write('100')

