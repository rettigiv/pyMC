from random import random
import sys

N = int(sys.argv[1])

outfile = open(('in%d.txt' % N), 'w')
outfile2 = open(('startVMD%d.gro' % N), 'w')

l = 1.0

outfile2.write('Starting position in GROMACS format\n')
outfile2.write('%d\n' % N)

for i in range(0, N):
    x = random()*l
    y = random()*l
    z = random()*l

    m = 2.0
    q = 1.0
    sig = 0
    eps = 0
    
    outfile.write('AT \t %3.3f \t %3.3f \t %3.3f \t %3.3f \t %3.3f \t %3.3f \t %3.3f\n' % (x, y, z, m, q, sig, eps))
    outfile2.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n' % (i+1, 'ATM', 'AT', i+1, x, y, z, 0, 0, 0))

outfile2.write('%8.4f %8.4f %8.4f\n' % (l, l, l))
outfile.close()
outfile2.close()
