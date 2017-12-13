import sys
import MDAnalysis
from random import random
from math import exp

class Molecule:
    def __init__(self, i, x, y, z, m, q, sigma, epsilon):
        self.i = i
        
        self.x = x
        self.y = y
        self.z = z

        self.vx = 0.0
        self.vy = 0.0
        self.vz = 0.0
        
        self.q = q
        self.sigma = sigma
        self.epsilon = epsilon



numMols = 0
maxChangeDist = 0.5
boxX = 1.0
boxY = 1.0
boxZ = 1.0

universe = MDAnalysis.Universe(('startVMD' + sys.argv[1][2:sys.argv[1].index('.')]+ '.gro'))
writer = MDAnalysis.Writer(('out' + sys.argv[1][2:sys.argv[1].index('.')] + '.trr'), 1)
outColor = open('color.txt','w')
E = 0.0
avgEn = 0.0
molecules = []
stepNum = 1

R = 1.9872156e-3
B = 0.0
T = 0.0
N = 0
nSteps = 0
stepNum = 1
densCenter = 0.0
densities = [0]*1000
Z = 0.0

numLower = 0
numHigherAcc = 0
numHigherRej = 0

def wrapDist(x1, x2):
    return min(abs(x1 - x2), abs(1-max(x1,x2) + min(x1,x2)))
    
def dist((x1,y1,z1), (x2,y2,z2)):
    return ((x2 - x1)**2.0 + (y2 - y1)**2.0 + (z2 - z1)**2.0)**0.5

        
def coloumb(mol1, mol2):
    colPer = 0
    for xDisp in range(-1,1):
        for yDisp in range(-1,1):
            for zDisp in range(-1,1):
                if(xDisp == 0 and yDisp == 0 and zDisp == 0 and mol1.i <= mol2.i):
                    continue
                d = dist((mol1.x, mol1.y, mol1.z), (mol2.x + xDisp, mol2.y + yDisp, mol2.z + zDisp))
                if(d > 0.5):
                    continue
                colPer += (138.935485)*(mol1.q * mol2.q) / d

    return colPer

def vdw(mol1, mol2):
    return 0.0


def PotentialEnergy(mol):
    pEnergy = 0.0
    for i in range(0, N):
        pEnergy += coloumb(mol, molecules[i])
        pEnergy += vdw(mol, molecules[i])
    return pEnergy


def KineticEnergy(mol):
    return 0.0
    return m*(vx**2.0 + vy**2.0 + vz**2.0)/2.0


def Hamiltonian():
    ener = 0.0
    for i in range(0,N):
        ener += PotentialEnergy(molecules[i])
        ener += KineticEnergy(molecules[i])

    return ener


def writeCoords():
    for i in range(0,N):
        universe.atoms[i].position = [molecules[i].x * 10.0, molecules[i].y*10.0, molecules[i].z*10.0]
    writer.write(universe.atoms)

def writeColor(i, c):
    outColor.write('%d %d\n' %(i, c))

        
def shiftMol(mol):
    mol.x = (mol.x + (maxChangeDist * 2 * random() - maxChangeDist)) % boxX
    mol.y = (mol.y + (maxChangeDist * 2 * random() - maxChangeDist)) % boxY
    mol.z = (mol.z + (maxChangeDist * 2 * random() - maxChangeDist)) % boxZ
    

def stepZ():
    global E
    global Z
    
    # move a random molecule
    for molInd in range(0,N):
        sMol = molecules[molInd]
        shiftMol(sMol)

    En = Hamiltonian()
    binNum = int((En - densCenter*0.4)/(densCenter*0.8) * 1000)

    if(binNum > 999):
        binNum = 999
    if(binNum < 0):
        binNum = 0

    densities[binNum] += 1
    Z += exp(-B * En)

def step():
    global E
    global avgEn
    global numLower
    global numHigherAcc
    global numHigherRej
    global Z
    
    # move a random molecule
    sMol = molecules[int(random()*N)]

    # color the particle to chosen
    writeColor(sMol.i, 50)
    writeCoords()
    
    smX = sMol.x
    smY = sMol.y
    smZ = sMol.z
    shiftMol(sMol)

    # now show shifted particle
    writeCoords()

    En = Hamiltonian()
    if(En < E):
        # downhill, keep it
        E = En
        numLower += 1
        writeColor(sMol.i, 100)
    else:
        # uphill, take it with boltzmann prob
        binNum = int((En - densCenter*0.4)/(densCenter*0.8) * 1000)
        if(binNum > 999):
            binNum = 999
        if(binNum < 0):
            binNum = 0

        prob = exp(-B * (En - E)) *densities[binNum] / Z
        if(random() < prob):
            # keep the change
            E = En
            writeColor(sMol.i, 100)
            numHigherAcc += 1
        else:
            # revert change
            shiftMol.x = smX
            shiftMol.y = smY
            shiftMol.z = smZ
            writeColor(sMol.i, 25)
            numHigherRej += 1

    #display colored mol then uncolor
    writeCoords()
    writeColor(sMol.i, 0)

    
    #store step energy for avg at end
    avgEn = avgEn + E


def handleInput():
    global N
    infile = open(sys.argv[1], 'r')

    molNum = 0
    for line in infile:
        args = line.split()
        x = float(args[1])
        y = float(args[2])
        z = float(args[3])

        m = float(args[4])*1.66e-27
        q = float(args[5])
        sig = float(args[6])
        eps = float(args[7])

        newMol = Molecule(molNum, x, y, z, m, q, sig, eps)
        molecules.append(newMol)
        molNum = molNum + 1

    N = molNum
        
    
def main():
    global E
    global writer
    global stepNum
    global densCenter

    # get run params
    nSteps = int(sys.argv[3])
    T = float(sys.argv[4])
    B = 1.0 / (R * T)

    # load up the input coords
    handleInput()

    # define starting energy
    E = Hamiltonian()
    print (E)

    densCenter = E
    for sNum in range(1,nSteps):
        stepZ()
    #writer.close()
    #writer = MDAnalysis.Writer('out.trr', N)
    #writer.write(universe.atoms)
    writeCoords()
    

    # loop through MC steps
    for stepNum in range(1, nSteps):
        step()
        writeCoords()

    #writer.close()
    print('kJ = %f\n' % (avgEn / nSteps))
    print('kJ norm = %f\n' % (avgEn / nSteps / N))
    print('Lower = %d\n' % (numLower))
    print('HigherAcc = %d\n' % (numHigherAcc))
    print('HigherRej = %d\n' % (numHigherRej))



main()






