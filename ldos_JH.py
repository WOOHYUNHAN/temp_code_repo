#!/usr/bin/env python
#############################################################################
# This program was written for calculating LDOS of NW                       #
# Program was designed by J.-S. Park.                                       #
# Latest update is done at 2011-07-25.                                      #
# This program is open source, but would not be distributed by the internet.#
# INPUT file? PARCHG files in the specific directory and EIGENVAL           #
#############################################################################

import os.path
import os
from math import *
from math import sqrt

def value(x,mu,sigma):
        if (abs(x-mu)/2/sigma > 4) :
                return 0
        else :
                return (1/sigma/sqrt(2* pi))* exp(-pow(x-mu,2) / 2/ pow(sigma,2))

def getMedian(a):
  a_len = len(a)
  if (a_len == 0): return None
  a_center = a_len / 2

  if (a_len % 2 == 1):
    return a[a_center]
  else:
    return (a[a_center - 1] + a[a_center]) / 2.0

def get(name,cut,energy,total,VASP):
        state = 0
        value = []
        sc = []
        posx = []
        posy = []
        posz = []

        number = 0
        grid = []

        f = open(name)
        buffer = f.readlines()
        length = len(buffer)

        for i in xrange(2,5) :
                sc.append(buffer[i].split())

        scx = float(sc[0][0])
        scy = float(sc[1][1])
        scz = float(sc[2][2])
        #print scx

        atom = buffer[5+VASP].split()
        for j in xrange(len(atom)) :
                number = number + int(atom[j])
        #print number

        for i in xrange(7+VASP,7+VASP+number) :
                temp = buffer[i].split()
                posx.append(float(temp[0]))
                posy.append(float(temp[1]))
                posz.append(float(temp[2]))

        posx.sort()
        posy.sort()
        posz.sort()

        centerx = scx * getMedian(posx)
        centery = scy * getMedian(posy)
        centerz = scz * getMedian(posz)

        grid = buffer[8+VASP+number].split()
        gridx = int(grid[0])
        gridy = int(grid[1])
        gridz = int(grid[2])

        nz = [] # in this case, nxy is nz. xy averaged. 

        for i in xrange(9+VASP+number,length) :
                temp = buffer[i].split()
                for j in xrange(len(temp)) :
                	value.append(float(temp[j]))

        for k in xrange(gridz) :
                nz.append(0.0)
                for ij in xrange(gridx*gridy):
                        nz[k] = nz[k]+ value[k*gridx*gridy + ij]
                #nxy[ij] = nxy[ij] / gridz / (scx * scy * scz)
                nz[k] = nz[k]/gridx/gridy / (scx * scy * scz)

        return nz

#print 'STRAT PROGRAM'
bandinfo = []
nameinfo = []

#print 'OPEN EIGENVAL'
f = open('EIGENVAL')
buffer = f.readlines()
length = len(buffer)

#print 'READ HEAD OF EIGENVAL'
temp = buffer[5].split()
nBAND= int(temp[2])
KPOINTS = int(temp[1])

##############################################################
# Input parameters                                           #
##############################################################
EMIN = 4.46
EMAX = 8.46 
total = 420
gridE = 500
sigma = 0.06
##############################################################

#print 'NOW READING ENERGY PART OF EIGENVAL FILE'
for j in xrange(0,KPOINTS):
        for i in xrange(0,nBAND):
                #print buffer[i + 8 + (nBAND+2)*j]
                temp = buffer[i + 8 + (nBAND+2)*j].split()
                #print temp[1 + 8 + (nBAND+2)*j]
                temp = float(temp[1])
                if (EMIN < temp):
                        if (EMAX > temp):
                                bandinfo.append([i+1,j+1,temp])

#print bandinfo

#print 'MAKING BANDINFO'
for i in xrange(len(bandinfo)):
        a = "%.4d" % bandinfo[i][0]
        b = "%.4d" % bandinfo[i][1]
        temp = 'PARCHG.' + a + '.' + b
        #print temp
        nameinfo.append(temp)

#print 'GETTING CHARGE DENSITY'
#Increase it after program is well made
N = []
for i in xrange(len(nameinfo)):
        N.append(get(nameinfo[i],0.1,bandinfo[i][2],total,1))

#print 'MAKING ENERGY GRID'
E = []
for i in xrange(gridE):
        E.append(EMIN+i*(EMAX-EMIN)/(gridE-1))

#print 'MAKING EMPTY MATRIX AND FILL IT BY DENSITY OF STATES'
result = []
for i in xrange(len(E)):
        result.append([0]*total)

#print 'CALCULATE MATRIX COMPONENTS'
for i in xrange(len(E)):
        for j in xrange(total):
                for k in xrange(len(bandinfo)):
                        result[i][j] = result[i][j] + value(E[i],bandinfo[k][2],sigma) * N[k][j]

output = open('LDOS.txt','w')
#print 'WRITING'
for i in xrange(len(E)):
        for j in xrange(total):
                tt = str(i) + " " + str(j) + " " + str(result[i][j])
                output.write(tt+'\n')

