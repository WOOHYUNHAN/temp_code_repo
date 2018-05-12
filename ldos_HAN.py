import os.path
import time
from math import *
from math import sqrt
from numpy.linalg import norm
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
#from numeric import *

def make_difference_CHGCAR(ref_file, target_file, output_name):
    f = open(ref_file)
    fbuffer = f.readlines()
    f.close()

    g = open(target_file)
    gbuffer = g.readlines()
    g.close()

    univf = float(fbuffer[1].split()[0])
    latt_af = np.array([float(fbuffer[2].split()[i]) for i in range(3)]) * univf
    latt_bf = np.array([float(fbuffer[3].split()[i]) for i in range(3)]) * univf
    latt_cf = np.array([float(fbuffer[4].split()[i]) for i in range(3)]) * univf

    univg = float(gbuffer[1].split()[0])
    latt_ag = np.array([float(gbuffer[2].split()[i]) for i in range(3)]) * univg
    latt_bg = np.array([float(gbuffer[3].split()[i]) for i in range(3)]) * univg
    latt_cg = np.array([float(gbuffer[4].split()[i]) for i in range(3)]) * univg

    #if not latt_af == latt_ag and latt_af == latt_ag and latt_af == latt_ag:
    #    print 'Two charge files have different cell size'
    #    return 0




def average_CHGCAR(input_file='CHGCAR', average_axis='z'):
    f = open(input_file)
    fbuffer = f.readlines()
    f.close()
    univ = float(fbuffer[1].split()[0])
    latt_a = np.array([float(fbuffer[2].split()[i]) for i in range(3)]) * univ
    latt_b = np.array([float(fbuffer[3].split()[i]) for i in range(3)]) * univ
    latt_c = np.array([float(fbuffer[4].split()[i]) for i in range(3)]) * univ
    a = np.linalg.norm(latt_a)
    b = np.linalg.norm(latt_b)
    c = np.linalg.norm(latt_c)
    print a, b, c

    num_data=5

    for i in range(10):
        temp = fbuffer[i].split()
        if temp[0] == 'Direct':
            start = i
            break

    temp_array = fbuffer[start-1].split()
    #print len(temp_array)
    atom = np.array([int(temp_array[k]) for k in range(len(temp_array))])
    total = np.sum(atom)

    grid_line = start + total + 1 + 1
    #print grid_line
    grid_data = [int(fbuffer[grid_line].split()[i]) for i in range(3)]
    NX, NY, NZ = grid_data[0], grid_data[1], grid_data[2]

    print NX, NY, NZ

    if NX*NY*NZ % num_data == 0:
        expected_line = ((NX*NY*NZ)/num_data)
        last_line = num_data
    else:
        expected_line = ((NX*NY*NZ)/num_data) + 1
        last_line = np.mod(NX*NY*NZ, num_data)
    value = []



    #print grid_line+1+expected_line

   
    for i in range(grid_line+1, grid_line+1+expected_line-1):
        for j in range(num_data):
            value.append(float(fbuffer[i].split()[j]))

    
    for j in range(last_line):
        value.append(float(fbuffer[grid_line+1+expected_line-1].split()[j]))
    #print value


    if average_axis == 'x':
        print "You choose x axis as standard axis, now calculation for x-y averaged values will operate"
        #X, a, b = NX, NY, NZ
    elif average_axis == 'y':
        print "You choose y axis as standard axis, now calculation for x-z averaged values will operate"
        #X, a, b = NY, NX, NZ
    elif average_axis == 'z':
        print "You choose z axis as standard axis, now calculation for x-y averaged values will operate"
        #X, a, b = NZ, NX, NY
    else:
        print "input axis is something wrong"

    averaged_value=[]
    for i in range(NZ):
        temp = 0.0
        for j in range(NX):
            for k in range(NY):
                index = (NX*NY) * i + NX*k + j
                #index = (NY*NZ)*(j) + NZ*(k) + (i)
                temp += value[index]
                #print temp
        temp = temp / (NX*NY)
        #print temp
        averaged_value.append(temp)

    print len(averaged_value)

    dz = c / NZ
    g = open('averaged_file.out','w')
    for i in range(NZ):
        line = str(dz*i)+'\t'+str(averaged_value[i])+'\n'
        g.write(line)
    g.close()






def plot_ldos(input_file='LDOS.txt'):
    EMIN=4.2376
    EMAX=8.2376

    gridE=500
    base=6.2376 # Fermi energy
    #total = 240
    total = 400
    #sc = 19.237
    sc = 21.2062244259443133


    fig = plt.figure()
    f = open(input_file)
    fbuffer = f.readlines()
    f.close()

    dos = np.zeros((gridE,total))

    for i in range(len(fbuffer)):
        temp = fbuffer[i].split()
        a = int(temp[0])
        b = int(temp[1])
        c = float(temp[2])
        if c < 0.0001:
            c = 0.0001
        dos[a,b]=c # gridE x total

    dx, dy = sc/total, (EMAX-EMIN)/gridE
    #print dx, dy
    x = np.arange(0, sc, dx)
    y = np.arange(EMIN-base, EMAX-base, dy)
    #print x
    X,Y = np.meshgrid(x,y)
    plt.pcolormesh(X,Y,dos,norm=LogNorm(vmin=0.05,vmax=3.0)) # Ca2N
    #plt.pcolormesh(X,Y,dos,norm=LogNorm(vmin=0.1,vmax=10.0)) # Ag
    #plt.pcolormesh(X,Y,dos, vmin=0.0001,vmax=16.206)
    plt.axis([0,sc-dx,EMIN-base,EMAX-dy-base])
    plt.plot([0,sc-dx],[0,0], color='white', linestyle='--', linewidth=4)
    plt.ylabel(r'E - Fermi level (eV)')
    plt.xlabel(r'z-axis length ($\AA$)')
    fig.savefig('LDOS.png')
    #plt.plot([13.535292,13.535292],[EMIN-base,EMAX-dy-base])
    #plt.plot([16.0862337,16.0862337],[EMIN-base,EMAX-dy-base])
    #print np.max(dos)
    #print np.min(dos)


#plot_ldos()
#plt.show()
average_CHGCAR(input_file='LOCPOT')
#make_difference_CHGCAR('CHGCAR_yesNi', 'CHGCAR_noNi', 'CHGCAR_diff')