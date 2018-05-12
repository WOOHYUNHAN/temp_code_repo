import os.path
import time
from math import *
from math import sqrt
from numpy.linalg import norm
import numpy as np
import matplotlib.pyplot as plt
#from numeric import *


def cal_distance(latt_a,latt_b,latt_c,pos):
    cart = np.dot(np.array([latt_a,latt_b,latt_c]), pos)
    #cart = pos[0] * np.array(latt_a) + pos[1] * np.array(latt_b) + pos[2] * np.array(latt_c)
    distance = np.linalg.norm(cart)
    return distance

def transform_from_radian_to_degree(input_angle):
	ia = float(input_angle)
	output_angle = (180*ia) / np.pi

	return output_angle

def get_atomic_data_from_CONTCAR(filename='CONTCAR'):
    ### read lattice vector and atomic information from CONTCAR
    f = open(filename)
    fbuffer = f.readlines()
    f.close()

    univ = float(fbuffer[1].split()[0])
    latt_a = np.array([float(fbuffer[2].split()[i]) for i in range(3)]) * univ
    latt_b = np.array([float(fbuffer[3].split()[i]) for i in range(3)]) * univ
    latt_c = np.array([float(fbuffer[4].split()[i]) for i in range(3)]) * univ

    vol = np.dot(latt_a, np.cross(latt_b,latt_c))

    recip_a = (2*np.pi/vol) * np.cross(latt_b, latt_c)
    recip_b = (2*np.pi/vol) * np.cross(latt_c, latt_a)
    recip_c = (2*np.pi/vol) * np.cross(latt_a, latt_b)

    atom_name = np.array([str(fbuffer[5].split()[i]) for i in range(len(fbuffer[5].split()))])
    atom_number = np.array([int(fbuffer[6].split()[i]) for i in range(len(fbuffer[5].split()))])


    #print np.linalg.norm(recip_a)/np.linalg.norm(recip_b)
    #print np.linalg.norm(recip_a)/np.linalg.norm(recip_c)
    #print recip_c

    return atom_name, atom_number, latt_a, latt_b, latt_c, recip_a, recip_b, recip_c

def cal_dos_band(Ef,nband,ngrid,smearing_type,sigma,area_type):
    f = open('FSband.BXSF')
    fbuffer = f.readlines()
    f.close()

    #if smearing_type == 0:
    #    print "we use gaussian smearing method with sigma of " + str(sigma) + " eV"
    #else:
    #    print "check smearing type!"

    nk = np.array([int(fbuffer[15].split()[i]) for i in range(len(fbuffer[15].split()))])

    nx, ny, nz = int(nk[0]), int(nk[1]), int(nk[2])

    #print nx,ny,nz


    start = 0
    for i in range(len(fbuffer)):
        temp = fbuffer[i].split()
        if len(temp)==0:
            pass
        elif temp[0] == 'BAND:' and temp[1] == str(nband):
            start = i
            break
    if start == 0:
        print "Cannot find proper bands"
        return 0

    #eigenval_array = np.zeros(nx*ny*nz)
    eigenval = []

    #print start
    for i in range(nx*ny*nz):
        temp = fbuffer[start+1+i].split()
        if len(temp)==0:
            pass
        elif temp[0] == 'BAND:' or temp[0] == 'END_BANDGRID_3D' :
            #print i
            #print temp[0] == ''
            break
        else:
            #print 'ok'
            for j in range(nz):
                eigenval.append(float(temp[j]))

    #print len(eigenval)
    #print nx*ny*nz

    eigenval = np.array(eigenval) 

    #maxvalue = np.max(eigenval) + 0.5
    #minvalue = np.min(eigenval) - 0.5

    #print np.max(eigenval)
    #print np.min(eigenval)

    maxvalue = 5.0
    minvalue = -5.0

    #print maxvalue
    #print minvalue
    #x = np.asarray(range(1000))*0.01
    #print x
    egrid = []
    spacing = (maxvalue - minvalue)/(ngrid-1)
    #print spacing
    for i in range(ngrid):
        egrid.append(minvalue+spacing*i)
    egrid=np.array(egrid)
    #print egrid
    #print eigenval[2]

    total = np.zeros(ngrid)
    for i in range(len(eigenval)):
        temp = eigenval[i] - egrid
        #print temp
        g_temp = (1.0/sigma/sqrt(2*pi))*np.exp(-0.5*((temp)/sigma)**2)
        total += g_temp

    #normalized = 2 / (np.sum(spacing * total)) # ISPIN = 1
    normalized = 2 / (np.sum(spacing * total)) * 0.5 # ISPIN 2
    #total = total * normalized
    total = total * normalized

    #print np.sum(spacing * total)

    area = 0
    for i in range(len(egrid)):
        if area_type == 1: # upper part of Fermi level
            if float(egrid[i]) >= 0:
                area += spacing * total[i]
        elif area_type == -1: # lower part of Fermi level
            if float(egrid[i]) <= 0:
                area += spacing * total[i]  
        else:
            print 'ERROR: area_type'          

    DOS_EF = 0
    for i in range(len(eigenval)):
        temp = eigenval[i] - 0
        #print temp
        g_temp = (1.0/sigma/sqrt(2*pi))*np.exp(-0.5*((temp)/sigma)**2)
        DOS_EF += g_temp

    print 'DOS for Ef = 0 is ', DOS_EF * normalized
    #print 'DOS for Ef = 0 is ', total[(np.abs(egrid-0.0)).argmin()]
    print 'selected area is =', area
    
    g = open('band_projected_DOS','w')
    for i in range(len(egrid)):
        line = str(egrid[i]) + '\t' + "%0.9f" %total[i] + '\n'
        g.write(line)

    g.close()
    plt.plot(egrid,total)
    
    #plt.axis([0,0.2,0,2])
    plt.show()
    return 0


    #print start



#cal_dos_band(0.16,36,0,0.005)
#g = open('test','w')
#for i in range(500):
#    initial = 0.0
#    final = 0.2
#    spacing = (final-initial)/500
#    energy = initial + i*spacing
#    #print str(energy) + '\t' + str(cal_dos_band(energy,36,0,0.01))
#    g.write(str(energy) + '\t' + str(cal_dos_band(energy,37,0,0.005)) + '\n')
#g.close()

def EIGENVAL_DOS(Ef,nband,ngrid,smearing_type,sigma,area_type):
    f = open('EIGENVAL')
    fbuffer = f.readlines()
    f.close()

    info = np.array([int(fbuffer[5].split()[i]) for i in range(len(fbuffer[5].split()))])

    x, nk, nb = int(info[0]), int(info[1]), int(info[2])

    weight = []
    eigenval = []
    for i in range(nk):
        line = 7 + (nb+2)*i
        w = float(fbuffer[line].split()[3])
        weight.append(w)
        eigenval.append(float(fbuffer[line+nband].split()[1])-Ef)

    eigenval = np.array(eigenval)
    weight = np.array(weight)

    #maxvalue = np.max(eigenval) + 0.5
    #minvalue = np.min(eigenval) - 0.5
    maxvalue = 5.0
    minvalue = -5.0
    #print np.max(eigenval) - np.min(eigenval)

    #print maxvalue
    #print minvalue
    #x = np.asarray(range(1000))*0.01
    #print x
    egrid = []
    spacing = (maxvalue - minvalue)/(ngrid-1)
    for i in range(ngrid):
        egrid.append(minvalue+spacing*i)
    egrid=np.array(egrid)

    total = np.zeros(ngrid)
    for i in range(nk):
        temp = eigenval[i] - egrid
        #print temp
        g_temp = (1.0/sigma/sqrt(2*pi))*np.exp(-0.5*((temp)/sigma)**2) * weight[i]
        total += g_temp

    print (np.sum(spacing * total))
    normalized = 2 / (np.sum(spacing * total))
    total = total * normalized

    print normalized

    area = 0
    for i in range(len(egrid)):
        if area_type == 1: # upper part of Fermi level
            if float(egrid[i]) >= 0:
                area += spacing * total[i]
        elif area_type == -1: # lower part of Fermi level
            if float(egrid[i]) <= 0:
                area += spacing * total[i]  
        else:
            print 'ERROR: area_type'   
    DOS_EF = 0
    for i in range(len(eigenval)):
        temp = eigenval[i] - 0
        #print temp
        g_temp = (1.0/sigma/sqrt(2*pi))*np.exp(-0.5*((temp)/sigma)**2)
        DOS_EF += g_temp

    print 'DOS for Ef = 0 is ', DOS_EF * normalized
    #print 'DOS for Ef = 0 is ', total[(np.abs(egrid-0.0)).argmin()]
    print 'selected area is =', area

    g = open('band_projected_DOS','w')
    for i in range(len(egrid)):
        line = str(egrid[i]) + '\t' + "%0.9f" %total[i] + '\n'
        g.write(line)

    g.close()

    plt.plot(egrid,total)
    #print np.sum(spacing * total)
    #plt.axis([0,0.2,0,2])
    plt.show()
    return 0



        









EIGENVAL_DOS(7.6896,101,10000,0,0.05,1)
#cal_dos_band(0.0, 107 ,10000,0,0.05, -1)
#g = open('test','w')
#for i in range(500):
#    initial = -1.0
#    final = 1.0
#    spacing = (final-initial)/500
#    energy = initial + i*spacing
#    #print energy
#    gvalue = np.exp(-0.5*((0.1-energy)/0.05)**2)
#    #print gvalue
#    g.write(str(energy) + '\t' + str(gvalue) + '\n')
#    #print str(energy) + '\t' + str(cal_dos_band(energy,36,0,0.01))
#    #g.write(str(energy) + '\t' + str(EIGENVAL_DOS(energy,36,0,0.05)) + '\n')
#g.close()