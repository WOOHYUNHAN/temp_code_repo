import os.path
import time
from math import *
from math import sqrt
from numpy.linalg import norm
import numpy as np
from sys import argv
#from numeric import *


def make_OUTCAR_and_CONTCAR_from_atomic_relaxation(input_file):


	f = open(input_file,'r')
	tempf = f.readlines()
	f.close()

	for i in range(len(tempf)):
		if str(tempf[i].split()[0]) == 'relax_done':
			standard_line = i

	Final_energy = float(tempf[standard_line-1].split()[3])
	Final_force = float(tempf[standard_line-1].split()[4])

	a_axis = np.array([float(tempf[standard_line+1].split()[i]) for i in range(3)])
	b_axis = np.array([float(tempf[standard_line+2].split()[i]) for i in range(3)])
	c_axis = np.array([float(tempf[standard_line+3].split()[i]) for i in range(3)])

	latt_vec = np.array([a_axis, b_axis, c_axis])

	stress = np.array([float(tempf[standard_line+4].split()[i]) for i in range(6)])

	recip_fac =  (2*np.pi) / np.dot(latt_vec[0], np.cross(latt_vec[1], latt_vec[2]))

	recip_vec = np.array([np.cross(latt_vec[1], latt_vec[2]), np.cross(latt_vec[2], latt_vec[0]), np.cross(latt_vec[0], latt_vec[1])]) * recip_fac

	atom_type = []
	position_cart = []
	position_dir = []
	force = []

	for i in range(standard_line+5, len(tempf)):
		tempff = tempf[i].split()
		atom_type.append(str(tempff[0]))
		position_cart.append([float(tempff[1]),float(tempff[2]),float(tempff[3])])
		position_dir.append(np.dot(np.linalg.inv(latt_vec.transpose()),np.array([float(tempff[1]),float(tempff[2]),float(tempff[3])])))
		force.append([float(tempff[4]),float(tempff[5]),float(tempff[6])])

	atom_type = np.array(atom_type)
	position_cart = np.array(position_cart)
	position_dir = np.array(position_dir)
	force = np.array(force)

	unique, indices = np.unique(atom_type, return_index=True)
	#print len(unique)
	#print unique, indices
	count = []

	for i in range(len(unique)):
		print len(unique) - 1
		if i == len(unique) - 1: 
			#print 'here'
			count.append(len(atom_type)-indices[i])
    	if i != len(unique) -1:
    		#print 'here2'
    		count.append(indices[i+1]-indices[i])
	
	'''
	Make CONTCAR
	'''

	g1 = open('CONTCAR', 'w')

	g1.write(str(Final_energy)+'\n')
	g1.write(str(1.00000000000)+'\n')

	for i in range(3):
		line = str(latt_vec[i][0]) + ' ' + str(latt_vec[i][1]) + ' ' + str(latt_vec[i][2]) + '\n'
		g1.write(line)

	line = ''
	for i in range(len(unique)):
		line += str(unique[i]) + ' '
	line += '\n'
	g1.write(line)

	line = ''
	for i in range(len(count)):
		line += str(count[i]) + ' '
	line += '\n'
	g1.write(line)


	g1.write('Direct'+'\n')

	for i in range(len(position_dir)):
		line = str(position_dir[i][0]) + ' ' + str(position_dir[i][1]) + ' ' + str(position_dir[i][2]) + '\n'
		g1.write(line)

	g1.close()

	'''
	Make OUTCAR
	'''
	'''
	  energy without entropy=      -73.33416531  energy(sigma->0) =      -73.33416531

    E-fermi :   0.0     XC(G=0): -11.8289     alpha+bet :-15.8746

  Total      -0.03514     0.06455    -0.06081    -0.03049     0.02155    -0.02658
  in kB     -56.30808   103.42654   -97.42605   -48.85257    34.52251   -42.58771
  external pressure =       16.77 kB  Pullay stress =        0.00 kB

  energy-cutoff  :        0.00
  volume of cell :       83.71
      direct lattice vectors                 reciprocal lattice vectors
     6.188022412  0.000000000  0.000000000     1.015378563  0.036071253 -0.036917548
    -0.238001779  6.699570566  0.000000000     0.000000000  0.937848963 -0.020808094
     0.071826217  0.044801846  2.019279839     0.000000000  0.000000000  3.111597107

  length of vectors
     6.188022412  6.703796731  2.021053508     1.016689563  0.938079770  0.938079770

 POSITION                                       TOTAL-FORCE (eV/Angst)
 -----------------------------------------------------------------------------------
      0.38048      0.49683     -0.38669        -0.000565     -0.008181     -0.004788
      5.02494      4.37338      0.81258         0.000790      0.002788      0.028374
      0.37899      4.51548      1.78852        -0.011923      0.006752      0.023863
      5.42764      2.97530      1.86313         0.005572      0.012735     -0.000901
      5.01193      1.38673      1.73187         0.018471      0.008542      0.013002
      2.59631      4.27157      0.60783        -0.001386      0.004131     -0.013822
      3.73270      5.19835      1.52392         0.004793     -0.023811     -0.014549
      3.39991      0.96976      0.66340        -0.000684     -0.014422     -0.002791
      2.55720      6.07258      0.56481        -0.007147     -0.001086     -0.020359
      1.66424      3.36117      1.69763         0.021780      0.013624     -0.003126
      1.93600      1.66163     -0.48451        -0.006612      0.013218      0.002590
      0.84527      5.92382      0.64950        -0.023089     -0.014290     -0.007491
 -----------------------------------------------------------------------------------
    total drift:                               -0.000000     -0.000000      0.000000



	'''
	h = open('OUTCAR','w')

	line = 'energy without entropy=      '+str(Final_energy)+ '  energy(sigma->0) =      ' + str(Final_energy)+'\n'
	h.write(line)

	line = '\n'
	h.write(line)

	line = 'E-fermi :   0.0     XC(G=0): -11.8289     alpha+bet :-15.8746' + '\n'
	h.write(line)

	line = '\n'
	h.write(line)

	eV2KB = 1600.0

	line = 'Total' + '\t' + str(stress[0]) + '\t' + str(stress[1]) + '\t' + str(stress[2]) + '\t' + str(stress[3]) + '\t' + str(stress[4]) + '\t' + str(stress[5]) + '\n'
	h.write(line)

	line = 'in kB' + '\t' + str(stress[0]*eV2KB) + '\t' + str(stress[1]*eV2KB) + '\t' + str(stress[2]*eV2KB) + '\t' + str(stress[3]*eV2KB) + '\t' + str(stress[4]*eV2KB) + '\t' + str(stress[5]*eV2KB) + '\n'
	h.write(line)

	line = 'external pressure =       '+ str((stress[0]+stress[1]+stress[2])*eV2KB/3.0) +'kB' +'   Pullay stress =        0.00 kB' + '\n'
	h.write(line)

	line = '\n'
	h.write(line)

	line = '  energy-cutoff  :        0.00' + '\n'
	h.write(line)

	line = 'volume of cell :       ' + str(np.dot(latt_vec[0], np.cross(latt_vec[1], latt_vec[2]))) + '\n'
	h.write(line)

	line = '      direct lattice vectors                 reciprocal lattice vectors' + '\n'
	h.write(line)


	for i in range(3):
		line = str(latt_vec[i][0]) + '\t' + str(latt_vec[i][1]) + '\t' + str(latt_vec[i][2]) + '\t' + str(recip_vec[i][0]) + '\t' + str(recip_vec[i][1]) + '\t' + str(recip_vec[i][2]) + '\n'
		h.write(line)

	line = '\n'
	h.write(line)

	line = 'length of vectors' + '\n'
	h.write(line)

	line = str(np.linalg.norm(latt_vec[0])) + '\t' + str(np.linalg.norm(latt_vec[1])) + '\t' + str(np.linalg.norm(latt_vec[2])) + '\t' + str(np.linalg.norm(recip_vec[0])) + '\t' + str(np.linalg.norm(recip_vec[1])) + '\t' + str(np.linalg.norm(recip_vec[2])) + '\n'
	h.write(line)

	line = '\n'
	h.write(line)

	line = 'POSITION                                       TOTAL-FORCE (eV/Angst)' + '\n'
	h.write(line)

	line = '-----------------------------------------------------------------------------------' + '\n'
	h.write(line)

	for i in range(len(position_cart)):
		line = str(position_cart[i][0]) + '\t' + str(position_cart[i][1]) + '\t' + str(position_cart[i][2]) + '\t' + str(force[i][0]) + '\t' + str(force[i][1]) + '\t' + str(force[i][2]) + '\n'
		h.write(line)

	line = '-----------------------------------------------------------------------------------' + '\n'
	h.write(line)

	line = 'total drift:                               -0.000000     -0.000000      0.000000' + '\n'
	h.write(line)



	h.close()

make_OUTCAR_and_CONTCAR_from_atomic_relaxation('out')