import os.path
import time
from math import *
from math import sqrt
from numpy.linalg import norm
import numpy as np
import matplotlib.pyplot as plt
#from numeric import *

def transform_from_radian_to_degree(input_angle):
	ia = float(input_angle)
	output_angle = (180*ia) / np.pi

	return output_angle


def calculation_angle(first_list, second_list, center_list):
	firstar = np.array(first_list)
	secondar = np.array(second_list)
	centerar = np.array(center_list)

	vector1 = firstar - centerar
	vector2 = secondar - centerar
	vector3 = secondar - firstar

	mag_vector1 = np.linalg.norm(vector1)
	mag_vector2 = np.linalg.norm(vector2)
	mag_vector3 = np.linalg.norm(vector3)

	angle = acos(((mag_vector1)*(mag_vector1)+(mag_vector2)*(mag_vector2)-(mag_vector3)*(mag_vector3))/(2*(mag_vector1)*(mag_vector2)))

	return angle






def calculation_distance(first_list, second_list):
	first_x = float(first_list[0])
	first_y = float(first_list[1])
	first_z = float(first_list[2])
	second_x = float(second_list[0])
	second_y = float(second_list[1])
	second_z = float(second_list[2])

	distance = sqrt((first_x-second_x)*(first_x-second_x)+(first_y-second_y)*(first_y-second_y)+(first_z-second_z)*(first_z-second_z))

	return distance


def tranform_direct_cartesian(aaxislist, baxislist, caxislist, atomlist1):
        atom_cart_list1=[]

        for i in range(0,3):
                temp1 = float(aaxislist[i])*float(atomlist1[0])
                temp2 = float(baxislist[i])*float(atomlist1[1])
                temp3 = float(caxislist[i])*float(atomlist1[2])
                temp = temp1 + temp2 + temp3
                atom_cart_list1.append(temp)

        return atom_cart_list1

def making_real_coordi(input_name, output_name):
	f = open(input_name)
	fbuffer = f.readlines()

	g = open(output_name, 'w')

	aaxislist = []
	baxislist = []
	caxislist = []

	for i in range(3):
		af = fbuffer[2].split()
		bf = fbuffer[3].split()
		cf = fbuffer[4].split()
		aaxislist.append(float(af[i]))
		baxislist.append(float(bf[i]))
		caxislist.append(float(cf[i]))

	total_atom = fbuffer[6].split()
	info = total_atom[0]

	for i in range(int(info)):
		starting_pos = i + 8
		atomlist = fbuffer[starting_pos].split()
		new = tranform_direct_cartesian(aaxislist,baxislist,caxislist,atomlist)
		if i < int(info) - 1:
			templine = str(new[0]) + '\t' + str(new[1]) + '\t' + str(new[2]) + '\n'
			g.write(templine)
		else:
			templine = str(new[0]) + '\t' + str(new[1]) + '\t' + str(new[2])
			g.write(templine)

	return aaxislist, baxislist, caxislist


def making_snap_shot(input_name,eigen_name, displacement):
	f = open(input_name)
	fbuffer = f.readlines()
	h = open(eigen_name)
	hbuffer = h.readlines()
	displacement_name = int(displacement*10)

	g = open('POSCAR_'+str(displacement_name), 'w')

	aaxislist = []
	baxislist = []
	caxislist = []

	for i in range(7):
		g.write(fbuffer[i])
		#print i
	
	for j in range(3):
		af = fbuffer[2].split()
		aaxislist.append(float(af[j]))
		bf = fbuffer[3].split()
		baxislist.append(float(bf[j]))
		cf = fbuffer[4].split()
		caxislist.append(float(cf[j]))
		#print j


	g.write('Cartesian'+'\n')

	total_atom = fbuffer[6].split()
	info = total_atom[0]
	#print info
	#print aaxislist

	for i in range(int(info)):
		starting_pos = i + 8
		atomlist = fbuffer[starting_pos].split()
		new = tranform_direct_cartesian(aaxislist,baxislist,caxislist,atomlist)
		new_array = np.array(new)
		eigen = []
		for j in range(3):
			hline = hbuffer[i].split()
			eigen.append(float(hline[j]))
		eigen_array = np.array(eigen)
		modi_array = float(displacement)*eigen_array + new_array
		if i < int(info) - 1:
			templine = str(modi_array[0]) + '\t' + str(modi_array[1]) + '\t' + str(modi_array[2]) + '\n'
			g.write(templine)
		else:
			templine = str(modi_array[0]) + '\t' + str(modi_array[1]) + '\t' + str(modi_array[2])
			g.write(templine)



making_snap_shot('POSCAR','eigenvector_tetra222.txt',5.4)
