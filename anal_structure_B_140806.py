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


def analyzing_one_strucutre_info(input_name, target_atom, bonding_cutoff):
	tlist = making_real_coordi(input_name, 'test.out')
	name = 'bonding_length_angle_of_' + str(target_atom) + '.out'
	out = open(name,'w')

	#print np.array(tlist[0])

	f = open("test.out")
	fbuffer = f.readlines()

	bc = float(bonding_cutoff)
	ta = int(target_atom)

	ta_list = []
	ta_f = fbuffer[ta-1].split()
	for s in range(3):
		ta_list.append(float(ta_f[s]))
		
	### distnace check
	
	bond_list = []
	distance_list = []
	coordinate_list = []
	angle_list = []
	deg_list = []

	first_line = 'bonding_with' + '\t' + 'bond_length[aongstrom]' + '\n'
	out.write(first_line)


	for i in range(len(fbuffer)):
		
		ref_f = fbuffer[i].split()
		ref_list = []

		for s in range(3):
			ref_list.append(float(ref_f[s]))

		ref_array = np.array(ref_list)
		
		for j in range(3):
			new_j = -1 + j
			for k in range(3):
				new_k = -1 + k
				for l in range(3):
					new_l = -1 + l
					new_aaxis_array = (new_j)*np.array(tlist[0])
					new_baxis_array = (new_k)*np.array(tlist[1])
					new_caxis_array = (new_l)*np.array(tlist[2])
					new_ref_array = ref_array + new_aaxis_array + new_baxis_array + new_caxis_array
					distance = calculation_distance(ta_list, new_ref_array)
					if distance == 0:
						continue
					elif distance <= bc:
						coordinate_list.append(new_ref_array)
						distance_list.append(float(distance))
						bond_list.append(i+1)
						temp_line = 'B(' + str(i+1) + ')' + '\t' + str(distance) + '\n'
						out.write(temp_line)
						break
					else:
						continue

	#print bond_list
	#print distance_list
	#print coordinate_list

	
	second_line = 'bond_order' + '\t' + 'bond_angle[radian]' + '\t' + 'bond_angle[degree]' + '\n'
	out.write(second_line)
	##  bond angle check

	for i in range(len(bond_list)):
		for j in range(len(bond_list)):
			if i==j:
				break
			else:
				angle_rad = calculation_angle(coordinate_list[i],coordinate_list[j],ta_list)
				angle_list.append(angle_rad)
				angle_deg = transform_from_radian_to_degree(angle_rad)
				deg_list.append(angle_deg)
				temp_line = 'B(' + str(i+1) + ')-B(' + str(ta) + ')-B(' + str(j+1) + ')' + '\t' + str(angle_rad) + '\t' + str(angle_deg) + '\n'
				out.write(temp_line)


	coordination_number = len(bond_list)

	return coordination_number, deg_list

def analyzing_coordination_number(filename,total_number, bonding_cutoff):
	tn = int(total_number)
	bc = float(bonding_cutoff)
	summ = 0
	total_coordi_list = []

	for i in range(tn):
		temp_coordinate = analyzing_one_strucutre_info(filename,i+1,bc)
		total_coordi_list.append(int(temp_coordinate[0]))

	for i in range(len(total_coordi_list)):
		summ += total_coordi_list[i]


	average_sum = float(summ)/tn

	print summ
	print average_sum

def analyzing_angle_distribution(filename,total_number, bonding_cutoff, angle_spacing):
	tn = int(total_number)
	bc = float(bonding_cutoff)
	asing = float(angle_spacing)

	total_angle_list = []

	for i in range(tn):
		temp_angle = analyzing_one_strucutre_info(filename,i+1,bc)
		for j in range(len(temp_angle[1])):
			total_angle_list.append(float(temp_angle[1][j]))

	total_angle_list = np.sort(total_angle_list)
	
	countlist = []
	xlist = []
	angle = 0
	xlist.append(angle)

	while angle < 180:
		angle += asing
		xlist.append(angle)
		#print angle

	for i in range(len(xlist)-1):
		temp_count = 0
		print '--------------'+str(i+1)+'--------------'
		for j in range(len(total_angle_list)):
			if total_angle_list[j] >= float(xlist[i+1]):
				break
			elif float(xlist[i+1]) > total_angle_list[j] >= float(xlist[i]):
				temp_count += 1
			else:
				continue
		countlist.append(temp_count)

	new_xlist = xlist[1:]
	#print total_angle_list
	#print len(new_xlist)
	#print len(countlist)

	plt.plot(new_xlist,countlist)
	plt.axis([0,180,0,18])
	plt.xlabel('Angle (deg)')
	plt.ylabel('Count')
	plt.show()

#analyzing_coordination_number('CONTCAR',28, 2.0)
#analyzing_angle_distribution('CONTCAR2',26, 2.2, 1)
making_real_coordi('heet.txt', 'test_heet')
