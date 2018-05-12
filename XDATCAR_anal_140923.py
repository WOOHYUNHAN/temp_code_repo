import os.path
import time
from math import *
from math import sqrt
from numpy.linalg import norm
import numpy as np
import matplotlib.pyplot as plt
#from numeric import *

def XDATCAR_cut(start_relax, end_relax):
        f=open('XDATCAR')
        fbuffer = f.readlines()
        number = fbuffer[6].split()
        total_number = int(number[0])
        name = 'XDATCAR_'+str(start_relax)+'_'+str(end_relax)
        g=open(name,'w')
        g.write(fbuffer[0])
        g.write(fbuffer[1])
        g.write(fbuffer[2])
        g.write(fbuffer[3])
        g.write(fbuffer[4])
        g.write(fbuffer[5])
        g.write(fbuffer[6])
        for i in range(start_relax, end_relax+1):
                for j in range(total_number+1):
                        target_line = (int(i)-1)*(total_number+1) + 7 + j
                        g.write(fbuffer[target_line])
        g.close()
        

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

def making_real_coordi(relax_number, output_name):
	f = open('XDACTCAR_30000_32000')
	fbuffer = f.readlines()
	rn = relax_number
	number=[]
	number = fbuffer[6].split()
	total_number = int(number[0])
	target_line = (int(rn)-1)*(total_number+1) + 7

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


	for i in range(1,total_number+1):
		starting_pos = i + target_line
		atomlist = fbuffer[starting_pos].split()
		new = tranform_direct_cartesian(aaxislist,baxislist,caxislist,atomlist)
		if i < total_number:
			templine = str(new[0]) + '\t' + str(new[1]) + '\t' + str(new[2]) + '\n'
			g.write(templine)
		else:
			templine = str(new[0]) + '\t' + str(new[1]) + '\t' + str(new[2])
			g.write(templine)

	return aaxislist, baxislist, caxislist, total_number

def test_direct_coordi(test):
        if -0.5 <= test <= 0.5:
                final = test
        elif test > 0.5:
                final = test - 1
        elif test < -0.5:
                final = test + 1
        return final


def analyzing_difference_direct_coordi(start_relax, end_relax, filename):
        f = open(filename)
        fbuffer = f.readlines()
        number = fbuffer[6].split()
        total_number = int(number[0])

        g = open('direct_coordi', 'w')

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

        for i in range(start_relax, end_relax+1):
                for j in range(1, total_number+1):
                        target_line1 = (int(i)-1)*(total_number+1) + 7 + j
                        target_line2 = (int(i+1)-1)*(total_number+1) + 7 + j
                        t1 = fbuffer[target_line1].split()
                        t2 = fbuffer[target_line2].split()
                        tx = float(t2[0])-float(t1[0])
                        ty = float(t2[1])-float(t1[1])
                        tz = float(t2[2])-float(t1[2])
                        tx = test_direct_coordi(tx)
                        ty = test_direct_coordi(ty)
                        tz = test_direct_coordi(tz)
                        if j < total_number:
                                templine = str(tx) + '\t' + str(ty) + '\t' + str(tz) + '\n'
                                g.write(templine)
                        else:
                                templine = str(tx) + '\t' + str(ty) + '\t' + str(tz) + '\n' + '\n'
                                g.write(templine)
        g.close()
        return aaxislist, baxislist, caxislist, total_number
        
def analyzing_atomic_distance(start_relax, end_relax, atom_number):
        a = analyzing_difference_direct_coordi(start_relax, end_relax,'XDATCAR_27000_32000')
        aaxislist = a[0]
        baxislist = a[1]
        caxislist = a[2]
        total_number = a[3]
        f = open('direct_coordi')
        fbuffer = f.readlines()
        name = 'distance_'+str(atom_number)+'_from_'+str(start_relax)+'_to_'+str(end_relax)
        g= open(name, 'w')
        for i in range(start_relax,end_relax+1):
                target_line = (int(i)-1)*(total_number+1) + atom_number - 1
                direct = fbuffer[target_line].split()
                trans_direct = tranform_direct_cartesian(aaxislist, baxislist, caxislist, direct)
                trans_direct = np.array(trans_direct)
                distance = np.linalg.norm(trans_direct)
                temp_line = str(i) + '\t' + str(distance)+'\n'
                g.write(temp_line)
                print 'relax number is ' + str(i)

        g.close()                        





        

def making_real_coordi_line(relax_number, atom_number):
	f = open('XDACTCAR_30000_34000')
	fbuffer = f.readlines()
	rn = relax_number
	number=[]
	number = fbuffer[6].split()
	total_number = int(number[0])
	target_line = (int(rn)-1)*(total_number+1) + 7

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

	starting_pos = atom_number + target_line
	atomlist = fbuffer[starting_pos].split()
	new = tranform_direct_cartesian(aaxislist,baxislist,caxislist,atomlist)
	return new




    

def analyzing_average_coord(start_relax, end_relax, bonding_cutoff):
	sr = start_relax
	er = end_relax
	bc = float(bonding_cutoff)

	final_data = []

	for i in range(sr,er+1):
		temp_data = []
		tlist = making_real_coordi(i,'test.out')
		total_number = tlist[3]
		f = open("test.out")
		fbuffer = f.readlines()
		for a in range(total_number):
			ta_list = []
			ta_f = fbuffer[a].split()
			
			for b in range(3):
				ta_list.append(float(ta_f[b]))
			
			index = 0
			for c in range(len(fbuffer)):
				ref_f = fbuffer[c].split()
				ref_list = []
				for d in range(3):
					ref_list.append(float(ref_f[d]))
				ref_array = np.array(ref_list)
				for j in range(-1,2):
					for k in range(-1,2):
						for l in range(-1,2):
							new_aaxis_array = (j)*np.array(tlist[0])
							new_baxis_array = (k)*np.array(tlist[1])
							new_caxis_array = (l)*np.array(tlist[2])
							new_ref_array = ref_array + new_aaxis_array + new_baxis_array + new_caxis_array
							distance = calculation_distance(ta_list, new_ref_array)
							if distance == 0:
								continue
							elif distance <= bc:
								index += 1
			temp_data.append(index)
		temp_array = np.array(temp_data)
		#print len(temp_array)
		avg = np.average(temp_array)
		final_data.append(avg)
		print 'relax_number is '+str(i)
	#print final_data
	g = open('average_coordi from '+str(sr)+' to '+str(er), 'w')
	for i in range(len(final_data)):
		temp_line = str(i+1) + '\t' + str(final_data[i]) + '\n'
		g.write(temp_line)


def making_set_information(start_relax, end_relax):
        #ico1 = [1, 23, 15, 25, 26, 16, 6, 20, 27, 13, 21, 12]
        #ico1 = [11, 2, 16, 3, 19, 8, 1, 24, 12, 23, 13, 21]
        #ico1 = [1,97,177,113,155,209,89,123,43,201,193,161]
        ico1 =  [222,112,120,216,38,80,32,160,174,44,72,208]
        #ico2 = [11, 17, 18, 28,8,22, 7,10, 19,3,14,5]
        #ico2 = [27, 5, 10, 15, 4, 20, 6, 26, 9, 14, 28, 22]
        #ico2 = [108, 52,148,132,59,139,80,175,24,40,223,80]
        ico2 = [161,89,177,97,1,185,85,121,9,57,145,17]
        dum1 = [4, 9]
        dum2 = [2, 24]
        f = open('XDATCAR_5000_15000_2x2x2')
        fbuffer = f.readlines()
        number = fbuffer[6].split()
        total_number = int(number[0])
        ico1list = np.array([0,0,0])
        ico2list = np.array([0,0,0])
        dum1list = np.array([0,0,0])
        dum2list = np.array([0,0,0])

        g = open('XDATCAR_set','w')
        g.write(fbuffer[0])
        g.write(fbuffer[1])
        g.write(fbuffer[2])
        g.write(fbuffer[3])
        g.write(fbuffer[4])
        g.write(fbuffer[5])
        g.write(str(4)+'\n')
        for i in range(start_relax, end_relax+1):
                for j in ico1:
                        templist = []
                        target_line = (int(i)-1)*(total_number+1) + 7 + j
                        temp = fbuffer[target_line].split()
                        templist.append(float(temp[0]))
                        templist.append(float(temp[1]))
                        templist.append(float(temp[2]))
                        ico1list = ico1list + templist
                for k in ico2:
                        templist = []
                        target_line = (int(i)-1)*(total_number+1) + 7 + k
                        temp = fbuffer[target_line].split()
                        templist.append(float(temp[0]))
                        templist.append(float(temp[1]))
                        templist.append(float(temp[2]))
                        templist = np.array(templist)
                        ico2list = ico2list + templist
                for l in dum1:
                        templist = []
                        target_line = (int(i)-1)*(total_number+1) + 7 + l
                        temp = fbuffer[target_line].split()
                        templist.append(float(temp[0]))
                        templist.append(float(temp[1]))
                        templist.append(float(temp[2]))
                        templist = np.array(templist)
                        dum1list = dum1list + templist
                for m in dum2:
                        templist = []
                        target_line = (int(i)-1)*(total_number+1) + 7 + m
                        temp = fbuffer[target_line].split()
                        templist.append(float(temp[0]))
                        templist.append(float(temp[1]))
                        templist.append(float(temp[2]))
                        templist = np.array(templist)
                        dum2list = dum2list + templist
                print 'relax = '+str(i)
                ico1list = ico1list / len(ico1)
                ico2list = ico2list / len(ico2)
                dum1list = dum1list / len(dum1)
                dum2list = dum2list / len(dum2)
                g.write('\n')
                g.write(str(ico1list[0])+'\t'+str(ico1list[1])+'\t'+str(ico1list[2])+'\n')
                g.write(str(ico2list[0])+'\t'+str(ico2list[1])+'\t'+str(ico2list[2])+'\n')
                g.write(str(dum1list[0])+'\t'+str(dum1list[1])+'\t'+str(dum1list[2])+'\n')
                g.write(str(dum2list[0])+'\t'+str(dum2list[1])+'\t'+str(dum2list[2])+'\n')
        g.close()

def analyzing_com_distance(start_relax, end_relax, atom_number):
        a = analyzing_difference_direct_coordi(start_relax, end_relax,'XDATCAR_set')
        aaxislist = a[0]
        baxislist = a[1]
        caxislist = a[2]
        total_number = a[3]
        f = open('direct_coordi')
        fbuffer = f.readlines()
        name = 'distance_'+str(atom_number)+'_from_'+str(start_relax)+'_to_'+str(end_relax)
        g= open(name, 'w')
        for i in range(start_relax,end_relax+1):
                target_line = (int(i)-1)*(total_number+1) + atom_number - 1
                direct = fbuffer[target_line].split()
                trans_direct = tranform_direct_cartesian(aaxislist, baxislist, caxislist, direct)
                trans_direct = np.array(trans_direct)
                distance = np.linalg.norm(trans_direct)
                temp_line = str(i) + '\t' + str(distance)+'\n'
                g.write(temp_line)
                print 'relax number is ' + str(i)
        g.close()

def get_rdf(input_name, relax_number, element1, element2, r_max, bin_spacing, multi):
        f= open(input_name)
        fbuffer = f.readlines()
        number = fbuffer[6].split()
        total_number = int(number[0])
        
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

        target1 = (int(relax_number)-1)*(total_number+1) + 7 + element1
        target2 = (int(relax_number)-1)*(total_number+1) + 7 + element2
        f1 = fbuffer[target1].split()
        f2 = fbuffer[target2].split()
        
        pos1 = []
        pos2 = []

        volume = np.dot(aaxislist, np.cross(baxislist, caxislist))

        for i in range(3):
                pos1.append(float(f1[i]))
                pos2.append(float(f2[i]))
        pos1 = np.array(pos1)
        pos2 = np.array(pos2)
        #print pos2
        
        multi_atoms = []
        nx, ny, nz = multi
        for x in range(nx):
                x = x - 1 
                for y in range(ny):
                        y = y - 1
                        for z in range(nz):
                                z = z - 1
                                multi_atoms.append(pos2 + np.array([x, y, z]))
        r = []
        #print len(multi_atoms)
        for i in range(len(multi_atoms)):
                diff = pos1 - multi_atoms[i]
                if np.linalg.norm(diff) == 0:
                        continue
                diff = diff - np.round(diff/multi)*multi
                diff = tranform_direct_cartesian(aaxislist, baxislist, caxislist, diff)
                diff = np.linalg.norm(diff)
                r.append(diff)
        #print r
        hist, bin_edges = np.histogram(r, bins = int(r_max / bin_spacing), range = (0, r_max), density = True)
        rdf = hist / 4. / np.pi / (bin_edges + bin_spacing / 2.)[:-1] ** 2 * volume / bin_spacing
        if (np.array(r) < r_max).any() :
                hist, bin_edges = np.histogram(r, bins = int(r_max / bin_spacing), range = (0, r_max), density = True)
                rdf = hist / 4. / np.pi / (bin_edges + bin_spacing / 2.)[:-1] ** 2 * volume / bin_spacing
                return bin_edges[:-1], rdf
        else :
                bin = np.linspace(0, r_max, num = r_max / bin_spacing , endpoint = False)
                rdf = np.zeros((len(bin)))
                return bin, rdf

def get_pair_corr(input_name, relax_number, set_list, r_max, bin_spacing, multi):
        bin = np.linspace(0, r_max, num = r_max / bin_spacing , endpoint = False)
        pair_corr = np.zeros((len(bin)))
        for i in set_list:
                for j in set_list:
                        element1=i
                        element2=j
                        r, rdf = get_rdf(input_name, relax_number, element1, element2, r_max, bin_spacing, multi)
                        pair_corr = pair_corr + rdf
        return r, pair_corr
                        
def cos_dist(array1, array2):

    array1_norm = np.linalg.norm(array1)
    array2_norm = np.linalg.norm(array2)

    dist  = 0.5 * (1. - (np.array(array1) * np.array(array2)).sum() / array1_norm / array2_norm)
    return dist


def making_PCF_XDAT(input_name,output_name, set_list, r_max, bin_spacing, multi, start_relax, end_relax):
        g = open(output_name,'w')
        reference_PCF = np.array([      0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,   54083.5325232 ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,  100688.08565189,
             0.        ,   65577.88300882,   97236.18074801,
        192249.8384691 ,  126710.13263082,  187917.58681456,
             0.        ,   91865.06985294,  181688.70954893,
         89840.55745189,   88853.30607114,   58588.15931788,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,   24625.14480951,
             0.        ,   72836.77692044,  102461.32789861,
         77804.86029873,       0.        ,   76726.12312032,
         17583.49571183,   58207.4352863 ,   43355.53902734,
         68893.73464683,   34211.73581032,       0.        ,
             0.        ,   33520.60744125,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,   48299.59360613,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,    6618.86555921,
             0.        ,   15259.88990212,       0.        ,
         12924.88589568,    6424.20372438,       0.        ,
             0.        ,       0.        ,    8366.12173612,
             0.        ,       0.        ,       0.        ,
             0.        ,   12189.62588958,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,   23160.48450547,   23030.73414722,
             0.        ,   22774.4831881 ,       0.        ,
         22522.48531318,       0.        ,       0.        ,
         14768.17265491,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,   20006.17076797,   19901.97213901,
         19798.58544559,    4924.00056847,       0.        ,
          4873.30334555,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,   11926.10255321,
             0.        ,   11809.03722433,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
         11357.59875375,   11302.99499266,       0.        ,
             0.        ,    8356.1442825 ,    5544.23545052,
             0.        ,       0.        ,       0.        ,
         10879.99666177,       0.        ,   10777.95739736,
             0.        ,       0.        ,       0.        ,
             0.        ,   15793.57833879,   30130.88170256,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,   10194.817992  ,
         15222.5586185 ,       0.        ,   10056.42840226,
         40043.69703882,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,    3621.65825323,
             0.        ,       0.        ,       0.        ,
          3558.05188108,       0.        ,       0.        ,
          4681.92046797,       0.        ,       0.        ,
             0.        ,    6901.40613793,       0.        ,
         10262.77103377,       0.        ,       0.        ,
             0.        ,    4483.34644385,    8928.37373116,
             0.        ,       0.        ,   13222.31791274,
         19749.43680739,       0.        ,   26110.60167257,
             0.        ,       0.        ,       0.        ,
             0.        ,    4261.29608982,    4243.54070875,
             0.        ,       0.        ,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,    3079.05622293,       0.        ,
             0.        ,       0.        ,       0.        ,
             0.        ,    3004.64069573,       0.        ,
             0.        ,       0.        ,       0.        ,
         11778.67734365,       0.        ])
        for i in range(start_relax, end_relax+1):
                a = get_pair_corr(input_name, i, set_list, r_max, bin_spacing, multi)
                fin_1 = a[1]-1.
                fin_2 = reference_PCF - 1.
                fin_1 = fin_1/np.linalg.norm(fin_1)
                fin_2 = fin_2/np.linalg.norm(fin_2)
                dist = cos_dist(fin_1,fin_2)
                line = str(i)+'\t'+str(dist)+'\n'
                g.write(line)
                print 'relax is ' + str(i)
        g.close()
                        
        
        
#analyzing_average_coord(1, 1, 2.0)
#analyzing_average_coord(1,3000, 2.0)
#making_real_coordi(10000,'test.out')
#making_set_information(1, 10000)
#analyzing_com_distance(1, 9999, 1)
#analyzing_com_distance(1, 9999, 2)
#analyzing_com_distance(1, 9999, 3)
#analyzing_com_distance(1, 9999, 4)
#XDATCAR_cut(13000, 14000)
#print get_rdf('XDATCAR_35000_35001', 1, 2, 5, 5, 0.1, np.array([2, 2, 2]))
#print get_pair_corr('XDATCAR_35000_35001', 1, [1, 23, 15, 25, 26, 16, 6, 20, 27, 13, 21, 12], 5, 0.01, np.array([2, 2, 2]))
making_PCF_XDAT('XDATCAR_13000_14000','PCF_XDAT_1', [11, 2, 16, 3, 19, 8, 1, 24, 12, 23, 13, 21], 5, 0.01, np.array([2, 2, 2]), 1, 1000)
making_PCF_XDAT('XDATCAR_13000_14000','PCF_XDAT_2', [27, 5, 10, 15, 4, 20, 6, 26, 9, 14, 28, 22], 5, 0.01, np.array([2, 2, 2]), 1, 1000)
