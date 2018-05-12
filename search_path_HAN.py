import os.path
import time
from math import *
from math import sqrt
from numpy.linalg import *
import numpy as np
import matplotlib.pyplot as plt
import cmath

def tranform_direct_cartesian(aaxislist, baxislist, caxislist, atomlist1):
        atom_cart_list1=[]

        for i in range(0,3):
                temp1 = float(aaxislist[i])*float(atomlist1[0])
                temp2 = float(baxislist[i])*float(atomlist1[1])
                temp3 = float(caxislist[i])*float(atomlist1[2])
                temp = temp1 + temp2 + temp3
                atom_cart_list1.append(temp)

        return atom_cart_list1

def assign_atom_number_fixed_unitcell(initial_structure, final_structure, iterations, initial_out_name, final_out_name):
    '''
    Conditions: fixed cell & monoatomic system
    '''
    f = open(initial_structure, 'r')
    g = open(final_structure, 'r')

    tempf = f.readlines()
    tempg = g.readlines()

    f.close()
    g.close()

    f_univ = float(tempf[1].split()[0])
    f_a =  np.array([float(tempf[2].split()[i]) for i in range(3)]) * f_univ
    f_b =  np.array([float(tempf[3].split()[i]) for i in range(3)]) * f_univ
    f_c =  np.array([float(tempf[4].split()[i]) for i in range(3)]) * f_univ

    g_univ = float(tempg[1].split()[0])
    g_a =  np.array([float(tempg[2].split()[i]) for i in range(3)]) * g_univ
    g_b =  np.array([float(tempg[3].split()[i]) for i in range(3)]) * g_univ
    g_c =  np.array([float(tempg[4].split()[i]) for i in range(3)]) * g_univ

    f_atoms = int(tempf[6].split()[0])
    g_atoms = int(tempg[6].split()[0])

    thres = 1e-5

    if np.linalg.norm(f_a - g_a) > thres or np.linalg.norm(f_b - g_b) > thres or np.linalg.norm(f_c - g_c) > thres:
        print 'ERROR: Please check for lattice constants of initial and final structures'
        print np.linalg.norm(f_a - g_a)
        print np.linalg.norm(f_b - g_b)
        print np.linalg.norm(f_c - g_c)        
        #return 0
    if f_atoms != g_atoms:
        print 'ERROR: Please check for the number of atoms in initial and final structures'
        print f_atoms
        print g_atoms
        return 0

    f_atom_list = []
    g_atom_list = []

    for i in range(f_atoms):
        temp_f_atom = np.array([float(tempf[i+8].split()[j]) for j in range(3)])
        temp_g_atom = np.array([float(tempg[i+8].split()[j]) for j in range(3)])
        f_atom_list.append(temp_f_atom)
        g_atom_list.append(temp_g_atom)

    total_distance = cal_total_distance(f_atom_list,g_atom_list,f_a,f_b,f_c)
    print total_distance

    for i in range(iterations):
        test = True
        while test:
            x = np.random.randint(f_atoms,size=1)[0]
            y = np.random.randint(f_atoms,size=1)[0]
            if x != y:
                test = False
        #print x, y
        temp_f_list = []
        for j in range(f_atoms):
            temp_f_list.append(f_atom_list[j])
        temp_f_list[x], temp_f_list[y] = temp_f_list[y], temp_f_list[x]
        #print temp_f_list
        newdistance = cal_total_distance(temp_f_list,g_atom_list,f_a,f_b,f_c)
        #print newdistance
        if newdistance < total_distance:
            total_distance = newdistance
            f_atom_list = temp_f_list
            print total_distance

    print 'final distance = ' + str(total_distance)

    h = open(initial_out_name,'w')
    h.write('B'+'\n')
    h.write('1.0000'+'\n')
    h.write(str(f_a[0])+' '+str(f_a[1])+' '+str(f_a[2])+'\n')
    h.write(str(f_b[0])+' '+str(f_b[1])+' '+str(f_b[2])+'\n')
    h.write(str(f_c[0])+' '+str(f_c[1])+' '+str(f_c[2])+'\n')
    h.write('B'+'\n')
    h.write(str(f_atoms)+'\n')
    h.write('Direct'+'\n')
    for i in range(f_atoms):
        temp = str(f_atom_list[i][0])+' '+str(f_atom_list[i][1])+' '+str(f_atom_list[i][2])+'\n'
        h.write(temp)

    h.close()

    t = open(final_out_name,'w')
    t.write('B'+'\n')
    t.write('1.0000'+'\n')
    t.write(str(g_a[0])+' '+str(g_a[1])+' '+str(g_a[2])+'\n')
    t.write(str(g_b[0])+' '+str(g_b[1])+' '+str(g_b[2])+'\n')
    t.write(str(g_c[0])+' '+str(g_c[1])+' '+str(g_c[2])+'\n')
    t.write('B'+'\n')
    t.write(str(g_atoms)+'\n')
    t.write('Direct'+'\n')
    for i in range(g_atoms):
        temp = str(g_atom_list[i][0])+' '+str(g_atom_list[i][1])+' '+str(g_atom_list[i][2])+'\n'
        t.write(temp)

    t.close()

def cal_total_distance(atom_list1, atom_list2, latt_a, latt_b, latt_c):
    n = len(atom_list1)
    latt = np.array([latt_a,latt_b,latt_c])
    #print latt
    distance =0
    for i in range(n):
        diff = atom_list1[i] - atom_list2[i]
        diff = np.array(diff) - np.around(diff)
        diff = np.abs(diff)
        car_diff = np.dot(latt.transpose(), diff.transpose())
        #print car_diff
        distance += np.linalg.norm(car_diff)

    return distance


def assign_atom_number_unitcell(initial_structure, final_structure, iterations):
    '''
    Conditions: monoatomic system
    '''
    f = open(initial_structure, 'r')
    g = open(final_structure, 'r')

    tempf = f.readlines()
    tempg = g.readlines()

    f.close()
    g.close()

    f_univ = float(tempf[1].split()[0])
    f_a =  np.array([float(tempf[2].split()[i]) for i in range(3)]) * f_univ
    f_b =  np.array([float(tempf[3].split()[i]) for i in range(3)]) * f_univ
    f_c =  np.array([float(tempf[4].split()[i]) for i in range(3)]) * f_univ

    g_univ = float(tempg[1].split()[0])
    g_a =  np.array([float(tempg[2].split()[i]) for i in range(3)]) * g_univ
    g_b =  np.array([float(tempg[3].split()[i]) for i in range(3)]) * g_univ
    g_c =  np.array([float(tempg[4].split()[i]) for i in range(3)]) * g_univ

    f_atoms = int(tempf[6].split()[0])
    g_atoms = int(tempg[6].split()[0])

    f_atom_list = []
    g_atom_list = []

    delta_f = np.array([float(tempf[8].split()[j]) for j in range(3)]) - np.array([0,0,0])
    delta_g = np.array([float(tempg[8].split()[j]) for j in range(3)]) - np.array([0,0,0])

    for i in range(f_atoms):
        temp_f_atom = np.array([float(tempf[i+8].split()[j]) for j in range(3)]) - delta_f
        temp_g_atom = np.array([float(tempg[i+8].split()[j]) for j in range(3)]) - delta_g
        f_atom_list.append(temp_f_atom)
        g_atom_list.append(temp_g_atom)

    f_latt = np.array([f_a, f_b, f_c])
    g_latt = np.array([g_a, g_b, g_c])

    total_distance = cal_total_distance2(f_atom_list,g_atom_list,f_latt,g_latt)
    print total_distance

    for i in range(iterations):
        test = True
        while test:
            x = np.random.randint(f_atoms,size=1)[0]
            y = np.random.randint(f_atoms,size=1)[0]
            if x != y:
                test = False
        #print x, y
        temp_f_list = []
        for j in range(f_atoms):
            temp_f_list.append(f_atom_list[j])
        temp_f_list[x], temp_f_list[y] = temp_f_list[y], temp_f_list[x]
        #print temp_f_list
        newdistance = cal_total_distance2(temp_f_list,g_atom_list,f_latt,g_latt)
        #print newdistance
        if newdistance < total_distance:
            total_distance = newdistance
            f_atom_list = temp_f_list
            print 'i = '+str(i)+' distance = '+str(total_distance)

    print 'final distance = ' + str(total_distance)

    h = open('POSCAR_modified_initial','w')
    h.write('B'+'\n')
    h.write('1.0000'+'\n')
    h.write(str(f_a[0])+' '+str(f_a[1])+' '+str(f_a[2])+'\n')
    h.write(str(f_b[0])+' '+str(f_b[1])+' '+str(f_b[2])+'\n')
    h.write(str(f_c[0])+' '+str(f_c[1])+' '+str(f_c[2])+'\n')
    h.write('B'+'\n')
    h.write(str(f_atoms)+'\n')
    h.write('Direct'+'\n')
    for i in range(f_atoms):
        temp = str(f_atom_list[i][0])+' '+str(f_atom_list[i][1])+' '+str(f_atom_list[i][2])+'\n'
        h.write(temp)

    h.close()

    t = open('POSCAR_modified_final','w')
    t.write('B'+'\n')
    t.write('1.0000'+'\n')
    t.write(str(g_a[0])+' '+str(g_a[1])+' '+str(g_a[2])+'\n')
    t.write(str(g_b[0])+' '+str(g_b[1])+' '+str(g_b[2])+'\n')
    t.write(str(g_c[0])+' '+str(g_c[1])+' '+str(g_c[2])+'\n')
    t.write('B'+'\n')
    t.write(str(g_atoms)+'\n')
    t.write('Direct'+'\n')
    for i in range(g_atoms):
        temp = str(g_atom_list[i][0])+' '+str(g_atom_list[i][1])+' '+str(g_atom_list[i][2])+'\n'
        t.write(temp)

    t.close()


def cal_total_distance2(atom_list1, atom_list2, latt1, latt2):
    n = len(atom_list1)
    latta = np.array(latt1)
    lattb = np.array(latt2)
    #print latt
    distance =0
    for i in range(n):
        cara = np.dot(latta.transpose(), (atom_list1[i]).transpose())
        carb = np.dot(lattb.transpose(), (atom_list2[i]).transpose())
        temp = []
        for a in range(-1,2):
            for b in range(-1,2):
                for c in range(-1,2):
                    add = np.array([a,b,c]) * latta
                    car_diff = cara - carb + add.transpose()
                    tempdistance = np.linalg.norm(car_diff)
                    temp.append(tempdistance)
        distance += min(temp)
    return distance

def replace_atom(POSCAR_ref, POSCAR_change, change_list):

    change_list = np.array(change_list) - 1

    f = open(POSCAR_ref, 'r')
    g = open(POSCAR_change, 'r')

    tempf = f.readlines()
    tempg = g.readlines()

    f.close()
    g.close()

    f_univ = float(tempf[1].split()[0])
    f_a =  np.array([float(tempf[2].split()[i]) for i in range(3)]) * f_univ
    f_b =  np.array([float(tempf[3].split()[i]) for i in range(3)]) * f_univ
    f_c =  np.array([float(tempf[4].split()[i]) for i in range(3)]) * f_univ

    g_univ = float(tempg[1].split()[0])
    g_a =  np.array([float(tempg[2].split()[i]) for i in range(3)]) * g_univ
    g_b =  np.array([float(tempg[3].split()[i]) for i in range(3)]) * g_univ
    g_c =  np.array([float(tempg[4].split()[i]) for i in range(3)]) * g_univ

    f_atoms = int(tempf[6].split()[0])
    g_atoms = int(tempg[6].split()[0])

    f_atom_list = []
    g_atom_list = []

    for i in range(f_atoms):
        temp_f_atom = np.array([float(tempf[i+8].split()[j]) for j in range(3)])
        temp_g_atom = np.array([float(tempg[i+8].split()[j]) for j in range(3)])
        f_atom_list.append(temp_f_atom)
        g_atom_list.append(temp_g_atom)


    h = open('POSCAR_replace','w')
    h.write('B'+'\n')
    h.write('1.0000'+'\n')
    h.write(str(f_a[0])+' '+str(f_a[1])+' '+str(f_a[2])+'\n')
    h.write(str(f_b[0])+' '+str(f_b[1])+' '+str(f_b[2])+'\n')
    h.write(str(f_c[0])+' '+str(f_c[1])+' '+str(f_c[2])+'\n')
    h.write('B'+'\n')
    h.write(str(f_atoms)+'\n')
    h.write('Direct'+'\n')
    for i in range(f_atoms):
        if i in change_list:
            temp = str(g_atom_list[i][0])+' '+str(g_atom_list[i][1])+' '+str(g_atom_list[i][2])+'\n'
        else:
            temp = str(f_atom_list[i][0])+' '+str(f_atom_list[i][1])+' '+str(f_atom_list[i][2])+'\n'
        h.write(temp)

    h.close()

#clist = [1,18,61,79,25,73,49,31,37,67,7,8,55,43,13,20,26,62,50,32,68,80,44,9,38,56,74,21,2,14,27,33,51,81,63,45]
#clist = [1,18,61,79,25,73,49,31,37,67,7,8,55,43,13,20,26,62,50,32,68,80,44,9]
#clist = [39,57,75,69,22,15,28,3,52,82,34,10,40,58,64,76,70,16,23,46,4,11,35]
#clist = [1,18,61,79,25,73,49,31,37,67,7,8,55,43,13,20,26,62,50,32,68,80,44,9,38,56,74,21,2,14,27,33,51,81,63,45,39,57,75,69,22,15,28,3,52,82,34,10]
#clist = [45,39,57,75,22,15,69,28,3,52,82,34,10,40,81,70,63]
clist = [11,35,4,46,70,16,65,41]
print len(clist)
replace_atom('POSCAR_AA', 'POSCAR_BB',clist)

#assign_atom_number_fixed_unitcell('POSCAR_BB', 'POSCAR_AA',50000, 'POSCAR_m_BB', 'POSCAR_m_AA') 


#assign_atom_number_unitcell('POSCAR_ii', 'POSCAR_ff', 1000)
