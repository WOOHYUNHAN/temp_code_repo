import os.path
import time
from math import *
from math import sqrt
import numpy as np
#from numeric import *

def tranform_direct_cartesian(aaxislist, baxislist, caxislist, atomlist1):
        atom_cart_list1=[]

        for i in range(0,3):
                temp1 = float(aaxislist[i])*float(atomlist1[0])
                temp2 = float(baxislist[i])*float(atomlist1[1])
                temp3 = float(caxislist[i])*float(atomlist1[2])
                temp = temp1 + temp2 + temp3
                atom_cart_list1.append(temp)

        return atom_cart_list1

def POSCAR_maker(input_name, relax_number, output_name):
    rn = relax_number
    f = open(input_name)
    fbuffer = f.readlines()
    atom=[]
    number=[]

    number = fbuffer[6].split()

    total_number = int(number[0])
    
    out = open(output_name, 'w')


    for i in xrange(1,8):
        temp_position = i - 1
        write_line = fbuffer[temp_position].split()
        temp_line = ''
        for j in xrange(len(write_line)):
            if j == len(write_line) -1 :
                temp_line += ' ' + write_line[j] + '\n'
            else:
                temp_line += ' ' + write_line[j]
        out.write(temp_line)

    direct_line = 'Direct'+'\n'
    out.write(direct_line)

    target_line = (int(rn)-1)*(total_number+1) + 7

    for i in xrange(1,total_number+1):
        if i == total_number:
            temp_position = target_line + i
            direct = fbuffer[temp_position].split()
            x = direct[0]
            y = direct[1]
            z = direct[2]
            templine = str(x)+' '+str(y)+' '+str(z)
            out.write(templine)
        else:
            temp_position = target_line + i
            direct = fbuffer[temp_position].split()
            x = direct[0]
            y = direct[1]
            z = direct[2]
            templine = str(x)+' '+str(y)+' '+str(z)+'\n'
            out.write(templine)


def POSCAR_maker_2(input_name, relax_number, output_name):
    rn = relax_number
    f = open(input_name)
    fbuffer = f.readlines()
    atom=[]
    number=[]

    number = fbuffer[6].split()

    total_number = int(number[0])
    
    out = open(output_name, 'w')

    for i in xrange(len(fbuffer)):
        temp = fbuffer[i].split()
        search_1 = temp[0]
        if search_1 == 'Direct':
                search_2 = temp[2]
                if search_2 == relax_number:
                        target_line = i

    for i in xrange(1,8):
        temp_position = i - 1
        write_line = fbuffer[temp_position].split()
        temp_line = ''
        for j in xrange(len(write_line)):
            if j == len(write_line) -1 :
                temp_line += ' ' + write_line[j] + '\n'
            else:
                temp_line += ' ' + write_line[j]
        out.write(temp_line)

    direct_line = 'Direct'+'\n'
    out.write(direct_line)


    for i in xrange(1,total_number+1):
        if i == total_number:
            temp_position = target_line + i
            direct = fbuffer[temp_position].split()
            x = direct[0]
            y = direct[1]
            z = direct[2]
            templine = str(x)+' '+str(y)+' '+str(z)
            out.write(templine)
        else:
            temp_position = target_line + i
            direct = fbuffer[temp_position].split()
            x = direct[0]
            y = direct[1]
            z = direct[2]
            templine = str(x)+' '+str(y)+' '+str(z)+'\n'
            out.write(templine)

def POSCAR_maker_3(input_name, relax_number, output_name):
    rn = relax_number
    f = open(input_name)
    fbuffer = f.readlines()
    atom=[]
    number=[]

    number = fbuffer[6].split()

    total_number = int(number[0])
    
    out = open(output_name, 'w')


    for i in xrange(1,8):
        temp_position = i - 1
        write_line = fbuffer[temp_position].split()
        temp_line = ''
        for j in xrange(len(write_line)):
            if j == len(write_line) -1 :
                temp_line += ' ' + write_line[j] + '\n'
            else:
                temp_line += ' ' + write_line[j]
        out.write(temp_line)

    direct_line = 'Direct'+'\n'
    out.write(direct_line)

    target_line = (int(rn)-1)*(total_number+8) + 7

    for i in xrange(1,total_number+1):
        if i == total_number:
            temp_position = target_line + i
            direct = fbuffer[temp_position].split()
            x = direct[0]
            y = direct[1]
            z = direct[2]
            templine = str(x)+' '+str(y)+' '+str(z)
            out.write(templine)
        else:
            temp_position = target_line + i
            direct = fbuffer[temp_position].split()
            x = direct[0]
            y = direct[1]
            z = direct[2]
            templine = str(x)+' '+str(y)+' '+str(z)+'\n'
            out.write(templine)

def POSCAR_maker_4(input_name, relax_number, output_name):
    rn = relax_number
    f = open(input_name)
    fbuffer = f.readlines()
    atom=[]
    number=[]

    number = fbuffer[6].split()

    total_number = int(number[0])
    
    out = open(output_name, 'w')

    for i in xrange(len(fbuffer)):
        temp = fbuffer[i].split()
        search_1 = temp[0]
        if search_1 == 'Direct':
                search_2 = temp[2]
                if search_2 == relax_number:
                        target_line = i

    for i in xrange(1,8):
        temp_position = target_line - 7 + i - 1
        write_line = fbuffer[temp_position].split()
        temp_line = ''
        for j in xrange(len(write_line)):
            if j == len(write_line) -1 :
                temp_line += ' ' + write_line[j] + '\n'
            else:
                temp_line += ' ' + write_line[j]
        out.write(temp_line)

    direct_line = 'Direct'+'\n'
    out.write(direct_line)


    for i in xrange(1,total_number+1):
        if i == total_number:
            temp_position = target_line + i
            direct = fbuffer[temp_position].split()
            x = direct[0]
            y = direct[1]
            z = direct[2]
            templine = str(x)+' '+str(y)+' '+str(z)
            out.write(templine)
        else:
            temp_position = target_line + i
            direct = fbuffer[temp_position].split()
            x = direct[0]
            y = direct[1]
            z = direct[2]
            templine = str(x)+' '+str(y)+' '+str(z)+'\n'
            out.write(templine)

def anal_npt_MD(input_name, output_name):
    f = open(input_name)

    fbuffer = f.readlines()
    atom=[]
    number=[]

    number = fbuffer[6].split()

    total_number = int(number[0])
    
    out = open(output_name, 'w')
    first_line = 'step'+'\t'+'a'+'\t'+'b'+'\t'+'c'+'\t'+'alpha(deg)'+'\t'+'beta(deg)'+'\t'+'gamma(deg)'+'\n'
    out.write(first_line)
    target_line = []
    index = 0
    step = 0

    for i in xrange(len(fbuffer)):
        temp = fbuffer[i].split()
        search_1 = temp[0]
        if search_1 == 'Direct':
            target_line.append(i)
            index += 1
    #print target_line
    print 'total relaxation is '+str(index)

    for i in target_line:
        step += 1
        latt_a = np.array([float(a) for a in fbuffer[i-5].split()])
        latt_b = np.array([float(b) for b in fbuffer[i-4].split()])
        latt_c = np.array([float(c) for c in fbuffer[i-3].split()])
        a = np.linalg.norm(latt_a)
        b = np.linalg.norm(latt_b)
        c = np.linalg.norm(latt_c)
        alpha = degrees(acos((np.inner(latt_b, latt_c))/(b*c)))
        beta = degrees(acos((np.inner(latt_a, latt_c))/(a*c)))
        gamma = degrees(acos((np.inner(latt_a, latt_b))/(a*b)))
        line = str(step)+'\t'+str(a)+'\t'+str(b)+'\t'+str(c)+'\t'+str(alpha)+'\t'+str(beta)+'\t'+str(gamma)+'\n'
        out.write(line)
    out.close()

            





#anal_npt_MD('XDATCAR', 'XDAT_results')
for i in range(1):
    number=13001 + i*1000
    name='CONTCAR_'+str(number)
    POSCAR_maker_3('XDATCAR',number,name)
#POSCAR_maker_3('XDATCAR',7000,'CONTCAR_7000')