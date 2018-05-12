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

def basic_info(input_name, output_name):
    f = open(input_name)
    fbuffer = f.readlines()
    atom=[]
    number=[]

    aaxislist = fbuffer[2].split()
    baxislist = fbuffer[3].split()
    caxislist = fbuffer[4].split()

    atom = fbuffer[5].split()
    number = fbuffer[6].split()

    total_number = int(number[0])

    out = open(output_name, 'w')
    firstline = str(total_number)+'\n'
    secondline = '\n'

    total_relax = (len(fbuffer) - 7 )/(total_number+1)

    for i in xrange(1, total_relax+1):
            target_line = (int(i)-1)*(total_number+1) + 7
            out.write(firstline)
            out.write(secondline)
            print 'subprocess is done---------------'+str(i)
            for j in xrange(1, total_number+1):
                    temp_position = target_line + j
                    direct = fbuffer[temp_position].split()
                    temp = tranform_direct_cartesian(aaxislist, baxislist, caxislist, direct)
                    x = temp[0]
                    y = temp[1]
                    z = temp[2]
                    templine = str(atom[0])+'\t'+str(x)+'\t'+str(y)+'\t'+str(z)+'\t'+'\n'
                    out.write(templine)

def basic_info_2(input_name, output_name):
    f = open(input_name)
    fbuffer = f.readlines()
    atom=[]
    number=[]

    aaxislist = fbuffer[2].split()
    baxislist = fbuffer[3].split()
    caxislist = fbuffer[4].split()

    atom = fbuffer[5].split()
    number = fbuffer[6].split()

    total_number = int(number[0])

    out = open(output_name, 'w')
    firstline = str(total_number)+'\n'
    secondline = '\n'

    #total_relax = (len(fbuffer) - 7 )/(total_number+1)

    for i in xrange(1, 4400):
            target_line = (int(i)-1)*(total_number+8) + 7
            out.write(firstline)
            out.write(secondline)
            print 'subprocess is done---------------'+str(i)
            for j in xrange(1, total_number+1):
                    temp_position = target_line + j
                    direct = fbuffer[temp_position].split()
                    temp = tranform_direct_cartesian(aaxislist, baxislist, caxislist, direct)
                    x = temp[0]
                    y = temp[1]
                    z = temp[2]
                    templine = str(atom[0])+'\t'+str(x)+'\t'+str(y)+'\t'+str(z)+'\t'+'\n'
                    out.write(templine)


def multiple_info(input_name, output_name):
    f = open(input_name)
    fbuffer = f.readlines()
    atom=[]
    number=[]

    aaxis = []
    baxis = []
    caxis = []

    aaxislist = fbuffer[2].split()
    baxislist = fbuffer[3].split()
    caxislist = fbuffer[4].split()
    for i in xrange(3):
            aaxis.append(float(aaxislist[i]))
            baxis.append(float(baxislist[i]))
            caxis.append(float(caxislist[i]))
    print aaxis
    print baxis
    print caxis
    atom = fbuffer[5].split()
    number = fbuffer[6].split()

    total_number = int(number[0])
    multiple_number = total_number * 8

    out = open(output_name, 'w')
    firstline = str(multiple_number)+'\n'
    secondline = '\n'

    total_relax = (len(fbuffer) - 7 )/(total_number+1)

    for i in xrange(1, 10000):
            target_line = (int(i)-1)*(total_number+1) + 7
            out.write(firstline)
            out.write(secondline)
            print 'subprocess is done---------------'+str(i)
            for j in xrange(1, total_number+1):
                    temp_position = target_line + j
                    direct = fbuffer[temp_position].split()
                    direct_m = []
                    for k in xrange(len(direct)):
                            direct_m.append(float(direct[k]))
                    direct_modified = np.array(direct_m)
                    #print direct
                    #print direct_modified
                    for a in xrange(2):
                            new_a = a
                            for b in xrange(2):
                                    new_b = b
                                    for c in xrange(2):
                                            new_c = c
                                            pre_a = new_a * np.array([1,0,0])
                                            pre_b = new_b * np.array([0,1,0])
                                            pre_c = new_c * np.array([0,0,1])
                                            #temp = tranform_direct_cartesian(aaxislist, baxislist, caxislist, direct_modified)
                                            final_list = direct_modified + pre_a + pre_b + pre_c                                            
                                            final_list = tranform_direct_cartesian(aaxislist, baxislist, caxislist, final_list)
                                            x = final_list[0]
                                            y = final_list[1]
                                            z = final_list[2]
                                            templine = str(atom[0])+'\t'+str(x)+'\t'+str(y)+'\t'+str(z)+'\t'+'\n'
                                            #templine = str(x)+'\t'+str(y)+'\t'+str(z)+'\t'+'\n'
                                            out.write(templine)


def multiple_info_2(input_name, output_name,size, start,finish):
    f = open(input_name)
    fbuffer = f.readlines()
    atom=[]
    number=[]

    nx = size[0]
    ny = size[1]
    nz = size[2]

    aaxis = []
    baxis = []
    caxis = []

    aaxislist = fbuffer[2].split()
    baxislist = fbuffer[3].split()
    caxislist = fbuffer[4].split()
    for i in xrange(3):
            aaxis.append(float(aaxislist[i]))
            baxis.append(float(baxislist[i]))
            caxis.append(float(caxislist[i]))
    print aaxis
    print baxis
    print caxis
    atom = fbuffer[5].split()
    number = fbuffer[6].split()

    total_number = int(number[0])
    multiple_number = total_number * nx * ny * nz

    out = open(output_name, 'w')
    firstline = str(multiple_number)+'\n'
    secondline = '\n'

    total_relax = (len(fbuffer) - 7 )/(total_number+1)

    for i in xrange(start, finish+1):
            target_line = (int(i)-1)*(total_number+8) + 7
            out.write(firstline)
            out.write(secondline)
            print 'subprocess is done---------------'+str(i)
            for j in xrange(1, total_number+1):
                    temp_position = target_line + j
                    direct = fbuffer[temp_position].split()
                    direct_m = []
                    for k in xrange(len(direct)):
                            direct_m.append(float(direct[k]))
                    direct_modified = np.array(direct_m)
                    #print direct
                    #print direct_modified
                    for a in xrange(nx):
                            new_a = a
                            for b in xrange(ny):
                                    new_b = b
                                    for c in xrange(nz):
                                            new_c = c
                                            pre_a = new_a * np.array([1,0,0])
                                            pre_b = new_b * np.array([0,1,0])
                                            pre_c = new_c * np.array([0,0,1])
                                            #temp = tranform_direct_cartesian(aaxislist, baxislist, caxislist, direct_modified)
                                            final_list = direct_modified + pre_a + pre_b + pre_c                                            
                                            final_list = tranform_direct_cartesian(aaxislist, baxislist, caxislist, final_list)
                                            x = final_list[0]
                                            y = final_list[1]
                                            z = final_list[2]
                                            templine = str(atom[0])+'\t'+str(x)+'\t'+str(y)+'\t'+str(z)+'\t'+'\n'
                                            #templine = str(x)+'\t'+str(y)+'\t'+str(z)+'\t'+'\n'
                                            out.write(templine)




def relax_info(input_name, output_name):
    f = open(input_name)
    fbuffer = f.readlines()
    f.close()
    atom=[]
    number=[]

    aaxis = []
    baxis = []
    caxis = []

    atom = fbuffer[5].split()
    number = fbuffer[6].split()

    total_number = int(number[0])
    multiple_number = total_number * 8

    
    total_relax = (len(fbuffer) - 7 )/(total_number+1)

    print total_relax
    
def xyz_from_OUTCAR(input_name, nx, ny, nz, output_name):
    f = open(input_name)
    fbuffer = f.readlines()
    f.close()
    for i in range(len(fbuffer)):
        temp = fbuffer[i].split()
        if temp == []:
            pass
        elif temp[0] == 'Dimension': 
            temp2 = fbuffer[i+2].split()
            #print temp2
            NIONS = int(temp2[11])
            #print NIONS
            break

    lattice_array = []
    position_array = []
    for i in range(len(fbuffer)):
        temp = fbuffer[i].split()
        if temp == []:
            pass
        elif temp[0] == 'POSITION': 
            #print i
            lattice_array.append(i-44)
            position_array.append(i)
            #break

    #print lattice_array
    out = open(output_name, 'w')

    #for i in range(301323,301324):
    for i in lattice_array:
        atomic_position = []
        #print i
        aaxislist = np.array([float((fbuffer[i+1].split())[j]) for j in range(3)])
        baxislist = np.array([float((fbuffer[i+2].split())[j]) for j in range(3)])
        caxislist = np.array([float((fbuffer[i+3].split())[j]) for j in range(3)])
        baseline1 = str(NIONS*nx*ny*nz)+'\n'
        baseline2 = '\n'        
        out.write(baseline1)
        out.write(baseline2)
        for k in range(NIONS):
            position = np.array([float((fbuffer[i+44+2+k].split())[j]) for j in range(3)])
            for a in range(nx):
                for b in range(ny):
                    for c in range(nz):
                        add_list = a * np.array([1,0,0]) + b * np.array([0,1,0]) + c * np.array([0,0,1])
                        add = tranform_direct_cartesian(aaxislist, baxislist, caxislist, add_list)
                        nposition = position + add
                        atomic_position.append(position)
                        templine = 'B' + '\t' + str(nposition[0]) + '\t' + str(nposition[1]) + '\t' + str(nposition[2]) + '\n'
                        out.write(templine)
        #print len(atomic_position)


        #print atomic_position

    out.close()




#xyz_from_OUTCAR('OUTCAR1', 2, 2, 2, 'output_name.xyz')
multiple_info_2('XDATCAR','alpha.xyz',[5,5,5],1,1)
