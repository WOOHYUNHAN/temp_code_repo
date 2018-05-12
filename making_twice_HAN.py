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

def twice_XDATCAR_maker(input_name, output_name,a_mul,b_mul,c_mul,sr,er):
        f = open(input_name)
        fbuffer = f.readlines()
        atom=[]
        number=[]
        am = int(a_mul)
        bm = int(b_mul)
        cm = int(c_mul)
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
        a_axis = np.array(aaxis)*am
        b_axis = np.array(baxis)*bm
        c_axis = np.array(caxis)*cm
        number = fbuffer[6].split()
        uni = fbuffer[1].split()

        total_number = int(number[0])
        after_total_number = total_number * (am) * (bm) * (cm)

        out = open(output_name, 'w')

        direct_line_1 = 'gamma_B'+'\n'
        out.write(direct_line_1)

        uni_number = float(uni[0])
        direct_line_2 = str(uni_number)+'\n'
        out.write(direct_line_2)
        line_a = str(a_axis[0])+'\t'+str(a_axis[1])+'\t'+str(a_axis[2])+'\t'+'\n'
        out.write(line_a)
        line_b = str(b_axis[0])+'\t'+str(b_axis[1])+'\t'+str(b_axis[2])+'\t'+'\n'
        out.write(line_b)
        line_c = str(c_axis[0])+'\t'+str(c_axis[1])+'\t'+str(c_axis[2])+'\t'+'\n'
        out.write(line_c)

        direct_line_3 = 'B'+'\n'
        out.write(direct_line_3)

        direct_line_4 = str(after_total_number)+'\n'
        out.write(direct_line_4)

        for i in range(sr,er+1):
                print 'relax = ' +str(i)
                out.write('\n')
                for j in range(1, total_number+1):
                        temp_position = (int(i)-1)*(total_number+1) + 7 + j
                        direct = fbuffer[temp_position].split()
                        direct_m = []
                        direct_m.append(float(direct[0])/float(am))
                        direct_m.append(float(direct[1])/float(bm))
                        direct_m.append(float(direct[2])/float(cm))
                        direct_modified = np.array(direct_m)
                        for a in range(am):
                                for b in range(bm):
                                        for c in range(cm):
                                                corr_a = np.array([1,0,0])/float(am)
                                                pre_a = a*corr_a
                                                corr_b = np.array([0,1,0])/float(bm)
                                                pre_b = b*corr_b
                                                corr_c = np.array([0,0,1])/float(cm)
                                                pre_c = c*corr_c
                                                final_list = direct_modified + pre_a + pre_b + pre_c
                                                #final_list = tranform_direct_cartesian(aaxislist, baxislist, caxislist, final_list)
                                                templine = str(final_list[0])+'\t'+str(final_list[1])+'\t'+str(final_list[2])+'\t'+'\n'
                                                out.write(templine)
        
        out.close()
        
        




def twice_POSCAR_maker(input_name, output_name,a_mul,b_mul,c_mul):
    f = open(input_name)
    fbuffer = f.readlines()
    atom=[]
    number=[]
    am = int(a_mul)
    bm = int(b_mul)
    cm = int(c_mul)

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
    a_axis = np.array(aaxis)*am
    b_axis = np.array(baxis)*bm
    c_axis = np.array(caxis)*cm

    print aaxislist
    print baxislist
    print caxislist
    
    number = fbuffer[6].split()
    uni = fbuffer[1].split()

    total_number = int(number[0])
    after_total_number = total_number * (am) * (bm) * (cm)
    
    out = open(output_name, 'w')

    direct_line_1 = 'gamma_B'+'\n'
    out.write(direct_line_1)

    uni_number = float(uni[0])
    direct_line_2 = str(uni_number)+'\n'
    out.write(direct_line_2)

    line_a = str(a_axis[0])+'\t'+str(a_axis[1])+'\t'+str(a_axis[2])+'\t'+'\n'
    out.write(line_a)
    line_b = str(b_axis[0])+'\t'+str(b_axis[1])+'\t'+str(b_axis[2])+'\t'+'\n'
    out.write(line_b)
    line_c = str(c_axis[0])+'\t'+str(c_axis[1])+'\t'+str(c_axis[2])+'\t'+'\n'
    out.write(line_c)
    
    direct_line_3 = 'B'+'\n'
    out.write(direct_line_3)

    direct_line_4 = str(after_total_number)+'\n'
    out.write(direct_line_4)

    direct_line_5 = 'Direct'+'\n'
    out.write(direct_line_5)


    for i in xrange(1,total_number+1):
            temp_position = 7 + i
            direct = fbuffer[temp_position].split()
            direct_m = []
            direct_m.append(float(direct[0])/float(am))
            direct_m.append(float(direct[1])/float(bm))
            direct_m.append(float(direct[2])/float(cm))
            direct_modified = np.array(direct_m)
            x = direct_modified[0]
            y = direct_modified[1]
            z = direct_modified[2]
            templine = str(x)+'\t'+str(y)+'\t'+str(z)+'\t'+'\n'
            out.write(templine)


 
    for i in xrange(1,total_number+1):
            temp_position = 7 + i
            direct = fbuffer[temp_position].split()
            direct_m = []
            direct_m.append(float(direct[0])/float(am))
            direct_m.append(float(direct[1])/float(bm))
            direct_m.append(float(direct[2])/float(cm))
            print direct_m
            #direct_m.append(
            direct_modified = np.array(direct_m)
            #print direct_modified
            for a in xrange(am):
                    for b in xrange(bm):
                            for c in xrange(cm):
                                    if int(a) == 0 and int(b) == 0 and int(c) ==0:
                                        pass
                                    else:
                                        corr_a = np.array([1,0,0])/float(am)
                                        pre_a = a * corr_a
                                    #print pre_a
                                        corr_b = np.array([0,1,0])/float(bm)
                                        pre_b = b * corr_b
                                    #print pre_b
                                        corr_c = np.array([0,0,1])/float(cm)
                                        pre_c = c * corr_c
                                    #print pre_c
                                        final_list = direct_modified + pre_a + pre_b + pre_c
                                    #final = tranform_direct_cartesian(aaxislist, baxislist, caxislist, final_list)
                                        x = final_list[0]
                                        y = final_list[1]
                                        z = final_list[2]
                                        templine = str(x)+'\t'+str(y)+'\t'+str(z)+'\t'+'\n'
                                        out.write(templine)
def multiple_POSCAR_general(input_name, output_name,a_mul,b_mul,c_mul):
    f = open(input_name)
    fbuffer = f.readlines()
    am = int(a_mul)
    bm = int(b_mul)
    cm = int(c_mul)

    aaxis = []
    baxis = []
    caxis = []    
    uni = fbuffer[1].split()
    uni_number = float(uni[0])
    aaxislist = fbuffer[2].split()
    baxislist = fbuffer[3].split()
    caxislist = fbuffer[4].split()

    for i in xrange(3):
            aaxis.append(float(aaxislist[i]))
            baxis.append(float(baxislist[i]))
            caxis.append(float(caxislist[i]))
    #print aaxis
    a_axis = np.array(aaxis)*am
    b_axis = np.array(baxis)*bm
    c_axis = np.array(caxis)*cm

    #print a_axis
    #print baxislist
    #print caxislist

    atomlist = fbuffer[5].split()
    #print atom
    numberlist = fbuffer[6].split()
    #print numberlist
    number = []
    for i in numberlist:
        number.append(int(i))
    #print number
    new_number = np.array(number) * am *bm * cm
    #print new_number

    out = open(output_name, 'w')
    mapfile = open('map.in','w')
    direct_line_1 = 'system'+'\n'
    out.write(direct_line_1)

    uni_number = float(uni[0])
    direct_line_2 = str(uni_number)+'\n'
    out.write(direct_line_2)

    line_a = str(a_axis[0])+'\t'+str(a_axis[1])+'\t'+str(a_axis[2])+'\t'+'\n'
    out.write(line_a)
    line_b = str(b_axis[0])+'\t'+str(b_axis[1])+'\t'+str(b_axis[2])+'\t'+'\n'
    out.write(line_b)
    line_c = str(c_axis[0])+'\t'+str(c_axis[1])+'\t'+str(c_axis[2])+'\t'+'\n'
    out.write(line_c)
    
    direct_line_3 = ''
    for i in atomlist:
        direct_line_3 += str(i)+' '
    direct_line_3 = direct_line_3 + '\n'
    out.write(direct_line_3)

    direct_line_4 = ''
    for i in new_number:
        direct_line_4 += str(i)+' '
    direct_line_4 = direct_line_4 + '\n'
    out.write(direct_line_4)

    direct_line_5 = 'Direct'+'\n'
    out.write(direct_line_5)

    starting_number = 0
    atomcount = 0
    for i, n in enumerate(atomlist):
        for j in range(1,int(number[i])+1):
            temp_position = 7 + starting_number + j
            #print temp_position
            direct = fbuffer[temp_position].split()
            direct_m = np.array([float(direct[s]) for s in range(3)], dtype=float)
            remaining_part = np.array([str(direct[s]) for s in range(3,len(direct))])
            direct_modified = direct_m / np.array([am,bm,cm])
            for a in range(am):
                for b in range(bm):
                    for c in range(cm):
                        atomcount += 1
                        corr = np.array([1.0,1.0,1.0])/np.array([am,bm,cm], dtype=float)
                        pre = corr * np.array([a,b,c])
                        final = direct_modified + pre
                        final_line = str(final[0])+'\t'+str(final[1])+'\t'+str(final[2])
                        for s in range(len(remaining_part)):
                            final_line += ' ' + str(remaining_part[s]) 
                        final_line += '\n'
                        out.write(final_line)
                        map_line = str(a)+' '+str(b)+' '+str(c)+' '+str(j-1)+' '+str(atomcount)+'\n'
                        mapfile.write(map_line)
                        #print pre
        starting_number += int(number[i])

    out.close()
    mapfile.close()

def replace_atom_position(POSCAR_original, POSCAR_replace, POSCAR_new, atomic_number):
    f = open(POSCAR_original)
    fbuffer = f.readlines()

    g = open(POSCAR_replace)
    gbuffer = g.readlines()

    fnumber = np.array([int(fbuffer[6].split()[i]) for i in range(len(fbuffer[6].split()))])
    gnumber = np.array([int(fbuffer[6].split()[i]) for i in range(len(fbuffer[6].split()))])

    #print np.sum(fnumber)
    if np.sum(fnumber) != np.sum(gnumber):
        print 'ERROR: the number of atoms in POSCAR files does not match'

    f.close()
    g.close()

    out = open(POSCAR_new, 'w')


    for i in range(len(fbuffer)):
        temp = fbuffer[i].split()
        if len(temp) != 0:
            if temp[0] == 'Direct' or temp[0] == 'direct':
                start_f = i
                break

    for i in range(len(gbuffer)):
        temp = gbuffer[i].split()
        if len(temp) != 0:
            if temp[0] == 'Direct' or temp[0] == 'direct':
                start_g = i
                break
    #print start_f, start_g

    for i in range(start_f+1):
        out.write(fbuffer[i])

    replace = 0

    for i in range(1,np.sum(fnumber)+1):
        if i in atomic_number:
            replace += 1
            out.write(gbuffer[start_g+i])
        else:
            out.write(fbuffer[start_f+i])


    out.close()
    print 'The number of replacements = ' + str(replace)



def test_direct_coordi(test):
        if -0.5 <= test <= 0.5:
                final = test
        elif test > 0.5:
                final = test - 1
        elif test < -0.5:
                final = test + 1
        return final

def calculate_distance(position_list, aaxislist, baxislist, caxislist, target_atom):
    relative_position = position_list - position_list[target_atom-1]
    distance = []
    for i in range(len(relative_position)):
        tempx= test_direct_coordi(float(relative_position[i][0]))
        tempy= test_direct_coordi(float(relative_position[i][1]))
        tempz= test_direct_coordi(float(relative_position[i][2]))
        temp = [tempx,tempy,tempz]
        cart = np.array(tranform_direct_cartesian(aaxislist, baxislist, caxislist, temp))
        dis = np.linalg.norm(cart)
        distance.append(dis)

    return np.array(distance)


def make_defective_sample(POSCAR_original, output_name, a_mul,b_mul,c_mul, vn):
    multiple_POSCAR_general(POSCAR_original, 'POSCAR_multipled',a_mul,b_mul,c_mul)
    f = open('POSCAR_multipled')
    fbuffer = f.readlines()
    f.close()
    numberlist = fbuffer[6].split()
    #print numberlist
    number = []
    for i in numberlist:
        number.append(int(i))
    tn = np.sum(np.array(number))


    aaxis = []
    baxis = []
    caxis = []    
    uni = fbuffer[1].split()
    uni_number = float(uni[0])
    aaxislist = fbuffer[2].split()
    baxislist = fbuffer[3].split()
    caxislist = fbuffer[4].split()

    for i in xrange(3):
            aaxis.append(float(aaxislist[i])*uni_number)
            baxis.append(float(baxislist[i])*uni_number)
            caxis.append(float(caxislist[i])*uni_number)

    position_list = []

    for i in range(tn):
        temp = np.array([float(fbuffer[i+8].split()[j]) for j in range(3)])
        position_list.append(temp)

    position_list = np.array(position_list)
    target_atom = np.random.randint(1,tn+1,size=1)
    print 'target atom is ' + str(target_atom)
    distance = calculate_distance(position_list, aaxis, baxis, caxis, target_atom)
    dis_sort = np.argsort(distance)

    #print dis_sort

    vacancy_list=[int(dis_sort[i]) for i in range(vn)]

    print vacancy_list
    #print distance[vacancy_list[0]], distance[vacancy_list[1]], distance[vacancy_list[2]]

    g = open(output_name,'w')

    for i in range(6):
        g.write(fbuffer[i])




    #sw = True
    #while sw is True:
    #    n = np.random.randint(1,tn+1,size=vn)
    #    if len(n) != len(np.unique(n)):
    #        sw = True
    #    else:
    #        sw = False
    #
    #print n

    g.write(str(tn-vn)+'\n'+'Direct'+'\n')

    for i in range(tn):
        if i in vacancy_list:
            pass
        else:
            g.write(fbuffer[i+8])

    g.close()




def shift_POSCAR_general(input_name,shift_atom_range,x_num,y_num,z_num):
    f = open(input_name)
    fbuffer = f.readlines()
    aaxis = []
    baxis = []
    caxis = []    
    uni = fbuffer[1].split()
    uni_number = float(uni[0])
    aaxislist = fbuffer[2].split()
    baxislist = fbuffer[3].split()
    caxislist = fbuffer[4].split()

    for i in xrange(3):
            aaxis.append(float(aaxislist[i]))
            baxis.append(float(baxislist[i]))
            caxis.append(float(caxislist[i]))
    #print aaxis
    a_axis = np.array(aaxis)*uni_number
    b_axis = np.array(baxis)*uni_number
    c_axis = np.array(caxis)*uni_number


    atomlist = fbuffer[5].split()
    #print atom
    numberlist = fbuffer[6].split()
    #print numberlist
    number = []
    for i in numberlist:
        number.append(int(i))
    #print number
    new_number = np.array(number) 

    name_count = 0

    for i in range(x_num):
        for j in range(y_num):
            for k in range(z_num):
                name_count += 1
                name = 'POSCAR-'+ '%03d' %(name_count) 
                #print name
                shift_x = 0 + (0.5/x_num)*i
                shift_y = 0 + (0.5/y_num)*j
                shift_z = 0 + (0.5/z_num)*k
                shift = np.array([shift_x,shift_y,shift_z])
                out = open(name, 'w')
                direct_line_1 = 'system'+'\n'
                out.write(direct_line_1)

                uni_number = float(uni[0])
                direct_line_2 = '1.0000000'+'\n'
                out.write(direct_line_2)

                line_a = str(a_axis[0])+'\t'+str(a_axis[1])+'\t'+str(a_axis[2])+'\t'+'\n'
                out.write(line_a)
                line_b = str(b_axis[0])+'\t'+str(b_axis[1])+'\t'+str(b_axis[2])+'\t'+'\n'
                out.write(line_b)
                line_c = str(c_axis[0])+'\t'+str(c_axis[1])+'\t'+str(c_axis[2])+'\t'+'\n'
                out.write(line_c)

                direct_line_3 = ''
                for r in atomlist:
                    direct_line_3 += str(r)+' '
                direct_line_3 = direct_line_3 + '\n'
                out.write(direct_line_3)

                direct_line_4 = ''
                for r in new_number:
                    direct_line_4 += str(r)+' '
                direct_line_4 = direct_line_4 + '\n'
                out.write(direct_line_4)

                direct_line_5 = 'Direct'+'\n'
                out.write(direct_line_5)

                for x in range(np.sum(new_number)):
                    temp_position = 8 + x
                    direct = fbuffer[temp_position].split()
                    direct_m = np.array([direct_coordi for direct_coordi in direct], dtype=float)
                    if x+1 in shift_atom_range:
                        direct_m = direct_m + shift
                    final_line = str(direct_m[0])+'\t'+str(direct_m[1])+'\t'+str(direct_m[2])+'\t'+'\n'
                    out.write(final_line)           
                out.close()
    


#twice_POSCAR_maker('CONTCAR_4','POSCAR5',1,1,5)
#twice_POSCAR_maker('CONTCAR_4','POSCAR6',1,1,6)
#twice_POSCAR_maker('CONTCAR_4','POSCAR7',1,1,7)
#twice_POSCAR_maker('CONTCAR_4','POSCAR8',1,1,8)
#twice_POSCAR_maker('POSCAR','POSCAR_1',2,2,2)
#multiple_POSCAR_general('CONTCAR', 'POSCAR_Ag_sub',6,6,1)
multiple_POSCAR_general('POSCAR', 'CONTCAR', 3,1,1)
#twice_XDATCAR_maker('XDATCAR_5000_15000', 'XDATCAR_5000_15000_2x2x2',2,2,2,1,10000)
#replace_atom_position('POSCAR_1', 'POSCAR_2','POSCAR_new', [46,38,4,36,26,14,28,42,56,52,40 ,30   ,39,3,35,25,13,27,41,55,51,45,37,29])
#make_defective_sample('CONTCAR','POSCAR',2,2,2,3)
#for i in range(5):
#    name = 'POSCAR_defec'+str(i)
#    make_defective_sample('CONTCAR',name,2,1,1,2)
