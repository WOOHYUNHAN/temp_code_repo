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



def tranform_direct_cartesian(aaxislist, baxislist, caxislist, atomlist1):
        a = np.array(aaxislist)
        b = np.array(baxislist)
        c = np.array(caxislist)

        cart = (float(atomlist1[0]) * a) + (float(atomlist1[1]) * b) + (float(atomlist1[2]) * c)

        return cart


def data_structure(input_name):
        f= open(input_name)
        fbuffer = f.readlines()

        f.close()

        for i in range(10):
            if str(fbuffer[i].split()[0]) == 'Direct':
                start = i
                break
        #print start
        uni = float(fbuffer[1].split()[0])

        atom = np.array([str(fbuffer[5].split()[i]) for i in range(len(fbuffer[5].split()))])
        num_atom = np.array([int(fbuffer[6].split()[i]) for i in range(len(atom))])

        aaxislist = np.array([float(fbuffer[2].split()[i]) for i in range(3)]) * uni
        baxislist = np.array([float(fbuffer[3].split()[i]) for i in range(3)]) * uni
        caxislist = np.array([float(fbuffer[4].split()[i]) for i in range(3)]) * uni

        lattice = np.array([aaxislist,baxislist,caxislist])


        position = []

        for element in range(1, np.sum(num_atom)+1):
            #print element
            pos1 = np.array([float(fbuffer[start + element].split()[i]) for i in range(3)])
            position.append(pos1)


        return lattice, np.array(position), atom, num_atom





def DOS(input_data_structure, input_name, PDOS_list, option=0):
    ###
    # option:Ef
    ###
    lattice, position, atom, num_atom = data_structure(input_data_structure)

    f= open(input_name)
    fbuffer = f.readlines()
    f.close()

    grid = int(fbuffer[5].split()[2]) 
    Ef = float(fbuffer[5].split()[3])

    if option == 0:
        corr = Ef # Ef
        VBM_position = 0
    elif option ==1:
        final_list, VT_position, VBM_position, CBM_position, orig_bg, corr_bg = LIPR('CONTCAR', 'PROCAR', -0.7, Ef)
        corr = Ef + VBM_position  # VBM
    else:
        print 'Option ERROR'

    energy = np.array([float(fbuffer[6+i].split()[0]) for i in range(grid)]) - corr

    #print energy

    DOS_list = []
    ortibal_DOS_list = []

    #### DOS list = [[[1_s],[1_p],[1_d]],[[2_s],[2_p],[2_d]],[[3_s],[3_p],[3_d]]...]

    for i in range(np.sum(num_atom)):
        start_line = 5 + (i+1) * (grid + 1)
        s = []
        p = []
        d = []
        p1 = []
        p2 = []
        p3 = []
        d1 = []
        d2 = []
        d3 = []
        d4 = []
        d5 = []
        #print start_line
        for j in range(grid):
            #print len(fbuffer[start_line+j+1].split())
            x =  np.array([float(fbuffer[start_line+j+1].split()[k]) for k in range(len(fbuffer[start_line+j+1].split()))])
            #print x[5:10]
            s.append(np.sum(x[1:2]))
            p.append(np.sum(x[2:5]))
            d.append(np.sum(x[5:10]))
            p1.append(np.sum(x[2:3]))
            p2.append(np.sum(x[3:4]))
            p3.append(np.sum(x[4:5]))
            d1.append(np.sum(x[5:6]))
            d2.append(np.sum(x[6:7]))
            d3.append(np.sum(x[7:8]))
            d4.append(np.sum(x[8:9]))
            d5.append(np.sum(x[9:10]))
        temp = np.array([s,p,d])
        temp2 = np.array([s,p1,p2,p3,d1,d2,d3,d4,d5])
        DOS_list.append(temp)
        ortibal_DOS_list.append(temp2)


    DOS_list = np.array(DOS_list)
    ortibal_DOS_list = np.array(ortibal_DOS_list)
    #print DOS_list[0][2][23]

    temp = []
    name = []

    temp.append(energy)

    for i in PDOS_list:
        file_name = 'PDOS_' + str(i) +'.out'
        first_line = 'Energy' + '\t' + str(i)+'_s' + '\t' + str(i)+'_p' + '\t' + str(i)+'_d'+'\n'
        g = open(file_name,'w')
        g.write(first_line)
        for j in range(grid):
            temp_line = str(energy[j])
            #for k in range(3): #s p d
            #    temp_line += '\t' + str(DOS_list[i-1][k][j])
            for k in range(9): # s px py pz d1 d2 d3 d4 d5
                temp_line += '\t' + str(ortibal_DOS_list[i-1][k][j])
            temp_line += '\n'
            g.write(temp_line)
        g.close()
        # array for output
        # temp = [[energy, i_s, i_p, i_d, name],[.....]..]
        temp.append(DOS_list[i-1][0])
        temp.append(DOS_list[i-1][1])
        temp.append(DOS_list[i-1][2])
        name.append('O_'+str(i)+'s')
        name.append('O_'+str(i)+'p')
        name.append('O_'+str(i)+'d')


    temp.append(name)
    projected_out_dos = temp
    #print projected_out_dos[7]
    #print len(projected_out_dos)

    
    #print len(num_atom)
    final_list = []
    temp = []
    temp.append(energy)
    name = []
    #### final list = [[[In_s],[In_p],[In_d]],[[Ga_s],[Ga_p],[Ga_d]],[[Zn_s],[Zn_p],[Zn_d]]...]
    for i in range(len(num_atom)):
        o_1, o_2 = np.sum(num_atom[:i]), np.sum(num_atom[:i+1])
        #print o_1, o_2
        s, p, d = np.zeros(grid), np.zeros(grid), np.zeros(grid)
        for j in range(o_1, o_2):
            #print j
            s += np.array(DOS_list[j][0])
            p += np.array(DOS_list[j][1])
            d += np.array(DOS_list[j][2])
        #print d
        temp2 = np.array([s,p,d])
        final_list.append(temp2)
        temp.append(s)
        temp.append(p)
        temp.append(d)
        name.append(str(atom[i])+'_s')
        name.append(str(atom[i])+'_p')
        name.append(str(atom[i])+'_d')
    final_list = np.array(final_list)
    #print name
    temp.append(name)

    total_out_dos = temp
    # total_out_dos = [energy, In_s ....., O_d, name]
    #print total_out_dos[13]
    #print len(total_out_dos)

    #print name

    out_name = 'DOS_total.out'
    g = open(out_name, 'w')
    first_line = 'Energy'
    for i in range(len(name)):
        first_line += '\t' + str(name[i])
    first_line += '\t' + 'total' + '\n'
    g.write(first_line)

    for i in range(grid):
        total = np.sum(np.array([float(final_list[x][y][i]) for x in range(len(num_atom)) for y in range(3)]))
        temp_line = str(energy[i]) 
        for j in range(len(atom)):
            for k in range(3):
                temp_line += '\t' + str(final_list[j][k][i])
        temp_line += '\t' + str(total) + '\n'
        g.write(temp_line)
    g.close()

    return Ef, VBM_position, total_out_dos, projected_out_dos

def calculate_carrierdensity_fromQE(file_name, ref_energy, chemical_pot, dim, Area, temperature):
    ref=float(ref_energy)
    cp=float(chemical_pot)
    k = 8.6173303/1e5 # Kb unit: eV/K
    kT = float(k*temperature)

    if ref > cp:
        #print "Do you want to calculate hole carrier density?"
        opt = 1
    else:
        #print "Do you want to calculate electron carrier density?"
        opt = 2



    f=open(file_name,'r')
    tempf = f.readlines()
    f.close()

    cell_vol = float(tempf[3].split()[4]) # unit = ang3




    energy = [] # unit = eV
    DOS = [] # unit = states/eV/unitcell

    for i in range(5,len(tempf)):
        tt = tempf[i].split()
        energy.append(float(tt[0]))
        DOS.append(float(tt[1]))

    dE = float(energy[1]) - float(energy[0])

    cal = []

    for i in range(len(energy)):
        if opt ==1:
            #print 'calculate hole carrier density'
            if float(energy[i]) <= ref_energy:
                if (float(energy[i])-cp)/kT > 700:
                    FD = 0
                else:
                    FD = 1.0 / (1.0 + exp((float(energy[i])-cp)/kT))
                #FD = 0
                cal.append(-float(DOS[i]) * (1.0-FD))
        elif opt ==2:
            #print 'calculate electron carrier density'
            if ref_energy <= float(energy[i]):
                if (float(energy[i])-cp)/kT > 700:
                    FD = 0
                else:
                    FD = 1.0 / (1.0 + exp((float(energy[i])-cp)/kT))
                #FD = 1
                cal.append(float(DOS[i]) * FD)
        else:
            print 'optin error'
    cal = np.array(cal) * dE
    density = np.sum(cal)
    #print density

    if dim == 2:
        #print '2D case'
        unit = Area
        return (density/unit) * 1e16
        print (density/unit) * 1e16 # unit = cm-2
    elif dim == 3:
        #print '3D case'
        unit = cell_vol
        return (density/unit) * 1e24
        print (density/unit) * 1e24 # unit = cm-3

    #print cell_vol

def cal_mobility_from_boltzwann(file_name1, file_name2, output_filename, ref_energy, area, height, temperature):
    '''
    conductivity = e * density * mobiilty
    '''
    f = open(file_name1,'r')
    tempf = f.readlines()
    f.close()

    e = 1.6 * 1e-19

    energy = []
    mobility_xx =[]
    mobility_xy =[]
    mobility_yy =[]
    mobility_xz =[]
    mobility_yz =[]
    mobility_zz =[]
    den = []

    for i in range(3, len(tempf)):
        tt = tempf[i].split()
        energy.append(float(tt[0])-ref_energy)
        density = calculate_carrierdensity_fromQE(file_name2, ref_energy, float(tt[0]), 2, area, temperature)
        den.append(float(density))
        if density != 0.0:
            mobility_xx.append((float(tt[2])*0.01*height*1e-8)/density/e)
            mobility_xy.append((float(tt[3])*0.01*height*1e-8)/density/e)
            mobility_yy.append((float(tt[4])*0.01*height*1e-8)/density/e)
            mobility_xz.append((float(tt[5])*0.01*height*1e-8)/density/e)
            mobility_yz.append((float(tt[6])*0.01*height*1e-8)/density/e)
            mobility_zz.append((float(tt[7])*0.01*height*1e-8)/density/e)
        else:
            mobility_xx.append(0)
            mobility_xy.append(0)
            mobility_yy.append(0)
            mobility_xz.append(0)
            mobility_yz.append(0)
            mobility_zz.append(0)

    g = open(output_filename,'w')

    for i in range(len(energy)):
        templine = str(energy[i]) + '\t' + str(den[i]) + '\t' + str(mobility_xx[i]) + '\t' + str(mobility_xy[i]) + '\t' + str(mobility_yy[i]) + '\t' + str(mobility_xz[i]) + '\t' + str(mobility_yz[i]) + '\t' + str(mobility_zz[i]) + '\n'
        g.write(templine)
    g.close()

   # for i in range(len(energy)):
   #     density = calculate_carrierdensity_fromQE(file_name2, Ef, float(energy[i]), area)





DOS('CONTCAR', 'DOSCAR', np.array([1]), option=0)
#calculate_carrierdensity_fromQE('P_boltzdos.dat', -1.6129, -2.6129, 2, 14.990499426689647539,  300)
#cal_mobility_from_boltzwann('P_elcond.dat', 'P_boltzdos.dat', 'mobility.dat', -1.3414,  23.515859788966510861167104424881, 5.1, 300) #green
#cal_mobility_from_boltzwann('P_elcond.dat', 'P_boltzdos.dat', 'mobility.dat', -1.222,  15.200478194132299509, 5.21, 300) #black
#cal_mobility_from_boltzwann('P_elcond.dat', 'P_boltzdos.dat', 'mobility.dat', -2.623295,  9.3057456915905426882280120614882, 4.8, 500) #blue
#cal_mobility_from_boltzwann('P_elcond.dat', 'P_boltzdos.dat', 'mobility.dat', -2.21,  12.96843206215056107600036989585, 0.0302*2, 300) #silicene