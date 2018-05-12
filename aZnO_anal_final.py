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

def get_Ef(input_name):
    f= open(input_name)
    fbuffer = f.readlines()
    f.close()

    grid = int(fbuffer[5].split()[2]) 
    Ef = float(fbuffer[5].split()[3])

    return Ef
        
def count_bonding(lattice, position, element1, element2, r_cutoff, multi):
        aaxislist = np.array(lattice[0])
        baxislist = np.array(lattice[1])
        caxislist = np.array(lattice[2])

        pos1 = position[element1-1]
        pos2 = position[element2-1]

        volume = np.dot(aaxislist, np.cross(baxislist, caxislist))


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
        count = 0
        for i in range(len(multi_atoms)):
                diff = pos1 - multi_atoms[i]
                if np.linalg.norm(diff) == 0:
                        continue
                diff = diff - np.round(diff/multi)*multi
                diff = tranform_direct_cartesian(aaxislist, baxislist, caxislist, diff)
                diff = np.linalg.norm(diff)

                if diff <= r_cutoff:
                    count += 1
        return count



def get_rdf(lattice, position, element1, element2, r_max, bin_spacing, multi):
        aaxislist = np.array(lattice[0])
        baxislist = np.array(lattice[1])
        caxislist = np.array(lattice[2])

        pos1 = position[element1-1]
        pos2 = position[element2-1]

        volume = np.dot(aaxislist, np.cross(baxislist, caxislist))


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
        #hist, bin_edges = np.histogram(r, bins = int(r_max / bin_spacing), range = (0, r_max), density = True)
        #rdf = hist / 4. / np.pi / (bin_edges + bin_spacing / 2.)[:-1] ** 2 * volume / bin_spacing
        if (np.array(r) < r_max).any() :
                hist, bin_edges = np.histogram(r, bins = int(r_max / bin_spacing), range = (0, r_max), density = True)
                rdf = hist / 4. / np.pi / (bin_edges + bin_spacing / 2.)[:-1] ** 2 * volume / bin_spacing
                return bin_edges[:-1], rdf
        else :
                bin = np.linspace(0, r_max, num = r_max / bin_spacing , endpoint = False)
                rdf = np.zeros((len(bin)))
                return bin, rdf

def get_pair_corr(input_name, r_max, bin_spacing, multi, option1, option2):
        #######
        # option1 and option2 are the name of atom that we want to calculate pcf function
        #######
        lattice, position, atom, num_atom = data_structure(input_name)
        bin = np.linspace(0, r_max, num = r_max / bin_spacing , endpoint = False)
        pair_corr = np.zeros((len(bin)))
        if option1 == 'Zn':
            s = np.where(atom == 'Zn')[0]
        elif option1 == 'O':
            s = np.where(atom == 'O')[0]        
        else:
            print 'option1 is error; you insert wrong atom'

        if option2 == 'Zn':
            o = np.where(atom == 'Zn')[0]
            o_1, o_2 = np.sum(num_atom[:o]) , np.sum(num_atom[:o+1])
        elif option2 == 'O':
            o = np.where(atom == 'O')[0]
            o_1, o_2 = np.sum(num_atom[:o]) , np.sum(num_atom[:o+1])
        elif option2 == 'metal':
            o = np.where(atom == 'O')[0]
            o_1, o_2 = 0, np.sum(num_atom[:o])
        elif option2 == 'All':
            o_1, o_2 = 0, np.sum(num_atom)     
        else:
            print 'option2 is error; you insert wrong atom'
        s_1, s_2 = np.sum(num_atom[:s]) , np.sum(num_atom[:s+1])
        
        print s_1, s_2, o_1, o_2
        for i in range(s_1+1, s_2 + 1):
                for j in range(o_1+1, o_2 + 1):
                        element1=i
                        element2=j
                        r, rdf = get_rdf(lattice, position, element1, element2, r_max, bin_spacing, multi)
                        pair_corr = pair_corr + rdf
        return r, pair_corr

def plot_pair_corr(r, pair_corr):
    plt.plot(r, pair_corr)
    plt.show()
                        
def cos_dist(array1, array2):

    array1_norm = np.linalg.norm(array1)
    array2_norm = np.linalg.norm(array2)

    dist  = 0.5 * (1. - (np.array(array1) * np.array(array2)).sum() / array1_norm / array2_norm)
    return dist

def average_Nc_bc(input_name, r_cutoff, multi, option1, option2):
        #######
        # option1 and option2 are the name of atom that we want to calculate pcf function
        #######
        lattice, position, atom, num_atom = data_structure(input_name)
        if option1 == 'Zn':
            s = np.where(atom == 'Zn')[0]
        elif option1 == 'O':
            s = np.where(atom == 'O')[0]        
        else:
            print 'option1 is error; you insert wrong atom'

        if option2 == 'Zn':
            o = np.where(atom == 'Zn')[0]
            o_1, o_2 = np.sum(num_atom[:o]) , np.sum(num_atom[:o+1])
        elif option2 == 'O':
            o = np.where(atom == 'O')[0]
            o_1, o_2 = np.sum(num_atom[:o]) , np.sum(num_atom[:o+1])
        elif option2 == 'metal':
            o = np.where(atom == 'O')[0]
            o_1, o_2 = 0, np.sum(num_atom[:o])
        elif option2 == 'All':
            o_1, o_2 = 0, np.sum(num_atom)     
        else:
            print 'option2 is error; you insert wrong atom'
        s_1, s_2 = np.sum(num_atom[:s]) , np.sum(num_atom[:s+1])
        
        #print s_1, s_2, o_1, o_2

        final_set = [[],[],[],[],[],[],[]]
        final = 0
        for i in range(s_1+1, s_2 + 1):
            temp = 0
            for j in range(o_1+1, o_2 + 1):
                count = count_bonding(lattice, position, i, j, r_cutoff, multi)
                temp += count
                final += count
            if temp == 0:
                final_set[0].append(i)
            elif temp == 1:
                final_set[1].append(i)
            elif temp == 2:
                final_set[2].append(i)
            elif temp == 3:
                final_set[3].append(i)
            elif temp == 4:
                final_set[4].append(i)
            elif temp == 5:
                final_set[5].append(i)
            elif temp == 6:
                final_set[6].append(i)
            else:
                print 'CRITICAL ERROR: coordination number of ' + str(i) + ' is out of range!!!'
                print 'Its coordination number is ' + str(temp)
        
        
        len_final_set = [len(final_set[0]),len(final_set[1]),len(final_set[2]),len(final_set[3]),len(final_set[4]),len(final_set[5]),len(final_set[6])]
        print final_set
        print str(option1) + ' : ' + str(len_final_set)
        final = final / float(num_atom[s])

        print np.array(len_final_set)

        return final

def Nc_oxygen_with_metal(input_name, r_cutoff_set, multi):
    lattice, position, atom, num_atom = data_structure(input_name)
    #print atom, num_atom

    #find oxygen atoms
    o = np.where(atom == 'O')[0]
    o_1, o_2 = np.sum(num_atom[:o]) , np.sum(num_atom[:o+1])

    final = [[],[],[],[],[],[]] # list for several coorination number s

    for i in range(o_1 +1, o_2 +1):
        count = 0
        for j in range(len(r_cutoff_set)):
            t_1, t_2 = np.sum(num_atom[:j]), np.sum(num_atom[:j+1])
            for k in range(t_1 +1, t_2 +1):
                temp = count_bonding(lattice, position, i, k, float(r_cutoff_set[j]), multi)
                count += temp
        if count == 0:
            final[0].append(i)
        elif count == 1:
            final[1].append(i)
        elif count == 2:
            final[2].append(i)
        elif count == 3:
            final[3].append(i)
        elif count == 4:
            final[4].append(i)
        elif count == 5:
            final[5].append(i)
        else:
            print 'CRITICAL ERROR: coordination number of ' + str(i) + ' is out of range!!!'
  

    #print final
    len_final = [len(final[0]),len(final[1]),len(final[2]),len(final[3]),len(final[4]),len(final[5])]
    print 'O : ' + str(len_final)
    #print 'Nc = 0: ' + str(len(final[0]))
    #print 'Nc = 1: ' + str(len(final[1]))
    #print 'Nc = 2: ' + str(len(final[2]))
    #print 'Nc = 3: ' + str(len(final[3]))
    #print 'Nc = 4: ' + str(len(final[4]))
    #print 'Nc = 5: ' + str(len(final[5]))

    return 0







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

    #### DOS list = [[[1_s],[1_p],[1_d]],[[2_s],[2_p],[2_d]],[[3_s],[3_p],[3_d]]...]

    for i in range(np.sum(num_atom)):
        start_line = 5 + (i+1) * (grid + 1)
        s = []
        p = []
        d = []
        #print start_line
        for j in range(grid):
            #print len(fbuffer[start_line+j+1].split())
            x =  np.array([float(fbuffer[start_line+j+1].split()[k]) for k in range(len(fbuffer[start_line+j+1].split()))])
            #print x[5:10]
            s.append(np.sum(x[1:2]))
            p.append(np.sum(x[2:5]))
            d.append(np.sum(x[5:10]))
        temp = np.array([s,p,d])
        DOS_list.append(temp)


    DOS_list = np.array(DOS_list)
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
            temp_line = str(energy[j]) + '\t' + str(DOS_list[i-1][0][j]) + '\t' + str(DOS_list[i-1][1][j]) + '\t' + str(DOS_list[i-1][2][j])+'\n'
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

def plot_DOS(final_list, plot_energy_range):
    fig = plt.figure()
    nkpt = len(final_list)
    #if option == 0:
    #    corr = 0 # Ef
    #elif option ==1:
    #    corr = 0 # DOS is alreday set to be VBM at zero
    #else:
    #    print 'Option ERROR'

    energy = np.array(final_list[0])

    for i in range(nkpt-2):
        ylist = final_list[i+1]
        name = str(final_list[nkpt-1][i])
        plt.plot(energy,ylist,'-', label=name)

    #handles, labels = fig.get_legend_handles_labels()
    #leg=fig.legend(handles[::-1], labels[::-1], frameon=False ,loc='lower left',numpoints=1,handletextpad=0.8,borderpad=0.1, ncol=2,labelspacing=0.1, handlelength=1, prop={'size':19})
    plt.legend()
    plt.axis(plot_energy_range)
    #plt.plot([plot_energy_range[0],plot_energy_range[1]],[criteria, criteria], '--', color='red')
    #plt.plot([VBM_position, VBM_position],[plot_energy_range[2],plot_energy_range[3]], color='blue')
    #plt.plot([CBM_position, CBM_position],[plot_energy_range[2],plot_energy_range[3]], color='blue')
    plt.plot([0,0], [plot_energy_range[2],plot_energy_range[3]], color='black')
    fig.savefig('totalDOS.png')
    #plt.show()

def plot_PDOS(final_list, plot_energy_range):
    fig = plt.figure()
    nkpt = len(final_list)
    #if option == 0:
    #    corr = 0 # Ef
    #elif option ==1:
    #    corr = 0 # DOS is alreday set to be VBM at zero
    #else:
    #    print 'Option ERROR'

    energy = np.array(final_list[0])

    for i in range(nkpt-2):
        ylist = final_list[i+1]
        name = str(final_list[nkpt-1][i])
        plt.plot(energy,ylist,'-', label=name)

    #handles, labels = fig.get_legend_handles_labels()
    #leg=fig.legend(handles[::-1], labels[::-1], frameon=False ,loc='lower left',numpoints=1,handletextpad=0.8,borderpad=0.1, ncol=2,labelspacing=0.1, handlelength=1, prop={'size':19})
    plt.legend()
    plt.axis(plot_energy_range)
    #plt.plot([plot_energy_range[0],plot_energy_range[1]],[criteria, criteria], '--', color='red')
    #plt.plot([VBM_position, VBM_position],[plot_energy_range[2],plot_energy_range[3]], color='blue')
    #plt.plot([CBM_position, CBM_position],[plot_energy_range[2],plot_energy_range[3]], color='blue')
    plt.plot([0,0], [plot_energy_range[2],plot_energy_range[3]], color='black')
    fig.savefig('PDOS.png')
    #plt.show()


def LIPR(input_data_structure, input_name, criteria, Ef):
    lattice, position, atom, num_atom = data_structure(input_data_structure)
    f= open(input_name)
    fbuffer = f.readlines()
    f.close()

    #print fbuffer[1].split()
    nkpt = int(fbuffer[1].split()[3]) 
    nband = int(fbuffer[1].split()[7])
    nion = int(fbuffer[1].split()[11])

    print nkpt, nband, nion

    len_band = 1 + 1 + 1 + nion + 1 + 1 + nion + nion + 1
    len_kpoint = 3 + len_band * nband
    a =  3 + len_kpoint * (18 - 1) + 2 + len_band * (334 - 1) + 2
    #print fbuffer[a].split()


    print nion


    final_list = []

    for i in range(nkpt):
        print 'Now k-point is...' + str(i+1) + '/' + str(nkpt)
        band_num = []
        energy = []
        LIPR_list = []
        for j in range(nband):
            line = 3 + len_kpoint * (i) + 2 + len_band * (j) + 2
            temp_array = np.array([float(fbuffer[line + k + 1].split()[10]) for k in range(nion)])
            sum = np.sum(temp_array)
            temp_array = temp_array / sum
            temp_array = temp_array * temp_array
            LIPR = log(np.sum(temp_array)) / log(nion)
            LIPR_list.append(LIPR)
            band_num.append(int(fbuffer[line-2].split()[1]))
            energy.append(float(fbuffer[line-2].split()[4]))
        LIPR_list = np.array(LIPR_list)
        band_num = np.array(band_num)
        energy = np.array(energy) - Ef

        final_list.append([band_num, energy, LIPR_list])

    final_list = np.array(final_list)
    #print len(final_list)

    #final_list < kpoint < band number & energy & LIPR (< means contains)

    # We find VBM position using the criteria of LIPR

    energy_gamma = np.absolute(final_list[0][1])
    VT = np.argmin(energy_gamma) + 1  #- 1 
    CB = VT + 1 

    print VT
    print CB

    temp_CB = np.array([float(final_list[i][1][CB-1]) for i in range(nkpt)])
    temp_VB = np.array([float(final_list[i][1][VT-1]) for i in range(nkpt)])

    E_CB = np.amin(temp_CB)
    E_VT = np.amax(temp_VB)

    #print E_CB, E_VT

    #temp = np.array([float(final_list[j][2][VT-1]) for j in range(nkpt)])
    #print temp

    init = VT
    for i in range(1000):
        temp = np.array([float(final_list[j][2][init-1]) for j in range(nkpt)]) - criteria
        if any (temp < 0):
            #print 'hhh'
            p = np.where(temp < 0)[0]
            #print init
            #print p
            pp = np.array([float(final_list[k][1][init-1]) for k in p])
            #print pp
            VBM_position = np.amax(pp)
            break
        else:
            init -= 1

    corr_bg = E_CB - VBM_position
    orig_bg = E_CB - E_VT
    VB_tail = E_VT - VBM_position
    CBM_position = E_CB
    VT_position = E_VT

    g = open('LIPR.out','w')

    first_line = str(VT_position)+'\t'+str(VBM_position)+'\t'+str(CBM_position)+'\n'
    g.write(first_line)

    for i in range(nkpt):
        for j in range(nband):
            templine = str(final_list[i][1][j]) + '\t' + str(final_list[i][2][j]) +'\n'
            g.write(templine)
    g.close()

    Energy_VT = np.array([float(final_list[i][1][VT-1]) for i in range(nkpt)])
    LIPR_VT = np.array([float(final_list[i][2][VT-1]) for i in range(nkpt)])
    inv_LIPR_VT = 1.0/LIPR_VT
    weight_VT = inv_LIPR_VT / np.sum(inv_LIPR_VT)
    VT_LIPR_av = np.sum(weight_VT*Energy_VT)
    print 'LIPR weight averaged energy = '+str(VT_LIPR_av)
    #print Energy_VT
    #print LIPR_VT

    print 'Band # of Valence Band Top is ' + str(VT) + ' and Eigenvalue (eV) is ' + str(E_VT)
    print 'Band # of Conduction Band Bottom is ' + str(CB) + ' and Eigenvalue (eV) is ' + str(E_CB)
    print 'Difference between VT and CB is ' + str(orig_bg)
    print 'VBM_position is ' + str(VBM_position) + ' corrected band gap by LIPR is ' + str(corr_bg)
    print 'valence band tail is ' + str(VB_tail)
 
    return final_list, VT_position, VBM_position, CBM_position, orig_bg, corr_bg

def plot_LIPR(final_list, plot_energy_range,criteria, VT_position, VBM_position, CBM_position, option=0):
    nkpt = len(final_list)
    fig = plt.figure()

    if option == 0:
        corr = 0 # Ef
        VBM_position = VBM_position
        CBM_position = CBM_position
        VT_position = VT_position
    elif option ==1:
        corr = 0 + VBM_position  # VBM
        VBM_position = VBM_position - corr
        CBM_position = CBM_position - corr
        VT_position = VT_position - corr
    else:
        print 'Option ERROR'

    for i in range(nkpt):
        xlist = np.array(final_list[i][1]) - corr
        ylist = final_list[i][2]
        plt.plot(xlist,ylist,linestyle="none", marker="o", ms=3, markerfacecolor="black", markeredgecolor = "black")

    plt.axis(plot_energy_range)
    plt.plot([plot_energy_range[0],plot_energy_range[1]],[criteria, criteria], '--', color='red')
    plt.plot([VBM_position, VBM_position],[plot_energy_range[2],plot_energy_range[3]], color='blue')
    plt.plot([CBM_position, CBM_position],[plot_energy_range[2],plot_energy_range[3]], color='blue')
    plt.plot([VT_position,VT_position], [plot_energy_range[2],plot_energy_range[3]], color='green')

    fig.savefig('LIPR.png')



#Ef = get_Ef('DOSCAR')
#final_list, VT_position, VBM_position, CBM_position, orig_bg, corr_bg = LIPR('CONTCAR', 'PROCAR', -0.7, Ef)
#print Ef + VBM_position
#plot_LIPR(final_list, [-5.0,5.0,-1.0,-0.2],-0.7, VT_position, VBM_position, CBM_position, option= 0)
#########################################
#print data_structure('CONTCAR')
#Ef, VBM_position, total_out_dos, projected_out_dos = DOS('CONTCAR', 'DOSCAR', np.array([1]), option=0)
#plot_DOS(total_out_dos, [-5.0,5.0,0,50])
#plt.show()
#plot_DOS(projected_out_dos, [-10.0,10.0,0,3])
#lattice, position, atom, num_atom = data_structure('CONTCAR')
#print get_rdf(lattice, position, 66, 85, 5, 0.1, np.array([2, 2, 2]))
#r, pair_corr = get_pair_corr('CONTCAR', 5, 0.05, np.array([2, 2, 2]), 'Zn', 'Sn')
#plot_pair_corr(r, pair_corr)
average_Nc_bc('CONTCAR', 2.3, np.array([2, 2, 2]), 'Zn', 'O')

#Nc_oxygen_with_metal('CONTCAR', [2.5,2.3,2.3], [2,2,2])
