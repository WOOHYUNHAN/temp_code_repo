import os.path
import time
from math import *
from math import sqrt
from numpy.linalg import norm
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from mpl_toolkits.mplot3d import Axes3D
#from numeric import *


def transform(latt_a, latt_b, latt_c, pos):
    latt = np.transpose(np.array([latt_a,latt_b,latt_c]))
    pos = np.transpose(np.array(pos))

    cart = np.dot(latt,pos)

    return cart

def cal_distance(latt_a,latt_b,latt_c,pos):
    cart = np.dot(np.array([latt_a,latt_b,latt_c]), pos)
    #cart = pos[0] * np.array(latt_a) + pos[1] * np.array(latt_b) + pos[2] * np.array(latt_c)
    distance = np.linalg.norm(cart)
    return distance

def transform_from_radian_to_degree(input_angle):
	ia = float(input_angle)
	output_angle = (180*ia) / np.pi

	return output_angle

def get_atomic_data_from_CONTCAR(filename='CONTCAR'):
    ### read lattice vector and atomic information from CONTCAR
    f = open(filename)
    fbuffer = f.readlines()
    f.close()

    univ = float(fbuffer[1].split()[0])
    latt_a = np.array([float(fbuffer[2].split()[i]) for i in range(3)]) * univ
    latt_b = np.array([float(fbuffer[3].split()[i]) for i in range(3)]) * univ
    latt_c = np.array([float(fbuffer[4].split()[i]) for i in range(3)]) * univ

    vol = np.dot(latt_a, np.cross(latt_b,latt_c))

    recip_a = (2*np.pi/vol) * np.cross(latt_b, latt_c)
    recip_b = (2*np.pi/vol) * np.cross(latt_c, latt_a)
    recip_c = (2*np.pi/vol) * np.cross(latt_a, latt_b)

    atom_name = np.array([str(fbuffer[5].split()[i]) for i in range(len(fbuffer[5].split()))])
    atom_number = np.array([int(fbuffer[6].split()[i]) for i in range(len(fbuffer[5].split()))])


    #print np.linalg.norm(recip_a)/np.linalg.norm(recip_b)
    #print np.linalg.norm(recip_a)/np.linalg.norm(recip_c)
    #print recip_c

    return atom_name, atom_number, latt_a, latt_b, latt_c, recip_a, recip_b, recip_c

def get_Ef_from_DOSCAR(filename='DOSCAR'):
    f= open(filename)
    fbuffer = f.readlines()
    f.close()

    grid = int(fbuffer[5].split()[2]) 
    Ef = float(fbuffer[5].split()[3])

    return Ef

def get_high_symmetry_from_KPOITNS(filename='KPOINTS'):
    ### read input kpoints and their name from KPOINTS
    f = open(filename)
    fbuffer = f.readlines()
    f.close()

    ngrid = int(fbuffer[1].split()[0])

    temp_KPOINTS=[]

    temp_name=[]

    for i in range(4, len(fbuffer)):
        temp=fbuffer[i].split()
        if temp != []:
            #print temp
            temp_kpoint = np.array([float(fbuffer[i].split()[j]) for j in range(3)])
            temp_KPOINTS.append(temp_kpoint)
            temp_name.append(str(temp[4]))


    #print temp_name
    #a =  np.sort(np.unique(temp_name,True)[1])
    #print a

    #KPOINTS=[]
    #name = []

    #for i in a:
    #    KPOINTS.append(temp_KPOINTS[i])
    #    name.append(temp_name[i])
    
    KPOINTS = temp_KPOINTS
    name = temp_name

    if not len(KPOINTS)/2.0 - len(KPOINTS)/2 == 0:
        print 'ERROR: The number of KPOITNS is not even number'


    #print ngrid
    #print len(KPOINTS)
    #print name
    return ngrid, KPOINTS, name

def gen_kpoints_full_list(ngrids, KPOINTS):
    full = []
    if len(KPOINTS)/2.0 - len(KPOINTS)/2 != 0:
        print 'the number of kpoints is not even; plz check KPOITNS file'
        return 0
    else:
        pass

    for i in range(len(KPOINTS)/2):
        delta = (KPOINTS[2*i+1] - KPOINTS[2*i])/float(ngrids-1)
        for j in range(ngrids):       
            temp = KPOINTS[2*i] + j*delta
            full.append(temp)
    return full
    

def get_eigenvalue_from_PROCAR(filename='PROCAR', option_prjected=0):
    ### read eigenvalues from KPOINTS
    f = open(filename)
    fbuffer = f.readlines()
    f.close()   
    ngrids, KPOINTS_list, name = get_high_symmetry_from_KPOITNS()
    #print len(KPOINTS_list)/2
    KPOINTS = gen_kpoints_full_list(ngrids, KPOINTS_list)
    atom_number = get_atomic_data_from_CONTCAR()[1]
    #print atom_number
    #print fbuffer[1].split()

    num_kpoints = int(fbuffer[1].split()[3])
    num_bands = int(fbuffer[1].split()[7])
    num_ions = int(fbuffer[1].split()[11])

    ####option####
    if option_prjected == 0:
        print 'Turn off projected option'
        projected=[]
    elif option_prjected ==1:
        print 'Turn on projected option'
        projected=[int(i) for i in range(1,num_ions+1)]
    else:
        print 'option error'
        return 0


    #print num_kpoints, num_bands, num_ions

     
    len_inter_band = 2 + (1+num_ions+1) + 1
    len_inter_kpoint = 1 + (len_inter_band*num_bands) + 2

    #print len_inter_band
    #print len_inter_kpoint

    #print fbuffer[(3 + len_inter_kpoint*5) + 2+ len_inter_band*5].split()
    total = []
    final = []
    for j in range(num_bands): # band
        total_ffinal = []
        ffinal = []
        for i in range(num_kpoints): # k-point
            total_unit = []
            unit = [] # unit = [band, kpoint, eigenvalue, normalized_portion]
            #kpoints=[float(fbuffer[3 + len_inter_kpoint*i].split()[y+3]) for y in range(3)]
            starting = (3 + len_inter_kpoint*i) + 2+ len_inter_band*j
            energy = float(fbuffer[starting].split()[4])
            #print starting
            unit.append(j+1)
            unit.append(KPOINTS[i])
            unit.append(energy)
            temp_tot = 0
            for k in range(len(atom_number)):
                portion = np.zeros(9)
                for l in range(atom_number[k]):
                    portion += np.array([float(fbuffer[starting+2+int(np.sum(atom_number[0:k]))+(1+l)].split()[x+1]) for x in range(9)])
                temp_tot += np.sum(portion)
                unit.append(portion)  #  unit = [band, k-point, eigenvalue, portion]
            total_unit.append(temp_tot)
            ffinal.append(unit) # final
            total_ffinal.append(total_unit) # sum of every portion each band, kpoints
        #print ffinal
        print 'band-'+str(j+1)+' is finished'
        total.append(total_ffinal)
        final.append(ffinal)

    if len(projected) == 0:
        projected_final = []
    else:
        print 'Calculate projected component'
        projected_final = []
        for j in range(num_bands):
            ffinal = []
            for i in range(num_kpoints):
                unit = []
                #kpoints=[float(fbuffer[3 + len_inter_kpoint*i].split()[y+3]) for y in range(3)]
                starting = (3 + len_inter_kpoint*i) + 2+ len_inter_band*j
                energy = float(fbuffer[starting].split()[4])
                unit.append(j+1)
                unit.append(KPOINTS[i])
                unit.append(energy)
                for l in projected:
                    portion = np.zeros(9)
                    portion += np.array([float(fbuffer[starting+2+int(l)].split()[x+1]) for x in range(9)])
                    unit.append(portion)
                ffinal.append(unit)
            projected_final.append(ffinal)
        print 'projected is finished'

    #print final
    #print '# of atom types= ' + str(len(final[0][0])-3) + ' # of atoms = ' + str(len(projected_final[0][0])-3)
    return num_kpoints, num_bands, total, final, projected_final

def make_band_file(option):
    ngrid, KPOINTS, name = get_high_symmetry_from_KPOITNS()
    Ef = get_Ef_from_DOSCAR()
    atom_name, atom_number, latt_a, latt_b, latt_c, recip_a, recip_b, recip_c = get_atomic_data_from_CONTCAR()
    num_kpoints, num_bands, total, final, projected_final = get_eigenvalue_from_PROCAR(filename='PROCAR', option_prjected=option)
    orbital = ['s','py','pz','px','dxy','dyz','dz2','dxz','dx2']

    if option==0:
        print 'SETTING 1: get band data depending on atom speicies'
        print '======================================'
        print str(atom_name)
        print str(atom_number)
        write_name = atom_name
        write_number = atom_number
        B_data = final
        output_name='band_Ef.dat'
    elif option==1:
        print 'SETTING 2: get band data for specific atom number'
        print 'specific atom data are written as below'
        print '======================================'
        projected_atom_list = np.array(['atom'+str(i) for i in range(1, np.sum(np.array(atom_number))+1)])
        print str(projected_atom_list)
        write_name = projected_atom_list
        write_number = projected_atom_list
        B_data = projected_final
        output_name='projected_band_Ef.dat'
    else:
        print 'option error'
        return 0




    length = 0
    kpoint_length = []
    high_symmetry = []
    #print len(KPOINTS)/2
    for i in range(len(KPOINTS)/2):
        high_symmetry.append(length)
        kpoint_length.append(length)
        #print diff*(ngrid-1)
        #print cal_distance(recip_a, recip_b, recip_c, diff*(ngrid-1))
        for j in range(ngrid-1):
            diff1 = KPOINTS[2*i] + ((KPOINTS[2*i+1] - KPOINTS[2*i]) / float(ngrid-1)) * (j)
            diff2 = KPOINTS[2*i] + ((KPOINTS[2*i+1] - KPOINTS[2*i]) / float(ngrid-1)) * (j+1)
            diff1_cart = transform(recip_a, recip_b, recip_c, diff1)
            diff2_cart = transform(recip_a, recip_b, recip_c, diff2)
            temp = np.linalg.norm(diff2_cart-diff1_cart)
            length += temp
            kpoint_length.append(length)
    high_symmetry.append(length)



    #print len(kpoint_length)
    if num_kpoints != len(kpoint_length):
        print 'ERROR: the number of kpoints is wrong'

    g = open(output_name, 'w')

    line = ''
    for i in range(len(write_name)):
        line += write_name[i] + '\t'
    line += '\n'
    g.write(line) 

    line = str(num_kpoints) + '\t' + str(num_bands) + '\t' + str(ngrid) + '\n'
    g.write(line)

    line = ''
    for i in range(len(high_symmetry)):
        line += str(high_symmetry[i]) + '\t'
    line += '\n'
    g.write(line)

    first_line ='#-band' + '\t' + 'k-point_x' + '\t' + 'k-point_y' + '\t'+ 'k-point_z' + '\t' + 'k_x' + '\t' + 'eigenvalue' + '\t'
    for i in range(len(write_name)):
        for k in range(9):
            name = write_name[i]+'_'+orbital[k]
            first_line += name + '\t'
    first_line += 'total' + '\n'
    g.write(first_line)

    for i in range(num_bands):
        for j in range(num_kpoints):
            temp_line = str(i+1) + '\t' + str(B_data[i][j][1][0]) + '\t' + str(B_data[i][j][1][1]) + '\t' + str(B_data[i][j][1][2]) + '\t' + str(kpoint_length[j]) + '\t' + str(B_data[i][j][2]-Ef) + '\t'
            tot = str(total[i][j][0])
            for k in range(len(write_number)):
                for l in range(9):
                    temp_line += str(B_data[i][j][3 + k][l]) + '\t'
            temp_line += str(tot) + '\n'
            g.write(temp_line)
    g.close()


def draw_Pz_band_HAN(atom):
    #atom: atom number or atom speicies ex) ['B','N']
    #orbital: s, px, py, pz, dxy, dyz, dz2, dxz, dx2
    #xrange(yrange): range of x axis (y axis)
    #option: 0 (total), 1 (projected)
    filename = 'band_Ef.dat'

    fig = plt.figure()
    f = open(filename)
    fbuffer = f.readlines()
    f.close()

    atom_name = np.array([str(fbuffer[0].split()[i]) for i in range(len(fbuffer[0].split()))])
    num_kpoints, num_bands, ngrid = int(fbuffer[1].split()[0]), int(fbuffer[1].split()[1]), int(fbuffer[1].split()[2])
    sym_point = np.array([float(fbuffer[2].split()[i]) for i in range(len(fbuffer[2].split()))])

    atom_list = []
    for i in range(len(atom)):
        s = np.where(atom_name == atom[i])[0][0]
        atom_list.append(s)
    print atom_list

    p_1, p_2 = 2 , 3

    plot = []
    for i in range(num_bands):
        x = [] # k_x
        y = [] # eigenvalues - Ef
        z = [] # portion
        for j in range(num_kpoints):
            starting = 4 + num_kpoints*i + j
            x.append(float(fbuffer[starting].split()[4]))
            y.append(float(fbuffer[starting].split()[5]))
            pz_portion = np.zeros(p_2-p_1)
            for l in atom_list:
                pz_portion += np.array([float(fbuffer[starting].split()[k]) for k in range(6+9*l+p_1, 6+9*l+p_2)])
            total = float(fbuffer[starting].split()[-1])
            z.append(np.sum(pz_portion)/total)
        plot.append([x,y,z])

    #print plot[0][2]

    energy_window = [-2,2]
    #bubble_size = 30

    plt.axis([np.min(sym_point) -0.01 , np.max(sym_point) + 0.01,energy_window[0],energy_window[1]])
    plt.plot([np.min(sym_point) -0.01 , np.max(sym_point) + 0.01], [0,0], linestyle='-', color='red')
    plt.ylabel(r'E - Fermi level (eV)')

    for i in range(num_bands):
        plt.plot(plot[i][0], plot[i][1],linestyle="-", color="black")
        plt.scatter(plot[i][0],plot[i][1], c= np.array(plot[i][2]),s=30,  cmap="coolwarm", edgecolors='none', vmin=0.0, vmax=1.0, alpha=0.5)
    
    for i in range(len(sym_point)):
        plt.plot([sym_point[i], sym_point[i]], energy_window, linestyle="-", color="red")
    fig.savefig('pz_band_Ef.png')
    plt.show()


    return 0

def draw_band_HAN(atom, orbital, option):
    #atom: atom number or atom speicies ex) ['B','N']
    #orbital: s, px, py, pz, dxy, dyz, dz2, dxz, dx2
    #xrange(yrange): range of x axis (y axis)
    #option: 0 (total), 1 (projected)
    if option ==0:
        filename = 'band_Ef.dat'
        print 'BAND STRUCTURE FOR ATOM TYPE'
    elif option ==1:
        filename = 'projected_band_Ef.dat'
        print 'BAND STRUCTURE FOR ATOM NUMBER'
    else:
        print 'ERROR: option is wrong'

    fig = plt.figure()
    f = open(filename)
    fbuffer = f.readlines()
    f.close()

    atom_name = np.array([str(fbuffer[0].split()[i]) for i in range(len(fbuffer[0].split()))])
    num_kpoints, num_bands, ngrid = int(fbuffer[1].split()[0]), int(fbuffer[1].split()[1]), int(fbuffer[1].split()[2])
    sym_point = np.array([float(fbuffer[2].split()[i]) for i in range(len(fbuffer[2].split()))])

    print 'Choose what you want to draw among below elements'
    print str(atom_name)

    atom_list=[]
    if 'projected' in atom_name:
        atom_list=[0]
    else:
        for i in range(len(atom)):
            s = np.where(atom_name == atom[i])[0][0]
            atom_list.append(s)
        #print atom_list

    if orbital == 's':
        p_1, p_2 = 0, 1
        color = 'red'
    elif orbital == 'p':
        p_1, p_2 = 1, 4
        color = 'blue'
    elif orbital == 'd':
        p_1, p_2 = 4, 9
        color = 'cyan'
    elif orbital == 'sp':
        p_1, p_2 = 0, 4
        color = 'green'
    elif orbital == 'spd':
        p_1, p_2 = 0, 9
        color = 'magenta'
    elif orbital == 'py':
        p_1, p_2 = 1, 2
        color = 'green'
    elif orbital == 'pz':
        p_1, p_2 = 2, 3
        color = 'red'  
    elif orbital == 'px':
        p_1, p_2 = 3, 4
        color = 'blue'
    elif orbital == 'dxy':
        p_1, p_2 = 4, 5
        color = 'red'
    elif orbital == 'dyz':
        p_1, p_2 = 5, 6
        color = 'red'
    elif orbital == 'dz2':
        p_1, p_2 = 6, 7
        color = 'red'
    elif orbital == 'dxz':
        p_1, p_2 = 7, 8
        color = 'red'
    elif orbital == 'dx2':
        p_1, p_2 = 8, 9
        color = 'red'

    else:
        print 'ERROR: input for orbital is strange'

    plot = []
    for i in range(num_bands):
        x = [] # k_x
        y = [] # eigenvalues - Ef
        z = [] # portion
        for j in range(num_kpoints):
            starting = 4 + num_kpoints*i + j
            x.append(float(fbuffer[starting].split()[4]))
            y.append(float(fbuffer[starting].split()[5]))
            temp_portion = np.zeros(p_2-p_1)
            for l in atom_list:
                temp_portion += np.array([float(fbuffer[starting].split()[k]) for k in range(6+9*l+p_1, 6+9*l+p_2)])
            total = float(fbuffer[starting].split()[-1])
            z.append(np.sum(temp_portion)/total)
        plot.append([x,y,z])

    #print plot[0][2]

    energy_window = [-8,8]
    bubble_size = 30

    plt.axis([np.min(sym_point) -0.01 , np.max(sym_point) + 0.01,energy_window[0],energy_window[1]])
    plt.plot([np.min(sym_point) -0.01 , np.max(sym_point) + 0.01], [0,0], linestyle='-', color='red')
    plt.ylabel(r'E - Fermi level (eV)')

    for i in range(num_bands):
        #plt.plot(plot[i][0], plot[i][1],linestyle="-", color="black")
        plt.scatter(plot[i][0],plot[i][1], bubble_size*np.array(plot[i][2]), color=color)
    
    for i in range(len(sym_point)):
        plt.plot([sym_point[i], sym_point[i]], energy_window, linestyle="-", color="red")
    fig.savefig('band_Ef.png')
    plt.show()


    return 0


def draw_band_HAN_simple(Ef,ymin,ymax):
    # Ef is fermi level
    # ymin/ymax is energy widnow with respect to Ef
    #ymin = ymin - Ef
    #ymax = ymax - Ef
    nk, k_array, name = get_high_symmetry_from_KPOITNS()
    #print len(k_array)
    atom_name, atom_number, latt_a, latt_b, latt_c, recip_a, recip_b, recip_c = get_atomic_data_from_CONTCAR()
    print recip_a, recip_b, recip_c
    temp = []
    for i in range(len(k_array)-1): # remove duplications
        if np.all(k_array[i] == k_array[i+1]):
            pass
        else:
            temp.append(i)
    temp.append(len(k_array)-1)
    #print temp
    special_k_array = [k_array[i] for i in temp]

    #print special_k_array

    f = open('EIGENVAL')
    fbuffer = f.readlines()
    f.close()

    info = np.array([int(fbuffer[5].split()[i]) for i in range(len(fbuffer[5].split()))])

    x, nk, nb = int(info[0]), int(info[1]), int(info[2])

    kpoints = []
    #weight = []
    eigenval = []
    for j in range(nb):
        eigenval.append([])

    for i in range(nk):
        line = 7 + (nb+2)*i
        #w = float(fbuffer[line].split()[3])]
        temp_k = np.array([float(fbuffer[line].split()[k]) for k in range(3)])
        kpoints.append(temp_k)
        #weight.append(w)
        for j in range(nb):
            eigenval[j].append(float(fbuffer[line+j+1].split()[1])-Ef)

    #print len(kpoints)
    #print len(eigenval[98])



    #print k_array
    n_section = len(k_array)/2
    #print n_section
    sspecial_x = []
    kk_x = []
    slength = 0
    length = 0
    inter_nk = nk/n_section
    for i in range(n_section):
        #slength += length
        sspecial_x.append(length)
        kk_x.append(length)
        for j in range(inter_nk-1):
            diff1 = k_array[2*i] + ((k_array[2*i+1] - k_array[2*i]) / float(inter_nk-1)) * (j)
            diff2 = k_array[2*i] + ((k_array[2*i+1] - k_array[2*i]) / float(inter_nk-1)) * (j+1)
            diff1_cart = transform(recip_a, recip_b, recip_c, diff1)
            diff2_cart = transform(recip_a, recip_b, recip_c, diff2)
            temp = np.linalg.norm(diff2_cart-diff1_cart)
            length += temp
            kk_x.append(length)
    sspecial_x.append(length)

    print sspecial_x
    print kk_x



    g=open('band_structure_Ef','w')

    line = str(nb) + '\t' + str(nk) + '\n'
    g.write(line)

    line =''
    for i in range(len(sspecial_x)):
        line += str(sspecial_x[i])+'\t'
    line += '\n'
    g.write(line)

    #print len(eigenval[1])
    
    for i in range(nk):
        line = str(kk_x[i]) + '\t'
        for j in range(nb):
            line += str(eigenval[j][i]) + '\t'
        line += '\n'
        g.write(line)

    g.close()

    n_vbm =1
    CBM = min(eigenval[n_vbm+1])
    VBM = max(eigenval[n_vbm])
    pos_CBM = np.argmin(eigenval[n_vbm+1])+1
    pos_VBM = np.argmax(eigenval[n_vbm])+1
    print 'CBM = '+str(CBM) + ' at k = ' + str(pos_CBM) + ' band = ' + str(n_vbm+2) 
    print 'VBM = '+str(VBM) + ' at k = ' + str(pos_VBM) + ' band = ' + str(n_vbm+1)
    print 'Gap = '+str(CBM-VBM)
    #plt.plot([k_x[pos_VBM],k_x[pos_CBM]], [VBM, CBM] , color='red',linewidth='3', linestyle='-')


    #print eigenval
    for i in range(nb):
        if i== 8 or i==9 or i==17 or i==18 :
            print 'ok'
            #plt.plot(kk_x,eigenval[i],color='red',linewidth='3')
            print len(eigenval[i])
    #    elif i== 34 or i== 35 or i== 69 or i==68 :
    #        print 'ok'
    #        plt.plot(kk_x,eigenval[i],color='blue',linewidth='3')
    #    else:
    #        plt.plot(kk_x,eigenval[i],color='black')
    for i in range(nb):
        plt.plot(kk_x,eigenval[i],color='black')
    for i in range(len(sspecial_x)):
        plt.plot([sspecial_x[i],sspecial_x[i]],[-100,100], color='red')
        #print special_x[i]
    xmin = min(sspecial_x)-0.1
    xmax = max(sspecial_x)+0.1
    plt.plot([xmin,xmax],[0,0],color='black',linestyle='--')

    #plt.plot([xmin,xmax],[0.4776,0.4776],color='red',linestyle='--')
    #plt.plot([xmin,xmax],[0.4593,0.45930],color='red',linestyle='--')
    #plt.plot([xmin,xmax],[0.25510,0.25510],color='red',linestyle='--')

    #plt.plot([xmin,xmax],[-0.26390,-0.26390],color='blue',linestyle='--')
    #plt.plot([xmin,xmax],[-0.5237,-0.5237],color='blue',linestyle='--')
    #plt.plot([xmin,xmax],[-0.5491,-0.5491],color='blue',linestyle='--')

    plt.axis([xmin,xmax,ymin,ymax])
    #plt.xais(off)
    plt.show()

    return special_x, k_x, eigenval

def draw_band_HAN_simple_polarized(Ef,ymin,ymax):
    # Ef is fermi level
    # ymin/ymax is energy widnow with respect to Ef
    #ymin = ymin - Ef
    #ymax = ymax - Ef
    nk, k_array, name = get_high_symmetry_from_KPOITNS()
    #print len(k_array)
    atom_name, atom_number, latt_a, latt_b, latt_c, recip_a, recip_b, recip_c = get_atomic_data_from_CONTCAR()
    #print recip_a, recip_b, recip_c
    temp = []
    for i in range(len(k_array)-1): # remove duplications
        if np.all(k_array[i] == k_array[i+1]):
            pass
        else:
            temp.append(i)
    temp.append(len(k_array)-1)
    #print temp
    special_k_array = [k_array[i] for i in temp]

    #print special_k_array

    f = open('EIGENVAL')
    fbuffer = f.readlines()
    f.close()

    info = np.array([int(fbuffer[5].split()[i]) for i in range(len(fbuffer[5].split()))])

    x, nk, nb = int(info[0]), int(info[1]), int(info[2])
    print nk, nb

    kpoints = []
    #weight = []
    eigenval = []
    eigenval_2 = []
    for j in range(nb):
        eigenval.append([])
        eigenval_2.append([])


    for i in range(nk):
        line = 7 + (nb+2)*i
        #w = float(fbuffer[line].split()[3])]
        temp_k = np.array([float(fbuffer[line].split()[k]) for k in range(3)])
        kpoints.append(temp_k)
        #weight.append(w)
        for j in range(nb):
            eigenval[j].append(float(fbuffer[line+j+1].split()[1])-Ef)
            eigenval_2[j].append(float(fbuffer[line+j+1].split()[2])-Ef)

    #print len(kpoints)
    #print len(eigenval[98])



    #print k_array
    n_section = len(k_array)/2
    #print n_section
    sspecial_x = []
    kk_x = []
    slength = 0
    length = 0
    inter_nk = nk/n_section
    for i in range(n_section):
        #slength += length
        sspecial_x.append(length)
        kk_x.append(length)
        for j in range(inter_nk-1):
            diff1 = k_array[2*i] + ((k_array[2*i+1] - k_array[2*i]) / float(inter_nk-1)) * (j)
            diff2 = k_array[2*i] + ((k_array[2*i+1] - k_array[2*i]) / float(inter_nk-1)) * (j+1)
            diff1_cart = transform(recip_a, recip_b, recip_c, diff1)
            diff2_cart = transform(recip_a, recip_b, recip_c, diff2)
            temp = np.linalg.norm(diff2_cart-diff1_cart)
            length += temp
            kk_x.append(length)
    sspecial_x.append(length)

    print sspecial_x
    print kk_x



    g=open('band_structure_Ef_spin1','w')
    g_2=open('band_structure_Ef_spin2','w')

    line = str(nb) + '\t' + str(nk) + '\n'
    g.write(line)
    g_2.write(line)

    line =''
    for i in range(len(sspecial_x)):
        line += str(sspecial_x[i])+'\t'
    line += '\n'
    g.write(line)
    g_2.write(line)

    #print len(eigenval[1])
    
    for i in range(nk):
        line = str(kk_x[i]) + '\t'
        line_2 = str(kk_x[i]) + '\t'
        for j in range(nb):
            line += str(eigenval[j][i]) + '\t'
            line_2 += str(eigenval_2[j][i]) + '\t'
        line += '\n'
        line_2 += '\n'
        g.write(line)
        g_2.write(line_2)

    g.close()
    g_2.close()

    #n_vbm =1
    #CBM = min(eigenval[n_vbm+1])
    #VBM = max(eigenval[n_vbm])
    #pos_CBM = np.argmin(eigenval[n_vbm+1])+1
    #pos_VBM = np.argmax(eigenval[n_vbm])+1
    #print 'CBM = '+str(CBM) + ' at k = ' + str(pos_CBM) + ' band = ' + str(n_vbm+2) 
    #print 'VBM = '+str(VBM) + ' at k = ' + str(pos_VBM) + ' band = ' + str(n_vbm+1)
    #print 'Gap = '+str(CBM-VBM)
    #plt.plot([k_x[pos_VBM],k_x[pos_CBM]], [VBM, CBM] , color='red',linewidth='3', linestyle='-')


    #print eigenval
    for i in range(nb):
        if i == 44 :
            plt.plot(kk_x,eigenval[i],color='black',linestyle='--')
            plt.plot(kk_x,eigenval_2[i],color='red',linestyle='--')
        elif i == 45:
            plt.plot(kk_x,eigenval[i],color='black',linestyle='-.')
            plt.plot(kk_x,eigenval_2[i],color='red',linestyle='-.')        
        else:
            plt.plot(kk_x,eigenval[i],color='black')
            plt.plot(kk_x,eigenval_2[i],color='red')
    for i in range(len(sspecial_x)):
        plt.plot([sspecial_x[i],sspecial_x[i]],[-100,100], color='red')
        #print special_x[i]
    xmin = min(sspecial_x)-0.1
    xmax = max(sspecial_x)+0.1
    plt.plot([xmin,xmax],[0,0],color='black',linestyle='--')

    #plt.plot([xmin,xmax],[0.4776,0.4776],color='red',linestyle='--')
    #plt.plot([xmin,xmax],[0.4593,0.45930],color='red',linestyle='--')
    #plt.plot([xmin,xmax],[0.25510,0.25510],color='red',linestyle='--')

    #plt.plot([xmin,xmax],[-0.26390,-0.26390],color='blue',linestyle='--')
    #plt.plot([xmin,xmax],[-0.5237,-0.5237],color='blue',linestyle='--')
    #plt.plot([xmin,xmax],[-0.5491,-0.5491],color='blue',linestyle='--')

    plt.axis([xmin,xmax,ymin,ymax])
    #plt.xais(off)
    plt.show()

    return special_x, k_x, eigenval

def plot_EIGENVAL(Ef,ymin,ymax):
    f = open('EIGENVAL')
    fbuffer = f.readlines()
    f.close()

    info = np.array([int(fbuffer[5].split()[i]) for i in range(len(fbuffer[5].split()))])

    x, nk, nb = int(info[0]), int(info[1]), int(info[2])

    kpoints = []
    #weight = []
    eigenval = []
    for j in range(nb):
        eigenval.append([])

    for i in range(nk):
        line = 7 + (nb+2)*i
        #w = float(fbuffer[line].split()[3])]
        #temp_k = np.array([float(fbuffer[line].split()[k]) for k in range(3)])
        kpoints.append(i+1)
        #weight.append(w)
        for j in range(nb):
            eigenval[j].append(float(fbuffer[line+j+1].split()[1])-Ef)

    print kpoints

    for i in range(nb):
        plt.plot(kpoints,eigenval[i],'-',color='black')
    
    #plt.axis([0,104,ymin,ymax])
    #plt.plot([-20000,20000],[0,0],color='red')
    #plt.xais(off)
    plt.show()

    #print len(kpoints)
    #print len(eigenval[98])
def shift_qpoint_from_convention_to_primitive(initialpoint):
    f = open('POSCAR_primitive')
    temp = f.readlines()
    f.close()
    univ1 = float(temp[1].split()[0])
    a1 = np.array([float(temp[2].split()[i]) for i in range(3)])*univ1
    b1 = np.array([float(temp[3].split()[i]) for i in range(3)])*univ1
    c1 = np.array([float(temp[4].split()[i]) for i in range(3)])*univ1
#print a1
#print b1
#print c1

    g = open('POSCAR_convention')
    temp = g.readlines()
    g.close()
    univ2 = float(temp[1].split()[0])
    a2 = np.array([float(temp[2].split()[i]) for i in range(3)])*univ2
    b2 = np.array([float(temp[3].split()[i]) for i in range(3)])*univ2
    c2 = np.array([float(temp[4].split()[i]) for i in range(3)])*univ2

    vol1 = np.inner(np.cross(a1,b1),c1)
    vol2 = np.inner(np.cross(a2,b2),c2)

#print np.inner(b2,c2)

    recip_a1 = 2*pi/vol1 * np.cross(b1,c1)
    recip_b1 = 2*pi/vol1 * np.cross(c1,a1)
    recip_c1 = 2*pi/vol1 * np.cross(a1,b1)

    recip_a2 = 2*pi/vol2 * np.cross(b2,c2)
    recip_b2 = 2*pi/vol2 * np.cross(c2,a2)
    recip_c2 = 2*pi/vol2 * np.cross(a2,b2)

#print np.inner(recip_b2,recip_c2)

    recip1 = np.matrix([recip_a1,recip_b1,recip_c1])
    recip2 = np.matrix([recip_a2,recip_b2,recip_c2])



    middle = np.transpose(recip2) * initialpoint
    #print middle
    final = np.linalg.inv(np.transpose(recip1)) * middle
    print str(float(final[0]))+ '\t' + str(float(final[1])) + '\t' +str(float(final[2]))

    return 0

def draw_band_HAN_QE(inputfile,Ef,ymin,ymax):
    f = open(inputfile)
    tempf = f.readlines()
    f.close()
    nbnd = int(tempf[0].split()[2])
    nks = int(tempf[0].split()[4])
    print nbnd, nks

    k_x = []
    band = []
    for i in range(nbnd):
        band.append([])

    dis = 0
    for i in range(nks):
        if i == 0:
            k_x.append(dis)
        else:
            x1 = np.array([float(tempf[2*(i)-1].split()[k]) for k in range(3)])
            x2 = np.array([float(tempf[2*(i+1)-1].split()[k]) for k in range(3)])
            cal_dis = np.linalg.norm(x2 - x1)
            dis += cal_dis
            k_x.append(dis)
        for j in range(nbnd):
            band[j].append(float(tempf[2*(i+1)].split()[j])-Ef)


    g = open('bandplot.QE','w')

    for i in range(nks):
        line = str(k_x[i])
        for j in range(nbnd):
            line += '\t' + str(band[j][i])
        line += '\n'
        g.write(line)
    g.close()


def Read_BerryCurvature(filename, direction, area, threshold):
    f = open(filename, 'r')
    tempf = f.readlines()
    f.close()

    curv = []

    for i in range(len(tempf)):
        #print tempf[i].split()
        curv.append(float(tempf[i].split()[direction]))

    curv = np.array(curv) * area

    positive_curv = (curv + np.abs(curv))  /2.0
    negative_curv = (curv - np.abs(curv))  /2.0

    for i in range(len(curv)):
        if abs(positive_curv[i]) < threshold:
            positive_curv[i] = 0.0
        if abs(negative_curv[i]) < threshold:
            negative_curv[i] = 0.0
        if abs(curv[i]) < threshold:
            curv[i] = 0.0        

    num_points = int(np.sqrt(len(curv)))
    #print num_points

    curv_matrix = curv.reshape((num_points,num_points)).T

    #print curv_matrix[0]

    g = open('out', 'w')

    for i in range(num_points):
        temp = ''
        for j in range(num_points):
            temp += str(curv_matrix[i][j]) + '\t'
        temp += '\n'
        g.write(temp)

    g.close()


    print np.sum(positive_curv)
    print np.sum(negative_curv)

    return 0




    




#Read_BerryCurvature('wannier90-kslice-curv.dat', 2, 2.23E-05, 0.5)
#get_eigenvalue_from_PROCAR()
#atom_name, atom_number, latt_a, latt_b, latt_c, recip_a, recip_b, recip_c = get_atomic_data_from_CONTCAR()
#print recip_a
#ngrids, KPOINTS_list, name = get_high_symmetry_from_KPOITNS()
#print len(KPOINTS_list)/2
#a = gen_kpoints_full_list(ngrids, KPOINTS_list)
#print a[49]
#print cal_distance(latt_a, latt_b, latt_c, np.array([0.5000000000000000,  0.2420979416751834,  0.4686769833746567]))
#make_band_file(option=0)
#draw_Pz_band_HAN(['B'])
#draw_band_HAN(['Si'], 's', option=0)
#draw_band_HAN(['atom39'], 'pz', option=1)
draw_band_HAN_simple(-1.4749,-1,1)
#draw_band_HAN_simple_polarized(-4.113,-1,1)
#plot_EIGENVAL(0, -10, 3)
#draw_band_HAN_QE('Pbands.dat',0.0,-2,2)

#shift_qpoint_from_convention_to_primitive(np.matrix([[0.0],[0.0],[0.0]]))
#shift_qpoint_from_convention_to_primitive(np.matrix([[0.5],[0.0],[0.0]]))
#shift_qpoint_from_convention_to_primitive(np.matrix([[0.33333],[0.33333],[0.0]]))
#shift_qpoint_from_convention_to_primitive(np.matrix([[0.0],[0.0],[0.0]]))
#shift_qpoint_from_convention_to_primitive(np.matrix([[0.0],[0.0],[0.5]]))
#shift_qpoint_from_convention_to_primitive(np.matrix([[0.5],[0.0],[0.5]]))
#shift_qpoint_from_convention_to_primitive(np.matrix([[0.33333],[0.33333],[0.5]]))
#shift_qpoint_from_convention_to_primitive(np.matrix([[0.0],[0.0],[0.5]]))

#shift_qpoint_from_convention_to_primitive(np.matrix([[0.0],[0.0],[0.0]]))
#shift_qpoint_from_convention_to_primitive(np.matrix([[0.5],[0.0],[0.0]]))
#shift_qpoint_from_convention_to_primitive(np.matrix([[0.5],[0.5],[0.0]]))
#shift_qpoint_from_convention_to_primitive(np.matrix([[0.0],[0.5],[0.0]]))
#shift_qpoint_from_convention_to_primitive(np.matrix([[0.0],[0.0],[0.0]]))