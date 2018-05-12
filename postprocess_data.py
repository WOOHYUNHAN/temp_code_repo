import os.path
import time
from math import *
from math import sqrt
from numpy.linalg import norm
import numpy as np
import re
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import matplotlib as mpl

keys = {

            "figure.figsize": np.array([10, 10]),
            "font.size": 12,
            "axes.labelsize": 11,
            "axes.titlesize": 12,
            "axes.linewidth": 0.8 ,    # edge linewidth
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 10,

            "grid.linewidth": 1,
            "lines.linewidth": 1.75,
            "patch.linewidth": .3,
            "lines.markersize": 7,
            "lines.markeredgewidth": 0,

            "xtick.major.width": 1,
            "ytick.major.width": 1,
            "xtick.minor.width": .5,
            "ytick.minor.width": .5,

            "xtick.major.pad": 7,
            "ytick.major.pad": 7,

"xtick.direction"      : "in",    # direction: in, out, or inout
"ytick.direction"      : "in",    # direction: in, out, or inout            
"xtick.top"            : True,   # draw ticks on the top side
"xtick.bottom"         : True ,  # draw ticks on the bottom side
"ytick.left"           : True ,  # draw ticks on the left side
"ytick.right"          : True , # draw ticks on the right side

"figure.subplot.left"    : 0.08 , # the left side of the subplots of the figure
"figure.subplot.right"   : 0.95 ,   # the right side of the subplots of the figure
"figure.subplot.bottom"  : 0.06 ,   # the bottom of the subplots of the figure
"figure.subplot.top"     : 0.97 ,   # the top of the subplots of the figure
"figure.subplot.wspace"  : 0.3  ,  # the amount of width reserved for blank space between subplots,
                                 # expressed as a fraction of the average axis width
"figure.subplot.hspace"  : 0.45 ,   # the amount of height reserved for white space between subplots,
                                 # expressed as a fraction of the average axis height

"font.family"         : "sans-serif",
}

# For more information about keys : http://matplotlib.org/users/customizing.html

mpl.rcParams.update(keys)


def cal_distance(latt_a,latt_b,latt_c,pos):
    cart = np.dot(np.array([latt_a,latt_b,latt_c]), pos)
    #cart = pos[0] * np.array(latt_a) + pos[1] * np.array(latt_b) + pos[2] * np.array(latt_c)
    distance = np.linalg.norm(cart)
    return distance


def print_line(start,end):
    for i in range(start,end+1):
        #line = '/home/users/hanwooh/POT_fitting_aenet_boron/12_morestructures/structure_set_B10_mB16_T52/'+ '%04d' %(i) +'.xsf'
        line = '/home/users/hanwooh/POT_fitting_aenet_carbon/01_initial_training_set/10_bctC4/MD_n20gpa_1500K_400steps/'+ '%04d' %(i) +'.xsf'
        print line

def read_xsf_info(start,end):
    g = open('xsf_info.out','w')
    g.write('number'+'\t'+'atom'+'\t'+'total_energy'+'\t'+'DFT_energy(eV/atom)'+'\t'+'MAF_x(eV/A)'+'\t'+'MAF_y(eV/A)'+'\t'+'MAF_z(eV/A)'+'\t'+'RMS(eV/A)'+'\t'+'vol(A3)'+'\n')
    for i in range(start, end+1):
        print i
        name= '%04d' %(i) +'.xsf'
        #name= 'structure' + '%04d' %(i) +'.xsf'
        f=open(name,'r')
        temp = f.readlines()
        f.close()
        a = np.array([float(temp[4].split()[h]) for h in range(3)])
        b = np.array([float(temp[5].split()[h]) for h in range(3)])
        c = np.array([float(temp[6].split()[h]) for h in range(3)])
        vol = np.inner(a,np.cross(b,c))
        tot_en = float(temp[0].split()[4])
        num = int(temp[8].split()[0])
        force = np.zeros(3)
        force2 = np.zeros(3)
        for j in range(num):
            temp_force = np.array([abs(float(temp[9+j].split()[k+4])) for k in range(3)])
            #temp_force2 = np.array([pow(float(temp[9+j].split()[k+4]),2) for k in range(3)])
            force += temp_force
            force2 += temp_force * temp_force
        RMS_f = sqrt(np.sum(force2)/num)
        force = force/num
        #final_force = np.linalg.norm(force)
        #print force
        line = str(i)+'\t'+str(num)+'\t'+str(tot_en)+'\t'+str(tot_en/num)+'\t'+str(force[0])+'\t'+str(force[1])+'\t'+str(force[2])+'\t'+str(RMS_f)+'\t'+str(vol)+'\n'
        g.write(line)
    g.close()
    return 0

def read_predict_info(start,end):
    g = open('predict_info.out','w')
    g.write('number'+'\t'+'atom'+'\t'+'total_energy'+'\t'+'ANN_energy(eV/atom)'+'\t'+'MAF_x(eV/A)'+'\t'+'MAF_y(eV/A)'+'\t'+'MAF_z(eV/A)'+'\t'+'RMS_f(eV/A)'+'\t'+'vol(A3)'+'\n')
    f = open('predict.out')
    temp = f.readlines()
    f.close()
    index=[]
    for i in range(len(temp)):
        if not temp[i].split() == []:
            #print i
            if temp[i].split()[0] == 'Number' and temp[i].split()[1] == 'of' and temp[i].split()[2] == 'atoms':
                index.append(i)
    #print len(index)
    for i in range(start,end+1):
        #print i
        line = index[i-1]
        num = int(temp[line].split()[4])
        tot_en = float(temp[line+13+num+1+2].split()[3])
        force = np.array([float(temp[line+13+num+1+2+2].split()[4+j]) for j in range(3)])
        RMS = float(temp[line+13+num+1+2+2+2].split()[3])
        final_force = np.linalg.norm(force)
        a = np.array([float(temp[line+5].split()[k+3]) for k in range(3)])
        b = np.array([float(temp[line+6].split()[k+3]) for k in range(3)])
        c = np.array([float(temp[line+7].split()[k+3]) for k in range(3)])
        vol = np.inner(a,np.cross(b,c))
        line = str(i)+'\t'+str(num)+'\t'+str(tot_en)+'\t'+str(tot_en/num)+'\t'+str(force[0])+'\t'+str(force[1])+'\t'+str(force[2])+'\t'+str(RMS)+'\t'+str(vol)+'\n'
        g.write(line)
    g.close()

def from_xsf_to_POSCAR(start,end, xsf_start, poscar_start, variations):
    for i in range(start, end+1):
        xsf_name = '%04d' %(i+xsf_start) +'.xsf'
        POSCAR_name = 'POSCAR_'+'%05d' %(i+poscar_start)
        f = open(xsf_name)
        temp = f.readlines()
        f.close()
        g = open(POSCAR_name, 'w')
        tot_en = float(temp[0].split()[4])
        num_atom = int(temp[8].split()[0])
        #print tot_en
        g.write(str(tot_en)+'\n')
        g.write(str(1.0000)+'\n')
        for j in range(3):
            g.write(temp[4+j])
        g.write('Si'+'\n')
        g.write(str(num_atom)+'\n')
        g.write('Cart'+'\n')
        for j in range(num_atom):
            random_move = (np.random.rand(3) - 1.0) * variations
            temp_f = temp[j+9].split()
            line = str(float(temp_f[1]) + random_move[0])+'\t'+str(float(temp_f[2]) + random_move[1])+'\t'+str(float(temp_f[3]) + random_move[2])+'\n'
            g.write(line)
        g.close()

    return 0

def from_xsf_to_xsf2(start,end, xsf_start, poscar_start, variations):
    for i in range(start, end+1):
        xsf_name = '%04d' %(i+xsf_start) +'.xsf'
        xsf2_name = '%06d' %(i+poscar_start) +'.xsf'
        f = open(xsf_name)
        temp = f.readlines()
        f.close()
        g = open(xsf2_name, 'w')
        tot_en = float(temp[0].split()[4])
        num_atom = int(temp[8].split()[0])
        #print tot_en
        for j in range(9):
            g.write(temp[j])
        for j in range(num_atom):
            random_move = (np.random.rand(3) - 1.0) * variations
            temp_f = temp[j+9].split()
            line = str(temp_f[0])+' ' + str(float(temp_f[1]) + random_move[0])+' '+str(float(temp_f[2]) + random_move[1])+' '+str(float(temp_f[3]) + random_move[2]) +' 0.0 0.0 0.0' +      '\n'
            g.write(line)
        g.close()

    return 0


def choose_uncorrected_data(input_DFT_file, input_ANN_file, energy_error):
    f = open(input_DFT_file, 'r')
    tempf = f.readlines()
    f.close()

    g = open(input_ANN_file, 'r')
    tempg = g.readlines()
    g.close()

    data_length = min(len(tempf), len(tempg))
    print data_length

    if len(tempf) != len(tempg):
        print 'DFT and ANN files are not matched up'
        print 'But Process will go assuming that minimum length is standard'
        print 'FYI: ANN length = ' + str(len(tempg)) + ' DFT length = ' + str(len(tempf)) + ' Min legnth = ' + str(data_length)
        #return 0

    f_number = []
    f_energy = []
    f_RMS = []

    for i in range(1, data_length):
        temp = tempf[i].split()
        f_number.append(int(temp[0]))
        f_energy.append(float(temp[3]))
        f_RMS.append(float(temp[7]))

    g_number = []
    g_energy = []
    g_RMS = []



    for i in range(1, data_length):
        temp = tempg[i].split()
        g_number.append(int(temp[0]))
        g_energy.append(float(temp[3]))
        g_RMS.append(float(temp[7]))




    delta_energy = np.absolute(np.array(f_energy) - np.array(g_energy))

    uncorrected_index = []

    for i in range(len(delta_energy)):
        if delta_energy[i] > energy_error:
            uncorrected_index.append(f_number[i])

    # print what you know
    print 'Information'
    print 'Enery error : ' + str(energy_error)
    num_unco =  '# of uncorrect  : ' + str(len(uncorrected_index)) + ' / ' + str(len(delta_energy))
    print num_unco
    mae = 'MAE E : ' + '%03d' %(np.average(delta_energy)*1000) + ' meV/atom'
    print mae

    # Draw nice and simple figure
    gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[1, 1])
    ax1 = plt.subplot(gs[0,0])

    for i in range(1, data_length):
        if i in uncorrected_index:
            ax1.plot(f_energy[i-1], g_energy[i-1], 'r.')
            #pass
        else:
            ax1.plot(f_energy[i-1], g_energy[i-1], 'k.')
            #pass
    ax1.plot([1000000000,1000000000],[1000000000,1000000000], 'r.', label=r'Uncorrect within '+str(energy_error)+ ' eV/atom')
    ax1.plot([1000000000,1000000000],[1000000000,1000000000], 'k.', label=r'Correct within '+str(energy_error)+ ' eV/atom')
    ax1.plot([1000000000,1000000000],[1000000000,1000000000], linestyle='none', label=num_unco)
    ax1.plot([1000000000,1000000000],[1000000000,1000000000], linestyle='none', label=mae)
    ax1.set_xlabel(r'DFT energy (eV/atom)')
    ax1.set_ylabel(r'ANN energy (eV/atom)')
    ax1.plot([min(np.min(f_energy), np.min(g_energy)),max(np.max(f_energy),np.max(g_energy))], [min(np.min(f_energy), np.min(g_energy)),max(np.max(f_energy),np.max(g_energy))], 'k--', linewidth=1)
    ax1.set_xlim(min(np.min(f_energy), np.min(g_energy)),max(np.max(f_energy),np.max(g_energy)))
    ax1.set_ylim(min(np.min(f_energy), np.min(g_energy)),max(np.max(f_energy),np.max(g_energy)))
    handles, labels = ax1.get_legend_handles_labels()
    leg=ax1.legend(handles[::-1], labels[::-1], frameon=True ,loc='upper left',numpoints=1,handletextpad=0.5,borderpad=0.05, ncol=1, labelspacing=0.3, handlelength=2, prop={'size':10})

    ax2 = plt.subplot(gs[0,1])

    for i in range(1, data_length):
        if i in uncorrected_index:
            ax2.plot(f_RMS[i-1], g_RMS[i-1], 'r.')
            #pass
        else:
            ax2.plot(f_RMS[i-1], g_RMS[i-1], 'k.')
            #pass
    
    ax2.set_xlabel(r'DFT force (eV/A)')
    ax2.set_ylabel(r'ANN force (eV/A)')
    ax2.plot([0,1], [0,1], 'k--', linewidth=1)
    ax2.set_xlim(0,1)
    ax2.set_ylim(0,1)
    leg=ax2.legend(handles[::-1], labels[::-1], frameon=True ,loc='upper left',numpoints=1,handletextpad=0.5,borderpad=0.05, ncol=1, labelspacing=0.3, handlelength=2, prop={'size':10})

    ax3 = plt.subplot(gs[1,0])

    x = np.linspace(0.0, 1.0, num=len(delta_energy))

    for i in range(1, data_length):
        ax3.plot(f_energy[i-1], g_energy[i-1], '.', linestyle='none', markerfacecolor=(0.0+x[i-1],0.0,1.0-x[i-1],0.5), markeredgecolor='white')

 
    ax3.set_xlabel(r'DFT energy (eV/atom)')
    ax3.set_ylabel(r'ANN energy (eV/atom)')
    ax3.plot([min(np.min(f_energy), np.min(g_energy)),max(np.max(f_energy),np.max(g_energy))], [min(np.min(f_energy), np.min(g_energy)),max(np.max(f_energy),np.max(g_energy))], 'k--', linewidth=1)
    ax3.set_xlim(min(np.min(f_energy), np.min(g_energy)),max(np.max(f_energy),np.max(g_energy)))
    ax3.set_ylim(min(np.min(f_energy), np.min(g_energy)),max(np.max(f_energy),np.max(g_energy)))
   
    ax4 = plt.subplot(gs[1,1])

    for i in range(1, data_length):
        ax4.plot(f_RMS[i-1], g_RMS[i-1], '.', linestyle='none', markerfacecolor=(0.0+x[i-1],0.0,1.0-x[i-1],0.5), markeredgecolor='white')
    
    ax4.set_xlabel(r'DFT force (eV/A)')
    ax4.set_ylabel(r'ANN force (eV/A)')
    ax4.plot([0,1], [0,1], 'k--', linewidth=1)
    ax4.set_xlim(0,1)
    ax4.set_ylim(0,1)
  
    plt.savefig("test.png")
    plt.show()
    ##########################################################

    return uncorrected_index

def make_training_set_from_error(uncorrected_index, num_var, move):
    out_count = 1
    for i in range(len(uncorrected_index)):
        index = int(uncorrected_index[i])
        for j in range(num_var):
            if j == 0:
                #from_xsf_to_POSCAR(0, 0, index, out_count, 0.0)
                from_xsf_to_xsf2(0, 0, index, out_count, 0.0)
            else:
                #from_xsf_to_POSCAR(0, 0, index, out_count, move)
                from_xsf_to_xsf2(0, 0, index, out_count, move)
            out_count += 1


def analysis_data(input_DFT_file, input_ANN_file, energy_error):
    f = open(input_DFT_file, 'r')
    tempf = f.readlines()
    f.close()

    g = open(input_ANN_file, 'r')
    tempg = g.readlines()
    g.close()

    if len(tempf) != len(tempg):
        print 'DFT and ANN files are not matched up'
        return 0

    f_number = []
    f_energy = []
    f_RMS = []

    for i in range(1, len(tempf)):
        temp = tempf[i].split()
        f_number.append(int(temp[0]))
        f_energy.append(float(temp[3]))
        f_RMS.append(float(temp[7]))

    g_number = []
    g_energy = []
    g_RMS = []

    for i in range(1, len(tempg)):
        temp = tempg[i].split()
        g_number.append(int(temp[0]))
        g_energy.append(float(temp[3]))
        g_RMS.append(float(temp[7]))


    delta_energy = np.absolute(np.array(f_energy) - np.array(g_energy))


    uncorrected_index = []

    for i in range(len(delta_energy)):
        if delta_energy[i] > energy_error:
            uncorrected_index.append(f_number[i])

    gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[1, 1])
    ax1 = plt.subplot(gs[0,0])

    x = np.linspace(0.0, 1.0, num=len(delta_energy))

    dft_correct = f_energy - np.min(f_energy)

    for i in f_number:
        ax1.plot(dft_correct[i-1], delta_energy[i-1], '.', linestyle='none', markerfacecolor=(0.0+x[i-1],0.0,1.0-x[i-1],0.5), markeredgecolor='white')


    ax1.set_xlabel(r'DFT energy (eV/atom)')
    ax1.set_ylabel(r'Energy difference (eV/atom)')
    ax1.plot([0,len(delta_energy)], [energy_error,energy_error], 'k--', linewidth=1)
    ax1.set_xlim(0,0.3)
    ax1.set_ylim(0,energy_error*2)

    ax2 = plt.subplot(gs[0,1])

    x = np.linspace(0.0, 1.0, num=len(delta_energy))

    for i in f_number:
        ax2.plot(f_RMS[i-1], delta_energy[i-1], '.', linestyle='none', markerfacecolor=(0.0+x[i-1],0.0,1.0-x[i-1],0.5), markeredgecolor='white')


    ax3 = plt.subplot(gs[1,0])

    for i in f_number:
        if i in uncorrected_index:
            ax3.plot(f_number[i-1], f_energy[i-1], 'r.')
            #pass
        else:
            ax3.plot(f_number[i-1], g_energy[i-1], 'k.')
            #pass

    ax3.set_xlim(0,len(f_number))
    ax3.set_ylim(min(np.min(f_energy), np.min(g_energy)),max(np.max(f_energy),np.max(g_energy)))
    plt.show()

def plot_bank_energy_from_csa_out(npop, natom):
    f = open('csa.out', 'r')
    tempf = f.readlines()
    f.close()

    bank_plot = []
    count = 0

    for i in range(len(tempf)):
        temp = tempf[i].split()
        if len(temp) != 0:
            if temp[0] == 'new' or temp[0] == 'old':
                count += 1
                num_line = (npop-1) / 10  + 1
                #print num_line
                bank = []
                for j in range(i+2, i+2+num_line):
                    temp2 = tempf[j].split()
                    for k in range(10):
                        bank.append(float(temp2[k]))
                bank = np.sort(np.array(bank) / natom)
                bank_plot.append([count, bank])

    #print bank_plot[0]
    gs = gridspec.GridSpec(1, 1, width_ratios=[1], height_ratios=[1])
    ax1 = plt.subplot(gs[0,0])

    num_update = len(bank_plot)

    print num_update
    print np.average(bank_plot[num_update-1][1])


    for i in range(num_update):
        for j in range(npop):
            ax1.plot(bank_plot[i][0], bank_plot[i][1][j], '.', linestyle='none', markerfacecolor='black', markeredgecolor='white')

    min_y = np.min([float(bank_plot[i][1][0]) for i in range(num_update)])

    ax1.set_xlabel(r'# of updates')
    ax1.set_ylabel(r'Energy (eV/atom)')
    #ax1.plot([0,len(delta_energy)], [energy_error,energy_error], 'k--', linewidth=1)
    ax1.set_xlim(1, num_update)
    ax1.set_ylim(min_y-0.01,-4.8)

    plt.savefig("bank_plot.png")
    plt.show()



    return 0
plot_bank_energy_from_csa_out(30,32)
#uncorrected_index =  choose_uncorrected_data('xsf_info_DFT_8th_180221.out', 'xsf_info_ANN_8th_180221.out', 0.005)
#make_training_set_from_error(uncorrected_index, 1, 0.005)
#analysis_data('xsf_info_DFT_8th.out', 'xsf_info_ANN_8th.out', 0.02)
#read_xsf_info(1,23200)
#read_predict_info(1,23200)
#print_line(1,300)
