import os.path
import time
from math import *
from math import sqrt
from numpy.linalg import norm
import numpy as np
import matplotlib.pyplot as plt
#from numeric import *

def make_bxsf_from_pyTB_2D(load_file, lattice_file, out_tag, Fermi_level, nkz):
    eigval = np.load(load_file)
    lattice = np.load(lattice_file)
    nband = len(eigval)
    nkx = int(np.sqrt(len(eigval[0])))
    nky = int(np.sqrt(len(eigval[0])))
    #print eigval[0][101*8-1]
    eigval_xy = []
    for i in range(nband):
        eigval_xy.append(np.reshape(eigval[i], (nkx,nky)))
    #print eigval_xy[7][100]
    nkz = int(nkz)
    datanum_oneline = int(6)
    latt_vec = np.array([[lattice[0][0], lattice[0][1], 0.0], [lattice[1][0], lattice[1][1], 0.0], [0.0, 0.0, 10.0]])
    recip_fac =  (2*np.pi) / np.dot(latt_vec[0], np.cross(latt_vec[1], latt_vec[2]))
    recip_vec = np.array([np.cross(latt_vec[1], latt_vec[2]), np.cross(latt_vec[2], latt_vec[0]), np.cross(latt_vec[0], latt_vec[1])]) * recip_fac
    recip_vec = recip_vec / recip_vec[0][0]
    #print recip_vec
    ######################################
    filename = 'FS_'+str(out_tag)+'_Ef_%0.2f.bxsf' % Fermi_level
    f = open(filename, 'w')
    f.write(' BEGIN_INFO'+'\n')
    f.write(' #'+'\n')
    f.write(' # this is a Band-XCRYSDEN-Structure-File'+'\n')
    f.write(' # aimed at Visualization of Fermi Surface'+'\n')
    f.write(' #'+'\n')
    f.write(' # Case:   B_fsup.bxsf'+'\n')
    f.write(' #'+'\n')
    f.write(' Fermi Energy:        '+str(Fermi_level)+'\n')
    f.write(' END_INFO'+'\n')
    f.write(' BEGIN_BLOCK_BANDGRID_3D'+'\n')
    f.write(' band_energies'+'\n')
    f.write(' BANDGRID_3D_BANDS'+'\n')
    f.write(' '+str(nband)+'\n')
    f.write(' '+str(nkx)+' '+str(nky)+' '+ str(nkz) +'\n')
    f.write(' '+'0.000000  0.000000  0.000000' + '\n')
    f.write(' '+str(recip_vec[0][0])+' ' + str(recip_vec[0][1])+ ' ' + str(recip_vec[0][2])+ '\n')
    f.write(' '+str(recip_vec[1][0])+' ' + str(recip_vec[1][1])+ ' ' + str(recip_vec[1][2])+ '\n')
    f.write(' '+str(recip_vec[2][0])+' ' + str(recip_vec[2][1])+ ' ' + str(recip_vec[2][2])+ '\n')
    for i in range(nband):
        f.write('BAND: '+ str(i+1)+'\n')
        total_line = int(nkx*nky*nkz) / datanum_oneline
        remaining = int(nkx*nky*nkz) % datanum_oneline
        print total_line, remaining
        index = 0
        for j in range(total_line):
            templine = ' '
            for k in range(datanum_oneline/nkz):
                for l in range(nkz):
                    templine += str(eigval[i][index+k]) + ' '
            templine += '\n'
            f.write(templine)
            index += datanum_oneline/nkz
        #templine = ''
        print index
        if remaining != 0:
            for l in range(remaining/nkz):
                templine = ' '
                for m in range(nkz):
                    templine += str(eigval[i][index]) + ' '
                index += 1
            print index
            templine += '\n'
            f.write(templine)

    f.write('END_BANDGRID_3D'+'\n')
    f.write('END_BLOCK_BANDGRID_3D'+'\n')

    f.close()

def read_reciprocal_lattice(filename):
    f = open(filename)
    tempf = f.readlines()
    f.close()

    a_axis = np.array([float(tempf[15].split()[i]) for i in range(3)])
    b_axis = np.array([float(tempf[16].split()[i]) for i in range(3)])
    c_axis = np.array([float(tempf[17].split()[i]) for i in range(3)])

    recip_lattice = np.array([a_axis, b_axis, c_axis])

    return recip_lattice

def read_bxsf_file_2D(filename, nth_band):
    f = open(filename)
    tempf = f.readlines()
    f.close()

    Ef = float(tempf[7].split()[2])

    nx, ny, nz = int(tempf[13].split()[0]), int(tempf[13].split()[1]), int(tempf[13].split()[2])

    band_line = []

    for i in range(len(tempf)):
        if not tempf[i] == []:
            if tempf[i].split()[0] == 'BAND:':
                band_line.append(i)
    start_line = band_line[nth_band-1]
    interval_line = band_line[1] - band_line[0]

    #print start_line, interval_line

    all_data = []
    for i in range(start_line+1, start_line+interval_line):
        ttemp = tempf[i].split()
        for i in range(len(ttemp)):
            all_data.append(float(ttemp[i])-Ef)
    #print all_data

    bgrid_2D = np.zeros((nx,ny))

    for i in range(nx):
        for j in range(ny):
            for k in range(1):
                #print i*(ny*nz)+j*(nz)+k
                bgrid_2D[i][j] = float(all_data[i*(ny*nz)+j*(nz)+k])


    #print bgrid_2D[0][0]

    return nx, ny, bgrid_2D

def approx_from_delta_to_gaussian(value, sigma):
    approx = (1.0/(sqrt(2*np.pi)*sigma)) * exp(-(value/(2*sigma))**2)
    return float(approx)



def cal_nesting_factor_singleq(filename, nth_band, sigma, q_point):
    nx, ny, bgrid_2D = read_bxsf_file_2D(filename, nth_band)

    unitx = 1.0 / (nx-1)
    unity = 1.0 / (ny-1)

    deltax = int(q_point[0] / unitx)
    deltay = int(q_point[1] / unity)

    #print deltax, deltay

    Final = 0

    for i in range(nx-1):
        for j in range(ny-1):
            x = approx_from_delta_to_gaussian(bgrid_2D[i][j], sigma) * approx_from_delta_to_gaussian(bgrid_2D[(i+deltax)%nx][(j+deltay)%ny], sigma)
            #print x
            Final += x

    return Final

def cal_nesting_factor_lineq(filename, nth_band, sigma, lineq_file, output_file_name):
    recip_lattice = read_reciprocal_lattice(filename)
    #print recip_lattice
    f = open(lineq_file, 'r')
    tempf = f.readlines()
    f.close()

    high_sym_point = []
    number_qpoints = []
    for i in range(len(tempf)):
        temp_q = np.array([float(tempf[i].split()[j]) for j in range(3)])
        high_sym_point.append(temp_q)
        number_qpoints.append(int(tempf[i].split()[3]))


    #print len(tempf)

    g = open(output_file_name, 'w')
    k_length = 0.0
    for i in range(len(high_sym_point)-1):
        for j in range(number_qpoints[i]):
            diff = (high_sym_point[i+1] - high_sym_point[i]) / number_qpoints[i]
            q_point = high_sym_point[i] + diff*j
            nesting = cal_nesting_factor_singleq(filename, nth_band, sigma, q_point)
            #print str(k_length) + '\t' + str(nesting)
            line = str(k_length) + '\t' + str(nesting) + '\n'
            g.write(line)
            cart_diff = np.dot(recip_lattice.T, diff)
            length_diff = np.linalg.norm(cart_diff)
            k_length += length_diff
    
    g.close()

    print 'job done'
    return 0

def cal_nesting_factor_areaq(filename, nth_band, sigma, qpoint_array, output_file_name):
    recip_lattice = read_reciprocal_lattice(filename)
    #print recip_lattice

    q_array = np.load(qpoint_array)

    q_array_direct =  np.dot(np.linalg.inv(recip_lattice.T), q_array.T).T

    nesting_array = []

    #print len(q_array_direct)

    for i in range(len(q_array_direct)):
        nesting = cal_nesting_factor_singleq(filename, nth_band, sigma, q_array_direct[i])
        nesting_array.append(nesting)
        print i

    nesting_array = np.array(nesting_array)

    np.save(output_file_name+'_nestingfunc', nesting_array)

    print 'job done'

    return 0

#make_bxsf_from_pyTB_2D('evals2_NN5_300.npy', 'lattice.npy', 'NN5_300', 2.40, 3)
#make_bxsf_from_pyTB_2D('evals2_NN5_300.npy', 'lattice.npy', 'NN5_300', 2.35, 3)
#make_bxsf_from_pyTB_2D('evals2_NN5_300.npy', 'lattice.npy', 'NN5_300', 2.30, 3)
#make_bxsf_from_pyTB_2D('evals2_NN5_300.npy', 'lattice.npy', 'NN5_300', 2.25, 3)
#make_bxsf_from_pyTB_2D('evals2_NN5_300.npy', 'lattice.npy', 'NN5_300', 2.20, 3)
#make_bxsf_from_pyTB_2D('evals2_NN5_300.npy', 'lattice.npy', 'NN5_300', 2.15, 3)
#make_bxsf_from_pyTB_2D('evals2_NN5_300.npy', 'lattice.npy', 'NN5_300', 2.10, 3)
#make_bxsf_from_pyTB_2D('evals2_NN5_300.npy', 'lattice.npy', 'NN5_300', 2.05, 3)
#make_bxsf_from_pyTB_2D('evals2_NN5_300.npy', 'lattice.npy', 'NN5_300', 2.00, 3)
#make_bxsf_from_pyTB_2D('evals2_NN5_300.npy', 'lattice.npy', 'NN5_300', 1.95, 3)
#make_bxsf_from_pyTB_2D('evals2_NN5_300.npy', 'lattice.npy', 'NN5_300', 1.90, 3)
#make_bxsf_from_pyTB_2D('evals2_NN5_300.npy', 'lattice.npy', 'NN5_300', 1.85, 3)
#make_bxsf_from_pyTB_2D('evals2_NN5_300.npy', 'lattice.npy', 'NN5_300', 1.80, 3)
#make_bxsf_from_pyTB_2D('evals2_NN5_300.npy', 'lattice.npy', 'NN5_300', 1.75, 3)
#make_bxsf_from_pyTB_2D('evals2_NN5_300.npy', 'lattice.npy', 'NN5_300', 1.70, 3)



#read_bxsf_file_2D('B_fsup.bxsf', 3)
band = 3
sigma = 0.05 # unit = eV 0.001 Ry = 0.013 eV

#cal_nesting_factor_areaq('FS_Ef_2.1.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'NNN_10_Ef_2.1')
#cal_nesting_factor_areaq('FS_Ef_2.0.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'NNN_10_Ef_2.0')
#cal_nesting_factor_areaq('FS_Ef_1.9.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'NNN_10_Ef_1.9')

#print cal_nesting_factor_singleq('FS_Ef_2.2.bxsf', band, sigma, [0.0, 0.0, 0.0])
#print cal_nesting_factor_singleq('FS_Ef_2.15.bxsf', band, sigma,[0.0, 0.0, 0.0])
#print cal_nesting_factor_singleq('FS_Ef_2.1.bxsf', band, sigma, [0.0, 0.0, 0.0])
#print cal_nesting_factor_singleq('FS_Ef_2.05.bxsf', band, sigma,[0.0, 0.0, 0.0])
#print cal_nesting_factor_singleq('FS_Ef_2.0.bxsf', band, sigma, [0.0, 0.0, 0.0])
#print cal_nesting_factor_singleq('FS_Ef_1.95.bxsf', band, sigma,[0.0, 0.0, 0.0])
#print cal_nesting_factor_singleq('FS_Ef_1.9.bxsf', band, sigma, [0.0, 0.0, 0.0])
#print cal_nesting_factor_singleq('FS_Ef_1.85.bxsf', band, sigma,[0.0, 0.0, 0.0])
#print cal_nesting_factor_singleq('FS_Ef_1.8.bxsf', band, sigma, [0.0, 0.0, 0.0])
#print cal_nesting_factor_singleq('FS_Ef_1.75.bxsf', band, sigma, [0.0, 0.0, 0.0])
#print cal_nesting_factor_singleq('FS_Ef_1.7.bxsf', band, sigma, [0.0, 0.0, 0.0])
#print cal_nesting_factor_singleq('FS_Ef_1.65.bxsf', band, sigma, [0.0, 0.0, 0.0])
#print cal_nesting_factor_singleq('FS_Ef_1.6.bxsf', band, sigma, [0.0, 0.0, 0.0])
#print cal_nesting_factor_singleq('FS_NN5_300_Ef_2.40.bxsf', band, sigma, [0.333,0.333,0.0])
#print cal_nesting_factor_singleq('FS_NN5_300_Ef_2.35.bxsf', band, sigma, [0.333,0.333,0.0])
#print cal_nesting_factor_singleq('FS_NN5_300_Ef_2.30.bxsf', band, sigma, [0.333,0.333,0.0])
#print cal_nesting_factor_singleq('FS_NN5_300_Ef_2.25.bxsf', band, sigma, [0.333,0.333,0.0])
print cal_nesting_factor_singleq('FS_NN20_300_Ef_2.20.bxsf', band, sigma, [0.333,0.333,0.0])
print cal_nesting_factor_singleq('FS_NN20_300_Ef_2.15.bxsf', band, sigma, [0.333,0.333,0.0])
print cal_nesting_factor_singleq('FS_NN20_300_Ef_2.10.bxsf', band, sigma, [0.333,0.333,0.0])
print cal_nesting_factor_singleq('FS_NN20_300_Ef_2.05.bxsf', band, sigma, [0.333,0.333,0.0])
print cal_nesting_factor_singleq('FS_NN20_300_Ef_2.00.bxsf', band, sigma, [0.333,0.333,0.0])
print cal_nesting_factor_singleq('FS_NN20_300_Ef_1.95.bxsf', band, sigma, [0.333,0.333,0.0])
print cal_nesting_factor_singleq('FS_NN20_300_Ef_1.90.bxsf', band, sigma, [0.333,0.333,0.0])
print cal_nesting_factor_singleq('FS_NN20_300_Ef_1.85.bxsf', band, sigma, [0.333,0.333,0.0])
print cal_nesting_factor_singleq('FS_NN20_300_Ef_1.80.bxsf', band, sigma, [0.333,0.333,0.0])
print cal_nesting_factor_singleq('FS_NN20_300_Ef_1.75.bxsf', band, sigma, [0.333,0.333,0.0])
print cal_nesting_factor_singleq('FS_NN20_300_Ef_1.70.bxsf', band, sigma, [0.333,0.333,0.0])

#cal_nesting_factor_areaq('FS_NN20_300_Ef_1.95.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN20_300_Ef_1.95')
#cal_nesting_factor_areaq('FS_NN10_300_Ef_1.90.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN10_300_Ef_1.90')
#cal_nesting_factor_areaq('FS_NN5_300_Ef_1.85.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN5_300_Ef_1.85')


#cal_nesting_factor_lineq('FS_NN20_Ef_2.20.bxsf', band, sigma, 'lineq.dat', 'nesting_line_NNN_20_Ef_2.20.out')
#cal_nesting_factor_lineq('FS_NN20_Ef_2.15.bxsf', band, sigma, 'lineq.dat', 'nesting_line_NNN_20_Ef_2.15.out')
#cal_nesting_factor_lineq('FS_NN20_Ef_2.10.bxsf', band, sigma, 'lineq.dat', 'nesting_line_NNN_20_Ef_2.10.out')
#cal_nesting_factor_lineq('FS_NN20_Ef_2.05.bxsf', band, sigma, 'lineq.dat', 'nesting_line_NNN_20_Ef_2.05.out')
#cal_nesting_factor_lineq('FS_NN20_Ef_2.00.bxsf', band, sigma, 'lineq.dat', 'nesting_line_NNN_20_Ef_2.00.out')
#cal_nesting_factor_lineq('FS_NN20_Ef_1.95.bxsf', band, sigma, 'lineq.dat', 'nesting_line_NNN_20_Ef_1.95.out')
#cal_nesting_factor_lineq('FS_NN20_Ef_1.90.bxsf', band, sigma, 'lineq.dat', 'nesting_line_NNN_20_Ef_1.90.out')
#cal_nesting_factor_lineq('FS_NN20_Ef_1.85.bxsf', band, sigma, 'lineq.dat', 'nesting_line_NNN_20_Ef_1.85.out')
#cal_nesting_factor_lineq('FS_NN20_Ef_1.80.bxsf', band, sigma, 'lineq.dat', 'nesting_line_NNN_20_Ef_1.80.out')
#cal_nesting_factor_lineq('FS_NN20_Ef_1.75.bxsf', band, sigma, 'lineq.dat', 'nesting_line_NNN_20_Ef_1.75.out')
#cal_nesting_factor_lineq('FS_NN20_Ef_1.70.bxsf', band, sigma, 'lineq.dat', 'nesting_line_NNN_20_Ef_1.70.out')

#cal_nesting_factor_lineq('FS_NN20_300_Ef_2.20.bxsf', band, sigma, 'lineq.dat', 'nesting_line_NNN_20_300_Ef_2.20.out')
#cal_nesting_factor_lineq('FS_NN20_300_Ef_2.15.bxsf', band, sigma, 'lineq.dat', 'nesting_line_NNN_20_300_Ef_2.15.out')
#cal_nesting_factor_lineq('FS_NN20_300_Ef_2.10.bxsf', band, sigma, 'lineq.dat', 'nesting_line_NNN_20_300_Ef_2.10.out')
#cal_nesting_factor_lineq('FS_NN20_300_Ef_2.05.bxsf', band, sigma, 'lineq.dat', 'nesting_line_NNN_20_300_Ef_2.05.out')
#cal_nesting_factor_lineq('FS_NN20_300_Ef_2.00.bxsf', band, sigma, 'lineq.dat', 'nesting_line_NNN_20_300_Ef_2.00.out')
#cal_nesting_factor_lineq('FS_NN20_300_Ef_1.95.bxsf', band, sigma, 'lineq.dat', 'nesting_line_NNN_20_300_Ef_1.95.out')
#cal_nesting_factor_lineq('FS_NN20_300_Ef_1.90.bxsf', band, sigma, 'lineq.dat', 'nesting_line_NNN_20_300_Ef_1.90.out')
#cal_nesting_factor_lineq('FS_NN20_300_Ef_1.85.bxsf', band, sigma, 'lineq.dat', 'nesting_line_NNN_20_300_Ef_1.85.out')
#cal_nesting_factor_lineq('FS_NN20_300_Ef_1.80.bxsf', band, sigma, 'lineq.dat', 'nesting_line_NNN_20_300_Ef_1.80.out')
#cal_nesting_factor_lineq('FS_NN20_300_Ef_1.75.bxsf', band, sigma, 'lineq.dat', 'nesting_line_NNN_20_300_Ef_1.75.out')
#cal_nesting_factor_lineq('FS_NN20_300_Ef_1.70.bxsf', band, sigma, 'lineq.dat', 'nesting_line_NNN_20_300_Ef_1.70.out')

#cal_nesting_factor_lineq('FS_NN15_Ef_2.20.bxsf', band, sigma, 'lineq.dat', 'nesting_line_NNN_15_Ef_2.20.out')
#cal_nesting_factor_lineq('FS_NN15_Ef_2.15.bxsf', band, sigma, 'lineq.dat', 'nesting_line_NNN_15_Ef_2.15.out')
#cal_nesting_factor_lineq('FS_NN15_Ef_2.10.bxsf', band, sigma, 'lineq.dat', 'nesting_line_NNN_15_Ef_2.10.out')
#cal_nesting_factor_lineq('FS_NN15_Ef_2.05.bxsf', band, sigma, 'lineq.dat', 'nesting_line_NNN_15_Ef_2.05.out')
#cal_nesting_factor_lineq('FS_NN15_Ef_2.00.bxsf', band, sigma, 'lineq.dat', 'nesting_line_NNN_15_Ef_2.00.out')
#cal_nesting_factor_lineq('FS_NN15_Ef_1.95.bxsf', band, sigma, 'lineq.dat', 'nesting_line_NNN_15_Ef_1.95.out')
#cal_nesting_factor_lineq('FS_NN15_Ef_1.90.bxsf', band, sigma, 'lineq.dat', 'nesting_line_NNN_15_Ef_1.90.out')
#cal_nesting_factor_lineq('FS_NN15_Ef_1.85.bxsf', band, sigma, 'lineq.dat', 'nesting_line_NNN_15_Ef_1.85.out')
#cal_nesting_factor_lineq('FS_NN15_Ef_1.80.bxsf', band, sigma, 'lineq.dat', 'nesting_line_NNN_15_Ef_1.80.out')
#cal_nesting_factor_lineq('FS_NN15_Ef_1.75.bxsf', band, sigma, 'lineq.dat', 'nesting_line_NNN_15_Ef_1.75.out')
#cal_nesting_factor_lineq('FS_NN15_Ef_1.70.bxsf', band, sigma, 'lineq.dat', 'nesting_line_NNN_15_Ef_1.70.out')

#cal_nesting_factor_lineq('FS_NN05_Ef_1.65.bxsf', band, sigma, 'lineq.dat', 'nesting_line_NNN_05_Ef_1.65.out')
#cal_nesting_factor_lineq('FS_NN05_Ef_1.60.bxsf', band, sigma, 'lineq.dat', 'nesting_line_NNN_05_Ef_1.60.out')
#cal_nesting_factor_lineq('FS_NN05_Ef_1.55.bxsf', band, sigma, 'lineq.dat', 'nesting_line_NNN_05_Ef_1.55.out')
#cal_nesting_factor_lineq('FS_NN05_Ef_1.50.bxsf', band, sigma, 'lineq.dat', 'nesting_line_NNN_05_Ef_1.50.out')

#cal_nesting_factor_areaq('FS_NN20_Ef_2.20.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_20_Ef_2.20')
#cal_nesting_factor_areaq('FS_NN20_Ef_2.15.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_20_Ef_2.15')
#cal_nesting_factor_areaq('FS_NN20_Ef_2.10.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_20_Ef_2.10')
#cal_nesting_factor_areaq('FS_NN20_Ef_2.05.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_20_Ef_2.05')
#cal_nesting_factor_areaq('FS_NN20_Ef_2.00.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_20_Ef_2.00')
#cal_nesting_factor_areaq('FS_NN20_Ef_1.95.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_20_Ef_1.95_sigma001')
#cal_nesting_factor_areaq('FS_NN20_Ef_1.90.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_20_Ef_1.90_sigma001 ')
#cal_nesting_factor_areaq('FS_NN20_Ef_1.85.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_20_Ef_1.85')
#cal_nesting_factor_areaq('FS_NN20_Ef_1.80.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_20_Ef_1.80')
#cal_nesting_factor_areaq('FS_NN20_Ef_1.75.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_20_Ef_1.75')
#cal_nesting_factor_areaq('FS_NN20_Ef_1.70.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_20_Ef_1.70')

#cal_nesting_factor_areaq('FS_NN15_Ef_2.20.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_15_Ef_2.20')
#cal_nesting_factor_areaq('FS_NN15_Ef_2.15.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_15_Ef_2.15')
#cal_nesting_factor_areaq('FS_NN15_Ef_2.10.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_15_Ef_2.10')
#cal_nesting_factor_areaq('FS_NN15_Ef_2.05.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_15_Ef_2.05')
#cal_nesting_factor_areaq('FS_NN15_Ef_2.00.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_15_Ef_2.00')
#cal_nesting_factor_areaq('FS_NN15_Ef_1.95.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_15_Ef_1.95')
#cal_nesting_factor_areaq('FS_NN15_Ef_1.90.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_15_Ef_1.90')
#cal_nesting_factor_areaq('FS_NN15_Ef_1.85.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_15_Ef_1.85')
#cal_nesting_factor_areaq('FS_NN15_Ef_1.80.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_15_Ef_1.80')
#cal_nesting_factor_areaq('FS_NN15_Ef_1.75.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_15_Ef_1.75')
#cal_nesting_factor_areaq('FS_NN15_Ef_1.70.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_15_Ef_1.70')
#
#cal_nesting_factor_areaq('FS_NN10_Ef_2.20.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_10_Ef_2.20')
#cal_nesting_factor_areaq('FS_NN10_Ef_2.15.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_10_Ef_2.15')
#cal_nesting_factor_areaq('FS_NN10_Ef_2.10.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_10_Ef_2.10')
#cal_nesting_factor_areaq('FS_NN10_Ef_2.05.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_10_Ef_2.05')
#cal_nesting_factor_areaq('FS_NN10_Ef_2.00.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_10_Ef_2.00')
#cal_nesting_factor_areaq('FS_NN10_Ef_1.95.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_10_Ef_1.95')
#cal_nesting_factor_areaq('FS_NN10_Ef_1.90.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_10_Ef_1.90')
#cal_nesting_factor_areaq('FS_NN10_Ef_1.85.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_10_Ef_1.85')
#cal_nesting_factor_areaq('FS_NN10_Ef_1.80.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_10_Ef_1.80')
#cal_nesting_factor_areaq('FS_NN10_Ef_1.75.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_10_Ef_1.75')
#cal_nesting_factor_areaq('FS_NN10_Ef_1.70.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_10_Ef_1.70')
#
#cal_nesting_factor_areaq('FS_NN05_Ef_2.20.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_05_Ef_2.20')
#cal_nesting_factor_areaq('FS_NN05_Ef_2.15.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_05_Ef_2.15')
#cal_nesting_factor_areaq('FS_NN05_Ef_2.10.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_05_Ef_2.10')
#cal_nesting_factor_areaq('FS_NN05_Ef_2.05.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_05_Ef_2.05')
#cal_nesting_factor_areaq('FS_NN05_Ef_2.00.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_05_Ef_2.00')
#cal_nesting_factor_areaq('FS_NN05_Ef_1.95.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_05_Ef_1.95')
#cal_nesting_factor_areaq('FS_NN05_Ef_1.90.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_05_Ef_1.90')
#cal_nesting_factor_areaq('FS_NN05_Ef_1.85.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_05_Ef_1.85')
#cal_nesting_factor_areaq('FS_NN05_Ef_1.80.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_05_Ef_1.80')
#cal_nesting_factor_areaq('FS_NN05_Ef_1.75.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_05_Ef_1.75')
#cal_nesting_factor_areaq('FS_NN05_Ef_1.70.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'nesting_area_NNN_05_Ef_1.70')



#cal_nesting_factor_areaq('B_fs_no_s18_spin.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'h0dope_s18_spin')
#cal_nesting_factor_areaq('B_fs_no_s18_nospin.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'h0dope_s18_nospin')
#cal_nesting_factor_areaq('B_fs_h1_s18_nospin.bxsf', band, sigma, 'h1dope_spin_kpoint.npy', 'h1dope_s18_nospin')
#cal_nesting_factor_areaq('B_fs_h2_s18_nospin.bxsf', band, sigma, 'h2dope_spin_kpoint.npy', 'h2dope_s18_nospin')
#cal_nesting_factor_areaq('B_fs_h3_s18_nospin.bxsf', band, sigma, 'h3dope_spin_kpoint.npy', 'h3dope_s18_nospin')


#cal_nesting_factor_areaq('B_fsup_no.bxsf', band, sigma, 'nodope_spin_kpoint.npy', 'h0dope_spin')
#cal_nesting_factor_areaq('B_fsup_h1.bxsf', band, sigma, 'h1dope_spin_kpoint.npy', 'h1dope_spin')
#cal_nesting_factor_areaq('B_fsup_h2.bxsf', band, sigma, 'h2dope_spin_kpoint.npy', 'h2dope_spin')
#cal_nesting_factor_areaq('B_fsup_h3.bxsf', band, sigma, 'h3dope_spin_kpoint.npy', 'h3dope_spin')

#print cal_nesting_factor_singleq('B_fsup_no_s22.bxsf', band, sigma, [0.333,0.333,0.0])
#print cal_nesting_factor_singleq('B_fsup_h1_s22.bxsf', band, sigma, [0.333,0.333,0.0])
#print cal_nesting_factor_singleq('B_fsup_h2_s22.bxsf', band, sigma, [0.333,0.333,0.0])
#print cal_nesting_factor_singleq('B_fsup_h3_s22.bxsf', band, sigma, [0.333,0.333,0.0])
#print cal_nesting_factor_singleq('B_fsup_h4_s22.bxsf', band, sigma, [0.333,0.333,0.0])
#print cal_nesting_factor_singleq('B_fsup_h5_s22.bxsf', band, sigma, [0.333,0.333,0.0])
#print cal_nesting_factor_singleq('B_fsup_h6_s22.bxsf', 3, sigma, [0.333,0.333,0.0])
#print cal_nesting_factor_singleq('B_fsup_h4.bxsf', band, sigma, [0.333,0.333,0.0])
#print cal_nesting_factor_singleq('B_fsup_h5.bxsf', band, sigma, [0.333,0.333,0.0])
#print cal_nesting_factor_singleq('B_fsup_h6.bxsf', band, sigma, [0.333,0.333,0.0])
#print cal_nesting_factor_singleq('B_fsup_h7.bxsf', band, sigma, [0.333,0.333,0.0])
#print cal_nesting_factor_singleq('B_fsup_h8.bxsf', band, sigma, [0.333,0.333,0.0])
#print cal_nesting_factor_singleq('B_fsup_h9.bxsf', band, sigma, [0.333,0.333,0.0])
#print ' '
#print cal_nesting_factor_singleq('B_fs_no.bxsf'  , band, sigma, [0.333,0.333,0.0])
#print cal_nesting_factor_singleq('B_fs_h1.bxsf'  , band, sigma, [0.333,0.333,0.0])
#print cal_nesting_factor_singleq('B_fs_h2.bxsf'  , band, sigma, [0.333,0.333,0.0])
#print cal_nesting_factor_singleq('B_fs_h3.bxsf'  , band, sigma, [0.333,0.333,0.0])
#print cal_nesting_factor_singleq('B_fs_h4.bxsf'  , band, sigma, [0.333,0.333,0.0])
#print cal_nesting_factor_singleq('B_fs_h5.bxsf'  , band, sigma, [0.333,0.333,0.0])
#print cal_nesting_factor_singleq('B_fs_h6.bxsf'  , band, sigma, [0.333,0.333,0.0])
#print cal_nesting_factor_singleq('B_fs_h7.bxsf'  , band, sigma, [0.333,0.333,0.0])
#print cal_nesting_factor_singleq('B_fs_h8.bxsf'  , band, sigma, [0.333,0.333,0.0])
#print cal_nesting_factor_singleq('B_fs_h9.bxsf'  , band, sigma, [0.333,0.333,0.0])

#print cal_nesting_factor_singleq('B_fsup_no_s24.bxsf', band, sigma, [0.333,0.333,0.0])
#print cal_nesting_factor_singleq('B_fsup_h1_s24.bxsf', band, sigma, [0.333,0.333,0.0])
#print cal_nesting_factor_singleq('B_fsup_h2_s24.bxsf', band, sigma, [0.333,0.333,0.0])
#print cal_nesting_factor_singleq('B_fsup_h3_s24.bxsf', band, sigma, [0.333,0.333,0.0])
#print cal_nesting_factor_singleq('B_fsup_h4_s24.bxsf', band, sigma, [0.333,0.333,0.0])
#print cal_nesting_factor_singleq('B_fsup_h5_s24.bxsf', band, sigma, [0.333,0.333,0.0])
#print cal_nesting_factor_singleq('B_fsup_h6_s24.bxsf', 2, sigma, [0.333,0.333,0.0])

#cal_nesting_factor_lineq('B_fsup_no.bxsf', band, sigma, 'lineq.dat', 'nesting_nodope_spin.out')
#cal_nesting_factor_lineq('B_fsup_h1.bxsf', band, sigma, 'lineq.dat', 'nesting_h1dope_spin.out')
#cal_nesting_factor_lineq('B_fsup_h2.bxsf', band, sigma, 'lineq.dat', 'nesting_h2dope_spin.out')
#cal_nesting_factor_lineq('B_fsup_h3.bxsf', band, sigma, 'lineq.dat', 'nesting_h3dope_spin.out')
#cal_nesting_factor_lineq('B_fsup_h4.bxsf', band, sigma, 'lineq.dat', 'nesting_nodope_spin.out')
#cal_nesting_factor_lineq('B_fsup_h5.bxsf', band, sigma, 'lineq.dat', 'nesting_h1dope_spin.out')
#cal_nesting_factor_lineq('B_fsup_h6.bxsf', band, sigma, 'lineq.dat', 'nesting_h2dope_spin.out')
#cal_nesting_factor_lineq('B_fsup_h7.bxsf', band, sigma, 'lineq.dat', 'nesting_h3dope_spin.out')
#cal_nesting_factor_lineq('B_fsup_h8.bxsf', band, sigma, 'lineq.dat', 'nesting_h2dope_spin.out')
#cal_nesting_factor_lineq('B_fsup_h9.bxsf', band, sigma, 'lineq.dat', 'nesting_h3dope_spin.out')
#
#
#cal_nesting_factor_lineq('B_fs_no.bxsf'  , band, sigma, 'lineq.dat', 'nesting_nodope_nospin.out')
#cal_nesting_factor_lineq('B_fs_h1.bxsf'  , band, sigma, 'lineq.dat', 'nesting_h1dope_nospin.out')
#cal_nesting_factor_lineq('B_fs_h2.bxsf'  , band, sigma, 'lineq.dat', 'nesting_h2dope_nospin.out')
#cal_nesting_factor_lineq('B_fs_h3.bxsf'  , band, sigma, 'lineq.dat', 'nesting_h3dope_nospin.out')
#cal_nesting_factor_lineq('B_fs_h4.bxsf'  , band, sigma, 'lineq.dat', 'nesting_h1dope_nospin.out')
#cal_nesting_factor_lineq('B_fs_h5.bxsf'  , band, sigma, 'lineq.dat', 'nesting_h2dope_nospin.out')
#cal_nesting_factor_lineq('B_fs_h6.bxsf'  , band, sigma, 'lineq.dat', 'nesting_h3dope_nospin.out')
#cal_nesting_factor_lineq('B_fs_h7.bxsf'  , band, sigma, 'lineq.dat', 'nesting_h1dope_nospin.out')
#cal_nesting_factor_lineq('B_fs_h8.bxsf'  , band, sigma, 'lineq.dat', 'nesting_h2dope_nospin.out')
#cal_nesting_factor_lineq('B_fs_h9.bxsf'  , band, sigma, 'lineq.dat', 'nesting_h3dope_nospin.out')#