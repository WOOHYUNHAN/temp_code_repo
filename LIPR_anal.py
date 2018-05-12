import os.path
import time
from math import *
from math import sqrt
from numpy.linalg import norm
import numpy as np
import matplotlib.pyplot as plt
#from numeric import *
        
def get_projected_dos(kpoint_num, band_num, atom_num, orbital, filename='PROCAR'):
    f = open(filename,'r')
    temp_f = f.readlines()
    band_list = []

    for i in range(len(temp_f)):
        tempf = temp_f[i].split()
        if len(tempf) == 0:
            pass
        else:
            if tempf[0] == 'k-point':
                if tempf[1] == str(kpoint_num):
                    temp_line = i
                    break
    #print temp_line
    specific_line = temp_line + 2 + (band_num-1)*(2+85+2+170+2) + 2 + atom_num
    atom_info = np.array(temp_f[specific_line].split())

    #print specific_line
    #print atom_info
    #print np.array(atom_info[1:3],dtype=float)

    if orbital == 's':
        output = np.sum(np.array(atom_info[1:2],dtype=float))
    elif orbital == 'p':
        output = np.sum(np.array(atom_info[2:5],dtype=float))
    elif orbital == 'd':
        output = np.sum(np.array(atom_info[5:9],dtype=float))
    else:
        print 'This is error'

    return output

def get_data(kpoint_num, band_num, start_atom, end_atom, orbital, filename='PROCAR'):
    tot_atom = end_atom - start_atom + 1
    out = []

    for i in range(start_atom, end_atom+1):
        temp = get_projected_dos(kpoint_num, band_num, i, orbital)
        print str(i)
        out.append(temp)

    out = np.array(out)
    out = out / np.sum(out)

    output_name = 'k_'+str(kpoint_num)+'_band_'+str(band_num)+'_from_'+str(start_atom)+'_to_'+str(end_atom)+'_'+str(orbital)+'orbital.out'
    g = open(output_name,'w')
    g.write('atom'+'\t'+'value'+'\n')

    for i in range(tot_atom):
        index = start_atom + i
        line = str(i+1) + '\t' + str(out[i]) + '\n'
        g.write(line)

    g.close()
    return out




#get_projected_dos(3, 355, 85, 'p')
get_data(1, 360, 37, 85, 'p', filename='PROCAR')
get_data(1, 368, 37, 85, 'p', filename='PROCAR')
get_data(1, 369, 37, 85, 'p', filename='PROCAR')
get_data(1, 375, 37, 85, 'p', filename='PROCAR')