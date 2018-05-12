'''
Created on 2014.  05.  14.
 
@author: nwan
@author: hanwooh
'''
import numpy as np 
import sys
import os, re
from os.path import isfile, exists


def get_info(dir_name, file_name):
    # get total energy from line of outcar
    '  free  energy   TOTEN  =       -53.472728 eV'
    '''
    FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)
    ---------------------------------------------------
    free  energy   TOTEN  =      -178.22800012 eV

    energy  without entropy=     -178.22769568  energy(sigma->0) =     -178.22789864
    enthalpy is  TOTEN    =      -153.45223930 eV   P V=       24.77576082
    '''
    f = open(str(dir_name)+'/'+str(file_name), 'r')
    outcar = f.readlines()
    f.close()
    reg = '(?<=free  energy   TOTEN)\s*=\s+[-+]?[0-9]*\.?[0-9]+'
    find = []
    for line in outcar:
        m = re.findall(reg, line)
        find += m
    tot_E = float(find[-1].replace('=',''))
    file = str(dir_name)+'/'+str(file_name)
    cmd1 = 'grep vol ' + file
    #print cmd1
    f = os.popen(cmd1)
    output = f.readlines()
    vol = float(output[-1].split()[-1])
    final_E = (tot_E / 28 )/ 27.211
    final_V = (vol / 28) * 6.7483344939499953831306552639989
    outline = str(final_V) + ' ' + str(final_E)
    return outline

print get_info('.', 'OUTCAR_n5gpa')
print get_info('.', 'OUTCAR_n3gpa')
print get_info('.', 'OUTCAR_0gpa')
print get_info('.', 'OUTCAR_3gpa')
print get_info('.', 'OUTCAR_5gpa')
print get_info('.', 'OUTCAR_10gpa')
print get_info('.', 'OUTCAR_15gpa')

