import os.path
import time
from math import *
from math import sqrt
from numpy.linalg import norm
import numpy as np
import matplotlib.pyplot as plt
#from numeric import *

def extract_electrostatic(filename, NIONS):
    f = open(filename)
    fbuffer = f.readlines()
    f.close()
    for i in range(len(fbuffer)):
        #print fbuffer[i].split()
        if len(fbuffer[i].split()) == 0:
            #print 'oh'
            pass
        else:
            if fbuffer[i].split()[0] == 'average' and fbuffer[i].split()[1] == '(electrostatic)':
                #print 'oh'
                line = i
    num_list = []
    potential_list = []

    add = (NIONS/5) + 1

    for i in range(line+3, line+3+add):
        temp = fbuffer[i].split()
        for j in range(len(temp)/2):
            num_list.append(int(temp[2*j]))
            potential_list.append(float(temp[2*j+1]))


    print num_list
    print potential_list

    g = open('electrostatic_pot', 'w')
    for i in range(len(num_list)):
        line = str(num_list[i]) + '\t' + str(potential_list[i]) + '\n'
        g.write(line)
    g.close()



extract_electrostatic('OUTCAR',28)