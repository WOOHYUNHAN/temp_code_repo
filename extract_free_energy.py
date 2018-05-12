import os.path
import time
from math import *
from math import sqrt
import numpy as np
#from numeric import *




def extract_free_energy(input_name):
        f = open('thermal_properties.yaml','r')
        g = open(input_name,'w')
        fbuffer = f.readlines()
        atom = fbuffer[8].split()
        atom = float(atom[1])
        for i in range(len(fbuffer)):
                test1 = fbuffer[i].split()
                #print test1
                #print len(test1)
                if len(test1) == 0:
                       pass
                else:
                       if test1[0] == 'free_energy:':
                               if test1[1] == 'kJ/mol':   # 1 J/K = 0.0000103641 eV
                                       pass
                               else: 
                                       #print 'yes'
                                       fren = (float(test1[1]) * 0.0000103641 * 1000)/atom
                                       test2 = fbuffer[i-1].split()
                                       temper = float(test2[2])
                                       test3 = fbuffer[i+2].split()
                                       heat = (float(test3[1])*0.0000103641)/atom
                                       test4 = fbuffer[i+1].split()
                                       entro = (float(test4[1])*0.0000103641)/atom
                                       #print fren, temper
                                       line = str(temper)+'\t'+str(entro)+'\t'+str(fren)+'\t'+str(heat)+'\n'
                                       g.write(line)
        f.close()
        g.close()
        
extract_free_energy('557gamma')

