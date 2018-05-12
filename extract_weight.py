import os.path
import time
from math import *
from math import sqrt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource, Normalize
import matplotlib as mpl
#from numeric import *



def extract_weight(input_name, output_name):
        f = open(input_name,'r')
        g = open(output_name,'w')
        fbuffer = f.readlines()
        kpoint_array = []
        weight_array = []
        for i in range(len(fbuffer)):
                test1 = fbuffer[i].split()
                #print test1
                #print len(test1)
                if len(test1) == 0:
                      pass
                else:
                      if test1[0] == 'nqs=':
                            weight_array.append(float(test1[1]))
                            test2 = fbuffer[i+1].split()
                            kpoint_array.append([float(test2[1]),float(test2[2]),float(test2[3])])
                      else:
                            pass
        #print kpoint_array
        g.write(str(len(weight_array))+'\n')
        line = ''
        for i in range(len(weight_array)):
              line += str(kpoint_array[i][0]) + ' '+ str(kpoint_array[i][1]) + ' '+str(kpoint_array[i][2]) + ' '+str(weight_array[i]) + '\n'
        g.write(line)
        for i in range(len(weight_array)):
            line = 'elph_dir/elph.inp_lambda.'+str(i+1)+'\n'
            g.write(line)
        f.close()
        g.close()
        print len(weight_array)
        print len(kpoint_array)

def extract_linewidth(filename, nbnd, nks):
    f = open(filename, 'r')
    tempf = f.readlines()
    f.close()

    dataperline = 6

    linec = nbnd / dataperline
    liner = nbnd % dataperline

    kpoint_array = []
    linewidth_array =  []

    for i in range(nks):
        startline = (linec + 1 + 1) * i + 1
        #print startline
        #print float(tempf[startline].split()[2])
        temp_k = [float(tempf[startline].split()[z]) for z in range(3)]
        kpoint_array.append(temp_k)
        linewidth_array.append([])
        for j in range(linec + 1):
            if j == linec:
                for k in range(liner):
                    linewidth_array[i].append((float(tempf[startline+j+1][10*k:10*(k+1)])))
            else:
                for k in range(dataperline):
                    linewidth_array[i].append((float(tempf[startline+j+1][10*k:10*(k+1)])))



    #print kpoint_array[7]
    #print linewidth_array[7]
    return np.array(kpoint_array), np.array(linewidth_array)

def cal_lambda(fliename1, fliename2, nbnd, nks, dos, outdata_name):
    thz2ry = 0.00030428427852
    cm2thz = 0.02998
    kpoint_array1, linewidth_array = extract_linewidth(fliename1, nbnd, nks) #unit = GHz
    kpoint_array2, freq_array = extract_linewidth(fliename2, nbnd, nks) # unit = THz

    lamb_mode = (linewidth_array * 0.001 * thz2ry) / (dos * pi) / (freq_array*cm2thz*thz2ry) / (freq_array*cm2thz*thz2ry)
    #lamb_mode = []
    #print freq_array[2]
    #for i in range(len(freq_array)):
    #    if freq_array[i] < 20: # cutoff = 20cm-1
    #        lamb_mode.append(0.0)
    #    else:
    #        temp = (linewidth_array[i] * 0.001 * thz2ry) / (dos * pi) / (freq_array[i]*cm2thz*thz2ry) / (freq_array[i]*cm2thz*thz2ry)
    #        lamb_mode.append(temp)
        #print linewidth_array
    #lamb_mode = np.array(lamb_mode)
    lamb = []
    for i in range(len(lamb_mode)):
        find_cutoff_freq = np.where(freq_array[i] < 20.0)[0]
        if len(find_cutoff_freq) != 0:
            #print zfind_zero_freq
            for j in range(len(find_cutoff_freq)):
                #print lamb_mode[i]
                lamb_mode[i][j] = 0.0
        temp = np.sum(lamb_mode[i])
        lamb.append(temp)
    lamb = np.array(lamb)

    #print kpoint_array1[5100]
    #print kpoint_array2[5100]
    #print linewidth_array[5100]
    #print freq_array[5100]
    #print lamb[5100]
    #print lamb[5101]
    print 'max lambda is ', max(lamb), 'min lambda is ', min(lamb)
    print 'All sum of lambda = ', np.sum(lamb)/len(freq_array)
    #print len(kpoint_array1)#.reshape((101,101))
    z =  lamb.reshape((101,101)).T
    #print lamb[101]
    #print kpoint_array1[3]
    #print z[1][0]


    cmap = plt.cm.gist_earth
    ls = LightSource(315, 45)
    rgb = ls.shade(z, cmap, vmin=0.0, vmax=7)
    fig, ax = plt.subplots(figsize=(8, 7))
    ax.imshow(rgb, interpolation='gaussian', aspect="auto")
    ax1 = fig.add_axes([0.05, 0.05, 0.9, 0.05])
    cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap, orientation='horizontal')
    plt.show()

    np.save(outdata_name+'_kpoint', kpoint_array1)
    np.save(outdata_name+'_lambda', lamb)


    return kpoint_array1, lamb


 

#extract_weight('B.q2r.out', 'kpoint1')
#extract_linewidth('elph.gamma.1', 27, 10000)
#extract_linewidth('B555.freq', 27, 10201)
cal_lambda('elph.gamma.1', 'B555.freq', 27, 10201, 49.669501, 'h1dope_nospin')

