#!/usr/bin/python
import numpy as np
import sys, os
from os.path import exists
import time
#from  gen_sym_pos mk_sym_pos, writePOSCAR


def delete_job(ID_list):
    if len(ID_list) == 0:
	line = 'nothing'
    else:
    	line = 'qdel '
    	for i in range(len(ID_list)):
        	line += str(ID_list[i]) + ' '
    		#os.popen(line)
        	#print line
    	os.popen(line)
    return line

def find_job_list(job_name, sleep_time):
    cmd = 'qstat > l1'
    os.popen(cmd)
    time.sleep(2)
    f = open('l1')
    f_buffer = f.readlines()
    f.close()

    initial_list = []
    for i in range(len(f_buffer)):
        ff = f_buffer[i].split()
        name = ff[1]
        ID = ff[0]
        if name == job_name:
            initial_list.append(ID)

    time.sleep(sleep_time)

    cmd = 'qstat > l2'
    os.popen(cmd)
    time.sleep(2)
    h = open('l2')
    h_buffer = h.readlines()
    h.close()

    final_list = []
    for i in range(len(h_buffer)):
        hh = h_buffer[i].split()
        name = hh[1]
        ID = hh[0]
        if name == job_name:
            final_list.append(ID)

    initial = np.array(initial_list)
    final = np.array(final_list)

    common = np.intersect1d(initial, final)

    return common

def main(job_name, sleep_time):
    while True:
        print 'program is finding stuck jobs'
        job_list = find_job_list(job_name, sleep_time)
        line = delete_job(job_list)
        print 'Program deletes '+ str(line)
        time.sleep(10)

main('version3_GF_Ge', 270)
#a = find_job_list('Selected_GF', 100)
#print a

#b = delete_job([])
#print b
