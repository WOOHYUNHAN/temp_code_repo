import os.path
import time
from math import *
from math import sqrt
#from numeric import *


def phonopy_setting(start,end):
    a = int(start)
    b = int(end)
    
    for i in range(a,b+1):
        if i < 10:
            filenumber = str(i)
            os.system("mkdir "+'disp-'+filenumber)
            print filenumber
            #os.system("cp " + 'disp-'+filenumber+'/vasprun.xml' + ' '+filenumber)
            os.system("cp INCAR "+'disp-'+filenumber)
            os.system("cp KPOINTS "+'disp-'+filenumber)
            os.system("cp POTCAR "+'disp-'+filenumber)
            os.system("cp sh "+'disp-'+filenumber)
            os.system("cp POSCAR-00"+str(i)+" "+'disp-'+filenumber+'/POSCAR')
            os.system("ln WAVECAR"+" "+'disp-'+filenumber+'/.')
        elif i < 100:
            filenumber = str(i)
            os.system("mkdir "+'disp-'+filenumber)
            print filenumber
            #os.system("cp " + 'disp-'+filenumber+'/vasprun.xml' + ' '+filenumber)
            os.system("cp INCAR "+'disp-'+filenumber)
            os.system("cp KPOINTS "+'disp-'+filenumber)
            os.system("cp POTCAR "+'disp-'+filenumber)
            os.system("cp sh "+'disp-'+filenumber)
            os.system("cp POSCAR-0"+str(i)+" "+'disp-'+filenumber+'/POSCAR')
            os.system("ln WAVECAR"+" "+'disp-'+filenumber+'/.')
        else:
            filenumber = str(i)
            os.system("mkdir "+'disp-'+filenumber)
            print filenumber
            #os.system("cp " + 'disp-'+filenumber+'/vasprun.xml' + ' '+filenumber)
            os.system("cp INCAR "+'disp-'+filenumber)
            os.system("cp KPOINTS "+'disp-'+filenumber)
            os.system("cp POTCAR "+'disp-'+filenumber)
            os.system("cp sh "+'disp-'+filenumber)
            os.system("cp POSCAR-"+str(i)+" "+'disp-'+filenumber+'/POSCAR')
            os.system("ln WAVECAR"+" "+'disp-'+filenumber+'/.')

    #os.system("rm PARCHG.0"+filenumber+".0001")

    
def submit_pbs_job(f_start, f_end):
    f=os.popen('pwd').readlines()
    #os.system('pwd')
    filehead='disp-'
    for i in range(f_end-f_start+1):
        filename= filehead+str(i+f_start)
        dir_name=f[0][:-1]+'/'+str(filename)
        #print dir_name
        cmd = 'cd ' + dir_name + ';qsub sh'
        os.system(cmd)
        print str(i+1)+'/'+str(f_end-f_start+1)
        time.sleep(5)

phonopy_setting(1, 8)
submit_pbs_job(1, 8)
