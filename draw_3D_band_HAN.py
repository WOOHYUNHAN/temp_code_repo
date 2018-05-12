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

def Read_kpoints_for_3Dband(input_file='KPOITNS'):
    f = open(input_file)
    fbuffer = f.readlines()
    f.close()

    nk = int(fbuffer[1].split()[0]) 

    kx = []
    ky = []

    for i in range(nk):
        tempf = fbuffer[i+3].split()
        kx.append(float(tempf[0]))
        ky.append(float(tempf[1]))

    kx = np.array(kx)
    ky = np.array(ky)

    return kx, ky

def Read_eigenvalue_for_3Dband(input_file, Ef):
    f = open(input_file)
    fbuffer = f.readlines()
    f.close()

    info = np.array([int(fbuffer[5].split()[i]) for i in range(len(fbuffer[5].split()))])

    x, nk, nb = int(info[0]), int(info[1]), int(info[2])

    EIGEN_array = np.zeros((nb,nk))

    #print EIGEN_array[23][2400]

    for i in range(nk):
        for j in range(nb):
            k_position = 7 + (i*(nb+2))
            band_position = j+1
            value = float(fbuffer[k_position+band_position].split()[1])
            EIGEN_array[j,i] = value - Ef

    return EIGEN_array

def Read_eigenvalue_for_3Dband_wannier90(nx, ny, Ef, num_wann):
    f = open('wannier90-kslice-bands.dat')
    fbuffer = f.readlines()
    f.close()

    g = open('wannier90-kslice-coord.dat')
    gbuffer = g.readlines()
    g.close()

    num_kpoints = nx*ny

    kx = []
    ky = []

    for i in range(nx):
      kx.append(float(gbuffer[i*ny].split()[0]))

    for i in range(ny):
      ky.append(float(gbuffer[i].split()[1]))




    bands=np.loadtxt('wannier90-kslice-bands.dat')
    EIGEN_array=bands.reshape((num_kpoints,num_wann)).T

    #EIGEN_array = np.zeros((num_wann,num_kpoints))
    #for i in range(num_kpoints):
    #  for j in range(num_wann):
    #    ftemp = float(fbuffer[num_kpoints*j + i].split()[0]) - Ef
    #    EIGEN_array[j,i] = ftemp

    return np.array(kx), np.array(ky), np.array(EIGEN_array)-Ef








from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib import rc
from scipy import interpolate

params = {'backend': 'ps',
      'text.latex.preamble': [r"\usepackage{upgreek}",
                              r"\usepackage{siunitx}",
                              r"\usepackage{amsmath}",
                              r"\usepackage{amstext}",],
      'axes.labelsize': 20,
      'axes.linewidth' : 2,
      'axes.titlesize' : 16,

      'lines.linewidth': 2,
      'lines.markersize' : 8,
      'lines.markeredgewidth' : 2,

      'text.fontsize':100,
      'legend.fontsize': 15,

      'axes.labelpad': 100,

      'xtick.major.width'    : 2,    # major tick width in points
      'xtick.minor.width'    : 1,    # minor tick width in points
      'ytick.major.width'    : 2,    # major tick width in points
      'ytick.minor.width'    : 1,    # minor tick width in points

      'xtick.major.pad'      : 25,      # distance to major tick label in points
      'xtick.minor.pad'      : 25,      # distance to the minor tick label in points
      'ytick.major.pad'      : 8,     # distance to major tick label in points
      'ytick.minor.pad'      : 8,     # distance to the minor tick label in points
      
      'xtick.labelsize': 16,
      'ytick.labelsize': 16,

      'figure.figsize': [8,8],
      'figure.dpi': 800,
      # 'text.usetex': True,
      'axes.unicode_minus': True,
      'ps.usedistiller' : 'xpdf'
      }     
rc('font',**{'family':'serif','serif':['Helvetica']})

fig = plt.figure(figsize=(8,8))
ax = fig.gca(projection='3d')


#kx, ky = Read_kpoints_for_3Dband('KPOINTS')

#print len(kx), len(ky)

#EIGEN_array = Read_eigenvalue_for_3Dband('EIGENVAL',-2.1097)


# Make data.
#nx, ny = (80, 80)
#x = np.linspace(-0.16, 0.16, nx)
#y = np.linspace(-0.10, 0.10, ny)
#xv, yv = np.meshgrid(x, y)

#print np.min(EIGEN_array[18])
#print np.max(EIGEN_array[17])


###############
###WANNIER90###
###############



nx, ny = (200, 200)
kx, ky, EIGEN_array = Read_eigenvalue_for_3Dband_wannier90(nx, ny, -1.20, 16)
x = np.linspace(np.min(kx), np.max(kx), nx)
y = np.linspace(np.min(ky), np.max(ky), ny)
xv, yv = np.meshgrid(x, y)

#print len(EIGEN_array[5])

a = EIGEN_array[5].reshape((nx,ny))
#af = interpolate.interp2d(xv, yv, a, kind='linear')
b = EIGEN_array[6].reshape((nx,ny))

#print np.max(a), np.min(b)
#print kx
surf = ax.plot_surface(xv,yv,a, cmap=cm.coolwarm, linewidth=1, antialiased=False, vmin=-1.0, vmax=1.0)
surf = ax.plot_surface(xv,yv,b, cmap=cm.coolwarm, linewidth=1, antialiased=False, vmin=-1.0, vmax=1.0)

#ax.scatter(kx, ky, EIGEN_array[5], s=10, c= '#2ca25f', edgecolor=(0,0,0,0))
#ax.scatter(kx, ky, EIGEN_array[6], s=10, c= '#2ca25f', edgecolor=(0,0,0,0))
#ax.scatter(kx, ky, EIGEN_array[20], s=10, c= '#08519C', edgecolor=(0,0,0,0))
#ax.scatter(kx, ky, EIGEN_array[41], s=10, c= '#08519C', edgecolor=(0,0,0,0))
#cbar=ax.contourf(xv, yv, b-a, zdir='z', offset=-1,cmap=cm.coolwarm, min=0.0, max=1.0)
#fig.colorbar(cbar, shrink=0.5, aspect=5)

#a = np.arange(6).reshape((2,3)).T
#print a
# Customize the z axis.
#ax.set_zlim(-1.0, 1.0)
#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

ax.set_xlabel(r'$k_x$', labelpad=10)
ax.set_ylabel(r'$k_y$')
ax.set_zlabel(r'E-$E_F$ (energy)')
#ax.w_xaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
#ax.w_yaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
#ax.w_zaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
#plt.rcParams['xtick.major.pad'] = 10
plt.tick_params(
axis='x',          # changes apply to the x-axis
which='both',      # both major and minor ticks are affected
bottom='on',      # ticks along the bottom edge are off
top='off',         # ticks along the top edge are off
labelbottom='on') # labels along the bottom edge are o
# Add a color bar which maps values to colors.
#fig.colorbar(surf, shrink=0.5, aspect=5)

plt.tick_params(
axis='y',          # changes apply to the x-axis
which='both',      # both major and minor ticks are affected
bottom='on',      # ticks along the bottom edge are off
top='off',         # ticks along the top edge are off
labelleft='on') # labels along the bottom edge are o
# Add a color bar which maps values to colors.
#fig.colorbar(surf, shrink=0.5, aspect=5)

plt.tick_params(
axis='z',          # changes apply to the x-axis
which='both',      # both major and minor ticks are affected
bottom='off',      # ticks along the bottom edge are off
top='off',         # ticks along the top edge are off
labelbottom='on') # labels along the bottom edge are o
# Add a color bar which maps values to colors.
#fig.colorbar(surf, shrink=0.5, aspect=5)

#fig.savefig('output.eps')
#fig.savefig('output.png')
#fig.savefig('output.pdf')

plt.show()