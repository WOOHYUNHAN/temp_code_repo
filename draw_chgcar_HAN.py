import os.path
import time
from math import *
from math import sqrt
from numpy.linalg import norm
import numpy as np
import matplotlib.pyplot as plt


def linear_criteria(temp_c, pos1_c, pos2_c, h, sigma):
    temp = False
    logx = False
    logy = False


    tempx = (np.sum(np.square(temp_c - pos1_c))-h*h)
    tempy = (np.sum(np.square(temp_c - pos2_c))-h*h)


    if tempx <= np.sum(np.square(pos1_c-pos2_c)):
        logx = True

    if tempy <= np.sum(np.square(pos1_c-pos2_c)):
        logy = True


    if logx == True and logy == True:
        fh = h
        temp = True
        #print 'ok'
    else:
        minval = min(tempx, tempy) 
        #fh = sqrt(minval + h*h)
        fh = 10000000000000000000000000

    # gauss broadening

    f = exp(-(fh*fh)/(2*sigma*sigma))


    return f,temp


def triangular_criteria(temp_c, pos1_c, pos2_c, pos3_c, h, sigma):
    temp = False
    logx = False
    logy = False
    logz = False

    tempx = (np.sum(np.square(temp_c - pos1_c))-h*h)
    tempy = (np.sum(np.square(temp_c - pos2_c))-h*h)
    tempz = (np.sum(np.square(temp_c - pos3_c))-h*h)


    if tempx <= max(np.sum(np.square(pos2_c-pos1_c)),np.sum(np.square(pos3_c-pos1_c))):
        logx = True

    if tempy <= max(np.sum(np.square(pos3_c-pos2_c)),np.sum(np.square(pos1_c-pos2_c))):
        logy = True

    if tempz <= max(np.sum(np.square(pos2_c-pos3_c)),np.sum(np.square(pos1_c-pos3_c))):
        logz = True


    if logx == True and logy == True and logz == True:
        fh = h
        temp = True
        #print 'ok'
    else:
        minval = min(tempx, tempy, tempz) 
        #fh = sqrt(minval + h*h)
        fh = 10000000000000000000000000

    # gauss broadening

    f = exp(-(fh*fh)/(2*sigma*sigma))


    return f,temp



def sum_charge_density_two_point(filename, sigma, atom1, atom2):
    f = open(filename, 'r')#
    fbuffer = f.readlines()
    f.close()


    line = 5

    #print fbuffer[6].split()
    #print len(fbuffer[6])
    univ = float(fbuffer[1].split()[0])

    axis_a = np.array([float(fbuffer[2].split()[i]) for i in range(3)]) * univ
    axis_b = np.array([float(fbuffer[3].split()[i]) for i in range(3)]) * univ
    axis_c = np.array([float(fbuffer[4].split()[i]) for i in range(3)]) * univ

    trans = np.array([axis_a,axis_b,axis_c])
    trans = np.transpose(trans)

    atoms_array = np.array([int(fbuffer[6].split()[i]) for i in range(len(fbuffer[6].split()))])
    num_atoms = np.sum(atoms_array)

    position1 = np.array([float(fbuffer[7+atom1].split()[i]) for i in range(3)])
    position2 = np.array([float(fbuffer[7+atom2].split()[i]) for i in range(3)])

    center = (position1+position2)/2 -np.array([0.5, 0.5, 0.5])

    #find equation of line



    pos1_c = np.dot(trans, position1-center)
    pos2_c = np.dot(trans, position2-center)

    v = pos2_c - pos1_c # directional vector for the line


    ##  Eq: ax + by + cz = d

    grida = int(fbuffer[7+num_atoms+2].split()[0])
    gridb = int(fbuffer[7+num_atoms+2].split()[1])
    gridc = int(fbuffer[7+num_atoms+2].split()[2])

    deltax, deltay, deltaz = 1.0/grida, 1.0/gridb, 1.0/gridc

    totalgrid = grida * gridb * gridc

    charge = 0
    in_line = 0
    out_line = 0



    for i in range(gridc):
        for j in range(gridb):
            for k in range(grida):
                temp = [deltax*k, deltay*j, deltaz*i]
                temp_c = np.dot(trans, temp-center)
                switch = False
                h = np.linalg.norm(np.cross((temp_c-pos1_c),v)) / np.linalg.norm(v)
                factor, check = linear_criteria(temp_c, pos1_c, pos2_c, h, sigma)
                if check == True:
                    in_line += 1
                else:
                    out_line += 1
                number = (gridb * grida * i) + (grida * j) + (k)
                col = int((float(number/line) - int(number/line))*5)
                row = int(number/line)
                #print col
                temp_charge = float(fbuffer[row+7+num_atoms+2].split()[col])
                charge += temp_charge * factor

    final = charge / totalgrid
    total = in_line+out_line

    print str(sigma)+ '\t' + str(final) + '\t'

    return final, in_line, out_line, float(in_line/total)


def sum_charge_density_three_point(filename, sigma, atom1, atom2, atom3):
    f = open(filename, 'r')#
    fbuffer = f.readlines()
    f.close()


    line = 5

    #print fbuffer[6].split()
    #print len(fbuffer[6])
    univ = float(fbuffer[1].split()[0])

    axis_a = np.array([float(fbuffer[2].split()[i]) for i in range(3)]) * univ
    axis_b = np.array([float(fbuffer[3].split()[i]) for i in range(3)]) * univ
    axis_c = np.array([float(fbuffer[4].split()[i]) for i in range(3)]) * univ

    trans = np.array([axis_a,axis_b,axis_c])
    trans = np.transpose(trans)

    atoms_array = np.array([int(fbuffer[6].split()[i]) for i in range(len(fbuffer[6].split()))])
    num_atoms = np.sum(atoms_array)

    position1 = np.array([float(fbuffer[7+atom1].split()[i]) for i in range(3)])
    position2 = np.array([float(fbuffer[7+atom2].split()[i]) for i in range(3)])
    position3 = np.array([float(fbuffer[7+atom3].split()[i]) for i in range(3)])

    center = (position1+position2+position3)/3 -np.array([0.5, 0.5, 0.5])

    #find equation of plane



    pos1_c = np.dot(trans, position1-center)
    pos2_c = np.dot(trans, position2-center)
    pos3_c = np.dot(trans, position3-center)

    a, b, c = np.cross(pos3_c - pos1_c, pos2_c - pos1_c)
    d = np.dot(np.array([a,b,c]), pos3_c)

    #print a, b, c, -d , -d/c

    ##  Eq: ax + by + cz = d

    grida = int(fbuffer[7+num_atoms+2].split()[0])
    gridb = int(fbuffer[7+num_atoms+2].split()[1])
    gridc = int(fbuffer[7+num_atoms+2].split()[2])

    deltax, deltay, deltaz = 1.0/grida, 1.0/gridb, 1.0/gridc

    totalgrid = grida * gridb * gridc

    charge = 0
    in_tri = 0
    out_tri = 0



    for i in range(gridc):
        for j in range(gridb):
            for k in range(grida):
                temp = [deltax*k, deltay*j, deltaz*i]
                temp_c = np.dot(trans, temp-center)
                switch = False
                h = abs(np.dot([a,b,c], temp_c) - d) / sqrt(a*a + b*b + c*c)
                factor, check = triangular_criteria(temp_c, pos1_c, pos2_c, pos3_c, h, sigma)
                if check == True:
                    in_tri += 1
                else:
                    out_tri += 1
                number = (gridb * grida * i) + (grida * j) + (k)
                col = int((float(number/line) - int(number/line))*5)
                row = int(number/line)
                #print col
                temp_charge = float(fbuffer[row+7+num_atoms+2].split()[col])
                charge += temp_charge * factor

    final = charge / totalgrid
    total = in_tri+out_tri

    print str(sigma)+ '\t' + str(final) + '\t'

    return final, in_tri, out_tri, float(in_tri/total)

def cal_charge_density_one_point(filename, sigma, point):
    f = open(filename, 'r')#
    fbuffer = f.readlines()
    f.close()


    line = 5

    #print fbuffer[6].split()
    #print len(fbuffer[6])
    univ = float(fbuffer[1].split()[0])

    axis_a = np.array([float(fbuffer[2].split()[i]) for i in range(3)]) * univ
    axis_b = np.array([float(fbuffer[3].split()[i]) for i in range(3)]) * univ
    axis_c = np.array([float(fbuffer[4].split()[i]) for i in range(3)]) * univ

    vol = np.dot(np.cross(axis_a,axis_b),axis_c)
    #print vol
    #sigma = sqrt(np.linalg.norm(axis_a)**2 + np.linalg.norm(axis_b)**2 + np.linalg.norm(axis_c)**2) / 60
    #sigma = 6* sigma
    #print sigma


    trans = np.array([axis_a,axis_b,axis_c])
    trans = np.transpose(trans)

    atoms_array = np.array([int(fbuffer[6].split()[i]) for i in range(len(fbuffer[6].split()))])
    num_atoms = np.sum(atoms_array)


    #position1 = np.array([float(fbuffer[7+atom1].split()[i]) for i in range(3)])
    #position2 = np.array([float(fbuffer[7+atom2].split()[i]) for i in range(3)])
    #position3 = np.array([float(fbuffer[7+atom3].split()[i]) for i in range(3)])


    point = [point[0]-int(point[0]), point[1]-int(point[1]), point[2]-int(point[2])]
    point = np.array(point)
    #point = (position1+position2+position3)/3.0
    #point = (position1+position2)/2.0

    center = point-np.array([0.5, 0.5, 0.5])
    point_c = np.dot(trans, point-center)
    #print point, point-center, point_c

    grida = int(fbuffer[7+num_atoms+2].split()[0])
    gridb = int(fbuffer[7+num_atoms+2].split()[1])
    gridc = int(fbuffer[7+num_atoms+2].split()[2])

    deltax, deltay, deltaz = 1.0/grida, 1.0/gridb, 1.0/gridc

    totalgrid = grida * gridb * gridc

    charge = 0


    for i in range(gridc):
        for j in range(gridb):
            for k in range(grida):
                temp = [deltax*k, deltay*j, deltaz*i]
                temp_c = np.dot(trans, temp-center)
                h = np.linalg.norm(temp_c - point_c)
                #print h
                factor = 1 #= exp(-(h*h)/(2*sigma*sigma))
                if h < sigma:
                    number = (gridb * grida * i) + (grida * j) + (k)
                    col = int((float(number/line) - int(number/line))*5)
                    row = int(number/line)
                    #print col
                    temp_charge = float(fbuffer[row+7+num_atoms+2].split()[col])
                    charge += temp_charge * factor

    final = charge /totalgrid


    print str(sigma)+ '\t' + str(final/vol) + '\t'

    return final, vol

# inter 3c = 7 8 9
# inter 2c = 1 4
# intra 3c = 1 2 3
# intra 2c = 1 2


#print '0gpa'
#sum_charge_density_three_point('CHGCAR_0gpa', 0.1, 7,8,9)
#sum_charge_density_three_point('CHGCAR_0gpa', 0.2, 7,8,9)
#sum_charge_density_three_point('CHGCAR_0gpa', 0.3, 7,8,9)
#sum_charge_density_three_point('CHGCAR_0gpa', 0.4, 7,8,9)
#sum_charge_density_three_point('CHGCAR_0gpa', 0.5, 7,8,9)
#sum_charge_density_three_point('CHGCAR_0gpa', 0.6, 7,8,9)
#sum_charge_density_three_point('CHGCAR_0gpa', 0.7, 7,8,9)
#sum_charge_density_three_point('CHGCAR_0gpa', 0.8, 7,8,9)
#sum_charge_density_three_point('CHGCAR_0gpa', 0.9, 7,8,9)
#sum_charge_density_three_point('CHGCAR_0gpa', 1.1, 7,8,9)
#sum_charge_density_three_point('CHGCAR_0gpa', 1.2, 7,8,9)
#sum_charge_density_three_point('CHGCAR_0gpa', 1.3, 7,8,9)
#sum_charge_density_three_point('CHGCAR_0gpa', 1.4, 7,8,9)
#sum_charge_density_three_point('CHGCAR_0gpa', 1.5, 7,8,9)
#
#print '20gpa'
#sum_charge_density_three_point('CHGCAR_20gpa', 0.1,7,8,9)
#sum_charge_density_three_point('CHGCAR_20gpa', 0.2,7,8,9)
#sum_charge_density_three_point('CHGCAR_20gpa', 0.3,7,8,9)
#sum_charge_density_three_point('CHGCAR_20gpa', 0.4,7,8,9)
#sum_charge_density_three_point('CHGCAR_20gpa', 0.5,7,8,9)
#sum_charge_density_three_point('CHGCAR_20gpa', 0.6,7,8,9)
#sum_charge_density_three_point('CHGCAR_20gpa', 0.7,7,8,9)
#sum_charge_density_three_point('CHGCAR_20gpa', 0.8,7,8,9)
#sum_charge_density_three_point('CHGCAR_20gpa', 0.9,7,8,9)
#sum_charge_density_three_point('CHGCAR_20gpa', 1.1,7,8,9)
#sum_charge_density_three_point('CHGCAR_20gpa', 1.2,7,8,9)
#sum_charge_density_three_point('CHGCAR_20gpa', 1.3,7,8,9)
#sum_charge_density_three_point('CHGCAR_20gpa', 1.4,7,8,9)
#sum_charge_density_three_point('CHGCAR_20gpa', 1.5,7,8,9)
#
#print '40gpa'
#sum_charge_density_three_point('CHGCAR_40gpa', 0.1,7,8,9)
#sum_charge_density_three_point('CHGCAR_40gpa', 0.2,7,8,9)
#sum_charge_density_three_point('CHGCAR_40gpa', 0.3,7,8,9)
#sum_charge_density_three_point('CHGCAR_40gpa', 0.4,7,8,9)
#sum_charge_density_three_point('CHGCAR_40gpa', 0.5,7,8,9)
#sum_charge_density_three_point('CHGCAR_40gpa', 0.6,7,8,9)
#sum_charge_density_three_point('CHGCAR_40gpa', 0.7,7,8,9)
#sum_charge_density_three_point('CHGCAR_40gpa', 0.8,7,8,9)
#sum_charge_density_three_point('CHGCAR_40gpa', 0.9,7,8,9)
#sum_charge_density_three_point('CHGCAR_40gpa', 1.1,7,8,9)
#sum_charge_density_three_point('CHGCAR_40gpa', 1.2,7,8,9)
#sum_charge_density_three_point('CHGCAR_40gpa', 1.3,7,8,9)
#sum_charge_density_three_point('CHGCAR_40gpa', 1.4,7,8,9)
#sum_charge_density_three_point('CHGCAR_40gpa', 1.5,7,8,9)
#
#print '50gpa'
#sum_charge_density_three_point('CHGCAR_50gpa', 0.1,7,8,9)
#sum_charge_density_three_point('CHGCAR_50gpa', 0.2,7,8,9)
#sum_charge_density_three_point('CHGCAR_50gpa', 0.3,7,8,9)
#sum_charge_density_three_point('CHGCAR_50gpa', 0.4,7,8,9)
#sum_charge_density_three_point('CHGCAR_50gpa', 0.5,7,8,9)
#sum_charge_density_three_point('CHGCAR_50gpa', 0.6,7,8,9)
#sum_charge_density_three_point('CHGCAR_50gpa', 0.7,7,8,9)
#sum_charge_density_three_point('CHGCAR_50gpa', 0.8,7,8,9)
#sum_charge_density_three_point('CHGCAR_50gpa', 0.9,7,8,9)
#sum_charge_density_three_point('CHGCAR_50gpa', 1.1,7,8,9)
#sum_charge_density_three_point('CHGCAR_50gpa', 1.2,7,8,9)
#sum_charge_density_three_point('CHGCAR_50gpa', 1.3,7,8,9)
#sum_charge_density_three_point('CHGCAR_50gpa', 1.4,7,8,9)
#sum_charge_density_three_point('CHGCAR_50gpa', 1.5,7,8,9)