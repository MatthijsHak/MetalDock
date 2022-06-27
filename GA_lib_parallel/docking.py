import os
import rdkit

import random
from random import seed

import numpy as np
import networkx as nx

from openbabel import pybel as py
from rdkit import Chem

import input_variables as iv
import variable_class as vc

def distance(x1,x2,y1,y2,z1,z2):

    d = np.sqrt( (x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2  )
    d = format(d, '.4f')

    return str(d)

# Translation Matrix - Translation is switched from column to row because we work with (1,4) coordinate matrix
def translation_matrix(matrix,dx,dy,dz):
    translation_matrix = np.array(([1,0,0,0],
                                   [0,1,0,0],
                                   [0,0,1,0],
                                   [dx,dy,dz,1]),dtype=np.float64)

    y_dim,x_dim = np.shape(matrix)
    extra_dim_arr = np.c_[matrix, np.ones(y_dim)]

    new_extra_dim_arr = np.dot(extra_dim_arr, translation_matrix)
    new_arr = np.delete(new_extra_dim_arr, (3), 1)

    return new_arr

# Rotation Matrices
def x_axis_rotation(matrix,theta):
    rotation_matrix = np.array(([1,0,0,0],
                                [0,np.cos(theta),-np.sin(theta),0],
                                [0,np.sin(theta),np.cos(theta),0],
                                [0,0,0,1]), dtype=np.float64)

    y_dim,x_dim = np.shape(matrix)
    extra_dim_arr = np.c_[matrix, np.ones(y_dim)]

    new_extra_dim_arr = np.dot(extra_dim_arr, rotation_matrix)
    new_arr = np.delete(new_extra_dim_arr, (3), 1)

    return new_arr

def y_axis_rotation(matrix,theta):
    rotation_matrix = np.array(([np.cos(theta),0,np.sin(theta),0],
                                [0,1,0,0],
                                [-np.sin(theta),0,np.cos(theta),0],
                                [0,0,0,1]),dtype=np.float64)

    y_dim,x_dim = np.shape(matrix)
    extra_dim_arr = np.c_[matrix, np.ones(y_dim)]

    new_extra_dim_arr = np.dot(extra_dim_arr, rotation_matrix)
    new_arr = np.delete(new_extra_dim_arr, (3), 1)

    return new_arr

def z_axis_rotation(matrix,theta):
    rotation_matrix = np.array(([np.cos(theta),-np.sin(theta),0,0],
                                [np.sin(theta),np.cos(theta),0,0],
                                [0,0,1,0],
                                [0,0,0,1]),dtype=np.float64)

    y_dim,x_dim = np.shape(matrix)
    extra_dim_arr = np.c_[matrix, np.ones(y_dim)]

    new_extra_dim_arr = np.dot(extra_dim_arr, rotation_matrix)
    new_arr = np.delete(new_extra_dim_arr, (3), 1)

    return new_arr

def randomize_translation_rotation(pdbqt_file):
    os.system('''awk '{if ($1 == "ATOM") print $0}' '''+pdbqt_file+''' > temp_1''')

    dock_site = open('temp_1','r')
    coord = [line.split() for line in dock_site]

    for i in range(0,len(coord)):
        del coord[i][0:6]
        del coord[i][3:7]

    coord_float = [[float(j) for j in i] for i in coord]
    coord_array = np.array(coord_float)

    # Random Coordinates
    random.seed()

    theta = random.uniform(0, 2*np.pi) # rotation angle

    box_size = list(iv.var.box_size.split(','))
    box_size = [float(x) for x in box_size]
    box_size = [(x*0.375) / 2 for x in box_size]

    dx = random.uniform(-box_size[0], box_size[0]) # translation x-direction
    dy = random.uniform(-box_size[1], box_size[1]) # translation y-direction
    dz = random.uniform(-box_size[2], box_size[2]) # translation z-direction

    random_coord = translation_matrix(coord_array,dx,dy,dz)
    random_coord = x_axis_rotation(random_coord,theta)
    random_coord = y_axis_rotation(random_coord,theta)
    random_coord = z_axis_rotation(random_coord,theta)

    output_array = [["%5.2f" % j for j in i] for i in random_coord]
    output_array = [[str(j) for j in i] for i in output_array]

    output = open('output_test','w')

    for elem in output_array:
        output.write('\t'.join(elem))
        output.write('\n')

    output.close()

    os.system("paste -d' 'temp_1 temp_1 output_test> temp_2")
    os.system(r'''awk '{ printf "%-4s %6s %2s %5s %1s %3s %11s %7s %7s %5s %5s %9s %-3s\n",$1,$2,$3,$4,$5,$6,$14,$15,$16,$10,$11,$12,$13}' temp_2 > temp_3''')

    new_file = open('temp_3','r')
    new_file_coord = [line.split() for line in new_file]

    old_file = open(''+pdbqt_file+'','r')
    old_file_coord = [line.split() for line in old_file]

    new_pdbqt = []
    k = 0

    for i in range(0,len(old_file_coord)):
        for j in range(0,len(new_file_coord)):
            if old_file_coord[i][0] == new_file_coord[j][0]:
                new_pdbqt.append(new_file_coord[k])
                k+=1
                break
            else:
                new_pdbqt.append(old_file_coord[i])
                break

    new_output = open('new_output_test','w')

    for elem in new_pdbqt:
        new_output.write('\t'.join(elem))
        new_output.write('\n')

    new_output.close()
    os.system(r'''awk '{ printf "%-4s %5s %4s %4s %1s %3s %11s %7s %7s %5s %5s %9s %-3s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' new_output_test > docking.pdbqt''')

    os.system("rm temp_1 temp_2 temp_3 output_test new_output_test")

