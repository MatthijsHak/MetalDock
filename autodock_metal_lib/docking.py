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

def create_ligand_pdbqt_file():
    #mol2 = next(py.readfile('xyz',''+iv.var.name_ligand+'_c.xyz'))
    #mol2.write('mol2',iv.var.name_ligand+'.mol2',overwrite=True)
    #mol2 = next(py.readfile('xyz','output.xyz'))
    #mol2.write('mol2',iv.var.name_ligand+'.mol2',overwrite=True)

    # Grep the correct part  of the itp file
    os.system("awk '/@<TRIPOS>ATOM/{flag=1; next} /@<TRIPOS>BOND/{flag=0} flag' "+iv.var.name_ligand+".mol2  > almost")

    # Create charge file if CM5
    os.system("awk '{if (NR!=1) {print}}' CM5_charges > new")
    os.system(r'''awk '{printf "%8s\n",$2}' new > new_charge''')

    # Insert extra column
    os.system("paste -d' 'test almost new_charge > there")
    #os.system("paste -d' 'test almost charges > there")

    # Switch Columns
    os.system(r'''awk '{ printf "%7s %-3s %14s %9s %9s %-5s %3s %5s %12s \n",$1,$2,$3,$4,$5,$6,$7,$8,$10}' there > correct''')

    # Delete previous stuff
    os.system("sed -n '1,/@<TRIPOS>ATOM/p;/@<TRIPOS>BOND/,$p' "+iv.var.name_ligand+".mol2 > ligand_almost")

    # Insert in ligand_par.itp
    os.system("sed '/@<TRIPOS>ATOM/ r correct' ligand_almost > "+iv.var.name_ligand+".mol2")
    os.system("rm new new_charge ligand_almost correct there almost")

    #os.system(os.environ['PYTHON_2']+''' '''+os.environ['MGLTOOLS']+'''/prepare_ligand4.py -l '''+iv.var.name_ligand+'''.mol2 -U \""" -C''')
    pdbqt = next(py.readfile('mol2',iv.var.name_ligand+'.mol2'))
    pdbqt.write('pdbqt',iv.var.name_ligand+'.pdbqt',overwrite=True)

def get_coordinates():
    os.system('''awk '$1 == "Ru" { print $0 }' ref.xyz > coordinates''')

    dock_site = open('coordinates','r')
    coord = [line.split() for line in dock_site]

    global dock_x, dock_y, dock_z

    dock_x = str(coord[0][1])
    dock_y = str(coord[0][2])
    dock_z = str(coord[0][3])

def users_coordinates():
    global dock_x, dock_y, dock_z

    dock_x = iv.var.dock_x
    dock_y = iv.var.dock_y
    dock_z = iv.var.dock_z


def prepare_receptor():
    os.system(os.environ['PYTHON_2']+' '+os.environ['MGLTOOLS']+'/prepare_receptor4.py -A check_hydrogens -r clean_'+iv.var.name_protein+'.pdb')

def add_to_dat_file():
    dat = open(''+iv.var.parameter_file+'', 'a')
    dat.write('nbp_r_eps '+iv.var.r_OA+'   '+iv.var.e_OA+' 12 6 OA '+iv.var.metal_symbol+'\n')
    dat.write('nbp_r_eps '+iv.var.r_SA+'   '+iv.var.e_SA+' 12 6 SA '+iv.var.metal_symbol+'\n')
    dat.write('nbp_r_eps '+iv.var.r_HD+'   '+iv.var.e_HD+' 12 6 HD '+iv.var.metal_symbol+'\n')
    dat.write('nbp_r_eps '+iv.var.r_NA+'   '+iv.var.e_NA+' 12 6 NA '+iv.var.metal_symbol+'\n')
    dat.write('nbp_r_eps '+iv.var.r_N+'   '+iv.var.e_N+' 12 6  N '+iv.var.metal_symbol+'\n')
    dat.close()

def create_gpf():
    os.system(os.environ['PYTHON_2']+" "+os.environ['MGLTOOLS']+"/prepare_gpf4.py -l "+iv.var.name_ligand+".pdbqt  -r clean_"+iv.var.name_protein+".pdbqt -p parameter_file="+iv.var.parameter_file+" -p npts='"+iv.var.box_size+"' -p gridcenter='"+dock_x+","+dock_y+","+dock_z+"'")
    gpf = open('clean_'+iv.var.name_protein+'.gpf', 'a')
    gpf.write('nbp_r_eps '+iv.var.r_OA+'   '+iv.var.e_OA+' 12 6 OA '+iv.var.metal_symbol+'\n')
    gpf.write('nbp_r_eps '+iv.var.r_SA+'   '+iv.var.e_SA+' 12 6 SA '+iv.var.metal_symbol+'\n')
    gpf.write('nbp_r_eps '+iv.var.r_HD+'   '+iv.var.e_HD+' 12 6 HD '+iv.var.metal_symbol+'\n')
    gpf.write('nbp_r_eps '+iv.var.r_NA+'   '+iv.var.e_NA+' 12 6 NA '+iv.var.metal_symbol+'\n')
    gpf.write('nbp_r_eps '+iv.var.r_N+'   '+iv.var.e_N+' 12 6  N '+iv.var.metal_symbol+'\n')
    gpf.close()

def autogrid():
    os.system(os.environ['AUTODOCK']+'/autogrid4 -p clean_'+iv.var.name_protein+'.gpf')

def create_dpf():
    os.system(os.environ['PYTHON_2']+" "+os.environ['MGLTOOLS']+"/prepare_dpf42.py -l "+iv.var.name_ligand+".pdbqt -r clean_"+iv.var.name_protein+".pdb -p parameter_file="+iv.var.parameter_file)

def autodock():
    os.system(os.environ['AUTODOCK']+'/autodock4 -p '+iv.var.name_ligand+'_clean_'+iv.var.name_protein+'.dpf')

def write_all_conformations():
     os.system(os.environ['PYTHON_2']+" "+os.environ['MGLTOOLS']+"/write_conformations_from_dlg.py -d "+iv.var.name_ligand+"_clean_"+iv.var.name_protein+".dlg")

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

    output_array = [["%.3f" % j for j in i] for i in random_coord]
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
    os.system(r'''awk '{ printf "%-4s %6s %3s %4s %1s %3s %11s %7s %7s %5s %5s %9s %-3s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' new_output_test > '''+iv.var.name_ligand+'''.pdbqt''')

    os.system("rm temp_1 temp_2 temp_3 output_test new_output_test")


