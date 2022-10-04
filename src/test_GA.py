#!/usr/bin/env python3 

import argparse
import glob,sys,os, shutil
import subprocess
import pygad
import numpy as np
import statistics

import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter

from openbabel import pybel as py
from openbabel import openbabel as ob

import docking as d
import input_variables as iv
import environment_variables


''' Each position in list equals the following paramater:
[0]  = R_OA  [1]  = e_OA
[2]  = R_SA  [3]  = e_SA
[4]  = R_HD  [5]  = e_HD
[6]  = R_NA  [7]  = e_NA
[8]  = R_N   [9]  = e_N
[10] = R_M  [11]  = e_M

''' 

standard_set = {'V' : [2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0],
                'Cr': [2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0],
                'Fe': [2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0],
                'Co': [2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0],
                'Ni': [2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0],
                'Cu': [2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0],
                'Y' : [2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0],
                'Mo': [2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0],
                'Ru': [2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0],
                'Rh': [2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0],
                'Pd': [2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0],
                'Re': [2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0],
                'Os': [2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0],
                'Ir': [2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0],
                'Pt': [2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0],
                'Au': [2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 5.0]
            }




# def docking_func(parameter_set, name_ligand, name_protein, energy, dock, npts):
#     os.system('cp '+os.environ['WORKING_DIR']+'/'+iv.var.parameter_file+' .')
#     os.system('cp '+os.environ['WORKING_DIR']+'/'+iv.var.name_ligand+'_c.xyz .')
#     create_ligand_pdbqt_file()

#     # insert solution for R and epsilon for H-bond
#     os.system(r'''awk '{ if ($2 == "'''+iv.var.metal_cap+'''" || $2 == "'''+iv.var.metal_symbol+'''") ($7 = '''+str(parameter_set[10])+''') && ($8 = '''+str(parameter_set[11])+'''); print $0}' '''+iv.var.parameter_file+''' > file_1''')
#     os.system(r'''awk '{ if ($2 == "'''+iv.var.metal_cap+'''" || $2 == "'''+iv.var.metal_symbol+r'''") printf"%-8s %-3s %7s %8s %8s %9s %4s %4s %2s %3s %3s %2s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12; else print $0}' file_1 > '''+iv.var.parameter_file)
   
#     #create_gpf():
#     os.system(os.environ['PYTHON_2']+" "+os.environ['MGLTOOLS']+"/prepare_gpf4.py -l "+name_ligand+".pdbqt  -r clean_"+name_protein+".pdbqt -p parameter_file="+iv.var.parameter_file+" -p npts='{},{},{}'".format(npts[0],npts[1],npts[2])+" -p gridcenter='{:.4},{:.4},{:.4}' ".format(dock[0],dock[1],dock[2])+" >  /dev/null 2>&1")
#     gpf = open('clean_'+name_protein+'.gpf', 'a')
#     gpf.write('nbp_r_eps {:>.4f}  {:>.4f}'.format(parameter_set[0],parameter_set[1])+'    12 6 OA '+iv.var.metal_symbol+'\n')
#     gpf.write('nbp_r_eps {:>.4f}  {:>.4f}'.format(parameter_set[2],parameter_set[3])+'    12 6 SA '+iv.var.metal_symbol+'\n')
#     gpf.write('nbp_r_eps {:>.4f}  {:>.4f}'.format(parameter_set[4],parameter_set[5])+'    12 6 HD '+iv.var.metal_symbol+'\n')
#     gpf.write('nbp_r_eps {:>.4f}  {:>.4f}'.format(parameter_set[6],parameter_set[7])+'    12 6 NA '+iv.var.metal_symbol+'\n')
#     gpf.write('nbp_r_eps {:>.4f}  {:>.4f}'.format(parameter_set[8],parameter_set[9])+'    12 6  N '+iv.var.metal_symbol+'\n')
#     gpf.close()

#     #autogrid()
#     os.system(os.environ['AUTODOCK']+'/autogrid4 -p clean_'+name_protein+'.gpf > /dev/null 2>&1')

#     #create_dpf()
#     d.write_dpf_file('clean_'+name_protein+'.gpf', name_ligand, 'clean_'+name_protein, iv.var.parameter_file, energy, random_pos=iv.var.random_position, SA=iv.var.docking_simulated_annealing, GA=iv.var.docking_genetic_algorithm)

#     #autodock()
#     os.system(os.environ['AUTODOCK']+'/autodock4 -p '+name_ligand+'_clean_'+name_protein+'.dpf >  /dev/null 2>&1')

#     #write_all_conformations()
#     os.system(os.environ['PYTHON_2']+" "+os.environ['MGLTOOLS']+"/write_conformations_from_dlg.py -d "+name_ligand+"_clean_"+name_protein+".dlg")

#     return

# def rmsd_func(name_ligand, n_prot, train=False, standard=False, test=False):
#     rmsd_avg = []
#     rmsd_list = []
#     avg_list = []
#     min_list = []
#     rmsd_print_list = []
    
#     if standard == True:
#         output = [f"---------------------------------------     PROTEIN {n_prot} STANDARD     ---------------------------------------\n"]
#     if test == True:
#         output = [f"---------------------------------------     PROTEIN {n_prot} TEST         ---------------------------------------\n"]
#     if train == True:
#         output = [f"-------------------------------------------     PROTEIN {n_prot}      --------------------------------------------\n"]

#     i = 1
#     while os.path.exists(name_ligand+"_%i.pdbqt" % i):
#         subprocess.call(os.environ['OBABEL']+" -ipdbqt "+name_ligand+"_{}.pdbqt".format(i)+" -oxyz "+name_ligand+"_{}.xyz".format(i)+" -d > "+name_ligand+"_{}.xyz".format(i), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) 

#         rmsd_non_rotate = float(subprocess.getoutput([os.environ['PYTHON_3']+' '+os.environ['LIB_DIR']+'/calculate_rmsd.py ref.xyz '+name_ligand+'_{}.xyz'.format(i)+' --reorder --rotation none --translation none']))
#         rmsd = rmsd_non_rotate

#         rmsd_print_list.append("RMSD for Conformation %i = %.4f\n"% (i, rmsd))
#         rmsd_list.append(rmsd)
#         i += 1
    

#     for j in range(0,len(rmsd_print_list)):
#         output.append(rmsd_print_list[j])

#         if standard == True:
#             std = open(os.environ['WORKING_DIR']+'/standard/standard_conformations','a')
#             std.write(rmsd_print_list[j])

#             protein = open(os.environ['WORKING_DIR']+f'/standard/protein_{n_prot}','a')
#             protein.write(rmsd_print_list[j])


#         if test == True:   
#             t = open(os.environ['WORKING_DIR']+'/test/test_conformations','a')
#             t.write(rmsd_print_list[j])
            
#             protein = open(os.environ['WORKING_DIR']+f'/test/protein_{n_prot}','a')
#             protein.write(rmsd_print_list[j])

#     avg_output = np.mean(rmsd_list)
#     avg_list.append(avg_output)
#     output.append(f"Average RMSD: {avg_output:.4}\n")

#     minimum_rmsd = min(rmsd_list)
#     min_list.append(minimum_rmsd)
#     output.append(f"Lowest RMSD: {minimum_rmsd:.4}\n")

#     stdv_rmsd = np.std(rmsd_list)
#     output.append(f"Standard Deviation RMSD: {stdv_rmsd:.4}\n")

#     var_rmsd = np.var(rmsd_list)
#     output.append(f"Variance RMSD: {var_rmsd:.4}\n")
#     output.append("-----------------------------------------------------------------------------------------------------------\n")

#     return avg_list, min_list, print(''.join(output))



if __name__=='__main__':

    iv.insert_arguments()

    os.environ['WORKING_DIR']=os.getcwd()

    if os.path.isdir(os.environ['WORKING_DIR']+'/standard'):
        shutil.rmtree(os.environ['WORKING_DIR']+'/standard', ignore_errors=True)
        os.mkdir(os.environ['WORKING_DIR']+'/standard')
    else:
        os.mkdir(os.environ['WORKING_DIR']+'/standard')

    if os.path.isdir(os.environ['WORKING_DIR']+'/test'):
        shutil.rmtree(os.environ['WORKING_DIR']+'/test', ignore_errors=True)
        os.mkdir(os.environ['WORKING_DIR']+'/test')
    else:
        os.mkdir(os.environ['WORKING_DIR']+'/test')
    
    
    # Make list of the protein numbers to iterate over
    dir_list = os.listdir(os.getcwd())
    dir_list = [str(i).replace('protein_','') for i in dir_list]
    dir_list = [int(i) for i in dir_list if convertible(i)]
    dir_list = sorted(dir_list)

    # Parameters 
    desired_output = 0
    parameter_set = [iv.var.r_OA, iv.var.e_OA, iv.var.r_SA, iv.var.e_SA, iv.var.r_HD, iv.var.e_HD, iv.var.r_NA, iv.var.e_NA, iv.var.r_N, iv.var.e_N, iv.var.r_M, iv.var.e_M]

    std_avg_list = []
    std_min_avg_list = []

    test_avg_list = []
    test_min_avg_list = []

    for n_prot in dir_list:

        os.chdir(os.environ['WORKING_DIR']+f'/protein_{n_prot}')
        os.environ['OUTPUT_DIR'] = os.environ['WORKING_DIR']+f'/protein_{n_prot}/output'

        for files in glob.glob("*.xyz"):
            file_list = files.split('_c.')
            name_ligand = file_list[0]

        for files in glob.glob("*.pdb"):
            file_list = files.split('.')
            name_protein = file_list[0]

        if iv.var.scale_factor != None:
            npts = box_size_func('*.sdf', 0.375, iv.var.scale_factor)
        
        if iv.var.box_size != None: 
            npts = iv.var.box_size
            
        ##### AutoDock #####
        os.chdir(os.environ['OUTPUT_DIR'])

        if os.path.isdir('docking') == False:
            os.mkdir('docking')
            os.chdir('docking')
        else:
            os.chdir('docking')

        e = open('energy','r')
        lines = [line.split() for line in e]
        energy = lines[0][4]

        dock = docking_centre('coordinates')

        os.system('cp '+os.environ['WORKING_DIR']+'/'+iv.var.parameter_file+' .')
        os.system('cp '+os.environ['WORKING_DIR']+'/'+iv.var.name_ligand+'_c.xyz .')
        create_ligand_pdbqt_file()

        # dock with standard parameter set 
        docking_func(standard_set[iv.var.metal_symbol], name_ligand, name_protein, energy, dock, npts)
        avg_list_std, min_list_std, output_std = rmsd_func(name_ligand, n_prot, standard=True)

        std_avg_list.append(avg_list_std)
        std_min_avg_list.append(min_list_std)

        # dock with optimized parameter_set 
        docking_func(parameter_set, name_ligand, name_protein, energy, dock, npts)
        avg_list, min_list, output = rmsd_func(name_ligand, n_prot,test=True)

        test_avg_list.append(avg_list)
        test_min_avg_list.append(min_list)


    std_avg_output = np.mean(np.array((std_avg_list)))
    std_min_avg_rmsd = np.mean(np.array(std_min_avg_list))

    standard_fitness = 1 / np.abs(std_avg_output - desired_output)

    test_avg_output = np.mean(np.array((test_avg_list)))
    test_min_avg_rmsd = np.mean(np.array(test_min_avg_list))

    test_fitness = 1 / np.abs(test_avg_output - desired_output)

    output = [f"------------------------------------------     RESULTS         --------------------------------------------\n"]
    output.append(f"Results standard parameters: RMSD average  {std_avg_output:2.4};  Fitness {standard_fitness:2.4}\n")
    output.append(f"Results test parameters:     RMSD average {test_avg_output:2.4};  Fitness {test_fitness:2.4}\n")
    output.append("-----------------------------------------------------------------------------------------------------------\n")

    print(''.join(output))

    if os.path.isdir(os.environ['WORKING_DIR']+'/figures'):
        shutil.rmtree(os.environ['WORKING_DIR']+'/figures', ignore_errors=True)
        os.mkdir(os.environ['WORKING_DIR']+'/figures')
    else:
        os.mkdir(os.environ['WORKING_DIR']+'/figures')

    os.chdir(os.environ['WORKING_DIR']+'/figures')
    
    plot_each_protein(dir_list)
    plot_total_conformations()