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


ob.obErrorLog.SetOutputLevel(0)

def convertible(v):
    try:
        int(v)
        return True
    except (TypeError, ValueError):
        return False

def is_float(string):
    try:
        float(string)
        return True
    except ValueError:
        return False


def remove_suffix(text, suffix):
    if text.endswith(suffix):
        return text[len(suffix):]
    return text

def create_ligand_pdbqt_file():
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

    return


def box_size_func(sdf_file, spacing, scale_factor):
    # Open SDF file
    sdf = glob.glob(sdf_file)
    sdf = ''.join(str(x) for x in sdf)
    sdf = open(sdf, 'r')
        
    # Extract coordinates 
    lines = [line.split() for line in sdf]
    del lines[:4]
        
    for j in range(0,len(lines)):
        del lines[j][4:]
        string = lines[j][3]
            
        if is_float(string) == True:
            del lines[j:]
            break

    coordinates = []
    elements = []
    x_axis = []
    y_axis = []
    z_axis = []
    
    for k in range(0,len(lines)):
        if lines[k][3] == ''+iv.var.metal_symbol+'':
            metal = lines[k][:3]
        
        coordinates.append(lines[k][:3])
        elements.append(lines[k][3])
        
        x_axis.append(float(lines[k][0]))
        y_axis.append(float(lines[k][1]))
        z_axis.append(float(lines[k][2]))
        
    # Shift axis to centre at metal
    metal = [float(i) for i in metal]
    metal = np.array(metal)
    
    x_max = np.max(x_axis-metal[0])
    x_min = np.min(x_axis-metal[0])
    
    y_max = np.max(y_axis-metal[1])
    y_min = np.min(y_axis-metal[1])

    z_max = np.max(z_axis-metal[2])
    z_min = np.min(z_axis-metal[2])


    x_axis = np.abs(np.max(x_axis-metal[0]) - np.min(x_axis-metal[0]))*scale_factor
    y_axis = np.abs(np.max(y_axis-metal[1]) - np.min(y_axis-metal[1]))*scale_factor
    z_axis = np.abs(np.max(z_axis-metal[2]) - np.min(z_axis-metal[2]))*scale_factor
    
    x_npts = (round(x_axis / spacing)) & (-2)
    y_npts = (round(y_axis / spacing)) & (-2)
    z_npts = (round(z_axis / spacing)) & (-2)

    max_side = max([x_npts,y_npts,z_npts])

    # Box dimensions in npts
    npts = [max_side, max_side, max_side]

    return npts

def docking_centre(coordinate_file):
    dock_site = open(coordinate_file,'r')
    coord = [line.split() for line in dock_site]

    if is_float(coord[0][0]) == True:
        dock = [coord[0][0], coord[0][1], coord[0][2]]
    else:
        dock = [coord[0][1], coord[0][2], coord[0][3]]

    return dock

def docking_func(parameter_set, name_ligand, name_protein, energy, dock, npts):
    os.system('cp '+os.environ['WORKING_DIR']+'/'+iv.var.parameter_file+' .')
    os.system('cp '+os.environ['WORKING_DIR']+'/'+iv.var.name_ligand+'_c.xyz .')
    create_ligand_pdbqt_file()

    # insert solution for R and epsilon for H-bond
    os.system(r'''awk '{ if ($2 == "'''+iv.var.metal_cap+'''" || $2 == "'''+iv.var.metal_symbol+'''") ($7 = '''+str(parameter_set[10])+''') && ($8 = '''+str(parameter_set[11])+'''); print $0}' '''+iv.var.parameter_file+''' > file_1''')
    os.system(r'''awk '{ if ($2 == "'''+iv.var.metal_cap+'''" || $2 == "'''+iv.var.metal_symbol+r'''") printf"%-8s %-3s %7s %8s %8s %9s %4s %4s %2s %3s %3s %2s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12; else print $0}' file_1 > '''+iv.var.parameter_file)
   
    #create_gpf():
    os.system(os.environ['PYTHON_2']+" "+os.environ['MGLTOOLS']+"/prepare_gpf4.py -l "+name_ligand+".pdbqt  -r clean_"+name_protein+".pdbqt -p parameter_file="+iv.var.parameter_file+" -p npts='{},{},{}'".format(npts[0],npts[1],npts[2])+" -p gridcenter='{:.4},{:.4},{:.4}' ".format(dock[0],dock[1],dock[2])+" >  /dev/null 2>&1")
    gpf = open('clean_'+name_protein+'.gpf', 'a')
    gpf.write('nbp_r_eps {:>.4f}  {:>.4f}'.format(parameter_set[0],parameter_set[1])+'    12 6 OA '+iv.var.metal_symbol+'\n')
    gpf.write('nbp_r_eps {:>.4f}  {:>.4f}'.format(parameter_set[2],parameter_set[3])+'    12 6 SA '+iv.var.metal_symbol+'\n')
    gpf.write('nbp_r_eps {:>.4f}  {:>.4f}'.format(parameter_set[4],parameter_set[5])+'    12 6 HD '+iv.var.metal_symbol+'\n')
    gpf.write('nbp_r_eps {:>.4f}  {:>.4f}'.format(parameter_set[6],parameter_set[7])+'    12 6 NA '+iv.var.metal_symbol+'\n')
    gpf.write('nbp_r_eps {:>.4f}  {:>.4f}'.format(parameter_set[8],parameter_set[9])+'    12 6  N '+iv.var.metal_symbol+'\n')
    gpf.close()

    #autogrid()
    os.system(os.environ['AUTODOCK']+'/autogrid4 -p clean_'+name_protein+'.gpf > /dev/null 2>&1')

    #create_dpf()
    d.write_dpf_file('clean_'+name_protein+'.gpf', name_ligand, 'clean_'+name_protein, iv.var.parameter_file, energy, random_pos=iv.var.random_position, SA=iv.var.docking_simulated_annealing, GA=iv.var.docking_genetic_algorithm)

    #autodock()
    os.system(os.environ['AUTODOCK']+'/autodock4 -p '+name_ligand+'_clean_'+name_protein+'.dpf >  /dev/null 2>&1')

    #write_all_conformations()
    os.system(os.environ['PYTHON_2']+" "+os.environ['MGLTOOLS']+"/write_conformations_from_dlg.py -d "+name_ligand+"_clean_"+name_protein+".dlg")

    return

def rmsd_func(name_ligand, n_prot, standard=False, test=False):
    rmsd_avg = []
    rmsd_list = []
    avg_list = []
    min_list = []
    rmsd_print_list = []
    
    if standard == True:
        output = [f"---------------------------------------     PROTEIN {n_prot} STANDARD     ---------------------------------------\n"]
    if test == True:
        output = [f"---------------------------------------     PROTEIN {n_prot} TEST         ---------------------------------------\n"]
            
    i = 1
    while os.path.exists(name_ligand+"_%i.pdbqt" % i):
        subprocess.call(os.environ['OBABEL']+" -ipdbqt "+name_ligand+"_{}.pdbqt".format(i)+" -oxyz "+name_ligand+"_{}.xyz".format(i)+" -d > "+name_ligand+"_{}.xyz".format(i), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) 

        rmsd_non_rotate = float(subprocess.getoutput([os.environ['PYTHON_3']+' '+os.environ['LIB_DIR']+'/calculate_rmsd.py ref.xyz '+name_ligand+'_{}.xyz'.format(i)+' --reorder --rotation none --translation none']))
        rmsd = rmsd_non_rotate

        rmsd_print_list.append("RMSD for Conformation %i = %.4f\n"% (i, rmsd))
        rmsd_list.append(rmsd)
        i += 1
    

    for j in range(0,len(rmsd_print_list)):
        output.append(rmsd_print_list[j])

        if standard == True:
            std = open(os.environ['WORKING_DIR']+'/standard/standard_conformations','a')
            std.write(rmsd_print_list[j])

            protein = open(os.environ['WORKING_DIR']+f'/standard/protein_{n_prot}','a')
            protein.write(rmsd_print_list[j])


        if test == True:   
            t = open(os.environ['WORKING_DIR']+'/test/test_conformations','a')
            t.write(rmsd_print_list[j])
            
            protein = open(os.environ['WORKING_DIR']+f'/test/protein_{n_prot}','a')
            protein.write(rmsd_print_list[j])

    avg_output = np.mean(rmsd_list)
    avg_list.append(avg_output)
    output.append(f"Average RMSD: {avg_output:.4}\n")

    minimum_rmsd = min(rmsd_list)
    min_list.append(minimum_rmsd)
    output.append(f"Lowest RMSD: {minimum_rmsd:.4}\n")

    stdv_rmsd = np.std(rmsd_list)
    output.append(f"Standard Deviation RMSD: {stdv_rmsd:.4}\n")

    var_rmsd = np.var(rmsd_list)
    output.append(f"Variance RMSD: {var_rmsd:.4}\n")
    output.append("-----------------------------------------------------------------------------------------------------------\n")

    return avg_list, min_list, print(''.join(output))


def plot_each_protein(protein_numbers):
    mpl.rcParams['figure.dpi'] = 400

    plt.rcParams.update({'font.size': 12})
    plt.rcParams['font.family'] = 'DeJavu Serif'
    plt.rcParams['font.serif'] = ['Palatino']

    standard_gen = []
    test_gen = []

    protein_number = protein_numbers

    standard_gen_dict = {}
    test_gen_dict = {}

    standard_protein_rmsd_dict = {}
    test_protein_rmsd_dict = {}

    for x in range(0,len(protein_number)):
        standard_gen_dict[x] = protein_number[x]
        standard_protein_rmsd_dict[protein_number[x]] = []
        
        standard_gen_dict[x] = os.environ['WORKING_DIR']+'/standard/protein_'+str(standard_gen_dict[x])
        
        protein = open(standard_gen_dict[x])
        protein_rmsd = [line.split() for line in protein]
        
        for i in range(0,len(protein_rmsd)):
            standard_gen.append(float(protein_rmsd[i][5]))
            standard_protein_rmsd_dict[protein_number[x]].append(float(protein_rmsd[i][5]))
            
        test_gen_dict[x] = protein_number[x]
        test_protein_rmsd_dict[protein_number[x]] = []
        
        test_gen_dict[x] = os.environ['WORKING_DIR']+'/test/protein_'+str(test_gen_dict[x])
        
        protein = open(test_gen_dict[x])
        protein_rmsd = [line.split() for line in protein]
        
        for i in range(0,len(protein_rmsd)):
            test_gen.append(float(protein_rmsd[i][5]))
            test_protein_rmsd_dict[protein_number[x]].append(float(protein_rmsd[i][5]))

        mu_first = statistics.mean(standard_protein_rmsd_dict[protein_number[x]])
        sigma_first = statistics.stdev(standard_protein_rmsd_dict[protein_number[x]])
        
        mu_last = statistics.mean(test_protein_rmsd_dict[protein_number[x]])
        sigma_last = statistics.stdev(test_protein_rmsd_dict[protein_number[x]])
        
        label_first = [f'Standard: $\mu$ = {mu_first:1.2f} $\sigma$ = {sigma_first:1.2f}']
        label_last = [f'GA Parameters: $\mu$ = {mu_last:1.2f} $\sigma$ = {sigma_last:1.2f}']

        bins = np.linspace(0, 15, 100)

        fig, ax = plt.subplots(1, figsize=(6,4), sharex=True)

        ax.hist(standard_protein_rmsd_dict[protein_number[x]], density=False, color='b', alpha=0.5, bins=bins, label=label_first)# density=False would make counts
        ax.hist(test_protein_rmsd_dict[protein_number[x]], density=False, color='r', alpha=0.5, bins =bins, label=label_last)

        ax.set_title('Protein {}'.format(protein_number[x]))

        ax.set_xticks(range(0,15,2))
        
        ax.xaxis.set_major_formatter(FormatStrFormatter('%i'))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%i'))

        ax.set_xlim(0,15)
        ax.set_ylabel('Conformations (N)')
        ax.set_xlabel('RMSD (Å)')

        plt.legend(loc='upper right')
    
        plt.tight_layout()
        plt.savefig(f"protein_"+str(protein_number[x])+".png", bbox_inches='tight')
        plt.close()

def plot_total_conformations():
    mpl.rcParams['figure.dpi'] = 400

    plt.rcParams.update({'font.size': 12})
    plt.rcParams['font.family'] = 'DeJavu Serif'
    plt.rcParams['font.serif'] = ['Palatino']

    standard_gen_all = []
    test_gen_all = []

    standard_gen = open(os.environ['WORKING_DIR']+'/standard/standard_conformations')
    standard_rmsd = [line.split() for line in standard_gen]
        
    for i in range(0,len(standard_rmsd)):
        standard_gen_all.append(float(standard_rmsd[i][5]))

    test_gen = open(os.environ['WORKING_DIR']+'/test/test_conformations')
    test_rmsd = [line.split() for line in test_gen]
        
    for i in range(0,len(test_rmsd)):
        test_gen_all.append(float(test_rmsd[i][5]))

    mu_first = statistics.mean(standard_gen_all)
    sigma_first = statistics.stdev(standard_gen_all)
    
    mu_last = statistics.mean(test_gen_all)
    sigma_last = statistics.stdev(test_gen_all)
    
    label_first = [f'Standard: $\mu$ = {mu_first:1.2f} $\sigma$ = {sigma_first:1.2f}']
    label_last = [f'GA Parameters: $\mu$ = {mu_last:1.2f} $\sigma$ = {sigma_last:1.2f}']

    bins = np.linspace(0, 15, 100)

    fig, ax = plt.subplots(1, figsize=(6,4), sharex=True)
    
    ax.hist(standard_gen_all, density=False, color='b', alpha=0.5, bins=bins, label=label_first)# density=False would make counts
    ax.hist(test_gen_all, density=False, color='r', alpha=0.5, bins =bins, label=label_last)
    
    ax.legend(loc="upper right")
    ax.set_title('All Conformations')
    
    ax.set_xlabel('RMSD (Å)')
    ax.set_ylabel('Conformations(N)')
    
    ax.set_xticks(range(0,15,2))

    ax.xaxis.set_major_formatter(FormatStrFormatter('%i'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%i'))

    plt.legend(loc='upper right')

    plt.tight_layout()
    plt.savefig(f"all_conformations.png", bbox_inches='tight')
    plt.close()


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