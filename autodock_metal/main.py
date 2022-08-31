#!/usr/bin/python3
import os
from os import path
import argparse
import subprocess


from openbabel import pybel as py
from openbabel import openbabel as ob


import calculate_rmsd as rmsd
import docking as d
import pdb_extraction as pdb
import quantum_calculation as q
import input_variables as iv
import environment_variables

ob.obErrorLog.SetOutputLevel(0)

if __name__=='__main__':

    iv.insert_arguments()

    os.environ['WORKING_DIR']=os.getcwd()
    os.environ['OUTPUT_DIR']=os.getcwd()+'/output'

    ###### Canonicalize Ligand ######
    if os.path.exists(iv.var.name_ligand+'_c.xyz') == False:
        os.system(os.environ['OBABEL']+' -imol2 '+iv.var.name_ligand+'.mol2 -oxyz '+iv.var.name_ligand+'_c.xyz --canonical > '+iv.var.name_ligand+'_c.xyz')

    if iv.var.scale_factor != None:
        npts = d.box_size_func('*.sdf', 0.375, iv.var.scale_factor)
        
    if iv.var.box_size != None: 
        npts = iv.var.box_size

    ###### Generate Output Dir #######
    if os.path.isdir('output') == False:
        os.mkdir('output')
        os.chdir('output')
    else:
        os.chdir('output')

    if os.path.isdir('file_prep') == False:
        os.mkdir('file_prep')
        os.chdir('file_prep')
    else:
        os.chdir('file_prep')

    ###### Create pdb files ###### 
    os.system('cp '+os.environ['WORKING_DIR']+'/'+iv.var.pdb_file_protein+' .')

    if os.path.exists('clean_'+iv.var.name_protein+'.pdb') == False:
        os.system("cp  "+os.environ['WORKING_DIR']+"/"+iv.var.name_protein+".pdb .")
        #pdb.clean_protein_pdb(iv.var.pdb_file_protein)
        pdb.protonate_pdb(iv.var.name_protein)
        pdb.clean_protein_pdb('pdb_prot.pdb')

    ###### Single Point ######
    os.chdir(os.environ['OUTPUT_DIR'])

    if os.path.isdir('single_point') == False:
        os.mkdir('single_point')
        os.chdir('single_point')
    else:
        os.chdir('single_point')

    os.system('cp '+os.environ['WORKING_DIR']+'/'+iv.var.name_ligand+'_c.xyz .')

    # If single point successful Skip otherwise Run Again#
    if os.path.isdir(os.environ['OUTPUT_DIR']+'/single_point/plams_workdir/plamsjob') == False:
        q.plams_single_point(iv.var.name_ligand+'_c.xyz')
        os.chdir(os.environ['OUTPUT_DIR']+'/single_point/plams_workdir/plamsjob')
        q.single_point_check('ams.log')
        os.system(os.environ['AMSBIN']+'/amsreport adf.rkf CM5 > CM5_charges')
        os.system("grep 'kcal/mol' ams.log > energy")

    else:
        os.chdir(os.environ['OUTPUT_DIR']+'/single_point/plams_workdir/plamsjob')
        q.single_point_check('ams.log')
        os.system(os.environ['AMSBIN']+'/amsreport adf.rkf CM5 > CM5_charges')
        os.system("grep 'kcal/mol' ams.log > energy")

    ##### AutoDock #####
    os.chdir(os.environ['OUTPUT_DIR'])

    if os.path.isdir('docking') == False:
        os.mkdir('docking')
        os.chdir('docking')
    else:
        os.chdir('docking')

    os.system('cp '+os.environ['WORKING_DIR']+'/'+iv.var.parameter_file+' .')

    # os.system(r'''awk '{ if ($2 == "'''+iv.var.metal_cap+'''" || $2 == "'''+iv.var.metal_symbol+'''") ($7 = '''+iv.var.r_M+''') && ($8 = '''+iv.var.e_M+'''); print $0}' '''+iv.var.parameter_file+''' > file_1''')
    # os.system(r'''awk '{ if ($2 == "'''+iv.var.metal_cap+'''" || $2 == "'''+iv.var.metal_symbol+r'''") printf "%-8s %-3s %7s %8s %8s %9s %4s %4s %2s %3s %3s %2s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12; else print $0}' file_1 > '''+iv.var.parameter_file)
    # os.system("rm file_1")

    os.system('cp '+os.environ['OUTPUT_DIR']+'/file_prep/clean_'+iv.var.pdb_file_protein+' .')

    os.system('cp '+os.environ['WORKING_DIR']+'/'+iv.var.name_ligand+'.mol2 .')
    os.system('cp '+os.environ['OUTPUT_DIR']+'/single_point/plams_workdir/plamsjob/CM5_charges .')
    os.system('cp '+os.environ['OUTPUT_DIR']+'/single_point/plams_workdir/plamsjob/energy .')

    e = open('energy','r')
    lines = [line.split() for line in e]
    energy = lines[0][4]

    os.system('cp '+os.environ['WORKING_DIR']+'/'+iv.var.name_ligand+'_c.xyz .')

    if os.path.isfile('ref.xyz') == False:
        os.system(os.environ['OBABEL']+" -ixyz "+iv.var.name_ligand+"_c.xyz -oxyz ref.xyz -d > ref.xyz")

    if iv.var.dock_x and iv.var.dock_y and iv.var.dock_z != None:
        dock = d.users_coordinates()
    if iv.var.dock_x and iv.var.dock_y and iv.var.dock_z == None:
        dock = d.get_coordinates()

    if iv.var.standard == True:
        parameter_set = [2.0, 10.0, 2.0, 10.0, 2.0, 10.0, 2.0, 10.0, 2.0, 10.0, 2.0, 10.0]
        d.docking_func(parameter_set, iv.var.name_ligand, iv.var.name_protein, energy, dock, npts)
    else:
        parameter_set = [iv.var.r_OA, iv.var.e_OA, iv.var.r_SA, iv.var.e_SA, iv.var.r_HD, iv.var.e_HD, iv.var.r_NA, iv.var.e_NA, iv.var.r_N, iv.var.e_N, iv.var.r_M, iv.var.e_M]
        d.docking_func(parameter_set, iv.var.name_ligand, iv.var.name_protein, energy, dock, npts)

    rmsd_list = []
    i = 1
    while os.path.exists(iv.var.name_ligand+"_%i.pdbqt" % i):
        os.system(os.environ['OBABEL']+" -ipdbqt "+iv.var.name_ligand+"_{}.pdbqt".format(i)+" -oxyz "+iv.var.name_ligand+"_{}.xyz".format(i)+" -d > "+iv.var.name_ligand+"_{}.xyz".format(i))

        rmsd_non_rotate = float(subprocess.getoutput([os.environ['PYTHON_3']+' '+os.environ['DOCK_LIB_DIR']+'/calculate_rmsd.py ref.xyz '+iv.var.name_ligand+'_{}.xyz'.format(i)+' -nh --reorder --rotation none --translation none']))
        rmsd = rmsd_non_rotate

        rmsd_list.append("RMSD for Conformation %i = %.4f"% (i, rmsd))
        i += 1

    for i in range(0,len(rmsd_list)):
        print(rmsd_list[i])
