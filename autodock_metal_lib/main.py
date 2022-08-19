#!/usr/bin/python3
import os
from os import path
import argparse
import subprocess


from openbabel import pybel as py
from openbabel import openbabel as ob


import calculate_rmsd as rmsd
import docking as dock
import pdb_extraction as pdb
import quantum_calculation as q
import input_variables as iv
import environment_variables

from resp_charges import resp_charges

ob.obErrorLog.SetOutputLevel(0)

if __name__=='__main__':

    iv.insert_arguments()

    os.environ['WORKING_DIR']=os.getcwd()
    os.environ['OUTPUT_DIR']=os.getcwd()+'/output'

    ###### Canonicalize Ligand ######
    if os.path.exists(iv.var.name_ligand+'_c.xyz') == False:
        os.system(os.environ['OBABEL']+' -imol2 '+iv.var.name_ligand+'.mol2 -oxyz '+iv.var.name_ligand+'_c.xyz --canonical > '+iv.var.name_ligand+'_c.xyz')

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
        #resp_charges('../../'+iv.var.name_ligand+'_c.xyz')
    else:
        os.chdir(os.environ['OUTPUT_DIR']+'/single_point/plams_workdir/plamsjob')
        q.single_point_check('ams.log')
        os.system(os.environ['AMSBIN']+'/amsreport adf.rkf CM5 > CM5_charges')
       #resp_charges('../../'+iv.var.name_ligand+'_c.xyz')

    ##### AutoDock #####
    os.chdir(os.environ['OUTPUT_DIR'])

    if os.path.isdir('docking') == False:
        os.mkdir('docking')
        os.chdir('docking')
    else:
        os.chdir('docking')

    #os.system('cp '+os.environ['WORKING_DIR']+'/'+iv.var.parameter_file+' .')

    #os.system(r'''awk '{ if ($2 == "'''+iv.var.metal_cap+'''" || $2 == "'''+iv.var.metal_symbol+'''") ($7 = '''+iv.var.r_Ru_Ru+''') && ($8 = '''+iv.var.e_Ru_Ru+'''); print $0}' '''+iv.var.parameter_file+''' > file_1''')
    #os.system(r'''awk '{ if ($2 == "'''+iv.var.metal_cap+'''" || $2 == "'''+iv.var.metal_symbol+r'''") printf "%-8s %-3s %7s %8s %8s %9s %4s %4s %2s %3s %3s %2s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12; else print $0}' file_1 > '''+iv.var.parameter_file)
    #os.system("rm file_1")

    os.system('cp '+os.environ['OUTPUT_DIR']+'/file_prep/clean_'+iv.var.pdb_file_protein+' .')

    os.system('cp '+os.environ['WORKING_DIR']+'/'+iv.var.name_ligand+'.mol2 .')
    os.system('cp '+os.environ['OUTPUT_DIR']+'/single_point/plams_workdir/plamsjob/CM5_charges .')

    os.system('cp '+os.environ['WORKING_DIR']+'/'+iv.var.name_ligand+'_c.xyz .')

    if os.path.isfile('ref.xyz') == False:
        os.system(os.environ['OBABEL']+" -ixyz "+iv.var.name_ligand+"_c.xyz -oxyz ref.xyz -d > ref.xyz")

    if iv.var.reference_docking == True:
        #dock.users_coordinates()
        dock.get_coordinates()
    else:
        dock.users_coordinates()

    if iv.var.rmsd == False or iv.var.rmsd == None:
        dock.docking()

    if iv.var.reference_docking == True:

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
