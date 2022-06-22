#!/usr/bin/python3
import os
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
        os.system(os.environ['OBABEL']+' -ixyz '+iv.var.name_ligand+'.xyz -oxyz '+iv.var.name_ligand+'_c.xyz --canonical > '+iv.var.name_ligand+'_c.xyz')

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
        os.system("cp  "+os.environ['WORKING_DIR']+"/clean_"+iv.var.name_protein+".pdb' .")
        #pdb.clean_protein_pdb(iv.var.pdb_file_protein)

    if os.path.exists('clean_'+iv.var.name_protein+'.pdb') == False:
        os.system("cp  "+os.environ['WORKING_DIR']+"/ref.pdb .")
        #pdb.get_ref_pdb(iv.var.pdb_file_protein)

    if iv.var.reference_docking == True:
        dock.get_coordinates()
    else:
        dock.users_coordinates()

    ##### GO ######
    os.chdir(os.environ['OUTPUT_DIR'])

    if os.path.isdir('gfnxtb') == False:
        os.mkdir('gfnxtb')
        os.chdir('gfnxtb')
    else:
        os.chdir('gfnxtb')

    os.system('mv '+os.environ['WORKING_DIR']+'/'+iv.var.name_ligand+'_c.xyz .')

    # If Geometry Converged Skip otherwise Run Again#
    if os.path.isdir(os.environ['OUTPUT_DIR']+'/gfnxtb/plams_workdir/plamsjob') == False:
        q.run_gfnxtb(iv.var.name_ligand+'_c.xyz')
        os.chdir(os.environ['OUTPUT_DIR']+'/gfnxtb/plams_workdir/plamsjob')
        q.gfnxtb_converged('ams.log')
        os.system(os.environ['AMSBIN']+'/amsreport dftb.rkf sdf > output.sdf')
    else:
        os.chdir(os.environ['OUTPUT_DIR']+'/gfnxtb/plams_workdir/plamsjob')
        q.gfnxtb_converged('ams.log')
        os.system(os.environ['AMSBIN']+'/amsreport dftb.rkf sdf > output.sdf')

    ###### Single Point ######
    os.chdir(os.environ['OUTPUT_DIR'])

    if os.path.isdir('single_point') == False:
        os.mkdir('single_point')
        os.chdir('single_point')
    else:
        os.chdir('single_point')

    os.system('cp '+os.environ['OUTPUT_DIR']+'/gfnxtb/plams_workdir/plamsjob/output.xyz .')
    os.system('mv output.xyz '+iv.var.name_ligand+'_c.xyz')

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

    os.system('cp '+os.environ['WORKING_DIR']+'/'+iv.var.parameter_file+' .')

    os.system(r'''awk '{ if ($2 == "RU" || $2 == "Ru") ($7 = '''+iv.var.r_Ru_Ru+''') && ($8 = '''+iv.var.e_Ru_Ru+'''); print $0}' '''+iv.var.parameter_file+''' > file_1''')
    os.system(r'''awk '{ if ($2 == "RU" || $2 == "Ru") printf "%-8s %-3s %7s %8s %8s %9s %4s %4s %2s %3s %3s %2s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12; else print $0}' file_1 > '''+iv.var.parameter_file)


    os.system('cp '+os.environ['OUTPUT_DIR']+'/gfnxtb/plams_workdir/plamsjob/output.xyz .')
    os.system('cp '+os.environ['OUTPUT_DIR']+'/gfnxtb/plams_workdir/plamsjob/output.sdf .')
    os.system('cp '+os.environ['OUTPUT_DIR']+'/single_point/plams_workdir/plamsjob/CM5_charges .')

    os.system('cp '+os.environ['OUTPUT_DIR']+'/file_prep/clean_'+iv.var.name_protein+'.pdb .')

    if iv.var.rmsd == False or iv.var.rmsd == None:
        os.system('cp '+os.environ['OUTPUT_DIR']+'/file_prep/ref.pdb .')

        #if os.path.exists(iv.var.name_ligand+'.pdbqt') == False:
        dock.create_ligand_pdbqt_file()

        #if os.path.exists('clean_'+iv.var.name_protein+'.pdb') == False:
        dock.prepare_receptor()

        dock.randomize_translation_rotation(iv.var.name_ligand+'.pdbqt')
        dock.add_to_dat_file()

        dock.create_gpf()
        dock.autogrid()

        dock.create_dpf()
        dock.autodock()

        dock.write_all_conformations()

    if iv.var.reference_docking == True:
        #if os.path.exists('ref.xyz') == False:
        #pdb = next(py.readfile('pdb','ref.pdb'))
        #pdb.write('xyz','ref.xyz',overwrite=True)

        rmsd_list = []

        os.system(os.environ['OBABEL']+" -isdf output.sdf -oxyz normalize.xyz -d > normalize.xyz")
        rmsd_normalize = float(subprocess.getoutput([os.environ['PYTHON_3']+' '+os.environ['DOCK_LIB_DIR']+'/calculate_rmsd.py ref.xyz normalize.xyz --reorder']))

        rmsd_list.append("RMSD between reference ligand and quantum optimized structure: %.4f" % rmsd_normalize)

        i = 1
        while os.path.exists(iv.var.name_ligand+"_%i.pdbqt" % i):
            os.system(os.environ['OBABEL']+" -ipdbqt "+iv.var.name_ligand+"_{}.pdbqt".format(i)+" -oxyz "+iv.var.name_ligand+"_{}.xyz".format(i)+" -d > "+iv.var.name_ligand+"_{}.xyz".format(i))

            rmsd_non_rotate = float(subprocess.getoutput([os.environ['PYTHON_3']+' '+os.environ['DOCK_LIB_DIR']+'/calculate_rmsd.py ref.xyz '+iv.var.name_ligand+'_{}.xyz'.format(i)+' -nh --reorder --rotation none --translation none']))
            rmsd = rmsd_non_rotate

            rmsd_list.append("RMSD for Conformation %i = %.4f"% (i, rmsd))
            i += 1

        for i in range(0,len(rmsd_list)):
            print(rmsd_list[i])
