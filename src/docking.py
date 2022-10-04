import os,sys
import subprocess
import environment_variables

import pdb_extraction as pdb
import prepare_dock as d

import adf_engine as adf 
import gaussian_engine as g
import psi4_engine as psi


from parser import Parser

def docking(input_file):
    par = Parser(input_file)

    input_dir = os.getcwd()
    output_dir = f'{input_dir}/output'

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

    if os.path.exists(par.name_ligand+'_c.xyz') == False:
        subprocess.call([os.environ['OBABEL']+f' -ixyz {input_dir}/'+par.xyz_file+' -oxyz '+par.name_ligand+'_c.xyz --canonical > '+par.name_ligand+'_c.xyz'],shell=True,  stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    if os.path.exists(par.name_ligand+'.mol2') == False:
        subprocess.call([os.environ['OBABEL']+' -ixyz '+par.name_ligand+'_c.xyz -omol2 '+par.name_ligand+'.mol2  > '+par.name_ligand+'.mol2'],shell=True,  stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    if os.path.exists(par.name_ligand+'.sdf') == False:
        subprocess.call([os.environ['OBABEL']+' -ixyz '+par.name_ligand+'_c.xyz -osdf '+par.name_ligand+'.sdf  > '+par.name_ligand+'.sdf'],shell=True,  stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    ###### Create pdb files ###### 
    os.system(f'cp {input_dir}/'+par.pdb_file+' .')

    if os.path.exists('clean_'+par.name_protein+'.pdb') == False:
        os.system(f"cp  {input_dir}/"+par.name_protein+".pdb .")
        pdb.protonate_pdb(par.pdb_file, par.pH)
        pdb.clean_protein_pdb(par.name_protein, par.pdb_file, par.clean_pdb)

    ###### Quantum Calculations ######
    os.chdir(f'{output_dir}')

    if os.path.isdir('QM') == False:
        os.mkdir('QM')
        os.chdir('QM')
    else:
        os.chdir('QM')

    if par.engine.lower() == 'adf':
        qm_dir = adf.adf_engine(par.name_ligand+'_c.xyz', par, output_dir)

    if par.engine.lower() == 'gaussian':
        qm_dir = g.gaussian_engine(par.name_ligand+'_c.xyz', par, output_dir)

    if par.engine.lower() == 'psi4':
        qm_dir = psi.psi4_engine(par.name_ligand+'_c.xyz', par, output_dir)
    
    # ## Geometry Optimization ##
    # if par.geom_opt == True:
    #     if os.path.isdir('geom_opt') == False:
    #         os.mkdir('geom_opt')
    #         os.chdir('geom_opt')
    #     else:
    #         os.chdir('geom_opt')

    #     subprocess.call([f'cp {output_dir}/file_prep/'+par.name_ligand+'_c.xyz .'], shell=True)

    #     # If Geometry Converged Skip otherwise Run Again#
    #     if os.path.isdir(f'{output_dir}/QM/geom_opt/plams_workdir/plamsjob') == False:
    #         q.geom_opt(par.name_ligand+'_c.xyz', par)
    #         os.chdir(f'{output_dir}/QM/geom_opt/plams_workdir/plamsjob')
    #         q.gfnxtb_converged('ams.log')
    #     else:
    #         os.chdir(f'{output_dir}/QM/geom_opt/plams_workdir/plamsjob')
    #         q.gfnxtb_converged('ams.log')


    # ## Single Point ##
    # os.chdir(f'{output_dir}/QM')

    # if os.path.isdir('single_point') == False:
    #     os.mkdir('single_point')
    #     os.chdir('single_point')
    # else:
    #     os.chdir('single_point')

    # if par.geom_opt == True:
    #     subprocess.call([f'cp {output_dir}/QM/geom_opt/plams_workdir/plamsjob/output.xyz '+par.name_ligand+'_c.xyz'], shell=True)
    # else:
    #     subprocess.call([f'cp {output_dir}/file_prep/'+par.name_ligand+'_c.xyz .'], shell=True)

    # # If single point successful Skip otherwise Run Again#
    # if os.path.isdir(f'{output_dir}/QM/single_point/plams_workdir/plamsjob') == False:
    #     q.single_point(par.name_ligand+'_c.xyz', par)
    #     os.chdir(f'{output_dir}/QM/single_point/plams_workdir/plamsjob')
    #     q.single_point_check('ams.log')
    #     subprocess.call([os.environ['AMSBIN']+'/amsreport adf.rkf CM5 > CM5_charges'], shell=True)
    #     subprocess.call(["grep 'kcal/mol' ams.log > energy"], shell=True)
    # else:
    #     os.chdir(f'{output_dir}/QM/single_point/plams_workdir/plamsjob')
    #     q.single_point_check('ams.log')
    #     subprocess.call([os.environ['AMSBIN']+'/amsreport adf.rkf CM5 > CM5_charges'], shell=True)
    #     subprocess.call(["grep 'kcal/mol' ams.log > energy"], shell=True)

      
    ##### AutoDock #####
    os.chdir(f'{output_dir}')

    if os.path.isdir('docking') == False:
        os.mkdir('docking')
        os.chdir('docking')
    else:
        os.chdir('docking')

    subprocess.call([f'cp {input_dir}/'+par.parameter_file+' .'], shell=True)

    subprocess.call([f'cp {output_dir}/file_prep/clean_'+par.name_protein+'.pdb .'], shell=True)

    subprocess.call([f'cp {output_dir}/file_prep/'+par.name_ligand+'.mol2 .'], shell=True)
    subprocess.call([f'cp {qm_dir}/CM5_charges .'], shell=True)
    subprocess.call([f'cp {qm_dir}/energy .'], shell=True)

    e = open('energy','r')
    lines = [line.split() for line in e]
    energy = lines[0][4]

    subprocess.call([f'cp {output_dir}/file_prep/'+par.name_ligand+'_c.xyz .'], shell=True)

    if par.rmsd == True:
        if os.path.isfile('ref.xyz') == False:
            subprocess.call([os.environ['OBABEL']+" -ixyz "+par.name_ligand+"_c.xyz -oxyz ref.xyz -d > ref.xyz"], shell=True)

    if par.dock_x and par.dock_y and par.dock_z != None:
        dock = d.users_coordinates(par.dock_x, par.dock_y, par.dock_z)
    else:
        dock = d.get_coordinates(par.metal_symbol)

    if par.standard == True:
        parameter_set = [2.0, 10.0, 2.0, 10.0, 2.0, 10.0, 2.0, 10.0, 2.0, 10.0, 2.0, 10.0]
        d.create_ligand_pdbqt_file(par.name_ligand)
        d.prepare_receptor(par.name_protein)
        d.docking_func(parameter_set, par.parameter_file, par.metal_symbol, par.name_ligand, par.name_protein, energy, dock, par.box_size, par.dock_algorithm, par.random_pos, par.ga_dock, par.sa_dock)
    else:
        d.create_ligand_pdbqt_file(par.name_ligand)
        d.prepare_receptor(par.name_protein)
        d.docking_func(par.parameter_set, par.parameter_file, par.metal_symbol, par.name_ligand, par.name_protein, energy, dock, par.box_size, par.dock_algorithm, par.random_pos, par.ga_dock, par.sa_dock)

    if par.rmsd == True:
        rmsd_list = []
        i = 1
        while os.path.exists(par.name_ligand+"_%i.pdbqt" % i):
            subprocess.call([os.environ['OBABEL']+" -ipdbqt "+par.name_ligand+"_{}.pdbqt".format(i)+" -oxyz "+par.name_ligand+"_{}.xyz".format(i)+" -d > "+par.name_ligand+"_{}.xyz".format(i)], shell=True)

            rmsd_non_rotate = float(subprocess.getoutput([os.environ['PYTHON_3']+' '+os.environ['DOCK_LIB_DIR']+'/calculate_rmsd.py ref.xyz '+par.name_ligand+'_{}.xyz'.format(i)+' -nh --reorder --rotation none --translation none']))
            rmsd = rmsd_non_rotate

            rmsd_list.append("RMSD for Conformation %i = %.4f"% (i, rmsd))
            i += 1

        for i in range(0,len(rmsd_list)):
            print(rmsd_list[i])


    ##### results #####
    os.chdir(f'{output_dir}')

    if os.path.isdir('results') == False:
        os.mkdir('results')
        os.chdir('results')
    else:
        os.chdir('results')

    subprocess.call([f'cp {output_dir}/docking/'+par.name_ligand+'_*.pdbqt .'], shell=True)
    subprocess.call([f'cp {output_dir}/docking/clean_'+par.pdb_file+' .'], shell=True)
 
    print("\nDOCKING SUCCESFULLY COMPLETED\n")
    print("THE PRINTED POSES AND PROTEIN CAN BE FOUND IN THE RESULTS DIRECTORY\n")
    return