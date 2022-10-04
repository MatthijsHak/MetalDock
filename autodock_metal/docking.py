import os, glob
import numpy as np

from random import seed

import numpy as np


from openbabel import pybel as py
from rdkit import Chem
import rdkit
import networkx as nx

'''
All the docking parameters were obtained by a genetic alogrithm
'''

def is_float(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

def create_ligand_pdbqt_file(name_ligand):
    # Grep the correct part  of the itp file
    os.system("awk '/@<TRIPOS>ATOM/{flag=1; next} /@<TRIPOS>BOND/{flag=0} flag' "+name_ligand+".mol2  > almost")

    # Create charge file if CM5
    os.system("awk '{if (NR!=1) {print}}' CM5_charges > new")
    os.system(r'''awk '{printf "%8s\n",$2}' new > new_charge''')

    # Insert extra column
    os.system("paste -d' 'test almost new_charge > there")
    #os.system("paste -d' 'test almost charges > there")

    # Switch Columns
    os.system(r'''awk '{ printf "%7s %-3s %14s %9s %9s %-5s %3s %5s %12s \n",$1,$2,$3,$4,$5,$6,$7,$8,$10}' there > correct''')

    # Delete previous stuff
    os.system("sed -n '1,/@<TRIPOS>ATOM/p;/@<TRIPOS>BOND/,$p' "+name_ligand+".mol2 > ligand_almost")

    # Insert in ligand_par.itp
    os.system("sed '/@<TRIPOS>ATOM/ r correct' ligand_almost > "+name_ligand+".mol2")
    os.system("rm new new_charge ligand_almost correct there almost")

    #os.system(os.environ['PYTHON_2']+''' '''+os.environ['MGLTOOLS']+'''/prepare_ligand4.py -l '''+iv.var.name_ligand+'''.mol2 -U \""" -C''')
    pdbqt = next(py.readfile('mol2',name_ligand+'.mol2'))
    pdbqt.write('pdbqt',name_ligand+'.pdbqt',overwrite=True)

def get_coordinates():
    os.system('''awk '$1 == "'''+iv.var.metal_symbol+r'''" { print $0 }' ref.xyz > coordinates''')

    dock_site = open('coordinates','r')
    coord = [line.split() for line in dock_site]

    dock_x = str(coord[0][1])
    dock_y = str(coord[0][2])
    dock_z = str(coord[0][3])

    dock = [dock_x, dock_y, dock_z]

    return dock

def users_coordinates():
    dock_x = iv.var.dock_x
    dock_y = iv.var.dock_y
    dock_z = iv.var.dock_z

    dock = [dock_x, dock_y, dock_z]

    return dock

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

def prepare_receptor(pdb_file):
    os.system(os.environ['PYTHON_2']+' '+os.environ['MGLTOOLS']+'/prepare_receptor4.py -A check_hydrogens -r clean_'+pdb_file)


def docking_func(parameter_set, name_ligand, name_protein, energy, dock, npts):
    os.system('cp '+os.environ['WORKING_DIR']+'/'+iv.var.parameter_file+' .')
    os.system('cp '+os.environ['WORKING_DIR']+'/'+iv.var.name_ligand+'.mol2 .')

    create_ligand_pdbqt_file(iv.var.name_ligand)
    prepare_receptor(iv.var.pdb_file_protein)

    # insert parameters for R and epsilon for H-bond
    os.system(r'''awk '{ if ($2 == "'''+iv.var.metal_cap+'''" || $2 == "'''+iv.var.metal_symbol+'''") ($7 = '''+str(parameter_set[10])+''') && ($8 = '''+str(parameter_set[11])+'''); print $0}' '''+iv.var.parameter_file+''' > file_1''')
    os.system(r'''awk '{ if ($2 == "'''+iv.var.metal_cap+'''" || $2 == "'''+iv.var.metal_symbol+r'''") printf"%-8s %-3s %7s %8s %8s %9s %4s %4s %2s %3s %3s %2s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12; else print $0}' file_1 > '''+iv.var.parameter_file)

    #create_gpf():
    os.system(os.environ['PYTHON_2']+" "+os.environ['MGLTOOLS']+"/prepare_gpf4.py -l "+name_ligand+".pdbqt  -r clean_"+name_protein+".pdbqt -p parameter_file="+iv.var.parameter_file+" -p npts='{},{},{}'".format(npts[0],npts[1],npts[2])+" -p gridcenter='{:.4},{:.4},{:.4}' ".format(dock[0],dock[1],dock[2]))
    gpf = open('clean_'+name_protein+'.gpf', 'a')
    gpf.write('nbp_r_eps {:>.4f}  {:>.4f}'.format(parameter_set[0],parameter_set[1])+'    12 6 OA '+iv.var.metal_symbol+'\n')
    gpf.write('nbp_r_eps {:>.4f}  {:>.4f}'.format(parameter_set[2],parameter_set[3])+'    12 6 SA '+iv.var.metal_symbol+'\n')
    gpf.write('nbp_r_eps {:>.4f}  {:>.4f}'.format(parameter_set[4],parameter_set[5])+'    12 6 HD '+iv.var.metal_symbol+'\n')
    gpf.write('nbp_r_eps {:>.4f}  {:>.4f}'.format(parameter_set[6],parameter_set[7])+'    12 6 NA '+iv.var.metal_symbol+'\n')
    gpf.write('nbp_r_eps {:>.4f}  {:>.4f}'.format(parameter_set[8],parameter_set[9])+'    12 6  N '+iv.var.metal_symbol+'\n')
    gpf.close()

    #autogrid()
    os.system(os.environ['AUTODOCK']+'/autogrid4 -p clean_'+name_protein+'.gpf')

    #create_dpf()
    write_dpf_file('clean_'+name_protein+'.gpf', name_ligand, 'clean_'+name_protein, iv.var.parameter_file, energy, random_pos=iv.var.random_position, SA=iv.var.docking_simulated_annealing, GA=iv.var.docking_genetic_algorithm)

    #autodock()
    os.system(os.environ['AUTODOCK']+'/autodock4 -p '+name_ligand+'_clean_'+name_protein+'.dpf')

    #write_all_conformations()
    os.system(os.environ['PYTHON_2']+" "+os.environ['MGLTOOLS']+"/write_conformations_from_dlg.py -d "+name_ligand+"_clean_"+name_protein+".dlg")

    return

def write_dpf_file(gpf_file, name_ligand, name_protein, parameter_file, energy_ligand, random_pos=False, GA=False, SA=False):
    gpf_file = open(gpf_file,'r')
    gpf_lines = [line.split() for line in gpf_file]

    ligand_type = gpf_lines[5]
    del ligand_type[0]
    del ligand_type[-4:]
    ligand_type_str = ' '.join(ligand_type)

    dpf_file = open(name_ligand+'_'+name_protein+'.dpf','w')
    dpf_file.write('autodock_parameter_version 4.2       # used by autodock to validate parameter set\n')
    dpf_file.write('parameter_file '+parameter_file+'    # parameter library filename\n')
    dpf_file.write('outlev 1                             # diagnostic output level\n')
    dpf_file.write('intelec                              # calculate internal electrostatics\n')
    dpf_file.write('seed pid time                        # seeds for random generator\n')
    dpf_file.write('ligand_types '+ligand_type_str+'     # atoms types in ligand\n')
    dpf_file.write('fld '+name_protein+'.maps.fld              # grid_data_file\n')
    for i in range(0,len(ligand_type)):
        dpf_file.write('map '+name_protein+'.'+ligand_type[i]+'.map                 # atom-specific affinity map\n')

    dpf_file.write('elecmap '+name_protein+'.e.map             # electrostatics map\n')
    dpf_file.write('desolvmap '+name_protein+'.d.map           # desolvation map\n\n')

    dpf_file.write('# Unbound Ligand Parameters\n')
    dpf_file.write('unbound_energy '+str(energy_ligand)+'                      # set the energy of the unbound state\n')
    dpf_file.write('move '+name_ligand+'.pdbqt                # small molecule\n')

    if random_pos == True:
        dpf_file.write('tran0 random  # initial coordinates/A or random\n')
        dpf_file.write('quaternion0 random                   # initial orientation\n')
        dpf_file.write('dihe0 random                         # initial dihedrals (relative) or random\n')

    if GA == True:
        dpf_file.write('# GA parameters\n')
        dpf_file.write('ga_pop_size '+iv.var.ga_pop_size+'                      # number of individuals in population\n')
        dpf_file.write('ga_num_evals '+iv.var.ga_num_evals+'                 # maximum number of energy evaluations\n')
        dpf_file.write('ga_num_generations '+iv.var.ga_num_generation+'             # maximum number of generations\n')
        dpf_file.write('ga_elitism '+iv.var.ga_elitism+'                         # number of top individuals to survive to next generation\n')
        dpf_file.write('ga_mutation_rate '+iv.var.ga_mutation_rate+'                # rate of gene mutation\n')
        dpf_file.write('ga_crossover_rate '+iv.var.ga_crossover_rate+'                # rate of crossover\n')
        dpf_file.write('ga_window_size '+iv.var.ga_window_size+'                    # number of preceding generation when deciding threshold for worst individual current population\n')
        dpf_file.write('ga_cauchy_alpha 0.0                  # Alpha parameter of Cauchy distribution\n')
        dpf_file.write('ga_cauchy_beta 1.0                   # Beta parameter Cauchy distribution\n')

        dpf_file.write('# Local Search Parameters\n')
        dpf_file.write('sw_max_its 300                       # iterations of Solis & Wets local search\n')
        dpf_file.write('sw_max_succ 4                        # consecutive successes before changing rho\n')
        dpf_file.write('sw_max_fail 4                        # consecutive failures before changing rho\n')
        dpf_file.write('sw_rho 1.0                           # size of local search space to sample\n')
        dpf_file.write('sw_lb_rho 0.01                       # lower bound on rho\n')
        dpf_file.write('ls_search_freq 0.06                  # probability of performing local search on individual\n')
        dpf_file.write('# Activate LGA\n')
        dpf_file.write('set_ga                               # set the above parameters for GA or LGA\n')
        dpf_file.write('set_psw1                             # set the above pseudo-Solis & Wets parameters\n')
        dpf_file.write('ga_run 10                            # do this many hybrid GA-LS runs\n')
    if SA == True:
        dpf_file.write('# SA Parameters\n')
        dpf_file.write('tstep 2.0\n')
        #dpf_file.write('e0max 0.0 10000                      # max initial energy; max number of retries\n')
        dpf_file.write('linear_schedule                      # linear_schedule or geometric_schedule\n')
        dpf_file.write('rt0 500                              # initial annealing temperature (absolute tmperature multiplied by gas constant\n')
        dpf_file.write('rtrf 0.9                             # annealing temperature reductin factor < 1 cools > 1 heats system\n')
        dpf_file.write('runs 50                              # number of docking runs\n')
        dpf_file.write('cycles 50                            # number of temperature reduction cycles\n')
        dpf_file.write('accs 30000                           # maximum number of accepted steps per cycle\n')
        dpf_file.write('rejs 30000                           # maximum number of rejected steps per cycle\n')
        dpf_file.write('select m                             # m selects the minimum state, 1 selects the last state during each cycle\n')
        dpf_file.write('trnrf 1.0                            # per cycle reduction factor for translation steps\n')
        dpf_file.write('quarf 1.0                            # per cycle reduction factor for orientation steps\n')
        dpf_file.write('dihrf 1.0                            # per cycle reduction factor for torsional dihedral steps\n')

        dpf_file.write('# Activate SA\n')
        dpf_file.write('simanneal 10                         # run this many SA docking\n')

    dpf_file.write('# Perform Analysis\n')
    dpf_file.write('rmsref '+name_ligand+'.pdbqt              # RMSD will be calculated with this file\n')
    dpf_file.write('rmsmode atype                        # Method to calculate the RMSD\n')
    dpf_file.write('rmstol 2.0                           # RMSD tolerance\n')
    dpf_file.write('analysis                             # perforem a ranked cluster analysis\n')

    return


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
    os.system(r'''awk '{ printf "%-4s %5s %4s %4s %1s %3s %11s %7s %7s %5s %5s %9s %-3s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' new_output_test > '''+iv.var.name_ligand+'''.pdbqt''')

    os.system("rm temp_1 temp_2 temp_3 output_test new_output_test")


