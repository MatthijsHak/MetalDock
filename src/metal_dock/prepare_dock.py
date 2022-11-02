import os, sys, glob, subprocess
import random

import numpy as np
from random import seed


def is_float(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

def check_pdbqt(pdbqt_file):
    count = 0
    with open(pdbqt_file) as fin:
        for line in fin:
            if 'MODEL' in line:
                count+=1
    if count > 1:
        print('LIGAND CANNOT BE RECOGNIZED AS ONE MOLECULE')
        print('REOPTIMIZE COMPOUND WITH DIFFERENT INITIAL SETTINGS OR DELETE NON-COVALENT LIGANDS\n')
        sys.exit()

def create_ligand_pdbqt_file(name_ligand):
    # Grep the correct part  of the itp file
    subprocess.call(["awk '/@<TRIPOS>ATOM/{flag=1; next} /@<TRIPOS>BOND/{flag=0} flag' "+name_ligand+".mol2  > almost"], shell=True)

    # Create charge file if CM5
    subprocess.call(["awk '{if (NR!=1) {print}}' CM5_charges > new"], shell=True)
    subprocess.call([r'''awk '{printf "%8s\n",$2}' new > new_charge'''], shell=True)

    # Insert extra column
    subprocess.call(["paste -d' 'test almost new_charge > there"], shell=True)

    # Switch Columns
    subprocess.call([r'''awk '{ printf "%7s %-3s %14s %9s %9s %-5s %3s %5s %12s \n",$1,$2,$3,$4,$5,$6,$7,$8,$10}' there > correct'''], shell=True)

    # Delete previous stuff
    subprocess.call(["sed -n '1,/@<TRIPOS>ATOM/p;/@<TRIPOS>BOND/,$p' "+name_ligand+".mol2 > ligand_almost"], shell=True)

    # Insert in ligand_par.itp
    subprocess.call(["sed '/@<TRIPOS>ATOM/ r correct' ligand_almost > "+name_ligand+".mol2"], shell=True)
    subprocess.call(["rm new new_charge ligand_almost correct there almost"], shell=True)

    subprocess.call([os.environ['OBABEL']+' -imol2 '+name_ligand+'.mol2 -opdbqt '+name_ligand+'.pdbqt  > '+name_ligand+'.pdbqt'],shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    check_pdbqt(name_ligand+'.pdbqt')
    return

def get_coordinates(xyz_file, metal_symbol):
    subprocess.call(['''awk '$1 == "'''+metal_symbol+r'''" { print $0 }' '''+xyz_file+''' > coordinates'''], shell=True)

    dock_site = open('coordinates','r')
    coord = [line.split() for line in dock_site]

    dock_x = str(coord[0][1])
    dock_y = str(coord[0][2])
    dock_z = str(coord[0][3])

    dock = [dock_x, dock_y, dock_z]

    return dock

def users_coordinates(dock_x, dock_y, dock_z):
    dock = [dock_x, dock_y, dock_z]
    return dock

def box_size_func(xyz_file, metal_symbol, spacing, scale_factor):
    # Open xyz file
    xyz = open(xyz_file, 'r')

    # Extract coordinates 
    lines = [line.split() for line in xyz]
    del lines[:2]

    coordinates = []
    elements = []
    x_axis = []
    y_axis = []
    z_axis = []

    for k in range(0,len(lines)):
        if lines[k][0] == ''+metal_symbol+'':
            metal = lines[k][1:4]

        coordinates.append(lines[k][:3])
        elements.append(lines[k][0])

        x_axis.append(float(lines[k][1]))
        y_axis.append(float(lines[k][2]))
        z_axis.append(float(lines[k][3]))
        
    # Shift axis to centre at metal
    metal = [float(i) for i in metal]
    metal = np.array(metal)

    x_dist = np.abs(np.max(x_axis-metal[0]) - np.min(x_axis-metal[0]))*scale_factor
    y_dist = np.abs(np.max(y_axis-metal[1]) - np.min(y_axis-metal[1]))*scale_factor
    z_dist = np.abs(np.max(z_axis-metal[2]) - np.min(z_axis-metal[2]))*scale_factor

    if x_dist > 20:
        x_dist = 20
    if y_dist > 20:
        y_dist = 20
    if z_dist > 20:
        z_dist = 20

    x_npts = (round(x_dist / spacing)) & (-2)
    y_npts = (round(y_dist / spacing)) & (-2)
    z_npts = (round(z_dist / spacing)) & (-2)

    max_side = max([x_npts,y_npts,z_npts])
    #print('Box size is {} {} {}'.format(max_side,max_side,max_side))

    return max_side

def prepare_receptor(name_protein):
    subprocess.call([os.environ['PYTHON_2']+' '+os.environ['MGLTOOLS']+'/prepare_receptor4.py -A check_hydrogens -r clean_'+name_protein+'.pdb'],shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    return

def docking_func(parameter_set, parameter_file, metal_symbol, name_ligand, name_protein, dock, box_size, num_poses, dock_algorithm, random_pos, ga_dock, sa_dock, energy=None):
    # Insert parameters for R and epsilon for H-bond
    subprocess.call([r'''awk '{ if ($2 == "'''+metal_symbol.upper()+'''" || $2 == "'''+metal_symbol+'''") ($7 = '''+str(parameter_set[10])+''') && ($8 = '''+str(parameter_set[11])+'''); print $0}' '''+parameter_file+''' > file_1'''], shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    subprocess.call([r'''awk '{ if ($2 == "'''+metal_symbol.upper()+'''" || $2 == "'''+metal_symbol+r'''") printf"%-8s %-3s %7s %8s %8s %9s %4s %4s %2s %3s %3s %2s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12; else print $0}' file_1 > '''+parameter_file], shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    subprocess.call(['rm file_1'], shell=True)

    #create_gpf():
    subprocess.call([os.environ['PYTHON_2']+" "+os.environ['MGLTOOLS']+"/prepare_gpf4.py -l "+name_ligand+".pdbqt  -r clean_"+name_protein+".pdbqt -p parameter_file="+parameter_file+" -p npts='{},{},{}'".format(box_size,box_size,box_size)+" -p gridcenter='{:.4},{:.4},{:.4}' ".format(dock[0],dock[1],dock[2])], shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    gpf = open('clean_'+name_protein+'.gpf', 'a')
    gpf.write('nbp_r_eps {:>.4f}  {:>.4f}'.format(parameter_set[0],parameter_set[1])+'    12 6 OA '+metal_symbol+'\n')
    gpf.write('nbp_r_eps {:>.4f}  {:>.4f}'.format(parameter_set[2],parameter_set[3])+'    12 6 SA '+metal_symbol+'\n')
    gpf.write('nbp_r_eps {:>.4f}  {:>.4f}'.format(parameter_set[4],parameter_set[5])+'    12 6 HD '+metal_symbol+'\n')
    gpf.write('nbp_r_eps {:>.4f}  {:>.4f}'.format(parameter_set[6],parameter_set[7])+'    12 6 NA '+metal_symbol+'\n')
    gpf.write('nbp_r_eps {:>.4f}  {:>.4f}'.format(parameter_set[8],parameter_set[9])+'    12 6  N '+metal_symbol+'\n')
    gpf.close()

    #autogrid()
    subprocess.call([os.environ['ROOT_DIR']+'/external/AutoDock/autogrid4 -p clean_'+name_protein+'.gpf'], shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    #create_dpf()
    write_dpf_file('clean_'+name_protein+'.gpf', name_ligand, 'clean_'+name_protein, parameter_file, num_poses, dock_algorithm, random_pos=random_pos, SA=sa_dock, GA=ga_dock, energy_ligand=energy)

    #autodock()
    subprocess.call([os.environ['ROOT_DIR']+'/external/AutoDock/autodock4 -p '+name_ligand+'_clean_'+name_protein+'.dpf'], shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    #write_all_conformations()
    subprocess.call([os.environ['PYTHON_2']+" "+os.environ['MGLTOOLS']+"/write_conformations_from_dlg.py -d "+name_ligand+"_clean_"+name_protein+".dlg"], shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    return

def write_dpf_file(gpf_file, name_ligand, name_protein, parameter_file, num_poses, dock_algorithm, random_pos=False, GA=False, SA=False, energy_ligand=None):
    gpf_file = open(gpf_file,'r')
    gpf_lines = [line.split() for line in gpf_file]

    ligand_type = gpf_lines[5]
    del ligand_type[0]
    del ligand_type[-4:]
    ligand_type_str = ' '.join(ligand_type)

    dpf_file = open(name_ligand+'_'+name_protein+'.dpf','w')
    dpf_file.write('autodock_parameter_version 4.2       # used by autodock to validate parameter set\n')
    dpf_file.write('parameter_file '+parameter_file+' # parameter library filename\n')
    dpf_file.write('outlev 1                             # diagnostic output level\n')
    dpf_file.write('intelec                              # calculate internal electrostatics\n')
    dpf_file.write('seed pid time                        # seeds for random generator\n')
    dpf_file.write('ligand_types '+ligand_type_str+'             # atoms types in ligand\n')
    dpf_file.write('fld '+name_protein+'.maps.fld              # grid_data_file\n')
    for i in range(0,len(ligand_type)):
        dpf_file.write('map '+name_protein+'.'+ligand_type[i]+'.map                 # atom-specific affinity map\n')

    dpf_file.write('elecmap '+name_protein+'.e.map             # electrostatics map\n')
    dpf_file.write('desolvmap '+name_protein+'.d.map           # desolvation map\n\n')

    dpf_file.write('# Unbound Ligand Parameters\n')
    if energy_ligand != None:
        dpf_file.write('unbound_energy '+str(energy_ligand)+'              # set the energy of the unbound state\n')
    
    dpf_file.write('move '+name_ligand+'.pdbqt                # small molecule\n')

    if random_pos == True:
        dpf_file.write('tran0 random                         # initial coordinates/A or random\n')
        dpf_file.write('quaternion0 random                   # initial orientation\n')
        dpf_file.write('dihe0 random                         # initial dihedrals (relative) or random\n')

    if GA == True:
        dpf_file.write('# GA parameters\n')
        dpf_file.write('ga_pop_size '+str(dock_algorithm[0])+'                      # number of individuals in population\n')
        dpf_file.write('ga_num_evals '+str(dock_algorithm[1])+'                 # maximum number of energy evaluations\n')
        dpf_file.write('ga_num_generations '+str(dock_algorithm[2])+'             # maximum number of generations\n')
        dpf_file.write('ga_elitism '+str(dock_algorithm[3])+'                         # number of top individuals to survive to next generation\n')
        dpf_file.write('ga_mutation_rate '+str(dock_algorithm[4])+'                 # rate of gene mutation\n')
        dpf_file.write('ga_crossover_rate '+str(dock_algorithm[5])+'                # rate of crossover\n')
        dpf_file.write('ga_window_size '+str(dock_algorithm[6])+'                    # number of preceding generation when deciding threshold for worst individual current population\n')
        dpf_file.write('ga_cauchy_alpha 0.0                  # Alpha parameter of Cauchy distribution\n')
        dpf_file.write('ga_cauchy_beta 1.0                   # Beta parameter Cauchy distribution\n')

        dpf_file.write('# Local Search Parameters\n')
        dpf_file.write('sw_max_its 300                       # iterations of Solis & Wets local search\n')
        dpf_file.write('sw_max_succ 4                        # consecutive successes before changing rho\n')
        dpf_file.write('sw_max_fail 4                        # consecutive failures before changing rho\n')
        dpf_file.write('sw_rho 1.0                           # size of local search space to sample\n')
        dpf_file.write('sw_lb_rho 0.01                       # lower bound on rho\n')
        dpf_file.write('ls_search_freq 0.06                  # probability of performing local search on individual\n')
        # dpf_file.write('do_local_only 20\n')
        dpf_file.write('# Activate LGA\n')
        dpf_file.write('set_ga                               # set the above parameters for GA or LGA\n')
        dpf_file.write('set_psw1                             # set the above pseudo-Solis & Wets parameters\n')
        dpf_file.write('ga_run '+str(num_poses)+'                             # do this many hybrid GA-LS runs\n')
    if SA == True:
        dpf_file.write('# SA Parameters\n')
        dpf_file.write('tstep 2.0\n')
        #dpf_file.write('e0max 0.0 10000                      # max initial energy; max number of retries\n')
        dpf_file.write('linear_schedule                      # linear_schedule or geometric_schedule\n')
        dpf_file.write('rt0 500                              # initial annealing temperature (absolute tmperature multiplied by gas constant\n')
        dpf_file.write('rtrf '+str(dock_algorithm[0])+'           # annealing temperature reductin factor < 1 cools > 1 heats system\n')
        dpf_file.write('runs '+str(dock_algorithm[1])+'           # number of docking runs\n')
        dpf_file.write('cycles '+str(dock_algorithm[2])+'         # number of temperature reduction cycles\n')
        dpf_file.write('accs 30000                           # maximum number of accepted steps per cycle\n')
        dpf_file.write('rejs 30000                           # maximum number of rejected steps per cycle\n')
        dpf_file.write('select m                             # m selects the minimum state, 1 selects the last state during each cycle\n')
        dpf_file.write('trnrf 1.0                            # per cycle reduction factor for translation steps\n')
        dpf_file.write('quarf 1.0                            # per cycle reduction factor for orientation steps\n')
        dpf_file.write('dihrf 1.0                            # per cycle reduction factor for torsional dihedral steps\n')

        dpf_file.write('# Activate SA\n')
        dpf_file.write('simanneal '+str(num_poses)+'                         # run this many SA docking\n')

    # dpf_file.write('# Perform Analysis\n')
    # dpf_file.write('rmsref '+name_ligand+'.pdbqt              # RMSD will be calculated with this file\n')
    # dpf_file.write('rmsmode atype                        # Method to calculate the RMSD\n')
    # dpf_file.write('rmstol 2.0                           # RMSD tolerance\n')
    dpf_file.write('analysis                             # perforem a ranked cluster analysis\n')
    return


def rmsd_func(name_ligand, n_prot, generation, directory, num_gen=None, train=False, standard=False, test=False):
    rmsd_avg = []
    rmsd_list = []
    avg_list = []
    min_list = []
    rmsd_print_list = []
    
    if standard == True:
        output = [f"---------------------------------------     PROTEIN {n_prot} STANDARD     ---------------------------------------\n"]
    if test == True:
        output = [f"---------------------------------------     PROTEIN {n_prot} TEST         ---------------------------------------\n"]
    if train == True:
        output = [f"-------------------------------------------     PROTEIN {n_prot}      --------------------------------------------\n"]

    i = 1
    while os.path.exists(name_ligand+"_%i.pdbqt" % i):
        subprocess.call(os.environ['OBABEL']+" -ipdbqt "+name_ligand+"_{}.pdbqt".format(i)+" -oxyz "+name_ligand+"_{}.xyz".format(i)+" -d > "+name_ligand+"_{}.xyz".format(i), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) 

        rmsd_non_rotate = float(subprocess.getoutput([os.environ['PYTHON_3']+' '+os.environ['ROOT_DIR']+'/metal_dock/calculate_rmsd.py ref.xyz '+name_ligand+'_{}.xyz'.format(i)+' --reorder --rotation none --translation none']))
        rmsd = rmsd_non_rotate

        rmsd_print_list.append("RMSD for Conformation %i = %.4f\n"% (i, rmsd))
        rmsd_list.append(rmsd)
        i += 1
    

    for j in range(0,len(rmsd_print_list)):
        output.append(rmsd_print_list[j])

        if standard == True:
            std = open(f'{directory}/standard/standard_conformations','a')
            std.write(rmsd_print_list[j])

            protein = open(f'{directory}/standard/protein_{n_prot}','a')
            protein.write(rmsd_print_list[j])


        if test == True:   
            t = open(f'{directory}/test/test_conformations','a')
            t.write(rmsd_print_list[j])
            
            protein = open(f'{directory}/test/protein_{n_prot}','a')
            protein.write(rmsd_print_list[j])

        if train == True:
            if generation == 0:
                first_gen = open(f'{directory}/first_gen/all_conf_first_gen','a')
                first_gen.write(rmsd_print_list[j])

                protein = open(f'{directory}/first_gen/protein_{n_prot}','a')
                protein.write(rmsd_print_list[j])

            if generation == num_gen+1:   
                last_gen = open(f'{directory}/last_gen/all_conf_last_gen','a')
                last_gen.write(rmsd_print_list[j])
            
                protein = open(f'{directory}/last_gen/protein_{n_prot}','a')
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


# def distance(x1,x2,y1,y2,z1,z2):

#     d = np.sqrt( (x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2  )
#     d = format(d, '.4f')

#     return str(d)

# # Translation Matrix - Translation is switched from column to row because we work with (1,4) coordinate matrix
# def translation_matrix(matrix,dx,dy,dz):
#     translation_matrix = np.array(([1,0,0,0],
#                                    [0,1,0,0],
#                                    [0,0,1,0],
#                                    [dx,dy,dz,1]),dtype=np.float64)

#     y_dim,x_dim = np.shape(matrix)
#     extra_dim_arr = np.c_[matrix, np.ones(y_dim)]

#     new_extra_dim_arr = np.dot(extra_dim_arr, translation_matrix)
#     new_arr = np.delete(new_extra_dim_arr, (3), 1)

#     return new_arr

# # Rotation Matrices
# def x_axis_rotation(matrix,theta):
#     rotation_matrix = np.array(([1,0,0,0],
#                                 [0,np.cos(theta),-np.sin(theta),0],
#                                 [0,np.sin(theta),np.cos(theta),0],
#                                 [0,0,0,1]), dtype=np.float64)

#     y_dim,x_dim = np.shape(matrix)
#     extra_dim_arr = np.c_[matrix, np.ones(y_dim)]

#     new_extra_dim_arr = np.dot(extra_dim_arr, rotation_matrix)
#     new_arr = np.delete(new_extra_dim_arr, (3), 1)

#     return new_arr

# def y_axis_rotation(matrix,theta):
#     rotation_matrix = np.array(([np.cos(theta),0,np.sin(theta),0],
#                                 [0,1,0,0],
#                                 [-np.sin(theta),0,np.cos(theta),0],
#                                 [0,0,0,1]),dtype=np.float64)

#     y_dim,x_dim = np.shape(matrix)
#     extra_dim_arr = np.c_[matrix, np.ones(y_dim)]

#     new_extra_dim_arr = np.dot(extra_dim_arr, rotation_matrix)
#     new_arr = np.delete(new_extra_dim_arr, (3), 1)

#     return new_arr

# def z_axis_rotation(matrix,theta):
#     rotation_matrix = np.array(([np.cos(theta),-np.sin(theta),0,0],
#                                 [np.sin(theta),np.cos(theta),0,0],
#                                 [0,0,1,0],
#                                 [0,0,0,1]),dtype=np.float64)

#     y_dim,x_dim = np.shape(matrix)
#     extra_dim_arr = np.c_[matrix, np.ones(y_dim)]

#     new_extra_dim_arr = np.dot(extra_dim_arr, rotation_matrix)
#     new_arr = np.delete(new_extra_dim_arr, (3), 1)

#     return new_arr

# def randomize_translation_rotation(pdbqt_file):
#     os.system('''awk '{if ($1 == "ATOM") print $0}' '''+pdbqt_file+''' > temp_1''')

#     dock_site = open('temp_1','r')
#     coord = [line.split() for line in dock_site]

#     for i in range(0,len(coord)):
#         del coord[i][0:6]
#         del coord[i][3:7]

#     coord_float = [[float(j) for j in i] for i in coord]
#     coord_array = np.array(coord_float)

#     # Random Coordinates
#     random.seed()

#     theta = random.uniform(0, 2*np.pi) # rotation angle

#     box_size = list(iv.var.box_size.split(','))
#     box_size = [float(x) for x in box_size]
#     box_size = [(x*0.375) / 2 for x in box_size]

#     dx = random.uniform(-box_size[0], box_size[0]) # translation x-direction
#     dy = random.uniform(-box_size[1], box_size[1]) # translation y-direction
#     dz = random.uniform(-box_size[2], box_size[2]) # translation z-direction

#     random_coord = translation_matrix(coord_array,dx,dy,dz)
#     random_coord = x_axis_rotation(random_coord,theta)
#     random_coord = y_axis_rotation(random_coord,theta)
#     random_coord = z_axis_rotation(random_coord,theta)

#     output_array = [["%.3f" % j for j in i] for i in random_coord]
#     output_array = [[str(j) for j in i] for i in output_array]

#     output = open('output_test','w')

#     for elem in output_array:
#         output.write('\t'.join(elem))
#         output.write('\n')

#     output.close()

#     os.system("paste -d' 'temp_1 temp_1 output_test> temp_2")
#     os.system(r'''awk '{ printf "%-4s %6s %2s %5s %1s %3s %11s %7s %7s %5s %5s %9s %-3s\n",$1,$2,$3,$4,$5,$6,$14,$15,$16,$10,$11,$12,$13}' temp_2 > temp_3''')

#     new_file = open('temp_3','r')
#     new_file_coord = [line.split() for line in new_file]

#     old_file = open(''+pdbqt_file+'','r')
#     old_file_coord = [line.split() for line in old_file]

#     new_pdbqt = []
#     k = 0

#     for i in range(0,len(old_file_coord)):
#         for j in range(0,len(new_file_coord)):
#             if old_file_coord[i][0] == new_file_coord[j][0]:
#                 new_pdbqt.append(new_file_coord[k])
#                 k+=1
#                 break
#             else:
#                 new_pdbqt.append(old_file_coord[i])
#                 break

#     new_output = open('new_output_test','w')

#     for elem in new_pdbqt:
#         new_output.write('\t'.join(elem))
#         new_output.write('\n')

#     new_output.close()
#     os.system(r'''awk '{ printf "%-4s %5s %4s %4s %1s %3s %11s %7s %7s %5s %5s %9s %-3s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' new_output_test > '''+iv.var.name_ligand+'''.pdbqt''')

#     os.system("rm temp_1 temp_2 temp_3 output_test new_output_test")


