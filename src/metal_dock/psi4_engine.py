import os, sys, subprocess
import psi4
import multiprocessing

def psi4_engine(xyz_file, var, output_dir):
    ## Geometry Optimization ##
    if os.getenv('SLURM_STEP_TASKS_PER_NODE')!= None:
        n_threads = int(os.getenv('SLURM_STEP_TASKS_PER_NODE'))
    elif os.getenv('PBS_NP') != None:
        n_threads = int(os.getenv('PBS_NP'))
    elif os.getenv('NSLOTS') != None:
        n_threads = int(os.getenv('NSLOTS'))
    else:
        n_threads = multiprocessing.cpu_count()
        
    print(n_threads)

    psi4.set_num_threads(n_threads)
    if var.geom_opt == True:
        if os.path.isdir('geom_opt') == False:
            os.mkdir('geom_opt')
            os.chdir('geom_opt')
        else:
            os.chdir('geom_opt')

        subprocess.call([f'cp {output_dir}/file_prep/'+xyz_file+' .'], shell=True)

        # If Geometry Converged Skip otherwise Run Again#
        #if os.path.exists(f'{output_dir}/QM/geom_opt/geom_opt.chk') == False:
        psi4_geom_opt(xyz_file, var)
        # gaussian_opt_converged('geom_opt.log')
        


def psi4_geom_opt(xyz_file, var):

    M = 2 * var.spin + 1 

    with open(xyz_file, 'r') as fin:
        xyz = fin.read()

        input_xyz = xyz.split('\n')[2:]
        insert =  str(var.charge)+' '+str(M)
        input_xyz.insert(0, insert)
        input_xyz = '\n'.join(input_xyz)

        mol = psi4.geometry(input_xyz)

        # mol = psi4.core.Molecule.from_string(xyz, dtype='xyz')

        if var.spin == 0:
            psi4.set_options({'reference' : 'RKS',
                            'basis' : var.basis_set,
                            })
            psi4.optimize('B3LYP-D3BJ', molecule=mol)

        if var.spin != 0:
            psi4.set_options({'reference' : 'UKS',
                            'basis' : var.basis_set,
                            })
            psi4.optimize('B3LYP-D3BJ', molecule=mol)

        

    mol.save_xyz_file('name.xyz',1)


