import os,sys, subprocess
import multiprocessing


from ase.io import read
from ase.calculators.gaussian import Gaussian, GaussianOptimizer


def gaussian_engine(xyz_file, var, output_dir):
    ## Geometry Optimization ##
    if var.geom_opt == True:
        if os.path.isdir('geom_opt') == False:
            os.mkdir('geom_opt')
            os.chdir('geom_opt')
        else:
            os.chdir('geom_opt')

        subprocess.call([f'cp {output_dir}/file_prep/'+xyz_file+' .'], shell=True)

        # If Geometry Converged Skip otherwise Run Again#
        if os.path.exists(f'{output_dir}/QM/geom_opt/geom_opt.chk') == False:
            gaussian_geom_opt(xyz_file, var)
            gaussian_opt_converged('geom_opt.log')
        
        if gaussian_opt_converged('geom_opt.log') == True:
            pass

        extract_xyz(xyz_file)

        ## Single Point ##
        os.chdir(f'{output_dir}/QM')

        if os.path.isdir('single_point') == False:
            os.mkdir('single_point')
            os.chdir('single_point')
        else:
            os.chdir('single_point')

        subprocess.call([f'cp {output_dir}/QM/geom_opt/output.xyz .'], shell=True)

        if os.path.exists(f'{output_dir}/QM/single_point/single_point.chk') == False:
            gaussian_sp('output.xyz', var)
            extract_CM5('single_point.log', 'output.xyz')
            energy = gaussian_extract_energy('single_point.log')
        
        if gaussian_sp_converged('single_point.log') == True:
            extract_CM5('single_point.log', 'output.xyz')
            energy = gaussian_extract_energy('single_point.log')
        
    else:
        ## Single Point ##
        os.chdir(f'{output_dir}/QM')

        if os.path.isdir('single_point') == False:
            os.mkdir('single_point')
            os.chdir('single_point')
        else:
            os.chdir('single_point')

        subprocess.call([f'cp {output_dir}/file_prep/'+xyz_file+' output.xyz'], shell=True)
        
        if os.path.exists(f'{output_dir}/QM/single_point/single_point.chk') == False:
            gaussian_sp('output.xyz', var)
            gaussian_sp_converged('single_point.log')
            gaussian_extract_CM5('single_point.log', 'output.xyz')
            energy = gaussian_extract_energy('single_point.log')
        
        if gaussian_sp_converged('single_point.log') == True:
            extract_CM5('single_point.log', 'output.xyz')
            energy = gaussian_extract_energy('single_point.log')

    return os.getcwd(), energy       

    

def gaussian_extract_energy(log_file):
    with open(log_file,'r') as fin:
        for line in fin:
            if line.startswith('SCF Done'):
                energy = line.split()[3]
                return energy


def gaussian_extract_CM5(log_file, xyz_file):
    mol = read(xyz_file)
    N = len(mol.positions) 

    subprocess.call(["grep -A"+str(N+1)+" 'Hirshfeld charges, spin densities, dipoles, and CM5 charges' "+log_file+" > charge_1"], shell=True)
    with open('charge_1','r') as fin:
        with open('CM5_charges','w') as fout:
            fin_lines = [line.split() for line in fin]

            for i in fin_lines[2:]:
                fout.write('{} {}\n'.format(i[1],i[6]))

    subprocess.call(['rm charge_1'], shell=True)
    return

            
def gaussian_opt_converged(log_file):
     with open(log_file) as log:
        if 'Optimization completed.' in log.read():
            print('GEOMETRY CONVERGED')
            return True
        else:
            print('GEOMETRY NOT CONVERGED - DELETE .chk AND RERUN')
            return sys.exit()

def gaussian_sp_converged(log_file):
     with open(log_file) as log:
        if 'SCF Done' in log.read():
            print('\nSINGLE POINT SUCCESSFULLY PERFORMED\n')
            return True
        else:
            print('\nSINGLE POINT NOT SUCCESSFUL\n')
            return sys.exit()

def gaussian_geom_opt(xyz_file, var):
    n_procs = var.ncpu 
    M = 2 * var.spin + 1 

    mol = read(xyz_file)
    # symbols = list(mol.symbols)

    s   = Gaussian(label='geom_opt',
                    nprocshared=n_procs,
                    mem='4GB',
                    chk='geom_opt.chk',
                    xc=var.functional,
                    charge=var.charge,
                    mult=M,
                    basis=var.basis_set,
                    pop='Hirshfeld',
                    SCRF='Solvent=Water',
                    EmpiricalDispersion=var.dispersion)

    opt = GaussianOptimizer(mol, s)
    opt.run(fmax='tight')
    opt.write('output.xyz')
    return 


def gaussian_sp(xyz_file, var):
    n_procs = var.ncpu
    M = 2 * var.spin + 1 

    mol = read(xyz_file)
    mol.calc = Gaussian(label='single_point',
                        nprocshared=n_procs,
                        mem='4GB',
                        chk='single_point.chk',
                        xc=var.functional,
                        charge=var.charge,
                        mult=M,
                        basis=var.basis_set,
                        pop='Hirshfeld',
                        SCRF='Solvent=Water',
                        scf='tight')

    mol.get_potential_energy()
    return