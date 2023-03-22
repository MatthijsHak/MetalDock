import os, sys, shutil, subprocess
import multiprocessing
import itertools

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

        shutil.copyfile(xyz_file, os.getcwd()+f'/{var.name_ligand}_c.xyz')

        # If Geometry Converged Skip otherwise Run Again#
        if os.path.exists(f'{output_dir}/QM/geom_opt/geom_opt.chk') == False:
            gaussian_geom_opt(xyz_file, var)
            gaussian_opt_converged('geom_opt.log')
        
        else: 
            gaussian_opt_converged('geom_opt.log') 

        ## Single Point ##
        os.chdir(f'{output_dir}/QM')

        if os.path.isdir('single_point') == False:
            os.mkdir('single_point')
            os.chdir('single_point')
        else:
            os.chdir('single_point')

        shutil.copyfile(f'{output_dir}/QM/geom_opt/output.xyz', os.getcwd()+'/output.xyz')

        if os.path.exists(f'{output_dir}/QM/single_point/single_point.chk') == False:
            gaussian_sp('output.xyz', var)
            gaussian_extract_CM5('single_point.log', 'output.xyz')
            energy = gaussian_extract_energy('single_point.log')
        
        if gaussian_sp_converged('single_point.log') == True:
            gaussian_extract_CM5('single_point.log', 'output.xyz')
            energy = gaussian_extract_energy('single_point.log')
        
    else:
        ## Single Point ##
        os.chdir(f'{output_dir}/QM')

        if os.path.isdir('single_point') == False:
            os.mkdir('single_point')
            os.chdir('single_point')
        else:
            os.chdir('single_point')

        shutil.copyfile(xyz_file, os.getcwd()+'/output.xyz')
        
        if os.path.exists(f'{output_dir}/QM/single_point/single_point.chk') == False:
            gaussian_sp('output.xyz', var)
            gaussian_sp_converged('single_point.log')
            gaussian_extract_CM5('single_point.log', 'output.xyz')
            energy = gaussian_extract_energy('single_point.log')
        
        if gaussian_sp_converged('single_point.log') == True:
            gaussian_extract_CM5('single_point.log', 'output.xyz')
            energy = gaussian_extract_energy('single_point.log')

    return os.getcwd(), energy       

    

def gaussian_extract_energy(log_file):
    with open(log_file,'r') as fin:
        for line in fin:
            if line.startswith(' SCF Done:') == True:
                energy = line.split()[4]
                return energy


def gaussian_extract_CM5(log_file, xyz_file):
    mol = read(xyz_file)
    N = len(mol.positions) 

    with open(log_file) as log:
        for line in log:
            if ' Hirshfeld charges, spin densities, dipoles, and CM5 charges' in line:
                fin_lines = list(itertools.islice(log, N+1))
                fin_lines = [line.strip().split() for line in fin_lines]
                with open('CM5_charges','w') as fout:
                    for i in fin_lines[1:]:
                        fout.write('{} {}\n'.format(i[1],i[7]))

    return
            
def gaussian_opt_converged(log_file):
     with open(log_file) as log:
        if 'Optimization completed.' in log.read():
            print('GEOMETRY CONVERGED')
            return True
        else:
            print('GEOMETRY NOT CONVERGED - DELETE .chk, .com & .log FILES TO RERUN')
            return sys.exit()

def gaussian_sp_converged(log_file):
     with open(log_file) as log:
        if 'SCF Done' in log.read():
            print('\nSINGLE POINT SUCCESSFULLY PERFORMED\n')
            return True
        else:
            print('\nSINGLE POINT NOT SUCCESSFUL - DELETE .chk, .com & .log FILES TO RERUN\n')
            return sys.exit()

def gaussian_geom_opt(xyz_file, var):
    M = 2 * var.spin + 1 

    mol = read(xyz_file)

    if var.solvent == '':
        s   = Gaussian(label='geom_opt',
                        nprocshared=var.ncpu ,
                        mem='4GB',
                        chk='geom_opt.chk',
                        xc=var.functional,
                        charge=var.charge,
                        mult=M,
                        basis=var.basis_set,
                        pop='Hirshfeld',
                        EmpiricalDispersion=var.dispersion,
                        integral='dkh')
    else:
        s   = Gaussian(label='geom_opt',
                nprocshared=var.ncpu ,
                mem='4GB',
                chk='geom_opt.chk',
                xc=var.functional,
                charge=var.charge,
                mult=M,
                basis=var.basis_set,
                pop='Hirshfeld',
                SCRF=f'PCM, solvent={var.solvent}',
                EmpiricalDispersion=var.dispersion,
                integral='dkh')

    opt = GaussianOptimizer(mol, s)
    opt.run(fmax='tight')
    mol.write('output.xyz')
    return 


def gaussian_sp(xyz_file, var):
    M = 2 * var.spin + 1 

    mol = read(xyz_file)
    if var.solvent == '':
        mol.calc = Gaussian(label='single_point',
                            nprocshared= var.ncpu,
                            mem='4GB',
                            chk='single_point.chk',
                            xc=var.functional,
                            charge=var.charge,
                            mult=M,
                            basis=var.basis_set,
                            pop='Hirshfeld',
                            scf='tight',
                            integral='dkh')
    else:
        mol.calc = Gaussian(label='single_point',
                    nprocshared= var.ncpu,
                    mem='4GB',
                    chk='single_point.chk',
                    xc=var.functional,
                    charge=var.charge,
                    mult=M,
                    basis=var.basis_set,
                    pop='Hirshfeld',
                    SCRF=f'PCM, solvent={var.solvent}',
                    scf='tight',
                    integral='dkh')

    mol.get_potential_energy()
    return