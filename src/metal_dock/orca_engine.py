import os, sys, subprocess
from ase.calculators.orca import ORCA
from ase.io import read, write

def orca_engine(xyz_file, var, output_dir):
    ## Geometry Optimization ##
    if var.geom_opt == True:
        if os.path.isdir('geom_opt') == False:
            os.mkdir('geom_opt')
            os.chdir('geom_opt')
        else:
            os.chdir('geom_opt')

        subprocess.call([f'cp {output_dir}/file_prep/'+xyz_file+' .'], shell=True)

        # If Geometry Converged Skip otherwise Run Again#
        if os.path.exists(f'{output_dir}/QM/geom_opt/geom.out') == False:
            orca_geom_opt(xyz_file, var)
            orca_opt_converged('geom.out')
            orca_extract_CM5('geom.out', xyz_file)
            energy = orca_extract_energy('geom.out')
        
        else:
            orca_opt_converged('geom.out')
            orca_extract_CM5('geom.out', xyz_file)
            energy = orca_extract_energy('geom.out')

    else:
        ## Single Point ##
        os.chdir(f'{output_dir}/QM')

        if os.path.isdir('single_point') == False:
            os.mkdir('single_point')
            os.chdir('single_point')
        else:
            os.chdir('single_point')

        subprocess.call([f'cp {output_dir}/file_prep/'+xyz_file+' .'], shell=True)

        # If Geometry Converged Skip otherwise Run Again#
        if os.path.exists(f'{output_dir}/QM/single_point/single_point.out') == False:
            orca_single_point(xyz_file, var)
            orca_sp_converged('single_point.out')
            orca_extract_CM5('single_point.out', xyz_file)
            energy = orca_extract_energy('single_point.out')

        if orca_sp_converged('single_point.out') == True:
            orca_extract_CM5('single_point.out', xyz_file)
            energy = orca_extract_energy('single_point.out')

    return os.getcwd(), energy 


def orca_extract_energy(log_file):
    with open(log_file,'r') as fin:
        for line in fin:
            if line.startswith('FINAL'):
                energy = line.split()[4]
                return energy
                
def orca_extract_CM5(log_file, xyz_file):
    mol = read(xyz_file)
    N = len(mol.positions) 

    subprocess.call(["grep -A"+str(N+6)+" 'HIRSHFELD ANALYSIS' "+log_file+" > charge_1"], shell=True)
    with open('charge_1','r') as fin:
        with open('CM5_charges','w') as fout:
            fin_lines = [line.split() for line in fin]
            fout.write('\n')
            for i in fin_lines[7:]:
                fout.write('{} {}\n'.format(i[1],i[2]))

    subprocess.call(['rm charge_1'], shell=True)
    return

def orca_opt_converged(log_file):
     with open(log_file) as log:
        if 'SUCCESS' in log.read():
            print('GEOMETRY CONVERGED')
            return True
        else:
            print('GEOMETRY NOT CONVERGED - DELETE geom.out')
            return sys.exit()

def orca_sp_converged(log_file):
     with open(log_file) as log:
        if 'SUCCESS' in log.read():
            print('GEOMETRY CONVERGED')
            return True
        else:
            print('GEOMETRY NOT CONVERGED - DELETE geom.out')
            return sys.exit()

def orca_geom_opt(xyz_file, var):
    M = 2 * var.spin + 1

    mol = read(xyz_file)
    mol.calc = ORCA(label='geom',
                    charge=var.charge,
                    mult=M,
                    orcasimpleinput='Opt '+var.functional+' '+var.basis_set+' '+var.dispersion+' CPCM(Water)',
                    orcablocks='%pal nprocs '+var.ncpu+' end % output Print[P_hirshfeld] 1 end',
                    )

    mol.get_potential_energy()
    mol.write('output.xyz')
    return

def orca_single_point(xyz_file, var):
    M = 2 * var.spin + 1

    mol = read(xyz_file)
    mol.calc = ORCA(label='single_point',
                    charge=var.charge,
                    mult=M,
                    orcasimpleinput=var.functional+' '+var.basis_set+' '+var.dispersion+' CPCM(Water)',
                    orcablocks='%pal nprocs '+var.ncpu+' end % output Print[P_hirshfeld] 1 end',
                    )

    mol.get_potential_energy()
    return