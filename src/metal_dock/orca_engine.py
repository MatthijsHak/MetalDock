import os, sys, csv, subprocess, shutil

from . import orca2CM5 as oc

from ase.calculators.orca import ORCA
from ase.io import read, write
from ase.optimize.lbfgs import LBFGS

def orca_engine(xyz_file, var, output_dir):
    ## Geometry Optimization ##
    if var.geom_opt == True:
        if os.path.isdir('geom_opt') == False:
            os.mkdir('geom_opt')
            os.chdir('geom_opt')
        else:
            os.chdir('geom_opt')

        shutil.copyfile(xyz_file, os.path.join(os.getcwd(),f'{var.name_ligand}_c.xyz'))

        # If Geometry Converged Skip otherwise Run Again#
        out_file = os.path.join(output_dir,'QM','geom_opt','geom.out')
        if os.path.exists(out_file) == False:
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
        os.chdir(os.path.join(output_dir,'QM'))

        if os.path.isdir('single_point') == False:
            os.mkdir('single_point')
            os.chdir('single_point')
        else:
            os.chdir('single_point')

        shutil.copyfile(xyz_file, os.path.join(os.getcwd(),f'{var.name_ligand}_c.xyz'))

        # If Geometry Converged Skip otherwise Run Again#
        out_file = os.path.join(output_dir,'QM','single_point','single_point.out')
        if os.path.exists(out_file) == False:
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
    # mol = read(xyz_file)
    # N = len(mol.positions) 

    a0,rd,pt = oc.LoadModel() 
    data = oc.GetLogFile(log_file, pt, rd)
    qcm5 = oc.HirshfeldToCM5(xyz_file, data, a0)
    qcm5.to_csv('CM5_charges.csv',index=False,float_format='%6.4f')

    with open('CM5_charges.csv', newline='') as f:
        reader = csv.reader(f)
        #! delete underneath
        data = list(reader)

        with open('CM5_charges','w') as fout:
            fout.write('\n')
            for i in data[1:]:
                fout.write('{} {}\n'.format(i[0],i[8]))
                
    return

def orca_opt_converged(log_file):
     with open(log_file) as log:
        if 'SUCCESS' in log.read():
            print('GEOMETRY CONVERGED')
            return True
        else:
            print('GEOMETRY NOT CONVERGED - VERIFY THE PROBLEM IN geom.out AND DELETE FILE BEFORE RUNNING AGAIN')
            return sys.exit()

def orca_sp_converged(log_file):
     with open(log_file) as log:
        if 'SUCCESS' in log.read():
            print('SINGLE POINT CONVERGED')
            return True
        else:
            print('SINGLE POINT NOT CONVERGED - VERIFY THE PROBLEM IN single_point.out AND DELETE FILE BEFORE RUNNING AGAIN')
            return sys.exit()

def orca_geom_opt(xyz_file, var):
    M = 2 * (var.spin*0.5) + 1

    mol = read(xyz_file)
    mol.calc = ORCA(label='geom',
                    charge=var.charge,
                    mult=M,
                    orcasimpleinput=f'{var.orcasimpleinput}',
                    orcablocks=f'%maxcore {str(var.memory)} %pal nprocs {str(var.ncpu)} end %output Print[P_hirshfeld] 1 end {var.orcablocks}',
                    )
    
    opt = LBFGS(mol)
    opt.run(fmax=0.05)

    mol.get_potential_energy()
    mol.write('output.xyz')

def orca_single_point(xyz_file, var):
    M = 2 * (var.spin*0.5) + 1

    mol = read(xyz_file)
    mol.calc = ORCA(label='single_point',
                    charge=var.charge,
                    mult=M,
                    orcasimpleinput=f'{var.orcasimpleinput}',
                    orcablocks=f'%pal nprocs {str(var.ncpu)} end %output Print[P_hirshfeld] 1 end {var.orcablocks}',
                    )

    mol.get_potential_energy()
    mol.write('output.xyz')

# def optimize_hydrogens(xyz_file, var):
#     M = 2 * (var.spin*0.5) + 1

#     mol = read(xyz_file)
#     mol.calc = ORCA(label='hydrogen_opt',
#                     charge=var.charge,
#                     mult=M,
#                     orcasimpleinput=f'Opt {var.orcasimpleinput}',
#                     orcablocks=f'%maxcore {str(var.memory)} %pal nprocs {str(var.ncpu)} end %geom optimizehydrogens true end')

#     mol.get_potential_energy()
#     mol.write(xyz_file)