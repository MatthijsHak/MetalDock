import os, sys, subprocess
import scm.plams as scm

def adf_engine(xyz_file, var, output_dir):
    ## Geometry Optimization ##
    if var.geom_opt == True:
        if os.path.isdir('geom_opt') == False:
            os.mkdir('geom_opt')
            os.chdir('geom_opt')
        else:
            os.chdir('geom_opt')

        subprocess.call([f'cp '+xyz_file+' .'], shell=True)

        # If Geometry Converged Skip otherwise Run Again#
        if os.path.isdir(f'{output_dir}/QM/geom_opt/plams_workdir/plamsjob') == False:
            adf_geom_opt(xyz_file, var)
            os.chdir(f'{output_dir}/QM/geom_opt/plams_workdir/plamsjob')
            adf_opt_converged('ams.log')
            subprocess.call([os.environ['AMSBIN']+'/amsreport adf.rkf CM5 > CM5_charges'], shell=True)
            energy = adf_extract_energy('ams.log')
        else:
            os.chdir(f'{output_dir}/QM/geom_opt/plams_workdir/plamsjob')
            adf_opt_converged('ams.log')
            subprocess.call([os.environ['AMSBIN']+'/amsreport adf.rkf CM5 > CM5_charges'], shell=True)
            energy = adf_extract_energy('ams.log')

    else:
        
        if os.path.isdir('single_point') == False:
            os.mkdir('single_point')
            os.chdir('single_point')
        else:
            os.chdir('single_point')


        # If single point successful Skip otherwise Run Again#
        if os.path.isdir(f'{output_dir}/QM/single_point/plams_workdir/plamsjob') == False:
            adf_sp(xyz_file, var)
            os.chdir(f'{output_dir}/QM/single_point/plams_workdir/plamsjob')
            subprocess.call([os.environ['AMSBIN']+'/amsreport adf.rkf CM5 > CM5_charges'], shell=True)
            energy = adf_extract_energy('ams.log')
        else:
            os.chdir(f'{output_dir}/QM/single_point/plams_workdir/plamsjob')
            subprocess.call([os.environ['AMSBIN']+'/amsreport adf.rkf CM5 > CM5_charges'], shell=True)
            energy = adf_extract_energy('ams.log')

    return os.getcwd(), energy


def adf_extract_energy(log_file):
    with open(log_file,'r') as fin:
        for line in fin:
            if 'kcal/mol' in line:
                energy = line.split()[4]
                return energy


def adf_opt_converged(ams_log):
    with open(ams_log) as log:
        if 'Geometry optimization converged' in log.read():
            return print('\nGEOMETRY CONVERGED\n')
        else:
            print('\nGEOMETRY NOT CONVERGED\n')
            return sys.exit()

def adf_sp_converged(ams_log):
    with open(ams_log) as log:
        if 'WARNING: partially occupied orbitals' in log.read() or 'ERROR' in log.read():
            print('\nSINGLE POINT NOT SUCCESSFUL\n')
            return sys.exit()
        else:
            return print('\nSINGLE POINT SUCCESSFULLY PERFORMED\n')

def adf_geom_opt(xyz_file, var):
    scm.init()

    #Molecule Structure
    m = scm.Molecule(xyz_file)
    m.properties.charge = ''+str(var.charge)+''

    #Settings
    s = scm.Settings()
    #AMS driver input
    s.input.ams.Task = 'GeometryOptimization'

    #ADF engine input
    s.input.adf.scf.iterations='500'
    # s.input.adf.AtomicChargesTypeForAMS='CM5'
    s.input.adf.basis.type=''+var.basis_set.upper()+''
    s.input.adf.basis.core='None'

    if var.functional_type.lower() == 'gga':
        s.input.adf.xc.gga=''+var.functional.upper()+''

    if var.functional_type.lower() == 'hybrid':
        s.input.adf.xc.hybrid=''+var.functional.upper()+''

    if var.dispersion != None:
        s.input.adf.xc.dispersion=''+var.dispersion.upper()+''

    if var.spin != 0:
        s.input.adf.unrestricted='yes'
        s.input.adf.spinpolarization=''+str(var.spin)+''

    s.input.adf.relativity.formalism='ZORA'
    s.input.adf.relativity.level='Scalar'

    s.input.adf.Solvation.Solv = "Name=Water"

    #Run Job
    j = scm.AMSJob(molecule=m, settings=s)
    result = j.run()

    scm.finish()

def adf_sp(xyz_file, var):
    scm.init()
    #Molecule Struture
    m = scm.Molecule(xyz_file)
    m.properties.charge = ''+str(var.charge)+''

    #Settings
    s = scm.Settings()
    #AMS driver input
    s.input.ams.Task='SinglePoint'

    #ADF engine input
    s.input.adf.scf.iterations='500'
    s.input.adf.AtomicChargesTypeForAMS='CM5'
    s.input.adf.basis.type=''+var.basis_set.upper()+''
    s.input.adf.basis.core='None'

    if var.functional_type.lower() == 'gga':
        s.input.adf.xc.gga=''+var.functional.upper()+''

    if var.functional_type.lower() == 'hybrid':
        s.input.adf.xc.hybrid=''+var.functional.upper()+''

    if var.dispersion != None:
        s.input.adf.xc.dispersion=''+var.dispersion.upper()+''

    if var.spin != 0:
        s.input.adf.unrestricted='yes'
        s.input.adf.spinpolarization=''+str(var.spin)+''

    s.input.adf.relativity.formalism='ZORA'
    s.input.adf.relativity.level='Scalar'

    s.input.adf.Solvation.Solv = "Name=Water"

    #Run Job
    j = scm.AMSJob(molecule=m, settings=s)
    result = j.run()

    scm.finish()