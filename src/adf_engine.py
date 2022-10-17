import os, sys, subprocess


def adf_engine(xyz_file, var, output_dir):

    from scm import plams

    ## Geometry Optimization ##
    if var.geom_opt == True:
        if os.path.isdir('geom_opt') == False:
            os.mkdir('geom_opt')
            os.chdir('geom_opt')
        else:
            os.chdir('geom_opt')

        subprocess.call([f'cp {output_dir}/file_prep/'+xyz_file+' .'], shell=True)

        # If Geometry Converged Skip otherwise Run Again#
        if os.path.isdir(f'{output_dir}/QM/geom_opt/plams_workdir/plamsjob') == False:
            adf_geom_opt(xyz_file, var)
            os.chdir(f'{output_dir}/QM/geom_opt/plams_workdir/plamsjob')
            geometry_optimized_check('ams.log')
        else:
            os.chdir(f'{output_dir}/QM/geom_opt/plams_workdir/plamsjob')
            geometry_optimized_check('ams.log')

    ## Single Point ##
    os.chdir(f'{output_dir}/QM')

    if os.path.isdir('single_point') == False:
        os.mkdir('single_point')
        os.chdir('single_point')
    else:
        os.chdir('single_point')

    if var.geom_opt == True:
        subprocess.call([f'cp {output_dir}/QM/geom_opt/plams_workdir/plamsjob/output.xyz '+xyz_file], shell=True)
    else:
        subprocess.call([f'cp {output_dir}/file_prep/'+xyz_file+' .'], shell=True)

    # If single point successful Skip otherwise Run Again#
    if os.path.isdir(f'{output_dir}/QM/single_point/plams_workdir/plamsjob') == False:
        adf_sp(xyz_file, var)
        os.chdir(f'{output_dir}/QM/single_point/plams_workdir/plamsjob')
        single_point_check('ams.log')
        subprocess.call([os.environ['AMSBIN']+'/amsreport adf.rkf CM5 > CM5_charges'], shell=True)
        subprocess.call(["grep 'kcal/mol' ams.log > energy"], shell=True)
    else:
        os.chdir(f'{output_dir}/QM/single_point/plams_workdir/plamsjob')
        single_point_check('ams.log')
        subprocess.call([os.environ['AMSBIN']+'/amsreport adf.rkf CM5 > CM5_charges'], shell=True)
        subprocess.call(["grep 'kcal/mol' ams.log > energy"], shell=True)

    return os.getcwd()


def adf_geom_opt(xyz_file, var):
    init()

    #Molecule Structure
    m = Molecule(xyz_file)
    m.properties.charge = ''+str(var.charge)+''

    #Settings
    s = Settings()
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
    j = AMSJob(molecule=m, settings=s)
    result = j.run()

    finish()

def geometry_optimized_check(ams_log):
    with open(ams_log) as log:
        if 'Geometry optimization converged' in log.read():
            return print('\nGEOMETRY CONVERGED\n')
        else:
            print('\nGEOMETRY NOT CONVERGED\n')
            return sys.exit()

def adf_sp(xyz_file, var):
    init()
    #Molecule Struture
    m = Molecule(xyz_file)
    m.properties.charge = ''+str(var.charge)+''

    #Settings
    s = Settings()
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
    j = AMSJob(molecule=m, settings=s)
    result = j.run()

    finish()

def single_point_check(ams_log):
    with open(ams_log) as log:
        if 'WARNING: partially occupied orbitals' in log.read() or 'ERROR' in log.read():
            print('\nSINGLE POINT NOT SUCCESSFUL\n')
            return sys.exit()
        else:
            return print('\nSINGLE POINT SUCCESSFULLY PERFORMED\n')
