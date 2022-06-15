from scm.plams import *
import input_variables as iv
import sys

def run_gfnxtb(xyz_file):
    init()

    #Molecule Structure
    m = Molecule(xyz_file)
    m.properties.charge = ''+iv.var.charge_ligand+''

    #Settings
    s = Settings()
    #AMS driver input
    s.input.ams.Task = 'GeometryOptimization'
    s.input.ams.Properties.NormalModes = 'Yes'

    #DFTB engine inpu
    s.input.DFTB.Model = 'GFN1-xTB'
    s.input.DFTB.DispersionCorrection = 'D3-BJ'
    s.input.DFTB.Solvation.Solvent = 'H2O'

    if iv.var.spin_ligand is not None:
        s.input.DFTB.Occupation.Strategy ='Aufbau'
        s.input.DFTB.UnpairedElectrons=''+iv.var.spin_ligand+''
    else:
        s.input.DFTB.Occupation.Strategy ='Auto'

    #Run Job
    j = AMSJob(molecule=m, settings=s)
    result = j.run()

    finish()

def gfnxtb_converged(ams_log):
    with open(ams_log) as log:
        if 'Geometry optimization successful!' in log.read():
            print('Geometry converged')
        else:
            print('Geometry not converged')
            sys.exit()

def plams_single_point(xyz_file):
    init()
    #Molecule Struture
    m = Molecule(xyz_file)
    m.properties.charge = ''+iv.var.charge_ligand+''

    #Settings
    s = Settings()
    #AMS driver input
    s.input.ams.Task='SinglePoint'

    #ADF engine input
    s.input.adf.AtomicChargesTypeForAMS='CM5'
    s.input.adf.basis.type=''+iv.var.basis_set.upper()+''
    s.input.adf.basis.core='None'

    if iv.var.gga_functional is not None:
        s.input.adf.xc.gga=''+iv.var.gga_functional.upper()+''

    if iv.var.hybrid_functional is not None:
        s.input.adf.xc.hybrid=''+iv.var.hybrid_functional.upper()+''

    if iv.var.dispersion_correction is not None:
        s.input.adf.xc.dispersion=''+iv.var.dispersion_correction.upper()+''

    if iv.var.spin_ligand is not None:
        s.input.adf.unrestricted='yes'
        s.input.adf.spinpolarization=''+iv.var.spin_ligand+''

    s.input.adf.relativity.formalism='ZORA'
    s.input.adf.relativity.level='Scalar'

    #Run Job
    j = AMSJob(molecule=m, settings=s)
    result = j.run()

    finish()

def single_point_check(ams_log):
    with open(ams_log) as log:
        if 'WARNING: partially occupied orbitals' in log.read():
            print('Single point not successful')
            sys.exit()
        else:
            print('Single point succesfully performed')
