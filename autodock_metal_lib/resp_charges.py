import numpy as np
import os, stat
import environment_variables

import subprocess
from subprocess import call

from vdw_surface import vdw_surface
from espfit import fit


def resp_charges(xyz_file):
    options = {}

    # Check options
     # RESP options have large case keys
    options = {k.upper(): v for k, v in sorted(options.items())}

    # VDW surface options
    if 'ESP' not in options:
        options['ESP'] = []
    if 'GRID' not in options:
        options['GRID'] = []
    if 'VDW_SCALE_FACTORS' not in options:
        options['VDW_SCALE_FACTORS'] = [1.4, 1.6, 1.8, 2.0]
    if 'VDW_POINT_DENSITY' not in options:
        options['VDW_POINT_DENSITY'] = 1.0

    # Hyperbolic restraint options
    if 'WEIGHT' not in options:
        options['WEIGHT'] = [1]*1
    if 'RESTRAINT' not in options:
        options['RESTRAINT'] = True
    if options['RESTRAINT']:
        if 'RESP_A' not in options:
            options['RESP_A'] = 0.0005
        if 'RESP_B' not in options:
            options['RESP_B'] = 0.1
        if 'IHFREE' not in options:
            options['IHFREE'] = True
        if 'TOLER' not in options:
            options['TOLER'] = 1e-5
        if 'MAX_IT' not in options:
            options['MAX_IT'] = 25

    # Constraint options
    if 'CONSTRAINT_CHARGE' not in options:
        options['CONSTRAINT_CHARGE'] = []
    if 'CONSTRAINT_GROUP' not in options:
        options['CONSTRAINT_GROUP'] = []


    data = {}

    # data same conformer
    data['natoms'] = []
    data['symbols'] = []
    data['mol_charge'] = 0

    # data eac conformer
    data['coordinates'] = []
    data['esp_values'] = []
    data['invr'] = []


    # number of atoms
    os.system("sed -n '1p' "+xyz_file+" > natoms")

    with open('natoms') as f:
        data['natoms'] = int(f.readline())

    # coordinates & symbols
    os.system("awk 'NR > 2 { print }' "+xyz_file+" > xyz_clean")

    xyz_clean = open('xyz_clean')
    lines_xyz = [line.split() for line in xyz_clean]

    coordinates = [[] for i in range(0,len(lines_xyz))]
    elements = [[] for i in range(0,len(lines_xyz))]

    for i in range(0,len(lines_xyz)):
        coordinates[i].append(lines_xyz[i][1])
        coordinates[i].append(lines_xyz[i][2])
        coordinates[i].append(lines_xyz[i][3])

        elements[i].append(lines_xyz[i][0])

    coordinates = [[float(j) for j in i] for i in coordinates]
    data['coordinates'] = np.array(coordinates)

    elements = [item for sublist in elements for item in sublist]
    data['symbols'] = [el.upper() for el in elements]

    # generate surface points 
    points = []

    for scale_factor in options['VDW_SCALE_FACTORS']:
        shell, radii = vdw_surface(data['coordinates'], data['symbols'], scale_factor)
        points.append(shell)

    points = np.concatenate(points)
    np.savetxt('grid.dat',points,fmt='%15.10f')

    # calculate esp values
    if os.path.exists('TAPE41') == True:
        os.system('rm TAPE41')

    os.system('cp '+os.environ['DOCK_LIB_DIR']+'/densf.sh .')
    os.system("sed '/Inline/ r grid.dat' densf.sh > run.sh")
    os.system("sh run.sh > densf_output")
    os.system('rm run.sh densf.sh')

    os.system(os.environ['AMSBIN']+'''/dmpkf TAPE41 SCF > grid_esp.dat''')

    os.system("awk 'NR > 3' grid_esp.dat > grid_esp_clean")

    grid_esp_clean = open('grid_esp_clean')
    grid_V = [line.split() for line in grid_esp_clean]
    grid_V = [item for sublist in grid_V for item in sublist]
    grid_V = [float(i) for i in grid_V]

    data['esp_values'].append(np.asarray(grid_V))

    # build a matrix of the inverse distance from each ESP point to each nucleus
    invr = np.zeros((len(points), len(coordinates)))
    for i in range(invr.shape[0]):
        for j in range(invr.shape[1]):
            invr[i, j] = 1/np.linalg.norm(points[i]-coordinates[j])

    data['invr'].append(invr)

    # calculate charges
    qf, labelf, notes = fit(options,data)

    with open("results_resp.out", "w") as f:
        f.write("Electrostatic potential parameters\n")
        if not options['GRID']:
            f.write("    van der Waals radii (Angstrom):\n")
            for i, j in radii.items():
                f.write("                               %8s%8.3f\n" %(i, j/scale_factor))
            f.write("    VDW scale factors:              ")
            for i in options["VDW_SCALE_FACTORS"]:
                f.write('%6.2f' %i)
            f.write('\n')
            f.write("    VDW point density:                %.3f\n" %(options["VDW_POINT_DENSITY"]))

        f.write("\nGrid information grid.dat:\n")
        f.write("    Number of grid points:            %d\n" %len(data['esp_values'][0]))

        f.write("\nRestraint\n")
        if options['RESTRAINT']:
            f.write("    Hyperbolic restraint to a charge of zero\n")
            if options['IHFREE']:
                    f.write("    Hydrogen atoms are not restrained\n")
            f.write("    resp_a:                           %.4f\n" %(options["RESP_A"]))
            f.write("    resp_b:                           %.4f\n" %(options["RESP_B"]))

            f.write("\nFit\n")
            f.write(notes)
            f.write("\nElectrostatic Potential Charges\n")
            f.write("   Center  Symbol")
            for i in labelf:
                f.write("%10s" %i)
            f.write("\n")
            for i in range(data['natoms']):
                f.write("  %5d    %s     " %(i+1, data['symbols'][i]))
                for j in qf:
                    f.write("%12.8f" %j[i])
                f.write("\n")
            f.write("Total Charge:    ")
            for i in qf:
                f.write("%12.8f"
                        %np.sum(i))
            f.write('\n')

        with open("charges","w") as output:
            for i in range(data['natoms']):
                output.write("%12.3f" %qf[1][i])
                output.write("\n")

    os.system("rm grid_esp_clean natoms xyz_clean")
    return
