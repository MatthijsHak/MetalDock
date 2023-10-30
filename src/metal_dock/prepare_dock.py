import random
import math
import os, sys, shutil, glob, subprocess

import itertools as it
import numpy as np
import networkx as nx
from networkx.algorithms.components.connected import connected_components

from collections import defaultdict 
from random import seed
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger

# Dictionary of all elements matched with their atomic masses.
mass_dict = {'H' : 1.008,'HE' : 4.003, 'LI' : 6.941, 'BE' : 9.012,\
            'B' : 10.811, 'C' : 12.011, 'N' : 14.007, 'O' : 15.999,\
            'F' : 18.998, 'NE' : 20.180, 'NA' : 22.990, 'MG' : 24.305,\
            'AL' : 26.982, 'SI' : 28.086, 'P' : 30.974, 'S' : 32.066,\
            'CL' : 35.453, 'AR' : 39.948, 'K' : 39.098, 'CA' : 40.078,\
            'SC' : 44.956, 'TI' : 47.867, 'V' : 50.942, 'CR' : 51.996,\
            'MN' : 54.938, 'FE' : 55.845, 'CO' : 58.933, 'NI' : 58.693,\
            'CU' : 63.546, 'ZN' : 65.38, 'GA' : 69.723, 'GE' : 72.631,\
            'AS' : 74.922, 'SE' : 78.971, 'BR' : 79.904, 'KR' : 84.798,\
            'RB' : 84.468, 'SR' : 87.62, 'Y' : 88.906, 'ZR' : 91.224,\
            'NB' : 92.906, 'MO' : 95.95, 'TC' : 98.907, 'RU' : 101.07,\
            'RH' : 102.906, 'PD' : 106.42, 'AG' : 107.868, 'CD' : 112.414,\
            'IN' : 114.818, 'SN' : 118.711, 'SB' : 121.760, 'TE' : 126.7,\
            'I' : 126.904, 'XE' : 131.294, 'CS' : 132.905, 'BA' : 137.328,\
            'LA' : 138.905, 'CE' : 140.116, 'PR' : 140.908, 'ND' : 144.243,\
            'PM' : 144.913, 'SM' : 150.36, 'EU' : 151.964, 'GD' : 157.25,\
            'TB' : 158.925, 'DY': 162.500, 'HO' : 164.930, 'ER' : 167.259,\
            'TM' : 168.934, 'YB' : 173.055, 'LU' : 174.967, 'HF' : 178.49,\
            'TA' : 180.948, 'W' : 183.84, 'RE' : 186.207, 'OS' : 190.23,\
            'IR' : 192.217, 'PT' : 195.085, 'AU' : 196.967, 'HG' : 200.592,\
            'TL' : 204.383, 'PB' : 207.2, 'BI' : 208.980, 'PO' : 208.982,\
            'AT' : 209.987, 'RN' : 222.081, 'FR' : 223.020, 'RA' : 226.025,\
            'AC' : 227.028, 'TH' : 232.038, 'PA' : 231.036, 'U' : 238.029,\
            'NP' : 237, 'PU' : 244, 'AM' : 243, 'CM' : 247, 'BK' : 247,\
            'CT' : 251, 'ES' : 252, 'FM' : 257, 'MD' : 258, 'NO' : 259,\
            'LR' : 262, 'RF' : 261, 'DB' : 262, 'SG' : 266, 'BH' : 264,\
            'HS' : 269, 'MT' : 268, 'DS' : 271, 'RG' : 272, 'CN' : 285,\
            'NH' : 284, 'FL' : 289, 'MC' : 288, 'LV' : 292, 'TS' : 294,\
            'OG' : 294}

def is_float(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

def center_of_mass(ligand):
    sum_mass = np.sum(ligand,axis=0)[3]

    centre = [sum([ligand[i][j]*ligand[i][3] for i in range(len(ligand))])  for j in range(3)]
    centre /= sum_mass
    return centre

def distance(a,b):
    d = np.sqrt( (a[0] - b[0])**2 + (a[1] - b[1])**2 + (a[2] - b[2])**2  )

    return d

def sphere_distance(r, a, b):
    dot = np.dot(a,b)
    cosine = dot / r**2 
    cosine_radian = cosine * (math.pi/180)
    d = r * np.arccos(cosine_radian)

    return d

def angle(a,b,c):
    a = np.array(a)
    b = np.array(b)
    c = np.array(c)
    ba = a - b
    bc = c - b

    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)
    angle = np.degrees(angle)

    return angle

def flatten(l):
    return [item for sublist in l for item in sublist]

def readXYZ(xyz_file, no_hydrogen=True):
    with open(xyz_file, 'r') as fin:
        xyz = [line.strip().split() for line in fin]
        xyz = xyz[2:]

        if no_hydrogen == True:
            coord = [[i[0], float(i[1]), float(i[2]), float(i[3])] for i in xyz if 'H' not in i]
        else:
            coord = [[float(i[1]), float(i[2]), float(i[3])] for i in xyz]
      
        return coord

def pdbqtToMol2(name_ligand):
    subprocess.call([os.environ['OBABEL']+f' -ipdbqt {name_ligand}.pdbqt -omol2 {name_ligand}.mol2 > {name_ligand}.mol2'],shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)
    return Chem.MolFromMol2File(f'{name_ligand}.mol2',sanitize=False)

def pdbqt_to_nx(mol2):
    G = nx.Graph()

    for atom in mol2.GetAtoms():
        G.add_node(atom.GetIdx(),
                   atomic_num=atom.GetAtomicNum(),
                   is_aromatic=atom.GetIsAromatic(),
                   atom_symbol=atom.GetSymbol())

    for bond in mol2.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(),
                   bond.GetEndAtomIdx(),
                   bond_type=bond.GetBondType())

    return G

def asCartesian(rthetaphi):
    #takes list rthetaphi (single coord) in radians
    r       = rthetaphi[0]
    theta   = rthetaphi[1]
    phi     = rthetaphi[2]
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return [x,y,z]

def asSpherical(xyz):
    #takes list xyz (single coord) to radians
    x       = xyz[0]
    y       = xyz[1]
    z       = xyz[2]
    r       =  np.sqrt(x*x + y*y + z*z)
    theta   =  np.arccos(z/r)
    phi     =  np.arctan2(y,x)
    return [r,theta,phi]

def return_idx_max_list(list):
    maxi = 0
    row = 0 
    for idx, x in enumerate(list):
        sum = 0
        for y in x:
            sum+= np.sum(y)    
        if sum > maxi:
            maxi = sum
            row = idx

    return row

def maximum_distance_sphere(n_ligands, atom_pos, metal_pos, radius=0.75):
    # set metal at origin
    atom_pos = [coord - metal_pos for coord in atom_pos]
    # get r and scale each ligand to unit sphere 
    new_atom_pos = [asSpherical(coord) for coord in atom_pos[1:]]
    r_pos = [new_atom_pos[coord][0] for coord in range(len(new_atom_pos))]
    scale_r = [radius / r_pos[coord] for coord in range(len(new_atom_pos))]

    for i in range(len(new_atom_pos)):
        new_atom_pos[i][0] = new_atom_pos[i][0] * scale_r[i]

    # return to Cartesian 
    new_atom_pos = [asCartesian(coord) for coord in new_atom_pos]

    pos_hydrogen = []
    dist_list = []
    new_xyz = []
    # generate equidistant points on sphere
    n_points = 10000

    alpha = (4 * np.pi * radius **2) / n_points
    d = np.sqrt(alpha)
    m_nu = int(np.round(np.pi/d))
    d_nu = np.pi/m_nu
    d_phi = alpha/d_nu
    
    for m in range (0,m_nu):
        nu = np.pi*(m+0.5)/m_nu
        m_phi = int(np.round(( 2 * np.pi *np.sin(nu))  / d_phi ))
        for n in range (0,m_phi):
            phi = 2*np.pi*n/m_phi
            x = radius*np.sin(nu)*np.cos(phi)
            y = radius*np.sin(nu)*np.sin(phi)
            z = radius*np.cos(nu)
            
            new_xyz = [x,y,z]

            sphere_dist = [sphere_distance(radius, atom, new_xyz)  for atom in new_atom_pos]
            euclid_dist = [distance(atom, new_xyz) for atom in new_atom_pos]

            dist = sphere_dist + euclid_dist
            dist_list.append(dist)

            pos_hydrogen.append([new_xyz[0],new_xyz[1],new_xyz[2]])

    # return row where distance is maximum
    row = return_idx_max_list(dist_list)
    # translate back to original coordinates
    new_pos_hydrogen = pos_hydrogen[row] + metal_pos
    return [new_pos_hydrogen[0], new_pos_hydrogen[1], new_pos_hydrogen[2]]


def delete_hydrogen(pdbqt_file):
    with open(pdbqt_file,'r') as fin:
        with open('output.pdbqt','w') as fout:
            for line in fin:
                if 'DD' in line:
                    pass
                else:
                    fout.write(line)
    shutil.move('output.pdbqt', pdbqt_file)
    return 
                
def create_ligand_pdbqt_file(par, name_ligand):
    with open(f'{name_ligand}.mol2','r') as fin_1:
        with open('CM5_charges','r') as fin_2:
            cm = [line.strip().split() for line in fin_2]
            cm = cm[1:]
            with open('output.mol2', 'w') as fout:
                atom_id = 0
                for line in fin_1:
                    if len(line.strip().split()) < 7:
                        fout.write(line)
                    else:
                        line = line.strip().split()
                        line[8] = cm[atom_id][1]
                        atom_id+=1
                        fout.write(f'     {line[0]:>2} {line[1]:<2}      {line[2]:>9} {line[3]:>9} {line[4]:>9} {line[5]:<5}   {line[6]:>1}  {line[7]:>4}       {line[8]:>6}\n')

    subprocess.call([os.environ['OBABEL']+' -imol2 output.mol2 -opdbqt '+name_ligand+'.pdbqt  > '+name_ligand+'.pdbqt'],shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    n_mdl = 0 
    with open(f'{name_ligand}.pdbqt','r') as fin:
        for line in fin:
            if 'MODEL' in line:
                n_mdl+=1
    
    if n_mdl == 1:
        one_model_file(par, name_ligand, f'{name_ligand}_c.xyz', f'{name_ligand}.pdbqt')
    else:
        multiple_model_file(par, name_ligand, f'{name_ligand}_c.xyz', f'{name_ligand}.pdbqt')
    return

def one_model_file(par, name_ligand, xyz_file, pdbqt_file):
    with open(pdbqt_file, 'r') as fin:
        new_lines = [line.strip().split() for line in fin]

    new_lines = [s for s in new_lines if "REMARK" not in s]

    for item in new_lines:
        try:
            if item[2] == par.metal_symbol.upper():
                metal_atom = item[1]
        except IndexError:
            pass

    new_pdbqt = []
    for line in new_lines:
        if 'MODEL' in line:
            pass
        elif 'ENDMDL' in line:
            pass
        else:
            new_pdbqt.append(line)

    if par.vacant_site == True:
        pos_hydrogen = find_pos_dummy(par, f'{name_ligand}_c.xyz', f'{name_ligand}.pdbqt',f'{name_ligand}')
        write_pdbqt(par, xyz_file, new_pdbqt, metal_atom, pos_hydrogen)
    else:
        write_pdbqt(par, xyz_file, new_pdbqt, metal_atom)
        
    return

def multiple_model_file(par, name_ligand, xyz_file, pdbqt_file):
    with open(pdbqt_file, 'r') as fin:
        new_lines = [line.strip().split() for line in fin]

    mdl = 1
    atom = 0
    mol_id = 0
    torsdof = 0

    new_pdbqt = []

    mol = globals()

    new_lines = [s for s in new_lines if "REMARK" not in s]
    
    for item in new_lines:
        try:
            if item[2] == par.metal_symbol.upper():
                metal_atom = item[1]
                metal_mdl = mdl 
        except IndexError:
            pass

        if 'MODEL' in item:
            mol[f'mol_{mdl}'] = []
            mol_id+=1

        if 'ENDMDL' not in item:
            mol[f'mol_{mdl}'].append(item)
        else:
            mol[f'mol_{mdl}'].append(item)
            mdl+=1

    # first insert the part that contains the metal
    for line in mol[f'mol_{metal_mdl}']:
        if 'TORSDOF' in line:
            torsdof+=int(line[1])
        elif 'ATOM' in line:
            atom+=1
            new_pdbqt.append(line)
        elif 'MODEL' in line:
            pass
        elif 'ENDMDL' in line:
            pass
        else:
            new_pdbqt.append(line)

    iter_over = list(range(mol_id))
    iter_over.remove(metal_mdl-1)

    for lig in iter_over:
        branch = False
        for line in mol[f'mol_{lig+1}']:
            if 'ATOM' in line and line[2] != 'H':
                branch_id = atom+int(line[1])
                break

        for line in mol[f'mol_{lig+1}']:
            if 'TORSDOF' in line:
                torsdof+=int(line[1])
            elif 'MODEL' in line:
                pass
            elif 'ROOT' in line:
                pass
            elif 'ATOM' in line:
                atom+=1
                if branch == False:
                    new_pdbqt.append(['BRANCH',f'{metal_atom}',f'{branch_id}'])
                    branch = True 
                new_pdbqt.append(line)
                
        new_pdbqt.append(['ENDBRANCH',f'{metal_atom}',f'{branch_id}'])
    new_pdbqt.append(['TORSDOF', str(torsdof)])

    if par.vacant_site == True:
        write_pdbqt(par, xyz_file, new_pdbqt, metal_atom)
        pos_hydrogen = find_pos_dummy(par, f'{name_ligand}_c.xyz', f'{name_ligand}.pdbqt',f'{name_ligand}')
        write_pdbqt(par, xyz_file, new_pdbqt, metal_atom, pos_hydrogen)
    else:
        write_pdbqt(par, xyz_file, new_pdbqt, metal_atom)

    return

def write_pdbqt(par, xyz_file, lines, metal_atom, pos_hydrogen=None):
    n_atoms = 1
    hydrogen_added = False
    with open(f'{par.name_ligand}.pdbqt', 'w') as fout:
        fout.write('ROOT\n')
        for line in lines:
            if pos_hydrogen != None and n_atoms == int(metal_atom)+1:
                fout.write(f'ATOM     {n_atoms:>2}  H   LIG A   1     {pos_hydrogen[0]:>7.3f} {pos_hydrogen[1]:>7.3f} {pos_hydrogen[2]:>7.3f}  0.00  0.00     0.000 DD\n')
                n_atoms+=1
                hydrogen_added = True

            if 'MODEL' in line:
                pass

            elif 'REMARK' in line:
                pass
            
            elif 'ENDROOT' in line:
                fout.write(f'ENDROOT\n')

            elif 'ATOM' in line:
                if len(line) < 13: 
                    fout.write(f'{line[0]}   {n_atoms:>4} {line[2]:>2}   LIG A   1     {line[5]:>7} {line[6]:>7} {line[7]:>7}  {line[8]:>4}  {line[9]:>4}    {line[10]:>6} {line[11]:<2}\n')
                    n_atoms+=1
                else:
                    fout.write(f'{line[0]}   {n_atoms:>4} {line[2]:>2}   LIG A   1     {line[6]:>7} {line[7]:>7} {line[8]:>7}  {line[9]:>4}  {line[10]:>4}    {line[11]:>6} {line[12]:<2}\n')
                    n_atoms+=1

            elif 'BRANCH' in line:
                if hydrogen_added == True:
                    if int(line[1]) > int(metal_atom):
                        line[1] = int(line[1])+1
                    if int(line[2]) > int(metal_atom):
                        line[2] = int(line[2])+1
                fout.write(f'{line[0]}   {line[1]}  {line[2]}\n')

            elif 'ENDBRANCH' in line:
                if hydrogen_added == True:
                    if int(line[1]) > int(metal_atom):
                        line[1] = int(line[1])+1
                    if int(line[2]) > int(metal_atom):
                        line[2] = int(line[2])+1
                fout.write(f'{line[0]}   {line[1]}  {line[2]}\n')

            elif 'TORSDOF' in line:
                fout.write(f'{line[0]} {line[1]}\n')
    return
 
def merge_common(lists): 
    neigh = defaultdict(set) 
    visited = set() 
    for each in lists: 
        for item in each: 
            neigh[item].update(each) 
    def comp(node, neigh = neigh, visited = visited, vis = visited.add): 
        nodes = set([node]) 
        next_node = nodes.pop 
        while nodes: 
            node = next_node() 
            vis(node) 
            nodes |= neigh[node] - visited 
            yield node 
    for node in neigh: 
        if node not in visited: 
            yield sorted(comp(node))

def find_pos_dummy(par, xyz_file, pdbqt_file, name_ligand):
    mol2 = pdbqtToMol2(name_ligand)
    G = pdbqt_to_nx(mol2)

    #get positions of the atoms from the pdbqt:
    coord = mol2.GetConformer(0).GetPositions()

    for k,v in G.nodes(data=True):
        if v['atom_symbol'].upper() == par.metal_symbol.upper():
            metal = k

    metal_pos = coord[metal]
    neighbor_list = sorted(list(G.neighbors(metal)))
    G.remove_node(metal)

    # Calculate atoms that belong to each ligand
    ligand_list = []
    dict_ligands = {}
    for n in neighbor_list:
        n_list = list(G.neighbors(n))
        n_list.append(n)
        ligand_list.append(sorted(n_list))

    ligand_list = list(merge_common(ligand_list))
    n_ligands = len(ligand_list)
    atom_pos = []
    atom_pos.append(metal_pos)
    G = pdbqt_to_nx(mol2)
    # Calculate centre of mass of ligands 
    for idx in range(n_ligands):
        mass = []
        for v in ligand_list[idx]:
            if v in neighbor_list: # ensure that the 
                xyz = list(coord[v])
                xyz.append(mass_dict[nx.get_node_attributes(G,'atom_symbol')[v].upper()])
                mass.append(xyz)

        atom_pos.append(center_of_mass((mass)))

    new_pos_hydrogen = maximum_distance_sphere(n_ligands, atom_pos, metal_pos)
    return new_pos_hydrogen

def get_coordinates(xyz_file, metal_symbol):
    dock_x = None
    
    with open(xyz_file, 'r') as fin:
        for _ in range(2):
            next(fin)
        for line in fin:
            if metal_symbol in line:
                 coordinates = line.strip().split()
                 dock_x = coordinates[1]
                 dock_y = coordinates[2]
                 dock_z = coordinates[3]
        if dock_x == None:
            raise ValueError('metal symbol not found in xyz file')
            
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

    return max_side

def prepare_receptor(name_protein):
    prepare_gpf4 = os.path.join(os.environ['MGLTOOLS'], 'prepare_receptor4.py')
    command = os.environ['PYTHON_3']+f' {prepare_gpf4} -U nphs -A None -r clean_{name_protein}.pdb'
    process = subprocess.Popen(
        command,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True
    )

    stdout, stderr = process.communicate()

    if process.returncode != 0:
        with open('prepare_receptor.out', 'w') as fout:
            fout.write(stdout)
            fout.write(stderr)
        print(f'ERROR DURING PREPARATION OF RECEPTOR, SEE /output/docking/prepare_receptor.out FOR DETAILS\n')
        print(f'IF PDB FILE IS NOT WRITTEN IN CORRECT PDB FORMAT, PLEASE EDIT MANUALLY\n')
        sys.exit()

def docking_func(par, parameter_set, name_ligand, name_protein, dock, box_size, energy=None):
    if os.path.exists(f'clean_{name_protein}.gpf'):
        os.remove(f'clean_{name_protein}.gpf')

    #create_gpf():
    prepare_gpf4 = os.path.join(os.environ['MGLTOOLS'], 'prepare_gpf4.py')
    subprocess.call([os.environ['PYTHON_3']+f" {prepare_gpf4} -l {name_ligand}.pdbqt  -r clean_{name_protein}.pdbqt -p parameter_file={par.parameter_file} -p npts='{box_size},{box_size},{box_size}' -p gridcenter='{dock[0]:.6},{dock[1]:.6},{dock[2]:.6}'"], shell=True)

    gpf = open(f'clean_{name_protein}.gpf', 'a')
    gpf.write(f'nbp_r_eps 0.25 23.2135   12 6  NA TZ\n')
    gpf.write(f'nbp_r_eps 2.10  3.8453   12 6  OA Zn\n')
    gpf.write(f'nbp_r_eps 2.25  7.5914   12 6  SA Zn\n')
    gpf.write(f'nbp_r_eps 1.00  0.0000   12 6  HD Zn\n')
    gpf.write(f'nbp_r_eps 2.00  0.0060   12 6  NA Zn\n')
    gpf.write(f'nbp_r_eps 2.00  0.2966   12 6  N  Zn\n')
    gpf.write(f'nbp_r_eps 2.20  {parameter_set[0]:>.4f}   12 10 NA {par.metal_symbol}\n')
    gpf.write(f'nbp_r_eps 2.25  {parameter_set[1]:>.4f}   12 10 OA {par.metal_symbol}\n')
    gpf.write(f'nbp_r_eps 2.30  {parameter_set[2]:>.4f}   12 10 SA {par.metal_symbol}\n')
    gpf.write(f'nbp_r_eps 1.00  {parameter_set[3]:>.4f}   12 6  HD {par.metal_symbol}\n')
    gpf.close()

    #autogrid()
    autogrid4 = os.path.join(os.environ['ROOT_DIR'],'external','AutoDock','autogrid4')
    if par.method.lower() == 'mc':
        subprocess.call([f'{autogrid4} -p clean_{name_protein}.gpf'], shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    else:
        subprocess.call([f'{autogrid4} -p clean_{name_protein}.gpf'], shell=True)

    #create_dpf()
    write_dpf_file(f'clean_{name_protein}.gpf', name_ligand, f'clean_{name_protein}', par.parameter_file, par.num_poses, par.dock_algorithm, random_pos=par.random_pos, SA=par.sa_dock, GA=par.ga_dock, energy_ligand=energy)

    #autodock()
    autodock4 = os.path.join(os.environ['ROOT_DIR'],'external','AutoDock','autodock4')
    if par.method.lower() == 'train' or par.method.lower() == 'mc':
        subprocess.call([f'{autodock4} -p {name_ligand}_clean_{name_protein}.dpf'], shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    else:
        subprocess.call([f'{autodock4} -p {name_ligand}_clean_{name_protein}.dpf'], shell=True)

    #write_all_conformations()
    write_conformations(name_ligand, name_protein)
    # write_conformations = os.path.join(os.environ['MGLTOOLS'], 'write_conformations_from_dlg.py')
    # subprocess.call([os.environ['PYTHON_3']+f" {write_conformations} -d {name_ligand}_clean_{name_protein}.dlg"], shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    return

def write_conformations(name_ligand, name_protein):
    with open(f'{name_ligand}_clean_{name_protein}.dlg', 'r') as fin:
        pose = []
        atom_id = 0
        mol_id = 1
        in_docked_block = False  
        docked_block = []  
        for line in fin:
            if 'DOCKED: ROOT' in line:
                in_docked_block = True
                docked_block.append(line)
            elif in_docked_block and 'TER' in line:
                in_docked_block = False
                docked_block.append(line)
                with open(f'{name_ligand}_{mol_id}.pdbqt', 'w') as output_file:
                    for block_line in docked_block:
                        cleaned_line = block_line.replace('DOCKED: ', '', 1)
                        output_file.write(cleaned_line)
                mol_id += 1
                docked_block = []  
            elif in_docked_block:
                docked_block.append(line)


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

    if GA == True and SA == False:
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
    if GA == False and SA == True:
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

    dpf_file.write('analysis                             # perforem a ranked cluster analysis\n')
    return

def rmsd_func(name_ligand, n_prot, directory, generation=None, num_gen=None, train=False, standard=False, test=False):
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
    while os.path.exists(os.path.join(directory,'docking',f'{name_ligand}_{i}.pdbqt')):
        pdbqt_file = os.path.join(output_dir,'docking',f'{name_ligand}_{i}.pdbqt')
        d.delete_hydrogen(pdbqt_file)
        subprocess.call([os.environ['OBABEL']+" -ipdbqt "+name_ligand+"_{}.pdbqt".format(i)+" -oxyz "+name_ligand+"_{}.xyz".format(i)+" -d > "+name_ligand+"_{}.xyz".format(i)], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        rmsd_func = os.path.join(os.environ['ROOT_DIR'],'metal_dock','calculate_rmsd.py')
        rmsd_non_rotate = float(subprocess.getoutput([os.environ['PYTHON_3']+f' {rmsd_func} ref.xyz '+name_ligand+'_{}.xyz'.format(i)+' -nh --reorder --rotation none --translation none']))
        rmsd = rmsd_non_rotate

        rmsd_list.append(rmsd)

        rmsd_list.append("RMSD for Conformation %i = %.4f"% (i, rmsd))
        i += 1

    for j in range(0,len(rmsd_print_list)):
        output.append(rmsd_print_list[j])

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
