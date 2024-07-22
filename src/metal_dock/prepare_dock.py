import random
import math
import os, sys, shutil, glob, subprocess

import itertools as it
import numpy as np
import networkx as nx
from networkx.algorithms.components.connected import connected_components

from .xyz2graph import MolGraph, to_networkx_graph
from openbabel import openbabel as ob

from scipy.spatial.distance import cdist

from collections import defaultdict 
from random import seed
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger

ob.obErrorLog.SetOutputLevel(0)

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

def extract_interacting_residues(pose, protein, cutoff=4.0):
    pose_xyz = []

    with open(pose, 'r') as fin:
        for _ in range(2):
            next(fin)
        for line in fin:
            split = line.strip().split()[1:4]
            pose_xyz.append([float(split[0]), float(split[1]), float(split[2])])

    interacting_residues = []
    interacting_residues_xyz = []
    # calculate all the residues that are within 4 angstroms of the ligand
    with open(protein,'r') as fin:
        for line in fin:
            if 'ATOM' in line or 'HETATM' in line:
                # calculate the distance between the pose_xyz and the 
                split = line.strip().split()
                interacting_residues.append((split[3], split[5]))
                interacting_residues_xyz.append([float(split[6]), float(split[7]), float(split[8])])

    distance_matrix = cdist(pose_xyz, interacting_residues_xyz,'euclidean')

    # Find indices of residues within 4 angstroms
    within_4_angstrom_indices = np.any(distance_matrix <= 4.0, axis=0)

    # Get the unique residues that are within 4 angstroms
    unique_residues_within_4_angstrom = set()
    for i, within in enumerate(within_4_angstrom_indices):
        if within:
            unique_residues_within_4_angstrom.add(interacting_residues[i])

    # Convert set to a list for further processing or output
    unique_residues_within_4_angstrom = list(unique_residues_within_4_angstrom)

    # order based on the id of the residues 
    unique_residues_within_4_angstrom = sorted(unique_residues_within_4_angstrom, key=lambda x: int(x[1]))

    return unique_residues_within_4_angstrom

def extract_dlg(dlg_file, var):
    # print the scores of the docking
    idx = 0
    scores = []
    with open(dlg_file,'r') as fin:
        for line in fin:
            if 'RMSD TABLE' in line or idx > 0:
                scores.append(line)
                idx += 1
                if idx == 19:
                    break

    binding_energy = []
    for score in scores:
        # extract the binding energy 
        split = score.strip().split()
        if len(split) == 7:
            binding_energy.append((int(split[2]), float(split[3])))

    # sort the binding eneryg based on the first element
    binding_energy = sorted(binding_energy, key=lambda x: x[0])

    # calculate the ligand efficiency
    ligand_efficiency = []
    # ligand efficiecny = binding energy / number of heavy atoms
    for be in binding_energy:
        ligand_efficiency.append(be[1] / var.n_heavy_atoms)

    return binding_energy, ligand_efficiency

def write_pose_to_pdb(xyz_file_name, pdb_file_name):
    graph = MolGraph()
    graph.read_xyz(xyz_file_name)
    
    G = to_networkx_graph(graph)

    nodes_sorted = sorted(G.nodes(data=True), key=lambda x: x[0])

    with open(pdb_file_name, 'w') as f:
        f.write('MODEL         1\n')
        for node in nodes_sorted:
            f.write(f"HETATM{node[0]+1:>5} {node[1]['element'].upper():>2}   UNL          {node[1]['xyz'][0]:>8.3f}{node[1]['xyz'][1]:>8.3f}{node[1]['xyz'][2]:>8.3f}  1.00  0.00          {node[1]['element']:>2}\n")
        for edge in G.edges():
            f.write(f"CONECT {edge[0]+1:>4} {edge[1]+1:>4}\n")
        f.write('ENDMDL\n')

# extract from the following bonds the implicit H, set subsequently the number of implicit Hs correct
def add_non_polar_hydrogens(name_xyz_H, name_pdb_H, name_xyz_noH, output_file_name):
    #=== Count the number of hydrogen atoms per heavy atom  ===#
    atoms = []
    bonds = []
    atom_constraints = [] 
    with open(name_pdb_H,'r') as fin:
        for line in fin:
            if 'HETATM' in line or 'ATOM' in line:
                split = line.strip().split()
                atoms.append((int(split[1]), split[-1]))

            elif 'CONECT' in line:
                split = line.strip().split()
                # add all combinations to the bonds list in form of tuples 
                atom1 = split[1]
                for atom2 in split[2:]:
                    bonds.append((int(atom1), int(atom2))) 

    # count the number of hydrogen atoms bound to each atom carbon atom 
    implicit_H = [0 for i in range(len(atoms))]
    for bond in bonds:
        atom1 = bond[0]
        atom2 = bond[1]

        # make sure it is a carbon atom 
        if 'C' == atoms[atom1 - 1][1] and 'H' == atoms[atom2 - 1][1]: 
            implicit_H[atom1 - 1] += 1

    # go through the xyz without the hydrogens and ensure that all indices are in the constraints list
    with open(name_xyz_noH,'r') as fin:
        for _ in range(2):
            next(fin)
        for idx, line in enumerate(fin):
            atom_constraints.append(idx)

    #=== Perform mapping between the docked and initial structure ===#
    # Step 1: Obtain graphs for struct_init and struct_dock
    mg_init = MolGraph()
    mg_init.read_xyz(name_xyz_H)
    G_init = to_networkx_graph(mg_init)

    mg_dock = MolGraph()
    mg_dock.read_xyz(name_xyz_noH)
    G_dock = to_networkx_graph(mg_dock)

    G_init_noH = G_init.copy()
    G_dock_noH = G_dock.copy()

    for node in list(G_init_noH.nodes):
        if G_init_noH.nodes[node]['element'] == 'H':
            G_init_noH.remove_node(node)

    for node in list(G_dock_noH.nodes):
        if G_dock_noH.nodes[node]['element'] == 'H':
            G_dock_noH.remove_node(node)

    gm = nx.algorithms.isomorphism.GraphMatcher(G_init_noH, G_dock_noH)
    gm.is_isomorphic()
    mapping = gm.mapping   

    # Reorder implicit_H using mapping
    implicit_H_mapped = [None] * len(implicit_H)  
    for current_index, new_index in mapping.items():
        implicit_H_mapped[new_index] = implicit_H[current_index]

    #=== Add and optimize non-polar hydrogen atoms ===#
    conv = ob.OBConversion()
    conv.SetInAndOutFormats('xyz','xyz')
    mol = ob.OBMol()
    conv.ReadFile(mol, name_xyz_noH)

    for i in list(range(1,len(atoms)+1)):
        # add hydrogens
        atom = mol.GetAtom(i)

        if implicit_H_mapped[i-1] == None:
            atom.SetImplicitHCount(0)
        else:
            atom.SetImplicitHCount(implicit_H_mapped[i-1])
        mol.AddHydrogens(atom)

    # Define constraints
    constraints = ob.OBFFConstraints()
    for i in atom_constraints:
        constraints.AddAtomConstraint(i)

    # Setup the force field with the constraints
    forcefield = ob.OBForceField.FindForceField("gaff")
    forcefield.Setup(mol)
    forcefield.SetConstraints(constraints)

    # Do a 500 steps conjugate gradient minimiazation
    # and save the coordinates to mol.
    forcefield.ConjugateGradients(1000, 1E-8)
    forcefield.GetCoordinates(mol)

    # Write the mol to a file
    conv.WriteFile(mol, output_file_name)

    # return the indices of the atoms that have not been added 
    return atom_constraints

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

def center_molecule(input_dir, xyz_file, metal_symbol):
    xyz = []
    atom_symbols = []

    with open(xyz_file, 'r') as f:
        for _ in range(2):
            next(f)
        for line in f:
            if not line.strip():  
                break 
            coords = line.split()[1:]
            xyz.append([float(coord) for coord in coords])
            atom_symbols.append(line.split()[0])

    # find the ruthenium atom and center the molecule around it
    ruthenium_index = atom_symbols.index(metal_symbol)
    ruthenium_coords = xyz[ruthenium_index]

    centered_xyz = np.array(xyz) - ruthenium_coords

    out_file = xyz_file.split('/')[-1][:-4]

    # write the centered coordinates to a new file
    with open(f'{out_file}_centered.xyz', 'w') as f:
        f.write(f'{len(centered_xyz)}\n')
        f.write('centered molecule\n')
        for atom, coords in zip(atom_symbols, centered_xyz):
            f.write(f'{atom} {coords[0]:10.7f} {coords[1]:10.7f} {coords[2]:10.7f}\n')

    return os.path.join(os.getcwd(), f'{out_file}_centered.xyz'  )

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


def delete_dummy_atom(pdbqt_file):
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
            if par.engine.lower() != 'gaussian':
                cm = cm[1:]
                print(cm)

            with open('output.mol2', 'w') as fout:
                atom_id = 0
                for line in fin_1:
                    if len(line.strip().split()) < 7 or 'Properties' in line:
                        fout.write(line)
                    else:
                        line = line.strip().split()
                        line[8] = cm[atom_id][1]
                        atom_id+=1
                        fout.write(f'     {line[0]:>2} {line[1]:<2}      {line[2]:>9} {line[3]:>9} {line[4]:>9} {line[5]:<5}   {line[6]:>1}  {line[7]:>4}       {line[8]:>6}\n')

    subprocess.call([os.environ['OBABEL']+' -imol2 output.mol2 -opdbqt '+name_ligand+'.pdbqt  > '+name_ligand+'.pdbqt'],shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    subprocess.call([os.environ['OBABEL']+' -imol2 output.mol2 -oxyz '+name_ligand+'_c.xyz  > '+name_ligand+'_c.xyz'],shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    n_mdl = 0 
    with open(f'{name_ligand}.pdbqt','r') as fin:
        for line in fin:
            if 'MODEL' in line:
                n_mdl+=1

    if n_mdl == 1:
        one_model_file(par, name_ligand, f'{name_ligand}_c.xyz', f'{name_ligand}.pdbqt')
    else:
        multiple_model_file(par, name_ligand, f'{name_ligand}_c.xyz', f'{name_ligand}.pdbqt')

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
        write_pdbqt(par, new_pdbqt, pos_hydrogen)
    else:
        write_pdbqt(par, new_pdbqt)
        
    return

def multiple_model_file(par, name_ligand, xyz_file, pdbqt_file):
    # rename and reorder the pdbqt file 
    reindex_atoms_branches(pdbqt_file, 'out.pdbqt')
    # I need to ensure that the atom numbering is the same for the graph as the pdbqt file 
    with open('out.pdbqt', 'r') as fin:
        new_lines = [line.strip().split() for line in fin]

    n_atoms = 0
    xyz_lines = []
    for line in new_lines:
        if 'ATOM' in line:
            n_atoms+=1
            xyz_lines.append([line[2], float(line[6]), float(line[7]), float(line[8])])

    # add to the beginning of the list and empty list and the number of atoms
    xyz_lines.insert(0, [n_atoms])

    with open('out.xyz', 'w') as f:
        for line in xyz_lines:
            if len(line) == 1:
                f.write(f'{line[0]}\n\n')
            else:
                # set the second letter to small letter if length of letter is 2
                if len(line[0]) == 2:
                    f.write(f'{line[0][0]}{line[0][1].lower()} {line[1]:.6f} {line[2]:.6f} {line[3]:.6f}\n')
                else:
                    f.write(f'{line[0]} {line[1]:.6f} {line[2]:.6f} {line[3]:.6f}\n')

    pdbqt_file = 'out.pdbqt'
    atom = 0

    new_lines = [s for s in new_lines if "REMARK" not in s]

    graph = MolGraph()
    graph.read_xyz('out.xyz')
    
    G = to_networkx_graph(graph)

    # Remove all nodes that have the element 'H' except for the 'H' atoms
    # that are connected to O, N, S, P
    remove_hydrogens(G)

    # find which atom is the metal atom from the graph
    for node in G.nodes(data=True):
        if node[1]['element'] == par.metal_symbol:
            metal_idx= node[0] # add 1 to the index to match the pdbqt file

    # generate the ligand indicies and obtain the hapticity interactions 
    hapt_atoms = []
    neighbor_list = sorted(list(G.neighbors(metal_idx)))

    ligands_list = []
    for node in neighbor_list:
        ligand = bfs(G, node, metal_idx)
        if ligand not in ligands_list:
            ligands_list.append(ligand)

    for n in neighbor_list:
        path_from_n = find_groups(G, n, neighbor_list, set())
        if sorted(path_from_n) not in hapt_atoms:
            hapt_atoms.append(sorted(path_from_n))

    hapticity_atoms = []
    for idx,ligand in enumerate(ligands_list):
        for connect in hapt_atoms: 
            if any(item in ligand for item in connect):
                hapticity_atoms.append(connect) 

    # extract all branches 
    branches, connections = extract_branches(new_lines)

    root = []
    branches_to_root = []
    root.append(['ROOT'])

    # identify the model in which the metal atom is located
    branches_metal, _ = identify_model(branches, connections, metal_idx+1, G)

    metal_branch = None
    for n in branches_metal.keys():
        if str(metal_idx) in branches_metal[n]['atoms']:
            for line in branches_metal[n]['lines']:
                if 'ATOM' in line:
                    root.append(line)
            metal_branch = n
    # remove the branch from the branches list
    branches_metal.pop(metal_branch)

    # The ligand_list and hapticity atoms that are already in the root can be removed 
    for line in root[1:]:
        for idx, (ligand, hapt_atoms) in enumerate(zip(ligands_list, hapticity_atoms)):
            if int(line[1]) in ligand or int(line[1]) in hapt_atoms:
                ligands_list.pop(idx)
                hapticity_atoms.pop(idx)

    # The ligand atoms need to be placed in the root branch
    for idx, (ligand, hapt_atoms) in enumerate(zip(ligands_list, hapticity_atoms)):
        ligand_graph = G.subgraph(ligand)
        # No hapticity ligand with small number of atoms
        # Add all atoms to the root branch
        if len(ligand) <= 3 and len(hapt_atoms) == 1:
            for atom in ligand:
                root.append(search_pdbqt_file(pdbqt_file, G.nodes[atom]['xyz']))

        # No hapticity ligand with larger number of atoms
        elif len(ligand) > 3 and len(hapt_atoms) == 1:
            print(ligand, hapt_atoms)
            # The ligand atoms need to be placed in a new branch
            # Add branch between metal atom and the interacting atom 
            # Search the ligand for more branches 
            branches, connections = extract_branches(new_lines)
            branches, connections = identify_model(branches, connections, hapt_atoms, G)
            branches_ligand, connections_ligand = remove_root_branches(branches, root, connections)

            # # all the atoms in the ligand are already added
            if branches_ligand == None and connections_ligand == None:
                continue

            branches_to_root = add_branches_to_root(branches_to_root, branches_ligand, connections_ligand, metal_idx,  G, root)

        # Hapticity ligand where all atoms are connected to the metal atom
        elif ligand == hapt_atoms:
            # Add to the root branch 
            for atom in hapt_atoms:
                root.append(search_pdbqt_file(pdbqt_file, G.nodes[atom]['xyz']))

        # Hapticity ligand with small number of atoms
        elif len(ligand) <= 3 and len(hapt_atoms) > 1:
            # The ligand atoms need to be placed in the root branch
            for atom in ligand:
                root.append(search_pdbqt_file(pdbqt_file, G.nodes[atom]['xyz']))

        # Hapticity ligand with larger number of atoms
        elif len(ligand) > 3 and len(hapt_atoms) > 1:
            # Add hapticity atoms to the root branch if they are not already added
            for atom in hapt_atoms:
                # get xyz of the atom
                xyz = G.nodes[atom]['xyz']
                # Check if the atom is already in the root branch
                atom_exists = False
                for line in root:
                    if 'ATOM' not in line:
                        continue
                    # Compare the coordinates with a tolerance to account for floating-point precision issues
                    if line[6] == f'{xyz[0]:.3f}' and line[7] == f'{xyz[1]:.3f}' and line[8] == f'{xyz[2]:.3f}':
                        atom_exists = True
                        break
                
                # If the atom is not in the root branch, add it
                if not atom_exists:
                    root.append(search_pdbqt_file(pdbqt_file, xyz))
   
            # extract all branches 
            branches, connections = extract_branches(new_lines)
            branches, connections = identify_model(branches, connections, hapt_atoms, G)
            branches_ligand, connections_ligand = remove_root_branches(branches, root, connections)

            # # all the atoms in the ligand are already added
            if branches_ligand == None and connections_ligand == None:
                continue
            
            branches_to_root = add_branches_to_root(branches_to_root, branches_ligand, connections_ligand, metal_idx,  G, root)

    root.append(['ENDROOT'])


    # append the branches to the root
    for branch in branches_to_root:
        root.append(branch)

    # root.append(final_branch[0])

    # find the metal_atom
    for item in root:
        try:
            if item[2] == par.metal_symbol.upper():
                metal_atom = item[1]
        except IndexError:
            pass

    # remove out.pdbqt and out.xyz
    os.remove('out.pdbqt')
    os.remove('out.xyz')

    if par.vacant_site == True:
        write_pdbqt(par, root)
        pos_hydrogen = find_pos_dummy(par, f'{name_ligand}_c.xyz', f'{name_ligand}.pdbqt',f'{name_ligand}')
        write_pdbqt(par, root, pos_hydrogen)
    else:
        write_pdbqt(par, root)

def add_branches_to_root(branches_to_root, branches_ligand, connections_ligand, metal_idx, G, root):
    final_branch = []
    branches_added = []
    end_branches = []
    # find the branch that is connected to the root branch 
    root_atoms = [int(line[1])-1 for line in root if 'ATOM' in line]
    # I need to iterate not over the hapt atoms, but over the atoms in the root 
    for root_atom in root_atoms:
        n = list(G.neighbors(root_atom)) 
        for neighbor in n:
            # if we find a branch start new branch
            end_branches = []
            # I need to know that the neighbor is not in the root branch 
            if neighbor not in root_atoms and neighbor != metal_idx and root_atom != metal_idx:
                root_branch_atom_xyz = G.nodes[neighbor]['xyz']
                branches_to_root.append(['BRANCH', root_atom+1, neighbor+1])
                final_branch.append(['ENDBRANCH', root_atom+1, neighbor+1])

                main_branch = None
                # verify which branch the xyz is located and add the branch to the root
                for branch in branches_ligand.keys():
                    if list(root_branch_atom_xyz) in branches_ligand[branch]['xyz']:
                        # add the branch to the root 
                        for line in branches_ligand[branch]['lines']:
                            if 'BRANCH' in line:
                                pass # atom_id_new_branch = [line[1], line[2]]
                            else:
                                branches_to_root.append(line)

                        main_branch = branch

                connected_branches = []
                for connection in connections_ligand:
                    if (connection[0] in branches_ligand[main_branch]['atoms'] 
                        or connection[1] in branches_ligand[main_branch]['atoms']):
                        connected_branches.append(connection)
                
                if connected_branches == []:
                    for endbranch in end_branches:
                        branches_to_root.append(endbranch)
                    continue
                else:
                    branches_to_root, end_branches = recursive_branches_write(
                        main_branch, connections_ligand, branches_ligand, 
                        end_branches, branches_added, branches_to_root
                    )

                    if end_branches:
                        for endbranch in end_branches[::-1]:
                            branches_to_root.append(endbranch)

    branches_to_root.append(final_branch[0])
    return branches_to_root

def recursive_branches_write(current_branch, connections, branches_ligand, end_branches, branches_added, root):
    current_branch_atoms = branches_ligand[current_branch]['atoms']
    next_branch = None
    current_connection = None

    # I need to keep track of the added branches up until a point where the branch ends
    # if next_branch is None:  # I need to verify if any of the branches in the branches_added has a connection to other connections if not then print the branch to the root
    for branch, line in branches_added[::-1]:
        if any(atom in connection for atom in branches_ligand[branch]['atoms'] for connection in connections):
            break
        else:
            root.append(line)
            branches_added.remove((branch, line))
            end_branches.remove(line)

    
    # Use the atom in the connection that is in the current branch to obtain the atom in the next branch
    for connection in connections:
        if connection[0] in current_branch_atoms:
            next_branch_atom = connection[1]
            current_connection = connection
            connections.remove(connection)
            break
        elif connection[1] in current_branch_atoms:
            next_branch_atom = connection[0]
            current_connection = connection
            connections.remove(connection)
            break

    # Find the next branch 
    # if current_connection:
    for branch in branches_ligand.keys():
        if next_branch_atom in branches_ligand[branch]['atoms']:
            next_branch = branch
            break

    # Write the next branch to the root
    # if current_connection and next_branch:
    root.append(['BRANCH', current_connection[0], current_connection[1]]) 
    end_branches.append(['ENDBRANCH', current_connection[0], current_connection[1]])
    branches_added.append((next_branch, ['ENDBRANCH', current_connection[0], current_connection[1]]))
    for line in branches_ligand[next_branch]['lines']:
        if 'ATOM' in line:
            root.append(line)
    

    next_branch_atoms = branches_ligand[next_branch]['atoms']
    connection_to_this_branch = []
    for connection in connections:
        if connection[0] in next_branch_atoms or connection[1] in next_branch_atoms:
            connection_to_this_branch.append(connection)


    if connection_to_this_branch:
        for connection in connection_to_this_branch:
            root, end_branches = recursive_branches_write(next_branch, connections, branches_ligand, end_branches, branches_added, root)

    # remove all tuples from list
    for branch in end_branches:
        if type(branch) == tuple:
            end_branches.remove(branch)
    return root, end_branches
                            
def identify_model(branches, connections, atoms, G):
    # go over the xyz coordinates and find the model that contains the hapticity atoms 
    # if hapticity atoms is list then iterate over the list
    if hasattr(atoms, '__iter__'):
        for atom in atoms:
            xyz_hapt = G.nodes[atom]['xyz']
            for model in branches.keys():
                for branch in branches[model].keys():
                    for line in branches[model][branch]['xyz']:
                        if xyz_hapt == tuple(line):
                            return branches[model], connections[model]
    else:
        xyz_hapt = G.nodes[atoms]['xyz']
        for model in branches.keys():
            for branch in branches[model].keys():
                for line in branches[model][branch]['xyz']:
                    if xyz_hapt == tuple(line):
                        return branches[model], connections[model]

def reindex_atoms_branches(input_file, output_file):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    atom_index_offset = 0
    branch_index_offset = 0

    model = {}
    mdl = 0

    # Iterate over the different text

    # For first model no changes 

    # For second model add number of atoms to atom index and each branch index

    # For third model add number of atoms to atom index and each branch index

    with open(output_file, 'w') as f:
        for line in lines:
            if 'MODEL' in line:
                mdl += 1
                model[mdl] = {'atom_index_offset': 0}

            if line.startswith("ATOM") and mdl == 1:
                f.write(line)
                model[mdl]['atom_index_offset'] += 1

            elif line.startswith("ATOM") and mdl > 1:
                atom_index_offset = model[mdl-1]['atom_index_offset']
                atom_index = int(line.split()[1])
                new_atom_index = atom_index + atom_index_offset
                f.write(f"ATOM {new_atom_index} {' '.join(line.split()[2:])}\n")
                model[mdl]['atom_index_offset'] += 1

            elif line.startswith("BRANCH") and mdl == 1:
                f.write(line)

            elif line.startswith("BRANCH") and mdl > 1:
                atom_index_offset = model[mdl-1]['atom_index_offset']
                branch_indices = [int(i) for i in line.split()[1:]]
                new_branch_indices = [str(i + atom_index_offset) for i in branch_indices]
                f.write(f"BRANCH {' '.join(new_branch_indices)}\n")

            elif line.startswith("ENDBRANCH") and mdl == 1:
                f.write(line)

            elif line.startswith("ENDBRANCH") and mdl > 1:
                atom_index_offset = model[mdl-1]['atom_index_offset']
                branch_indices = [int(i) for i in line.split()[1:]]
                new_branch_indices = [str(i + atom_index_offset) for i in branch_indices]
                f.write(f"ENDBRANCH {' '.join(new_branch_indices)}\n")

            elif 'REMARK' in line:
                pass
            else:
                f.write(line)

def remove_root_branches(branches, root, connections_ligand):
    delete_branches = set()
    # remove the branches that are already in the root based on the xyz coordinates
    for branch in branches.keys():
        for line in branches[branch]['lines']:
            for root_line in root:
                if 'ATOM' in line and 'ATOM' in root_line:
                # find if the xyz of the branch is in the root
                    if [float(root_line[6]), float(root_line[7]), float(root_line[8])] in branches[branch]['xyz']:
                        # Check if key exists before deleting
                        delete_branches.add(branch)
                        break

    for branch in delete_branches:
        del branches[branch]

    # also remove the connections that are already in the root based on the atom indices
    atom_indices_root = []
    for line in root:
        if 'ATOM' in line:
            atom_indices_root.append(line[1])

    delete_connections = []
    for atom in atom_indices_root:
        for connections in connections_ligand:
            if atom in connections:
                delete_connections.append(connections)

    for connection in delete_connections:
        connections_ligand.remove(connection)

    return branches, connections_ligand

def extract_branches(ligand_pdbqt_lines):
    model_idx = 0
    branch_idx = 0
    branches = {}
    connections = {}

    current_root = False
    current_branch = False

    for line in ligand_pdbqt_lines:
        if 'MODEL' in line:
            model_idx += 1
            branches[model_idx] = {}
            connections[model_idx] = []

        if 'ROOT' in line:
            current_root = True
            branches[model_idx][branch_idx] = {'lines': [], 'atoms': [], 'xyz': []}
            branches[model_idx][branch_idx]['lines'].append(line)
            
        elif 'ENDROOT' in line:
            branches[model_idx][branch_idx]['lines'].append(line)
            current_root = False
            branch_idx += 1

        elif current_root:
            branches[model_idx][branch_idx]['lines'].append(line)
            branches[model_idx][branch_idx]['atoms'].append(line[1])
            branches[model_idx][branch_idx]['xyz'].append([float(line[6]), float(line[7]), float(line[8])])

        elif 'BRANCH' in line: #and not current_branch:
            current_branch = True
            branch_idx += 1
            branches[model_idx][branch_idx] = {'lines': [], 'atoms': [], 'xyz': []}
            branches[model_idx][branch_idx]['lines'].append(line)
            connections[model_idx].append([line[1], line[2]])
            # branches[branch_idx]['atoms'].append(line[1])
            # branches[branch_idx]['xyz'].append([float(line[5]), float(line[6]), float(line[7])])

        elif 'ENDBRANCH' in line:
            branches[model_idx][branch_idx]['lines'].append(line)
            current_branch = False

        elif current_branch:
            branches[model_idx][branch_idx]['lines'].append(line)
            branches[model_idx][branch_idx]['atoms'].append(line[1])
            branches[model_idx][branch_idx]['xyz'].append([float(line[6]), float(line[7]), float(line[8])])

    return branches, connections

def remove_hydrogens(graph):
    nodes_to_remove = []
    
    # Iterate over all nodes with their attributes
    for node, data in graph.nodes(data=True):
        if data.get('element') == 'H':
            keep_hydrogen = False
            # Check if this hydrogen is connected to any O, N, S, P
            for neighbor in graph.neighbors(node):
                if graph.nodes[neighbor].get('element') in ['O', 'N', 'S', 'P']:
                    keep_hydrogen = True
                    break
            if not keep_hydrogen:
                nodes_to_remove.append(node)
    
    # Remove the nodes outside the iteration to avoid modification during iteration
    graph.remove_nodes_from(nodes_to_remove)

def search_pdbqt_file(pdbqt_file, xyz):
    '''
    Search the pdbqt file and return the line of the atom where the xyz coordinates
    match the input coordinates
    '''
    with open(pdbqt_file, 'r') as fin:
        for line in fin:
            if 'ATOM' in line or 'HETATM' in line:
                split = line.strip().split()
                if compare_decimals([float(split[6]), float(split[7]), float(split[8])], xyz):
                    return split 

def compare_decimals(set1, set2):
    rounded_set1 = tuple(round(num, 3) for num in set1)
    rounded_set2 = tuple(round(num, 3) for num in set2)
    return rounded_set1 == rounded_set2

def bfs(G, starting_node, metal_idx):
    ligand_atoms = []
    queue = [starting_node]

    while queue:
        atom = queue.pop(0)

        if atom not in ligand_atoms:
            ligand_atoms.append(atom)

        for i in G.neighbors(atom):
            if i not in ligand_atoms and i != metal_idx:
                queue.append(i)
                
    return sorted(ligand_atoms)

def find_groups(graph, current_node, metal_neighbors, visited):
    idx_group = []  # Create a new group for each starting node

    def recursive_search(node):
        if node in metal_neighbors and node not in visited and node not in idx_group:
            visited.add(node)
            idx_group.append(node)
            for neighbor in graph.neighbors(node):
                recursive_search(neighbor)

    recursive_search(current_node)
    return idx_group

def write_pdbqt(par, lines, pos_hydrogen=None):
    # I need to keep track of the old index and the new index such that I can shift the branch ids 
    mapping = {}

    if pos_hydrogen != None:
        n_atoms = 2
    else:
        n_atoms = 1

    # create mapping
    for line in lines:
        if 'ATOM' in line:
            mapping[int(line[1])] = n_atoms
            n_atoms+=1

    n_atoms = 1

    torsdof = 0
    with open(f'{par.name_ligand}.pdbqt', 'w') as fout:
        fout.write('ROOT\n')
        for line in lines:
            if pos_hydrogen != None and n_atoms == 1:
                fout.write(f'ATOM     {n_atoms:>2}  H   LIG A   1     {pos_hydrogen[0]:>7.3f} {pos_hydrogen[1]:>7.3f} {pos_hydrogen[2]:>7.3f}  0.00  0.00     0.000 DD\n')
                n_atoms+=1

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
                elif len(line[2]) == 3:
                    fout.write(f'{line[0]}   {n_atoms:>4} {line[2]:>3}   LIG A   1    {line[6]:>7} {line[7]:>7} {line[8]:>7}  {line[9]:>4}  {line[10]:>4}    {line[11]:>6} {line[12]:<2}\n')
                    n_atoms+=1
                elif len(line[2]) == 4:
                    fout.write(f'{line[0]}   {n_atoms:>4} {line[2]:>4}   LIG A   1    {line[6]:>7} {line[7]:>7} {line[8]:>7}  {line[9]:>4}  {line[10]:>4}    {line[11]:>6} {line[12]:<2}\n')
                    n_atoms+=1
                else:
                    fout.write(f'{line[0]}   {n_atoms:>4} {line[2]:>2}   LIG A   1     {line[6]:>7} {line[7]:>7} {line[8]:>7}  {line[9]:>4}  {line[10]:>4}    {line[11]:>6} {line[12]:<2}\n')
                    n_atoms+=1

            elif 'BRANCH' in line:
                fout.write(f'{line[0]}   {mapping[int(line[1])]}  {mapping[int(line[2])]}\n')
                torsdof+=1

            elif 'ENDBRANCH' in line:
                fout.write(f'{line[0]}   {mapping[int(line[1])]}  {mapping[int(line[2])]}\n')

        
        fout.write(f'TORSDOF {torsdof}\n')
 
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

    return [x_npts, y_npts, z_npts]

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

def docking_func(par, name_ligand, name_protein, dock, box_size, energy=None):
    if os.path.exists(f'clean_{name_protein}.gpf'):
        os.remove(f'clean_{name_protein}.gpf')

    #create_gpf():
    prepare_gpf4 = os.path.join(os.environ['MGLTOOLS'], 'prepare_gpf4.py')
    subprocess.call([os.environ['PYTHON_3']+f" {prepare_gpf4} -l {name_ligand}.pdbqt  -r clean_{name_protein}.pdbqt -p parameter_file={par.parameter_file} -p npts='{box_size[0]},{box_size[1]},{box_size[2]}' -p gridcenter='{dock[0]:.6},{dock[1]:.6},{dock[2]:.6}'"], shell=True)

    gpf = open(f'clean_{name_protein}.gpf', 'a')
    gpf.write(f'nbp_r_eps 0.25 23.2135   12 6  NA TZ\n')
    gpf.write(f'nbp_r_eps 2.10  3.8453   12 6  OA Zn\n')
    gpf.write(f'nbp_r_eps 2.25  7.5914   12 6  SA Zn\n')
    gpf.write(f'nbp_r_eps 1.00  0.0000   12 6  HD Zn\n')
    gpf.write(f'nbp_r_eps 2.00  0.0060   12 6  NA Zn\n')
    gpf.write(f'nbp_r_eps 2.00  0.2966   12 6  N  Zn\n')
    if par.internal_param == False:
        gpf.write(f'nbp_r_eps 2.20  {par.parameter_set[0]:>.4f}   12 10 NA {par.metal_symbol}\n')
        gpf.write(f'nbp_r_eps 2.25  {par.parameter_set[1]:>.4f}   12 10 OA {par.metal_symbol}\n')
        gpf.write(f'nbp_r_eps 2.30  {par.parameter_set[2]:>.4f}   12 10 SA {par.metal_symbol}\n')
        gpf.write(f'nbp_r_eps 1.00  {par.parameter_set[3]:>.4f}   12 6  HD {par.metal_symbol}\n')
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

    # dpf_file.write('# Unbound Ligand Parameters\n')
    # if energy_ligand != None:
        # dpf_file.write('unbound_energy '+str(energy_ligand)+'              # set the energy of the unbound state\n')
    
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
