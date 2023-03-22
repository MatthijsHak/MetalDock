import os, sys
import numpy as np
import networkx as nx 
import networkx.algorithms.isomorphism as iso

from rdkit import Chem
from rdkit.Chem import AllChem
import rdkit.Chem.rdMolTransforms
import rdkit.Geometry

from scipy.optimize import minimize

from typing import Tuple

import mdtraj as md


def xyz_to_nx(name_ligand: str) -> Tuple[nx.Graph, AllChem.rdchem.Mol]:
    """ Creates Graph Object from Mol File

    Arguments 
    ---------
    name_ligand : str

    Output
    ---------
    Graph : nx.Graph
    Mol : rdkit.Chem.rdchem.Mol
    """

    # mol_1 = next(py.readfile('xyz',f'{xyz_file}'))
    # mol_1.write('mol', f'{name_ligand}.mol',overwrite=True)

    Mol = Chem.MolFromMolFile(f'{name_ligand}.mol',removeHs=False, sanitize=False)
    G = nx.Graph()

    for idx, atom in enumerate(Mol.GetAtoms()):
        coord = Mol.GetConformer().GetAtomPosition(idx)
        G.add_node(atom.GetIdx(),
                atomic_num=atom.GetAtomicNum(),
                is_aromatic=atom.GetIsAromatic(),
                atom_symbol=atom.GetSymbol(),
                x_coord=coord.x,
                y_coord=coord.y,
                z_coord=coord.z,
                atom_pos=idx+1)

    for bond in Mol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(),
                bond.GetEndAtomIdx(),
                bond_type=bond.GetBondType())

    return G, Mol


def squared_residuals(new_pos_hydrogen, atom_coords, target_distances):
    distances = np.linalg.norm(new_pos_hydrogen - atom_coords, axis=1)
    return np.sum((distances - target_distances)**2)


if __name__ == '__main__':

    init_graph, init_mol = xyz_to_nx('init')
    pose_graph, pose_mol = xyz_to_nx('pose')

    GM = nx.algorithms.isomorphism.GraphMatcher(init_graph,pose_graph)
    
    for subgraph in GM.subgraph_isomorphisms_iter():
        sub_dict = subgraph

    pose_atom_order = [[nx.get_node_attributes(init_graph,'atom_pos')[key], nx.get_node_attributes(init_graph,'atom_pos')[value]] for key,value in sub_dict.items()]

    pose_coords = []
    # Initalize a list for the new coordinates
    for idx, pose_atom in enumerate(pose_mol.GetAtoms()):      
       pose_coords.append(pose_mol.GetConformer().GetAtomPosition(idx))
    
    pose_coords = np.array(pose_coords)

    # Initialize a list for the new positions
    init_hydrogens = []
    pose_hydrogens = []
    placed_hydrogens = []
    is_hydrogen_placed = False
    for i, init_atom in enumerate(init_mol.GetAtoms()):
        if init_atom.GetSymbol() == "H":
            # Get the coordinates of the first atom
            hydrogen_atom = np.array(init_mol.GetConformer().GetAtomPosition(i))
            init_hydrogens.append(hydrogen_atom)
            # Initialize an array to store the target distances
            target_distances = []

            # Loop over all atoms
            for j in pose_atom_order:
                # Get the position of the current atom
                atom_pos = np.array(init_mol.GetConformer().GetAtomPosition(j[0]-1))
                # Subtract the position of the first atom to get the difference vector
                dist_init = np.linalg.norm(atom_pos - hydrogen_atom)
                # Store the difference vector
                target_distances.append(dist_init)
            
            # # If a hydrogen is placed consider them as well 
            # if is_hydrogen_placed == True:
            #     atom_pos = np.array(init_mol.GetConformer().GetAtomPosition(i-1))
            #     placed_hydrogens.append(atom_pos)
            #     for pos in placed_hydrogens:
            #         dist_init = np.linalg.norm(pos - hydrogen_atom)
            #         target_distances.append(dist_init)

            #     pose_coords = np.vstack((pose_coords, atom_pos))

            target_distances = np.array(target_distances)
            optimized_pos = minimize(squared_residuals, hydrogen_atom, args=(pose_coords,target_distances), method='BFGS',tol=1e-10)
            pose_hydrogens.append(optimized_pos.x)
            # is_hydrogen_placed = True

    # Write hydrogens to new xyz file
    with open('pose.xyz','r') as fin:
        fin_lines = fin.readlines()
        fin_lines[0] = f"{init_mol.GetNumAtoms()}\n"
        with open('test_output.xyz','w+') as fout:
            fout.writelines(fin_lines)
            for coord in pose_hydrogens:
                fout.write(f"H        {coord[0]:>9.5f}      {coord[1]:>9.5f}      {coord[2]:>9.5f}\n")