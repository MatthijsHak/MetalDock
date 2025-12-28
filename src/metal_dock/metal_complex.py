import os

import subprocess as sp
import numpy as np
import networkx as nx

from scipy.optimize import differential_evolution

from src.metal_dock.adf_engine import ADFEngine
from src.metal_dock.gaussian_engine import GaussianEngine
from src.metal_dock.orca_engine import ORCAEngine
from src.metal_dock.xyz2graph import MolGraph, to_networkx_graph
from src.metal_dock.logger import MetalDockLogger


class MetalComplex:
    def __init__(self, par):
        self.par = par
        self.logger = MetalDockLogger()
        self.qm_engine = ADFEngine(par, self) if par.engine.lower() == 'adf' else \
                        GaussianEngine(par, self) if par.engine.lower() == 'gaussian' else \
                        ORCAEngine(par, self) if par.engine.lower() == 'orca' else None

    def canonicalize_ligand(self):
        """
        This function canonicalizes the ligand.
        """
        self.logger.info("CANONICALIZING LIGAND STRUCTURE\n")
        output_xyz_file = self.par.output_dir / 'file_prep'/ f'{self.par.name_ligand}_c.xyz'
        sp.call([os.environ['OBABEL']+f' -ixyz {self.par.xyz_file} -oxyz {output_xyz_file} --canonical > {output_xyz_file}'],shell=True, stdout=sp.PIPE, stderr=sp.PIPE)

    def create_mol_graph(self, 
                         run_type: str):
        """
        This function creates the molecular graph from the xyz file.

        Args:
            run_type (str): The type of run.
        """
        mg = MolGraph()
        if self.par.score_only:
            mg.read_xyz(self.par.input_dir / f'{self.par.name_ligand}.xyz')
        else:
            mg.read_xyz(self.par.output_dir / 'QM' / run_type / 'output.xyz')
        self.graph = to_networkx_graph(mg)

    def add_charges_to_graph(self):
        """
        This function adds the charges to the graph.
        """
        for node, data in self.graph.nodes(data=True):
            element = data.get("element")
            if element:
                node_key = f'{element.upper()}{node+1}'
                if node_key in self.charges:
                    self.graph.nodes[node]['charge'] = self.charges[node_key]

    def create_ligand_pdbqt_file(self):
        """
        This function creates the pdbqt file for the ligand.
        """
        # Start at the metal atom
        for node in self.graph.nodes():
            if self.graph.nodes[node]['element'] == self.par.metal_symbol:
                self.metal_atom = node
                break

        ligand_graphs = self._obtain_ligand_subgraphs()

        self.rotatable_bonds = self._find_rotatable_bonds(ligand_graphs)

        #! I need to check fi the branches always exclude hydrogens 
        branches, branch_connections = self._create_ligand_branches(ligand_graphs)

        # Calculate atom index mapping and distances to hydrogen
        self._initialize_atom_index_mapping()
        self._convert_elements_to_autodock_elements()
        self._add_old_to_pdbqt_mapping(branches)    
        self.write_pdbqt(branches, branch_connections, N_ligands=len(ligand_graphs))

    def write_pdbqt(self, 
                    branches: dict, 
                    branch_connections: dict, 
                    N_ligands: int):
        """
        This function writes the pdbqt file.

        Args:
            branches (dict): A dictionary containing the branches.
            branch_connections (dict): A dictionary containing the branch connections.
            N_ligands (int): The number of ligands.
            pos_hydrogen (list): A list containing the position of the hydrogen atom.
        """

        # create the pdbqt file
        pdbqt_path = self.par.output_dir / 'file_prep' / f'{self.par.name_ligand}.pdbqt'
        with open(pdbqt_path, 'w') as fout:
            branch_atoms = branches['ROOT']
            # First write the root
            fout.write(f'ROOT\n')
            for atom in branch_atoms:
                element = self.graph.nodes[atom]['element']
                if element == 'DD':
                    element = 'H'

                autodock_element = self.atom_index_mapping[atom]['autodock_element']
                xyz = self.graph.nodes[atom]['xyz']
                charge = self.graph.nodes[atom]['charge']
                pdbqt_index = self.atom_index_mapping[atom]['pdbqt_index'] + 1
                fout.write(f'ATOM   {pdbqt_index:>4} {element:<2}    LIG A   1    {xyz[0]:>7.3f} {xyz[1]:>7.3f} {xyz[2]:>7.3f}  0.00  0.00    {float(charge):>6.3f} {autodock_element:<2}\n')

            fout.write('ENDROOT\n')
            # Then extract all branches for each ligand
            for ligand_idx in range(N_ligands):
                # Obtain all branches that include 'ligand_{ligand_idx}' in their keys
                ligand_branches = {key: value for key, value in branches.items() if f'ligand_{ligand_idx}' in key}

                # Sort the branches based on the branch index
                ligand_branches = sorted(ligand_branches.items(), key=lambda x: int(x[0].split('_')[1]))
                endbranch_lines = []  # Collect ENDBRANCH lines

                for branch_name, branch_atoms in ligand_branches:
                        # Process each branch as needed
                        connection_info = next((v for k, v in branch_connections.items() if v['to_name'] == branch_name), None)
                        if connection_info:
                            from_atom = connection_info['from_atom']
                            from_atom = self.atom_index_mapping[from_atom]['pdbqt_index'] + 1
            
                            to_atom = connection_info['to_atom'] 
                            to_atom = self.atom_index_mapping[to_atom]['pdbqt_index'] + 1

                            fout.write(f'BRANCH {from_atom} {to_atom}\n')
                        
                        for atom in branch_atoms:
                            element = self.graph.nodes[atom]['element']
                            autodock_element = self.atom_index_mapping[atom]['autodock_element']
                            xyz = self.graph.nodes[atom]['xyz']
                            charge = self.graph.nodes[atom]['charge']
                            pdbqt_index = self.atom_index_mapping[atom]['pdbqt_index'] + 1 

                            fout.write(f'ATOM     {pdbqt_index:>2} {element:<2}    LIG A   1    {xyz[0]:>7.3f} {xyz[1]:>7.3f} {xyz[2]:>7.3f}  0.00  0.00    {float(charge):>6.3f} {autodock_element:<2}\n')
                        
                        if connection_info:
                            endbranch_lines.append(f'ENDBRANCH {from_atom} {to_atom}\n')

                # Write ENDBRANCH lines in reverse order
                for line in reversed(endbranch_lines):
                    fout.write(line)

    def _initialize_atom_index_mapping(self):
        """
        This function initializes the atom index mapping.
        """
        
        self.atom_index_mapping = {}
        # add extra position for the dummy atom at node position -1 
        if self.par.vacant_site == True:
            self.atom_index_mapping[-1] = {'pdbqt_index': None, 'autodock_element': None}
            self._add_dummy_atom_to_graph()
            
        for node in self.graph.nodes():
            self.atom_index_mapping[node] = {'pdbqt_index': None, 'autodock_element': None}

    def _convert_elements_to_autodock_elements(self):
        """
        This function converts the element to the autodock element.
        """
        # loop over the graph and convert the element to the autodock element
        for node in self.graph.nodes():
            # identify the element
            element = self.graph.nodes[node]['element']
            if element == 'H':
                # verify if it is bound to a O, N, or S
                for neighbor in self.graph.neighbors(node):
                    if self.graph.nodes[neighbor]['element'] in ['O', 'N', 'S']:
                        autodock_element = 'HD'
                        break
                else:
                    autodock_element = 'H'

            elif element == 'C':
                # Check if part of a ring
                if any(node in cycle for cycle in nx.cycle_basis(self.graph)):
                    # Obtain the bond order of the bonds to this element
                    for neighbor in self.graph.neighbors(node):
                        if self.graph.nodes[neighbor]['element'] == 'H':
                            pass
                        # Check the bond order for all elements
                        if self.graph.edges[node, neighbor]['bond_order'] > 1.10:
                            autodock_element = 'A'  # Aromatic or special element
                        else:
                            autodock_element = element.upper()  # Default to uppercase of the element
                else:
                    autodock_element = 'C'

            elif element == 'O':
                autodock_element = 'OA'

            elif element == 'N':
                # verif if it has three neigbors 
                if len(list(self.graph.neighbors(node))) == 3:
                    autodock_element = 'N'
                else:
                    autodock_element = 'NA'

            else:
                autodock_element = element

            self.atom_index_mapping[node]['autodock_element'] = autodock_element

    def _add_old_to_pdbqt_mapping(self, 
                                  branches: dict):
        """
        This function adds the old index to the pdbqt index mapping.

        Args:
            branches (dict): A dictionary containing the branches.
        """
        pdbqt_index = 0  # Start new index from 1 for PDBQT file

        if self.par.vacant_site == True:
            # add the dummy atom to the root after the metal atom
            branches['ROOT'].insert(branches['ROOT'].index(self.metal_atom) + 1, -1)

        # Iterate over the ROOT branch first
        for atom in branches['ROOT']:
            self.atom_index_mapping[atom]['pdbqt_index'] = pdbqt_index
            pdbqt_index += 1

        # obtain all branches except the ROOT
        ligand_branches = {key: value for key, value in branches.items() if key != 'ROOT'}

        # sort the branches based on the ligand index
        ligand_branches = sorted(ligand_branches.items(), key=lambda x: int(x[0].split('_')[1]))

        # iterate over all branches and print the atom index mapping
        for branch_name, branch_atoms in ligand_branches:
            for atom in branch_atoms:
                self.atom_index_mapping[atom]['pdbqt_index'] = pdbqt_index
                pdbqt_index += 1

    def _add_dummy_atom_to_graph(self):
        """
        This function adds a dummy atom to the graph.
        It calculates the center of mass for ligand atoms in the same ring interacting with the metal atom.
        """
        # Determine all atoms bound to the metal atom
        metal_atom_neighbors = list(self.graph.neighbors(self.metal_atom))

        # Initialize a list to store positions of ligand atoms or their COM
        ligand_positions = []

        # Get all cycles (rings) in the graph
        cycles = nx.cycle_basis(self.graph)

        # Iterate over each neighbor
        for neighbor in metal_atom_neighbors:
            # Check if the neighbor is part of any ring
            neighbor_in_cycle = False
            for cycle in cycles:
                if neighbor in cycle:
                    # Find all other neighbors of the metal atom in the same cycle
                    ring_neighbors = [n for n in metal_atom_neighbors if n in cycle]
                    if len(ring_neighbors) > 1:
                        # Calculate the center of mass for these neighbors
                        positions = [self.graph.nodes[n]['xyz'] for n in ring_neighbors]
                        center_of_mass = np.mean(positions, axis=0)
                        ligand_positions.append(center_of_mass)
                        neighbor_in_cycle = True
                        break
                    
            if not neighbor_in_cycle:
                # If the neighbor is not part of any ring, add its position
                ligand_positions.append(self.graph.nodes[neighbor]['xyz'])

        hydrogen_position = self._find_max_distance_position(ligand_positions, self.graph.nodes[self.metal_atom]['xyz'])

        # add to the mapping 
        self.graph.add_node(-1, element='DD', xyz=hydrogen_position, charge=0)
        self.atom_index_mapping[-1] = {'pdbqt_index': None, 'autodock_element': 'DD'}

    @staticmethod
    def spherical_to_cartesian(radius: float, 
                               theta: float, 
                               phi: float, 
                               center: list):
        """
        Convert spherical to cartesian coordinates.

        Args:
            radius (float): The radius of the sphere.
            theta (float): The angle in radians.
            phi (float): The angle in radians.
            center (list): The center of the sphere.

        Returns:
            list: The cartesian coordinates.
        """
        sin_theta = np.sin(theta)
        x = radius * sin_theta * np.cos(phi) + center[0]
        y = radius * sin_theta * np.sin(phi) + center[1]
        z = radius * np.cos(theta) + center[2]
        return np.array([x, y, z])

    def _find_max_distance_position(self, 
                                   ligand_positions: list, 
                                   metal_pos: list):
        """
        Find the position on a unit sphere that maximizes distance from ligand positions.

        Args:
            ligand_positions (list): A list of ligand positions.
            metal_pos (list): The position of the metal atom.

        Returns:
            list: The position of the maximum distance position.
        """
        # Convert inputs to numpy arrays for efficient calculations
        metal_pos = np.asarray(metal_pos)
        ligand_positions = np.asarray(ligand_positions) - metal_pos  # Shift to origin

        def objective_function(params):
            """Negative sum of distances to ligand positions (to maximize)."""
            theta, phi = params
            test_point = self.spherical_to_cartesian(1, theta, phi, [0, 0, 0])  # Centered at origin
            distances = np.linalg.norm(test_point - ligand_positions, axis=1)
            return -np.sum(distances)  # Negative for maximization

        # Global optimization using differential evolution
        result = differential_evolution(
            objective_function, 
            bounds=[(0, np.pi), (0, 2 * np.pi)],  
            maxiter=5000,  # Reduce iterations
            tol=1e-8
        )

        # Convert optimal spherical coordinates back to cartesian
        theta_opt, phi_opt = result.x
        max_distance_position = self.spherical_to_cartesian(1, theta_opt, phi_opt, metal_pos)
        
        return max_distance_position

    def _create_ligand_branches(self, 
                                ligand_graphs: list):
        """
        This function creates the branches for the ligand graphs.

        Args:
            ligand_graphs (list): A list of ligand graphs.
            metal_atom (int): The index of the metal atom in the graph.

        Returns:
            tuple: A tuple containing the branches and branch connections.
        """
        # Initialize branches with the metal atom as the root
        visited = set()
        branches = {'ROOT': [self.metal_atom]}
        branch_connections = {}  # Dictionary to store branch connections
        ligand_idx = 0
        connection_idx = 0  # To keep track of connection indices

        # Iterate over each ligand graph
        for ligand_graph in ligand_graphs:
            # Start at the metal atom for each ligand
            current_branch = 'ROOT'
            current_path = [self.metal_atom]

            def dfs(atom, current_branch, visited):
                nonlocal connection_idx
                visited.add(atom)  # Mark the atom as visited
                for neighbor in ligand_graph.neighbors(atom):
                    if neighbor not in visited:
                        current_path.append(neighbor)
                        if self._is_rotatable_bond(atom, neighbor) and self._is_not_in_ring(atom, neighbor) and self._is_dihedral_bond(atom, neighbor):
                            # Create a new branch if the bond is rotatable
                            new_branch_name = f"Branch_{len(branches)}_ligand_{ligand_idx}"
                            branches[new_branch_name] = [neighbor]
                            # Store the connection information
                            connection_name = f'connection_{connection_idx}'
                            branch_connections[connection_name] = {
                                'from_name': current_branch,
                                'to_name': new_branch_name,
                                'from_atom': atom,
                                'to_atom': neighbor
                            }
                            connection_idx += 1
                            dfs(neighbor, new_branch_name, visited)
                        else:
                            # Add to the current branch if the bond is not rotatable
                            branches[current_branch].append(neighbor)
                            dfs(neighbor, current_branch, visited)
                        current_path.pop()

            # Start DFS from the metal atom
            dfs(self.metal_atom, current_branch, set())
            ligand_idx += 1

        return branches, branch_connections  # Return both branches and connections 

    def _is_dihedral_bond(self, 
                          atom1: int, 
                          atom2: int):
        """
        This function checks if the bond between two atoms is the middle bond of a dihedral angle.
        Hydrogens are not considered in the neighbors.
        """
        # Get neighbors excluding hydrogens
        atom_1_neighbors = [n for n in self.graph.neighbors(atom1) if self.graph.nodes[n]['element'] != 'H']
        atom_2_neighbors = [n for n in self.graph.neighbors(atom2) if self.graph.nodes[n]['element'] != 'H']

        # Check if both atoms have at least one non-hydrogen neighbor
        if len(atom_1_neighbors) >= 2 and len(atom_2_neighbors) >= 2:
            return True
        else:
            return False

    def _is_not_in_ring(self, 
                        atom1: int, 
                        atom2: int):
        """
        This function checks if the bond between two atoms is not in the same ring.
        Returns True if the atoms are not in the same ring, False otherwise.
        """
        for cycle in nx.cycle_basis(self.graph):
            # Check if both atoms are in the same cycle
            if atom1 in cycle and atom2 in cycle:
                return False
        return True

    def _is_rotatable_bond(self, 
                            atom1: int, 
                            atom2: int):
        """
        This function checks if the bond between two atoms is rotatable.

        Args:
            atom1 (int): The index of the first atom.
            atom2 (int): The index of the second atom.

        Returns:
            bool: True if the bond is rotatable, False otherwise.
        """
        return (atom1, atom2) in self.rotatable_bonds or (atom2, atom1) in self.rotatable_bonds

    def _find_rotatable_bonds(self, 
                              ligand_graphs: list):
        """
        This function finds all the rotatable bonds in the ligand graphs.

        Args:
            ligand_graphs (list): A list of ligand graphs.

        Returns:
            set: A set of rotatable bonds.
        """
        rotatable_bonds = set()
        for ligand_graph in ligand_graphs:
            rotatable_bonds.update(self._generate_proper_dihedrals(ligand_graph))

        # set should be order independent
        rotatable_bonds = {(min(a, b), max(a, b)) for a, b in rotatable_bonds}
        return rotatable_bonds

    def _generate_proper_dihedrals(self, 
                                  ligand_graph: nx.Graph):
        """
        This function generates all the proper dihedrals in the ligand graph.

        Args:
            ligand_graph (nx.Graph): The ligand graph.

        Returns:
            set: A set of proper dihedrals.
        """
        proper_dihedrals = set()

        for atom in ligand_graph.nodes():
            dihedrals = self._find_dihedrals(ligand_graph, atom, 3)

            for dihedral in dihedrals:
                # Check bond order between the 2nd and 3rd elements
                atom2, atom3 = dihedral[1], dihedral[2]
                if self._is_valid_bond_order(atom2, atom3) and not (atom2, atom3) in proper_dihedrals:
                    proper_dihedrals.add((atom2, atom3))

        return proper_dihedrals

    def _find_dihedrals(self, 
                        graph: nx.Graph, 
                        node: int, 
                        depth: int):
        """
        Generates the dihedrals for a given node in a graph.

        Args:
            graph (nx.Graph): The ligand graph.
            node (int): The index of the node.
            depth (int): The depth of the dihedral.

        Returns:
            list: A list of dihedrals.
        """
        if depth == 0:
            return[[node]]
        
        dihedral = []
        for neighbor in graph.neighbors(node):
            for path in self._find_dihedrals(graph,neighbor,depth-1):
                if node not in path:
                    dihedral.append([node]+path)
        return dihedral

    def _is_valid_bond_order(self, 
                             atom1: int, 
                             atom2: int):
        """
        This function checks if the bond order between two atoms is valid. 
        This is done by checking if the bond order is within the range of 0.8  to 1.20.
        Args:
            atom1 (int): The index of the first atom.
            atom2 (int): The index of the second atom.

        Returns:
            bool: True if the bond order is valid, False otherwise.
        """
        # Retrieve the bond order for the edge between atom1 and atom2
        edge_data = self.graph.get_edge_data(atom1, atom2)
        if edge_data is None:
            return False  # or handle it as needed

        bond_order = edge_data.get('bond_order', 0)
        # Check if the bond order is within the specified range
        return 0.8  < bond_order < 1.2

    def _obtain_ligand_subgraphs(self):
        """
        This function obtains the ligand subgraphs from the metal complex graph.

        Args:
            metal_atom (int): The index of the metal atom in the graph.

        Returns:
            list: A list of ligand subgraphs.
        """
        # Get all the atoms that are bonded to the metal atom
        metal_bonded_atoms = list(self.graph.neighbors(self.metal_atom))

        # To keep track of visited atoms
        visited = set()

        # Obtain subgraphs for each ligand
        ligand_subgraphs = []
        for atom in metal_bonded_atoms:
            if atom not in visited:
                # Perform a BFS/DFS to get all connected atoms of this ligand
                subgraph_atoms = self._get_ligand_subgraph(atom, visited)
                # add the metal atom to the subgraph
                subgraph_atoms.add(self.metal_atom)
                ligand_subgraphs.append(self.graph.subgraph(subgraph_atoms).copy())

        return ligand_subgraphs

    def _get_ligand_subgraph(self, 
                             start_atom: int, 
                             visited: set):
        """
        Breadth First Search to obtain the ligand subgraph.

        Args:
            start_atom (int): The index of the starting atom in the graph.
            visited (set): A set of visited atoms.
            metal_atom (int): The index of the metal atom in the graph.

        Returns:
            set: A set of atoms that are part of the ligand subgraph.
        """
        # Use a queue for BFS
        queue = [start_atom]
        subgraph_atoms = set()

        while queue:
            current_atom = queue.pop(0)
            if current_atom not in visited:
                visited.add(current_atom)
                subgraph_atoms.add(current_atom)
                # Add neighbors to the queue if they are not the metal atom
                for neighbor in self.graph.neighbors(current_atom):
                    # Ensure we do not traverse back to the metal atom
                    if neighbor != self.metal_atom and neighbor not in visited:
                        queue.append(neighbor)

        return subgraph_atoms
