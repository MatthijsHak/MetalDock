#!/usr/bin/python

# import packages
from __future__ import division
from __future__ import with_statement
from __future__ import absolute_import
import sys
import os
import math as mt
import bisect
from functools import partial
import itertools
from operator import itemgetter
from itertools import groupby
from io import open

# test input arguments and save dictionary
def optparser(argv=sys.argv[1:]):
	
	''' Returns a dictionary "args" with bias parameter file name, gpf file name,
		dpf file name and/or dpf directory name
	'''
	
	# usage: input line
	usage = '''
	Usage:

	python prepare_bias.py -b bias_parameter_file [-g gpf] [-d dpf] [-D dpfs_dir] [-m map_file]

		-b : input file for bias (x y z Energy Radius Type)
		-g : AutoDock4 grid parameter file
		-d : AutoDock4 docking parameter file
		-D : modify all dpf files (and ligand pdbqts if necessary) in
		     dpfs_dir directory
		-m : specific input map file
				
	Requires -b and at least one of the other arguments (-g, -d, -D or -m)
	'''
	
	# initiate dictionary with arguments
	args = dict()
	
	# add bias parameter file name to dictionary
	if '-b' in argv:
		bias_file_index = argv.index('-b') + 1
		args.update({'bias_file' : argv[bias_file_index]})
	else:
		print 'Error: -b argument required indicating bias parameter file'
		print usage
		sys.exit(2)		
	
	# add gpf file name to dictionary
	if '-g' in argv:
		gpf_index = argv.index('-g') + 1
		args.update({'gpf_file' : argv[gpf_index]})
	
	# add dpf file name to dictionary
	if '-d' in argv:
		dpf_index = argv.index('-d') + 1
		args.update({'dpf_file' : argv[dpf_index]})
	
	# add dpfs directory name to dictionary
	if '-D' in argv:
		dir_index = argv.index('-D') + 1
		args.update({'dpfs_directoy' : argv[dir_index]})

	# add map file name to dictionary
	if '-m' in argv:
		map_index = argv.index('-m') + 1
		args.update({'map_file' : argv[map_index]})
		
	# test number of arguments
	if len(argv) < 4:
		print 'Error: invalid number of arguments\n'
		print 'Requires -b and at least one of the other arguments (-g, -d, -D or -m)'
		print usage
		sys.exit(2)
	
	# output dictionary with bias_file_name, gpf_file_name, dpf_file_name and/or dpf_directory_name
	return args


# requires input bias parameter file in current directory
# get coordinates, energy, radius and type of bias sites
class Bias_Sites(object):
	
	''' Returns a list of dictionaries "BS" with coordinates, energy,
		radius and type for each bias site
	'''
	
	def __init__(self, bias_file):

		self.biassites = self.parse_file(bias_file)
	
	# parse bias parameter file and get coordinates, energy, radius and type of bias sites
	def parse_file(self, bias_file):

		# open bias parameter file
		with open(bias_file) as f:
			
			# initiate BS list
			BS = []

			# readlines puts each line of the file as an element of the list
			for i, line in enumerate(f.readlines()):
				
				# line has a list with each space or tab seaprated string of the line as element
				line = line.split()
				
				# each line needs 6 columns
				if ((len(line) != 6)):
					print 'Error: invalid number of columns in bias parameter file'
					sys.exit(2)
				# only operate in numeric lines (skip title line)
				elif (((line[0].split('.')[0]).strip("-").isnumeric())):
					# create a dictionary (bs) with bias site data
					bs = {	'x':float(line[0]),
							'y':float(line[1]),
							'z':float(line[2]),
							'radius':float(line[4]),
							'int_type':str(line[5].lstrip('\n')) }
					# energy bias definition (negative number)
					if float(line[3]) < 0:
						bs.update( {'DG':float(line[3]) } )
					# PFP bias definition (positive number)
					else:
						bs.update( {'PFP':float(line[3])} )
					
					# BS is then a list of dictionaries, one dictionary per bias site
					# with coordinates, energy, radius and type of the site
					BS.append(bs)
		
		# return a list of dictionaries, one dictionary per bias site
		# with coordinates, energy, radius and type of the site
		return(BS)


# requires gpf and map files in current directory
# extract grid parameters: spacing, center, points in each direction, total points, coordinate limits
# extract ligand atom types, map names, receptor name
class Grid(object):

	''' Returns a dictionary "grid" with grid parameters (spacing, ...)
		and ligand atom types, map names and receptor name
	'''

	# execute the functions defined below	
	def __init__(self, gpf_file):
		
		self.grid = self.parse_file(gpf_file)
		self.grid = self.additional(self.grid)
	
	# get spacing, no. of points, center of the grid + ligand atom types and map names
	def parse_file(self, gpf_file):

		# initiate dictionary
		grid = dict()

		# get ligand atom types and map names from gpf
		# extract one of the energy maps as a reference to get spacing, no. of points, ...
		with open(gpf_file) as f:
			line = f.readline().split()
			count_maps = 0
			while line:
				# get ligang atom types
				# consider line with ligand types having a comment (#) or not
				if line[0] == 'ligand_types':
					if '#' in line:
						last_index = line.index('#')
						ligand_types = line[1:last_index]
					else:
						ligand_types = line[1:]
				# get map names for each atom type
				elif line[0] == 'map':
					grid.update({'map_' + ligand_types[count_maps] : line[1]})
					count_maps += 1
					ref_map = line[1]
				line = f.readline().split()
		
		grid.update({'ligand_types' : ligand_types})
		grid.update({'ref_map' : ref_map})

		if not ref_map in [f for f in os.listdir('.') if f.endswith('.map')]:
			raise ValueError('There should be a map file "' + ref_map + '" in the current directory according to gpf file "' + gpf_file + '".')
		
		# open reference energy map to get spacing, no. of points, ...
		with open(ref_map) as f:
			
			# iterate from 0 to 5 (6 times, the variable "_" is not used)
			for _ in range(6):
				
				# read the first 6 lines of the energy map
				# for each line generate a list "line" with each space separated string as an element
				line = f.readline().split()
				
				# add to "grid" dictionary the spacing 
				if line[0] == 'SPACING':
					grid.update( {'spacing' : float(line[1]) })
				
				# add to "grid" dictionary the no. of points in each direction
				# the real no. of points is always odd (NELEMENT + 1 in each direction)
				if line[0] == 'NELEMENTS':
					grid.update( {'NX' : int(line[1]) + 1 })
					grid.update( {'NY' : int(line[2]) + 1 })
					grid.update( {'NZ' : int(line[3]) + 1 })

				# add to "grid" dictionary the coordinates of the center of the grid
				if line[0]=='CENTER':
					grid.update( {'CX' : float(line[1]) })
					grid.update( {'CY' : float(line[2]) })
					grid.update( {'CZ' : float(line[3]) })

				# add to "grid" dictionary the file name of the receptor
				if line[0] == 'MACROMOLECULE':
					grid.update( {'file' : line[1] })
					maps_basename = grid['file'].split('.pdbqt')[0]
					grid.update( {'basename' : maps_basename})
		
		# return the dictionary with the grid parameters
		return(grid)
	
	# get additional grid parameters
	def additional(self, grid):
		
		# coordinate limits of the grid (M = max, m = min)
		grid.update({'Mx' : grid['CX'] + (grid['NX'] - 1) * grid['spacing'] / 2 })
		grid.update({'My' : grid['CY'] + (grid['NY'] - 1) * grid['spacing'] / 2 })
		grid.update({'Mz' : grid['CZ'] + (grid['NZ'] - 1) * grid['spacing'] / 2 })
		
		grid.update({'mx' : grid['CX'] - (grid['NX'] - 1) * grid['spacing'] / 2 })
		grid.update({'my' : grid['CY'] - (grid['NY'] - 1) * grid['spacing'] / 2 })
		grid.update({'mz' : grid['CZ'] - (grid['NZ'] - 1) * grid['spacing'] / 2 })
		
		# total no. of grid points
		self.grid.update( {'total_points': self.grid['NX'] * self.grid['NY'] * self.grid['NZ'] })

		# return the dictionary with more grid parameters
		return(grid)


# requires map file in current directory
# extract grid parameters: spacing, center, points in each direction, total points, coordinate limits
class Grid_Only(object):

	''' Returns a dictionary "grid" with grid parameters (spacing, ...)
	'''

	# execute the functions defined below	
	def __init__(self, map_file):
		
		self.grid = self.parse_file(map_file)
		self.grid = self.additional(self.grid)
	
	# get spacing, no. of points, center of the grid
	def parse_file(self, map_file):

		# initiate dictionary
		grid = dict()

		if not map_file in [f for f in os.listdir('.') if f.endswith('.map')]:
			raise ValueError('There should be a map file "' + map_file + '" in the current directory."')
		
		# open reference energy map to get spacing, no. of points, ...
		with open(map_file) as f:
			
			# iterate from 0 to 5 (6 times, the variable "_" is not used)
			for _ in range(6):
				
				# read the first 6 lines of the energy map
				# for each line generate a list "line" with each space separated string as an element
				line = f.readline().split()
				
				# add to "grid" dictionary the spacing 
				if line[0] == 'SPACING':
					grid.update( {'spacing' : float(line[1]) })
				
				# add to "grid" dictionary the no. of points in each direction
				# the real no. of points is always odd (NELEMENT + 1 in each direction)
				if line[0] == 'NELEMENTS':
					grid.update( {'NX' : int(line[1]) + 1 })
					grid.update( {'NY' : int(line[2]) + 1 })
					grid.update( {'NZ' : int(line[3]) + 1 })

				# add to "grid" dictionary the coordinates of the center of the grid
				if line[0]=='CENTER':
					grid.update( {'CX' : float(line[1]) })
					grid.update( {'CY' : float(line[2]) })
					grid.update( {'CZ' : float(line[3]) })

				# add to "grid" dictionary the file name of the receptor
				if line[0] == 'MACROMOLECULE':
					grid.update( {'file' : line[1] })
					maps_basename = grid['file'].split('.pdbqt')[0]
					grid.update( {'basename' : maps_basename})
		
		# return the dictionary with the grid parameters
		return(grid)
	
	# get additional grid parameters
	def additional(self, grid):
		
		# coordinate limits of the grid (M = max, m = min)
		grid.update({'Mx' : grid['CX'] + (grid['NX'] - 1) * grid['spacing'] / 2 })
		grid.update({'My' : grid['CY'] + (grid['NY'] - 1) * grid['spacing'] / 2 })
		grid.update({'Mz' : grid['CZ'] + (grid['NZ'] - 1) * grid['spacing'] / 2 })
		
		grid.update({'mx' : grid['CX'] - (grid['NX'] - 1) * grid['spacing'] / 2 })
		grid.update({'my' : grid['CY'] - (grid['NY'] - 1) * grid['spacing'] / 2 })
		grid.update({'mz' : grid['CZ'] - (grid['NZ'] - 1) * grid['spacing'] / 2 })
		
		# total no. of grid points
		self.grid.update( {'total_points': self.grid['NX'] * self.grid['NY'] * self.grid['NZ'] })

		# return the dictionary with more grid parameters
		return(grid)


# requires grid dictionary and maps in current directory
# create extra map for dummy atom in center of aromatic rings
# (energy = 0 for all points at this stage)
class Extra_Map(object):

	''' Writes map for dummy atom in center of aromatic rings
		(energy = 0 for all points at this stage)
	'''

	def __init__(self, grid):

		self.extramap = self.generate_map(grid)
	
	def generate_map(self, grid):

		with open(grid['basename'] + '.AC.map', 'w') as f:
			with open(grid['ref_map']) as g:
				for _ in range(6):
					line = g.readline()
					f.write(line)						
				for _ in range(grid['NX'] * grid['NY'] * grid['NZ']):
					f.write(u'0.000\n')

		return True


# requires grid and bias sites dictionaries
# get cube representation for the bias sites
# translocate bias site to the grid points
class Cubes(object):
	
	''' Returns bias sites coordinates and props
		once translocated to points in the grid
	'''
	
	def __init__(self, biassites, grid):
		
		cube_set = [self.get_cube(bs, grid) for bs in biassites]

		self.cube_set = [cube for cube in cube_set]
	
	# cube = bias site representation
	def get_cube(self, bs, grid):

		# no. of points in half the side of the cube
		# (twice the size of the bias site to allow smooth energy modifications through space)
		points_from_center = mt.ceil(2 * bs['radius'] / grid['spacing'])

        # no. of points in the side of the cube
        # (+ 1 to have the same no. of points on each side of the central point)
		side_points = points_from_center * 2  + 1

        # total no. of points in the cube
		total_points = side_points ** 3		

        # no. of points between center of the bias site and lowest extreme of the grid
        # (does not include extreme point)
		cnx = round( (bs['x'] - grid['mx']) / grid['spacing'] )
		cny = round( (bs['y'] - grid['my']) / grid['spacing'] )
		cnz = round( (bs['z'] - grid['mz']) / grid['spacing'] )

        # new coordinates for center of the cube coincident with a grid point (translocation)
		cube_center_x = cnx * grid['spacing'] + grid['mx']
		cube_center_y = cny * grid['spacing'] + grid['my']
		cube_center_z = cnz * grid['spacing'] + grid['mz']
		
		# line number in the grid file where the energy for the center of the cube is stored
		# (point of the grid) - 6 lines of metadata + parse of file (x -> y -> z)
		center_point_map = 6 + cnx + (grid['NX']) * cny + (grid['NX']) * (grid['NY']) * cnz

		# define extreme coordinates (minimum) of the cube (matching grid points)
		min_cube_x = cube_center_x - points_from_center * grid['spacing']
		min_cube_y = cube_center_y - points_from_center * grid['spacing']
		min_cube_z = cube_center_z - points_from_center * grid['spacing']

		# define extreme coordinates (maximum) of the cube (matching grid points)
		max_cube_x = cube_center_x + points_from_center * grid['spacing']
		max_cube_y = cube_center_y + points_from_center * grid['spacing']
		max_cube_z = cube_center_z + points_from_center * grid['spacing']

		# no. of points between extreme of the cube and extreme of the grid (minimum)
		min_nx = cnx - points_from_center
		min_ny = cny - points_from_center
		min_nz = cnz - points_from_center

		# check if the cube totally fits inside the grid space
		# if not, remove the points that lie outside and recalculate params
		if grid['Mx'] < max_cube_x:
			total_points = total_points - (max_cube_x - grid['Mx']) / grid['spacing']
			max_cube_x = grid['Mx']
		if grid['My'] < max_cube_y:
			total_points = total_points - (max_cube_y - grid['My']) / grid['spacing']
			max_cube_y = grid['My']
		if grid['Mz'] < max_cube_z:
			total_points = total_points - (max_cube_z - grid['Mz']) / grid['spacing']
			max_cube_z = grid['Mz']
		if grid['mx'] > min_cube_x:
			total_points = total_points - (grid['mx'] - min_cube_x) / grid['spacing']
			min_cube_x = grid['mx']
			min_nx = cnx - (cube_center_x - min_cube_x) / grid['spacing']
		if grid['my'] > min_cube_y:
			total_points = total_points - (grid['my'] - min_cube_y) / grid['spacing']
			min_cube_y = grid['my']
			min_ny = cny - (cube_center_y - min_cube_y) / grid['spacing']
		if grid['mz'] > min_cube_z:
			total_points = total_points - (grid['mz'] - min_cube_z) / grid['spacing']
			min_cube_z = grid['mz']
			min_nz = cnz - (cube_center_z - min_cube_z) / grid['spacing']

		# distance between original center of the bias site and
		# new center of the cube coincident with the grid point
		distance_to_grid = ( (bs['x'] - cube_center_x)**2 + (bs['y'] - cube_center_y)**2 + (bs['z'] - cube_center_z)**2 ) ** (1/2) 

		# dictionary with cube parameters
		cube_props = {'grid_limit' : {'NX': grid['NX'], 'NY': grid['NY'], 'NZ': grid['NZ']},
					 'bias_site_center': {'x': bs['x'], 'y': bs['y'], 'z': bs['z']},
					 'cube_center': {'x': cube_center_x, 'y': cube_center_y, 'z': cube_center_z},
					 'min_cube_coords': {'x': min_cube_x, 'y': min_cube_y, 'z': min_cube_z},
					 'max_cube_coords': {'x': max_cube_x, 'y': max_cube_y, 'z':max_cube_z},
					 'points_cubeMin_gridMin': {'x': min_nx, 'y': min_ny, 'z': min_nz },
					 'center_point_map': int(center_point_map),
					 'spacing': grid['spacing'],
					 'total_points': int(total_points),
					 'radius': bs['radius'],
					 'distance_to_grid': distance_to_grid,
					 'int_type': bs['int_type']
					 }

		# save energy or pfp in dictionary
		if 'DG' in bs:
			cube_props.update( {'DG': bs['DG']} )
		elif 'PFP' in bs:
			cube_props.update( {'PFP': bs['PFP']} )
		else:
			print 'Energy / PFP not found for the bias site ', cube_props['bias_site_center']
		
		# return cube parameters for representing the bias site (translocated to the grid points)
		return cube_props


# requires cubes and grid dictionaries
# modifies the corresponding energy maps
class Insert_Bias(object):

	''' routine that inserts Gaussian modifications
		to the maps according to type of bias
	'''

	def __init__(self, cubes, grid, files):
	
		cubes_acc = []
		cubes_don = []
		cubes_aro = []
		cubes_map = []
		modifications_by_cubes_acc = []
		modifications_by_cubes_don = []
		modifications_by_cubes_aro = []
		modifications_by_cubes_map = []
		
		for cube in cubes:
			
			# store cubes and energy modifications for acceptor biases
			if cube['int_type'] == 'acc' :
				cubes_acc.append(cube)
				modifications_by_cubes_acc.append(self.cube_modifications(cube))

			# store cubes and energy modifications for donor biases
			elif cube['int_type'] == 'don' :
				cubes_don.append(cube)
				modifications_by_cubes_don.append(self.cube_modifications(cube))

			# store cubes and energy modifications for aromatic biases
			elif cube['int_type'] == 'aro' :
				cubes_aro.append(cube)
				modifications_by_cubes_aro.append(self.cube_modifications(cube))

			# store cubes and energy modifications for specific map biases
			elif cube['int_type'] == 'map' :
				cubes_map.append(cube)
				modifications_by_cubes_map.append(self.cube_modifications(cube))

		# for acceptor biases
		if modifications_by_cubes_acc:
			# list of energy modifications (#point, energy)
			modifications_lst_acc = self.merge_boxes(modifications_by_cubes_acc)
			# modify the original acceptor maps with the bias information
			if 'map_OA' in grid:
				self.insert_bias(cubes_acc, modifications_lst_acc, grid['map_OA'], grid['map_OA'].split('.map')[0] + '.biased.map')
			if 'map_NA' in grid:
				self.insert_bias(cubes_acc, modifications_lst_acc, grid['map_NA'], grid['map_NA'].split('.map')[0] + '.biased.map')

		# for donor biases
		if modifications_by_cubes_don:
			# list of energy modifications (#point, energy)
			modifications_lst_don = self.merge_boxes(modifications_by_cubes_don)
			# modify the original donor map with the bias information
			if 'map_HD' in grid:
				self.insert_bias(cubes_don, modifications_lst_don, grid['map_HD'], grid['map_HD'].split('.map')[0] + '.biased.map')

		# for aromatic biases
		if modifications_by_cubes_aro:
			# list of energy modifications (#point, energy)
			modifications_lst_arom = self.merge_boxes(modifications_by_cubes_aro)
			# modify the dummy map with the bias information
			self.insert_bias(cubes_aro, modifications_lst_arom, grid['basename'] + '.AC.map', grid['basename'] + '.AC.biased.map')
			# remove the map with 0 energy created for this sake
			os.remove(grid['basename'] + '.AC.map')

		# for specific map biases
		if modifications_by_cubes_map:
			# list of energy modifications (#point, energy)
			modifications_lst_map = self.merge_boxes(modifications_by_cubes_map)
			# modify the specific map with the bias information
			self.insert_bias(cubes_map, modifications_lst_map, files['map_file'], files['map_file'].split('.map')[0] + '.biased.map')
			
	# returns a list of tuples: map positions in the map file (line number), their energy reward (dE)
	def cube_modifications(self, cube):

		# initiate empty list
		modifications = []
		
		# energy at 298.15 K (kcal/mol)
		kT = 0.59248

		# extreme coordinates (minimum) of the cube (matching grid points)
		x = cube['min_cube_coords']['x']
		y = cube['min_cube_coords']['y']
		z = cube['min_cube_coords']['z']

		# no. of points between extreme of the cube and extreme of the grid (minimum)
		nx = cube['points_cubeMin_gridMin']['x']
		ny = cube['points_cubeMin_gridMin']['y']
		nz = cube['points_cubeMin_gridMin']['z']

		# iterate times equal to total no. of points in the cube
		for _ in range(cube['total_points']):

			# if x > extreme coordinates (maximum) of the cube (matching grid points)
			# move one y point and start again parsing x points
			if x > cube['max_cube_coords']['x']:
				x =  cube['min_cube_coords']['x']
				y += cube['spacing']
				
				nx = cube['points_cubeMin_gridMin']['x']
				ny += 1
			
			# if y > extreme coordinates (maximum) of the cube (matching grid points)
			# move one z point and start again parsing x,y points
			if y > cube['max_cube_coords']['y']:
				y =  cube['min_cube_coords']['y']
				z += cube['spacing']
				
				ny = cube['points_cubeMin_gridMin']['y']
				nz += 1

			# distance from center of bias site
			distance = ( (x - cube['bias_site_center']['x'])**2 + (y - cube['bias_site_center']['y'])**2 + (z - cube['bias_site_center']['z'])**2 ) ** (1/2)
			
			# if the energy input for the bias site was PFP (positive number)
			if 'PFP' in cube:
				
				# define energy reward (gaussian function for distance from bias site center)
				dE = -1 * kT * mt.log(cube['PFP']) * mt.exp(-1 * (distance ** 2) / (cube['radius'] ** 2))
			
			# if the energy input for the bias site was DG (negative number)
			elif 'DG' in cube:
				
				# define energy reward (gaussian function for distance from bias site center)
				dE = cube['DG'] * mt.exp(-1 * (distance ** 2) / (cube['radius'] ** 2))	

			# if the energy reward is stronger than 0.01 kcal/mol
			if dE < -0.01:
				
				# get line from map file to change the energy in that point (parse only inside the cube)
				map_position = 6 + nx + (cube['grid_limit']['NX']) * ny + (cube['grid_limit']['NX']) * (cube['grid_limit']['NY']) * nz
			
				# store line from map file and energy modification in "modifications"
				modifications.append((map_position, dE))

			# move one x point in grid
			x += cube['spacing']
			nx += 1
		
		# returns a list of map positions in the map file with their energy reward (dE)
		return modifications

	# if the same point of the grid has 2 or more modifications due to different
	# bias sites, only keep the most favorable one
	def merge_boxes(self, modifications):
		
		if len(modifications) > 0:
			
			# transform into "list of tuples" with [point index, energy]
			flat_list = list(itertools.chain.from_iterable(modifications))
			
			# sort list first by point index and then by energy
			flat_list.sort(key=itemgetter(0,1))
			
			# save the lowest energy for each point
			# generate a diccionary {point index: lowest energy}
			lst_unique=[(key,) + tuple(elem for _, elem in group) for key, group in groupby(flat_list, lambda pair: pair[0])]
			energy_dictionary = {}
			for i in lst_unique:
				energy_dictionary[i[0]] = i[1]
		
		# return dictionary {point index: lowest energy modification}
		return energy_dictionary

	# modify the corresponding map
	def insert_bias(self, cube_set, modifications, map_in, map_out):
		
		# initialize dictionary
		Eori = {}
		
		# open map
		with open(map_in) as f:
			with open(map_out, 'w') as g:
				for index_map_in, original_E in enumerate(f.readlines()):
					# find point in dictionary of modifications
					if index_map_in in modifications:
						# define original energy
						Eori[index_map_in] = float(original_E.rstrip('\n'))
						# define new energy with bias modification
						E = Eori[index_map_in] + modifications[index_map_in]
						# write new energy value
						g.write(u'{:.3f}'.format(E) + '\n')
					else:
						# keep original energy value
						g.write(original_E)
						
		# output with bias information (grid_modif_"ATOMYTPE".map)
		dat_modif = 'grid_modif_' + map_in.split('.')[-2] +'.dat'
		lines = []
		with open(dat_modif, 'w') as p:
			p.write(u"\t".join(['x', 'y', 'z', 'point', 'dist', 'E_ori', 'E_modif', 'E_delta', 'radius'])+'\n')
			for cube in cube_set:
				line = ["%0.3f" % cube['cube_center']['x'],
					    "%0.3f" % cube['cube_center']['y'], 
						"%0.3f" % cube['cube_center']['z'],
						str(cube['center_point_map']), 
						"%0.3f" % cube['distance_to_grid'],
						"%0.3f" % Eori[cube['center_point_map']],
						"%0.3f" % (float(Eori[cube['center_point_map']]) + float(modifications[cube['center_point_map']])), 
						"%0.3f" % float(modifications[cube['center_point_map']]), 
						"%0.3f" % cube['radius']]
				lines.append("\t".join(line))
			for line in lines:
				p.write(unicode(line) + '\n')


# requires dpf file in current directory
# parses dpf to extract info and generate a new modified dpf
class DPF(object):

	def __init__(self, dpf_file):
		
		self.dpf = self.parse_file(dpf_file)
		
	# parse dpf and extract ligand types, map names, ligand name, ...
	def parse_file(self, dpf_file):

		# initiate dictionary
		dpf_props = dict()

		# parse dpf and extract ligand types, map names, ligand name, ...
		with open(dpf_file) as f:
			line = f.readline().split()
			count_maps = 0
			while line:
				if line[0] == 'ligand_types':
					if '#' in line:
						dpf_props.update({'#' : 'yes'})
						last_index = line.index('#')
						ligand_types = line[1:last_index]
					else:
						dpf_props.update({'#' : 'no'})
						ligand_types = line[1:]
					dpf_props.update({'no_ligand_types' : len(ligand_types)})
				elif line[0] == 'map':
					dpf_props.update({ligand_types[count_maps] : line[1]})
					count_maps += 1
				elif line[0] == 'move':
					dpf_props.update({'ligand_name' : line[1]})
					basename = line[1].split('.pdbqt')[0]
					dpf_props.update({'ligand_basename' : basename})
				line = f.readline().split()
			
			acceptors = []
			acceptor_maps = []
			if 'OA' in dpf_props:
				acceptors.append('OA')
				acceptor_maps += dpf_props['OA'].split()
			if 'NA' in dpf_props:
				acceptors.append('NA')
				acceptor_maps += dpf_props['NA'].split()
			dpf_props.update({'acceptors' : acceptors})
			dpf_props.update({'acceptor_maps' : acceptor_maps})
			
			donors = []
			donor_maps = []
			if 'HD' in dpf_props:
				donors.append('HD')
				donor_maps += dpf_props['HD'].split()
			dpf_props.update({'donors' : donors})
			dpf_props.update({'donor_maps' : donor_maps})
			
			if 'A' in dpf_props:
				map_basename = dpf_props['A'].split('A.map')[0]
			elif 'C' in dpf_props:
				map_basename = dpf_props['C'].split('C.map')[0]
			elif 'OA' in dpf_props:
				map_basename = dpf_props['OA'].split('OA.map')[0]
			elif 'N' in dpf_props:
				map_basename = dpf_props['N'].split('N.map')[0]
			dpf_props.update({'map_basename' : map_basename})
		
		return dpf_props


# create a new dpf with new maps, dummy type, parameter file, ...
class New_DPF(object):

	def __init__(self, dpf_props, dpf_file, has_dummy, map_file):
		
		self.new_dpf = self.create_dpf(dpf_props, dpf_file, has_dummy, map_file, acc_bias, don_bias, aro_bias, map_bias)
		
	# create a new dpf with new maps, dummy type, parameter file, ...		
	def create_dpf(self, dpf_props, dpf_file, has_dummy, map_file, acc_bias, don_bias, aro_bias, map_bias):

		# new dpf name: basename.biased.dpf
		new_dpf_file = dpf_file.split('.dpf')[0] + '.biased.dpf'

		# start writing new dpf
		with open(new_dpf_file, 'w') as f:

			# if there is aromatic bias and aromatic atoms in the ligand
			# add first line to call parameter file with dummy atom parameters
			if aro_bias and has_dummy:
				f.write(u'parameter_file ad4_arom_params.dat\n')

			# open original dpf
			with open(dpf_file) as g:
				line = g.readline()
				line_lst = line.split()
				
				map_count = 0
				while line:
					
					# add dummy atom type at the end of ligand_types list if
					# there is aromatic bias and aromatic atoms in the ligand
					if line_lst[0] == 'ligand_types' and aro_bias and has_dummy:

						if dpf_props['#'] == 'yes':
							last_index = line_lst.index('#')
							new_line_lst = line_lst[0:last_index]
							new_line_lst.append('AC')
							new_line = ' '.join(str(e) for e in new_line_lst)
							f.write(unicode(new_line) + '\n')
						
						elif dpf_props['#'] == 'no':
							f.write(line + ' AC\n')
					
					# modify map names if corresponds	
					elif line_lst[0] == 'map':
						map_count += 1

						# modify specific map name if corresponds
						if map_bias:
							if line_lst[1] == map_file:
								basename = line_lst[1].split('.map')[0]
								f.write('map ' + basename + '.biased.map\n')
							else:
								f.write(line)
						
						# modify acceptor / donor map names if corresponds	
						else:
							
							if 'OA' in dpf_props:
								if line_lst[1] == dpf_props['OA']:
									if acc_bias:
										basename = line_lst[1].split('.map')[0]
										f.write('map ' + basename + '.biased.map\n')
									else:
										f.write(line)
							
							if 'NA' in dpf_props:
								if line_lst[1] == dpf_props['NA']:
									if acc_bias:
										basename = line_lst[1].split('.map')[0]
										f.write('map ' + basename + '.biased.map\n')
									else:
										f.write(line)
							
							if 'HD' in dpf_props:
								if line_lst[1] == dpf_props['HD']:
									if don_bias:
										basename = line_lst[1].split('.map')[0]
										f.write('map ' + basename + '.biased.map\n')
									else:
										f.write(line)
							
							if line_lst[1] not in dpf_props['acceptor_maps'] and line_lst[1] not in dpf_props['donor_maps']:
								f.write(line)
					
					# add dummy map at the end if corresponds
					elif map_count == dpf_props['no_ligand_types'] and aro_bias and has_dummy:
						f.write('map ' + dpf_props['map_basename'] + 'AC.biased.map\n')
						f.write(line)
						map_count = 0
					
					# change ligand name to dock if dummy atoms were added to its structure	
					elif line_lst[0] == 'move' and aro_bias and has_dummy:
						f.write('move ' + dpf_props['ligand_basename'] + '.dum.pdbqt\n')
					
					# otherwise keep lines from original dpf
					else:
						f.write(line)
					
					line = g.readline()
					line_lst = line.split()
		
		return dpf_props


# Adds dummy atom to each aromatic ring center if these two conditions are satisfied:
# - Open Babel recognize the ring as aromatic.
# - The ring has at least one aromatic carbon according to the input PDBQT.
class Dummies(object):
	
	def __init__(self, ligand_filename):
	
		self.add_dummies = self.add_dummy(ligand_filename)
	
	# function to get center of mass
	def center_of_mass(self, Vx, Vy, Vz, mass):
		cgx = np.sum(Vx*mass)/np.sum(mass)
		cgy = np.sum(Vy*mass)/np.sum(mass)
		cgz = np.sum(Vz*mass)/np.sum(mass)
		return (cgx, cgy, cgz)

	# pdbqt parser
	class PDBQT(object):
		
		def __init__(self, line):
			self._parse_common(line)     # useful for PDB and PDBQT
			self._parse_specific(line)   # differs in PDBQT

		def getline(self):
			txt = self._print_common()    # no \n; PDB + PDBQT
			txt += self._print_specific() # differs in PDBQT
			return txt
			
		def _parse_common(self, line):
			"""Common to PDB and PDBQT formats"""
			self.keyword     = line      [ 0: 6]     # ATOM or HETATM
			self.serial      = int(line  [ 6:11])    # atom id
			#                            [11:12]
			self.name        = line      [12:16]     # atom name
			self.altLoc      = line      [16:17]     # Alternate location
			self.resName     = line      [17:20]     # Residue name
			#                            [20:21] 
			self.chain       = line      [21:22]     # chain
			self.resNum      = int(line  [22:26])    # Residue number
			self.icode       = line      [26:27]     # ???
			#                            [27:30]
			self.x           = float(line[30:38])    # X
			self.y           = float(line[38:46])    # Y
			self.z           = float(line[46:54])    # Z
			self.occupancy   = float(line[54:60])    # Occupancy
			self.bfact       = float(line[60:66])    # Temperature factor

		def _parse_specific(self, line):
			""" PDBQT characters [68:79] """
			self.charge      = float(line[68:76])   # Charge
			self.atype       = line      [77:79]    # Atom type
			#self.atype = self.atype.strip().upper()
			self.atype = self.atype.strip()

		def _print_common(self):
			""" Characters [0:68]"""
			linestr = ''
			linestr += '%6s' % (self.keyword)
			linestr += '%5d' % (self.serial)
			linestr += ' ' 
			linestr += '%4s' % (self.name)
			linestr += '%1s' % (self.altLoc) 
			linestr += '%3s' % (self.resName)
			linestr += ' ' 
			linestr += '%1s' % (self.chain)
			linestr += '%4d' % (self.resNum)
			linestr += '%1s' % (self.icode)
			linestr += ' ' * 3 
			linestr += '%8.3f' % (self.x)
			linestr += '%8.3f' % (self.y)
			linestr += '%8.3f' % (self.z)
			linestr += '%6.2f' % (self.occupancy)
			linestr += '%6.2f' % (self.bfact)
			return linestr

		def _print_specific(self):
			""" PDBQT characters [68:79] """
			linestr =  ' ' * 2                      # [66:68]
			linestr += '%8.3f' % (self.charge)      # [68:76]
			linestr += ' ' * 1                      # [76:77]
			linestr += '%-2s' % (self.atype)        # [77:79]
			#linestr += '\n'
			return linestr

	# get total number of atoms from a pdbqt file
	def NumberOfAtoms(self, pdbqtFile):
		natoms=0
		with open(pdbqtFile, 'r') as f:
			# iterate over file lines
			for line in f:
				# count atom if line starts with ATOM or HETATM
				if line.startswith('ATOM  ') or line.startswith('HETATM'):
					natoms+=1
		return natoms

	# store if block has aromatic rings according to input pdbqt
	class Block(object):
		""" Class to provide print block and define if it aromatic
		according to the pdbqt input file
		"""
		coordinates = []
		# prints whole pdbqt block
		def printBlock(self):
			txt = ""
			if len(self.coordinates) > 0:
				for line in self.coordinates:
					atom = Dummies.PDBQT(line)
					txt += atom.getline() + "\n"
			return str(txt)
			
		# defines if block is aromatic (i.e. it has at least one A atom)
		def isaromatic(self):
			AList = []
			if len(self.coordinates) > 0:
				for line in self.coordinates:
					atom = Dummies.PDBQT(line)
					AList.append(atom.atype)
				if "A" in AList:
								return True
				else:
					return False
			else:
				return False

	# find aromatic rings (according to obabel) and return a list of their centers of mass
	def return_COMs(self, inputstring):
		""" return one or more centers of mass for aromatic ring in a unique block
		"""
		# set a new object to be manipulated
		mol = openbabel.OBMol()
		obConversion = openbabel.OBConversion()
		# input and output be the same format (pdbqt)
		obConversion.SetInAndOutFormats("pdbqt", "pdbqt")
		obConversion.ReadString(mol, inputstring)
		# set the output variable
		out = []
		# iterate over the rings present in the molecule (block)
		for r in mol.GetSSSR():
			# save if aromatic
			if r.IsAromatic():
				Vx = np.array([])
				Vy = np.array([])
				Vz = np.array([])
				mass = np.array([])
				# iterate over the atoms
				for obatom in openbabel.OBMolAtomIter(mol):
					# only atoms defined in the present ring
					if  r.IsInRing(obatom.GetIdx()):
						Vxcurrent=np.array([obatom.GetX()])
						Vx=np.concatenate((Vx,Vxcurrent))
						Vycurrent=np.array([obatom.GetY()])
						Vy=np.concatenate((Vy,Vycurrent))
						Vzcurrent=np.array([obatom.GetZ()])
						Vz=np.concatenate((Vz,Vzcurrent))
						mass_current=np.array([obatom.GetAtomicMass()])	
						mass=np.concatenate((mass,mass_current))
				# return aromatic center to the current ring
				center = self.center_of_mass(Vx,Vy,Vz,mass)
				# save center of mass in a list
				out.append(center)
		return out

	# format center of mass as a pdbqt line
	def return_COM_format(self, resname, resnum, serial, COMlist, chain, keyword):
		''' Returns pdbqt lines for centers of mass
		'''
		# example of pdbqt line format
		lineexample="ATOM      1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00     0.000 AC"
		# parse example pdbqt line
		atom = self.PDBQT(lineexample)
		# out variable
		out = []
		# iterate over the different centers of mass that were found
		for center in COMlist:
			new_line = atom.getline()
			add_line = self.PDBQT(new_line)
			# get atom number
			add_line.serial = serial
			# get coordinates
			add_line.x = center[0]
			add_line.y = center[1]
			add_line.z = center[2]
			# get atom name, type, residue name, number, etc
			add_line.name = str('CM ')
			add_line.resNum = resnum
			add_line.atype = 'AC'
			add_line.resName = resname
			add_line.chain = chain
			add_line.keyword = keyword
			# print pdbqt line for COM
			out.append(str(add_line.getline()))
			serial+=1
		return out

	# main function
	def parse_tree(self, pdbqtFile):
		# get number of atoms in the pdbqt file
		natoms = self.NumberOfAtoms(pdbqtFile)
		# variable to store the different blocks of the pdbqt file
		blockstring = []
		# out variable
		out = ""
		# aromatic flag
		arom_flag = False
		# open pdbqt input file
		with open(pdbqtFile, 'r') as f:
			# iterate over file lines
			for line in f:
				# if line starts with ATOM or HETATM
				if line.startswith('ATOM  ') or line.startswith('HETATM'):
					# parse line
					atom = self.PDBQT(line)
					# get keyword (ATOM or HETATM), residue name, number and chain
					keyword = atom.keyword
					resname = atom.resName
					resnum = atom.resNum
					chain = atom.chain
					# store line in block
					blockstring.append(line)
				# when block ends
				else:
					# if it was an actual block (e.g. not initial remarks)
					if len(blockstring) > 0:
						# define a new block
						my_block = self.Block()
						# save the pdbqt lines (with coordinates) of the new block
						my_block.coordinates = blockstring
						# if block has A atom according to input pdbqt
						if my_block.isaromatic() is True:
							# get center(s) of mass for aromatic rings in this block (obabel defines aromaticity at this point)
							COMlist = self.return_COMs(my_block.printBlock())
							# if there is at least one aromatic ring found, turn on aromatic flag
							if COMlist:
								arom_flag = True
							# format center(s) of mass as pdbqt lines
							new_lines = self.return_COM_format(resname,resnum,natoms+1,COMlist,chain,keyword)
							# update block with COM dummy atoms
							my_block.coordinates = my_block.coordinates + new_lines
							# update the total number of atoms
							natoms = natoms + int(len(COMlist))
						# print block to out variable
						out += (my_block.printBlock())
					# print root, branch, endbranch, endroot lines to out variable
					out+=line
					# reset blockstring
					blockstring = []
		# remove final empty line
		out = os.linesep.join([s for s in out.splitlines() if s])
		# if there is at least one aromatic ring, print new pdbqt file
		if arom_flag:
			return out
		# else return False
		else:
			print 'WARNING. No aromatic rings found in', pdbqtFile, '/ Skipping file.'
			return arom_flag

	# modify pdbqt adding a dummy atom in the center of aromatic rings
	# using the above defined functions
	def add_dummy(self, input_filename):
		
		has_dummy = False
		
		# define output name
		ligand_basename = input_filename.split('.pdbqt')[0]
		output = ligand_basename + ".dum.pdbqt"
		
		# get output pdbqt with dummy atoms in center of aromatic rings
		coordinates = self.parse_tree(input_filename)
		
		# write output only if it has aromatic rings
		if coordinates:
			out = open(output, 'w')
			out.write(coordinates)
			out.close()
			has_dummy = True
		
		return has_dummy


######################################################
###################### Main ##########################
######################################################

# test input arguments and save them in a dict
files = optparser()

# get coordinates, energy, radius and type of the bias sites
biases = Bias_Sites(files['bias_file'])

# establish which type of biases should be applied
acc_bias = False
don_bias = False
aro_bias = False
map_bias = False
for i in range(0,len(biases.biassites)):
	if biases.biassites[i]['int_type'] == 'acc' and not acc_bias:
		acc_bias = True
	elif biases.biassites[i]['int_type'] == 'don' and not don_bias:
		don_bias = True
	elif biases.biassites[i]['int_type'] == 'aro' and not aro_bias:
		aro_bias = True
	elif biases.biassites[i]['int_type'] == 'map' and not map_bias:
		map_bias = True

# if -m invoked, generate new modified map from original specific map and bias info
if 'map_file' in files:
	
	# specific map modification may not be combined with other biases in the same run
	if acc_bias or don_bias or aro_bias:
		print 'Specific map modification (-m) may not be combined with other biases in the same run.\n'
		print 'The bias parameter file',files['bias_file'],'should only have "map" type biases.'
		sys.exit(2)
	
	# get grid parameters: spacing, center, no. of points, ...
	# get ligand atom types, map names, receptor name, ...
	mygrid = Grid_Only(files['map_file'])

	# get cube representation for the bias sites (translocated to the grid points)
	cube_set = Cubes(biases.biassites, mygrid.grid)
	
	# generate new map with modified energies
	modifications = Insert_Bias(cube_set.cube_set, mygrid.grid, files)
	
	# WARNING if no dpf was indicated
	if 'dpf_file' not in files:
		map_out = files['map_file'].split('.map')[0]
		print '>>>>>>>:',files['map_file'],'has been modified to',map_out,'.biased.map'
		print 'WARNING: No DPF file was indicated'
		print 'WARNING: Remember to modify the map line in the DPF file accordingly.'

# if -g invoked, generate new modified maps from gpf and bias info
if 'gpf_file' in files:

	# get grid parameters: spacing, center, no. of points, ...
	# get ligand atom types, map names, receptor name, ...
	mygrid = Grid(files['gpf_file'])

	# generate extra map for aromatic bias (with zero energy at this point)
	if aro_bias:
		Extra_Map(mygrid.grid)

	# get cube representation for the bias sites (translocated to the grid points)
	cube_set = Cubes(biases.biassites, mygrid.grid)
	
	# generate new maps with modified energies
	modifications = Insert_Bias(cube_set.cube_set, mygrid.grid, files)

# if -d invoked, generate new dpf and, if necessary, new ligand pdbqt with dummies
if 'dpf_file' in files:

	# parse dpf file to get properties
	mydpf = DPF(files['dpf_file'])
	
	# if there is a specific map bias
	if map_bias:
		# generate new dpf file to read the new map
		mynewdpf = New_DPF(mydpf.dpf, files['dpf_file'], False, files['map_file'])

	# if there is aromatic bias
	elif aro_bias:
		
		import openbabel, pybel
		from itertools import chain
		from os.path import basename
		import numpy as np
		
		# add dummy to ligand if necessary
		mydum = Dummies(mydpf.dpf['ligand_name'])
		# generate new dpf file with or without dummy atom according to ligand
		mynewdpf = New_DPF(mydpf.dpf, files['dpf_file'], mydum.add_dummies, False)
	
	# if there is not aromatic bias
	else:	
		# generate new dpf file without dummy atom
		mynewdpf = New_DPF(mydpf.dpf, files['dpf_file'], False, False)
		
# if -D invoked, generate new dpfs and, if necessary, new ligand pdbqts with dummies
# for all the dpfs/pdbqts in the specified directory
if 'dpfs_directoy' in files:
	
	# change to specified directory
	os.chdir(files['dpfs_directoy'])
	# list dpfs
	dpf_files = [f for f in os.listdir('.') if f.endswith('.dpf')]
	if len(dpf_files) == 0:
		raise ValueError('There should be at least one dpf file in the -D directory.')
	else:
		
		if aro_bias:
			import openbabel, pybel
			from itertools import chain
			from os.path import basename
			import numpy as np
		
		# generate new dpf files and add dummies to ligands if necessary
		for i in dpf_files:
			mydpf = DPF(i)
			
			# if there is aromatic bias
			if aro_bias:
				if os.path.isfile(mydpf.dpf['ligand_name']):
					# add dummy to ligand if necessary
					mydum = Dummies(mydpf.dpf['ligand_name'])
					# generate new dpf file with or without dummy atom according to ligand
					mynewdpf = New_DPF(mydpf.dpf, i, mydum.add_dummies, False)
				else:
					raise ValueError('Could not add dummy atom to ligand structure.'
					' There should be a ligand pdbqt file named "' + mydpf.dpf['ligand_name'] +
					'" in the -D directory according to "' + i + '".')
			# if there is not aromatic bias
			else:
				# generate new dpf file without dummy atom
				mynewdpf = New_DPF(mydpf.dpf, i, False, False)

# create parameter file for dummy atom if necessary
if aro_bias:
	f = open('ad4_arom_params.dat', 'w')
	f.write(u'# $Id: prepare_bias.py,v 1.1.2.1 2018/12/18 01:02:17 annao Exp $\n'
			'# \n'
			'# AutoDock \n'
			'# \n'
			'# Copyright (C) 1989-2007,  Garrett M. Morris, David S. Goodsell, Ruth Huey, Arthur J. Olson, \n'
			'# All Rights Reserved.\n'
			'# \n'
			'# AutoDock is a Trade Mark of The Scripps Research Institute.\n'
			'# \n'
			'# This program is free software; you can redistribute it and/or\n'
			'# modify it under the terms of the GNU General Public License\n'
			'# as published by the Free Software Foundation; either version 2\n'
			'# of the License, or (at your option) any later version.\n'
			'# \n'
			'# This program is distributed in the hope that it will be useful,\n'
			'# but WITHOUT ANY WARRANTY; without even the implied warranty of\n'
			'# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n'
			'# GNU General Public License for more details.\n'
			'# \n'
			'# You should have received a copy of the GNU General Public License\n'
			'# along with this program; if not, write to the Free Software\n'
			'# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.\n'
			'\n'
			'# AutoDock Linear Free Energy Model Coefficients and Energetic Parameters\n'
			'#                   Version 4.1 Bound\n'
			'#                    $Revision: 1.1.2.1 $\n'
			'\n'
			'# FE_unbound_model is used to specify how the internal energy of the\n'
			'# ligand should be treated when estimating the free energy of binding,\n'
			'# and can be set to one of the following strings:\n'
			'#   unbound_same_as_bound, extended, or compact\n'
			'# unbound_same_as_bound -- this assumes the internal energy of the ligand is the\n'
			'#                          same before and after binding.\n'
			'# extended -- this assumes the internal energy of the ligand is that of an \n'
			'#             extended conformation when unbound.\n'
			'# compact -- this assumes the internal energy of the ligand is that of a \n'
			'#            compact conformation when unbound.\n'
			'#FE_unbound_model unbound_same_as_bound\n'
			'\n'
			'# AutoDock 4 free energy coefficients with respect to original (AD2) energetic parameters\n'
			'#  This model assumes that the bound and unbound conformations are the same.\n'
			'#  See Table 3 in Huey,Morris,Olson&Goodsell (2007) J Comput Chem 28: 1145-1152.\n'
			'#\n'
			'#               Free Energy Coefficient\n'
			'#               ------\n'
			'FE_coeff_vdW    0.1662\n'
			'FE_coeff_hbond  0.1209\n'
			'FE_coeff_estat  0.1406\n'
			'FE_coeff_desolv 0.1322\n'
			'FE_coeff_tors   0.2983\n'
			'\n'
			'# AutoDock 4 Energy Parameters\n'
			'\n'
			'# - Atomic solvation volumes and parameters\n'
			'# - Unweighted vdW and Unweighted H-bond Well Depths\n'
			'#\n'
			'# - Atom Types\n'
			'# - Rii = sum of vdW radii of two like atoms (in Angstrom)\n'
			'# - epsii = vdW well depth (in Kcal/mol)\n'
			'# - vol = atomic solvation volume (in Angstrom^3)\n'
			'# - solpar = atomic solvation parameter\n'
			'# - Rij_hb = H-bond radius of the heteroatom in contact with a hydrogen (in Angstrom)\n'
			'# - epsij_hb = well depth of H-bond (in Kcal/mol)\n'
			'# - hbond = integer indicating type of H-bonding atom (0=no H-bond)\n'
			'# - rec_index = initialised to -1, but later on holds count of how many of this atom type are in receptor\n'
			'# - map_index = initialised to -1, but later on holds the index of the AutoGrid map\n'
			'# - bond_index = used in AutoDock to detect bonds; see "mdist.h", enum {C,N,O,H,XX,P,S}\n'
			'#\n'
			'# - To obtain the Rij value for non H-bonding atoms, calculate the \n'
			'#        arithmetic mean of the Rii values for the two atom types.\n'
			'#        Rij = (Rii + Rjj) / 2\n'
			'#\n'
			'# - To obtain the epsij value for non H-bonding atoms, calculate the \n'
			'#        geometric mean of the epsii values for the two atom types.\n'
			'#        epsij = sqrt( epsii * epsjj )\n'
			'#\n'
			'# - Note that the Rij_hb value is non-zero for heteroatoms only, and zero for H atoms;\n'
			'#        to obtain the length of an H-bond, look up Rij_hb for the heteroatom only; \n'
			'#        this is combined with the Rii value for H in the receptor, in AutoGrid.\n'
			'#        For example, the Rij_hb for OA-HD H-bonds will be (1.9 + 1.0) Angstrom, \n'
			'#        and the weighted epsij_hb will be 5.0 kcal/mol * FE_coeff_hbond.\n'
			'#\n'
			'#        Atom   Rii                             Rij_hb       rec_index\n'
			'#        Type         epsii           solpar         epsij_hb    map_index\n'
			'#                            vol                          hbond     bond_index\n'
			'#        --     ----  -----  -------  --------  ---  ---  -  --  -- --\n'
			'atom_par H      2.00  0.020   0.0000   0.00051  0.0  0.0  0  -1  -1  3	# Non H-bonding Hydrogen\n'
			'atom_par HD     2.00  0.020   0.0000   0.00051  0.0  0.0  2  -1  -1  3	# Donor 1 H-bond Hydrogen\n'
			'atom_par HS     2.00  0.020   0.0000   0.00051  0.0  0.0  1  -1  -1  3	# Donor S Spherical Hydrogen\n'
			'atom_par C      4.00  0.150  33.5103  -0.00143  0.0  0.0  0  -1  -1  0	# Non H-bonding Aliphatic Carbon\n'
			'atom_par A      4.00  0.150  33.5103  -0.00052  0.0  0.0  0  -1  -1  0	# Non H-bonding Aromatic Carbon\n'
			'atom_par N      3.50  0.160  22.4493  -0.00162  0.0  0.0  0  -1  -1  1	# Non H-bonding Nitrogen\n'
			'atom_par NA     3.50  0.160  22.4493  -0.00162  1.9  5.0  4  -1  -1  1	# Acceptor 1 H-bond Nitrogen\n'
			'atom_par NS     3.50  0.160  22.4493  -0.00162  1.9  5.0  3  -1  -1  1	# Acceptor S Spherical Nitrogen\n'
			'atom_par OA     3.20  0.200  17.1573  -0.00251  1.9  5.0  5  -1  -1  2	# Acceptor 2 H-bonds Oxygen\n'
			'atom_par OS     3.20  0.200  17.1573  -0.00251  1.9  5.0  3  -1  -1  2	# Acceptor S Spherical Oxygen\n'
			'atom_par F      3.09  0.080  15.4480  -0.00110  0.0  0.0  0  -1  -1  4	# Non H-bonding Fluorine\n'
			'atom_par Mg     1.30  0.875   1.5600  -0.00110  0.0  0.0  0  -1  -1  4	# Non H-bonding Magnesium\n'
			'atom_par MG     1.30  0.875   1.5600  -0.00110  0.0  0.0  0  -1  -1  4	# Non H-bonding Magnesium\n'
			'atom_par P      4.20  0.200  38.7924  -0.00110  0.0  0.0  0  -1  -1  5	# Non H-bonding Phosphorus\n'
			'atom_par SA     4.00  0.200  33.5103  -0.00214  2.5  1.0  5  -1  -1  6	# Acceptor 2 H-bonds Sulphur\n'
			'atom_par S      4.00  0.200  33.5103  -0.00214  0.0  0.0  0  -1  -1  6	# Non H-bonding Sulphur\n'
			'atom_par Cl     4.09  0.276  35.8235  -0.00110  0.0  0.0  0  -1  -1  4	# Non H-bonding Chlorine\n'
			'atom_par CL     4.09  0.276  35.8235  -0.00110  0.0  0.0  0  -1  -1  4	# Non H-bonding Chlorine\n'
			'atom_par Ca     1.98  0.550   2.7700  -0.00110  0.0  0.0  0  -1  -1  4	# Non H-bonding Calcium\n'
			'atom_par CA     1.98  0.550   2.7700  -0.00110  0.0  0.0  0  -1  -1  4	# Non H-bonding Calcium\n'
			'atom_par Mn     1.30  0.875   2.1400  -0.00110  0.0  0.0  0  -1  -1  4	# Non H-bonding Manganese\n'
			'atom_par MN     1.30  0.875   2.1400  -0.00110  0.0  0.0  0  -1  -1  4	# Non H-bonding Manganese\n'
			'atom_par Fe     1.30  0.010   1.8400  -0.00110  0.0  0.0  0  -1  -1  4	# Non H-bonding Iron\n'
			'atom_par FE     1.30  0.010   1.8400  -0.00110  0.0  0.0  0  -1  -1  4	# Non H-bonding Iron\n'
			'atom_par Zn     1.48  0.550   1.7000  -0.00110  0.0  0.0  0  -1  -1  4	# Non H-bonding Zinc\n'
			'atom_par ZN     1.48  0.550   1.7000  -0.00110  0.0  0.0  0  -1  -1  4	# Non H-bonding Zinc\n'
			'atom_par Br     4.33  0.389  42.5661  -0.00110  0.0  0.0  0  -1  -1  4	# Non H-bonding Bromine\n'
			'atom_par BR     4.33  0.389  42.5661  -0.00110  0.0  0.0  0  -1  -1  4	# Non H-bonding Bromine\n'
			'atom_par I      4.72  0.550  55.0585  -0.00110  0.0  0.0  0  -1  -1  4	# Non H-bonding Iodine\n'
			'atom_par Z      4.00  0.150  33.5103  -0.00143  0.0  0.0  0  -1  -1  0  # Non H-bonding covalent map\n'
			'atom_par G      4.00  0.150  33.5103  -0.00143  0.0  0.0  0  -1  -1  0	# Ring closure Glue Aliphatic Carbon  # SF\n'
			'atom_par GA     4.00  0.150  33.5103  -0.00052  0.0  0.0  0  -1  -1  0	# Ring closure Glue Aromatic Carbon   # SF\n'
			'atom_par J      4.00  0.150  33.5103  -0.00143  0.0  0.0  0  -1  -1  0	# Ring closure Glue Aliphatic Carbon  # SF\n'
			'atom_par Q      4.00  0.150  33.5103  -0.00143  0.0  0.0  0  -1  -1  0	# Ring closure Glue Aliphatic Carbon  # SF\n'
			'atom_par AC     0.00  0.000   0.0000   0.00000  0.0  0.0  0  -1  -1  4	# dummy atom in center of aromatic rings\n')
