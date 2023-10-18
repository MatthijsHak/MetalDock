#!/usr/bin/python

# import packages
import sys

# test input arguments
def optparser(argv=sys.argv[1:]):
	
	''' Returns "arg" with dlg file name
	'''
	
	# usage: input line
	usage = '''
	Usage:

	python pdb_poses.py docking_log_file
	'''
	
	# check number of arguments and return DLG file name
	if len(argv) != 1:
		print 'Error: invalid number of arguments'
		print usage
		sys.exit(2)
	else:
		arg = argv[0]
	
	return arg

# test input arguments and save the input dlg file name
input_file = optparser()

# open dlg
with open(input_file) as dlg:
	# read line
	line = dlg.readline()
	# create list with line elements
	line_lst = line.split()
	# initiate cluster number
	cluster_number = 0
	# parse DLG
	while line:
		# if line is not empty
		if len(line_lst) !=0:
			# create a new pdb file each time MODEL line is found
			if line_lst[0] == 'MODEL':
				cluster_number += 1
				pdb = open('rank%s.pdb' % cluster_number, 'w')
			# write ATOM lines to output PDB
			elif line_lst[0] == 'ATOM':
				pdb.write(line)
		# go to next line in DLG
		line = dlg.readline()
		line_lst = line.split()
