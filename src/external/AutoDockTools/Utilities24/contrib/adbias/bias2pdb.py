#!/usr/bin/python

# import packages
import sys
import os


# test input arguments
def optparser(argv=sys.argv[1:]):
	
	''' Returns "arg" with bias file name
	'''
	
	# usage: input line
	usage = '''
	Usage:

	python bias2pdb.py bias_parameter_file

		<< bias_parameter_file : input file for bias (x y z Energy Radius Type) >>
	'''
	
	# check number of arguments and return bias file name
	if len(argv) != 1:
		print 'Error: invalid number of arguments'
		print usage
		sys.exit(2)
	else:
		arg = argv[0]
	
	return arg


# requires input bias file in current directory
# get coordinates, energy, radius and type of bias sites
class Bias_Sites():
	
	''' Returns a list of dictionaries "BS" with coordinates, energy,
		radius and type for each bias site
	'''
	
	def __init__(self, bias_file):

		self.biassites = self.parse_file(bias_file)
	
	# parse bias file and get coordinates, energy, radius and type of bias sites
	def parse_file(self, bias_file):

		# open bias file
		with open(bias_file) as f:
			
			# initiate BS list
			BS = []

			# readlines puts each line of the file as an element of the list
			for i, line in enumerate(f.readlines()):
				
				# line has a list with each space or tab seaprated string of the line as element
				line = line.split()
				
				# each line needs 6 columns
				if ((len(line) != 6)):
					print 'Error: invalid number of columns in bias file'
					sys.exit(2)
				# only operate in numeric lines (skip title line)
				elif (((line[0].split(u'.')[0]).strip(u"-").isnumeric())):
					# create a dictionary (bs) with bias site data
					bs = {	'x':float(line[0]),
							'y':float(line[1]),
							'z':float(line[2]),
							'radius':float(line[4]),
							'int_type':unicode(line[5].lstrip('\n')) }
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


# get a line with PDB format
def write_PDB_format(resname,resid,name,atom_id,coordinates):
	x = coordinates[0]
	y = coordinates[1]
	z = coordinates[2]
	pdb_line = '{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}'.format('ATOM',atom_id,name,' ',resname,'X',resid,' ',x,y,z,1.00,0.00,'','')
	return(pdb_line)


###################### Main ##########################

# test input arguments and save the input bias file name
input_file = optparser()

# get coordinates, energy, radius and type of the bias sites
biases = Bias_Sites(input_file)

# open a file to store the pdb
f = open('bias_sites.pdb','w')
# write a pdb line for each bias site, depending its type
count = 1
for i in xrange(0,len(biases.biassites)):
	# get site coordinates
	coordinates = [biases.biassites[i]['x'],biases.biassites[i]['y'],biases.biassites[i]['z']]
	# write HB acceptor site
	if biases.biassites[i]['int_type'] == 'acc':
		pdb_line = write_PDB_format('ACC',count,'N',count,coordinates)
		count = count + 1
		f.write(pdb_line + '\n')
	# write HB donor site
	elif biases.biassites[i]['int_type'] == 'don':
		pdb_line = write_PDB_format('DON',count,'H',count,coordinates)
		count = count + 1
		f.write(pdb_line + '\n')
	# write aromatic site
	elif biases.biassites[i]['int_type'] == 'aro':
		pdb_line = write_PDB_format('ARO',count,'C',count,coordinates)
		count = count + 1
		f.write(pdb_line + '\n')
	# write unspecified site (specific map modification)
	elif biases.biassites[i]['int_type'] == 'map':
		pdb_line = write_PDB_format('MAP',count,'M',count,coordinates)
		count = count + 1
		f.write(pdb_line + '\n')
	# return error msg if unsupported type of bias is found
	else:
		print 'ERROR >>> Invalid type of bias found in input file.\n      >>> Type must be "acc", "don", "aro" or "map".'
		os.remove('bias_sites.pdb')
		sys.exit(2)
