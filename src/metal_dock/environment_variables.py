#!/usr/bin/env python
import os,sys

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))[:-11] # Project Root
os.environ['ROOT_DIR']=ROOT_DIR

os.environ['LIB_DIR']=sys.executable.split('/bin/',1)[0]

os.environ['OBABEL']=os.environ['LIB_DIR']+'/bin/obabel'
os.environ['PDB2PQR']=os.environ['LIB_DIR']+'/bin/pdb2pqr30'
os.environ['MGLTOOLS']=os.environ['LIB_DIR']+'/MGLToolsPckgs/AutoDockTools/Utilities24'
os.environ['PYTHON_2']=os.environ['LIB_DIR']+'/bin/python2.7'
os.environ['PYTHON_3']=os.environ['LIB_DIR']+'/bin/python3'