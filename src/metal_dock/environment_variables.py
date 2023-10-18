#!/usr/bin/env python
import os,sys,subprocess

def find_command_path(command):
    try:
        if sys.platform.startswith('win'):
            # On Windows, use the 'where' command
            result = subprocess.check_output(['where', command], universal_newlines=True)
        else:
            # On Unix-based systems, use the 'which' command
            result = subprocess.check_output(['which', command], universal_newlines=True)

        # Remove any leading/trailing whitespace and return the path
        return result.strip()
    except subprocess.CalledProcessError:
        print(f"Error: The command '{command}' was not found in your system's PATH.")
        sys.exit(1)


ROOT_DIR = os.path.dirname(os.path.abspath(__file__))[:-11] # Project Root
os.environ['ROOT_DIR']=ROOT_DIR

#os.environ['LIB_DIR']=sys.executable.split( os.path.sep+'bin'+os.path.sep,1)[0]

os.environ['OBABEL']=find_command_path('obabel')
#os.path.join(os.environ['LIB_DIR'],'bin','obabel')
os.environ['PDB2PQR']=find_command_path('pdb2pqr30')
#os.path.join(os.environ['LIB_DIR'],'bin','pdb2pqr30')
os.environ['MGLTOOLS']=os.path.join(os.environ['ROOT_DIR'],'external','AutoDockTools','Utilities24')
os.environ['PYTHON_2']=find_command_path('python2.7')
#os.path.join(os.environ['LIB_DIR'],'bin','python2.7')
os.environ['PYTHON_3']=find_command_path('python3')
#os.path.join(os.environ['LIB_DIR'],'bin','python3')


#obabel_path = find_command_path('obabel')
#print(obabel_path)
#
#pdb2pqr30 = find_command_path('pdb2pqr30')
#print(pdb2pqr30)
#
#python2_path = find_command_path('python2.7')
#print(python2_path)
#
#python3_path = find_command_path('python3')
#print(python3_path)
#
#print(os.environ['OBABEL'])
#print(os.environ['PDB2PQR'])
#print(os.environ['PYTHON_2'])
#print(os.environ['PYTHON_3'])
#sys.exit()
