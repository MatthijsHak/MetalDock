import os
import shutil
import sys
import subprocess as sp
from src.metal_dock.logger import MetalDockLogger

class Protein:
    def __init__(self, par):
        self.par = par
        self.logger = MetalDockLogger()
        self.protonate_pdb()
        self.clean_protein_pdb()

    def protonate_pdb(self):
        """
        Protonate the protein.
        """
        if not (self.par.output_dir / 'file_prep' / f'{self.par.name_protein}.pdb').exists():
            input_pdb = self.par.pdb_file
            output_pdb_protonated = self.par.output_dir / 'file_prep' / f'{self.par.name_protein}.pdb'

            if self.par.clean_pdb == True:
                sp.call(os.environ['PDB2PQR']+f' --noopt --pdb-output {output_pdb_protonated} --with-ph {str(self.par.pH)} --drop-water {input_pdb} {output_pdb_protonated}', shell=True,  stdout=sp.PIPE, stderr=sp.PIPE)
            else:
                sp.call(os.environ['PDB2PQR']+f' --noopt --pdb-output {output_pdb_protonated} --with-ph {str(self.par.pH)} {input_pdb} {output_pdb_protonated}', shell=True,  stdout=sp.PIPE, stderr=sp.PIPE)

    def clean_protein_pdb(self):
        """
        Remove cofactors or ligands from the pdb file 
        """
        pdb_input_path = self.par.output_dir / 'file_prep' / f'{self.par.name_protein}.pdb'
        pdb_output_path = self.par.output_dir / 'file_prep' / f'clean_{self.par.name_protein}.pdb'

        if self.par.clean_pdb == True:
            with open(pdb_input_path, 'r') as fin:
                with open(pdb_output_path, 'w') as fout:
                    for line in fin:
                        if 'HETATM'  not in line:
                            fout.write(line)
        else:
            shutil.move('pdb_prot.pdb', pdb_output_path)

        self.pdb_file = pdb_output_path

    def create_pdbqt_file(self):
        """
        Create the pdbqt file for the protein.
        """
        if not (self.par.output_dir / 'file_prep' / f'clean_{self.par.name_protein}.pdbqt').exists():
            pdb_path = self.par.output_dir / 'file_prep' / f'clean_{self.par.name_protein}.pdb'
            pdbqt_path = self.par.output_dir / 'file_prep' / f'clean_{self.par.name_protein}.pdbqt'
            prepare_gpf4 = os.path.join(os.environ['MGLTOOLS'], 'prepare_receptor4.py')
            command = os.environ['PYTHON_3']+f' {prepare_gpf4} -U nphs -A None -r {pdb_path} -o {pdbqt_path}'

            process = sp.Popen(
                command,
                shell=True,
                stdout=sp.PIPE,
                stderr=sp.PIPE,
                universal_newlines=True
            )
            stdout, stderr = process.communicate()

            if process.returncode != 0:
                out_file = self.par.output_dir / 'docking' / 'prepare_receptor.out'
                with open(out_file, 'w') as fout:
                    fout.write(stdout)
                    fout.write(stderr)
                self.logger.error('ERROR DURING PREPARATION OF RECEPTOR, SEE /output/docking/prepare_receptor.out FOR DETAILS')
                self.logger.error('IF PDB FILE IS NOT WRITTEN IN CORRECT PDB FORMAT, PLEASE EDIT MANUALLY')
                sys.exit()
        else:
            self.logger.info('RECEPTOR PDBQT FILE ALREADY EXISTS\n')
