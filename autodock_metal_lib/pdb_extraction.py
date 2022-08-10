import os
import input_variables as iv

# Get clean.pdb
def clean_protein_pdb(pdb_file):
    os.system('''sed -n '1,/\bTER\b/p;/CONECT/,$p' '''+pdb_file+''' > clean_1''')
    os.system("grep -v HETATM clean_1 > clean_2")
    os.system('''grep -v -e "ANISOU" -e "AALA" -e "BALA" -e "ACYS" -e "BCYS" -e "AASP" -e "BASP" -e "AGLU" -e "BGLU" -e "APHE" -e "BPHE" -e "AGLY" -e "BGLY" -e "AHIS" -e "BHIS" -e "AILE" -e "BILE" -e "ALYS" -e "BLYS" -e "ALEU" -e "BLEU" -e "AMET" -e "BMET" -e "AASN" -e "BASN" -e "APRO" -e "BPRO" -e "AGLN" -e "BGLN" -e "AARG" -e "BARG" -e "CARG" -e "ASER" -e "BSER" -e "ATHR" -e "BTHR" -e "AVAL" -e "BVAL" -e "ATRP" -e "BTRP" -e "ATYR" -e "BTYR" -e "APYL" -e "BPYL" -e "ASEC" -e "BSEC" clean_2 > clean_'''+iv.var.name_protein+'''.pdb''')
    os.system("rm clean_1 clean_2")

def protonate_pdb(pdb_file):
    os.system(os.environ['PYTHON_BIN']+'/pdb2pqr30 --pdb-output pdb_prot.pdb --pH 7.4 --drop-water '+pdb_file+' pdb_prot.pdb')


def get_ref_pdb(pdb_file):
    os.system("grep -v HOH "+pdb_file+" > protein_noH2O.pdb")
    os.system("sed -n '/^HETATM/p' protein_noH2O.pdb > LIG_1.pdb")
    os.system("grep -v SO4 LIG_1.pdb > LIG_2.pdb")

    fin = open('LIG_2.pdb')
    line = [line.split() for line in fin]

    TAG = ''

    for count, i in enumerate(range(0,len(line))):
        if ""+iv.var.metal_symbol+"" in line[i][-1] or ""+iv.var.metal_cap+"" in line[i][-1]:
            TAG = line[i][-7]
            break

    os.system(r'''awk '{if ($4 == "'''+TAG+'''" || $5 == "'''+TAG+'''" || $6 == "'''+TAG+'''") { print }}' LIG_2.pdb > ref.pdb''')
    os.system("rm protein_noH2O.pdb LIG_1.pdb LIG_2.pdb")
