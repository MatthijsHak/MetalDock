import os,sys
import subprocess

def protonate_pdb(pdb_file, pH):
    subprocess.call(os.environ['PYTHON_BIN']+'/pdb2pqr30 --pdb-output pdb_prot.pdb --pH '+str(pH)+' --drop-water '+pdb_file+' pdb_prot.pdb', shell=True,  stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return



# Get clean.pdb
def clean_protein_pdb(name_protein, pdb_file, clean_pdb=True):
    if clean_pdb == True:
        subprocess.call(["sed -n '1,/\bTER\b/p;/CONECT/,$p' pdb_prot.pdb > clean_1"], shell=True)
        subprocess.call(['grep -v "HETATM" clean_1 > clean_2'], shell=True)
        subprocess.call(["mv clean_2 "+pdb_file],shell=True)

    subprocess.call(['''grep -v -e "ANISOU" -e "AALA" -e "BALA" -e "ACYS" -e "BCYS" -e "AASP" -e "BASP" -e "AGLU" -e "BGLU" -e "APHE" -e "BPHE" -e "AGLY" -e "BGLY" -e "AHIS" -e "BHIS" -e "AILE" -e "BILE" -e "ALYS" -e "BLYS" -e "ALEU" -e "BLEU" -e "AMET" -e "BMET" -e "AASN" -e "BASN" -e "APRO" -e "BPRO" -e "AGLN" -e "BGLN" -e "AARG" -e "BARG" -e "CARG" -e "ASER" -e "BSER" -e "ATHR" -e "BTHR" -e "AVAL" -e "BVAL" -e "ATRP" -e "BTRP" -e "ATYR" -e "BTYR" -e "APYL" -e "BPYL" -e "ASEC" -e "BSEC" '''+pdb_file+''' > clean_'''+name_protein+'''.pdb'''],shell=True)
    subprocess.call(["rm clean_1"], shell=True)
    return