#! /usr/bin/python

# Determine the list of residues whose CA atoms are closer than 20 Å

import argparse
import re

from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBParser import PDBParser

#function for nice atom printing    
def atom_id(at):
    res = at.get_parent()
    chain = res.get_parent()
    return "{} {}{}.{}".format (res.get_resname(), chain.id, res.id[1], at.id)
  

parser = argparse.ArgumentParser (
    prog='exercise3', 
    description='List of atoms per residue'
)

parser.add_argument(
    '--res', 
    action='store', 
    dest='res',
    help='Residue Number as A220',
    required=True
)

parser.add_argument('pdb_file',help='Input PDB')

args = parser.parse_args()

for k, v in vars(args).items():
    print ('{:10}:'.format(k), v)

parser = PDBParser(PERMISSIVE=1)

print ('Parsing', args.pdb_file)
# load structure from PDB file

st = parser.get_structure('1UBQ', args.pdb_file)

# parsing input into chain + res number
m = re.match('([A-z]*)([0-9]*)', args.res)
chain_id = m[1]
res_num = int(m[2])

# Getting atoms
if chain_id in st[0]:
    if res_num in st[0][chain_id]:
        residue = st[0][chain_id][res_num]
        for at in residue.get_atoms():
            print(atom_id(at), at.get_coord())
    