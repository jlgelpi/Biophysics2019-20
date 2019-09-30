#! /usr/bin/python

# Determine the list of residues whose CA atoms are closer than 20 Å

import argparse

from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBParser import PDBParser

    
def atom_id(at):
    res = at.get_parent()
    chain = res.get_parent()
    return "{} {}{}.{}".format (res.get_resname(), chain.id, res.id[1], at.id)
  

parser = argparse.ArgumentParser (
    prog='exercise2', 
    description='CA pairs within max dist'
)

parser.add_argument(
    '--maxdist', 
    action='store', 
    dest='max_dist',
    default=20,
    type=float,
    help='Max contact distance'
)

parser.add_argument('pdb_file',help='Input PDB')

args = parser.parse_args()

for k, v in vars(args).items():
    print ('{:10}:'.format(k), v)

parser = PDBParser(PERMISSIVE=1)

print ('Parsing', args.pdb_file)
# load structure from PDB file

st = parser.get_structure('1UBQ', args.pdb_file)

ca_atoms=[]

for at in st.get_atoms():
    if at.id == 'CA':
        ca_atoms.append(at)

print (len(ca_atoms), 'CA Atoms found')

# Preparing search
nbsearch = NeighborSearch(ca_atoms)

at_pairs =  nbsearch.search_all(args.max_dist)



for at1, at2 in sorted(at_pairs, key=lambda at_pair: at_pair[0].serial_number):
    print (atom_id(at1), ":", atom_id(at2), at1-at2)
    

