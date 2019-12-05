#! /usr/bin/python3

"""
    Determine the list of residues whose CA atoms are closer than 20 Å
"""

import argparse

from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBParser import PDBParser


def atom_id(at):
    """ Function to build a friendly representation 
        of an atom id like ARG A23.CA 
    """
    res = at.get_parent()
    chain = res.get_parent()
    return "{} {}{}.{}".format (res.get_resname(), chain.id, res.id[1], at.id)
  
parser = argparse.ArgumentParser (
    prog='exercise2', 
    description='Getting CA pairs within max dist',
    usage='exercise2.py [options] pdb_file [> output_file]'
)

parser.add_argument(
    '--maxdist', 
    action='store', 
    dest='max_dist',
    default=20,
    type=float,
    help='Max contact distance (A)'
)

parser.add_argument('pdb_file',help='Input PDB', type=open)

args = parser.parse_args()

for k, v in vars(args).items():
    print ('{:10}:'.format(k), v)

print("PDB.filename:", args.pdb_file.name)

parser = PDBParser(PERMISSIVE=1)

print ('Parsing', args.pdb_file)

# load structure from PDB file of PDB ifle handler
st = parser.get_structure('STR', args.pdb_file.name)

# collecting CA atoms 
ca_atoms=[]

for at in st.get_atoms():
    if at.id == 'CA':
        ca_atoms.append(at)

print (len(ca_atoms), 'CA Atoms found')

# Preparing search
nbsearch = NeighborSearch(ca_atoms)

at_pairs =  nbsearch.search_all(args.max_dist)

# Output sorted by atom,serial_number, nbsearch returns ordered pairs
# Redirect the output with > output_list
for at1, at2 in sorted(at_pairs, key=lambda at_pair: at_pair[0].serial_number):
    print (atom_id(at1), ":", atom_id(at2), at1-at2)
