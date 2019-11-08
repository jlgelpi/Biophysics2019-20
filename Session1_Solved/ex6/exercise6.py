#! /usr/bin/python3

"""
    Generate a list of backbone connectivity.
"""

import argparse

from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBParser import PDBParser

PEP_BOND_ATS = ['N', 'C']
COVLNK = 2.5 

def atom_id(at):
    """ Function to build a friendly representation 
        of an atom id like ARG A23.CA 
    """
    res = at.get_parent()
    chain = res.get_parent()
    return "{} {}{}.{}".format (res.get_resname(), chain.id, res.id[1], at.id)
  
parser = argparse.ArgumentParser (
    prog='exercise6', 
    description='Generate a list of backbone connectivity',
    usage='exercise6.py [options] pdb_file [> output_file]'
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

# collecting atom candidates
bck_atoms=[]

for at in st.get_atoms():
    if at.id in PEP_BOND_ATS :
        bck_atoms.append(at)

print (len(bck_atoms), 'candidate atoms found')

# Preparing search
nbsearch = NeighborSearch(bck_atoms)

at_pairs =  nbsearch.search_all(COVLNK)


# Output sorted by atom,serial_number, nbsearch returns ordered pairs
# Redirect the output with > output_list
for at1, at2 in sorted(at_pairs, key=lambda at_pair: at_pair[0].serial_number):
    # Discard same residue
    if at1.get_parent() == at2.get_parent():
        continue
    print ('{:11} : {:11} {:8.3f}'.format(atom_id(at1),atom_id(at2), at1-at2))
