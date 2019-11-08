#! /usr/bin/python3

"""
    Determine all possible hydrogen bonds (Polar atoms at less than 3.5 Å).
"""

import argparse

from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBParser import PDBParser

POLAR_ELEMENTS = ['N', 'O', 'S']

# Optional test andclassify donor-acceptor pairs
HB_DONORS = ['ND1','ND2','NE1','NE2','NH1', 'NH2', 'NZ','N', 'OH', 'SG']
HB_ACCEPTORS = ['O', 'OH', 'SG', 'SD','OD1','OD2','OE1','OE2']

def atom_id(at):
    """ Function to build a friendly representation 
        of an atom id like ARG A23.CA 
    """
    res = at.get_parent()
    chain = res.get_parent()
    return "{} {}{}.{}".format (res.get_resname(), chain.id, res.id[1], at.id)
  
def check_hb (at1, at2):
    """ Check whether the HB is a valid one and classifies atoms """
    ad=False
    da=False
    if at1.id in HB_DONORS and at2.id in HB_ACCEPTORS:
        da = True
    if at1.id in HB_ACCEPTORS and at2.id in HB_DONORS:
        ad = True
    return da,ad
    
parser = argparse.ArgumentParser (
    prog='exercise4', 
    description='Determine all possible hydrogen bonds',
    usage='exercise4.py [options] pdb_file [> output_file]'
)

parser.add_argument(
    '--check', 
    action='store_true', 
    help='Check HB pairs and classify'
)

parser.add_argument(
    '--hb_max_dist', 
    action='store', 
    type=float,
    default=3.5,
    help='Define HB max dist (default: 3.5A)'
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

# collecting Polar atoms 
polar_atoms=[]

for at in st.get_atoms():
    if at.element in POLAR_ELEMENTS:
        polar_atoms.append(at)

print (len(polar_atoms), 'polar Atoms found')

# Preparing search
nbsearch = NeighborSearch(polar_atoms)

at_pairs =  nbsearch.search_all(args.hb_max_dist)


# Output sorted by atom,serial_number, nbsearch returns ordered pairs
# Redirect the output with > output_list
for at1, at2 in sorted(at_pairs, key=lambda at_pair: at_pair[0].serial_number):
    ad = True
    da = True
    if args.check:
        da,ad = check_hb (at1, at2)
    if da and ad: # classification not requested or not predictable
        print ('{:11}  : {:11}  {:8.3f}'.format(atom_id(at1),atom_id(at2), at1-at2))
    elif da: # Donor - Acceptor
        print ('{:11}D : {:11}A {:8.3f}'.format(atom_id(at1),atom_id(at2), at1-at2))
    elif ad : # Acceptor - Donor
        print ('{:11}A : {:11}D {:8.3f}'.format(atom_id(at1),atom_id(at2), at1-at2))
    # both false imply a discarded HB

# Simpler code without the --check option
#for at1, at2 in sorted(at_pairs, key=lambda at_pair: at_pair[0].serial_number):
#    print ('{:11}  : {:11}  {:8.3f}'.format(atom_id(at1),atom_id(at2), at1-at2))
