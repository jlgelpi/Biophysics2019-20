#! /usr/bin/python3
"""
    Generate a list of all CA atoms of given residue type with coordinates
"""

import argparse

from Bio.PDB.PDBParser import PDBParser

def atom_id(at):
    """ Function to build a friendly representation 
        of an atom id like ARG A23.CA 
    """
    res = at.get_parent()
    chain = res.get_parent()
    return "{} {}{}.{}".format (res.get_resname(), chain.id, res.id[1], at.id)

parser = argparse.ArgumentParser (
    prog='exercise5', 
    description='Generate a list of all CA atoms of given residue type with coordinates',
    usage='exercise5.py [options] pdb_file [> output_file]'
)


parser.add_argument(
    '--restype', 
    action='store', 
    dest='res_type',
    help='Residue type',
    required=True
)

parser.add_argument('pdb_file',help='Input PDB', type=open)

args = parser.parse_args()

print ("Arguments")

for k, v in vars(args).items():
    print ('{:10}:'.format(k), v)

print()
print ("PDB file:", args.pdb_file.name)

parser = PDBParser(PERMISSIVE=1)

print()
print ('Parsing', args.pdb_file.name)

# load structure from PDB file of PDB ifle handler
st = parser.get_structure('STR', args.pdb_file)

print()
print ("Atom list")

for at in st[0].get_atoms():
    if at.get_parent().get_resname() == args.res_type:
        print ('{:11} {}'.format(atom_id(at), at.get_coord()))
