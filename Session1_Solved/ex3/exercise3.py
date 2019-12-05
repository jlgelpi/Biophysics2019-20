#! /usr/bin/python3
"""
    Generate a list of all atoms for a given residue number
"""

import argparse

import sys

from Bio.PDB.PDBParser import PDBParser

def atom_id(at):
    """ Function to build a friendly representation 
        of an atom id like ARG A23.CA 
    """
    res = at.get_parent()
    chain = res.get_parent()
    return "{} {}{}.{}".format (res.get_resname(), chain.id, res.id[1], at.id)

parser = argparse.ArgumentParser (
    prog='exercise3', 
    description='Generate a list of all atoms for a given residue number',
    usage='exercise3.py [options] pdb_file [> output_file]'
)

# Assumed normal residue numbering, may fail in odd formatted PDB files

parser.add_argument(
    '--resnum', 
    action='store', 
    dest='res_num',
    type=int,
    help='Residue Number'
)

# Optional: For multichain proteint we may want to decide on the chain as well

parser.add_argument(
    '--chain', 
    action='store', 
    dest='chain_id',
    default='A',
    help='Chain id (default A)'
)

# Optional: Combine chain and number in single entry

parser.add_argument(
    '--residue', 
    action='store', 
    dest='res_id',
    help='Residue id (like A23) (optional, overrides --chain and --resnum)'
)

parser.add_argument('pdb_file',help='Input PDB', type=open)

args = parser.parse_args()

print ("Arguments")

for k, v in vars(args).items():
    print ('{:10}:'.format(k), v)

# Defining chain and residue number if --residue has been used

if args.res_id:
    chain_id = args.res_id[0:1]
    res_num = int(args.res_id[1:])
else:
    chain_id = args.chain_id
    res_num =args.res_num
print()
print ("PDB file:", args.pdb_file.name)
print ("Selected Chain: {}, Residue number: {}".format(chain_id,res_num))
# Check whether the input is complete sys.exit print on the std.err and exits
if not chain_id or not res_num:
    sys.exit("ERROR: missing either chain id or residue number")

parser = PDBParser(PERMISSIVE=1)

print()
print ('Parsing', args.pdb_file.name)

# load structure from PDB file of PDB ifle handler
st = parser.get_structure('STR', args.pdb_file)

# Checking residue exists
if chain_id not in st[0] or res_num not in st[0][chain_id]:
    sys.exit("ERROR: non existing chain or residue")

print()
print ("Atom list")

for at in st[0][chain_id][res_num].get_atoms():
    print (atom_id(at))
