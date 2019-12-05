#! /usr/bin/python3
"""
    Print all distances between the atoms of two given residues
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
    '--res1', 
    action='store', 
    dest='res1',
    help='Residue 1 (as A23)',
    required=True
)

parser.add_argument(
    '--res2', 
    action='store', 
    dest='res2',
    help='Residue 2 (as B25)',
    required=True
)


parser.add_argument('pdb_file',help='Input PDB', type=open)

args = parser.parse_args()

print ("Arguments")

for k, v in vars(args).items():
    print ('{:10}:'.format(k), v)

# Defining chains and residue numbers

chain_id1 = args.res1[0:1]
res_num1 = int(args.res1[1:])

chain_id2 = args.res2[0:1]
res_num2 = int(args.res2[1:])

print()
print ("PDB file:", args.pdb_file.name)
print ("Selected Residue 1 Chain: {}, Residue number: {}".format(chain_id1,res_num1))
print ("Selected Residue 2 Chain: {}, Residue number: {}".format(chain_id2,res_num2))

# Check whether the input is complete sys.exit print on the std.err and exits
if not chain_id1 or not res_num1:
    sys.exit("ERROR: unknown either chain id or residue 1 number")

if not chain_id2 or not res_num2:
    sys.exit("ERROR: unknown either chain id or residue 2 number")

parser = PDBParser(PERMISSIVE=1)

print()
print ('Parsing', args.pdb_file.name)

# load structure from PDB file of PDB ifle handler
st = parser.get_structure('STR', args.pdb_file)

# Checking residues exist and are different
if chain_id1 not in st[0] or res_num1 not in st[0][chain_id1]:
    sys.exit("ERROR: non existing chain or residue")
if chain_id2 not in st[0] or res_num2 not in st[0][chain_id2]:
    sys.exit("ERROR: non existing chain or residue")
if (chain_id1 == chain_id2) and (res_num1 == res_num2):
    sys.exit("ERROR: identical residues")
print()
print ("Atom list")

for at1 in st[0][chain_id1][res_num1].get_atoms():
    for at2 in st[0][chain_id2][res_num2].get_atoms():
        print ('{:11}-{:11}: {:8.3f} (A)'.format(atom_id(at1), atom_id(at2), at2-at1))


