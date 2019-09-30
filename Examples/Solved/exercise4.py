#! /usr/bin/python

# Determine all possible hydrogen bonds (Polar atoms at less than 3.5 Å).

import argparse

from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBParser import PDBParser

HBLNK = 3.5

POLAR = ['N','O','S']

def atom_id(at):
    res = at.get_parent()
    chain = res.get_parent()
    return "{} {}{}.{}".format (res.get_resname(), chain.id, res.id[1], at.id)
  

parser = argparse.ArgumentParser (
    prog='exercise4', 
    description='HBonds list'
)


parser.add_argument('pdb_file',help='Input PDB')

parser.add_argument('--hetatm', action='store_true', help='Include non protein atoms')

args = parser.parse_args()

for k, v in vars(args).items():
    print ('{:10}:'.format(k), v)

parser = PDBParser(PERMISSIVE=1)

print ('Parsing', args.pdb_file)
# load structure from PDB file

st = parser.get_structure('1UBQ', args.pdb_file)

polar_ats = []

for at in st.get_atoms():
    res = at.get_parent()
    if res.id[0].startswith('H_') or res.id[0].startswith('W') and not args.hetatm:
        continue
    if at.element in POLAR:
        polar_ats.append(at)

print (len(polar_ats), 'POLAR Atoms found')

# Preparing search
nbsearch = NeighborSearch(polar_ats)

at_pairs =  nbsearch.search_all(HBLNK)

hbs = {}

for at1, at2 in sorted(at_pairs, key=lambda at_pair: at_pair[0].serial_number):
    res1 = at1.get_parent()
    res2 = at2.get_parent()
    # remove atom pairs from the same residue and next in sequence
    if res2.id[1] - res1.id[1] < 2:
        continue
    # using the contact with the shortest distance between residues
    dist = at1 - at2
    if res1 not in hbs:
        hbs[res1] = {}
    if res2 not in hbs[res1]:
        hbs[res1][res2]= (at1, at2, dist)
    else:
        if dist < hbs[res1][res2][2]:
            hbs[res1][res2] = (at1, at2, dist)

for res1 in hbs:
    for res2 in hbs[res1]:
        print (
            atom_id(hbs[res1][res2][0]), ":", atom_id(hbs[res1][res2][1]), 
            hbs[res1][res2][2]
        )
