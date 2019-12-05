#! /usr/bin/python3

""" Obtain a PDB file including only a list of chains from the original structure """

import argparse

from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBParser import PDBParser

parser = argparse.ArgumentParser (
    prog='exercise8', 
    description='Obtain a PDB file including only a list of chains from the original structure',    
)

parser.add_argument(
    '--chains', 
    action='store', 
    dest='chains',
    help='Accepted chain list (comma sep)',
    required=True
)

parser.add_argument('pdb_file',help='Input PDB', type=open)
parser.add_argument('output_pdb_file',help='Output PDB', type=argparse.FileType('w'))

args = parser.parse_args()

for k, v in vars(args).items():
    print ('{:10}:'.format(k), v)

print("PDB.filename:", args.pdb_file.name)


parser = PDBParser()

st = parser.get_structure('STR', args.pdb_file)

print ('{} chains found: {}'.format(len(st[0]), ','.join([ch.id for ch in st[0]])))

chains_ok = args.chains.split(',')
if not chains_ok:
    sys.exit('ERROR: No chains selected')

deleted = 0
for chn in st[0]:
    if chn.id not in chains_ok:
        st[0].detach_child(chn.id)
        deleted += 1 
        print ("Removing chain", chn.id)

if not deleted:
    print ("WARNING: no chains deleted")
    
pdbio = PDBIO()

pdbio.set_structure(st)
pdbio.save(args.output_pdb_file)

