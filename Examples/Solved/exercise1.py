#! /usr/bin/python

# Build a command line (using argparse) with the following parameters

import argparse

parser = argparse.ArgumentParser (
    prog='exercise1', 
    description='Simple command line'
)

parser.add_argument(
    '--maxdist', 
    action='store', 
    dest='max_dist',
    help='Max contact distance'
)

parser.add_argument(
    '--atom_list',
    action='store',
    dest='at_list',
    help='Atom list, comma separated'
)
parser.add_argument(
    '-o',
    action='store',
    dest='output_path',
    help='Output file name'
)

parser.add_argument('pdb_file',help='Input PDB')

args = parser.parse_args()

for k, v in vars(args).items():
    print ('{:10}:'.format(k), v)
