#! /usr/bin/python3

"""
    Build a command line (using argparse) with the following parameters
      a. PDB file name (required)
      b. Max contact distance (float, mandatory)
      c. Atom list (string)
      d. Output file name (string)
    
"""
    

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
# 
# optional argument
# identical funcionality can be obtained by output redirections as in
# exercise1 pdb_file > output
#
parser.add_argument(
    '-o',
    action='store',
    type=argparse.FileType('w'),
    dest='output_path',
    help='Output file name'
)

parser.add_argument('pdb_file',type=open, help='Input PDB')

args = parser.parse_args()

for k, v in vars(args).items():
    print ('{:10}:'.format(k), v)
