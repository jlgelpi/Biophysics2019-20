#! /usr/bin/python3

"""
    Initial setup a structure for Energy evaluation
"""
import argparse
import sys
import os

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NACCESS import NACCESS_atomic
from Bio.PDB.PDBIO import PDBIO, Select

from residue_library import ResiduesDataLib
from forcefield import VdwParamset
import energies as en

NACCESS_BINARY = '/home/gelpi/DEVEL/BioPhysics/2019-20/code2019/soft/NACCESS/naccess'

parse_cmd = argparse.ArgumentParser (
    prog='structure_setup', 
    description='basic structure setup'
)

parse_cmd.add_argument(
    '--rlib',
    action='store', 
    dest='reslib_file',
    default='data/aaLib.lib', 
    help='Residue Library'
)
parse_cmd.add_argument(
    '--vdw',
    action='store', 
    dest='vdwprm_file',
    default='data/vdwprm', 
    help='Vdw parameters'
)

parse_cmd.add_argument(
    '--dist',
    action='store',
    dest='cutoff_dist',
    default=8.0,
    type=float,
    help='Cutoff distance for determining the interface'
)
parse_cmd.add_argument('pdb_file',help='Input PDB', type=open)

args = parse_cmd.parse_args()

print("PDB.filename:", args.pdb_file.name)
print("Residue Lib.:", args.reslib_file)
print("PDB.filename:", args.vdwprm_file)
print("Distance:", args.cutoff_dist)

# Loading Libraries
# loading residue library from data/aaLib.lib
residue_library = ResiduesDataLib(args.reslib_file)

# loading VdW parameters
ff_params = VdwParamset(args.vdwprm_file)

parser = PDBParser(PERMISSIVE=1)
print ('Parsing', args.pdb_file)
# load structure from PDB file of PDB ifle handler
st = parser.get_structure('STR', args.pdb_file.name)

# assign data types, and charges from libraries
# We will use the xtra attribute in Bio.PDB.Atom to hold the new data
# Possible errors on N-term and C-Term atoms
# Possible errors on HIS alternative forms

en.add_atom_parameters(st, residue_library, ff_params)

# Calculating surfaces
# The specific PATH to naccess script (in soft) is needed
# ASA goes to .xtra field directly

srf = NACCESS_atomic(st[0],naccess_binary=NACCESS_BINARY)

# Prepare surfaces for the separate chains
# Alternatively the twp PDB files can be prepared outside and parsed here
io = PDBIO()
st_chains={}

class SelectChain(Select):
    def __init__(self, chid):
        self.id = chid
    
    def accept_chain(self, chain):
        if chain.id == self.id:
            return 1
        else:
            return 0
        
for ch in st[0]:
    io.set_structure(st)
    io.save ('tmp.pdb',SelectChain(ch.id))
    st_chains[ch.id] = parser.get_structure('stA','tmp.pdb')
    en.add_atom_parameters(st_chains[ch.id], residue_library, ff_params)
    srfA = NACCESS_atomic(st_chains[ch.id][0], naccess_binary=NACCESS_BINARY)
os.remove('tmp.pdb')

## Interface residues

interface = en.get_interface(st, args.cutoff_dist)

## Energies
elec = {}
elec_ala = {}

vdw = {}
vdw_ala = {}

solvAB = {}
solvAB_ala = {}

solvA = {}
solvA_ala = {}

chids = []
for ch in st[0]:
    chids.append(ch.id)
    
totalIntElec = 0.
totalIntVdw = 0.
totalSolv = 0.
totalSolvMon = {'A':0.,'B':0.}
total = 0.

for ch in st[0]:
    for res in interface[ch.id]:
        elec[res], elec_ala[res], vdw[res], vdw_ala[res] = en.calc_int_energies(st[0],res)
        solvAB[res], solvAB_ala[res] = en.calc_solvation(st[0],res)
        solvA[res], solvA_ala[res] = en.calc_solvation(
            st_chains[ch.id], 
            st_chains[ch.id][0][ch.id][res.id[1]]
        )
        totalIntElec += elec[res]
        totalIntVdw += vdw[res]
        totalSolv += solvAB[res]
        totalSolvMon[ch.id] += solvA[res]
        total += elec[res] + vdw[res] + solvAB[res] - solvA[res]
print('{:20}: {:11.4f}'.format('Total Elec Int:',totalIntElec))
print('{:20}: {:11.4f}'.format('Total Vdw Int:',totalIntVdw))
print('{:20}: {:11.4f}'.format('Total Solv AB:',totalSolv))
print('{:20}: {:11.4f}'.format('Total Solv A:',totalSolvMon['A']))
print('{:20}: {:11.4f}'.format('Total Solv B:',totalSolvMon['B']))
print('{:20}: {:11.4f}'.format('DGintAB-A-B:',total))


print("Ala Scanning")
for ch in st[0]:
    for res in interface[ch.id]:
        print(
            '{:11} {:11.4f}{:11.4f}{:11.4f}{:11.4f}{:11.4f}'.format(
                en.residue_id(res), 
                totalIntElec - elec[res] + elec_ala[res], 
                totalIntVdw - vdw[res] + vdw_ala[res],
                totalSolv - solvAB[res] + solvAB_ala[res],
                totalSolvMon[ch.id] - solvA[res] + solvA_ala[res],
                - elec[res] + elec_ala[res] - vdw[res] + vdw_ala[res] -solvAB[res] +\
                    solvAB_ala[res] -solvA[res] + solvA_ala[res]
            )
        )
            