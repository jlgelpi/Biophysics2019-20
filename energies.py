''' Function to evaluate energies'''

from Bio.PDB.NeighborSearch import NeighborSearch
import math

#Possible Atom names that correspond to Ala atoms"
ala_atoms = {'N', 'H', 'CA', 'HA', 'C', 'O', 'CB', 'HB', 'HB1', 'HB2', 'HB3', 'HA1', 'HA2', 'HA3'}

def residue_id(res):
    '''Returns readable residue id'''
    return '{} {}{}'.format(res.get_resname(), res.get_parent().id, res.id[1])

def atom_id(at):
    '''Returns readable atom id'''
    return '{}.{}'.format(residue_id(at.get_parent()), at.id)

def add_atom_parameters(st, res_lib, ff_params):
    ''' Adds parameters from libraries to atom .xtra field
        For not recognized atoms, issues a warning and put default parameters
    '''
    for at in st.get_atoms():
        resname = at.get_parent().get_resname()
        params = res_lib.get_params(resname, at.id)
        if not params:
            print("WARNING: residue/atom pair not in library ("+atom_id(at) + ')')
            at.xtra['atom_type'] = at.element
            at.xtra['charge'] = 0
        else:
            at.xtra['atom_type'] = params.at_type
            at.xtra['charge'] = params.charge
        at.xtra['vdw'] = ff_params.at_types[at.xtra['atom_type']]

def get_interface(st, dist):
    ''' Detects interface residues within a distance(dist)
        Assumes two chains, i.e. a unique interface set per chain.
    '''
    select_ats = []
    for at in st.get_atoms():
        # Skip Hydrogens to reduce time
        if at.element != 'H':
            select_ats.append(at)
    nbsearch = NeighborSearch(select_ats)
    interface = {}
    # Sets are more efficient than lists. Use sets when order is not relevant
    for ch in st[0]:
        interface[ch.id] = set()

    for at1, at2 in nbsearch.search_all(dist):
        #Only different chains
        res1 = at1.get_parent()
        ch1 = res1.get_parent()
        res2 = at2.get_parent()
        ch2 = res2.get_parent()
        if ch1 != ch2:
            interface[ch1.id].add(res1)
            interface[ch2.id].add(res2)
    return interface

def calc_int_energies(st, res):
    '''Returns interaction energies (residue against other chains)
        for all atoms and for Ala atoms
    '''
    elec = 0.
    elec_ala = 0.
    vdw = 0.
    vdw_ala = 0.

    for at1 in res.get_atoms():
        for at2 in st.get_atoms():
        # skip same chain atom pairs
            if at2.get_parent().get_parent() != res.get_parent():
                r = at1 - at2
                e = elec_int(at1, at2, r)
                elec += e
                if at1.id in ala_atoms: #GLY are included implicitly
                    elec_ala += e
                e = vdw_int(at1, at2, r)
                vdw += e
                if at1.id in ala_atoms: #GLY are included implicitly
                    vdw_ala += e
    return elec, elec_ala, vdw, vdw_ala

def MH_diel(r):
    '''Mehler-Solmajer dielectric'''
    return 86.9525 / (1 - 7.7839 * math.exp(-0.3153 * r)) - 8.5525

def elec_int(at1, at2, r):
    '''Electrostatic interaction energy between two atoms at r distance'''
    return 332.16 * at1.xtra['charge'] * at2.xtra['charge'] / MH_diel(r) / r

def vdw_int(at1, at2, r):
    '''Vdw interaction energy between two atoms'''
    eps12 = math.sqrt(at1.xtra['vdw'].eps * at2.xtra['vdw'].eps)
    sig12_2 = at1.xtra['vdw'].sig * at2.xtra['vdw'].sig
    return 4 * eps12 * (sig12_2**6/r**12 - sig12_2**3/r**6)

def calc_solvation(st, res):
    '''Solvation energy based on ASA'''
    solv = 0.
    solv_ala = 0.
    for at in res.get_atoms():
        if 'EXP_NACCESS' not in at.xtra:
            continue
        s = float(at.xtra['EXP_NACCESS'])* at.xtra['vdw'].fsrf
        solv += s
        if at.id in ala_atoms:
            solv_ala += s
    return solv, solv_ala
