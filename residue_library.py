"""
 Residue library management
 Attypes and partial charges
 data from CMIP reslib formatted file
"""

import sys

class ResiduesDataLib():
    def __init__(self, fname):
        self.residue_data = {}
        try:
            fh = open(fname, "r")
        except OSError:
            print ("#ERROR while loading library file (", fname, ")")
            sys.exit(2)
        for line in fh:
            if line[0] == '#':
                continue
            data = line.split()
            r = Residue(data)
            self.residue_data[r.id] = r
        self.nres = len(self.residue_data)

    def get_params (self, resid, atid):
        atom_id = resid + ':' + atid
        if atom_id in self.residue_data:
            return self.residue_data[atom_id]
        else:
            print ("WARNING: atom not found in library (", atom_id, ')')
            return None

class Residue():
    def __init__(self, data):
        self .id = data[0] + ':' + data[1]
        self.at_type = data[2]
        self.charge  = float(data[3])
