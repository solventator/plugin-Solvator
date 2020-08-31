from __future__ import division

from olexFunctions import OlexFunctions
OV = OlexFunctions()

import numpy as np

class Solvator_atom():
    def __init__(self, symbol, label, sfac, coords, fvar = 1, occupancy = 1, uiso = 0.03, u  = None, coord_type = 'frac'):
        
        self.symbol = symbol
        self.label =  label
        self.sfac = sfac
        self.coords = coords
        self.fvar = fvar
        self.occupancy = occupancy
        self.uiso = uiso
        self.u = u
        self.peak_height = None

        self.coord_type = coord_type
        self.olex_id = None
        self.solavtor_id = None
        name = self.symbol + str(self.label)
        self.name = name
        
    def set_occupancy(self,new_occupancy):
        self.occupancy = new_occupancy
        
    def change_occupancy(self,change_in_occupancy):
        self.occupancy += change_in_occupancy
        
    def make_anisotropic(self):
        if not u:
            u = [self.uiso, self.uiso, self.uiso, self.uiso, self.uiso, self.uiso]
        else:
            print "Atom is already anisotropic: %s%s" % (atom.symbol, atom.label)
    
    def make_isotropic(self):
        self.uiso = (self.u[0] + self.u[1] + self.u[2])/3
        self.u = None
    
    def shift_slightly(self):
        self.coords[0][0]+=0.0001
        self.coords[1][0]+=0.0001
        self.coords[2][0]+=0.0001
        