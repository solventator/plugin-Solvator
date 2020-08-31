from __future__ import division

from olexFunctions import OlexFunctions
OV = OlexFunctions()

import numpy as np
import guest
import random

class Q_cluster():
    
    def __init__(self, peak_list, electron_count, centroid):
        
        try:
            self.peak_list = peak_list
            self.electron_count = electron_count
            self.centroid = centroid
            self.model = []
            self.solution = None
            self.starting_positions = []
            self.multiplicity = 1
        except:
            print "Could not create new cluster"
        
    def add_molecule(self, guest_name, centroid, p_path, SFAC_BY_SYMBOL, CART_TO_FRAC):

        print guest_name, centroid, p_path, SFAC_BY_SYMBOL
        molecule = guest.Guest(guest_name, 0.5, p_path, SFAC_BY_SYMBOL, CART_TO_FRAC, centroid) # always add things at 0.5 occupancy
        self.solution.append(molecule)    

    def delete_molecule(self):

        if len(self.solution) > 1:
            random_number = random.randrange(len(self.solution))
            del self.solution[random_number]
        else:
            print "Cannot delete only molecule in solution"
            print "Number of molecules in temp void solution: %s" % len(self.solution)

        return self.solution