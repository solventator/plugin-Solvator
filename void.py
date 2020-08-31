from __future__ import division

from olexFunctions import OlexFunctions
OV = OlexFunctions()

import numpy as np
import guest
import random

class Void():
    def __init__(self, number, volume, electrons):
        
        self.number = number
        self.cluster_centroids = None
        self.volume = volume
        self.electrons = electrons
        self.model = []
        self.solution = []
        self.multiplicity = 1