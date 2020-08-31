from __future__ import division

from olexFunctions import OlexFunctions
OV = OlexFunctions()
import sys
import re
import os
import numpy as np
import parser
import solvator_atom
import void
import copy

class Host():

    def __init__(self, filename, FVAR_floats):
        
        self.filename = filename
        resfile = filename + ".res"
        atom_pattern1 = re.compile(r'^[A-Z]+[a-z]*[0-9]+[A-Za-z]*\s+[0-9]+(\s+-*[0-9]+\.[0-9]+){3}.*')
        atom_pattern2 = re.compile(r'(\s+-*[0-9]+\.[0-9]+)+')
        HKLF_pattern = re.compile(r'^HKLF')
        atom_lines = []
        self.uniq_atoms = []
        self.voids = None
        self.uniq_voids = []
        self.clusters = []
        self.total_missing_electrons = None
        self.total_void_volume = None
        self.average_uiso = 0

        with open(resfile, 'r') as fh:
            for line in fh.readlines():
                if re.match(atom_pattern1, line) or re.match(atom_pattern2,line):
                    atom_lines.append(line)        
                if re.match(HKLF_pattern, line): # stop when HKLF instruction is reached
                    break

        for i, line in enumerate(atom_lines):
            if re.match(r'^\s*[A-Z]', line): # only consider lines that start with an atomic symbol
                line = line.split()
                if line[-1] == "=":
                    next_line = atom_lines[i+1].split() # joins lines that end with an equal sign to the following line.
                    line = line[:-1] + next_line
                atom_name = line[0].strip()
                atom_symbol = re.sub(r'([A-Z]+[a-z]*)\d*.*', r'\1', atom_name).lower().capitalize()
                atom_label = re.sub(r'[A-Z]+[a-z]*(\d*.*)', r'\1', atom_name)
                atom_sfac = int(line[1].strip())
                atom_x = float(line[2].strip())
                atom_y = float(line[3].strip())
                atom_z = float(line[4].strip())
                atom_coords = np.array([[atom_x],[atom_y],[atom_z]])
                line[5] = line[5].strip()
                split_point = re.search("\.", line[5]).start()
                atom_fvar = int(line[5][:split_point-1])
                # This next section takes into account the free variables and occupancies of each of the atoms to calculate the overall occupancy.
                if atom_fvar > 0:
                    atom_occupancy = float(line[5][split_point-1:].strip())*FVAR_floats[atom_fvar-1]
                else:
                    atom_occupancy = float(line[5][split_point-1:].strip())*(1-FVAR_floats[abs(atom_fvar)-1])
                #print len(line)
                if len(line) == 7:
                    atom_u = None
                    atom_uiso = float(line[6].strip())
                elif len(line) == 12:
                    atom_u = [float(line[6].strip()), float(line[7].strip()), float(line[8].strip()), float(line[9].strip()), float(line[10].strip()), float(line[11].strip())]
                    atom_uiso =  (atom_u[0] + atom_u[1] + atom_u[2])/3 # approximation
                
                thisatom = solvator_atom.Solvator_atom(atom_symbol, atom_label, atom_sfac, atom_coords, atom_fvar, atom_occupancy, atom_uiso, atom_u, 'frac')
                self.uniq_atoms.append(thisatom)
        
        print "Number of unique atoms: %s" % len(self.uniq_atoms)
                
        for atom in self.uniq_atoms:
            index = self.uniq_atoms.index(atom)
            if atom.uiso < 0: 
                for i in range(1,4):
                    if self.uniq_atoms[index-i].uiso > 0:
                        #print "This is a hydrogen riding on atom", atoms[index-i][0]
                        atom.uiso = atom.uiso * self.uniq_atoms[index-i].uiso # Calculates the value of the Uiso for the hydrogens based on their parent atom

        
    def get_voids(self):
        
        ciffile = self.filename + ".cif"
        
        if os.path.exists(ciffile):

            with open(ciffile, 'r') as fhand:
                line_num = 0
                start_line = 0

                lines = fhand.readlines()
                for line in lines:
                    line_num += 1
                    if "_void_content" in line: # search for Squeeze/mask analysis in CIF
                        start_line = line_num
                
                if start_line == 0:
                    print "No voids found"
                
                search_phrase = re.compile(r'^[\s]*([0-9]+)[\s]+(-*[0-9]+\.[0-9]+)[\s]+(-*[0-9]+\.[0-9]+)[\s]+(-*[0-9]+\.[0-9]+)[\s]+([0-9]+\.*[0-9]*)[\s]+([0-9]+\.*[0-9]*)[\s]+.*')

                self.voids = []
                while re.match(search_phrase,lines[start_line]):
                    m = lines[start_line].split()
                    void_number = int(m[0].strip())
                    #void_centroid = np.array([[float(m[1])],[float(m[2])],[float(m[3])]])
                    void_volume = float(m[4].strip())
                    void_content = float(m[5].strip())
                    thisvoid = void.Void(void_number, void_volume, void_content)
                    self.voids.append(thisvoid)
                    start_line +=1
                            
        self.total_missing_electrons = 0
        self.total_void_volume = 0
        
        for m_void in self.voids:
            self.total_missing_electrons += m_void.electrons
            self.total_void_volume += m_void.volume
            
            if self.uniq_voids == []:
                self.uniq_voids.append(copy.deepcopy(m_void))
            else:
                print m_void.number, m_void.volume, m_void.electrons

                for uniq_void in self.uniq_voids:
                    if (m_void.volume == uniq_void.volume) and (m_void.electrons == uniq_void.electrons):
                        uniq_void.multiplicity += 1
                    else:
                        self.uniq_voids.append(copy.deepcopy(m_void))
                    
        total_electrons = 0
        total_volume = 0
        for uniq_void in self.uniq_voids:
            total_electrons += uniq_void.electrons
            total_volume += uniq_void.volume
        
        self.total_uniq_missing_electrons = total_electrons
        self.total_uniq_void_volume = total_volume