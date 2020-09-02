from __future__ import division
from olexFunctions import OlexFunctions
OV = OlexFunctions()

import numpy as np
import solvator_atom
import copy
import re
import olx
import olexex
import os
from shutil import copyfile

class Guest():
    
    def __init__(self, solvent_name, occupancy, p_path, SFAC_by_symbol, CART_TO_FRAC, cluster_centroid = np.array([[0],[0],[0]])):
    #def __init__(self, solvent_name, occupancy, p_path, SFAC_by_symbol, cluster_centroid = np.array([[0],[0],[0]])):

        filename = os.path.join(p_path, "solvents/" ) + str(solvent_name) +".guest"
        self.atoms = [] # these atom coordinates are the original Cartesian coordinates and read from the .guest file and should remain unchanged.
        self.transformed_atoms = [] # these atom coordinates are fractional coordinates.
        self.occupancy = occupancy
        self.fixed_occupancy = False
        
        with open(filename) as fh:
            header = next(fh).split()
            self.name = header[0]
            self.electrons = float(header[1])
            self.volume = float(header[2])
            self.original_restraints = []
            self.original_hfix = []
            if next(fh).split()[0] == "restraints":
                line = next(fh)
                while line.split()[0].strip() != "end":
                    self.original_restraints.append(line.strip())
                    line = next(fh)
            if next(fh).split()[0] == "hydrogens":
                line = next(fh)
                while line.split()[0].strip() != "end":
                    self.original_hfix.append(line.strip())
                    line = next(fh)                    
            
            number = 0
            for line in fh:
                number += 1
                line = line.split()
                #atom_symbol = line[0]
                #atom_label = "%sJ" % number
                r = re.compile("([a-zA-Z]+)([0-9]*)([A-Za-z]*)")
                atom_symbol = r.match(line[0]).group(1)
                if r.match(line[0]).group(2):
                    atom_number = r.match(line[0]).group(2)
                else: 
                    atom_number = number
                if r.match(line[0]).group(3):
                    atom_suffix = r.match(line[0]).group(3)
                    atom_label = str(atom_number) + atom_suffix
                else:
                    atom_label = str(atom_number)
                
                # Check all atoms have SFAC number. If not, create new one sequentially.
                if atom_symbol not in SFAC_by_symbol:
                    length_of_dict = len(SFAC_by_symbol)
                    SFAC_by_symbol[atom_symbol] = (length_of_dict, 4) # adds the number 4 to the unit line. Perhaps update to count actual total number of atoms in file?
                    #SFAC_BY_NUMBER[length_of_dict] = (atom_symbol, 4)
                    print 'Element', atom.symbol, 'appended to SFAC dictionaries'
                atom_sfac = SFAC_by_symbol[atom_symbol][0]
                x_cart = float(line[1])
                y_cart = float(line[2])
                z_cart = float(line[3])
                coords_cart = np.array([[x_cart],[y_cart],[z_cart]])

                atom = solvator_atom.Solvator_atom(atom_symbol, atom_label, atom_sfac, coords_cart, None, self.occupancy,0.03, None, 'cart') 
                self.atoms.append(atom)
        
        # Create the first set of fractional coordinates and translate such that the centroid of the molecule is coincident with
        # the centroid of the cluster
        self.transformed_atoms = copy.deepcopy(self.atoms)
        for atom in self.transformed_atoms:
            atom.coords = np.dot(CART_TO_FRAC, atom.coords)
        current_centroid = self.calculate_centroid()
        self.translate(cluster_centroid - current_centroid)
        
        max_distance = 0
        for atom1 in self.atoms:
            for atom2 in self.atoms:
                difference = atom1.coords - atom2.coords
                distance = np.linalg.norm(difference)
                if distance > max_distance:
                    max_distance = distance
        
        self.max_dimension = max_distance


    def find_highest_atom_number_in_model(self):

        highest_atom_number = 0

        for atom in olexex.OlexRefinementModel()._atoms:
            #print atom['aunit_id'], atom['label']
            atom_label = re.match('[A-Za-z]+([0-9]+)[A-Za-z\']*', atom['label'] )
            if not atom_label:
                continue
            else:
                number = int(atom_label.group(1))
                if number > highest_atom_number:
                    highest_atom_number = number

        return highest_atom_number
        
    def add_to_olex_model(self, part):
        atom_number = self.find_highest_atom_number_in_model()            
        
        solvator_number = 0
        for atom in self.transformed_atoms:
            solvator_number += 1
            atom.solvator_id = "SOLV" + str(solvator_number)        
        
        #for i, atom in enumerate(olexex.OlexRefinementModel()._atoms):
        #    print i, olx.xf.au.GetAtomName(atom['aunit_id']), olx.xf.au.GetAtomCrd(atom['aunit_id'])
        print "Adding %s to the olex model" % self.name
        #use atom's old name, atom's new name and then replace in restraints?-------------?
        for atom in self.transformed_atoms:
            number_atoms_before = len(olexex.OlexRefinementModel()._atoms)
            atom_number +=1
            atom.olex_id = str(atom.symbol) + str(atom_number)
            #name = str(atom.symbol) + str(atom_number) +"J"
            olx.xf.au.NewAtom(atom.olex_id, float(atom.coords[0][0]),float(atom.coords[1][0]),float(atom.coords[2][0]))
            number_atoms_after = len(olexex.OlexRefinementModel()._atoms)
            if number_atoms_after == number_atoms_before:
                while number_atoms_after == number_atoms_before:
                    atom.shift_slightly()
                    olx.xf.au.NewAtom(atom.olex_id, float(atom.coords[0][0]),float(atom.coords[1][0]),float(atom.coords[2][0]))
                    number_atoms_after = len(olexex.OlexRefinementModel()._atoms)                
                print("WARNING! Atom %s%s has been shifted in order to add it successfully.") % (atom.symbol, atom.label)

        for atom in olexex.OlexRefinementModel()._atoms:
            for guest_atom in self.transformed_atoms:
                if(re.match(guest_atom.olex_id, atom['label'])):
                    olx.xf.au.SetAtomOccu(atom['aunit_id'], guest_atom.occupancy)
                    olx.xf.au.SetAtomU(atom['aunit_id'], guest_atom.uiso)
                    olx.xf.au.SetAtomPart(atom['aunit_id'], part)
            
        OV.cmd('fuse')
        if self.original_restraints == []:
            print "No restraints for this molecule"
        else:
        # Find and replace solvator_id with olex_id in each restraint. This has to be done in stages
        # just in case the guest atom's original name overlaps with any of the names of the host atoms.
            self.restraints = copy.deepcopy(self.original_restraints)
            for n in range(len(self.original_restraints)):
                for atom in self.transformed_atoms:
                    self.restraints[n] = self.restraints[n].replace(atom.name, atom.solvator_id)
                for atom in self.transformed_atoms:
                    self.restraints[n] = self.restraints[n].replace(atom.solvator_id, atom.olex_id)
        
        for restraint in self.restraints:
            OV.cmd(restraint)
            
    
    def remove_from_olex_model(self):
        for atom in olexex.OlexRefinementModel()._atoms:
            for guest_atom in self.transformed_atoms:
                if(re.match(guest_atom.olex_id, atom['label'])):
                    OV.cmd('kill %s' % (guest_atom.olex_id))
        
    def update_coordinates_from_olex_model(self):
        
        for atom in olexex.OlexRefinementModel()._atoms:
            self.occupancy = float(olx.xf.au.GetAtomOccu(atom['aunit_id'])) # so the occupancy of the molecule will match that of the last atom! :)
            for guest_atom in self.transformed_atoms:
                if(re.match(guest_atom.olex_id, atom['label'])):
                    guest_atom.coords[0][0], guest_atom.coords[1][0], guest_atom.coords[2][0] = olx.xf.au.GetAtomCrd(atom['aunit_id']).split(',')
                    guest_atom.uiso = float(olx.xf.au.GetAtomUiso(atom['aunit_id']))
                    guest_atom.occupancy = float(olx.xf.au.GetAtomOccu(atom['aunit_id']))
                    
        
        return self
    
    
    def add_to_olex_model_via_resfile(self):
        
        # We can only ever reap the resfile, otherwise, we will be opeining a file were the filename of the res
        # does not match the hkl filename.
        
        atom_number = self.find_highest_atom_number_in_model()
        
        solvator_number = 0
        for atom in self.transformed_atoms:
            atom_number +=1
            solvator_number += 1
            atom.solvator_id = "SOLV" + str(solvator_number)
            atom.olex_id = atom.symbol + str(atom_number)
        
        resfile = os.path.join(OV.FilePath(), OV.FileName()) + ".res"
        original_resfile = os.path.join(OV.FilePath(), OV.FileName()) + "_original.res"
        
        with open(original_resfile, 'r') as infile, open(resfile, 'w') as outfile:
            #atom_line = re.compile(r'^([a-zA-Z]+[0-9]+[A-Za-z]*\s+[0-9]+)\s+(-*)([0-9]\.[0-9]+)\s+(-*)([0-9]\.[0-9]+)\s+(-*)([0-9]\.[0-9]+)(.*)')
            for line in infile.readlines():
                #print line
                if re.match(r'^FVAR.*', line):
                    line = line.strip().split()
                    number_of_FVAR = len(line)
                    for item in line:
                        outfile.write("%s " % item)
                    outfile.write("%.4f\n" % self.occupancy)
                #elif re.match(atom_line, line):
                    #outfile.write("%s  %s1%.6f  %s1%.6f  %s1%.6f%s\n" % (atom_line.match(line).group(1), atom_line.match(line).group(2),float(atom_line.match(line).group(3)), atom_line.match(line).group(4), float(atom_line.match(line).group(5)), atom_line.match(line).group(6), float(atom_line.match(line).group(7)), atom_line.match(line).group(8)))
                elif re.match(r'HKLF.*', line):
                    outfile.write("PART -%s\n" % number_of_FVAR)
                    for atom in self.transformed_atoms:
                        outfile.write(("%s %s %.6f %.6f %.6f %s1.000000 %.6f\n") % (atom.olex_id, atom.sfac, atom.coords[0], atom.coords[1], atom.coords[2], number_of_FVAR,  atom.uiso))#, number_of_FVAR))
                        solvator_number += 1
                    outfile.write("PART 0\n")
                    outfile.write(line)
                else:
                    outfile.write(line)
        
        olx.Atreap(resfile)
        print "Number of atoms after writing resfile: %s" % len(olexex.OlexRefinementModel()._atoms)
        
        OV.cmd('fuse')
        if self.original_restraints == []:
            print "No restraints for this molecule"
        else:
        # Find and replace solvator_id with olex_id in each restraint. This has to be done in stages
        # just in case the guest atom's original name overlaps with any of the names of the host atoms.
            self.restraints = copy.deepcopy(self.original_restraints)
            for n in range(len(self.original_restraints)):
                for atom in self.transformed_atoms:
                    self.restraints[n] = self.restraints[n].replace(atom.name, atom.solvator_id)
                for atom in self.transformed_atoms:
                    self.restraints[n] = self.restraints[n].replace(atom.solvator_id, atom.olex_id)
        
        for restraint in self.restraints:
            OV.cmd(restraint)

        R_current = float(str(olx.CalcR()).split(',')[0])

        return R_current

            
    def add_hfix(self):
        
        if self.original_hfix == []:
            print "No HFIX supplied for this molecule"
        else:
            self.hfix = copy.deepcopy(self.original_hfix)
        # Find and replace solvator_id with olex_id in each HFIX. This has to be done in stages
        # just in case the guest atom's original name overlaps with any of the names of the host atoms.
            for n in range(len(self.hfix)):
                for atom in self.transformed_atoms:
                    #print atom.name, atom.solvator_id
                    self.hfix[n] = self.hfix[n].replace(atom.name, atom.solvator_id)
                for atom in self.transformed_atoms:
                    #print atom.solvator_id, atom.olex_id
                    self.hfix[n] = self.hfix[n].replace(atom.solvator_id, atom.olex_id)
        
        for hfix in self.hfix:
            OV.cmd(hfix)

    def calculate_centroid(self):

        x_centroid = 0
        y_centroid = 0
        z_centroid = 0

        for atom in self.transformed_atoms:
            x_centroid += atom.coords[0]
            y_centroid += atom.coords[1]
            z_centroid += atom.coords[2]
        x_centroid = x_centroid/len(self.transformed_atoms)
        y_centroid = y_centroid/len(self.transformed_atoms)
        z_centroid = z_centroid/len(self.transformed_atoms)
        centroid = np.array([x_centroid,y_centroid,z_centroid])

        return centroid

    def change_occupancy(self, value):
        
        self.occupancy += value
        
        if self.occupancy < 0.05:
            self.occupancy = 0.05
        if self.occupancy > 0.95:
            self.occupancy = 0.95

        for atom in self.transformed_atoms:
            atom.occupancy = self.occupancy
    
    def set_occupancy(self, value):
        
        self.occupancy = value
        
        if self.occupancy < 0.05:
            self.occupancy = 0.05
        if self.occupancy > 0.95:
            self.occupancy = 0.95

        for atom in self.transformed_atoms:
            atom.occupancy = self.occupancy    

    def translate(self, translation):
        
        for atom in self.transformed_atoms:
            atom.coords += translation

    def rotate(self, parameters, cart_to_frac, frac_to_cart, centre_of_rotation):
        
        """""
        This method should rotate the molecule around the centroid given. 
        """""
        
        phi = parameters[0]
        theta = parameters[1]
        psi = parameters[2]
        
        #print "phi, theta, psi are: %.1f, %.1f, %.1f degrees." % (np.degrees(phi), np.degrees(theta), np.degrees(psi))

        a1 = np.cos(phi)
        b1 = -np.sin(phi)
        c1 = np.sin(phi)
        d1 = np.cos(phi)
        e1 = np.cos(theta)
        f1 = np.sin(theta)
        g1 = -np.sin(theta)
        h1 = np.cos(theta)
        i1 = np.cos(psi)
        j1 = -np.sin(psi)
        k1 = np.sin(psi)
        l1 = np.cos(psi)

        # rotation matrix for rotation around angles phi, theta, psi around x, y, and z axes respectively
        rotator = np.array([
                        [e1*i1           , e1*j1,            f1    ],
                        [a1*k1 + b1*g1*i1, a1*l1 + b1*g1*j1, b1*h1 ],
                        [c1*k1 + d1*g1*i1, c1*l1 + d1*g1*j1, d1*h1 ]
                          ])

        for atom in self.transformed_atoms:

            atom.coords -= centre_of_rotation
            atom.coords = np.dot(frac_to_cart, atom.coords)
            atom.coords = np.dot(rotator, atom.coords)
            atom.coords = np.dot(cart_to_frac, atom.coords)
            atom.coords += centre_of_rotation

    # Perhaps later we might need to rotate individually        
    #rotator_x = np.array([
                       #[1, 0,         0             ],
                       #[0, np.cos(self.phi), -np.sin(self.phi)],
                       #[0, np.sin(self.phi),  np.cos(self.phi)]
                      #])
    
    #rotator_y = np.array([
                        #[ np.cos(self.theta), 0, np.sin(self.theta)],
                        #[0             , 1, 0            ],
                        #[-np.sin(self.theta), 0, np.cos(self.theta)]
                       #])  
                       
    #rotator_z = np.array([
                       #[np.cos(self.psi), -np.sin(self.psi), 0],
                       #[np.sin(self.psi),  np.cos(self.psi), 0],
                       #[0,0,1]
                      #])  
    def refine_using_olex_model(self, fixed_position = False, fixed_occupancy = False, fixed_iso = False):
        
        #resfile = os.path.join(OV.FilePath(), OV.FileName()) + ".res"
        original_resfile = os.path.join(OV.FilePath(), OV.FileName()) + "_original.res"

        #copyfile(original_resfile, resfile)
        
        olx.Atreap(original_resfile)
        
        
        #print "Page 3: Number of atoms after copying original_resfile to resfile: %s" % (len(olexex.OlexRefinementModel()._atoms))
        #for i, atom in enumerate(olexex.OlexRefinementModel()._atoms):
            #print i, olx.xf.au.GetAtomName(atom['aunit_id']), olx.xf.au.GetAtomCrd(atom['aunit_id'])
        
        if fixed_position == True:
            fp = 1
        else:
            fp = ""
        
        if fixed_occupancy == True:
            fo = 1
        else:
            fo = ""        
        
        if fixed_iso == True:
            fi = 1
        else:
            fi = ""
            
        self.add_to_olex_model()
        
        print "Number of atoms after adding_to_olex_model: %s" % len(olexex.OlexRefinementModel()._atoms)

        self.add_restraints()

        OV.cmd('refine 10 5')

        R_current = float(str(olx.CalcR()).split(',')[0])

        self.update_coordinates_from_olex_model()
        self.remove_from_olex_model()

        return [R_current, self]    