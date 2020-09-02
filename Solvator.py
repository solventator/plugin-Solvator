from __future__ import division
from olexFunctions import OlexFunctions

OV = OlexFunctions()

import os
import htmlTools
import olex
import olx
import gui

import time
debug = bool(OV.GetParam("olex2.debug", False))

instance_path = OV.DataDir()

try:
  from_outside = False
  p_path = os.path.dirname(os.path.abspath(__file__))
except:
  from_outside = True
  p_path = os.path.dirname(os.path.abspath("__file__"))

l = open(os.sep.join([p_path, 'def.txt'])).readlines()
d = {}
for line in l:
  line = line.strip()
  if not line or line.startswith("#"):
    continue
  d[line.split("=")[0].strip()] = line.split("=")[1].strip()

p_name = d['p_name']
p_htm = d['p_htm']
p_img = eval(d['p_img'])
p_scope = d['p_scope']

OV.SetVar('Solvator_plugin_path', p_path)

from PluginTools import PluginTools as PT

# My imports
import sys
import cctbx_olex_adapter
import re
import random
import numpy as np
import time
import copy
import host
import guest
import solvator_atom
import q_cluster
from cctbx import maptbx, miller, uctbx
import parser
import scipy.optimize

#from common_solvents import common_solvents
import olexex
from shutil import copyfile
from itertools import combinations, combinations_with_replacement, permutations
from operator import itemgetter

class Solvator(PT):

  def __init__(self):
    super(Solvator, self).__init__()
    self.p_name = p_name
    self.p_path = p_path
    self.p_scope = p_scope
    self.p_htm = p_htm
    self.p_img = p_img
    self.deal_with_phil(operation='read')
    self.print_version_date()
    if not from_outside:
      self.setup_gui()
    # END Generated =======================================

    OV.registerFunction(self.save_to_phil,True,"solvator")    
    OV.registerFunction(self.restore_original_resfile,True,"solvator")
    OV.registerFunction(self.solve,True,"solvator")
    OV.registerFunction(self.open_solvent_files_location,True,"solvator")
    
    print "\nWelcome to the SOLVATOR!"

  def open_solvent_files_location(self):
    
    olx.Shell(os.path.join(p_path, "solvents/"))
  
  def save_to_phil(self):

    self.deal_with_phil(operation='save')

  def read_res(self,resfile):
    
    """"
    This method is designed to parse the resfile and extract the parameters that are used by the SOLVATOR module.
    
    """
  
    try:
      with open(resfile, 'r') as fh:
  
        global SYMMETRY_GENERATORS
        global LATT
        SYMMETRY_GENERATORS = []
        ZERR_pattern = re.compile(r'^ZERR.*')
        SYMM_pattern = re.compile(r'^SYMM(.*)')
        FVAR_pattern = re.compile(r'FVAR.*')
        SFAC_pattern = re.compile(r'SFAC.*')
        UNIT_pattern = re.compile(r'UNIT.*')
        #SFAC_BY_NUMBER = {}
        SFAC_BY_SYMBOL = {}
        FVAR = [1] # Just a dummy value for FVAR[0] so that the FVAR matches the actual position on the line.
  
        for line in fh.readlines():
         
          if re.match('^CELL.*', line):
            CELL = line.split()
            wavelength = float(CELL[1].strip())
            a = float(CELL[2].strip())
            b = float(CELL[3].strip())
            c = float(CELL[4].strip())
            ALPHA = float(CELL[5].strip())
            BETA = float(CELL[6].strip())
            GAMMA = float(CELL[7].strip())
  
            if (a == b) and (a == c) and ALPHA == 90.0 and BETA == 90.0 and GAMMA == 90.0:
              crystal_system = "cubic"
            elif (a == b) and (a != c) and ALPHA == 90.0 and BETA == 90.0 and GAMMA == 90.0:
              crystal_system = "tetragonal"
            elif (a == b) and (a == c) and ALPHA == 90.0 and BETA == 90 and GAMMA == 120.0:
              crystal_system = "trigonal"
            elif ALPHA == 90 and BETA == 90.0 and GAMMA == 90.0:
              crystal_system = "orthorhombic"
            elif ALPHA == 90 and GAMMA == 90.0:
              crystal_system = "monoclinic"
            else:
              crystal_system = "triclinic"
          
          
          m = SYMM_pattern.match(line)
          if m:
            symmetry_operator = m.group(1)
            symmetry_operator.replace(" ","")
            symmetry_operator = symmetry_operator.strip()
            symmetry_operator = symmetry_operator.split(",")
            SYMMETRY_GENERATORS.append(symmetry_operator)
          
          if re.match(SFAC_pattern, line):
            elements = line.split()
          
          if re.match(UNIT_pattern, line):
            units = line.split()
          if re.match(r'FVAR.*', line):
            FVAR += line.split()[1:]
          if re.match(ZERR_pattern, line):
            Z = int(line.split()[1])
          if re.match(r'LATT.*', line):
            LATT = int(line.split()[1])
          

        volume = float(olx.xf.au.GetCellVolume().strip())
        FVAR_FLOATS = [float(item) for item in FVAR]
        OSF = FVAR_FLOATS[1]**2        
        #generate SFAC as dictionary with numbers or symbols as keys:
        #print elements, units
        for i in range(1,len(elements)):
          #SFAC_BY_NUMBER[i] = (elements[i],int(units[i]))
          SFAC_BY_SYMBOL[elements[i]] = (i, int(units[i]))
        
        ALPHA_RAD = np.deg2rad(ALPHA)
        BETA_RAD = np.deg2rad(BETA)
        GAMMA_RAD = np.deg2rad(GAMMA)
        
        SINALPHA = np.sin(ALPHA_RAD)
        SINBETA = np.sin(BETA_RAD)
        SINGAMMA = np.sin(GAMMA_RAD)
        COSALPHA = np.cos(ALPHA_RAD)
        COSBETA = np.cos(BETA_RAD)
        COSGAMMA = np.cos(GAMMA_RAD)

        global CART_TO_FRAC
        
        CART_TO_FRAC = np.zeros((3, 3))
        CART_TO_FRAC[0, 0] = 1.0 / a
        CART_TO_FRAC[0, 1] = -COSGAMMA / ( a *  SINGAMMA)
        CART_TO_FRAC[0, 2] =  b * c * ( COSALPHA *  COSGAMMA -  COSBETA) / ( volume *  SINGAMMA)
        CART_TO_FRAC[1, 1] = 1.0 / ( b *  SINGAMMA)
        CART_TO_FRAC[1, 2] =  a * c * ( COSBETA *  COSGAMMA -  COSALPHA) / ( volume *  SINGAMMA)
        CART_TO_FRAC[2, 2] =  a * b * SINGAMMA /  volume
        
        global FRAC_TO_CART
        
        FRAC_TO_CART = np.zeros((3,3))
        FRAC_TO_CART[0,0] = a
        FRAC_TO_CART[0,1] = b * COSGAMMA
        FRAC_TO_CART[0,2] = c * COSBETA
        FRAC_TO_CART[1,1] = b * SINGAMMA
        FRAC_TO_CART[1,2] = c * (COSALPHA - COSBETA * COSGAMMA)/SINGAMMA
        FRAC_TO_CART[2,2] = volume/(a * b * SINGAMMA)


        return (crystal_system,
                OSF,
                FVAR_FLOATS,
                SFAC_BY_SYMBOL
                )
  
    except:
      print "There was a problem reading the res file."

  
  def calculate_distance(self, atom1, atom2):#):
    
    #Calculates the distance between two atoms of the SOLVATOR_atom
    
    atom1_cart = np.dot(FRAC_TO_CART, atom1.coords)
    atom2_cart = np.dot(FRAC_TO_CART, atom2.coords)
    
    difference = atom2_cart - atom1_cart
    distance = np.linalg.norm(difference)
    
    return distance
  
  def calculate_plane(self, atom1, atom2, atom3):#):
    
    # Will only work for atoms with fractional coordinates. 
    coords1 = np.dot(FRAC_TO_CART, atom1.coords)
    coords2 = np.dot(FRAC_TO_CART, atom2.coords)
    coords3 = np.dot(FRAC_TO_CART, atom3.coords)
    
    p1 = np.array([coords1[0][0],coords1[1][0],coords1[2][0]])
    p2 = np.array([coords2[0][0],coords2[1][0],coords2[2][0]])
    p3 = np.array([coords3[0][0],coords3[1][0],coords3[2][0]])
    
    # These two vectors are in the plane
    v1 = p3 - p1
    v2 = p2 - p1
    
    # the cross product is a vector normal to the plane
    cp = np.cross(v1, v2)
    a, b, c = cp
    
    # This evaluates a * x3 + b * y3 + c * z3 which equals d
    d = np.dot(cp, p3)
    
    print('The equation is {0}x + {1}y + {2}z = {3}'.format(a, b, c, d))
    
    return (a,b,c,d)
  
  def calculate_normal(self, atom1, atom2, atom3):
    
    # Will calculate the normal to the plave of three atoms. Will only work for atoms with fractional coordinates. 
    
    coords1 = np.dot(FRAC_TO_CART, atom1.coords)
    coords2 = np.dot(FRAC_TO_CART, atom2.coords)
    coords3 = np.dot(FRAC_TO_CART, atom3.coords)
    
    p1 = np.array([coords1[0][0],coords1[1][0],coords1[2][0]])
    p2 = np.array([coords2[0][0],coords2[1][0],coords2[2][0]])
    p3 = np.array([coords3[0][0],coords3[1][0],coords3[2][0]])
    
    # These two vectors are in the plane
    v1 = p3 - p1
    v2 = p2 - p1
    
    # the cross product is a vector normal to the plane
    cp = np.cross(v1, v2)
    a, b, c = cp
    #print "cp:", cp
    
    # This evaluates a * x3 + b * y3 + c * z3 which equals d
    
    #print('The normal is {0}x + {1}y + {2}z'.format(a, b, c))
    
    return cp
  
  def calculate_rotation(self,x,y,z,a,b,c,u,v,w,theta):
    L = u **2 + v **2 + w **2
    sqrt_L = np.sqrt(L)
    #u,v,w, = vector of intersection of two planes
    #a,b,c = point on intersection
    #x,y,z = coordinates of point to be rotated.
    x1 = ((a * (v**2 + w**2) - u * ( b*v + c*w - u*x - v*y - w*z) ) * (1 - np.cos(theta)) + L * x * np.cos(theta) + sqrt_L * (-c*v + b*w - w*y + v*z) * np.sin(theta))/L
    y1 = ((b * (u**2 + w**2) - v * ( a*u + c*w - u*x - v*y - w*z) ) * (1 - np.cos(theta)) + L * y * np.cos(theta) + sqrt_L * ( c*u - a*w + w*x - u*z) * np.sin(theta))/L
    z1 = ((c * (u**2 + v**2) - w * ( a*u + b*v - u*x - v*y - w*z) ) * (1 - np.cos(theta)) + L * z * np.cos(theta) + sqrt_L * (-b*u + a*v - v*x + u*y) * np.sin(theta))/L
    
    return np.array([[x1],[y1],[z1]])  
  
  
  def calculate_angle_between_planes(self, plane1, plane2):
    
    #The planes should be a tuple (a,b,c,d) of three (for a vector) or four numbers (for a plane) representing ax + by + cz = d
    
    a1 = plane1[0]
    b1 = plane1[1]
    c1 = plane1[2]
    a2 = plane2[0]
    b2 = plane2[1]
    c2 = plane2[2]

    d = ( a1 * a2 + b1 * b2 + c1 * c2 ) 
    e1 = np.sqrt( a1 * a1 + b1 * b1 + c1 * c1)
    e2 = np.sqrt( a2 * a2 + b2 * b2 + c2 * c2) 
    d = d / (e1 * e2)
    if (np.isclose(1, abs(d), 0.001)):
      #print "Planes are coplanar"
      return 0
    else:
      angle = np.arccos(d)
      #print("Angle between planes %.1f") % np.degrees(angle) 

      return angle
  
  def calculate_angle_between_atoms(self, atom1, atom2, atom3):
    
    # This calculates the angle at atom2
    
    a = self.calculate_distance(atom1, atom3)
    b = self.calculate_distance(atom2, atom3)
    c = self.calculate_distance(atom1, atom2)
    
    #print atom1.coords, atom2.coords, atom3.coords
    
    value = (b**2 + c**2 - a**2)/(2*b*c)
    if value > 1.0:
      "Impossible value %.12f %s%s - %s%s - %s%s" % (value, atom1.symbol, atom1.label, atom2.symbol, atom2.label, atom3.symbol, atom3.label)
      value = 1
    elif value < -1.0:
      "Impossible value %.12f %s%s - %s%s - %s%s" % (value, atom1.symbol, atom1.label, atom2.symbol, atom2.label, atom3.symbol, atom3.label)
      value = -1 
    angle = np.arccos(value)
    #print "%s%s - %s%s - %s%s: %.1f" % (atom1.symbol, atom1.label, atom2.symbol, atom2.label, atom3.symbol, atom3.label, np.degrees(angle))
      
    return angle
  
  def get_vector(self, atom1, atom2):
    
    coords1 = np.dot(FRAC_TO_CART, atom1.coords)
    coords2 = np.dot(FRAC_TO_CART, atom2.coords)
    [x],[y],[z] = coords2 - coords1
    
    return x,y,z
  
  def is_chemically_sensible(self, guest_atoms, host_atoms):
    #print "Searching for chemically sensible model..."
    model_atoms = []
    
    model_atoms += self.expand_atoms_with_symmetry(guest_atoms)
    model_atom = self.move_all_atoms_into_cell(model_atoms)
    
    for model_atom in model_atoms:
      for host_atom in host_atoms:
        length = self.calculate_distance(model_atom, host_atom)
        if length < short_contact:
          #print "Found a distance of %.2f" % (length)
          return False
        else:
          continue
  
    #print "...done"
    return True
    
  def make_distance_matrix(self, atom_list):#):
    
    # This creates a distance matrix of all the atoms in the list.
    
    distance_matrix = np.zeros([len(atom_list), len(atom_list)])
    for i, peak_i in enumerate(atom_list):
      for j, peak_j in enumerate(atom_list):
        distance = self.calculate_distance(peak_i, peak_j)
        distance_matrix[i][j] = distance
    
    return distance_matrix
  
  def expand_atoms_with_symmetry(self, uniq_atoms): 
  
    """""Lattice type: 1=P, 2=I, 3=rhombohedral obverse on hexagonal axes, 4=F, 5=A, 6=B, 7=C. 
    N must be made negative if the structure is non-centrosymmetric. 
    
    Fdd2 	  0.125 0.125 0.500   	MOVE 0.25 0.25 1 -1
    I41 	  0.500 0.250 0.500   	MOVE 1 0.5 1 -1
    I4122 	  0.500 0.250 0.125   	MOVE 1 0.5 0.25 -1
    I41md 	  0.500 0.250 0.500   	MOVE 1 0.5 1 -1
    I41cd 	  0.500 0.250 0.500   	MOVE 1 0.5 1 -1
    I42d 	  0.500 0.250 0.125   	MOVE 1 0.5 0.25 -1
    F4132 	  0.125 0.125 0.125   	MOVE 0.25 0.25 0.25 -1
    """
    #print "Symmetry generators:", SYMMETRY_GENERATORS
    atoms = copy.deepcopy(uniq_atoms)
    
    if abs(LATT) == 1: # P-centred
      pass
    
    elif abs(LATT) == 2: # I-centred
      new_atoms = []
      for atom in uniq_atoms:
        new_atom = copy.deepcopy(atom)
        new_x = new_atom.coords[0] + 0.5
        new_y = new_atom.coords[1] + 0.5
        new_z = new_atom.coords[2] + 0.5
        new_atom.coords = (np.array([new_x,new_y,new_z]))
        new_atoms.append(new_atom) 
      atoms += new_atoms
    
    elif abs(LATT) == 3: # Rhombohedral obverse on hexagonal axes"
      new_atoms = []
      for atom in uniq_atoms:
        new_atom1 = copy.deepcopy(atom)
        new_atom2 = copy.deepcopy(atom)
        new_x1 = new_atom1.coords[0] + 0.666666
        new_y1 = new_atom1.coords[1] + 0.333333
        new_z1 = new_atom1.coords[2] + 0.333333
        new_x2 = new_atom2.coords[0] + 0.333333
        new_y2 = new_atom2.coords[1] + 0.666666
        new_z2 = new_atom2.coords[2] + 0.666666
        new_atom1.coords = (np.array([new_x1,new_y1,new_z1]))
        new_atom2.coords = (np.array([new_x2,new_y2,new_z2]))
        new_atoms.append(new_atom1)
        new_atoms.append(new_atom2) 
      atoms += new_atoms
    
    elif abs(LATT) == 4: # F-centred
      new_atoms = []
      for atom in uniq_atoms:
        new_atom1 = copy.deepcopy(atom)
        new_atom2 = copy.deepcopy(atom)
        new_atom3 = copy.deepcopy(atom)
        new_x1 = new_atom1.coords[0]
        new_y1 = new_atom1.coords[1] + 0.5
        new_z1 = new_atom1.coords[2] + 0.5
        new_x2 = new_atom1.coords[0] + 0.5
        new_y2 = new_atom1.coords[1]
        new_z2 = new_atom1.coords[2] + 0.5
        new_x3 = new_atom1.coords[0] + 0.5
        new_y3 = new_atom1.coords[1] + 0.5
        new_z3 = new_atom1.coords[2]
        new_atom1.coords = (np.array([new_x1,new_y1,new_z1]))
        new_atom2.coords = (np.array([new_x2,new_y2,new_z2]))
        new_atom3.coords = (np.array([new_x3,new_y3,new_z3]))
        new_atoms.append(new_atom1)
        new_atoms.append(new_atom2) 
        new_atoms.append(new_atom3) 
      atoms += new_atoms
    
    elif abs(LATT) == 5: # A-centred
      new_atoms = []
      for atom in uniq_atoms:
        new_atom = copy.deepcopy(atom)
        new_x = new_atom.coords[0]
        new_y = new_atom.coords[1] + 0.5
        new_z = new_atom.coords[2] + 0.5
        new_atom.coords = (np.array([new_x,new_y,new_z]))
        new_atoms.append(new_atom) 
      atoms += new_atoms
    
    elif abs(LATT) == 6: # B-centred
      new_atoms = []
      for atom in uniq_atoms:
        new_atom = copy.deepcopy(atom)
        new_x = new_atom.coords[0] + 0.5
        new_y = new_atom.coords[1]
        new_z = new_atom.coords[2] + 0.5
        new_atom.coords = (np.array([new_x,new_y,new_z]))
        new_atoms.append(new_atom) 
      atoms += new_atoms
    
    elif abs(LATT) == 7: # C-centred
      new_atoms = []
      for atom in uniq_atoms:
        new_atom = copy.deepcopy(atom)
        new_x = new_atom.coords[0] + 0.5
        new_y = new_atom.coords[1] + 0.5
        new_z = new_atom.coords[2]
        new_atom.coords = np.array([new_x,new_y,new_z])
        new_atoms.append(new_atom) 
      atoms += new_atoms
  
    for operator in SYMMETRY_GENERATORS:
      new_atoms = []
      x_equation = parser.expr(operator[0]).compile()
      y_equation = parser.expr(operator[1]).compile()
      z_equation = parser.expr(operator[2]).compile()
      for atom in atoms:
        new_atom = copy.deepcopy(atom)
        X = atom.coords[0]
        Y = atom.coords[1]
        Z = atom.coords[2]
        #print X, Y, Z, "x,y,z"
        new_x = eval(x_equation)
        new_y = eval(y_equation)
        new_z = eval(z_equation)
        new_atom.coords = np.array([new_x,new_y,new_z])
        #new_atom = self.find_closest_supercell_atom(atom, new_atom)
        new_atoms.append(new_atom)
  
      atoms += new_atoms
        
    
    if LATT > 0:  #Symmetry expansion for centrosymmetric structures.
      #But centrosymmetric space groups are not always inverted through 0 0 0, so watch out!
      #print "Structure is centrosymmetric. Adding atoms at (-x,-y,-z)"
  
      new_atoms = []
  
      for atom in atoms:
        new_atom = copy.deepcopy(atom)
        new_atom.coords = -new_atom.coords
        #new_atom = self.find_closest_supercell_atom(atom, new_atom)
        new_atoms.append(new_atom)
      
      atoms += new_atoms
      
    
    return atoms
    
  def move_all_atoms_into_cell(self, atoms):
    
    #This section moves all the atoms into the unit cell
    for atom in atoms:  
      x_int = int(atom.coords[0])
      atom.coords[0] = atom.coords[0] - x_int        
      y_int = int(atom.coords[1])
      atom.coords[1] = atom.coords[1] - y_int        
      z_int = int(atom.coords[2])
      atom.coords[2] = atom.coords[2] - z_int        
      if atom.coords[0] < 0:
        #print "a", atom.coords
        atom.coords[0] = 1 + atom.coords[0]
        #print "b", atom.coords
      if atom.coords[1] < 0:
        atom.coords[1] = 1 + atom.coords[1]
      if atom.coords[2] < 0:
        atom.coords[2] = 1 + atom.coords[2]
    
    return atoms
  
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

  
  def find_closest_supercell_atom(self, original_atom, super_atom):
    
    print "Finding closest supercell atom:"
    super_atoms = self.supercell_expand([super_atom], exclude_original=True)
    
    shortest_distance = self.calculate_distance(original_atom, super_atom)
    closest_atom = super_atom
    
    for super_atom in super_atoms:
      distance = self.calculate_distance(original_atom, super_atom)
      print distance
      if distance < shortest_distance:
        shortest_distance = distance
        closest_atom = super_atom
    
    return closest_atom

  def supercell_expand(self, atoms, exclude_original = False): 
    
    #print "Starting supercell expansion (excluding original: %s)" % exclude_original
    atoms = self.move_all_atoms_into_cell(atoms)
    new_atoms = []
    cells = [(-1, -1, -1), (-1, -1, 0), (-1, -1, 1),
             (-1,  0, -1), (-1,  0, 0), (-1,  0, 1),
             (-1,  1, -1), (-1,  1, 0), (-1,  1, 1),
             ( 0, -1, -1), ( 0, -1, 0), ( 0, -1, 1),
             ( 0,  0, -1), ( 0,  0, 0), ( 0,  0, 1),
             ( 0,  1, -1), ( 0,  1, 0), ( 0,  1, 1),
             ( 1, -1, -1), ( 1, -1, 0), ( 1, -1, 1),
             ( 1,  0, -1), ( 1,  0, 0), ( 1,  0, 1),
             ( 1,  1, -1), ( 1,  1, 0), ( 1,  1, 1)]
    cells = [np.array([[x],[y],[z]]) for x,y,z in cells]
    
    if exclude_original:
      cells.pop(13) # removes 000 from list
    
    i = 0
    for atom in atoms:
      for thing in cells:
        new_atom = copy.deepcopy(atom)
        #new_atom.label = new_atom.label.split('[', 1)[0]
        #new_atom.label = new_atom.label + "[" +str(thing[0][0])+str(thing[1][0]) + str(thing[2][0]) + "]"
        new_atom.coords += thing
        new_atoms.append(new_atom)
        #print "Old coordinates: (%.2f,%.2f,%.2f). New coordinates: (%.2f,%.2f,%.2f)" % (atom.coords[0][0], atom.coords[1][0], atom.coords[2][0], new_atom.coords[0][0], new_atom.coords[1][0], new_atom.coords[2][0])
        i += 1
            
    #print "Completed expansion to 3 x 3 supercell (%s new atoms:)" % i
    #print ["%s%s" % (atom.symbol, atom.label) for atom in new_atoms]
    return new_atoms
  
  def find_largest_distances(self, atom_list, number_of_clusters):
    
    # For a list of atoms, this find the n atoms with the highest combined distance between them. These atoms are chosen as starting
    # points for finding clusters (n = number of clusters)
    
    if number_of_clusters == 0:
      raise Exception("No Q-peak clusters found")
    
    if number_of_clusters == 1:
      return (0,) # the comma is necessary to make sure that the value is returned as a tuple and not an integer
    
    else:
      distance_matrix = self.make_distance_matrix(atom_list)#)
  
      # This creates a list of all possible combinations in a given range of indexes, 
      # e.g. [0,1,2,3] gives [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)] for groups of two
  
      lst = [i for i in range(len(atom_list))]
      my_combinations = list(combinations(lst,number_of_clusters))
      #This figures out which combination created above has the largest total distance between the peaks
      max_distance = 0
      for n, myset in enumerate(my_combinations):
        distance = 0
        pairs = list(combinations(myset,2))
        #print pairs
        for pair in pairs:
          distance += distance_matrix[pair]
          if distance > max_distance:
            max_distance = distance
            index = myset
            #print "updating_max_distance", max_distance
  
      print "The %s Q-peaks with the highest combined total distance (%.2f) are:" % (number_of_clusters, max_distance)
      for n in range(number_of_clusters):
        print "%s%s" % (atom_list[index[n]].symbol, atom_list[index[n]].label)

      return index
  
  def find_largest_distance_between_atoms(self, atom_list):

    # This creates a list of all possible combinations in a given range of indexes, 
    # e.g. [0,1,2,3] gives [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)] for groups of two

    atom_pairs = list(combinations(atom_list,2))
    #This figures out which combination created above has the largest distance between the peaks
    max_distance = 0
    for atom_pair in atom_pairs:
      distance = self.calculate_distance(atom_pair[0], atom_pair[1])
      if distance > max_distance:
        max_distance = distance

    return max_distance
  
  def find_smallest_distance_between_atoms(self, atom_list):

    # This creates a list of all possible combinations in a given range of indexes, 
    # e.g. [0,1,2,3] gives [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)] for groups of two

    atom_pairs = list(combinations(atom_list,2))
    #This figures out which combination created above has the smallest distance between the peaks
    min_distance = 1000
    for atom_pair in atom_pairs:
      distance = self.calculate_distance(atom_pair[0], atom_pair[1])
      if distance < min_distance:
        min_distance = distance

    return min_distance
  
  def merge_clusters(self, atom_list1, atom_list2, min_distance):

    expanded_peaks1 = self.expand_atoms_with_symmetry(atom_list1)
    superexpanded_peaks1 = self.supercell_expand(expanded_peaks1)
    expanded_peaks2 = self.expand_atoms_with_symmetry(atom_list2)
    superexpanded_peaks2 = self.supercell_expand(expanded_peaks2)
    
    
    for atom1 in superexpanded_peaks1:
      for atom2 in superexpanded_peaks2:
        distance = self.calculate_distance(atom1, atom2)
        if distance < min_distance:
          print "%s%s is only %.2f from %s%s so number of clusters should be reduced" % (atom1.symbol, atom1.name, distance, atom2.symbol, atom2.name)
          OV.Refresh()
          return True
    
    else:
      return False
    
  def calculate_centroid(self, atom_list):

    number_of_atoms = len(atom_list)
    
    if len(atom_list) == 0:
      print "Cannot calculate centroid for no atoms"
      return None
    
    try:
      x_centroid = 0
      y_centroid = 0
      z_centroid = 0
    
      for atom in atom_list:
        x_centroid += atom.coords[0]
        y_centroid += atom.coords[1]
        z_centroid += atom.coords[2]
      x_centroid = x_centroid/number_of_atoms
      y_centroid = y_centroid/number_of_atoms
      z_centroid = z_centroid/number_of_atoms
      
      centroid = np.array([x_centroid,y_centroid,z_centroid])      
      
      return centroid
    
    except:
      print "There was a problem calculating the centroid"
      return None
    
  def place_in_void(self, peak, host_atoms):
    if peak.peak_height > minimum_peak_height:
      place_in_void = True
      for atom in host_atoms:
        length = self.calculate_distance(atom, peak)
        #print "%s%s - %s%s: %.2f" %(peak.symbol, peak.label, atom.symbol, atom.label, length)
        if length < short_contact:
          print "%s%s has been rejected because it is only %.2f from %s%s" % (peak.symbol, peak.label, length, atom.symbol, atom.label)
          return False
      return True  
  
  def compress_Q_peaks(self, peaks):
    
    # This method takes a list of peaks  and arranges them so they are as close to each other as possible, taking into
    # account symmetry generated atoms. 
    
    compressed_q_peaks = []
    compressed_q_peaks.append(peaks[0]) # add the first q-peak to the compressed peaks      
    
    for i in range(1,len(peaks)):
      print "Finding ideal position for %s%s" % (peaks[i].symbol, peaks[i].label)
      expanded_Q_peaks = self.expand_atoms_with_symmetry([peaks[i]])
      supercell_Q_peaks = self.supercell_expand(expanded_Q_peaks)
      #print i, len(supercell_Q_peaks)
      minimum = 1000
      for peak in supercell_Q_peaks: # run through all n * 27 symmetry equivalent peaks where n is number of q-peaks and 27 is 3^3 (i.e. a 3*3*3 supercell)
        #print len(compressed_q_peaks)
        for compressed_peak in compressed_q_peaks: # check all peaks in the original
          length = self.calculate_distance(peak, compressed_peak)
          #print "length %.2f, %.2f, %.2f, %.2f: %.2f, %.2f, %.2f" % (length, compressed_peak.coords[0][0], compressed_peak.coords[1][0], compressed_peak.coords[2][0], peak.coords[0][0], peak.coords[1][0], peak.coords[2][0] )
          if length < minimum:
            minimum = length
            #print "(updating minimum %s for peak %s%s to %s%s)" % (minimum, peak.symbol, peak.label, compressed_peak.symbol, compressed_peak.label) 
            nearest_peak = peak
      compressed_q_peaks.append(nearest_peak)
      print "Added %s%s at a minimum distance of %.2f" % (nearest_peak.symbol, nearest_peak.label, minimum)
    return compressed_q_peaks  
  
  def find_q_peak_clusters(self, peaks, number_of_clusters):

    # This method searches for clusters in the q-peaks.
    
    if len(peaks) <= number_of_clusters:
      number_of_clusters = len(peaks)
    print "Number of clusters to find: %s" % number_of_clusters # Calculated from the total number of molecules that will fit in the void.

    indices = self.find_largest_distances(peaks, number_of_clusters)
    
    print "Clustering..."
    OV.Refresh()
    
    assign1 = []
    k = 1

    for i in range(25): # number of iterations

      if k == 1: # simply select the first n q-peaks 
        centroids = []
        label = 1
        for index in indices:
          centroid = solvator_atom.Solvator_atom("X", label, 1, peaks[index].coords) # define centroids as atoms so easy to expand, etc.
          label +=1
          centroids.append(centroid)

      else: #calculate centroid for each cluster according to the list in assign1

        cluster_peak_lists = []
        for c in range(number_of_clusters):
          peak_list = []
          for i, index in enumerate(assign2):
            if index == c:
              #print "Adding %s%s to Cluster %s" % (compressed_q_peaks[i].symbol, compressed_q_peaks[i].label, c+1)
              peak_list.append(peaks[i])
          cluster_peak_lists.append(peak_list)          
        
        centroids = []
        for peak_list in cluster_peak_lists:
          centroid_coords = self.calculate_centroid(peak_list)
          centroid = solvator_atom.Solvator_atom("X", i, 1, centroid_coords, coord_type='frac')
          centroids.append(centroid)

      # Calculate the distance of each peak to each centroid and decide on the smaller distance. 

      assign2 = [] # a list of indices assigning each entry in compressed_q_peak to a centroid.
      for peak in peaks:
        minimum = 10000000
        for i in range(number_of_clusters):
          distance = self.calculate_distance(centroids[i], peak)
          if distance < minimum:
            minimum = distance
            index = i
        # print "minimum distance %.2f, index %s" % (minimum, index)
        assign2.append(index)

      if assign1 == assign2:
        print "Reached convergence after %s iterations" % k
        break
      else:
        #print "assign1:", assign1
        #print "assign2:", assign2
        assign1 = assign2
        k +=1
    
    return cluster_peak_lists

  def get_peaks_in_voids(self, q_peaks, host_atoms):
    
    # This method sorts q-peaks that may be part of the structure, from those that are in a void, using some user-defined
    # criteria of minimum peak height, and the short contact (the distance from which a q-peak must lie from the main structure)
    
    peaks_in_void = []
    
    print "Number of Q-peaks to analyse: %s" % len(q_peaks)

    if len(q_peaks) == 0:
      raise Exception ("There are no Q-peaks to cluster.")
    q_peaks = self.move_all_atoms_into_cell(q_peaks)

    for peak in q_peaks:
      if self.place_in_void(peak, host_atoms):
        peaks_in_void.append(peak)

    print "Number of Q-peaks to cluster: %s" % len(peaks_in_void)
    if len(peaks_in_void) == 0:
      raise Exception ("There are no Q-peaks in the void.")
    
    return peaks_in_void  
  
  def get_superclusters(self, host_atoms, q_peaks, number_of_clusters, guests):
    
    expanded_host_atoms = self.supercell_expand(host_atoms)
    peaks_in_void = self.get_peaks_in_voids(q_peaks, expanded_host_atoms)
    compressed_q_peaks = self.compress_Q_peaks(peaks_in_void)
    
    if number_of_clusters == 1:
      cluster_peak_lists = [compressed_q_peaks]
    else:
    
      too_many_clusters = True
    
      while too_many_clusters == True:
        
        print "Looking for %s clusters" % number_of_clusters
        cluster_peak_lists = self.find_q_peak_clusters(compressed_q_peaks, number_of_clusters)
          
        too_many_clusters = False
        
        cluster_pairs = combinations(cluster_peak_lists,2)
        
        for pair in cluster_pairs:
          if self.merge_clusters(pair[0], pair[1], short_contact):
            too_many_clusters = True
            number_of_clusters -=1
            break
          else:
            continue
        
      
    superclusters = []
    
    for number, peak_list in enumerate(cluster_peak_lists):
      print "\nCluster %s (of %s)" % (number + 1, len(cluster_peak_lists))
      print "Number of independent peaks in this cluster: %s" % len(peak_list)
      original_peaks = len(peak_list)
      supercluster = self.grow_atoms(peak_list)
      print "Number of peaks in this supercluster: %s" % len(peak_list)
      multiplicity = np.ceil(float(len(supercluster))/original_peaks)

      centroid_coords = self.calculate_centroid(supercluster)
      centroid = solvator_atom.Solvator_atom("X", number, 1, centroid_coords)
        
      total_e_density = 0
      print "Finding cluster centroids (taking into account symmetry-generated atoms):"
      for peak in supercluster: # This is a very crude method of estimating the electron density for a particular cluster. Better to use density maps. Ask Horst.
        print "%s%s [%.2f, %.4f, %.4f]" % (peak.symbol, peak.label, peak.coords[0][0], peak.coords[1][0], peak.coords[2][0])
        total_e_density += peak.peak_height
        #print "%s%s %.2f, %.2f, %.2f, (%.1f e)" % (peak.symbol, peak.label, peak.coords[0], peak.coords[1], peak.coords[2], peak.peak_height)
      print "Centroid %s (of group of %s peaks): (%.2f, %.2f, %.2f). Total electron density: %.1f" % (number+1, len(supercluster), centroid.coords[0], centroid.coords[1], centroid.coords[2], total_e_density)        
      
      new_cluster = q_cluster.Q_cluster(supercluster, total_e_density, centroid)
      new_cluster.multiplicity = multiplicity
      print "Multiplicity of this cluster: %s"  % new_cluster.multiplicity
      superclusters.append(new_cluster)
    
      number +=1
    
    return superclusters
  
  def grow_atoms(self, list_of_atoms):
    """""
    for each cluster, grab nearby peaks - call it a supercluster.
    Then we have a list of superclusters
    """""
    expanded_Q_peaks = self.expand_atoms_with_symmetry(list_of_atoms)
    supercell_Q_peaks = self.supercell_expand(expanded_Q_peaks, exclude_original=False)
    supercluster = list_of_atoms

    duplicates = []
    while len(duplicates) < len(supercluster):
      for i, super_peak in enumerate(supercell_Q_peaks):
        for peak in supercluster:
          if np.array_equal(peak.coords, super_peak.coords):
            duplicates.append(i)
            
    duplicates.sort(reverse=True)
    for d in duplicates:
      supercell_Q_peaks.pop(d)   
      
    #print "Have removed duplicates"
    #print "Number of supercell_Q_peaks: %s" % len(supercell_Q_peaks)
    
    number_of_new_peaks_added = len(supercluster)
    
    while number_of_new_peaks_added > 0 and (len(supercluster) < 50):
      
      length_of_supercluster_at_start = len(supercluster)
      old_number_of_peaks_to_be_added = 0
      indices = []
      
      for peak in supercluster:
        for i, super_peak in enumerate(supercell_Q_peaks):
          distance = self.calculate_distance(peak, super_peak)
          if (distance < OV.GetParam("solvator.solve.short_contact")) and (i not in indices):
            indices.append(i)
            print "%s: Distance %.2f. Peak coords (%s%s): %.2f, %.2f, %.2f Superpeak coords (%s%s): %.2f, %.2f, %.2f" % (i, distance, peak.symbol, peak.label, peak.coords[0][0], peak.coords[1][0], peak.coords[2][0], super_peak.symbol, super_peak.label, super_peak.coords[0][0], super_peak.coords[1][0], super_peak.coords[2][0])

      indices.sort(reverse = True)
      #print indices
      
      if len(indices) > 0:
        print "Found %s more peaks nearby" % (len(indices))
      
      for i in indices:
        peak_to_add = supercell_Q_peaks.pop(i)
        #print "Length of supercell_Q_peaks: %s" % len(supercell_Q_peaks)
        supercluster.append(peak_to_add)
      
      length_of_supercluster_at_end = len(supercluster)
      #print "Length of supercluster after adding new peaks: %s" % length_of_supercluster_at_end
      number_of_new_peaks_added = length_of_supercluster_at_end - length_of_supercluster_at_start
      OV.Refresh()
  
    return supercluster

  
  def add_list_of_molecules_to_olex_model(self, molecule_list, multiplicity, add_restraints = True, add_hydrogens = False):

    print "Using multiplicity of %s" % multiplicity
    if multiplicity == 1:
      sign = " "
    elif multiplicity >= 2:
      sign = "-"
    
    site_occu = 1/multiplicity

    # We can only ever reap the resfile, otherwise, we will be opeining a file were the filename of the res
    # does not match the hkl filename.
    resfile = os.path.join(OV.FilePath(), OV.FileName()) + ".res"
    original_resfile = os.path.join(OV.FilePath(), OV.FileName()) + "_original.res"
    copyfile(original_resfile, resfile)
    olx.Atreap(resfile)
    
    atom_number = self.find_highest_atom_number_in_model()
    solvator_number = 0
    
    if len(molecule_list) > 2:
      SUMP = True
    else:
      SUMP = False
    
    for molecule in molecule_list:
      #print "Adding %s to model (with atom number starting %s)" % (molecule.name, atom_number + 1)
      #molecule.set_occupancy(molecule.occupancy/len(molecule_list))
      #print "Occupancy of this %s: %.2f" % (molecule.name, molecule.occupancy)
    
      for atom in molecule.transformed_atoms:
        atom_number +=1
        solvator_number += 1
        atom.solvator_id = "SOLV" + str(solvator_number)
        atom.olex_id = atom.symbol + str(atom_number)

    number_of_PARTS = 0
    FVAR_line_count = 0
    with open(original_resfile, 'r') as fh:
      for line in fh.readlines():
        if re.match(r'^PART.*', line):
          line = line.strip().split()
          PART_number = abs(int(line[1]))
          if PART_number > number_of_PARTS:
            number_of_PARTS = PART_number
    
    #print "Number of parts: %s" % number_of_PARTS
    
    number_of_PARTS +=1
    
    FVAR_line_number = 0
    with open(original_resfile, 'r') as infile, open(resfile, 'w') as outfile:
      #atom_line = re.compile(r'^([a-zA-Z]+[0-9]+[A-Za-z]*\s+[0-9]+)\s+(-*)([0-9]\.[0-9]+)\s+(-*)([0-9]\.[0-9]+)\s+(-*)([0-9]\.[0-9]+)(.*)')
      for line in infile.readlines():
        #print line
        if re.match(r'^FVAR.*', line):
          line = line.strip().split()
          number_of_FVAR = len(line)
          for item in line:
            outfile.write("%s " % item)
          if len(molecule_list) == 2:
            outfile.write("%.4f " % molecule_list[0].occupancy)
          elif len(molecule_list) >= 3:
            for molecule in molecule_list:
              outfile.write("%.4f " % molecule.occupancy)
          outfile.write("\n")
          if SUMP:
            outfile.write("SUMP 1.0 0.001 ")
            for i in range(number_of_FVAR, number_of_FVAR + len(molecule_list)):
              outfile.write("1 %s " % i)
            outfile.write("\n")
        #elif re.match(atom_line, line):
          #outfile.write("%s  %s1%.6f  %s1%.6f  %s1%.6f%s\n" % (atom_line.match(line).group(1), atom_line.match(line).group(2),float(atom_line.match(line).group(3)), atom_line.match(line).group(4), float(atom_line.match(line).group(5)), atom_line.match(line).group(6), float(atom_line.match(line).group(7)), atom_line.match(line).group(8)))
        elif re.match(r'HKLF.*', line):
          #number_of_FVAR -=1
          
          if len(molecule_list) == 1:
            outfile.write("PART %s%s\n" % (sign, number_of_PARTS))
            for atom in molecule.transformed_atoms:
              outfile.write(("%s %s %.6f %.6f %.6f 1%.6f %.6f\n") % (atom.olex_id, atom.sfac, atom.coords[0], atom.coords[1], atom.coords[2],  site_occu, atom.uiso))#, number_of_FVAR))
            outfile.write("PART 0\n")
          
          elif len(molecule_list) == 2:
            outfile.write("PART %s%s %s%.6f\n" % (sign, number_of_PARTS, number_of_FVAR, site_occu))
            for atom in molecule_list[0].transformed_atoms:
              outfile.write(("%s %s %.6f %.6f %.6f %s%.6f %.6f\n") % (atom.olex_id, atom.sfac, atom.coords[0], atom.coords[1], atom.coords[2], number_of_FVAR,  site_occu, atom.uiso))#, number_of_FVAR))
            number_of_PARTS +=1
            outfile.write("PART %s%s -%s%.6f\n" % (sign, number_of_PARTS, number_of_FVAR, site_occu))
            for atom in molecule_list[1].transformed_atoms:
              outfile.write(("%s %s %.6f %.6f %.6f -%s%.6f %.6f\n") % (atom.olex_id, atom.sfac, atom.coords[0], atom.coords[1], atom.coords[2], number_of_FVAR,  site_occu, atom.uiso))#, number_of_FVAR))
            outfile.write("PART 0\n")
              
          elif len(molecule_list) >= 3:
            for molecule in molecule_list:
              outfile.write("PART %s%s %s%.6f\n" % (sign, number_of_PARTS, number_of_FVAR, site_occu))
              for atom in molecule.transformed_atoms:
                outfile.write(("%s %s %.6f %.6f %.6f %s%.6f %.6f\n") % (atom.olex_id, atom.sfac, atom.coords[0], atom.coords[1], atom.coords[2], number_of_FVAR,  site_occu, atom.uiso))#, number_of_FVAR))
              number_of_FVAR += 1
              number_of_PARTS +=1
            outfile.write("PART 0\n")
            
          outfile.write(line)
        else:
          outfile.write(line)

    olx.Atreap(resfile)
    #print "Number of atoms after writing resfile: %s" % len(olexex.OlexRefinementModel()._atoms)

    OV.cmd('fuse')
    if add_restraints:
      for molecule in molecule_list:
        if molecule.original_restraints == []:
          print "No restraints for this %s molecule" % molecule.name
        else:
        # Find and replace solvator_id with olex_id in each restraint. This has to be done in stages
        # just in case the guest atom's original name overlaps with any of the names of the host atoms.
          molecule.restraints = copy.deepcopy(molecule.original_restraints)
          for n in range(len(molecule.original_restraints)):
            for atom in molecule.transformed_atoms:
              print "replacing %s with %s" % (atom.name, atom.solvator_id)
              molecule.restraints[n] = molecule.restraints[n].replace(atom.name, atom.solvator_id)
            for atom in molecule.transformed_atoms:
              print "replacing %s with %s" % (atom.solvator_id, atom.olex_id)
              molecule.restraints[n] = molecule.restraints[n].replace(atom.solvator_id, atom.olex_id)
  
        for restraint in molecule.restraints:
          OV.cmd(restraint)

    R_current = float(str(olx.CalcR()).split(',')[0])

    return R_current

  def get_starting_positions(self, supercluster, guest_molecules, host_atoms, R_host):
    
    resfile = os.path.join(OV.FilePath(), OV.FileName()) + ".res"
    original_resfile = os.path.join(OV.FilePath(), OV.FileName()) + "_original.res"
    
    lowest_R_so_far = 100
    q_peaks = supercluster.peak_list
    solutions = []
    print "Atoms in this supercluster:"
    for atom in q_peaks:
      print "%s%s" % (atom.symbol, atom.label)
    
    for molecule in guest_molecules:
      
      #Sort the qpeaks list such that we can skip over starting duplicates
      
      q_peaks = sorted(q_peaks, key=lambda x: x.label, reverse=False)
      print "Sorted peaks:"
      for atom in q_peaks:
        print "%s%s" % (atom.symbol, atom.label)
      
      number= len(q_peaks)
      
      # Now what we do depends on how many atoms are in the molecule. 
      
      if len(molecule.transformed_atoms) == 1:
        first_guest_atom = molecule.transformed_atoms[0]
        
        for i in range(number):
          #print "i:%s (%s%s)" % (i, q_peaks[i].symbol, q_peaks[i].label)
          if (i > 0) and (q_peaks[i].label == q_peaks[i-1].label): # just skip over atoms with the same label
            continue
    
          translation = q_peaks[i].coords - first_guest_atom.coords
          molecule.translate(translation)
          [a],[b],[c], = np.dot(FRAC_TO_CART, q_peaks[i].coords)
          
          if self.is_chemically_sensible(molecule.transformed_atoms,
                                         host_atoms): 
            oriented_molecule = copy.deepcopy(molecule)
            solutions.append(oriented_molecule)        
      
      if len(molecule.transformed_atoms) == 2:
        first_guest_atom = molecule.transformed_atoms[0]
        second_guest_atom = molecule.transformed_atoms[1]
        guest_distance1_2  = self.calculate_distance(first_guest_atom, second_guest_atom)
        
        print "Bond length in guest: %.2f" % (guest_distance1_2)
        
        for i in range(number):
          #print "i:%s (%s%s)" % (i, q_peaks[i].symbol, q_peaks[i].label)
          if (i > 0) and (q_peaks[i].label == q_peaks[i-1].label): # just skip over atoms with the same label
            continue
    
          translation = q_peaks[i].coords - first_guest_atom.coords
          molecule.translate(translation)

          
          for j in [x for x in range(number) if x != i]:
            q_peak_distance = self.calculate_distance(q_peaks[i], q_peaks[j])
            if np.isclose(q_peak_distance, guest_distance1_2, atol = 0.4):
              #print "Within range: j:%s (%s%s) is %.3f Angstroms from i(%s) (%s%s) (should be %.2f)" % (j, q_peaks[j].symbol, q_peaks[j].label, q_peak_distance, i, q_peaks[i].symbol, q_peaks[i].label, guest_distance1_2)
              second_q_peak = q_peaks[j]
              q_peak_vector1_2 = self.get_vector(q_peaks[j], q_peaks[i])
              guest_vector1_2 = self.get_vector(second_guest_atom, first_guest_atom)
              
              rotation_angle = self.calculate_angle_between_atoms(q_peaks[j], q_peaks[i], second_guest_atom)
    
              u,v,w = np.cross(q_peak_vector1_2, guest_vector1_2)
              
              for atom in molecule.transformed_atoms:
                [x],[y],[z] = np.dot(FRAC_TO_CART, atom.coords )
                cart_coords = self.calculate_rotation(x, y, z, a, b, c, u, v, w, -rotation_angle)
                atom.coords = np.dot(CART_TO_FRAC, cart_coords)
              
              if self.is_chemically_sensible(molecule.transformed_atoms,
                                             host_atoms): 
                oriented_molecule = copy.deepcopy(molecule)
                solutions.append(oriented_molecule)
      
      if len(molecule.transformed_atoms) >= 3:
      
        first_guest_atom = molecule.transformed_atoms[0]
        second_guest_atom = molecule.transformed_atoms[1]
        third_guest_atom = molecule.transformed_atoms[2]
    
        guest_distance1_2  = self.calculate_distance(first_guest_atom, second_guest_atom)
        guest_distance1_3 = self.calculate_distance(first_guest_atom, third_guest_atom)
        guest_distance2_3 = self.calculate_distance(second_guest_atom, third_guest_atom)
        guest_angle = self.calculate_angle_between_atoms(third_guest_atom, second_guest_atom, first_guest_atom)
        
        print "Guest distances: %.2f, %.2f, %.2f" % (guest_distance1_2, guest_distance1_3, guest_distance2_3)
        
        """"
        Find nearby q-peaks
        """
    
        for i in range(number):
          #print "i:%s (%s%s)" % (i, q_peaks[i].symbol, q_peaks[i].label)
          if (i > 0) and (q_peaks[i].label == q_peaks[i-1].label): # just skip over atoms with the same label
            continue
    
          translation = q_peaks[i].coords - first_guest_atom.coords
          molecule.translate(translation)
          [a],[b],[c], = np.dot(FRAC_TO_CART, q_peaks[i].coords)
          
          for j in [x for x in range(number) if x != i]:
            q_peak_distance = self.calculate_distance(q_peaks[i], q_peaks[j])
            if np.isclose(q_peak_distance, guest_distance1_2, atol = 0.4):
              #print "Within range: j:%s (%s%s) is %.3f Angstroms from i(%s) (%s%s) (should be %.2f)" % (j, q_peaks[j].symbol, q_peaks[j].label, q_peak_distance, i, q_peaks[i].symbol, q_peaks[i].label, guest_distance1_2)
              second_q_peak = q_peaks[j]
              q_peak_vector1_2 = self.get_vector(q_peaks[j], q_peaks[i])
              guest_vector1_2 = self.get_vector(second_guest_atom, first_guest_atom)
              
              rotation_angle = self.calculate_angle_between_atoms(q_peaks[j], q_peaks[i], second_guest_atom)
    
              u,v,w = np.cross(q_peak_vector1_2, guest_vector1_2)
              
              for atom in molecule.transformed_atoms:
                [x],[y],[z] = np.dot(FRAC_TO_CART, atom.coords )
                cart_coords = self.calculate_rotation(x, y, z, a, b, c, u, v, w, -rotation_angle)
                atom.coords = np.dot(CART_TO_FRAC, cart_coords)
            
              for k in [y for y in range(number) if (y != i and y!= j)]:
                q_peak_distance = self.calculate_distance(q_peaks[j], q_peaks[k])
                q_peak_angle = self.calculate_angle_between_atoms(q_peaks[k], q_peaks[j], q_peaks[i])
    
                if (np.isclose(q_peak_distance, guest_distance2_3, atol = 0.4)) and np.isclose(q_peak_angle,guest_angle,atol = 0.52):
                  #print "Within range: k:%s (%s%s) is %.3f Angstroms from j:%s (%s%s) (should be %.2f)" % (k, q_peaks[k].symbol, q_peaks[k].label, q_peak_distance, j, q_peaks[j].symbol, q_peaks[j].label, guest_distance2_3)
                  third_q_peak = q_peaks[k]
                  
                  normal1 = self.calculate_normal(q_peaks[i], q_peaks[j], q_peaks[k])
                  normal2 = self.calculate_normal(q_peaks[i], q_peaks[j], third_guest_atom)
                  rotation_angle = self.calculate_angle_between_planes(normal1, normal2)
                  
                  u,v,w = np.cross(normal1, normal2)
                  
                  for atom in molecule.transformed_atoms:
                    [x],[y],[z] = np.dot(FRAC_TO_CART, atom.coords )
                    cart_coords = self.calculate_rotation(x, y, z, a, b, c, u, v, w, -rotation_angle)
                    atom.coords = np.dot(CART_TO_FRAC, cart_coords)
  
                  if self.is_chemically_sensible(molecule.transformed_atoms,
                                                 host_atoms): 
                    oriented_molecule = copy.deepcopy(molecule)
                    solutions.append(oriented_molecule)
      
      return solutions

  def restore_original_resfile(self):
    
    resfile = os.path.join(OV.FilePath(), OV.FileName()) + ".res"
    restore_resfile = os.path.join(OV.FilePath(), OV.FileName()) + "_restore.res"
    
    copyfile(restore_resfile, resfile)
    
    olx.Atreap(resfile)
    
  
  def calculate_solvent_numbers(self, total_volume, total_electrons, guests_used, Z):
    model = []
    solvent_matrix = np.zeros((2, len(guests_used)))
    solvent_max_numbers = []
    solvent_min_numbers = []
    print "Total void volume %.1f A3, total missing electrons %.1f" % (total_volume, total_electrons)
    total_volume = total_volume/Z
    total_electrons = total_electrons/Z
    
    for solvent in guests_used:
      #print solvent.name, guests_used.index(solvent)
      solvent_matrix[0,guests_used.index(solvent)] = solvent.electrons
      #solvent_expanded[1,guests_used.index(solvent)] = solvent.volume
      max_solvent_by_volume = (total_volume/solvent.volume)
      max_solvent_by_electrons = (total_electrons/solvent.electrons)
      print solvent.name, ":"
      print "by volume, %.2f:" % max_solvent_by_volume
      print "by electrons, %.2f:" % max_solvent_by_electrons
      solvent_max_numbers.append(min(max_solvent_by_volume,max_solvent_by_electrons))
      solvent_min_numbers.append(0)

    #print "Max number of solvent molecules:", solvent_max_numbers
    #print "solvent matrix", solvent_matrix
    desired = np.array([total_electrons,total_volume])
    #print "desired solution:", desired

    #solution, residuals, rank, s = np.linalg.lstsq(solvent_matrix, desired,rcond=None)
    #Use linear least squares to find the best number of a particular solvent to fit in the void
    solution = scipy.optimize.lsq_linear(solvent_matrix,desired,bounds=(solvent_min_numbers,solvent_max_numbers)).x

    #Populate the model list with tuples containing the solvent name and number expected, e.g. ((benzene, 3), (toluene, 2))
    for i, solvent in enumerate(guests_used):
      number_of_molecules = solution[i]
      if not np.isclose(number_of_molecules, 0.0, atol = 0.1): # ignore solvent numbers less than 10% calculated occupancy
        print "Proposed number of", solvent.name, "molecules:", number_of_molecules
        model.append((solvent, number_of_molecules))
      else:
        print "(Probably no", solvent.name, ")"

    return model
  
  def get_Q_peaks(self):
    
    resfile = os.path.join(OV.FilePath(), OV.FileName()) + ".res"
    q_peak_pattern = re.compile(r'^Q[0-9]+\s+[0-9]+(\s+-*[0-9]+(\.[0-9]+)*){3}.*')
    q_peak_lines = []
    q_peaks = []
    
    with open(resfile, 'r') as fh:
        for line in fh.readlines():
            if re.match(q_peak_pattern, line):
              q_peak_lines.append(line)
    
    for i, line in enumerate(q_peak_lines):
        line = line.split()
        atom_name = line[0].strip()
        atom_symbol = re.sub(r'([A-Z]+[a-z]*)\d*.*', r'\1', atom_name).lower().capitalize()
        atom_label = re.sub(r'[A-Z]+[a-z]*(\d*.*)', r'\1', atom_name)
        atom_sfac = int(line[1].strip())
        atom_x = float(line[2].strip())
        atom_y = float(line[3].strip())
        atom_z = float(line[4].strip())
        atom_coords = np.array([[atom_x],[atom_y],[atom_z]])
       
        this_q_peak = solvator_atom.Solvator_atom(atom_symbol, atom_label, atom_sfac, atom_coords, 1, 1, 0.03, None, 'frac')
        this_q_peak.peak_height = float(line[7].strip())
        q_peaks.append(this_q_peak)
    
    print len(q_peaks)
    
    return q_peaks

  def solve(self):
    
    time_start = time.time()
    
    print "\n\nStarting SOLVATOR...\n"
    

    OV.Refresh()
    
    resfile = os.path.join(OV.FilePath(), OV.FileName()) + ".res"
    original_resfile = os.path.join(OV.FilePath(), OV.FileName()) + "_original.res"
    restore_file = os.path.join(OV.FilePath(), OV.FileName()) + "_restore.res"

    #This next section reads the res file via "read_res" and assigns "global-ish" variables, which are used in various methods called by the SOLVATOR.
    CRYSTAL_SYSTEM, \
    OSF,\
    FVAR_FLOATS, \
    SFAC_BY_SYMBOL = self.read_res(resfile)

    time_data_start = time.time()
    
    guests_used = []
    
    solvent_directory = os.path.join(p_path, "solvents/")
    
    print p_path
    print solvent_directory
    
    file_list = os.listdir(solvent_directory)
    
    search_pattern = re.compile(r'(.*).guest')
    
    available_solvents = []
    
    for solvent_file in file_list:
      m = search_pattern.match(solvent_file)
      if m:
        solvent_name = m.group(1)
        available_solvents.append(solvent_name)
    
    print available_solvents
    
    for solvent_name in available_solvents:
      print solvent_name
      used = OV.GetParam('solvator.solvent' + '.' + str(solvent_name))
      if used == True:
        print "used"
        thisguest = guest.Guest(solvent_name, 0.5, p_path, SFAC_BY_SYMBOL, CART_TO_FRAC)
        guests_used.append(thisguest)
    OV.Refresh()
       
    if not guests_used:
      raise Exception("Please select some solvent molecules")
    else:
      print "These are the solvents suggested by you:\n"
      for item in guests_used:
        print item.name
      print "\n"    

    mc_trials = OV.GetParam("solvator.solve.trials")
    mc_fineness = OV.GetParam("solvator.solve.fineness")
    
    global short_contact
    short_contact = OV.GetParam("solvator.solve.short_contact")
    global minimum_peak_height
    minimum_peak_height = OV.GetParam("solvator.solve.minimum_peak_height")
    cutoff = OV.GetParam("solvator.solve.cutoff")
    use_host_symmetry = OV.GetParam("solvator.solve.use_host_symmetry")
    max_number_molecules_in_solution = OV.GetParam("solvator.solve.max_molecules_in_solution")
    OV.SetParam('solvator.solve.stop_mc', False)
    
    max_tries_chemical_sense = 1000

    time_data_end = time.time()
    
    print "Time taken for data manipulations: %.2f s" % (time_data_end - time_data_start)


    time_setup_start = time.time()

    """
    This section runs through a number of refinements and makes sure that
    there are enough peaks for the SOLVATOR algorithm to work with
    """
    
    OV.SetParam('snum.refinement.use_solvent_mask',False)
    OV.cmd('refine 10 5')
    OV.cmd('kill $Q')

    
    print 'Using CalcSolv to extract voids'
    OV.SetParam('snum.refinement.use_solvent_mask',True)
    OV.SetParam('snum.masks.update_cif', True)
    OV.AddIns('ABIN')
    OV.cmd('refine 10 5')
    OV.SetParam('snum.refinement.recompute_mask_before_refinement',True)
    OV.cmd('refine 10 5') 
    R_masked = float(str(olx.CalcR()).split(',')[0])
    
    copyfile(resfile, original_resfile)
    copyfile(original_resfile, restore_file)
    
    with open(original_resfile, "r") as f:
        lines = f.readlines()
    with open(original_resfile, "w") as f:
        for line in lines:
            if line.strip("\n") != "ABIN":
                f.write(line)
    
    myhost = host.Host(OV.FileName(), FVAR_FLOATS)
    
    #print myhost.average_uiso
    host_atoms = self.expand_atoms_with_symmetry(myhost.uniq_atoms)
    host_atoms = self.move_all_atoms_into_cell(host_atoms)
    
    total_uiso = 0
    total_atoms = 0
    for atom in host_atoms:
      if atom.uiso < 0:
        pass
      else:
        total_uiso += atom.uiso
        total_atoms += 1
  
    average_uiso = total_uiso/total_atoms
    uiso_ratio = OV.GetParam("solvator.solve.uiso_ratio")
    
    print "Number of host atoms in cell:", len(host_atoms)
        
    myhost.get_voids() # We can do this now because we have just run the mask algorithm

    OV.SetParam('snum.refinement.use_solvent_mask',False)
    OV.DelIns('ABIN')
    OV.cmd('refine 10 %s' % OV.GetParam("solvator.solve.q_peaks"))    # need to get lots of Q-peaks in order to calculate sensible void clusters
    
    #R_host = self.calculate_R(host_fcf,reflections)
    R_host = float(str(olx.CalcR()).split(',')[0])
    print "R-factor for host alone: %.2f%%" % (R_host*100)
    OV.Refresh()
    
    time_setup_end = time.time()
    
    cluster_time = time_setup_end - time_setup_start
    
    q_peaks = self.get_Q_peaks()
    print "Total missing electrons:", myhost.total_missing_electrons
    print "Total void volume", myhost.total_void_volume
    
    if myhost.total_missing_electrons < 1:
      raise Exception("No electron density found. Please check that the data are not already masked.")
    Z = float(olx.xf.au.GetZ())
    
    #number_of_uniq_voids = len(myhost.uniq_voids)
    total_number_of_voids = len(myhost.voids)
    #print "Number of unique voids: %s" % number_of_uniq_voids
    #print "Number of total voids:  %s" % total_number_of_voids
    print "Z:", float(olx.xf.au.GetZ())
    
    print "VOIDS:"
    print "Number \t Volume \t Electrons \t Multiplicity"
    total_solvents = 0.0
    suggested_guests = []
    #for void in myhost.voids:
      #print "%s\t %.2f \t %.1f \t %s\n" % (void.number, void.volume, void.electrons, void.multiplicity)
    
      #void.model = self.calculate_solvent_numbers(void.volume, void.electrons, guests_used, Z)
      #print void.model
      #total_molecules_in_void = int(round(sum(x[1] for x in void.model))) 
      #if len(myhost.uniq_voids) == len(myhost.voids):
      #  total_molecules_in_void = total_molecules_in_void/Z
    #total_solvents += total_molecules_in_void
      #for solvent in void.model:
        #if (solvent[1] > 0.1):
          #suggested_guests.append(solvent[0])      
    
    model = self.calculate_solvent_numbers(myhost.total_void_volume, myhost.total_missing_electrons, guests_used, Z)
    for solvent in model:
      if (solvent[1] > 0.1):
        suggested_guests.append(solvent[0])
        total_solvents += solvent[1]
    
    print total_solvents
    total_solvents = int(np.ceil(total_solvents))

    
    suggested_guests = list(set(list(suggested_guests))) # removes duplicates
    print "Suggested guests for your structure:"
    for solvent in suggested_guests:
      print solvent.name
    
    print "Total solvent molecules:", total_solvents
    total_solvents += 2 #(because better to look for too many rather than too few clusters)

    superclusters = self.get_superclusters(host_atoms, q_peaks, total_solvents, suggested_guests)
    superclusters.sort(key=lambda x: x.electron_count, reverse=True)
    
    for c, supercluster in enumerate(superclusters):
  
      print "Finding the best starting positions of molecules to Q-peaks for cluster %s (of %s clusters)" % (c+1, len(superclusters))
  
      
      # This next line produces a list of suitable starting solutions for the molecules. 
      supercluster.starting_positions = self.get_starting_positions(supercluster, suggested_guests, host_atoms, R_host)
      
      
      if supercluster.starting_positions:
        print "Number of starting solutions for cluster %s: %s" % (c+1, len(supercluster.starting_positions))
          
      else:
        print "There were no suitable starting positions for Cluster %s" % (c+1)
        continue # on to next cluster
      
      if supercluster.starting_positions:
        #print supercluster.starting_positions
        
        max_number = min(max_number_molecules_in_solution, len(supercluster.starting_positions))
        best_R_factor = R_host
        for i in range(1,max_number+1):
          my_combinations = []
          if i ==1:
            print supercluster.starting_positions
            my_list = combinations(supercluster.starting_positions, i)
          else:
            print reduced_list
            my_list = combinations(reduced_list, i)
          my_combinations += my_list
        
          copy_combinations = copy.deepcopy(my_combinations)
          print "Number of molecules in these combinations %s" % i
          results = []
          print "R-factor to beat for cluster %s: %.2f" % (c+1, best_R_factor*100)
          for n, combination in enumerate(copy_combinations):
            R_factor = self.add_list_of_molecules_to_olex_model(combination, supercluster.multiplicity, add_restraints = False)
            print "R_factor: %.2f" % (R_factor*100)
            result = [combination, R_factor]
            results.append(result)
          
          results.sort(key = itemgetter(1))
          print results[0][1]
          sorted_combinations = [result[0] for result in results]
          if i == 1:
            reduced_list = []
            number = min(len(sorted_combinations), 20)
            for j in range(0,number):
              reduced_list += sorted_combinations[j]
          #print sorted_combinations[0]
          else:
            print len(sorted_combinations)
          
          n = min(3,len(sorted_combinations))
          for combination in sorted_combinations[0:n]:
            reject_solution = False
            self.add_list_of_molecules_to_olex_model(combination, supercluster.multiplicity)
            OV.cmd('refine 5 5')
            R_factor = float(str(olx.CalcR()).split(',')[0])
        
            for molecule in combination:
              molecule.update_coordinates_from_olex_model()
              for atom in molecule.transformed_atoms:
                if atom.uiso > (uiso_ratio*average_uiso):
                  print "Rejecting because of Uiso of %s%s. Average Uiso of structure: %.3f, this atom: %.3f" % (atom.symbol, atom.label, average_uiso, atom.uiso)
                  reject_solution = True
                  break
            
            if not reject_solution:
              print "Not rejected:", R_factor
              if R_factor < best_R_factor:# - 0.05*(len(combination)-1)): # A penalty to favour solutions with fewer molecules.
                print "We've got an improved combination here. R factor: %.2f" % (R_factor*100)
                best_R_factor = R_factor
                supercluster.solution = copy.deepcopy(combination)

      if supercluster.solution:
        print "The best solution for cluster %s has %s molecules" % (c+1, len(supercluster.solution))
        self.add_list_of_molecules_to_olex_model(supercluster.solution, multiplicity = supercluster.multiplicity)
        for molecule in supercluster.solution:
          OV.cmd('refine 5 5')
          if molecule.original_hfix != []:
            molecule.add_hfix()
            OV.cmd('refine 5 5')
          best_R_factor = float(str(olx.CalcR()).split(',')[0])
        copyfile(resfile,original_resfile)
      else:
        print "No solution for cluster %s" % (c+1)
        copyfile(original_resfile, resfile)
        olx.Atreap(resfile)
    
    R_factor = float(str(olx.CalcR()).split(',')[0])
    
    if R_host <= R_factor:
      print "The program has not produced a better solution than for the host alone. Sorry about that!\n"
      self.restore_original_resfile()

    time_end = time.time()
    time_taken = time_end - time_start
    
    print ("\n\nEnd of SOLVATOR")
    print "This process took {:.2f}".format(time_taken/60) , "minutes"  


Solvator_instance = Solvator()
