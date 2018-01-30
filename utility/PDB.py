#!/bin/python3 -u

# TODO Clean up and reference atom and residue.py
# To fix...
# Let run with only one or two of the three input file types
# Multithread if you want to waste some time with learning...
# Import vector class (copy into own file)

import math
import subprocess
import struct
from vec3 import vec3

class PDBFile:
  def __init__(self):
    self.chains = dict()
    self.molecules = []
    self.atoms = []
    self.filename = ""
    self.nAtoms = 0
  def CheckMolecules(self, topFile):
    if len(self.chains) < 1:
      if not self.AssembleMolecules():
        return
    if not topFile:
      print("Missing topology file.")
      return
    error = False
    # Else have data we need
    missing = []
    for k,v in self.chains.items():
      for m in v:
        if m.name not in topFile.residues:
          error = True
          if m.name not in missing:
            print("Missing topology data for residue", m.name)
            missing.append(m.name)
          break
        # Have the data
        r = topFile.residues[m.name]
        mol = set([n.name for n in m.atoms])
        res = set([n.name for n in r.atoms.values()])
        a = mol - res
        b = res - mol
        if len(a) > 0 or len(b) > 0:
          print("Molecule",m.id,m.name,"in chain",m.chain,"has problems...")
          error = True
          if len(a) > 0:
            print("\tExtra atoms not found in topology file: ", a)
          if len(b) > 0:
            print("\tAtoms missing from pdb file: ", b)
    if not error:
      print("Molecules in pdb file match topology definitions.")
    else:
      print("Problems with analysis... (See above)")
      print("Possibly patches which this program ignores for now.")
  def CheckBonds(self, topFile, minimum = 1.0):
    if len(self.chains) < 1:
      if not self.AssembleMolecules():
        return
    if not topFile:
#      print("Missing topology file.")
      pass
    print("Checking for bonds less than",minimum)
    print("Atom indices will start at 1.")
    for i in range(len(self.atoms)):
      a = self.atoms[i]
      for j in range(i,len(self.atoms)):
        if i == j:
          continue
        b = self.atoms[j]
        distance = a.Distance(b)
        if distance < minimum:
          print("Atom",a.id,"and",b.id,"are too damn close","{0: <6.3f}".format(distance),"Angstroms")
  def AssembleMolecules(self, verbose = False):
    if(self.nAtoms < 1):
      print("No atoms in PDB file.")
      return
    if verbose:
      print("Assembling molecules from",len(self.atoms),"atoms.")
    for i in self.atoms:
      resname = str(i.resname).strip() # Strip surrounding white space
      chain = i.chain
      rid = int(i.resID)
      if rid < 1:
        print("Error: Atom residue index is out of range...")
        return False
      if chain not in self.chains:
        self.chains[chain] = []
      while rid > len(self.chains[chain]):
        self.chains[chain].append(Molecule())
      self.chains[chain][rid-1].name = resname
      self.chains[chain][rid-1].id = rid
      self.chains[chain][rid-1].chain = chain
      self.chains[chain][rid-1].atoms.append(i)
    # Clean up
    for k,v in self.chains.items():
      index = 0
      while index < len(v):
        if len(v[index].atoms) == 0:
          v.pop(index)
          index -= 1
        index += 1
#      print("Chain",k,"reduced to",len(v),"molecules.")
    return True
  def CutWater(self, radius):
    min_vec = [0.,0.,0.]
    max_vec = [0.,0.,0.]
    while i < len(self.atoms):
      a = self.atoms[i]
      resname = str(a.resname).strip()
      if "TIP" not in resname:
        min_vec[0] = min(a.x, min_vec[0])
        min_vec[1] = min(a.y, min_vec[1])
        min_vec[2] = min(a.z, min_vec[2])
        max_vec[0] = max(a.x, max_vec[0])
        max_vec[1] = max(a.y, max_vec[1])
        max_vec[2] = max(a.z, max_vec[2])
      i += 1
    print("Max (not water):",max_vec)
    print("Min (not water):",min_vec)
    
    nDel = 0
    while i < len(self.atoms):
      a = self.atoms[i]
      resname = str(a.resname).strip()
      if "TIP" in resname:
        if x > min_vec[0] and x < max_vec[0]:
          if y > min_vec[1] and y < max_vec[1]:
            if z > min_vec[2] and z < max_vec[2]:
              # Intruder
              del self.atoms[i]
              nDel += 1
              i -= 1
      i += 1
    print("Deleted",nDel,"interior water atoms.")
    # Now catch stragglers...
    id1 = 0; id2 = 0; id3 = 0
    while i < len(self.atoms):
      a = self.atoms[i]
      resname = str(a.resname).strip()
      if "TIP" in resname:
        if id1 == 0:
          id1 = a.resID
        if i > 0:
          break
    print("Code incomplete...")
    return False
    
  def WritePDBFromChains(self,fn):
    with open(fn, 'w') as f:
      # Write main chains
      for key in sorted(self.chains.keys()):
        if "ION" in key:
          continue
        for mol in self.chains[key]:
          for a in mol.atoms:
            f.write(a.PDBLine())
      # Write ions last
      if "ION" in self.chains.keys():
        for mol in self.chains["ION"]:
          for a in mol.atoms:
            f.write(a.PDBLine())
    f.close()
        
class PDBAtom:
  pass
class Molecule:
  def __init__(self):
    self.name = ""
    self.chain = ""
    self.id = 0
    self.atoms = []
    self.nAtoms = 0
    self.mass = 0.
    self.CoM = vec3()
  def CenterOfMass(self):
    self.mass = 0.
    for a in self.atoms:
      if len(a.element) > 2:
        a.element = a.element[:2].strip()
      if "C" in a.element:
        self.mass += 12.011
      elif "H" in a.element:
        self.mass += 1.0079
      elif "N" in a.element:
        self.mass += 14.007
      elif "O" in a.element:
        self.mass += 15.999
      elif "NA" in a.element:
        self.mass += 22.990
      elif "CL" in a.element:
        self.mass += 35.453
      elif "P" in a.element:
        self.mass += 30.974
      elif "S" in a.element:
        self.mass += 32.065
      elif "K" in a.element:
        self.mass += 39.098
    self.CoM = vec3()
    for a in self.atoms:
      mass = 0.
      if "C" in a.element:
        mass += 12.011
      elif "H" in a.element:
        mass += 1.0079
      elif "N" in a.element:
        mass += 14.007
      elif "O" in a.element:
        mass += 15.999
      elif "NA" in a.element:
        mass += 22.990
      elif "CL" in a.element:
        mass += 35.453
      elif "P" in a.element:
        mass += 30.974
      elif "S" in a.element:
        mass += 32.065
      elif "K" in a.element:
        mass += 39.098
      self.CoM += (vec3([a.x, a.y, a.z])*mass/self.mass)
    return self.CoM

class Residue:
  def __init__(self):
    self.name = ""
    self.charge = 0.
    self.atoms = dict()
  def CheckCharge(self, verbose = False):
    charge = 0.
    for i in self.atoms.values():
      charge += i.charge
    dQ = abs(self.charge - charge)
    if dQ > 1e-4:
      print("Residue",self.name,"has bad topology data.")
      print("\tCharge is off by", dQ)

lastID = 1
class Atom:
  def ReadLine(self, l):
    global lastID
    if l[:4] != "ATOM":
      return False
    try:
      self.id = int(l[4:11])
    except:
      try:
        self.id = int(l[4:11], 16)
      except:
        try:
          if '*' in l[4:11]:
            self.id = lastID
        except:
          print("Error reading ATOM entry:",l)
          print("\tSkipping...")
          return False
    self.name = l[11:16].strip()
    self.resname = l[16:21].strip()
    self.chainID = l[21].strip()
    self.resID = int(l[22:26])
    self.x = float(l[30:38])
    self.y = float(l[38:46])
    self.z = float(l[46:54])
    self.occupancy = float(l[54:60])
    self.tempFactor = float(l[60:66])
    self.chain = l[72:76].strip()
    self.element = l[76:].strip()
    if self.element == "":
      self.element = " "
    lastID += 1
    return True
  def Distance(self, atom):
    vectorTo = [atom.x - self.x,
    atom.y - self.y, atom.z - self.z]
    return math.sqrt(vectorTo[0]**2 + vectorTo[1]**2 + vectorTo[2]**2)
  def PDBLine(self):
    # Specific output from VMD 1.9.1 for easy diff use
    if len(self.name) < 4:
      return "ATOM {0:>6}  {1:<4}{2:<4}{3}{4:>4}    \
{5:>8.3f}{6:>8.3f}{7:>8.3f}{8:6.2f}{9:6.2f}      {10:<4}{11:>2}\n".format(
self.id,self.name,self.resname,self.chainID,self.resID,self.x,self.y,self.z,self.occupancy,self.tempFactor,
self.chain,self.element)
    else:
      return "ATOM {0:>6} {1:<5}{2:<4}{3}{4:>4}    \
{5:>8.3f}{6:>8.3f}{7:>8.3f}{8:6.2f}{9:6.2f}      {10:<4}{11:>2}\n".format(
self.id,self.name,self.resname,self.chainID,self.resID,self.x,self.y,self.z,self.occupancy,self.tempFactor,
self.chain,self.element)
  
def ReadPDB(filename):
  try:
    f = open(filename, 'r')
  except:
    print("Failed to open",filename,"for reading.")
    return False
  p = PDBFile()
  p.filename = filename
  global lastID
  lastID = 0
  for l in f:
    if l[:4] == "ATOM":
      a = Atom()
      if a.ReadLine(l):
        p.atoms.append(a)
        continue
      try:
        a.id = int(l[4:11])
      except:
        try:
          a.id = int(l[4:11],16)
        except:
          try:
            if '*' in l[4:11]:
              a.id = lastID
              lastID += 1
          except:
            print("Error reading ATOM entry:",l)
            print("\tSkipping...")
            continue
      a.name = l[11:17].strip()
      a.resname = l[17:21].strip()
      a.chain = l[21:22].strip()
      a.resid = int(l[22:26])
      a.x = float(l[30:38])
      a.y = float(l[38:46])
      a.z = float(l[46:54])
      p.atoms.append(a)
    elif l[:3] == "END":
      break
  f.close()
  p.nAtoms = len(p.atoms)
  return p

