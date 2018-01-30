# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 15:04:12 2017

@author: Rory
Atom class for use in setting up MD simulations for LAMMPS and NAMD
"""

from vec3 import vec3
from ..data.elementList import elementList

class Atom:
  def __init__(self):
    self.name = "Unk"
    self.index = 0 # Index for atom and bonding lists
    self.atomicNumber = 0
    self.elementName = "Unk"
    self.formalCharge = 0
    self.partialCharge = 0
    self.atomTypeName = "Unk"
    self.residueIndex = 0
    self.residueName = "Unk"
    self.residueType = "Unk" # Not necessary?
    self.chainIndex = 0
    self.chainName = "Unk"
    self.position = Vec3()
    self.velocity = Vec3()
    self.occupancy = 1.          # PDB/NAMD
    self.temperatureFactor = 1.  # PDB/NAMD
  def distance(self, B):
    return self.position.distance(B.position)
  def readXYZLine(self, string, index = 0):
    l = string.split()
    self.name = l[0]
    self.position.x = float(l[1])
    self.position.x = float(l[2])
    self.position.x = float(l[3])
    self.index = index
  def readLineMS55(self, string):
    pass
  def readLineVMD(self, string):
    pass
  def writeLineMS55(self):
    pass
  def writeLineVMD(self):
    pass