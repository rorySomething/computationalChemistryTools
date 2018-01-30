# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 22:36:54 2017

@author: Rory
A full simulation cell consisting of multiple molecules/residues
"""
# TODO Incomplete implementation

from mdResidue import Residue
from mdAtom import Atom
from vec3 import Vec3
from math import sin,cos,asin,acos,radians,degrees

class System:
  def __init__(self):
    self.cell = [0.,0.,0.,90.,90.,90.]
    self.atoms = [] # Atom data in list
    self.residues = [] # Residues in list
    self.chains = []
  def addXYZ(self, filename):
    try:
      f = open(filename, 'r')
    except:
      print("ERROR: Failed to open",filename)
      return False
    atoms = int(f.readline())
    title = f.readline()
    for a in range(atoms):
      self.atoms.append(Atom())
      self.atoms[-1].readXYZLine(f.readline(), len(self.atoms)-1)
    f.close()
  def addPDB(self, filename):
    pass
  def writeXYZ(self, filename):
    pass
  def writePDB(self, filename):
    pass
  def writePSF(self, filename):
    pass
  def writeLAMMPSData(self, filename):
    pass
  def mergeTopologyData(self, filename):
    pass
  def addShells(self, atomNameList):
    pass
  # A=X,B,C or C=Z,A,B
  def getCellVectors(self, AequalsX = True):
    a = self.cell[0]
    b = self.cell[1]
    c = self.cell[2]
    alpha = radians(self.cell[3])
    beta = radians(self.cell[4])
    gamma = radians(self.cell[5])
    sa = sin(alpha); sb = sin(beta); sg = sin(gamma)
    ca = cos(alpha); cb = cos(beta); cg = cos(gamma)
    if AequalsX:
      #TODO Fix this
      A = Vec3([a, 0., 0.])
      B = Vec3([
          sg*b,
          cg*b,
          0. ])
      C = Vec3([
          sb*c,
          sa*cb*c,
          ca*cb*c ])
      return (A,B,C)
    else:
      C = Vec3([0.,0.,c])
      A = Vec3([cb*a, 0., sb*a])
      B = Vec3([cg*ca*b,
                sg*b])
      return (A,B,C)
  def getVolume(self):
    cell = self.getCellVectors()
    # Triple scalar product
    return cell[0].cross(cell[1]).dot(cell[2])