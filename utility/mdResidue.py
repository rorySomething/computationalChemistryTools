# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 15:14:18 2017

@author: Rory
"""

# TODO Incomplete implementation

from mdAtom import Atom

class Residue:
  def __init__(self):
    self.name = "Unk"
    self.index = 0
    self.type = "Unk"
    self.chainIndex = 0
    self.chainName = "Unk"
    self.atoms = [] # Atom index list
    self.bonds = [] # pairs of atom indices
    self.angles = [] # triplets of atom indices
    self.dihedrals = [] # 4 atom indices
    self.impropers = [] # center,A,B,C