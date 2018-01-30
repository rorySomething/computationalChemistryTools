"""
Handle PSF files related to NAMD use
"""
class PSFFile:
  def __init__(self):
    self.atoms = []
    # index chain resid resname atomname atomtype charge mass 0
    self.bonds = []
    self.angles = []
    self.dihedrals = []
    self.impropers = []
    self.filename = ""
    self.titleLines = []
    self.nAtoms = 0
  def CheckCharges(self, topFile):
    # First check sum
    charge = 0.
    fullTest = True
    missing = []
    for i in self.atoms:
      charge += float(i[6])
      resname = i[3].strip()
      if resname in topFile.residues: # If residue is known
        res = topFile.residues[resname]
        if i[4] in res.atoms: # If atom type is in residue
          atom = res.atoms[i[4]]
          if abs(float(i[6]) - float(atom.charge)) > 1e-4: # Test charge match with type
            print("Error: Charges don't match for atom",i[4]," of type",i[5],"in residue", i[3])
            print("\t",float(i[6]),"vs.",atom.charge)
        else:
          print("Error: Atom type",i[4],"is not defined in its assigned residue", i[3])
      else:
        if resname not in missing:
          missing.append(resname)
          print("Can't find topology data for residue",resname)
          print(sorted(topFile.residues.keys()))
        fullTest = False
    if abs(charge) > 1e-5:
      print("WARNING: The system is charged",charge)
    else:
      print("Total system charge is neglegable",charge)
    if fullTest:
      print("PSF charges match Topology charges.")
    else:
      print("Could not ensure that PSF charges match Topology charges.")
      
  # Make sure no atoms are missing
  def CheckMolecules(self, topFile):
    pass
class PSFAtom:
  pass
class TOPFile:
  def __init__(self):
    self.residues = dict()
    self.filenames = []
    self.nResidues = 0
  def CheckCharges(self):
    for i in self.residues.values():
      i.CheckCharge(True)
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
      
def ReadPSF(filename, onlyAtoms = False):
  try:
    f = open(filename, 'r')
  except:
    print("Failed to open",filename,"for reading.")
    return False
  p = PSFFile()
  p.filename = filename
  for l in f:
    if "!NTITLE" in l:
      p.titleLines = []
      for i in range(int(l.split()[0])):
        p.titleLines.append(f.readline())
    elif "!NATOM" in l:
      p.nAtoms = int(l.split()[0])
      for i in range(p.nAtoms):
        p.atoms.append(f.readline().split())
      if onlyAtoms:
        break
    elif "!NBOND" in l:
      p.nBonds = int(l.split()[0])
      bonds = 0
      while bonds < p.nBonds:
        a = f.readline().split()
        j = 0
        while j < (len(a)-1):
          # -1 to directly index from 0
          p.bonds.append( ( int(a[j]) - 1, int(a[j+1]) - 1 ) )
          j += 2
          bonds += 1
    elif "!NTHETA" in l:
      p.nAngles = int(l.split()[0])
      angles = 0
      while angles < p.nAngles:
        a = f.readline().split()
        j = 0
        while j < (len(a)-2):
          # -1 to directly index from 0
          p.angles.append( ( int(a[j]) - 1, int(a[j+1]) - 1, int(a[j+2]) - 1 ) )
          j += 3
          angles += 1
    elif "!NPHI" in l:
      p.nDihedrals = int(l.split()[0])
      dihed = 0
      while dihed < p.nDihedrals:
        a = f.readline().split()
        j = 0
        while j < (len(a)-3):
          # -1 to directly index from 0
          p.dihedrals.append( ( int(a[j]) - 1, int(a[j+1]) - 1,
          int(a[j+2]) - 1, int(a[j+3]) - 1 ) )
          j += 4
          dihed += 1
    elif "!NIMPHI" in l:
      p.nImpropers = int(l.split()[0])
      imp = 0
      while imp < p.nImpropers:
        a = f.readline().split()
        j = 0
        while j < (len(a)-3):
          # -1 to directly index from 0
          p.impropers.append( ( int(a[j]) - 1, int(a[j+1]) - 1,
          int(a[j+2]) - 1, int(a[j+3]) - 1 ) )
          j += 4
          imp += 1
    elif "!NDON" in l:
      p.nDonors = int(l.split()[0])
      if p.nDonors > 0:
        print("Ignoring Donors")
    elif "!NACC" in l:
      p.nAcceptors = int(l.split()[0])
      if p.nAcceptors > 0:
        print("Ignoring Acceptors")
    elif "!NNB" in l:
      p.nNonbonding = int(l.split()[0])
      if p.nNonbonding > 0:
        print("Ignoring Nonbonded")
    elif "!NGRP" in l:
      p.nGRP = int(l.split()[0])
      # Don't know what this is
  f.close()
  return p
