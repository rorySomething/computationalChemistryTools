
# TODO Clean up

from residue import *

class TOPAtom:
  def __init__(self):
    self.name = ""
    self.type = ""
    self.charge = 0.

class TOPFile:
  def __init__(self):
    self.residues = dict()
    self.patches = dict()
    self.filenames = []
    self.nResidues = 0
    self.atomTypes = dict() # MASS, #, Name, mass, Element ! Comment
    self.declares = []
    self.defaultPatches = []
  def CheckCharges(self):
    for i in self.residues.values():
      i.CheckCharge(True)
  def ReadTop(self, filename):
    self.filenames.append(filename)
    with open(filename, 'r') as f:
      reading = False; patch = False
      r = Residue()
      lineNumber = 0
      for l in f:
        lineNumber += 1
        if "MASS" in l:
          l = l.split()
          if l[2] in self.atomTypes.keys():
            print("Warning overwriting atom type", l[2])
          self.atomTypes[l[2]] = (
            float(l[3]), l[4], int(l[1])) # mass, element, index
        elif "DECL" in l[:6]:
          self.declares.append(l.split()[1])
        elif "DEFA" in l[:6]:
          l = l.split()
          self.defaultPatches.append(l[2]) # FIRST
          self.defaultPatches.append(l[4]) # LAST
        elif "RESI" in l or "PRES" in l:
          if "PRES" in l:
            patch = True
          else:
            patch = False
          l = l.split()
          comment = False
          for i in l:
            if '!' in i:
              comment = True
              break
            elif 'RESI' in i or 'PRES' in i:
              break
            elif 'DELE' in i:
              comment = True # Ignore
              break
          if comment:
            continue
          if reading:
            if patch:
              if r.name in self.patches:
                print("Multiple definitions of",r.name)
              self.patches[r.name] = r
              r = Residue()
            else:
              if r.name in self.residues:
                print("Multiple definitions of",r.name)
              self.residues[r.name] = r
              r = Residue()
          else:
            reading = True
          r.name = l[1]
          try:
            r.charge = float(l[2])
          except:
            r.charge = 0.0
            print("Error reading residue charge, set to 0.0 :>",l[2])
          print("Residue:",r.name, "q =", r.charge)
        elif reading:
          if "ATOM" in l.strip()[:5]:
            l = l.split()
            comment = False
            for i in l:
              if '!' in i:
                comment = True
                break
              elif 'ATOM' in i:
                break
              elif 'DELE' in i:
                comment = True # Ignore
                break
            if comment:
              continue
            a = TOPAtom()
            a.name = l[1]
            a.type = l[2]
            try:
              a.charge = float(l[3])
            except:
              print(l)
              print("Poorly formatted topology file around line",lineNumber)
            if a.name in r.atoms:
              print("Multiple definitions of",a.name,"in residue",r.name)
            r.atoms[a.name] = a
          elif "DELE" in l.strip()[:5]:
            l = l.split()
            if l[1] == "ATOM":
              r.delete.append(l[2])
            else:
              r.delete.append(l[1])
          elif "BOND" in l or "DOUBLE" in l:
            l = l.split()[1:]
            for i in range(0, len(l)-1, 2):
              r.bonds.append((l[i],l[i+1]))
          elif "IMPR" in l:
            l = l.split()[1:]
            for i in range(0, len(l)-4, 4):
              r.impropers.append(l[i:i+3])
          elif "DONOR" in l:
            l = l.split()
            r.donors.append(l[1:])
          elif "ACCEPTOR" in l:
            l = l.split()
            r.acceptors.append(l[1:])
          elif "IC" in l[:4]:
            r.internalCoordinates.append(l.split()[1:])
      # Catch last residue
      if r.name in self.residues:
        print("Multiple definitions of",r.name)
      self.residues[r.name] = r
    self.nResidues = len(self.residues)
    f.close()
