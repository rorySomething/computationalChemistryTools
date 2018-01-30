"""
Classes for manipulating VASP I/O
Python 3
"""

import vec3

class KPOINTS:
  x = y = z = w = 0
  bands = []
  def Print(self):
    print(self.x,self.y,self.z,self.w)
    for i in range(len(self.bands)):
      print(i+1,"\t",self.bands[i])
  def __init__(self):
    self.x = self.y = self.z = 0
    self.bands = []
  def Coord(self):
    a = [self.x,self.y,self.z,self.w]
    return a

class EIGENVAL:
  l1 = []
  l2 = []
  a3 = 0
  CAR = ""
  title = ""
  valence = kpts = bands = 0
  kpoints = []
  polarized = False
  def Print(self):
    print(self.l1,"\n",self.l2,"\n",self.a3,"\n",
      self.CAR,"\n",self.title,"\n",self.valence,self.kpts,self.bands,
      "\nPolarized:",self.polarized)
    for i in range(len(self.kpoints)):
      print("kpoint:",'{:>4}'.format(i+1),'  ({0:>7.6f}, {1:>7.6f}, {2:>7.6f}, {3:>7.6f})'.format(*self.kpoints[i].Coord()))
    #  self.kpoints[i].Print()
    print("Total Bands:",self.bands)
    print("Kpoints:",len(self.kpoints))
    #for i in self.kpoints:
    #  print("Bands:",'{:>4}'.format(len(i.bands))," ",end='\t')

class OUTCAR:
  def __init__(self):
    self.valence = 0
    self.kpts = 0
    self.bands = 0
    self.polarized = False
    self.kpoints = []
    self.eFermi = 0
    self.XC = 0
    self.alphaBeta = 0

class CONTCAR:
  """
    Apply scale to lattice automatically
  """
  def __init__(self, filename = None, applyScale = True):
    self.title = ""
    self.scale = 1.0
    self.lattice = []
    self.atomTypes = []
    self.atomCounts = []
    self.lattice = []
    self.direct = True
    self.dynamics = False
    self.atoms = []
    self.vectors = []
    # TODO Add error checking
    if filename:
      print("Debug: Opening POSCAR-type:",filename)
      with open(filename, 'r') as f:
        self.title = f.readline()
        self.scale = float(f.readline())
        self.lattice = [
            [float(x)*self.scale for x in f.readline().split()],
            [float(x)*self.scale for x in f.readline().split()],
            [float(x)*self.scale for x in f.readline().split()] ]
        l = f.readline().split()
        labeled = False
        try:
          i = int(l[0])
        except:
          labeled = True
        if labeled:
          self.atomTypes = l[:]
          self.atomCounts = [int(x) for x in f.readline().split()]
        else:
          self.atomTypes = [str(x) for x in range(len(l))]
          self.atomCounts = [int(x) for x in f.readline().split()]
        l = f.readline().lower()
        if l[:6] == "select":
          self.dynamics = True
        elif l == "direct":
          self.direct = True
        if self.dynamics:
          self.direct = (f.readline().lower() == "direct")
        # Now read coordinates
        numbers = []
        i = 0
        for c in self.atomCounts:
          for x in range(c):
            numbers.append(i)
          i += 1 # index of atom type next set of atoms is +1
        # Actual coordiate reading now
        for i in numbers:
          l = f.readline().split()
          self.atoms.append([]) # New atom, just a list
          self.atoms[-1].append([float(x) for x in l[:3]]) # xyz
          self.atoms[-1].append([self.atomTypes[i]]) # Name
          if self.dynamics:
            self.atoms[-1].append([l[3:]]) # list of frozen coordinates
          else:
            self.atoms[-1].append(['T', 'T', 'T']) # unfrozen
        # Blank line
        f.readline()
        # Vectors?
        l = f.readline()
        if not l:
          # EOF, no vectors
          self.vectors = [[0.,0.,0.] for x in range(len(self.atoms))]
        else:
          # Read in vectors
          for i in range(len(self.atoms)):
            self.vectors.append([float(x) for x in l.split()])
            l = f.readline() # Kinda out of normal order because of earlier check

