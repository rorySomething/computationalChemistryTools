#!/bin/python3

from math import sqrt
import math

class CalcData:
  def __init__(self):
    self.complete = False

class vec3:
  x = y = z = 0.
  def __init__(self, x = [0.,0.,0.]):
    self.x = float(x[0]); self.y = float(x[1]); self.z = float(x[2])
  def Distance(self, v):
    x = self.x - v.x; y = self.y - v.y; z = self.z - v.z
    return sqrt(x*x+y*y+z*z)
  def __add__(self,other):
    return vec3([self.x+other.x,self.y+other.y,self.z+other.z])
  def __sub__(self,other):
    return vec3([self.x-other.x,self.y-other.y,self.z-other.z])
  def __mul__(self,other):
    return vec3([self.x*float(other),self.y*float(other),self.z*float(other)])
  __rmul__ = __mul__
  def __div__(self,other):
    return vec3([self.x/float(other),self.y/float(other),self.z/float(other)])
  def __eq__(self,other):
    return self.x == other.x and self.y == other.y and self.z == other.z
  def __str__(self):
    return "{0:< 8.3f}, {1:< 8.3f}, {2:< 8.3f}".format(self.x,self.y,self.z)
  def __neg__(self):
    return vec3([-self.x, -self.y, -self.z])
  def Length(self):
    return sqrt(self.x*self.x+self.y*self.y+self.z*self.z)
  def Normalize(self):
    l = self.Length()
    if l == 0.:
      return vec3([0.,0.,0.])
    return vec3([self.x/l, self.y/l, self.z/l])
  def Dot(self, other):
    return self.x*other.x + self.y*other.y + self.z*other.z
  def Cross(self, other):
    return vec3([
    self.y*other.z-self.z*other.y,
    self.x*other.z-self.z*other.x,
    self.x*other.y-self.y*other.x
    ])    
  def __truediv__(self,other):
    return vec3([self.x/float(other),self.y/float(other),self.z/float(other)])
  """Assuming this is in cartesian coordinates, convert to fractional with
  coordinate vectors A, B and C"""
  def Fractional(self, A,B,C):
      la = A.Length(); lb = A.Length(); lc = A.Length()
      #volume = A.Dot(B.Cross(C)) # Triple Scalar Product
      # volume = sqrt(1 - cosa**2 - cosb**2 - cosg**2 + 2*ca*cb*cg)
      # 
      #
      cosA = B.Normalize().Dot(C.Normalize())
      cosB = A.Normalize().Dot(C.Normalize())
      cosG = A.Normalize().Dot(B.Normalize())
      sinA = sqrt(1. - cosA**2)
      sinB = sqrt(1. - cosB**2)
      sinG = sqrt(1. - cosG**2)
      asg = la*sinG
      volume = sqrt(1. - cosA**2 - cosB**2 - cosG**2 + 2*cosA*cosB*cosG)
      
      r1 = vec3([ 1./la, -cosG/(la*sinG), (cosA*cosG-cosB)/(la*sinG*volume) ])
      r2 = vec3([ 0.,    sinG/lb,         (cosB*cosG-cosA)/(lb*sinG*volume) ])
      r3 = vec3([ 0.,    0.,                          sinG/(lc*volume) ])
      # Cartesian to fractional from C++ code reference
      # frac = transpose(conversion_matrix) * xyz
      # conversion_matrix = {
      #   ( 1./Length(A), -cos(gamma)/(Length(A)*sin(gamma)), (cosa*cosg-cosb)/(la*vol*sing) ),
      #   ( 0., (1./lb)*sing, (cosb * cosg - cosa)/(lb*vol*sing) ),
      #   ( 0., 0., sing/(c*vol) )
      # }
      # Volume = sqrt(1 - cosa**2 - cosb**2 - cosg**2 + 2*cosa*cosb*cosg)
      #   ??? Triple product?
      #     = | (axb) dot c |
      #     = l(axb) * lc * costheta
      #
      # Matrix * vector = v.x*column1, v.y*column2, v.z*column3 ;  sum rows
      #  = Dot(vector, column)
      # CosA = 
      #
      return vec3([ self.Dot(r1), self.Dot(r2), self.Dot(r3) ])

    
class Atom:
  def __init__(self, pos = vec3([0.,0.,0.]), element = 0, name = ''):
    self.p = pos
    self.element = element
    self.name = name
  def SortedNeighbors(self, atoms):
    d = []
    l = []
    n = 0
    for i in atoms:
      d.append(self.p.Distance(i.p))
      l.append(n)
      n+=1
    swap = True
    while swap:
      swap = False
      for i in range(len(l)-1):
        a = l[i]; b = l[i+1]
        if d[b] < d[a]:
          swap = True
          l[i] = b; l[i+1] = a
    return l
  def Distances(self, l, atoms):
    d = []
    for i in range(len(l)):
      d.append(self.p.Distance(atoms[l[i]].p))
    return d
  def Shells(self, l, atoms, tol):
    d = self.Distances(l, atoms)
    last = d[1]
    sums = 0
    n = 0
    s = []
    for i in range(1,len(d)):
      if (d[i]-last) > tol:
        s.append(sums/n)
        s.append(n); sums = d[i]; n = 1
      else:
        n+=1; sums += d[i]
      last = d[i]
    if n != 0:
      s.append(sums/n)
      s.append(n);
    return s

types = []
atoms = []
title = ''
scale = 1.0
direct = ''
a = b = c = vec3([0.,0.,0.])

class CrystalCell:
  pass

def ReadCONTCAR(fn = "CONTCAR", verbose = False):
  try:
    f = open(fn, 'r')
  except:
    print("Can't open",fn,"for reading.")
    return False
  cell = CrystalCell()
  cell.atoms = []
  cell.types = []
  cell.counts = []
  cell.scale = 0.0
  cell.direct = False
  cell.lattice = []
  cell.title = f.readline() # Title
  cell.scale = float(f.readline()) # Lattice scale
  cell.labeled = False
  l = f.readline().split() # A, lattice vector
  cell.lattice.append(vec3(l)) #vec3(float(l[0]), float(l[1]), float(l[2]))
  l = f.readline().split() # B
  cell.lattice.append(vec3(l)) #vec3(float(l[0]), float(l[1]), float(l[2]))
  l = f.readline().split() # C
  cell.lattice.append(vec3(l)) #vec3(float(l[0]), float(l[1]), float(l[2]))
  natoms = f.readline().split() # Atom count list or name list
  try:
    n = int(natoms[0])
  except:
    # Elements are labeled?
    labels = True
    cell.labeled = True
  if labels:
    l = f.readline().split()
    for i in l:
      cell.counts.append(int(i))
    for i in natoms:
      cell.types.append(i)
  else:
    for i in natoms:
      cell.counts.append(int(i))
      cell.types.append(i)
  total = 0
  for i in cell.counts:
    total += i
  if verbose:
    print("Found {0} types and {1} total atoms:".format(len(cell.counts), total), cell.counts)
  cell.direct = (f.readline().lower() == "direct")
  for i in range(len(cell.counts)):
    for j in range(cell.counts[i]):
      l = f.readline().split()
      if not cell.direct:
        l = l[:3]
      if labels:
        cell.atoms.append(Atom(vec3(l),i,natoms[i]))
      else:
        cell.atoms.append(Atom(vec3(l),i,str(i+1)))
  if verbose:
    print("Coordinates stored:",len(cell.atoms))
  #Ignore velocities
  f.close()
  return cell
  

def PrintOUTCAR(directory):
  c = CalcData()
  try:
    f = open(directory+"OUTCAR", 'r')
  except:
    print("No "+directory+"OUTCAR")
    return c
  c.Ngeom = 0; c.Nelec = 0
  while True:
    l = f.readline()
    if not l:
      break
    if "Elapsed time" in l:
      c.Elapsed = float(l.split()[3])
      c.complete = True
    elif "energy  without entropy" in l:
      l = l.split()
      c.lastEnergy = float(l[3])
      c.lastEnergySigma = float(l[6])
      c.Ngeom += 1
    elif "energy without entropy" in l:
      c.Nelec += 1
    elif "NKPTS =" in l:
      c.NKPTS = int(l.split()[3])
    elif "E-fermi" in l:
      l = l.split()
      c.eFermi = float(l[2])
      c.XC = float(l[4])
      if len(l[6]) > 1:
        c.alphaBeta = float(l[6][1:])
      else:
        c.alphaBeta = float(l[7])
    elif "EDIFF  =" in l:
      c.EDIFF = float(l.split()[2])
    elif "ENCUT  =" in l:
      c.ENCUT = float(l.split()[2])
  f.close()
  if not c.complete:
    print("Calculation incomplete.\n\
Steps Elec., Geom. {0},{1}".format(c.Nelec, c.Ngeom))
    return c
  print("ENCUT: {0:<10.1f} EDIFF: {1}".format(c.ENCUT, c.EDIFF))
  log = math.log10(c.EDIFF)
  precision = 0
  if log < 0:
    precision = -math.floor(log)
  print("Final Energy: {0: .{prec}f}\tSigma->0: {1: .{prec}f}\n\
E-fermi: {6: .4f}\tAlpha+Beta: {7: .4f}\n\
Steps Elec., Geom. {2},{3}\n\
NKPTS: {5:3}\n\
Elapsed: {4: .3f}".format(c.lastEnergy, c.lastEnergySigma, c.Nelec, c.Ngeom,
c.Elapsed, c.NKPTS, c.eFermi+c.alphaBeta, c.alphaBeta, prec=precision))
  return c

  
def ReadOUTCAR(directory):
  c = CalcData()
  try:
    f = open(directory+"OUTCAR", 'r')
  except:
    print("No "+directory+"OUTCAR")
    return c
  c.Ngeom = 0; c.Nelec = 0
  while True:
    l = f.readline()
    if not l:
      break
    if "Elapsed time" in l:
      c.Elapsed = float(l.split()[3])
      c.complete = True
    elif "energy  without entropy" in l:
      l = l.split()
      c.lastEnergy = float(l[3])
      c.lastEnergySigma = float(l[6])
      c.Ngeom += 1
    elif "energy without entropy" in l:
      c.Nelec += 1
    elif "NKPTS =" in l:
      c.NKPTS = int(l.split()[3])
    elif "E-fermi" in l:
      l = l.split()
      c.eFermi = float(l[2])
      c.XC = float(l[4])
      if len(l[6]) > 1:
        c.alphaBeta = float(l[6][1:])
      else:
        c.alphaBeta = float(l[7])
    elif "EDIFF  =" in l:
      c.EDIFF = float(l.split()[2])
    elif "ENCUT  =" in l:
      c.ENCUT = float(l.split()[2])
  f.close()
  return c


def GetMPGrid(directory = './'):
  fn = directory+"KPOINTS"
  mp = "N/A"
  with open(fn, 'r') as f:
    l = f.readline()
    if "Automatic" in l:
      # Okay
      pass
    else:
      print(fn+" bad format")
    l = f.readline() # 0
    l = f.readline()
    if "Gamma" in l:
      mp = "G"
    elif "Monkhorst" in l:
      mp = "MP"
    else:
      print(fn+" bad format?")
      print("\t"+l)
    l = f.readline()
    mp += " {0}x{1}x{2}".format(*l.split())
  f.close()
  return mp


def GetData(directory):
  complete = False; contcar = False; poscar = False
  contcar = ReadCONTCAR(directory+"CONTCAR")
  if not contcar:
    print("No CONTCAR")
  poscar = ReadCONTCAR(directory+"POSCAR")
  if not poscar:
    print("No POSCAR")
  if contcar:
    print("CONTCAR:",len(contcar.atoms),"atoms")
    a = contcar.lattice[0]
    b = contcar.lattice[1]
    c = contcar.lattice[2]
    print("Lattice Vectors:\n\
    a: {0.x: .4f}, {0.y: .4f}, {0.z: .4f}\n\
    b: {1.x: .4f}, {1.y: .4f}, {1.z: .4f}\n\
    c: {2.x: .4f}, {2.y: .4f}, {2.z: .4f}".format(a,b,c))
    print("Lattice Parameters")
    la = a.Length(); lb = b.Length(); lc = c.Length()
    alpha = b.Dot(c)/lb/lc; beta = a.Dot(c)/la/lc; gamma = a.Dot(b)/la/lb;
    print("a, b, c {0:.4f}, {1:.4f}, {2:.4f}".format(la, lb, lc))
    print("alpha, beta, gamma {0:.3f}, {1:.3f}, {2:.3f}".format(
    math.degrees(math.acos(alpha)), math.degrees(math.acos(beta)),
    math.degrees(math.acos(gamma))))
    print("Volume: {0:.3f} Ang^3".format((a.Dot(b.Cross(c)))))
  elif poscar:
    print("POSCAR:",len(poscar.atoms),"atoms")
    a = poscar.lattice[0]
    b = poscar.lattice[1]
    c = poscar.lattice[2]
    print("Lattice Vectors:\n\
    a: {0.x: .4f}, {0.y: .4f}, {0.z: .4f}\n\
    b: {1.x: .4f}, {1.y: .4f}, {1.z: .4f}\n\
    c: {2.x: .4f}, {2.y: .4f}, {2.z: .4f}".format(a,b,c))
    print("Lattice Parameters")
    la = a.Length(); lb = b.Length(); lc = c.Length()
    alpha = b.Dot(c)/lb/lc; beta = a.Dot(c)/la/lc; gamma = a.Dot(b)/la/lb;
    print("a, b, c {0:.4f}, {1:.4f}, {2:.4f}".format(la, lb, lc))
    print("alpha, beta, gamma {0:.3f}, {1:.3f}, {2:.3f}".format(
    math.degrees(math.acos(alpha)), math.degrees(math.acos(beta)),
    math.degrees(math.acos(gamma))))
    print("Volume: {0:.3f} Ang^3".format((a.Dot(b.Cross(c)))))
  print("KPOINTS: "+GetMPGrid(directory))
  if PrintOUTCAR(directory).complete:
    print("Successful calculation.")
    
def ReadData(directory):
  complete = False; contcar = False; poscar = False
  contcar = ReadCONTCAR(directory+"CONTCAR")
  a = vec3(); b = vec3(); c = vec3()
  la = 0.; lb = 0.; lc = 0.;
  alpha = 0.; beta = 0.; gamma = 0.
  if not contcar:
    print("No CONTCAR")
  poscar = ReadCONTCAR(directory+"POSCAR")
  if not poscar:
    print("No POSCAR")
  if contcar:
    atoms = len(contcar.atoms)
    a = contcar.lattice[0]
    b = contcar.lattice[1]
    c = contcar.lattice[2]
    la = a.Length(); lb = b.Length(); lc = c.Length()
    alpha = b.Dot(c)/lb/lc; beta = a.Dot(c)/la/lc; gamma = a.Dot(b)/la/lb;
  elif poscar:
    atoms = len(poscar.atoms)
    a = poscar.lattice[0]
    b = poscar.lattice[1]
    c = poscar.lattice[2]
    la = a.Length(); lb = b.Length(); lc = c.Length()
    alpha = b.Dot(c)/lb/lc; beta = a.Dot(c)/la/lc; gamma = a.Dot(b)/la/lb;
  calc = ReadOUTCAR(directory)
  calc.directory = directory
  alpha = math.degrees(math.acos(alpha))
  beta = math.degrees(math.acos(beta))
  gamma = math.degrees(math.acos(gamma))
  calc.A = a; calc.a = la; calc.B = b; calc.b = lb; calc.C = c; calc.c = lc
  calc.alpha = alpha; calc.beta = beta; calc.gamma = gamma
  calc.volume = a.Dot(b.Cross(c))
  calc.KPOINTS = GetMPGrid(directory)
  calc.Natoms = atoms
  return calc
  
def PrintList(filename = "Vasp_Runs.txt", data = None):
  if not data:
    return
  with open(filename, 'w') as f:
    f.write("\
Directory  Energy       NKPTS  E-fermi  Alpha+  Elapsed Time  a  b  c  alpha  beta  gamma  NSCF  Ngeom\n\
           without                  Beta                                                     steps\n\
           entropy  KPOINTS\n")
    for d in data:
      if d.complete:
        #         Dir    Energy    NKPT    E-f       A+B       Elapsed
        f.write("{0:9}  {1:>11.5f}  {2:5}  {3: .4f}  {4: .4f}  {5:.3f}\
    {6: .4f}  {7: .4f}  {8: .4f}  {9: .3f}  {10: .3f}  {11: .3f}\
    {12}  {13} {14}\n".format(d.directory, d.lastEnergy, d.NKPTS, d.eFermi+d.alphaBeta,
    d.alphaBeta, d.Elapsed, d.a, d.b, d.c, d.alpha, d.beta, d.gamma, d.Nelec, d.Ngeom, d.KPOINTS))
      else:
        f.write("{0:9}  Incomplete  NSCF: {1}  Ngeom: {2} KPOINTS: {3}\n".format(d.directory, d.Nelec, d.Ngeom, d.KPOINTS))
  f.close()

def SystemInfo(args):
  i = 1
  directories = []
  bList = False
  while i < len(args):
    a = args[i]
    if len(a) > 1 and a[0] == '-':
      if a[1] == 'h' or a[1:3] == '-h':
        print('Reads a VASP POSCAR or CONTCAR and OUTCAR and outputs general info')
        print('Default behavior is to read from current directory')
        print('Otherwise reads from directories provided')
        return
      elif len(a) == 2 and a[1] == 'l':
        bList = True
    else:
      directories.append(a)
    i += 1
  if len(directories) == 0:
    directories.append(".")
  # Start doing things
  if bList:
    data = []
    for d in directories:
      if d[-1] != "/":
        d += "/"
      print(d)
      data.append(ReadData(d))
    PrintList("Vasp_Runs.txt", data)
  else:
    for d in directories:
      if d[-1] != "/":
        d += "/"
      print(d)
      GetData(d)
  print("Good Day.")

if __name__ == "__main__":
  import sys
  SystemInfo(sys.argv)

