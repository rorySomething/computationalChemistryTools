#!/bin/python3

from math import sqrt
import math

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

def StructureInfo(args):
  i = 1
  read = []
  crystal = ""
  defect = False
  prec = False
  p = 0
  scale = []
  outscale = []
  while i < len(args):
    a = args[i]
    if len(a) > 1 and a[0] == '-':
      if a[1] == 'h' or a[1:3] == '-h':
        print('Reads a VASP POSCAR or CONTCAR and outputs general info')
        print('Default behavior is to read from CONTCAR')
        print('Can give a list of POS/CONTCARs to read.')
        print('-p (exponent) optional to control output precision')
        print('\te.g. -p 4 prints 100.0001')
        print('\t     -p 6 prints 100.000127')
        print('\tDefault is 10^-4 for coordinates and 10^-3 for angles')
        return
      elif '-p' in a and len(args) > i+1:
        prec = True
        p = int(args[i+1])
        i += 1
    else:
        read.append(a)
    i += 1
  if len(read) == 0:
    read.append("CONTCAR")
  # Start doing things
  for r in read:
    print("Read Cell: "+r)
    cell = ReadCONTCAR(r, True)
    print(len(cell.atoms),"atoms")
    a = cell.lattice[0]
    b = cell.lattice[1]
    c = cell.lattice[2]
    print("Lattice Vectors:\n\
    a: {0.x}, {0.y}, {0.z}\n\
    b: {1.x}, {1.y}, {1.z}\n\
    c: {2.x}, {2.y}, {2.z}".format(a,b,c))
    print("Lattice Parameters")
    la = a.Length(); lb = b.Length(); lc = c.Length()
    alpha = b.Dot(c)/lb/lc; beta = a.Dot(c)/la/lc; gamma = a.Dot(b)/la/lb;
    if prec:
      print(("a, b, c {0:."+str(p)+"f}, {1:."+
      str(p)+"f}, {2:."+str(p)+"f}").format(la, lb, lc))
    else:
      print("a, b, c {0:.4f}, {1:.4f}, {2:.4f}".format(la, lb, lc))
    if prec:
      print(("alpha, beta, gamma {0:."+str(p)+"f}, {1:."+str(p)+"f}, {2:."+str(p)+"f}").format(
      math.degrees(math.acos(alpha)), math.degrees(math.acos(beta)),
      math.degrees(math.acos(gamma))))
    else:
      print("alpha, beta, gamma {0:.3f}, {1:.3f}, {2:.3f}".format(
      math.degrees(math.acos(alpha)), math.degrees(math.acos(beta)),
      math.degrees(math.acos(gamma))))
    if prec:
      print(("Volume: {0:."+str(p)+"f} Ang^3").format((a.Dot(b.Cross(c)))))
    else:
      print("Volume: {0:.3f} Ang^3".format((a.Dot(b.Cross(c)))))
  print("Good Day!")

if __name__ == "__main__":
  import sys
  StructureInfo(sys.argv)

