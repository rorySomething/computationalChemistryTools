#!/bin/python3

# Rory Vander Valk
# Python3 - rotate.py

from math import sqrt

# Rory Vander Valk
# Python3 - vec3.py

import math

class vec3:
  def __init__(self, x = [0.,0.,0.]):
    if not hasattr(x, '__getitem__') or len(x) < 3:
      # Not a vector, ignore it
      self.x,self.y,self.z = 0.,0.,0.
    else:
      self.x = float(x[0]); self.y = float(x[1]); self.z = float(x[2])
  def lengthSqr(self):
    # return self.dot(self)
    return self.x*self.x + self.y*self.y + self.z*self.z
  def length(self):
    # return math.sqrt(self.lengthSqr())
    return math.sqrt(self.x*self.x+self.y*self.y+self.z*self.z)
  def vectorTo(self, target):
    return vec3([
    target.x - self.x,
    target.y - self.y,
    target.z - self.z])
  def distanceSqr(self, other):
    # return self.vectorTo(other).lengthSqr()
    x = self.x - other.x; y = self.y - other.y; z = self.z - other.z
    return x*x+y*y+z*z
  def distance(self, other):
    return math.sqrt(self.distanceSqr(other))
  def __add__(self,other):
    return vec3([self.x+other.x,self.y+other.y,self.z+other.z])
  __radd__ = __add__
  def __sub__(self,other):
    return vec3([self.x-other.x,self.y-other.y,self.z-other.z])
  # TODO This may be wrong...
  # A -= B should be A = A-B so __rsub__ = __sub__ ???
  def __rsub__(self, x):
    return vec3([-self.x+x.x,-self.y+x.y,-self.z+x.z])
  # TODO Add a hasattr check here?
  def __mul__(self,other):
    return vec3([self.x*float(other),self.y*float(other),self.z*float(other)])
  __rmul__ = __mul__
  def __div__(self,other):
    return vec3([self.x/float(other),self.y/float(other),self.z/float(other)])
  def __eq__(self,other):
    return self.x == other.x and self.y == other.y and self.z == other.z
  def __str__(self):
    return "{0:< 8.3e}, {1:< 8.3e}, {2:< 8.3e}".format(self.x,self.y,self.z)
  def __neg__(self):
    return vec3([-self.x, -self.y, -self.z])
  def getList(self):
    return [self.x, self.y, self.z]
  def vectorMultiply(self, b):
    return vec3([self.x*b.x, self.y*b.y, self.z*b.z])
  def normalize(self):
    l = self.length()
    # TODO Better if we let a program fail silently
    #      by returning [1e-30, 0., 0.] when length=0.
    if l == 0.:
      return vec3([0.,0.,0.])
    return vec3([self.x/l, self.y/l, self.z/l])
  def dot(self, b):
    return (self.x*b.x + self.y*b.y + self.z*b.z)
  # self x B -> AxB -> right-handed
  def cross(self, b):
    return vec3([
    self.y*b.z-self.z*b.y,
    -(self.x*b.z-self.z*b.x),
    self.x*b.y-self.y*b.x])
  #
  # [ i  j  k  ]
  # [ ax ay az ]
  # [ bx by bz ]
  #
  # { i*(ay*bz-az*by) - j*(ax*bz-az*bx) + k*(ax*by-bx*ay) }
  # i = (ay*bz-az*by); j = -(ax*bz-az*bx); k = (ax*by-bx*ay)
  def angleDegrees(self, b):
    return math.degrees( math.acos( self.dot(b)/(self.length()*b.length()) ) )
  def __truediv__(self,other):
    return vec3([self.x/float(other),self.y/float(other),self.z/float(other)])
  """Assuming this is in cartesian coordinates, convert to fractional with
  coordinate vectors A, B and C"""
  def getFractional(self, A,B,C):
      la = A.length(); lb = B.length(); lc = C.length()
      #volume = A.Dot(B.Cross(C)) # Triple Scalar Product
      # volume = math.sqrt(1 - cosa**2 - cosb**2 - cosg**2 + 2*ca*cb*cg)
      # 
      #
      cosA = B.normalize().dot(C.normalize())
      cosB = A.normalize().dot(C.normalize())
      cosG = A.normalize().dot(B.normalize())
      sinA = math.sqrt(1. - cosA**2)
      sinB = math.sqrt(1. - cosB**2)
      sinG = math.sqrt(1. - cosG**2)
      asg = la*sinG
      volume = math.sqrt(1. - cosA**2 - cosB**2 - cosG**2 + 2*cosA*cosB*cosG)
      
      r1 = vec3([ 1./la, -cosG/(la*sinG), (cosA*cosG-cosB)/(la*sinG*volume) ])
      r2 = vec3([ 0.,    1./(lb* sinG),         (cosB*cosG-cosA)/(lb*sinG*volume) ])
      r3 = vec3([ 0.,    0.,                          sinG/(lc*volume) ])
      # Cartesian to fractional from C++ code reference
      # frac = transpose(conversion_matrix) * xyz
      # conversion_matrix = {
      #   ( 1./Length(A), -cos(gamma)/(Length(A)*sin(gamma)), (cosa*cosg-cosb)/(la*vol*sing) ),
      #   ( 0., 1./(lb*sing), (cosb * cosg - cosa)/(lb*vol*sing) ),
      #   ( 0., 0., sing/(c*vol) )
      # }
      # Volume = math.sqrt(1 - cosa**2 - cosb**2 - cosg**2 + 2*cosa*cosb*cosg)
      #   ??? Triple product?
      #     = | (axb) dot c |
      #     = l(axb) * lc * costheta
      #
      # Matrix * vector = v.x*column1, v.y*column2, v.z*column3 ;  sum rows
      #  = Dot(vector, column)
      # CosA = 
      #
      return vec3([ self.dot(r1), self.dot(r2), self.dot(r3) ])
  """Assuming this is in fractional coordinates, convert to cartesian with
     coordinate vectors A, B and C"""
  def getCartesian(self, A,B,C):
    return self.x*A + self.y*B + self.z*C
  xyz = property(getList, None)

"""
Return the position rotated about a given axis
position - vec3 initial position
vector - vec3 rotation axis, will be normalized
angleDegree - rotation angle
origin - vec3 rotate about this origin (offset position before rotation)
         Need to offset to rotate about an object's center of mass vs. center of entire system or bottom corner (0,0,0)
"""
def rotatePosition(position, vector, angleDegrees, origin = vec3([0.,0.,0.])):
  pos = position - origin
  # Rotate about center 2 atoms
  v = vector.normalize()
  rotationMatrix = []
  radians = math.radians(angleDegrees)
  cosA = math.cos(radians)
  sinA = math.sin(radians)
  ncA = 1. - cosA
  # wikipedia
  #  cosA = cos(theta) ; sinA = sin(theta)
  #  ncA = (1. - cosA)
  #  x2 = v.x*v.x ; y2 = v.y*v.y ; z2 = v.z*v.z
  #  xy = v.x*v.y ; etc...
  #  [ cosA + x2 ncA    ,  xy ncA - z sinA  ,  xz ncA + y sinA ]
  #  [ yx ncA + z sinA  ,  cosA + y2 ncA    ,  yz ncA - x sinA ]
  #  [ zx ncA - y sinA  ,  zy ncA + x sinA  ,  cosA + z2 ncA   ] = cosA + z2 cosA - z2 vs. z2 + cosA - z2 cosA
  #  xy = v.x * v.y * ncA
  #  xz = v.x * v.x * ncA
  #  yz = v.y * v.z * ncA
  x2 = v.x**2;  y2 = v.y**2;  z2 = v.z**2
  xy = v.x*v.y * ncA;  zS = v.z * sinA
  xz = v.x*v.z * ncA;  yS = v.y * sinA
  yz = v.y*v.z * ncA;  xS = v.x * sinA
  rotationMatrix.append( vec3([ cosA + x2*ncA,  xy - zS,        xz + yS]))
  rotationMatrix.append( vec3([ xy + zS,        cosA + y2*ncA,  yz - xS]))
  rotationMatrix.append( vec3([ xz - yS,        yz + xS,        cosA + z2*ncA]))
  r = vec3()
  r.x = rotationMatrix[0].dot(pos)
  r.y = rotationMatrix[1].dot(pos)
  r.z = rotationMatrix[2].dot(pos)
  r = r + origin
  return r


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

"""
Read in a VASP ~5 POSCAR/CONTCAR
Returns a CrystalCell structure
CrystalCell.atoms   Atom(.p .element .name)
CrystalCell.lattice [vec3, vec3, vec3]
CrystalCell.types
CrystalCell.counts
CrystalCell.scale
CrystalCell.direct
CrystalCell.title
CrystalCell.labeled
TODO This is an older file version
     This will fail if the file includes line like SELECTIVE DYNAMICS
"""
def ReadCONTCAR(fn = "CONTCAR", verbose = False):
  try:
    f = open(fn, 'r')
  except:
    print("Can't open",fn,"for reading.")
    return
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
      
def WritePOSCAR(cell, fn = "SUPERCELL"):
  try:
    f = open(fn, 'w')
  except:
    print("Can't open",fn,"for writing.")
    return
  f.write(cell.title[:-1]+'\n')
  f.write("{0:>21.16f}\n".format(float(cell.scale)))
  f.write(" {0.x:>22.16f}{0.y:>22.16f}{0.z:>22.16f}\n".format(cell.lattice[0]))
  f.write(" {0.x:>22.16f}{0.y:>22.16f}{0.z:>22.16f}\n".format(cell.lattice[1]))
  f.write(" {0.x:>22.16f}{0.y:>22.16f}{0.z:>22.16f}\n".format(cell.lattice[2]))
  if cell.labeled:
    for i in range(len(cell.types)):
      f.write("{0:>5}".format(cell.types[i]))
    f.write('\n')
  for i in range(len(cell.counts)):
    f.write("{0:>5}".format(cell.counts[i]))
  direct = "Direct"
  if not cell.direct:
    direct = "Direct" # No selective dynamics yet
  f.write("\n"+direct+"\n")
  if len(cell.counts) > 1:
    # Sort atoms by type
    atomsList = cell.atoms[:]
    cell.atoms.clear()
    if cell.labeled:
      for i in cell.types:
        j = 0
        while j < len(atomsList):
          if atomsList[j].name == i:
            cell.atoms.append(atomsList[j])
            del atomsList[j]
          else:
            j += 1
        
  for l in cell.atoms:
    f.write("{0.x:>20.16f}{0.y:>20.16f}{0.z:>20.16f}\n".format(l.p))
  f.write('\n')
  for l in cell.atoms:
    f.write("{0:>16.8E}{1:>16.8E}{2:>16.8E}\n".format(0.0,0.0,0.0))
  f.close()

def Rotate(args):
  i = 1
  read = "CONTCAR"
  angle = 38.213
  vector = vec3([0., 0., 1.]) # z-axis
  verbose = False
  while i < len(args):
    a = args[i]
    print(i,a)
    if len(a) == 2 and a[0] == '-':
      if a[1] == 'a':
        angle = float(args[i+1])
        i += 1
      elif a[1] == 'd':
        verbose = True
      elif a[1] == 'v':
        if len(args[i+1].split(',')) == 3:
          vector = vec3([float(x) for x in args[i+1].split(',')])
          i += 1
        else:
          vector = vec3([float(args[i+1]), float(args[i+2]), float(args[i+3])])
          i += 3
      elif a[1] == 'h' or a[1:3] == '-h':
        print('Rotate a crystal structure')
        print('Rotates about the cell center')
        print('Default behavior is to read from CONTCAR to create cell')
        print('Options: -a, -v, -d')
        print('\t-a - angle in degrees  e.g. -a 41.457')
        print('\t-v - vector to rotate about  e.g. -v 0.0,0.0,1.0')
        print('\t\tComma separated, no spaces...')
        print('\t-d - print extra debug text')
        return
    else:
      read = a
    i += 1
  # Start doing things
  if verbose:
    print("Read Cell")
  cell = ReadCONTCAR(read, verbose)
  if verbose:
    print(len(cell.atoms),"atoms")
  a = cell.lattice[0]*cell.scale
  b = cell.lattice[1]*cell.scale
  c = cell.lattice[2]*cell.scale
  if verbose:
    print("Lattice Vectors:\n\
    a: {0.x}, {0.y}, {0.z}\n\
    b: {1.x}, {1.y}, {1.z}\n\
    c: {2.x}, {2.y}, {2.z}".format(a,b,c))
  rotA = rotatePosition(a, vector, angle)
  rotB = rotatePosition(b, vector, angle)
  rotC = rotatePosition(c, vector, angle)
  center = a/2 + b/2 + c/2
  if verbose:
    print("Cell center: {0.x}, {0.y}, {0.z}".format(center))
  # Rotate atoms
  # Need to rotate atoms? or just the vectors?
  for atom in cell.atoms:
    # Atoms [.p .element .name]
    print('f',atom.p)
    frac = atom.p
    cart = frac.x*a + frac.y*b + frac.z*c
    print('f->c',cart)
    print('c->f',cart.getFractional(a,b,c))
    rotated = rotatePosition(cart, vector, angle)
    print("R",rotated)
    rotated = rotated.getFractional(rotA, rotB, rotC)
    print("Rf",rotated)
    atom.p = rotated # Update with new rotated coordinates
  cell.lattice = [rotA, rotB, rotC] # Update to rotated cell lattices
  output = 'rotated-'+read.strip('\/.')
  if verbose:
    print("Writing to", output)
  WritePOSCAR(cell, output)
  print("We're done.\nGood Day")

if __name__ == "__main__":
  import sys
  Rotate(sys.argv)

