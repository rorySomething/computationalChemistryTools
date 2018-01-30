#!/bin/python3

from math import sqrt

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
  l = f.readline().lower()
  if l[0] == "s":
    if verbose:
      print("Selective Dynamics")
    l = f.readline().lower()
  cell.direct = (l == "direct")
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

"""Returns a list [ [atom_1, atom_2, distance] ]"""
def GetPeriodicList(atoms, A, B, C):
  data = []
  for i,a in enumerate(atoms):
    for j in range(i+1,len(atoms)):
      b = atoms[j]
      dist = 1e10
      for ij in range(3):
        J = -1+ij
        for ik in range(3):
          K = -1+ik
          for il in range(3):
            L = -1+il
            v = J*A + K*B + L*C
            d1 = a.p.Distance(b.p + v)
            dist = min(d1,dist)
      data.append([i, j, dist])
  return data
"""Returns a list with two lists [ [[partner],[distance]] ]"""
def GetSortedList(li, atoms):
  data = [([],[]) for x in range(len(atoms))]
  for i in li:
    id1 = i[0]
    id2 = i[1]
    dist = i[2]
    
    data[id1][0].append(id2)
    data[id1][1].append(dist)
    
    data[id2][0].append(id1)
    data[id2][1].append(dist)
  # Data rearranged
  for i in data:
    swap = True
    while swap:
      swap = False
      for j in range(len(i[1])-1):
        a = i[1][j]; b = i[1][j+1]
        if b < a:
          swap = True
          i[1][j] = b; i[1][j+1] = a
          c = i[0][j]
          i[0][j] = i[0][j+1]; i[0][j+1] = c
  return data
  
"""Returns a list with a list of shells [ [(dist, count), (dist, count), ...] ]
Assumes the input list is sorted nearest to furthest"""
def GetShells(dis, atoms, tolerance):
  shells = [[] for x in dis]
  for i,x in enumerate(dis):
    indices = x[0]
    distances = x[1]
    last = distances[0]
    sums = distances[0]
    n = 1
    average = sums
    # Use a running average for the shells instead of looking at distance to next...
    for j in range(1,len(distances)):
      if (distances[j]-average) > tolerance: # Hit a new shell
        shells[i].append((sums/n, n))
        sums = distances[j]; n = 1
        average = sums
      else:
        n+=1; sums += distances[j]
        average = sums/float(n)
      last = distances[j]
    shells[i].append((average, n))
  return shells

def ToCartesian(atoms,A,B,C):
  for i in range(len(atoms)):
    v = atoms[i].p
    atoms[i].p = A*v.x + B*v.y + C*v.z

def PrintHelp():
  print("This script reads in two VASP CONTCAR/POSCAR to print out RMSD...")
  print("\tDefaults to compare POSCAR and CONTCAR in current folder")
  print("\tRMSD is calculated using cartesian coordinates, not fractional,")
  print("\ttherefore lattice mismatch should be accounted for...")
  print("Options: -a FILENAME -b FILENAME")
  print("\t-i Atom index starting from 0 of desired atoms to examine")
  print("\t-d Reports deviation of each atom individually along with system average")

def RMSD(args):
  a = "POSCAR"
  b = "CONTCAR"
  detail = False
  index = -1
  i = 1
  while i < len(args):
    if "-" == args[i][0] and len(args[i]) == 2:
      # Reading an option
      if args[i][1] == 'i':
        index = int(args[i+1])
        i += 2
      elif args[i][1] == 'd':
        detail = True
        i += 1
      elif args[i][1] == 'a':
        a = args[i+1]
        i += 2
      elif args[i][1] == 'b':
        b = args[i+1]
        i += 2
      else:
        PrintHelp()
        return
    else:
      PrintHelp()
      return
  A = ReadCONTCAR(a)
  B = ReadCONTCAR(b)
  if A.direct:
    lA = A.lattice[0] ; lB = A.lattice[1] ; lC = A.lattice[2];
    for a in A.atoms:
      a.p = a.p.x*lA + a.p.y*lB + a.p.z*lC
  if B.direct:
    lA = B.lattice[0] ; lB = B.lattice[1] ; lC = B.lattice[2];
    for a in B.atoms:
      a.p = a.p.x*lA + a.p.y*lB + a.p.z*lC
  if len(A.atoms) != len(B.atoms) or len(A.counts) != len(B.counts):
    print("Systems are not equivalent")
  bnAtoms = len(B.atoms)
  count = 0
  devTotal = 0.
  for i,a in enumerate(A.atoms):
    if index >= 0 and i != index:
      continue
    if i >= bnAtoms:
      break
    b = B.atoms[i]
    dp = (a.p - b.p).Length()
    if detail:
      print("{0:<2} {1:<3} Deviation: {2:<6.3f} A".format(
      a.name, i+1, dp))
    count += 1
    devTotal += dp
  if count == 0:
    print("No atoms...")
  else:
    print("RMSD: {0:>12.5e} A".format(devTotal/count))
  print("Good Day")
  return

if __name__ == "__main__":
  import sys
  RMSD(sys.argv)
