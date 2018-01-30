#!/bin/python3

from math import sqrt
from math import log
from time import clock
import numpy as np

class vec3:
  def __init__(self, x):
    self.x = float(x[0]); self.y = float(x[1]); self.z = float(x[2])
  def distanceSquared(self, v):
    x = self.x - v.x; y = self.y - v.y; z = self.z - v.z
    return x*x + y*y + z*z
  def distance(self, v):
    return sqrt(self.distanceSquared(v))
  def lengthSquared(self):
    return self.dot(self)
  def length(self):
    return sqrt(self.dot(self))
  def __mul__(self, x):
    return vec3([self.x*float(x),self.y*float(x),self.z*float(x)])
  def __rmul__(self, x):
    return self.__mul__(x)
  def __add__(self, x):
    return vec3([self.x+x.x,self.y+x.y,self.z+x.z])
  def __sub__(self, x):
    return vec3([self.x-x.x,self.y-x.y,self.z-x.z])
  def __rsub__(self, x):
    return vec3([-self.x+x.x,-self.y+x.y,-self.z+x.z])
  def __radd__(self, x):
    return self.__add__(x)
  def __str__(self):
    return "[{0}, {1}, {2}]".format(self.x, self.y, self.z)
  def getList(self):
    return [self.x,self.y,self.z]
  def normalize(self):
    l = self.length()
    return vec3([self.x/l, self.y/l, self.z/l])
  def dot(self, b):
    return self.x*b.x + self.y*b.y + self.z*b.z
  def vectorMultiply(self, b):
    return vec3([self.x*b.x, self.y*b.y, self.z*b.z])
  def Cross(self, b):
    return vec3([self.y*b.z-self.z*b.y, self.x*b.z-self.z*b.x, self.x*b.y-self.y*b.x])
    #
    # [ i  j  k  ]
    # [ ax ay az ]
    # [ bx by bz ]
    #
    # { i*(ay*bz-az*by) - j*(ax*bz-az*bx) + k*(ax*by-bx*ay) }
    # i = (ay*bz-az*by); j = (ax*bz-az*bx); k = (ax*by-bx*ay)
  xyz = property(getList, None)
    
class Atom:
  def __init__(self, pos, element = 0, name = ''):
    self.p = pos
    self.element = element
    self.name = name
  def toFractional(self, A,B,C):
    cart = np.array(self.p.xyz)
    m = np.matrix([A.xyz, B.xyz, C.xyz])
    frac = cart * m.I
    self.p = vec3(frac.tolist()[0])
  def toCartesian(self, A,B,C):
    frac = np.array(self.p.xyz)
    m = np.matrix([A.xyz, B.xyz, C.xyz])
    cart = frac * m
    self.p = vec3(cart.tolist()[0])
  def SortedNeighbors(self, atoms):
    d = []
    l = []
    n = 0
    for i in atoms:
      d.append(self.p.distance(i.p))
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
  def SortedNeighborsPeriodic(self, atoms, A, B, C):
    d = []
    l = []
    lAtom = []
    n = 0
    for i in atoms:
      dist = 1e10
      for ij in range(3):
        J = -1+ij
        for ik in range(3):
          K = -1+ik
          for il in range(3):
            L = -1+il
            v = J*A + K*B + L*C
            d1 = self.p.distance(i.p + v)
            dist = min(d1,dist)
      d.append(dist)
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
  def DistancesPeriodic(self, l, atoms, A, B, C):
    d = []
    for i in range(len(l)):
      dist = 1e10
      for ij in range(3):
        J = -1+ij
        for ik in range(3):
          K = -1+ik
          for il in range(3):
            L = -1+il
            v = J*A + K*B + L*C
            d1 = self.p.distance(atoms[l[i]].p + v)
            dist = min(d1,dist)
      d.append(dist)
    return d
  def ShellsPeriodic(self, l, atoms, tol, A, B, C):
    d = self.DistancesPeriodic(l, atoms, A, B, C)
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
    
  def Distances(self, l, atoms):
    d = []
    for i in range(len(l)):
      d.append(self.p.distance(atoms[l[i]].p))
    return d
    
  def Shells(self, l, atoms, tol):
    d = self.Distances(l, atoms)
    last = d[1]
    sums = d[1]
    n = 1
    average = sums
    s = []
    for i in range(1,len(d)):
      if (d[i]-average) > tol:
        s.append(average)
        s.append(n); sums = d[i]; n = 1
        average = sums
      else:
        n+=1; sums += d[i]
        average = sums/float(n)
      last = d[i]
    if n != 0:
      s.append(average)
      s.append(n);
    return s

""" Returns atoms[], types[], counts[], A,B,C, title, labels=True|False, direct="Direct"|"Cartesian"
"""
def ReadCONTCAR(filename, verbose = True):
  try:
    f = open(filename, 'r')
  except:
    print("Can't open file", filename)
    return
  types = []
  counts = []
  atoms = []
  A, B, C = 0, 0, 0
  title = f.readline()
  labels = False
  scale = float(f.readline())
  l = f.readline().split()
  A = vec3(l) * scale  #vec3(float(l[0]), float(l[1]), float(l[2]))
  l = f.readline().split()
  B = vec3(l) * scale  #vec3(float(l[0]), float(l[1]), float(l[2]))
  l = f.readline().split()
  C = vec3(l) * scale  #vec3(float(l[0]), float(l[1]), float(l[2]))
  types = f.readline().split()
  total = 0
  try:
    n = int(types[0])
  except:
    print("Elements are labeled?")
    labels = True
  if verbose:
    print("Title:",title)
    print("Lattice:\n\t{0:4.3f}\n\t{1.x} {1.y} {1.z}\
    \n\t{2.x} {2.y} {2.z}\n\t{3.x} {3.y} {3.z}".format(scale,A,B,C))
    if labels:
      print("Labeled atoms.")
    else:
      print("VASP 4.6 format.")
  if labels:
    l = f.readline().split()
    for i in l:
      counts.append(int(i))
  else:
    for i in natoms:
      counts.append(int(i))
  for i in counts:
    total += i
  l = [x for x in types]
  print("Found {0} types and {1} total atoms:".format(len(types), total), l)
  direct = f.readline()
  if "Selective" in direct: # Selective Dynamics followed by Cartesian or Recip. Flag
    direct = f.readline()
  n = 0
  for i in range(len(counts)):
    for j in range(counts[i]):
      l = f.readline().split()
      if labels:
        atoms.append(Atom(vec3(l),1,types[i]))
      else:
        atoms.append(Atom(vec3(l),1,str(i+1)))
  #Ignore velocities
  f.close()
  return (atoms[:], types[:], counts[:], A, B, C, title, labels, direct)

class CrystalCell:
  def __init__(self):
    self.atoms = []
    self.types = []
    self.counts = []
    self.scale = 1.
    self.direct = False
    self.selectiveDynamics = False
    self.lattice = [None, None, None]
    self.title = ""
    self.labeled = False
  def addAtom(self, atom):
    if atom.name not in self.types:
      self.types.append(atom.name)
      self.counts.append(0)
    self.counts[self.types.index(atom.name)] += 1
    self.atoms.append(atom)
  def addAtoms(self, atomList):
    for a in atomList:
      self.addAtom(a)
  def setA(self, A):
    #TODO Add a test for vec3 type vs. [x,y,z] type
    self.lattice[0] = A
  def setB(self, B):
    self.lattice[1] = B
  def setC(self, C):
    self.lattice[2] = C
  def checkLabels(self):
    self.labeled = True
    if len(self.atoms) == 0:
      self.labeled = False
    for t in self.types:
      if not t.isalpha():
        self.labeled = False  

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

"""Returns a list [ [atom_1, atom_2, distance] ]
Works in cartesian coordinates???"""
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
            d1 = a.p.distance(b.p + v)
            dist = min(d1,dist)
      data.append([i, j, dist])
  print(data[:217])
  return data
"""Returns a list with two lists [ [partner],[distance] ]"""
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

def checkMean(atoms,distances,neighborList,shells):
  shellDict = {}
  for i,a in enumerate(atoms):
    firstShell = shells[i][0]
    if a.name in shellDict.keys():
      shellDict[a.name][0] += float(firstShell[0])
      shellDict[a.name][1] += int(firstShell[1])
      shellDict[a.name][2] += 1
    else:
      shellDict[a.name] = []
      shellDict[a.name].append(float(firstShell[0])) # Distance
      shellDict[a.name].append(int(firstShell[1])) # Neighbor Count
      shellDict[a.name].append(1) # Count
  for k,v in shellDict.items():
    print("Atom {0} mean 1st shell {1: 5.3f}:{2:>2}".format(
      k, v[0]/float(v[2]), v[1]/int(v[2])))
    v.append(0.) # Distance Variance
    v.append(0.) # Count Variance
  for i,a in enumerate(atoms):
    firstShell = shells[i][0]
    totalCount = shellDict[a.name][2]
    meanDist = shellDict[a.name][0] / totalCount
    meanNeighbors = shellDict[a.name][1] / totalCount
    shellDict[a.name][-2] += ( (firstShell[0] - meanDist)**2 ) / totalCount # Dist Variance
    shellDict[a.name][-1] += ( (firstShell[1] - meanNeighbors)**2 ) / totalCount # Count Variance
  for k,v in shellDict.items():
    print("Atom {0} 1st shell Mean:Std.Dev. {1}:{2} Distance - {3}:{4} Neighbors".format(
      k, v[0]/v[2], sqrt(v[-2]), v[1]/v[2], sqrt(v[-1])))
    
  

def PrintHelp():
  print("This script will trim down and center a lattice structure.\n\
  Reads VASP CONTCAR/POSCAR, select atom index to center around and\n\
  a total number of atoms to trim down to...\n\
  ATOMS ARE INDEXED FROM 0\n\
  centerAndTrim.py index total filename")
  print("Options:\n\
  \t-c 0,0,0 Coordinates to center the atoms to\n\
  \t-l a,b,c New lattice constants")
  print("Example:\n\t\
  centerAndTrim.py 0 15 CONTCAR\n\t\tTrim down to 15 atoms around and including the first\n\
  \t\tindexed from 0")

def CenterAndTrim(args):
  fn = "CONTCAR"
  index = -1
  fatoms = -1
  center = [0.,0.,0.]
  flattice = None
  i = 1
  if len(args) < 3:
    PrintHelp()
    return
  while i < len(args):
    arg = args[i]
    if "-" == arg[0]:
      if arg[1] == 'c' and len(arg) == 2:
        # Specifying center?
        arg2 = args[i+1]
        if len(arg2.split()) == 3:
          try:
            cent = [float(x) for x in args.split()]
            center = cent
            i += 2
            continue
          except:
            print("Bad input arguments? -c x,y,z is expected")
            raise
        else:
          try:
            center[0] = float(arg2)
            center[1] = float(args[i+2])
            center[2] = float(args[i+3])
            i += 4
            continue
          except:
            print("Bad input arguments? -c x,y,z is expected or -c x y z")
            raise
      elif arg[1] == 'l' and len(arg) == 2:
        # Specifying final lattice?
        arg2 = args[i+1]
        flattice = [0.,0.,0.]
        if len(arg2.split(',')) == 3:
          try:
            lat = [float(x) for x in arg2.split(',')]
            flattice = lat
            i += 2
            continue
          except:
            print("Bad input arguments? -l a,b,c is expected e.g. -l 5.2,5.2,4.7")
            raise
        else:
          try:
            flattice[0] = float(arg2)
            flattice[1] = float(args[i+2])
            flattice[2] = float(args[i+3])
            i += 4
            continue
          except:
            print("Bad input arguments? -l a,b,c is expected or -l a b c    e.g. -l 5.2 5.2 4.7")
            raise
      elif (arg[1] == 'h' or '-h' in arg) and len(arg) < 7: # -h, -help, --help
        PrintHelp()
        return
    else:
      if index < 0:
        index = int(arg)
        i += 1
      elif fatoms < 0:
        fatoms = int(arg)
        i += 1
      else:
        fn = arg
        i += 1
  center = vec3(center)
  atoms, types, counts, A,B,C,title,labels,direct = ReadCONTCAR(fn)
  print(len(atoms),"atoms\nLattice:")
  print(A,B,C,sep='\n')
  #ToCartesian(atoms,A,B,C)
  for a in atoms:
    a.toCartesian(A,B,C)
  print("Getting periodic neighbor list")
  li = GetPeriodicList(atoms, A, B, C)
  dist = GetSortedList(li, atoms)
  #s = GetShells(dis, atoms, tolerance)
  for a in atoms:
    a.toFractional(A,B,C)
  translation = center - atoms[index].p # Vector to translate all to the center
  # dist = [ [ [neighbor], [distance] ] <- per atom index ]
  atomsToKeep = [index]
  if fatoms > 1:
    atomsToKeep.extend(dist[index][0][:fatoms-1])
  print("Atoms to keep",",".join([str(x) for x in atomsToKeep]))
  fatoms = []
  cell = CrystalCell()
  ascale = 1.; bscale = 1.; cscale = 1.
  nA = A; nB = B; nC = C
  if flattice:
    nA = A.normalize() * flattice[0]
    nB = B.normalize() * flattice[1]
    nC = C.normalize() * flattice[2]
  for i in atomsToKeep:
    a = atoms[i]
    a.p = (a.p + translation)
    a.toCartesian(A,B,C)
    a.toFractional(nA,nB,nC)
    cell.addAtom(a)
  cell.setA(nA); cell.setB(nB); cell.setC(nC);
  cell.checkLabels()
  cell.direct = "Direct" in direct
  cell.title = title.strip('\n')+" trimmed"
  WritePOSCAR(cell, fn + "Trimmed")
  print("Good Day")
  return

if __name__ == "__main__":
  import sys
  CenterAndTrim(sys.argv)

