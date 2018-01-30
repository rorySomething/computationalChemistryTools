#!/bin/python3

from math import sqrt
from math import log

class vec3:
  x = y = z = 0.
  def __init__(self, x):
    self.x = float(x[0]); self.y = float(x[1]); self.z = float(x[2])
  def Distance(self, v):
    x = self.x - v.x; y = self.y - v.y; z = self.z - v.z
    return sqrt(x*x+y*y+z*z)
  def __mul__(self, x):
    return vec3([self.x*float(x),self.y*float(x),self.z*float(x)])
  def __rmul__(self, x):
    return self.__mul__(x)
  def __add__(self, x):
    return vec3([self.x+x.x,self.y+x.y,self.z+x.z])
  def __radd__(self, x):
    return self.__add__(x)
  def __str__(self):
    return "[{0}, {1}, {2}]".format(self.x, self.y, self.z)
  def Dot(self, b):
    return self.x*b.x + self.y*b.y + self.z*b.z
  def Cross(self, b):
    return vec3([self.y*b.z-self.z*b.y, self.x*b.z-self.z*b.x, self.x*b.y-self.y*b.x])
  #
  # [ i  j  k  ]
  # [ ax ay az ]
  # [ bx by bz ]
  #
  # { i*(ay*bz-az*by) - j*(ax*bz-az*bx) + k*(ax*by-bx*ay) }
  # i = (ay*bz-az*by); j = (ax*bz-az*bx); k = (ax*by-bx*ay)
    
class Atom:
  def __init__(self):
    self.p = vec3()
    self.element = 0
    self.name = ''
  def __init__(self, pos, element, name):
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
            d1 = self.p.Distance(i.p + v)
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
            d1 = self.p.Distance(atoms[l[i]].p + v)
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
      d.append(self.p.Distance(atoms[l[i]].p))
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
      
types = []
atoms = []
A = 0 ; B = 0 ; C = 0
def ReadCONTCAR(filename, verbose = True):
  try:
    f = open(filename, 'r')
  except:
    print("Can't open file", filename)
    return
  types.clear()
  atoms.clear()
  global A
  global B
  global C
  title = f.readline()
  labels = False
  scale = float(f.readline())
  l = f.readline().split()
  A = vec3(l) #vec3(float(l[0]), float(l[1]), float(l[2]))
  l = f.readline().split()
  B = vec3(l) #vec3(float(l[0]), float(l[1]), float(l[2]))
  l = f.readline().split()
  C = vec3(l) #vec3(float(l[0]), float(l[1]), float(l[2]))
  natoms = f.readline().split()
  total = 0
  try:
    n = int(natoms[0])
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
      types.append(int(i))
  else:
    for i in natoms:
      types.append(int(i))
  for i in types:
    total += i
  l = [x for x in types]
  print("Found {0} types and {1} total atoms:".format(len(types), total), l)
  direct = f.readline()
  if "Selective" in direct: # Selective Dynamics followed by Cartesian or Recip. Flag
    direct = f.readline()
  n = 0
  for i in range(len(types)):
    for j in range(types[i]):
      l = f.readline().split()
      if labels:
        atoms.append(Atom(vec3(l),1,natoms[i]))
      else:
        atoms.append(Atom(vec3(l),1,str(i+1)))
  #Ignore velocities
  f.close()

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
  print("This script reads VASP CONTCAR/POSCAR to print out nearest neighbor info...")
  print("Options: -t tolerance to call two bond lengths similar")
  print("\t-i Atom index starting from 0 of desired atom to examine")
  print("\t-c Print 'cleaner' output")

def Nearest_Neighbors(args):
  fn = "CONTCAR"
  tolerance = 0.03
  clean = False
  index = -1
  i = 1
  while i < len(args):
    if "-" == args[i][0]:
      # Reading an option
      if args[i][1] == 'i':
        index = int(args[i+1])
        i += 2
      elif args[i][1] == 't':
        tolerance = float(args[i+1])
        i += 2
      elif args[i][1] == 'c':
        clean = True
        i += 1
      else:
        PrintHelp()
        return
    else:
     fn = args[i]
     i += 1

  ReadCONTCAR(fn)
  print(len(atoms),"atoms")
  ToCartesian(atoms,A,B,C)
  with open("Nearest_Neighbor_Shells.txt",'w') as out:
    print(A,B,C,sep='\n')
    if index < 0:
      print("Calculating distances... (with periodicity)")
      li = GetPeriodicList(atoms, A, B, C)
      print("Sorting pairs...")
      dis = GetSortedList(li, atoms)
      print("Getting shells...")
      print("Distance grouping tolerance {0: <8.2e} Angstroms.".format(tolerance))
      s = GetShells(dis, atoms, tolerance)
      if clean:
        decimals = max(int(round(log(1./tolerance, 10)+0.5)),0) # Add 0.5 to guarantee round up
        digits = decimals + 4
      for i,a in enumerate(atoms):
        if clean:
          print("Shells: {0}".format(", ".join(
          ["({0:.{3}f}, {1})".format(x,y,digits,decimals) for x,y in s[i]])))
        else:
          print("Shells:", s[i])
#        out.write("Atoms {0}\nShells: {1}\n".format(dis[i][0],
#        ["%0.3f, %i" % (x,y) for x,y in s[i]]))
        out.write("Atoms {0}\nShells: {1}\n".format(dis[i][0],
        ", ".join(["{0:0.3f}:{1}".format(x,y) for x,y in s[i]])))
      checkMean(atoms,dis,li,s)
    else:
      for i,a in enumerate(atoms):
        if index >= 0 and i != index:
          continue
        l = a.SortedNeighborsPeriodic(atoms, A, B, C)
        d = a.DistancesPeriodic(l,atoms,A,B,C)
        s = a.ShellsPeriodic(l,atoms, tolerance,A,B,C)
        print(l, "\n\t", d,"\nShells:", s)
        out.write("Atoms {0}\nShells: {1}\n".format(l, s))
      checkMean(atoms,dis,li,s)
  print("Good Day")

if __name__ == "__main__":
  import sys
  Nearest_Neighbors(sys.argv)

