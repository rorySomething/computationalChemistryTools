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
  def __init__(self):
    self.p = vec3()
    self.element = 0
    self.name = ''
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
def LatticeCheck(cell1, scale, cell2, verbose = False):
  A = [x*cell1.scale for x in cell1.lattice]; B = [x*cell2.scale for x in cell2.lattice]
  for i in range(len(scale)):
    A[i] = A[i]/scale[i]
  if (A[0]-B[0]).Length() < 0.0001:
    if (A[1]-B[1]).Length() < 0.0001:
      if (A[2]-B[2]).Length() < 0.0001:
        return True
      else:
        if verbose:
          print("Lattice mismatch, C vectors don't match.")
        return False
    else:
      if verbose:
        print("Lattice mismatch, B vectors don't match.")
      return False
  else:
    if verbose:
      print("Lattice mismatch, A vectors don't match.")
      print(A[0],B[0],":",(A[0]-B[0]).Length())
    return False
def GetSuperCell(cell, scale = None, crystal = None, outscale = None):
  if not outscale:
    return cell
  outcell = CrystalCell()
  outcell.title = "supercell from " + cell.title
  outcell.lattice = []
  outcell.atoms = []
  outcell.scale = cell.scale
  outcell.direct = cell.direct
  outcell.labeled = cell.labeled
  a = cell.lattice[0]; b = cell.lattice[1]; c = cell.lattice[2]
  print('L:',[str(x) for x in cell.lattice],'A,B,C:',a,b,c)
  outcell.counts = []
  outcell.types = []
  for i in range(len(cell.types)):
    outcell.counts.append(0)
    outcell.types.append(cell.types[i])
  if crystal and scale:
    ca = crystal.lattice[0]; cb = crystal.lattice[1]; cc = crystal.lattice[2]
    outcell.lattice.append(a/scale[0]*outscale[0])
    outcell.lattice.append(b/scale[1]*outscale[1])
    outcell.lattice.append(c/scale[2]*outscale[2])
    cartAtoms = [x for x in cell.atoms] # Leave fractional
    for x in cartAtoms:
#      x.p = a*x.p.x + b*x.p.y + c*x.p.z # To cartesian coordinates
       x.p.x = x.p.x * scale[0] / outscale[0]
       x.p.y = x.p.y * scale[1] / outscale[1]
       x.p.z = x.p.z * scale[2] / outscale[2]
    for i in range(outscale[0]):
      for j in range(outscale[1]):
        for k in range(outscale[2]): # Add all new crystal atoms
          if i < scale[0] and j < scale[1] and k < scale[2]:
            continue
          for l in crystal.atoms:
#            p = l.p.x * ca + l.p.y * cb + l.p.z * cc # To cartesian
            x = l.p.x / outscale[0] + i/ outscale[0]
            y = l.p.y / outscale[1] + j/ outscale[1]
            z = l.p.z / outscale[2] + k/ outscale[2]
            p = vec3([x,y,z])
            m = Atom(p, l.element, l.name)
#            p += ca * i + cb * j + cc * k # Move by periodicity
#            m = Atom(p, l.element, l.name)
            cartAtoms.append(m) # Add to new system
    a = outcell.lattice[0]
    b = outcell.lattice[1]
    c = outcell.lattice[2]
    for i in cartAtoms: # Change to fractional of new system
 #     outcell.atoms.append(Atom(i.p.Fractional(a,b,c), i.element, i.name))
      outcell.atoms.append(Atom(i.p, i.element, i.name))
      for j in range(len(outcell.types)):
        if i.name == outcell.types[j]:
          outcell.counts[j] += 1
  else:
    # Just make a supercell
    outcell.lattice.append(a*outscale[0])
    outcell.lattice.append(b*outscale[1])
    outcell.lattice.append(c*outscale[2])
    for i in range(outscale[0]):
      for j in range(outscale[1]):
        for k in range(outscale[2]): # Add all crystal atoms over supercell
          for l in cell.atoms:
            x = l.p.x / outscale[0] + i/ outscale[0]
            y = l.p.y / outscale[1] + j/ outscale[1]
            z = l.p.z / outscale[2] + k/ outscale[2]
            p = vec3([x,y,z])
            m = Atom(p, l.element, l.name)
            outcell.atoms.append(m) # Add to new system
            outcell.counts[l.element] += 1
  return outcell
    
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
  for l in cell.atoms:
    f.write("{0.x:>20.16f}{0.y:>20.16f}{0.z:>20.16f}\n".format(l.p))
  f.write('\n')
  for l in cell.atoms:
    f.write("{0:>16.8E}{1:>16.8E}{2:>16.8E}\n".format(0.0,0.0,0.0))
  f.close()

def Reflect(args):
  i = 1
  read = "CONTCAR"
  out = ""
  indices = []
  tolerance = 0.05
  floor = False
  while i < len(args):
    a = args[i]
    if len(a) > 1 and a[0] == '-':
      if a[1] == 'f':
        read = args[i+1]
        i += 1
      elif a[1] == 'o':
        out = args[i+1]
        i += 1
      elif a[1] == 't':
        tolerance = float(args[i+1])
        i += 1
      elif a[1] == 'n':
        floor = True
      elif a[1] == 'h' or a[1:3] == '-h':
        print('Reflects atoms along lattice vectors if within a distance of the walls.')
        print('\tUsage: reflect.py 0,1,5,7,9')
        print('Default behavior is to read from CONTCAR to create cell')
        print('Default tolerance is 0.05')
        print('e.g. If a specified atom is at fractional coordinates 0.0 0.5 0.96')
        print('\tIt will be moved to position 0.0 0.5 -0.04')
        print('Intended for use with supercell.py after a defect pushes atoms through a crystal lattice.')
        print('Output is a POSCAR file named REFLECT')
        print('Input atom indices from 0 separated by commas.')
        print('Options: -f, -o, -c')
        print('\t-f - Input file name (default is CONTCAR)')
        print('\t-o - Output file name (default is REFLECT)')
        print('\t-t - Fractional coordinate tolerance (default is 0.05)')
        print('\t-n - Floor all of the coordinates, only reflect coordinates > 1.0 - tolerance')
        print('\t\tThis is the best option to work nicely with supercell.py and the initial bulk')
        return
    else:
      if len(indices) == 0:
        try:
          indices.append(int(a))
        except:
          indices = a.split(',')
    i += 1
  if not out:
    out = 'REFLECT'
  # Start doing things
  print("Read Cell")
  cell = ReadCONTCAR(read, True)
  print(len(cell.atoms),"atoms")
  if len(indices) == 0:
    for i in range(len(cell.atoms)):
      indices.append(i)
  for i in indices:
    A = cell.atoms[int(i)]
    if A.p.x < tolerance and not floor:
      A.p.x += 1.0
    elif A.p.x > 1.0 - tolerance:
      A.p.x -= 1.0
    if A.p.y < tolerance and not floor:
      A.p.y += 1.0
    elif A.p.y > 1.0 - tolerance:
      A.p.y -= 1.0
    if A.p.z < tolerance and not floor:
      A.p.z += 1.0
    elif A.p.z > 1.0 - tolerance:
      A.p.z -= 1.0
  WritePOSCAR(cell, out)
  print("Good Day")

if __name__ == "__main__":
  import sys
  Reflect(sys.argv)

