
# TODO clean up more to share between GAMESS, VASP, NAMD, & LAMMPS scripts

from .vec3 import vec3

lastID = 1
class Atom:
  def __init__(self, pos = vec3([0.,0.,0.]), atomicNumber = 0., name = ""):
    # GAMESS
    # pos, atomicNumber, name
    # PDB
    # index, name, resname, resid, chain, x,y,z
    self.p = vec3(pos.getList())
    self.atomicNumber = atomicNumber
    self.element = ""
    self.name = name
    self.id = 0
    self.resName = ""
    self.resID = 0
    self.typeName = ""
    self.charge = 0.
    self.chainID = ""
    self.chain = " "
    self.occupancy = 1.0
    self.tempFactor = 1.0
  def PDBLine(self):
    # Specific output from VMD 1.9.1 for easy diff use
    if len(self.name) < 4:
      return "ATOM {0:>6}  {1:<4}{2:<4}{3}{4:>4}    \
{5:>8.3f}{6:>8.3f}{7:>8.3f}{8:6.2f}{9:6.2f}      {10:<4}{11:>2}\n".format(
self.id,self.name,self.resName,self.chainID,self.resID,self.p.x,self.p.y,self.p.z,self.occupancy,self.tempFactor,
self.chain,self.element)
    else:
      return "ATOM {0:>6} {1:<5}{2:<4}{3}{4:>4}    \
{5:>8.3f}{6:>8.3f}{7:>8.3f}{8:6.2f}{9:6.2f}      {10:<4}{11:>2}\n".format(
self.id,self.name,self.resName,self.chainID,self.resID,self.p.x,self.p.y,self.p.z,self.occupancy,self.tempFactor,
self.chain,self.element)
  
  def readPDBLine(self, line):
    global lastID
    if l[:4] != "ATOM":
      return False
    try:
      self.id = int(l[4:11])
    except:
      try:
        self.id = int(l[4:11], 16)
      except:
        try:
          if '*' in l[4:11]:
            self.id = lastID
        except:
          print("Error reading ATOM entry:",l)
          print("\tSkipping...")
          return False
    self.name = l[11:16].strip()
    self.resName = l[16:21].strip()
    self.chainID = l[21].strip()
    self.resID = int(l[22:26])
    self.x = float(l[30:38])
    self.y = float(l[38:46])
    self.z = float(l[46:54])
    self.occupancy = float(l[54:60])
    self.tempFactor = float(l[60:66])
    self.chain = l[72:76].strip()
    self.element = l[76:].strip()
    if self.element == "":
      self.element = " "
    lastID += 1
    return True

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
  def distance(self, b):
    return self.p.distance(b.p)
    
  def distances(self, l, atoms):
    d = []
    for i in range(len(l)):
      d.append(self.p.distance(atoms[l[i]].p))
    return d
    
  def Shells(self, l, atoms, tol):
    d = self.distances(l, atoms)
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
