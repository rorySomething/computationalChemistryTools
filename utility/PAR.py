# Reads angles and such from CHARMM parameter files
###TODO
# optimize searches by sorting bonds, angle, etc by typename

class Parameters:
  def __init__(self):
    self.filenames = []
    self.atoms = dict()   # Name:Mass
    self.bonds = []       # [A,B,k,r]
    self.angles = []      # [A,B,C,k,theta]
    self.dihedrals = []   # [A,B,C,D,k,g,phi]
    self.impropers = []   # [A,B,C,D,k,phi]
    self.nonbonding = []  # [A, 0.0, sigma, r0]
  def ReadPAR(self, filename):
    self.filenames.append(filename)
    with open(filename, 'r') as f:
      bonds = False
      angles = False
      dihedrals = False
      impropers = False
      LJ = False
      for l in f:
        l = l.strip()
        if not l or l[0] == '!':
          continue # Comment
        elif "END" in l[:3]:
          break
        elif "HBOND" in l[:5]:
          continue
        elif "MASS" in l[:4]:
          l = l.split()
          if len(l) < 4:
            continue
          self.atoms[l[2]] = float(l[3])
        if "BONDS" in l[:6]:
          bonds = True
          angles = False
          dihedrals = False
          impropers = False
          LJ = False
        elif "ANGLES" in l[:10]:
          bonds = False
          angles = True
          dihedrals = False
          impropers = False
          LJ = False
        elif "DIHEDRALS" in l[:12]:
          bonds = False
          angles = False
          dihedrals = True
          impropers = False
          LJ = False
        elif "IMPROPER" in l[:12]:
          bonds = False
          angles = False
          dihedrals = False
          impropers = True
          LJ = False
        elif "NONBONDED" in l[:12]:
          bonds = False
          angles = False
          dihedrals = False
          impropers = False
          LJ = True
        if bonds:
          l = l.split()
          if len(l) < 4:
            continue
          self.bonds.append(l[:4]) # A,B,k,r
        elif angles:
          l = l.split()
          if len(l) < 5:
            continue
          self.angles.append(l[:5]) # A,B,C,k,theta
        elif dihedrals:
          l = l.split()
          if len(l) < 7:
            continue
          self.dihedrals.append(l[:7]) # A,B,C,D,k,g,phi
        elif impropers:
          l = l.split()
          if len(l) < 7:
            continue
          self.impropers.append(l[:7]) # A,B,C,D,k,g,phi <- phase g,phi should be 0
        elif LJ:
          l = l.split()
          if l[0] == 'cutnb' or len(l) < 4:
            continue
          self.nonbonding.append(l[:4]) # A, 0.0, sigma, r0
    f.close()
    print("Bonds",len(self.bonds))
    print("Angles",len(self.angles))
    print("Dihedrals",len(self.dihedrals))
  
  # Type names A,B
  def getBondData(self,A,B):
    for b in self.bonds:
      if (A == b[0] and B == b[1]) or (A == b[1] and B == b[0]):
        return b
    return False
  
  # Type names A,B
  def getBondLength(self,A,B):
    b = self.getBondData(A,B)
    if b:
      return b[-1]
    else:
      return b # False
  
  def getBondForceConstant(self,A,B):
    b = self.getBondData(A,B)
    if b:
      return b[-2]
    else:
      return b # False
  
  # Type names A,B,C
  def getAngleData(self,A,B,C):
    for a in self.angles:
      if B == a[1] and (A == a[0] and C == a[2]) or (A == a[2] and C == a[0]):
        return a
    return False
  
  # Type names A,B,C
  def getAngle(self,A,B,C):
    a = self.getAngleData(A,B,C)
    if a:
      return a[-1]
    else:
      return a # False
  
  def getAngleForceConstant(self,A,B,C):
    a = self.getAngleData(A,B,C)
    if a:
      return a[-2]
    else:
      return a # False
  
  
  # Type names A,B,C,D
  def getDihedralData(self,A,B,C,D):
    for a in self.dihedrals:
      if (A == a[0] and B == a[1] and C == a[2] and D == a[3]) or (A == a[3] and B == a[2] and C == a[1] and D == a[0]):
        return a
    return False
  
  # Type names A,B,C,D
  def getDihedralAngle(self,A,B,C,D):
    a = self.getDihedralData(A,B,C,D)
    if a:
      return a[-1]
    else:
      return a # False
  
  def getDihedralForceConstant(self,A,B,C,D):
    a = self.getDihedralData(A,B,C,D)
    if a:
      return a[-3]
    else:
      return a # False
  def getDihedralMult(self,A,B,C,D):
    a = self.getDihedralData(A,B,C,D)
    if a:
      return a[-2]
    else:
      return a # False


  # Type names A,B,C,D
  # A is center, bonded to B,C and D
  # Angle is between planes ABC and BCD
  def getImproperList(self,A):
    imps = []
    for a in self.impropers:
      if A in a[:4]:
        imps.append(a[:])
    return imps[:]

  def getImproperData(self,A,B,C,D):
    for a in self.impropers:
      if (A == a[0] and B == a[1] and C == a[2] and D == a[3]) or (A == a[3] and B == a[2] and C == a[1] and D == a[0]):
        return a
    return False
  
  # Type names A,B,C,D
  def getImproperAngle(self,A,B,C,D):
    a = self.getImproperData(A,B,C,D)
    if a:
      return a[-3]
    else:
      return a # False
  
  # SHOULD BE 0
  def getImproperPhase(self,A,B,C,D):
    a = self.getImproperData(A,B,C,D)
    if a:
      return a[-1]
    else:
      return a # False
    
  # SHOULD BE 0
  def getImproperMult(self,A,B,C,D):
    a = self.getImproperData(A,B,C,D)
    if a:
      return a[-2]
    else:
      return a # False