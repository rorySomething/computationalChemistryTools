#!/bin/python3

import os
import math
class vec3:
  x = y = z = 0.
  def __init__(self, x = [0,0,0]):
    self.x = float(x[0]); self.y = float(x[1]); self.z = float(x[2])
  def Distance(self, v):
    x = self.x - v.x; y = self.y - v.y; z = self.z - v.z
    return math.sqrt(x*x+y*y+z*z)
  def Length(self):
    return math.sqrt(self.x*self.x+self.y*self.y+self.z*self.z)
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
  def AngleDegrees(self, b):
    return math.degrees( math.acos( self.Dot(b)/(self.Length()*b.Length()) ) )
    
class KPOINT:
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

class EIGENFILE:
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

ef = EIGENFILE()

efermi = XC = alpha_beta = 0

def ReadOUTCAR():
  try:
    f = open("OUTCAR")
  except:
    print("Can't open OUTCAR file")
    return
  slic = len(" E-fermi :")
  for l in f:
    if " E-fermi :" == l[:slic]:
      print("Found fermi level info")
      a = l.split()
      efermi = float(a[2])
      XC = float(a[4])
      alpha_beta = float(a[6][1:])
      print("E-fermi: {0: 5.3f}, XC(G=0): {1: 5.3f}, alpha+beta: {2: 5.3f}".
      format(efermi, XC, alpha_beta))
      break
  f.close()

def ReadKPoints(f, bands):
  f.readline() # Blank line
  a = f.readline().split()
  k = KPOINT()
  k.bands.clear()
  print(a)
  k.x = float(a[0])
  k.y = float(a[1])
  k.z = float(a[2])
  k.w = float(a[3])
  if ef.polarized:
    k.bands = [[],[]]
  for i in range(bands):
    if not ef.polarized:
      k.bands.append(float(f.readline().split()[1]))
    else:
      a = f.readline().split()
      k.bands[0].append(float(a[1]))
      k.bands[1].append(float(a[2]))
  return k

def Read_Eigenval():
  try:
    f = open("EIGENVAL")
  except:
    print("Can't open EIGENVAL file")
    return
  ef.l1 = f.readline().split()
  sp = ef.l1[3]
  ef.polarized = (2 == int(sp))
  ef.l2 = f.readline().split()
  ef.a3 = float(f.readline())
  ef.CAR = f.readline()
  ef.title = f.readline()
  a = f.readline().split()
  ef.valence = int(a[0])
  ef.kpts = int(a[1])
  ef.bands = int(a[2])
  for i in range(ef.kpts):
    ef.kpoints.append(ReadKPoints(f, ef.bands))
  f.close()
  ef.Print() # To confirm good read
  odd = ef.valence % 2 == 1
  homo = int(ef.valence/2)
  lumo = int(homo+1)
  print("VB, CB: {0}, {1}".format(homo,lumo))
  if odd:
    print("Odd number of electrons")
    if not ef.polarized:
      print("Calculation is closed shell with an odd number of electrons?")
      pass
    hAlpha = int((ef.valence-1)/2)+1
    lAlpha = hAlpha + 1
    hBeta = lAlpha
    lBeta = hBeta + 1
    print("Alpha VB, CB: {0}, {1}".format(hAlpha,lAlpha))
  homo -= 1; lumo -= 1 # Index from 0
  hmax = -1e10; hmin = 1e10
  hmaxkp = hminkp = 0
  lmax = -1e10; lmin = 1e10
  lmaxkp = lminkp = 0
  if odd:
    hAlpha -=1; lAlpha -= 1; hBeta -= 1; lBeta -= 1; #Index
    hamax = -1e10; hamin = 1e10; hamaxkp = haminkp = 0
    hbmax = -1e10; hbmin = 1e10; hbmaxkp = hbminkp = 0
    lamax = -1e10; lamin = 1e10; lamaxkp = laminkp = 0
    lbmax = -1e10; lbmin = 1e10; lbmaxkp = lbminkp = 0
  for i in range(len(ef.kpoints)):
    if not ef.polarized:
      vb = ef.kpoints[i].bands[homo]
      cb = ef.kpoints[i].bands[lumo]
    else:
      vb = ef.kpoints[i].bands[0][homo]
      cb = ef.kpoints[i].bands[0][lumo]
    if odd and ef.polarized:
      vba = ef.kpoints[i].bands[0][hAlpha]
      vbb = ef.kpoints[i].bands[1][hBeta]
      cba = ef.kpoints[i].bands[0][lAlpha]
      cbb = ef.kpoints[i].bands[1][lBeta]
      #Alpha
      if vba > hamax:
        hamax = vba; hamaxkp = i
      if vba < hamin:
        hamin = vba; haminkp = i
      if cba > lamax:
        lamax = cba; lamaxkp = i
      if cba < lamin:
        lamin = cba; laminkp = i
      # Beta
      if vbb > hbmax:
        hbmax = vbb; hbmaxkp = i
      if vbb < hbmin:
        hbmin = vbb; hbminkp = i
      if cbb > lbmax:
        lbmax = cbb; lbmaxkp = i
      if cbb < lbmin:
        lbmin = cbb; lbminkp = i
    # End of odd else
    if vb > hmax:
      hmax = vb; hmaxkp = i
    if vb < hmin:
      hmin = vb; hminkp = i
    if cb > lmax:
      lmax = cb; lmaxkp = i
    if cb < lmin:
      lmin = cb; lminkp = i
  if odd and ef.polarized:
    if hamaxkp == laminkp:
      print("Direct band gap (alpha) at", ef.kpoints[hamaxkp].Coord())
    else:
      print("Indirect band gap (alpha)")
    print("Gap (alpha):",'{:5.3f}'.format(lamin-hamax),"eV")
    if hbmaxkp == lbminkp:
      print("Direct band gap (beta) at", ef.kpoints[hbmaxkp].Coord())
    else:
      print("Indirect band gap (beta)")
    print("Gap  (beta):",'{:5.3f}'.format(lbmin-hbmax),"eV")
    print("Forbidden alpha -> beta:",'{:5.3f}'.format(lbmin-hamax),"eV")
    print("Forbidden beta -> alpha:",'{:5.3f}'.format(lamin-hbmax),"eV")
    print("Alpha VBE, CBE, Gap, VBW, CBW: \
    {0:5.3f}, {1:5.3f}, {2:5.3f}, {3:5.3f}, {4:5.3f} eV".format(hamax,lamin,
    lamin-hamax,hamax-hamin,lamax-lamin))
    print("Beta  VBE, CBE, Gap, VBW, CBW: \
    {0:5.3f}, {1:5.3f}, {2:5.3f}, {3:5.3f}, {4:5.3f} eV".format(hbmax,lbmin,
    lbmin-hbmax,hbmax-hbmin,lbmax-lbmin))
  else:
    if hmaxkp == lminkp:
      print("Direct band gap at", ef.kpoints[hmaxkp].Coord())
    else:
      print("Indirect band gap")
    print("Gap:",'{:5.3f}'.format(lmin-hmax),"eV")
    print("VBE, CBE, Gap, VBW, CBW: \
    {0:5.3f}, {1:5.3f}, {2:5.3f}, {3:5.3f}, {4:5.3f} eV".format(hmax,lmin,
    lmin-hmax,hmax-hmin,lmax-lmin))
  print("Read Outcar?")
  ReadOUTCAR()

class CalcInfo:
  def __init__(self):
    self.name = ""          # Identifier
    self.complete = False       # Run completed
    self.lastEnergy = 0.        # Energy without entropy
    self.title = ""             # Title of calc.
    self.atomCounts = []        # Atom sets from CONTCAR
    self.nAtoms = 0             # Total Atoms
    self.atoms = []             # Atoms
    self.Emp = ""               # Empirical Formula
    self.jobID = ""             # Job number
    self.elapsedTime = 0.       # Run time (Elapsed in OUTCAR)
    self.A = vec3()             # Lattice vector
    self.B = vec3()             # Lattice vector
    self.C = vec3()             # Lattice vector
    self.kpoints = []           # k-points and band data
    self.nKpoints = 0           # NKPT in OUTCAR
    self.polarized = False      # Spin polarized calculation ?
    self.bands = []             # Band data
    self.nBands = 0             # Bands computed
    self.bandStruct = False     # Band structure calculation ?
    self.bandComplete = False   # Band structure calculation completed ?
    self.nElectrons = 0         # NELECT
    self.iVBE = 0               # Valence band index
    self.iCBE = 0               # Conduction band index
    self.directGap = False      # Direct band gap ?
    self.bandGap = 0.           # Band gap energy
    self.directBandGap = 0.     # Direct gap energy
    self.dir = ""               # Calculation directory

def GetPoscarData(poscar, calcInfo):
  sum_ = 0
  with open(poscar, 'r') as f:
    calcInfo.title = f.readline() # Title
    f.readline() # Scale
    f.readline() # A
    f.readline() # B
    f.readline() # C
    a = f.readline().split() # Names?
    try:
      a[0] = int(a[0])
    except:
      # Caught atom names
      a = f.readline().split() # Names?
    calcInfo.atomCounts = []
    for i in a:
      calcInfo.atomCounts.append(int(i))
      sum_ += int(i)
  calcInfo.nAtoms = sum_
  
'''Return a string "directory_name Empirical_Formula Atom_Count Title"'''
def GetPotcarData(poscar, calcInfo):
  potcar = poscar[:-4]+'T'+poscar[-3:] # For calculations with 1 POSCAR and 1 POTCAR
  #potcar[-4] = 'T' # POSCAR to POTCAR
  GetPoscarData(poscar, calcInfo)
  calcInfo.atoms = []
  nAtoms = len(calcInfo.atomCounts)
  with open(potcar, 'r') as f:
    while True:
      l = f.readline()
      if not l:
        break
      if "VRHFIN" in l:
        a = l.split()[1][1:]
        if a[-1] == ':':
          a = a[:-1]
        calcInfo.atoms.append(a)
      if len(calcInfo.atoms) == nAtoms:
        break
  calcInfo.Emp = ""
  for i in range(len(calcInfo.atoms)):
    calcInfo.Emp += calcInfo.atoms[i]+str(calcInfo.atomCounts[i])
  i = min( max(0, poscar.rfind('/')+1), len(poscar) )
  calcInfo.jobID = poscar[i:]
  
def IsCalculationComplete(outcar, calcInfo = None):
  complete = False
  try:
    with open(outcar, 'r') as f:
      while True:
        l = f.readline()
        if not l:
          break
        if "Elapsed time" in l:
          complete = True
          if calcInfo:
            calcInfo.elapsedTime = float(l.split()[3])
    f.close()
  except:
    return False
  return complete

def GetOutcarData(outcar, calcInfo):
  calcInfo.lastEnergy = 0.
  calcInfo.complete = False
  try:
    f = open(outcar, 'r')
  except:
    try:
      f = open(outcar+'.out', 'r')
    except:
      print("Can't open",outcar,"(.out) for reading...")
      return
  while True:
    l = f.readline()
    if not l:
      break
    if "energy  without" in l:
      calcInfo.lastEnergy = float(l.split()[3])
    if "direct lattice vectors" in l:
      #calcInfo.A = vec3(f.readline().split()[:3])
      #calcInfo.B = vec3(f.readline().split()[:3])
      #calcInfo.C = vec3(f.readline().split()[:3])
      a = f.readline()
      calcInfo.A = vec3([a[:16],a[16:29],a[29:42]])
      a = f.readline()
      calcInfo.B = vec3([a[:16],a[16:29],a[29:42]])
      a = f.readline()
      calcInfo.C = vec3([a[:16],a[16:29],a[29:42]])
  f.close()
  calcInfo.complete = IsCalculationComplete(outcar, calcInfo)
  
def ReadKPoints(f, calc):
  l = f.readline() # Blank line
  l = f.readline().split() # x y z weight
  calc.kpoints.append(l)
  if calc.polarized:
    calc.bands.append([[],[]])
  else:
    calc.bands.append([])
  for i in range(calc.nBands):
    if not calc.polarized:
      l = f.readline().split()
      if not len(l) == 2:
        print("Incomplete Eignevalues file?",calc.dir,calc.title)
        return False
      calc.bands[-1].append(float(l[1]))
    else:
      l = f.readline()
      if not len(l) == 3:
        print("Incomplete Eignevalues file?",calc.dir,calc.title)
        return False
      calc.bands[-1][0].append(float(l[1]))
      calc.bands[-1][1].append(float(l[2]))
  return True
  
def GetBandData(rootname, calc):
  calc.bandStruct = False
  calc.bandComplete = False
  calc.nElectrons = 0
  calc.kpoints = []
  calc.nKpoints = 0
  calc.bands = []
  calc.nBands = 0
  calc.polarized = False
  calc.iVBE = 0
  calc.iCBE = 0
  try:
    f = open(rootname+"KPOINTS", 'r')
  except:
    print("Can't open",rootname+"KPOINTS for reading.")
    return
  l = f.readline()
  if "explicit k-points" in l.lower():
    calc.bandStruct = True
  f.close()
  if not calc.bandStruct:
    return
  try:
    f = open(rootname+"EIGENVAL", 'r')
  except:
    print("Can't open",rootname+"EIGENVAL for reading.")
    return
  l = f.readline().split() # 4 numbers... ? ? ? ISPIN
  if not len(l) == 4:
    print("Incomplete Eignevalues file?",rootname+"EIGENVAL")
    return
  calc.polarized = (2 == int(l[3]))
  f.readline() # 5 numbers... ???
  f.readline() # float
  f.readline() # CAR
  f.readline() # title
  l = f.readline().split() # 3 numbers
  if not len(l) == 3:
    print("Incomplete Eignevalues file?",rootname+"EIGENVAL")
    return
  calc.nElectrons = int(l[0])
  calc.nKpoints = int(l[1])
  calc.nBands = int(l[2])
  for i in range(calc.nKpoints):
    if not ReadKPoints(f, calc):
      return
  f.close()
  # Look through bands
  odd = (calc.nElectrons % 2 == 1)
  if not odd:
    calc.iVBE = int(calc.nElectrons/2)
    calc.iCBE = int(calc.iVBE+1)
    homo = calc.iVBE - 1; lumo = calc.iCBE - 1; # Index from 0
    hmax = -1e10; hmin = 1e10
    hmaxkp = hminkp = 0
    lmax = -1e10; lmin = 1e10
    lmaxkp = lminkp = 0
    gaps = []
    for i in range(calc.nKpoints):
      vb  = calc.bands[i][homo]
      cb  = calc.bands[i][lumo]
      gaps.append(cb - vb)
      if vb > hmax:
        hmax = vb; hmaxkp = i
      if vb < hmin:
        hmin = vb; hminkp = i
      if cb > lmax:
        lmax = cb; lmaxkp = i
      if cb < lmin:
        lmin = cb; lminkp = i
    calc.directGap = (hmaxkp == lminkp)
    calc.bandGap = lmin - hmax
    calc.directBandGap = calc.bandGap
    if not calc.directGap:
      calc.directBandGap = min(gaps)
    calc.bandComplete = True
  else: # Odd number of electrons
    return
    if not ef.polarized: # Copy above with elevated band edge indexes
      print("Calculation is closed shell with an odd number of electrons?")
      hAlpha = int((ef.valence-1)/2)+1
      lAlpha = hAlpha + 1
      hBeta = lAlpha
      lBeta = hBeta + 1
      pass
    else: # Calculation done right, get Alpha and Beta gaps
      hAlpha -=1; lAlpha -= 1; hBeta -= 1; lBeta -= 1; #Index
      hamax = -1e10; hamin = 1e10; hamaxkp = haminkp = 0
      hbmax = -1e10; hbmin = 1e10; hbmaxkp = hbminkp = 0
      lamax = -1e10; lamin = 1e10; lamaxkp = laminkp = 0
      lbmax = -1e10; lbmin = 1e10; lbmaxkp = lbminkp = 0
      for i in range(len(ef.kpoints)):
        if not ef.polarized:
          vb = ef.kpoints[i].bands[homo]
          cb = ef.kpoints[i].bands[lumo]
        else:
          vb = ef.kpoints[i].bands[0][homo]
          cb = ef.kpoints[i].bands[0][lumo]
        if odd and ef.polarized:
          vba = ef.kpoints[i].bands[0][hAlpha]
          vbb = ef.kpoints[i].bands[1][hBeta]
          cba = ef.kpoints[i].bands[0][lAlpha]
          cbb = ef.kpoints[i].bands[1][lBeta]
          #Alpha
          if vba > hamax:
            hamax = vba; hamaxkp = i
          if vba < hamin:
            hamin = vba; haminkp = i
          if cba > lamax:
            lamax = cba; lamaxkp = i
          if cba < lamin:
            lamin = cba; laminkp = i
          # Beta
          if vbb > hbmax:
            hbmax = vbb; hbmaxkp = i
          if vbb < hbmin:
            hbmin = vbb; hbminkp = i
          if cbb > lbmax:
            lbmax = cbb; lbmaxkp = i
          if cbb < lbmin:
            lbmin = cbb; lbminkp = i
        # End of odd else
        if vb > hmax:
          hmax = vb; hmaxkp = i
        if vb < hmin:
          hmin = vb; hminkp = i
        if cb > lmax:
          lmax = cb; lmaxkp = i
        if cb < lmin:
          lmin = cb; lminkp = i
      if odd and ef.polarized:
        if hamaxkp == laminkp:
          print("Direct band gap (alpha) at", ef.kpoints[hamaxkp].Coord())
        else:
          print("Indirect band gap (alpha)")
        print("Gap (alpha):",'{:5.3f}'.format(lamin-hamax),"eV")
        if hbmaxkp == lbminkp:
          print("Direct band gap (beta) at", ef.kpoints[hbmaxkp].Coord())
        else:
          print("Indirect band gap (beta)")
        print("Gap  (beta):",'{:5.3f}'.format(lbmin-hbmax),"eV")
        print("Forbidden alpha -> beta:",'{:5.3f}'.format(lbmin-hamax),"eV")
        print("Forbidden beta -> alpha:",'{:5.3f}'.format(lamin-hbmax),"eV")
        print("Alpha VBE, CBE, Gap, VBW, CBW: \
        {0:5.3f}, {1:5.3f}, {2:5.3f}, {3:5.3f}, {4:5.3f} eV".format(hamax,lamin,
        lamin-hamax,hamax-hamin,lamax-lamin))
        print("Beta  VBE, CBE, Gap, VBW, CBW: \
        {0:5.3f}, {1:5.3f}, {2:5.3f}, {3:5.3f}, {4:5.3f} eV".format(hbmax,lbmin,
        lbmin-hbmax,hbmax-hbmin,lbmax-lbmin))
      else:
        if hmaxkp == lminkp:
          print("Direct band gap at", ef.kpoints[hmaxkp].Coord())
        else:
          print("Indirect band gap")
        print("Gap:",'{:5.3f}'.format(lmin-hmax),"eV")
        print("VBE, CBE, Gap, VBW, CBW: \
        {0:5.3f}, {1:5.3f}, {2:5.3f}, {3:5.3f}, {4:5.3f} eV".format(hmax,lmin,
        lmin-hmax,hmax-hmin,lmax-lmin))

def ParseTime(elapsed):
  sec = elapsed % 60.
  minu = elapsed // 60 % 60
  hour = elapsed // 3600 % 24
  days = elapsed // (3600 * 24) % 30
  months = elapsed // (3600 * 24 * 30)
  return "{0:>3} months {1:>2} days {2:>2} hours {3:>2} minutes {4:>4.1f} seconds".format(months,days,hour,minu,sec)
  
# Write calculation data to html        
def WriteHTML(calcName, calcInfo):
  with open(calcName+".html", 'w') as f:
    f.write("<html>\n<body>\n")
    f.write('<h1>'+calcName+'</h1>\n')
    f.write('<p>\n')
    f.write("Last Known Whereabouts: "+calcInfo.dir+'</br>\n')
    f.write("Title: "+calcInfo.title+'</br>\n')
    f.write("ID (GoVASP): "+str(calcInfo.jobID)+'</br>\n')
    f.write("Atoms: "+str(calcInfo.nAtoms)+" Atom sets:"+str(calcInfo.atomCounts)+'</br>\n')
    f.write("System: "+calcInfo.Emp+'</br>\n')
    f.write("A: {0}</br>\nB: {1}</br>\nC: {2}</br>\n".format(calcInfo.A,calcInfo.B,calcInfo.C))
    f.write("k-points: "+str(calcInfo.nKpoints)+'</br>\n')
    if not calcInfo.complete:
      f.write("</p><h2>CALCULATION IS INCOMPLETE</h2><p>\n")
    else:
      f.write("Calculation complete in: "+ParseTime(calcInfo.elapsedTime)+'</br>\n')
      f.write("Energy without entropy: {0:> 15.8e} eV</br>\n".format(calcInfo.lastEnergy))
      if calcInfo.polarized:
        f.write("ISPIN = 2 (Spin Polarized)</br>\n")
      else:
        f.write("ISPIN = 1 (Non-Polarized)</br>\n")
      f.write(str(calcInfo.nElectrons)+" Electrons in "+str(calcInfo.nBands)+" Bands</br>\n")
      if calcInfo.bandStruct:
        if not calcInfo.bandComplete:
          f.write("</p><h2>INCOMPLETE BAND STRUCTURE CALCULATION</h2><p>\n")
        else:
          f.write("Completed band structure calculation</br>\n")
          f.write("VB, CB: "+str(calcInfo.iVBE)+", "+str(calcInfo.iCBE)+"</br>\n")
          if calcInfo.directGap:
            f.write("Direct Band Gap: {0:> 6.3f}</br>\n".format(calcInfo.bandGap))
          else:
            f.write("Indirect Band Gap: {0:> 6.3f}</br>\n".format(calcInfo.bandGap))
            f.write("<tab>Direct Band Gap: {0:> 6.3f}</br>\n".format(calcInfo.directBandGap))
    
    f.write('</p>\n</body>\n</html>\n')
  f.close()
  
# Read in old index and output new index file
def AddToIndex(htmlFile, linkName, htmlIndex, complete):
  # Read links from old index file
  linx = []; flinx = []; incomp = False
  try:
    with open(htmlIndex, 'r') as f:
      if "<html" not in f.readline():
        print("Bad header in index file?")
        f.close()
        return
      else:
        while True:
          l = f.readline()
          if not l:
            print("Error: Premature end of index file")
            break
          if "</html>" in l:
            break # Good EOF
          if "<h4> Incomplete" in l:
            incomp = True
          elif "<a href=" in l:
            if incomp:
              flinx.append(l)
            else:
              linx.append(l)
    f.close()
  except:
    pass
  i = 0
  attempts = 1
  ohtml = htmlFile
  if complete:
    while i < len(linx):
      if htmlFile in linx[i]:
        htmlFile = ohtml[:-5]+str(attempts)+ohtml[-5:]
        attempts += 1
        i = 0
      i += 1
  else:
    while i < len(flinx):
      if htmlFile in flinx[i]:
        htmlFile = ohtml[:-5]+str(attempts)+ohtml[-5:]
        attempts += 1
        i = 0
      i += 1
  if complete:
    linx.append('<a href="'+htmlFile+'" target="frame-1">'+linkName+' </a></br>\n')
  else:
    flinx.append('<a href="'+htmlFile+'" target="frame-1">'+linkName+' </a></br>\n')
  # Rewrite entire index file...
  with open(htmlIndex, 'w') as f:
    f.write('<html>\n<body>\n')
    for i in sorted(linx):
      f.write(i)
    f.write("<h4> Incomplete Calculations </h4>\n")
    for i in sorted(flinx):
      f.write(i)
    f.write('</body>\n</html>\n')
  f.close()
  return htmlFile[:-5]

# Read in old index and output new alphabetized index file
def RewriteIndex(htmlIndex = "DirIndex.html"):
  # Read links from old index file
  linx = []; flinx = []; incomp = False
  try:
    with open(htmlIndex, 'r') as f:
      if "<html" not in f.readline():
        print("Bad header in index file?")
        f.close()
        return
      else:
        while True:
          l = f.readline()
          if not l:
            print("Error: Premature end of index file")
            break
          if "</html>" in l:
            break # Good EOF
          if "<h4> Incomplete" in l:
            incomp = True
          elif "<a href=" in l:
            if incomp:
              flinx.append(l)
            else:
              linx.append(l)
    f.close()
  except:
    print("Can't rewrite non-existant index file...")
    return
  # Sort
  linx.sort(); flinx.sort()
  # Rewrite entire index file...
  with open(htmlIndex, 'w') as f:
    f.write('<html>\n<body>\n')
    # Divide right side
    f.write('<div style="border: 1px solid; float: left; overflow-y:scroll;overflow-x:scroll width: 20%; height:100%;" scrolling="yes" id="index">\n')
    for i in linx:
      f.write(i)
    f.write("<h4> Incomplete Calculations </h4>\n")
    for i in flinx:
      f.write(i)
    f.write('</div>\n')
    # Divide left side
    f.write('<div style="border: 1px solid; float: right; width: 75%; height:100%;position:relative;top:0px;left:0px;" id="frame">\n')
    # Frame buffer
    f.write('<iframe name="frame-1" width="100%" height="100%" style="border:none"></iframe>\n')
    f.write('</div>\n')
    f.write('</body>\n</html>\n')
  f.close()

"""Output calculation data to html?
Have a list down the left side and
a frame view to the right with the info."""
def CalcToHTML(calcName, calcInfo, htmlindex = "DirIndex.html"):
  # Add new file to index file
  calcName = AddToIndex(calcName+".html", calcName, htmlindex, calcInfo.complete)
#  with open(htmlindex, 'r+') as f:
#    if "<html" not in f.readline():
#      break
#    pos = f.tell()
#    while True:
#      l = f.readline()
#      if not l:
#        print("Error: Premature end of index file")
#        break
#      if "</html>" in l:
#        break # Good EOF
#      if "Last File" in l:
#        # Found last entry
#        f.seek(pos) # Insert new link
#        f.write('<a href="'+calcName+'.html">'+calcName+' </a>\n')
#        break
#      pos = f.tell()
#  f.close()
  # Create new html file
  WriteHTML(calcName, calcInfo)

'''Remove leading path with last /
/usr/local/bin -> bin
'''
def StripPath(filename):
  i = filename.rfind('/') #max(filename.rfind('/'),filename.rfind('\\'))
  i = max(0,i+1)
  return filename[i:]

def ReadDirectory(args):
  d = "."
  readpot = False
  readEnergy = False
  readBand = False
  printFormat = False
  i = 1
  jobs = []
  while i < len(args):
    a = args[i]
    if a[0] == '-':
      if a[1] == 'P':
        readpot = True
        print("Reading POTCAR and POSCAR to index empirical formula.")
      elif a[1] == 'E':
        readEnergy = True
        print('Reading OUTCAR for final "energy  without entropy" entry.')
      elif a[1] == 'J':
        print("Only reading from list of jobs:")
        jobs = args[i+1].split(',')
        print(jobs)
        i += 1
      elif a[1:5] == 'band':
        readBand = True
        print('Reading EIGENVAL for bandgaps, also check KPOINTS for explicit k-points.')
      elif a[1] == 'A':
        printFormat = True
        readBand = True
        readpot = True
      elif a[1] == 'h' or a[1:3] == '-h':
        print('Reads through a directory to read through VASP calculations. (ignores links)')
        print('Results are printed to IndexData.txt or JobsData.txt.')
        print('Options: -P, -E, -J, -band, -h')
        print('\t-P Read POTCAR to determine atoms in calculation and print Formula')
        print('\t-E Read OUTCAR to determine final energy or if calculation ran at all')
        print('\t-J Only read listed directories separated by a comma (e.g. 0,24,546,SiC)')
        print('\t   This will print out to JobsData.txt')
        print('\t-band Read EIGENVAL and KPOINTS files to give the band gap if complete')
        return
    else:
      d = a
    i+=1
  # Start doing stuff
  calc = CalcInfo()
  i = -1
  for root, dirs, files in os.walk(d):
    if len(jobs) == 0:
      for name in files:
        if "POSCAR" in name:
          calc.dir = root
          if readEnergy:
            GetOutcarData(root+'/'+name[:-6]+'OUTCAR', calc)
          if readBand:
            GetBandData(root+'/'+name[:-6], calc)
          if not readpot:
            #print(root,"/",name,": ",GetAtomCount(root+'/'+name), sep="")
            GetPoscarData(root+'/'+name, calc)
            with open("Atoms{0}.txt".format(calc.nAtoms), 'a') as f:
              f.write(root+"/"+name+'\n')
          else:
            GetPotcarData(root+'/'+name, calc)
          with open("IndexData.txt", 'a') as f:
            out = ""
            if readpot:
              # Title contains \n
              out = root+'\t'+calc.Emp+'\t'+str(calc.nAtoms)
              if readEnergy:
                if calc.complete:
                  out += "\tE= "+str(calc.lastEnergy)
                else:
                  out += "\tINCOMPLETE"
              if readBand:
                if not calc.bandStruct:
                  out += "\tNoBand\tN/A"
                elif not calc.bandComplete:
                  out += "\tIncomplete\tN/A"
                elif calc.directGap:
                  out += "\tDirect Gap\t{0: <5.2f}".format(calc.directBandGap)
                else:
                  out += "\tIndirect Gap\t{0: <5.2f}".format(calc.bandGap)
              out += '\t'+calc.title
            else:
              out = StripPath(root)
              if readEnergy:
                if calc.complete:
                  out += "\tE= "+str(calc.lastEnergy)
                else:
                  out += "\tINCOMPLETE"
                if readBand:
                  if not calc.bandStruct:
                    out += "\tNoBand\tN/A"
                  elif not calc.bandComplete:
                    out += "\tIncomplete\tN/A"
                  elif calc.directGap:
                    out += "\tDirect Gap\t{0: <5.2f}".format(calc.directBandGap)
                  else:
                    out += "\tIndirect Gap\t{0: <5.2f}".format(calc.bandGap)
              out += '\t'+calc.title+'\t'+str(calc.nAtoms)+'\n'
            f.write(out)
            if printFormat:
              with open("JobsData-A.txt", 'a') as f:
                bPrint = True
                count = [0,0,0,0]
                for i in range(len(calc.atoms)):
                  if calc.atoms[i] == "Al":
                    count[0] = calc.atomCounts[i]
                  elif calc.atoms[i] == "Ga":
                    count[1] = calc.atomCounts[i]
                  elif calc.atoms[i] == "In":
                    count[2] = calc.atomCounts[i]
                  elif calc.atoms[i] == "N":
                    count[3] = calc.atomCounts[i]
                    continue
                  else:
                    bPrint = False
                if bPrint:
                  out = StripPath(root) + '\t' + calc.Emp
                  out += '\t{0}\t{1}\t{2}\t{3}'.format(count[0],count[1],count[2],count[3])
                  if calc.bandComplete:
                    out += '\t{0: <5.2f}'.format(calc.bandGap)
                    if calc.directGap:
                      out += '\tDirect'
                    else:
                      out += '\tIndirect'
                  out += '\n'
                  f.write(out)
            if calc.complete:
              CalcToHTML(calc.Emp+" "+str(i),calc)
            else:
              CalcToHTML(calc.Emp+" "+str(i)+" Incomplete",calc)
            i += 1
    else: # Read only from listed jobs
      calc.dir = root
      if StripPath(root) in jobs:
        for name in files:
          if "POSCAR" in name:
            i += 1
            if readEnergy:
              GetOutcarData(root+'/'+name[:-6]+'OUTCAR', calc)
            if readBand:
              GetBandData(root+'/'+name[:-6], calc)
            if not readpot:
              #print(root,"/",name,": ",GetAtomCount(root+'/'+name), sep="")
              GetPoscarData(root+'/'+name, calc)
              with open("Atoms{0}.txt".format(calc.nAtoms), 'a') as f:
                f.write(root+"/"+name+'\n')
            else:
              GetPotcarData(root+'/'+name, calc)
            with open("JobsData.txt", 'a') as f:
              out = ""
              if readpot:
                # Title contains \n
                out += root+'\t'+calc.Emp+'\t'+str(calc.nAtoms)
                if readEnergy:
                  if calc.complete:
                    out += "\tE= "+str(calc.lastEnergy)
                  else:
                    out += "\tINCOMPLETE"
                if readBand:
                  if not calc.bandStruct:
                    out += "\tNoBand\tN/A"
                  elif not calc.bandComplete:
                    out += "\tIncomplete\tN/A"
                  elif calc.directGap:
                    out += "\tDirect Gap\t{0: <5.2f}".format(calc.directBandGap)
                  else:
                    out += "\tIndirect Gap\t{0: <5.2f}".format(calc.bandGap)
                out += '\t'+calc.title
              else:
                out = StripPath(root)
                if readEnergy:
                  if calc.complete:
                    out += "\tE= "+str(calc.lastEnergy)
                  else:
                    out += "\tINCOMPLETE"
                if readBand:
                  if not calc.bandStruct:
                    out += "\tNoBand\tN/A"
                  elif not calc.bandComplete:
                    out += "\tIncomplete\tN/A"
                  elif calc.directGap:
                    out += "\tDirect Gap\t{0: <5.2f}".format(calc.directBandGap)
                  else:
                    out += "\tIndirect Gap\t{0: <5.2f}".format(calc.bandGap)
                out += '\t'+calc.title+'\t'+str(calc.nAtoms)+'\n'
              f.write(out)
              if printFormat:
                with open("JobsData-A.txt", 'a') as f:
                  bPrint = True
                  count = [0,0,0,0]
                  for i in range(len(calc.atoms)):
                    if calc.atoms[i] == "Al":
                      count[0] = calc.atomCounts[i]
                    elif calc.atoms[i] == "Ga":
                      count[1] = calc.atomCounts[i]
                    elif calc.atoms[i] == "In":
                      count[2] = calc.atomCounts[i]
                    elif calc.atoms[i] == "N":
                      count[3] = calc.atomCounts[i]
                      continue
                    else:
                      bPrint = False
                  if bPrint:
                    out = StripPath(root) + '\t' + calc.Emp
                    out += '\t{0}\t{1}\t{2}\t{3}'.format(count[0],count[1],count[2],count[3])
                    if calc.bandComplete:
                      out += '\t{0: <5.2f}'.format(calc.bandGap)
                      if calc.directGap:
                        out += '\tDirect'
                      else:
                        out += '\tIndirect'
                    out += '\n'
                    f.write(out)
            if calc.complete:
              CalcToHTML(calc.Emp+" "+str(i),calc)
            else:
              CalcToHTML(calc.Emp+" "+str(i)+" Incomplete",calc)
              i += 1
  RewriteIndex()
  print("Good Day")
  
if __name__ == "__main__":
  import sys
  ReadDirectory(sys.argv)

