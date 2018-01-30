#!/bin/python3

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
  efermi = False
  XC = False
  alpha_beta = False
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
  return [efermi, XC, alpha_beta]

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
  x = ReadOUTCAR()
  print("Vacuum adjusted\n\tVBE, E-fermi, CBE: {0:5.3f}, {1:5.3f}, {2:5.3f}"
  .format(hmax + x[2], x[0] + x[2], lmin + x[2]))

if __name__ == "__main__":
  Read_Eigenval()

