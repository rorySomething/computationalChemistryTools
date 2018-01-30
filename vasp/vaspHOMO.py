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
    
class OUTCARFILE:
  def __init__(self):
    self.valence = 0
    self.kpts = 0
    self.bands = 0
    self.polarized = False
    self.kpoints = []
    self.eFermi = 0
    self.XC = 0
    self.alphaBeta = 0

ef = EIGENFILE()

efermi = XC = alpha_beta = 0

def BandAverageOUTCAR(ot):
  print("TODO: Check that the same band is being averaged when using occupancy for HOMO as with OUTCAR reading")
  homo = 0
  lumo = 0
  odd = ot.valence % 2 == 1
  homo = ot.valence//2 - 1 # -1 (index from zero)
  lumo = homo+1
  print("HOMO", homo+1)
  if not ot.polarized:
    print("Not ready to handle non-polarized")
    return
  else:
    homoSet = False
    hA = 0; hB = 0; lA = 0; lB = 0
    totalhA = 0; MA = 0; mA = 0
    totallA = 0; lMA = 0; lmA = 0
    for i in range(ot.kpts):
      kpt = ot.kpoints[i]
      eA = 0
      for j in range(ot.bands):
        if j <= homo: #not homoSet and kpt.bands[0][j][1] > 0.5: # bands[alpha, beta][bandIndex][energy, occupancy]
          #                                         # Occupied (>0.5)
          hA = j ; eA = kpt.bands[0][j][0]
        else:
          # Found Lumo
          #homoSet = True
          totalhA += eA
          leA = kpt.bands[0][j][0]
          totallA += leA
          if i == 0:
            MA = mA = eA
            lMA = lmA = leA
          if eA > MA:
            MA = eA
          if eA < mA:
            mA = eA
          if leA > lMA:
            lMA = leA
          if leA < lmA:
            lmA = leA
          break
    print("Alpha HOMO {0}".format(hA+1))
    print("\tHOMO band energy mean: {0:.3f} eV".format(totalhA/ot.kpts+ot.alphaBeta))
    print("\t\tHOMO range {0:.3f} - {1:.3f} eV".format(mA+ot.alphaBeta, MA+ot.alphaBeta))
    print("\tLUMO band energy mean: {0:.3f} eV".format(totallA/ot.kpts+ot.alphaBeta))
    print("\t\tLUMO range {0:.3f} - {1:.3f} eV".format(lmA+ot.alphaBeta, lMA+ot.alphaBeta))
    homoSet = False
    hA = 0; hB = 0; lA = 0; lB = 0
    totalhA = 0; MA = 0; mA = 0
    totallA = 0; lMA = 0; lmA = 0
    for i in range(ot.kpts):
      kpt = ot.kpoints[i]
      eA = 0
      for j in range(ot.bands):
        if j <= homo: # not homoSet and kpt.bands[1][j][1] > 0.5:
          #                # Occupied (>0.5)
          hA = j ; eA = kpt.bands[1][j][0]
        else:
          # Found Lumo
          homoSet = True
          totalhA += eA
          leA = kpt.bands[1][j][0]
          totallA += leA
          if i == 0:
            MA = mA = eA
            lMA = lmA = leA
          if eA > MA:
            MA = eA
          if eA < mA:
            mA = eA
          if leA > lMA:
            lMA = leA
          if leA < lmA:
            lmA = leA
          break
    print("Beta HOMO {0}".format(hA+1))
    print("\tHOMO band energy mean: {0:.3f} eV".format(totalhA/ot.kpts+ot.alphaBeta))
    print("\t\tHOMO range {0:.3f} - {1:.3f} eV".format(mA+ot.alphaBeta, MA+ot.alphaBeta))
    print("\tLUMO band energy mean: {0:.3f} eV".format(totallA/ot.kpts+ot.alphaBeta))
    print("\t\tLUMO range {0:.3f} - {1:.3f} eV".format(lmA+ot.alphaBeta, lMA+ot.alphaBeta))
  return

def ReadOUTCARKPoints(f, ot):
  if ot.polarized:
    # Read after E-fermi :
    # Blank Line
    # Blank Line
    # spin component 1
    # *Then repeated below
    # Blank Line
    # k-point     1 : x y z
    #  band No.   band energies   occupation
    # ...
    #  band No.   band energies   occupation
    # Blank
    # spin component 2
    # *Then repeat below
    # Blank Line
    # k-point     1 : x y z
    #  band No.   band energies   occupation
    # ...
    #  band No.   band energies   occupation
    f.readline() # Blank line
    f.readline() # Blank line
    a = f.readline()
    if not "spin component" in a:
      print("Error reading k-point data")
      return False
    for i in range(ot.kpts):
      f.readline() # Blank line
      a = f.readline().split()
      k = KPOINT()
      k.x = float(a[3])
      k.y = float(a[4])
      k.z = float(a[5])
      k.bands = [[],[]]
      f.readline() # band No. energy occupancy header
      for j in range(ot.bands):
        a = f.readline().split()
        # bandIndex  Energy  Occ.
        # Append (Energy, occupancy)
        k.bands[0].append((float(a[1]), float(a[2])))
      ot.kpoints.append(k)
    # Read beta spin
    f.readline() # Blank
    a = f.readline()
    if not "spin component" in a:
      print("Error reading k-point data")
      return False
    for i in range(ot.kpts):
      f.readline() # Blank line
      a = f.readline().split()
#      k = KPOINT()
#      k.x = float(a[3])
#      k.y = float(a[4])
#      k.z = float(a[5])
#      k.bands = [[],[]]
      f.readline() # band No. energy occupancy header
      for j in range(ot.bands):
        a = f.readline().split()
        # bandIndex  Energy  Occ.
        # Append (Energy, occupancy)
        ot.kpoints[i].bands[1].append((float(a[1]), float(a[2])))
  else:
    print("Not ready to read non-polarized data...")
    return False
  return True
def GetAlphaBeta(directory):
  apb = 0.
  with open(directory+"OUTCAR", 'r') as f:
    for l in f:
      if "E-fermi :" in l:
        a = l.split()
        aB = a[-1]
        if aB[0] == ":":
          aB = aB[1:]
        apb = float(aB)
        break
  f.close()
  return apb
def ReadOUTCAR(directory):
  try:
    f = open(directory+"OUTCAR", 'r')
  except:
    print("Can't open OUTCAR file")
    return
  ot = OUTCARFILE()
  slic = len(" E-fermi :")
  efermi = False
  XC = False
  alpha_beta = False
  for l in f:
    if "NKPTS" in l:
      l = l.split()
      ot.kpts = int(l[3])
      ot.bands = int(l[-1]) # NBANDS
    elif "ISPIN" in l:
      l = l.split()
      ot.polarized = (int(l[2]) == 2)
    elif "E-fermi :" in l:
      ot.kpoints.clear()
      a = l.split()
      ot.eFermi = float(a[2])
      ot.XC = float(a[4])
      aB = a[-1]
      if aB[0] == ":":
        aB = aB[1:]
      ot.alphaBeta = float(aB)
      ReadOUTCARKPoints(f, ot)
    elif "NELECT" in l:
      ot.valence = int(float(l.split()[2]))
  f.close()
  print("E-fermi: {0:6.3f} Alpha+Beta: {1:7.3f}".format(ot.eFermi+ot.alphaBeta, ot.alphaBeta))
  BandAverageOUTCAR(ot)

def BandAverage(ef, apb = 0.):
  odd = ef.valence % 2 == 1
  homo = ef.valence//2 - 1 # -1 (index from zero)
  lumo = homo+1
  if not ef.polarized:
    bandAverages = []
    bandMax = []
    bandMin = []
    for i in range(ef.bands):
      total = 0
      bandMax.append(ef.kpoints[0].bands[i] + apb)
      bandMin.append(ef.kpoints[0].bands[i] + apb)
      for j in range(len(ef.kpoints)):
        e = ef.kpoints[j].bands[i] + apb
        total += e
        if e > bandMax[-1]:
          bandMax[-1] = e
        if e < bandMin[-1]:
          bandMin[-1] = e
      bandAverages.append(total/ef.kpts)
    if odd:
      print("Alpha HOMO {0}".format(homo+1))
      print("\tHOMO band energy mean: {0:.3f} eV".format(bandAverages[homo+1]))
      print("\t\tHOMO range {0:.3f} - {1:.3f} eV".format(bandMin[homo+1], bandMax[homo+1]))
      print("\tLUMO band energy mean: {0:.3f} eV".format(bandAverages[lumo+1]))
      print("\t\tLUMO range {0:.3f} - {1:.3f} eV".format(bandMin[lumo+1], bandMax[lumo+1]))
      print("\tLUMO+1 range {0:.3f} - {1:.3f} eV".format(bandMin[lumo+3], bandMax[lumo+3]))
      print("Beta")
      print("\tHOMO band energy mean: {0:.3f} eV".format(bandAverages[homo]))
      print("\t\tHOMO range {0:.3f} - {1:.3f} eV".format(bandMin[homo], bandMax[homo]))
      print("\tLUMO band energy mean: {0:.3f} eV".format(bandAverages[lumo]))
      print("\t\tLUMO range {0:.3f} - {1:.3f} eV".format(bandMin[lumo], bandMax[lumo]))
      print("\t\tLUMO range {0:.3f} - {1:.3f} eV".format(bandMin[lumo+2], bandMax[lumo+2]))
    #-----------Even Occupancy--------------------------------------------------
    else:
      print("HOMO {1} band energy mean: {0:.3f} eV".format(bandAverages[homo], homo+1))
      print("\tHOMO range {0:.3f} - {1:.3f} eV".format(bandMin[homo], bandMax[homo]))
      print("LUMO band energy mean: {0:.3f} eV".format(bandAverages[lumo]))
      print("\tLUMO range {0:.3f} - {1:.3f} eV".format(bandMin[lumo], bandMax[lumo]))
      print("\tLUMO range {0:.3f} - {1:.3f} eV".format(bandMin[lumo+1], bandMax[lumo+1]))
  #-----------Spin Polarized----------------------------------------------------
  else:
    bandAveragesA = [] ; bandAveragesB = []
    bandMaxA = [] ; bandMaxB = []
    bandMinA = [] ; bandMinB = []
    for i in range(ef.bands):
      totalA = 0 ; totalB = 0
      bandMaxA.append(ef.kpoints[0].bands[0][i] + apb)
      bandMinA.append(ef.kpoints[0].bands[0][i] + apb)
      bandMaxB.append(ef.kpoints[0].bands[1][i] + apb)
      bandMinB.append(ef.kpoints[0].bands[1][i] + apb)
      for j in range(len(ef.kpoints)):
        eA = ef.kpoints[j].bands[0][i] + apb
        eB = ef.kpoints[j].bands[1][i] + apb
        totalA += eA ; totalB += eB
        if eA > bandMaxA[-1]:
          bandMaxA[-1] = eA
        if eA < bandMinA[-1]:
          bandMinA[-1] = eA
        if eB > bandMaxB[-1]:
          bandMaxB[-1] = eB
        if eB < bandMinB[-1]:
          bandMinB[-1] = eB
      bandAveragesA.append(totalA/ef.kpts)
      bandAveragesB.append(totalB/ef.kpts)
    if odd:
      print("Doublet state")
      print("HFOCO (Highest fully occupied crystal orbital)")
      print("Alpha States HFOCO {0}".format(homo+1))
      # State  VBMin     VBE      CBE     CBMax     VBA    CBA
      print("{0:^10}{1:^10}{2:^10}{3:^10}{4:^10}{5:^10}{6:^10}".format("","VBMin","VBE","CBE","CBMax","VBA","CBA"))
      print("{0:^10}{1:>10.3f}{2:>10.3f}{3:>10.3f}{4:>10.3f}{5:>10.3f}{6:>10.3f}".format("HFOCO",bandMinA[homo], bandMaxA[homo],bandMinA[lumo], bandMaxA[lumo], bandAveragesA[homo],bandAveragesA[lumo]))
      print("{0:^10}{1:>10.3f}{2:>10.3f}{3:>10.3f}{4:>10.3f}{5:>10.3f}{6:>10.3f}".format(" Beta",bandMinB[homo], bandMaxB[homo],bandMinB[lumo], bandMaxB[lumo], bandAveragesB[homo],bandAveragesB[lumo]))
      print("{0:^10}{1:>10.3f}{2:>10.3f}{3:>10.3f}{4:>10.3f}{5:>10.3f}{6:>10.3f}".format("HOCO",bandMinA[lumo], bandMaxA[lumo],bandMinA[lumo+1], bandMaxA[lumo+1], bandAveragesA[lumo],bandAveragesA[lumo+1]))
      print("{0:^10}{1:>10.3f}{2:>10.3f}{3:>10.3f}{4:>10.3f}{5:>10.3f}{6:>10.3f}".format("Beta",bandMinB[lumo], bandMaxB[lumo],bandMinB[lumo+1], bandMaxB[lumo+1], bandAveragesB[lumo],bandAveragesB[lumo+1]))
      print("{0:^10}{1:>10.3f}{2:>10.3f}{3:>10.3f}{4:>10.3f}{5:>10.3f}{6:>10.3f}".format("LUMO+1",bandMinA[lumo+1], bandMaxA[lumo+1],bandMinA[lumo+2], bandMaxA[lumo+2], bandAveragesA[lumo+1],bandAveragesA[lumo+2]))
    #-----------Even Occupancy--------------------------------------------------
    else:
      print("Alpha HOMO {0}".format(homo+1))
      print("\tHOMO band energy mean: {0:.3f} eV".format(bandAveragesA[homo]))
      print("\t\tHOMO range {0:.3f} - {1:.3f} eV".format(bandMinA[homo], bandMaxA[homo]))
      print("\tLUMO band energy mean: {0:.3f} eV".format(bandAveragesA[lumo]))
      print("\t\tLUMO range {0:.3f} - {1:.3f} eV".format(bandMinA[lumo], bandMaxA[lumo]))
      print("\tLUMO+1 range {0:.3f} - {1:.3f} eV".format(bandMinA[lumo+1], bandMaxA[lumo+1]))
      print("Beta")
      print("\tHOMO band energy mean: {0:.3f} eV".format(bandAveragesB[homo]))
      print("\t\tHOMO range {0:.3f} - {1:.3f} eV".format(bandMinB[homo], bandMaxB[homo]))
      print("\tLUMO band energy mean: {0:.3f} eV".format(bandAveragesB[lumo]))
      print("\t\tLUMO range {0:.3f} - {1:.3f} eV".format(bandMinB[lumo], bandMaxB[lumo]))
      print("\tLUMO+1 range {0:.3f} - {1:.3f} eV".format(bandMinB[lumo+1], bandMaxB[lumo+1]))

def ReadKPoints(f, polarized, bands):
  # Blank Line
  # x y z weight
  # bandIndex energyUp energyDown
  # bandIndex energyUp energyDown
  # ...
  # bandIndex energyUp energyDown
  f.readline() # Blank line
  a = f.readline().split()
  k = KPOINT()
  k.bands.clear()
  k.x = float(a[0])
  k.y = float(a[1])
  k.z = float(a[2])
  k.w = float(a[3])
  if polarized:
    k.bands = [[],[]]
  for i in range(bands):
    if not polarized:
      k.bands.append(float(f.readline().split()[1]))
    else:
      a = f.readline().split()
      k.bands[0].append(float(a[1]))
      k.bands[1].append(float(a[2]))
  return k

def ReadEIGENVAL(directory, apb = 0.):
  try:
    f = open(directory+"EIGENVAL", 'r')
  except:
    print("Can't open EIGENVAL file")
    return
  ef = EIGENFILE()
  ef.l1 = f.readline().split()
  # 9  9  100   2
  # ?  ?  ?     ISPIN
  sp = ef.l1[3]
  ef.polarized = (2 == int(sp))
  ef.l2 = f.readline().split()
  # 0.1817149E+02  0.5468620E-09  0.5468620E-09  0.5468620E-09  0.5000000E-15
  # ? ? ? ? ?
  ef.a3 = float(f.readline())
  # 1.0000000000000000E-004
  # ?
  ef.CAR = f.readline()
  # CAR
  ef.title = f.readline()
  # Si9
  a = f.readline().split()
  # 36   1305     32
  # valence  kpoints   bands
  ef.valence = int(a[0])
  ef.kpts = int(a[1])
  ef.bands = int(a[2])
  print("Valence: {0:<5} KPTS: {1:<5} Bands: {2:<5}".format(ef.valence, ef.kpts, ef.bands))
  # Loop to read kpoints
  for i in range(ef.kpts):
    ef.kpoints.append(ReadKPoints(f, ef.polarized, ef.bands))
  f.close()
#  ef.Print() # To confirm good read
  BandAverage(ef, apb)
  return
  
  ###
  ###
  print("Passed Return???")
  ###
  ### OLD Below
  ###
  odd = ef.valence % 2 == 1
  homo = ef.valence//2
  lumo = homo+1
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
  import sys
  args = sys.argv
  directory = "./"
  if len(args) > 1:
    directory = args[1]
    if "/" != directory[-1]:
      directory += "/"
  apb = GetAlphaBeta(directory)
  print("Reading EIGENVAL")
  ReadEIGENVAL(directory, apb)
  print("Reading OUTCAR")
  ReadOUTCAR(directory)

