#!/bin/python3

from math import sqrt
import math
import subprocess

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
class vec3:
  x = y = z = 0.
  def __init__(self):
    self.x = self.y = self.z = 0.
  def __init__(self, x, y, z):
    self.x = x; self.y = y; self.z = z
  def __init__(self, x):
    self.x = float(x[0]); self.y = float(x[1]); self.z = float(x[2])    
  def Distance(self, v):
    x = self.x - v.x; y = self.y - v.y; z = self.z - v.z
    return sqrt(x*x+y*y+z*z)
  def __sub__(self, a):
    return vec3([self.x-a.x, self.y-a.y, self.z-a.z])
  def __str__(self):
    return "{0:< 6.3f},{1:< 6.3f},{2:< 6.3f}".format(self.x,self.y,self.z)
  def __neg__(self):
    return vec3([-self.x, -self.y, -self.z])
  def Length(self):
    return sqrt(self.x*self.x+self.y*self.y+self.z*self.z)
  def Normalize(self):
    l = self.Length()
    if l == 0.:
      return vec3([0.,0.,0.])
    return vec3([self.x/l, self.y/l, self.z/l])
  def __add__(self,other):
    return vec3([self.x+other.x,self.y+other.y,self.z+other.z])
  def __mul__(self,other):
    return vec3([self.x*float(other),self.y*float(other),self.z*float(other)])
  def __truediv__(self,other):
    return vec3([self.x/float(other),self.y/float(other),self.z/float(other)])

def OUTCAR_ReadKPoint(f, points, bandData, nBands):
  f.readline() # Blank line
  a = f.readline().split()
  if len(a) != 6:
    return False
  points.append(vec3(a[3:])) # k-point x,y,z
  f.readline() # "  band No.  band energies     occupation "
  bandData.append( [points[-1], [], []] )
    # Append array with (k-point, [band energy], [occupancy])
  for i in range(nBands):
    a = f.readline().split()
    bandData[-1][1].append(float(a[1]))  # Append band energy
    bandData[-1][2].append(float(a[2]))  # Append occupancy
  return True
  
def printBandEdges(nElect, nBands, vacAlpha, vacBeta):
  if nBands < 1 or nElect < 1 or len(vacAlpha) < 1:
    return False
  polarized = len(vacBeta) > 0
  odd = (nElect % 2 == 1)
  homo = int(nElect/2)
  lumo = int(homo+1)
#  print("VB, CB: {0}, {1}".format(homo,lumo))
  if odd:
#    print("Odd number of electrons")
#    if not ef.polarized:
#      print("Calculation is closed shell with an odd number of electrons?")
#      pass
    hAlpha = int((nElect-1)/2)+1
    lAlpha = hAlpha + 1
    hBeta = lAlpha
    lBeta = hBeta + 1
    #print("Alpha VB, CB: {0}, {1}".format(hAlpha,lAlpha))
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
    print("HOMO/LUMO_a:",hAlpha,lAlpha)
    print("HOMO/LUMO_B:",hBeta,lBeta)
  kpoints = int(len(vacAlpha) / nBands)
  for i in range(kpoints):
    index = i * nBands
    if not polarized:
      vb = vacAlpha[index + homo]
      cb = vacAlpha[index + lumo]
    else:
      vb = vacAlpha[index + homo]
      cb = vacAlpha[index + lumo]
    if odd and polarized:
      vba = vacAlpha[index + hAlpha]
      vbb = vacBeta[index + hBeta]
      cba = vacAlpha[index + lAlpha]
      cbb = vacBeta[index + lBeta]
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
  if odd and polarized:
    if hamaxkp == laminkp:
      print("Direct band gap (alpha)")# at", ef.kpoints[hamaxkp].Coord())
    else:
      print("Indirect band gap (alpha)")
    print("Gap (alpha):",'{:5.3f}'.format(lamin-hamax),"eV")
    if hbmaxkp == lbminkp:
      print("Direct band gap (beta)")# at", ef.kpoints[hbmaxkp].Coord())
    else:
      print("Indirect band gap (beta)")
    print("Gap  (beta):",'{:5.3f}'.format(lbmin-hbmax),"eV")
#    print("Forbidden alpha -> beta:",'{:5.3f}'.format(lbmin-hamax),"eV")
#    print("Forbidden beta -> alpha:",'{:5.3f}'.format(lamin-hbmax),"eV")
    print("Alpha VBE, CBE, Gap, VBW, CBW: \
    {0:5.3f}, {1:5.3f}, {2:5.3f}, {3:5.3f}, {4:5.3f} eV".format(hamax,lamin,
    lamin-hamax,hamax-hamin,lamax-lamin))
    print("Beta  VBE, CBE, Gap, VBW, CBW: \
    {0:5.3f}, {1:5.3f}, {2:5.3f}, {3:5.3f}, {4:5.3f} eV".format(hbmax,lbmin,
    lbmin-hbmax,hbmax-hbmin,lbmax-lbmin))
  else:
    if hmaxkp == lminkp:
      print("Direct band gap")# at", ef.kpoints[hmaxkp].Coord())
    else:
      print("Indirect band gap")
    print("Gap:",'{:5.3f}'.format(lmin-hmax),"eV")
    print("VBE, CBE, Gap, VBW, CBW: \
    {0:5.3f}, {1:5.3f}, {2:5.3f}, {3:5.3f}, {4:5.3f} eV".format(hmax,lmin,
    lmin-hmax,hmax-hmin,lmax-lmin))
#  print("Read Outcar?")
#  x = ReadOUTCAR()

def writeGNUPlotCtlbasic(nBands):
  with open('gnuctl.txt', 'w') as f:
    f.write("set key top\n")
    f.write('set xtics ("X" 0, "G" 10, "L" 20, "W" 30, "G" 40, "W" 41, "X" 51, "L" 61)\n')
    f.write('plot ')
    for i in range(nBands):
      if i == nBands-1:
        f.write('"Band_Data.txt" using 1:{0} notitle with linespoints\n'.format(i+2))
      else:
        f.write('"Band_Data.txt" using 1:{0} notitle with linespoints,\\\n'.format(i+2))
    f.write('pause -1')
  f.close()

def writeGNUPlotCtl(nBands, homo, lumo, reference = None):
  with open('gnuctl.txt', 'w') as f:
    f.write("# Warning this is set up for a specific set of k-points\n")
    f.write("# May be generalized in the future...\n")
    f.write("set terminal pngcairo size 800,600 enhanced font 'Arial,10'\n")
    f.write('set output "Bands.png"\n')
    f.write('set title "Band structure calculation"\n')
    f.write("set key top center horizontal spacing 1.01\n")
    f.write('set xtics ("X" 0, "G" 10, "L" 20, "W" 30, "G" 40, "W" 41, "X" 51, "L" 61)\n')
    f.write('set mytics 5\n')
    f.write('show mytics\n')
    f.write('set grid mytics xtics\n')
    f.write('set pointsize 0.2\n')
    f.write('# Draw left side\n')
    f.write('set xrange [ 0 : 40 ]\n')
    f.write('set yrange [ -20 : 10 ]\n')
    f.write('set xlabel "k-point"\n')
    f.write('set ylabel "Energy [eV]"\n')
    f.write('set multiplot\n')
    f.write('set size 0.66,1\n')
    f.write('set origin 0.0,0.0\n')
    f.write('set lmargin 5\nset rmargin 0\n')
    f.write('plot ')
    for i in range(nBands):
      if i == nBands-1:  # Last entry
        if not reference:
          f.write('"Band_Data_a.txt" using 1:{0} notitle with linespoints\n'.format(i+2))
        else:
          f.write('"Band_Data_a.txt" using 1:{0} notitle with linespoints,\\\n'.format(i+2))
          f.write('"-" using 1:2 title "E-fermi ref." with lines lt rgb "black"\n')
          f.write('# k-point-index  E-fermi\n')
          for i in range(7):
            if i == 6:
              f.write("{0:<2} {1:< 6.3f}\n".format(i*10 + 1, reference))
            else:
              f.write("{0:<2} {1:< 6.3f}\n".format(i*10, reference))
          f.write("end\n")
      elif i == homo:
        f.write('"Band_Data_a.txt" using 1:{0} with linespoints title "HOMO" lw 2 lt rgb "red",\\\n'.format(i+2))
      elif i == lumo:
        f.write('"Band_Data_a.txt" using 1:{0} with linespoints title "LUMO" lw 2 lt rgb "blue",\\\n'.format(i+2))
      else:
        f.write('"Band_Data_a.txt" using 1:{0} notitle with linespoints,\\\n'.format(i+2))
    # Draw right side
    f.write('# Draw right side\n')
    f.write('set title " "\n')
    f.write('set xrange [ 41 : 61 ]\n')
    f.write('set ylabel ""   # Remove axis label\n')
    f.write('set format y "" # Remove axis values\n')
    f.write('set size 0.34,1\n')
    f.write('set origin 0.66,0.0\n')
    f.write('set lmargin 3\nset rmargin 5\n')
    f.write('plot ')
    for i in range(nBands):
      if i == nBands-1:
        if not reference:
          f.write('"Band_Data_a.txt" using 1:{0} notitle with linespoints\n'.format(i+2))
        else:
          f.write('"Band_Data_a.txt" using 1:{0} notitle with linespoints,\\\n'.format(i+2))
          f.write('"-" using 1:2 title "E-fermi ref." with lines lt rgb "black"\n')
          f.write('# k-point-index  E-fermi\n')
          for i in range(7):
            if i == 6:
              f.write("{0:<2} {1:< 6.3f}\n".format(i*10 + 1, reference))
            else:
              f.write("{0:<2} {1:< 6.3f}\n".format(i*10, reference))
          f.write("end\n")
      elif i == homo:
        f.write('"Band_Data_a.txt" using 1:{0} with linespoints title "HOMO" lw 2 lt rgb "red",\\\n'.format(i+2))
      elif i == lumo:
        f.write('"Band_Data_a.txt" using 1:{0} with linespoints title "LUMO" lw 2 lt rgb "blue",\\\n'.format(i+2))
      else:
        f.write('"Band_Data_a.txt" using 1:{0} notitle with linespoints,\\\n'.format(i+2))
    f.write('set nomultiplot\n')
  f.close()
  print("gnuplot load 'gnuctl.txt'")
  
def printBandInfo(nElect, nBands, spinPolarized, bandData, betaBandData):
  print("Searching for highest fully occupied band.")
  for i in range(nBands):
    occ = 0.
    occb = None
    for kpt in bandData:
      occ += kpt[2][i]
    for kpt in betaBandData:
      if not occb:
        occb = kpt[2][i]
      else:
        occb += kpt[2][i]
    print("Band",i+1,"Occupancy","Alpha: {0:< 4.1f}".format(occ/float(len(bandData))), end = " ")
    if occb:
      print("Beta: {0:< 4.1f}".format(occb/float(len(betaBandData))))
    else:
      print("") # Missing new line if not printing beta

def debugPrintBands(bandData, vacAlpha, betaBandData, vacBeta):
  for i,(x,y) in enumerate(zip(bandData, betaBandData)):
    print("Kpoint",i)
    if len(x[1]) != len(y[1]):
      print("ERROR unequal number of bands read...")
      continue
    if len(x[1]) != len(x[2]):
      print("ERROR missing occupancies of spin comp. 1...")
      continue
    if len(y[1]) != len(y[2]):
      print("ERROR missing occupancies of spin comp. 2...")
      continue
    print(len(x[1]),"bands")
    for j in range(len(x[1])):
      print("Band: {0:<2}, A: {1:< 6.3f} {2:< 4.1f}, B: {3:< 6.3f} {4:< 4.1f}".
      format(j+1, x[1][j], x[2][j], y[1][j], y[2][j]))

def readKPTFile(filename):
  points = []
  with open(filename, 'r') as f:
    while True:
      l = f.readline()
      if not l:
        break
      if l.lstrip()[0] == '#': # Strip whitespace
        continue # Comment line
      a = l.split()
      if len(a) < 3 or (len(a[0])==1 and a[0]=='-'):
        points.append("-")
      else:
        points.append(vec3(l.split()[:3]))
  f.close()
  return points[:]

"""Give list of kpts(vec3) and
number of search points per line"""
def getTotalKPoints(points, linePoints):
  i = 0 ; start = 1 ; restart = False
  nPoints = min(1, len(points))
  while i < (len(points)-1):
    p0 = points[i]; p1 = points[i+1]
    if p0 == "-" or p1 == "-":
      i += 1
      start = 0
      restart = True
      continue
    elif restart:
      restart = False
    else:
      start = 1
    # start = 1 avoids double counting each listed kpt
    # nPoints is initially == 1, for the 1st point
    # Then continuously add linePoints from p0 to p1
    # If there is a break from '-' entry
    # start = 0, to acount for the 1st leading point
    # of a new set of points separate from the previous set
    #   If that makes sense
    for j in range(start,linePoints+1):
      nPoints += 1
    i += 1
  return nPoints

def Setup(args):
  print("Assuming you used vasp_runbands script, we'll run through the finished\n\
  calculations to assemble band structure data to one file you can graph...")
  print("WARNING: Currently outputs graphing for a specific band structure calculation.")
  print("Expecting commandline arguments:\n\
  # of k-points per calculation, # of k-points per line search,\n\
  # of k-points along which you wanted calculated\n\
  Can give -r to give a kpts file to read from (as with vasp_runbands.py).\n\
  \tNot yet used to define kpts in graph\n\
  e.g. vaspAssembleBands.py 10 10 -r kpts-fd-3m.txt")
  if len(args) < 4:
    print("Not enough command line arguments")
  debug=False
  max_kpoints = False
  line_kpoints = False
  totalPoints = 0
  i = 1
  while i < len(args):
    arg = args[i]
    if len (arg) == 2 and arg[0] == '-':
      if '-r' in arg:
        # Read input file
        infile = args[i+1]
        points = readKPTFile(infile)
        i += 1
      elif '-d' in arg:
        debug = True
      else:
        print("Unrecognized argument:",arg)
    else:
      if not max_kpoints:
        max_kpoints = int(arg)
      elif not line_kpoints:
        line_kpoints = int(arg)
      else:
        nPoints = int(arg)
        totalPoints = ((nPoints-1) * line_kpoints) + 1
    i += 1
  if totalPoints == 0:
    totalPoints = getTotalKPoints(points, line_kpoints)
  # Collect data
  reference_E_fermi = -6.060 # eV
  points = []; bandData = []
  eFermi = []; XC = []; alphaBeta = []
  nDirectories = math.ceil(totalPoints/max_kpoints)
  print("Reading through",nDirectories," directories.")
  nBands = 0
  spinPolarized = False; nElect = 0
  betaPoints = []; betaBandData = []
  vacAlpha = []; vacBeta = []
  startPoints = 0
  for i in range(nDirectories):
    filename = "{0:0={nrange}}/OUTCAR".format(i, nrange=math.floor(math.log(nDirectories)))
    try:
      f = open(filename, 'r')
    except:
      print("Failed to open",filename,"to read.")
      break
    target = " E-fermi :"
    target2 = "NBANDS="
    tLen = len(target)
    tLen2 = len(target2)
    while True:
      l = f.readline()
      if not l:
        break
      if l[:tLen] == target:
        print("Found E-fermi entry")
        a = l.split()
        eFermi.append(float(a[2]))
        XC.append(float(a[4]))
        if len(a) == 8: # Ends as 'alpha+bet : -#.####
          alphaBeta.append(float(a[-1]))
        else: # Ends as 'alpha+bet :-##.####
          alphaBeta.append(float(a[6][1:]))
        f.readline() # Blank
        if spinPolarized:
          f.readline() # Blank
          f.readline() # spin component 1
        while OUTCAR_ReadKPoint(f, points, bandData, nBands):
          # Just keep reading them....
          pass
        if spinPolarized:
          # Failed upon reading spin component 2
          while OUTCAR_ReadKPoint(f, betaPoints, betaBandData, nBands):
            # Just keep reading them....
            pass
      if target2 in l:
        a = l.split()
        nBands = int(a[14])
        print(nBands,"Bands.")
      if "ISPIN" in l:
        sp = int(l.split()[2])
        if sp == 2:
          spinPolarized = True
        elif sp == 1:
          spinPolarized = False
        else:
          print("Unrecognized ISPIN setting: ",l)
      if "NELECT" in l:
        nElect = int(float(l.split()[2]))
    for i in range(startPoints, len(bandData)): # kpt index - i
      for j in range(len(bandData[i][1])): # band energy index - j
        #print(i, j, bandData[i][1][j], alphaBeta[-1], bandData[i][1][j] + alphaBeta[-1])
        vacAlpha.append(bandData[i][1][j] + alphaBeta[-1])
        if spinPolarized:
          vacBeta.append(betaBandData[i][1][j] + alphaBeta[-1])
    startPoints = len(bandData)
  # Should have all band data...
  if debug:
    debugPrintBands(bandData, vacAlpha, betaBandData, vacBeta)
  print("Outputting data for graphing...")
  with open("Band_Data.txt",'w') as o:
    o.write("#Band Structure data from VASP line search...\n#")
    for i in points:
      o.write(str(i)+" ")
    o.write('\n#Kpoint Bands >')
    for i in range(nBands):
      o.write(str(i+1) + " ")
    o.write('\n')
    if spinPolarized:
      for p in range(len(bandData)): # for i,b in enumerate(bandData):
        o.write("{0:<3} ".format(p))
        index = p * nBands
        for i in range(nBands): # in b[1]:
          o.write("{0:< 10.3f}".format(vacBeta[index+i]))#j))
        o.write("\n")
    else:
      for i,b in enumerate(bandData):
        o.write("{0:<3} ".format(i))
        for j in b[1]:
          o.write("{0:< 10.3f}".format(j))
        o.write("\n")
  o.close()
  print("Alpha+bet",alphaBeta)
  with open("Band_Data_a.txt",'w') as o:
    o.write("#Band Structure data from VASP line search...\n#")
    for i in points:
      o.write(str(i)+" ")
    o.write('\n#Kpoint Bands >')
    for i in range(nBands):
      o.write(str(i+1) + " ")
    o.write('\n')
    for p in range(len(bandData)):
      o.write("{0:<3} ".format(p)) # K-point index from 0
      index = p * nBands
      for i in range(nBands):
        o.write("{0:< 10.3f}".format(vacAlpha[index+i]))
      o.write("\n")
  o.close()
  homo = 17; lumo = 18
  homo = nElect//2
  if nElect%2 == 1:
    homo += 1 # Odd electrons
  homo -= 1 # Index from 0
  lumo = homo+1
  writeGNUPlotCtl(nBands, homo, lumo, reference_E_fermi)
  print("E-fermi:", end=" ")
  for i,e in enumerate(eFermi):
    print(e + alphaBeta[i], end=" ")
  print("")
  printBandEdges(nElect, nBands, vacAlpha, vacBeta)
  printBandInfo(nElect, nBands, spinPolarized, bandData, betaBandData)
  print("Good Day!")

if __name__ == "__main__":
  import sys
  Setup(sys.argv)

