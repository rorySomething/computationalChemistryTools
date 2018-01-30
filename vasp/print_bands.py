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

def Print_Bands():
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
  try:
    f = open("Bands.txt", 'w')
  except:
    print("Can't open Bands.txt for writing")
    return
  for i in range(ef.bands):
    s = "{0}".format(i+1)
    for j in ef.kpoints:
      if not ef.polarized:
        s = s + " {0: >6.3f} ".format(j.bands[i])
      else:
        s = s + " {0: >6.3f} ".format(j.bands[0][i])
    s = s + '\n'
    print(s)
    f.write(s)
  #print("Read Outcar?")
  #ReadOUTCAR()
  f.close()

if __name__ == "__main__":
  Print_Bands()

