#!/bin/python3

from math import sqrt
import math
import subprocess

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
  def __add__(self,other):
    return vec3([self.x+other.x,self.y+other.y,self.z+other.z])
  def __mul__(self,other):
    return vec3([self.x*float(other),self.y*float(other),self.z*float(other)])
  def __truediv__(self,other):
    return vec3([self.x/float(other),self.y/float(other),self.z/float(other)])

  
def WriteRunScript(total):
  print("Writing a linux script to run all of the jobs...")
  print("Assuming you have an alias vasp to run the job automatically.")
  try:
    f = open("runBandLines", 'w')
  except:
    print("Failed to open runBandLines for writing")
    return
  f.write("#!/bin/tcsh\n")
  f.write("source ~/.cshrc\n")
  for i in range(total):
    directory = "{0:0={nrange}}".format(i, nrange=math.floor(math.log(total)))
    if i == 0:
      f.write("cd "+directory+"\n")
    else:
      f.write("cd ../"+directory+"\n")
    f.write("vasp\n")
  f.write("cd ..\n")
  f.close()
  subprocess.call(["chmod", "+x", "runBandLines"], stdout=subprocess.DEVNULL)
  # Bash
  try:
    f = open("bashRunBandLines", 'w')
  except:
    print("Failed to open bashRunBandLines for writing")
    return
  f.write("#!/bin/bash\n")
  f.write("shopt -s expand_aliases\n")
  f.write("source ~/.bashrc\n")
  for i in range(total):
    directory = "{0:0={nrange}}".format(i, nrange=math.floor(math.log(total)))
    if i == 0:
      f.write("cd "+directory+"\n")
    else:
      f.write("cd ../"+directory+"\n")
    f.write("vasp\n")
  f.write("cd ..\n")
  f.close()
  subprocess.call(["chmod", "+x", "bashRunBandLines"], stdout=subprocess.DEVNULL)

def GetBandINCARCopy():
  incar = "# Band structure from vasp_runbands.py\n"
  with open('INCAR', 'r') as i:
    sigma = False
    while True:
      l = i.readline()
      if not l:
        break # EOF
      if "NSW" in l:
        incar += "NSW = 0\n"
      elif "ICHARG" in l:
        incar += "ICHARG = 11\n"
      elif "LWAVE" in l:
        incar += "LWAVE = .FALSE.\n"
      elif "LCHARG" in l:
        incar += "LCHARG = .FALSE.\n"
      elif "LAECHG" in l:
        incar += "LAECHG = .FALSE.\n"
      elif "ISMEAR" in l:
        incar += "ISMEAR = 2\n"
        if "SIGMA" in l and not sigma:
          a = l.split()
          for n,s in enumerate(a):
            if "SIGMA" in s:
              incar += "SIGMA = " + a[n+2] + "\n"
              sigma = True
              break
          if not sigma:
            print("Error reading SIGMA entry")
            incar += "SIGMA = 0.01\n"
            sigma = True
      elif "SIGMA" in l and not sigma:
        a = l.split()
        for n,s in enumerate(a):
          if "SIGMA" in s:
            incar += "SIGMA = " + a[n+2] + "\n"
            sigma = True
            break
        if not sigma:
          print("Error reading SIGMA entry")
          incar += "SIGMA = 0.01\n"
          sigma = True
      else:
        incar += (l)
    if not sigma:
      incar += "SIGMA = 0.01\n"
  i.close()
  return incar

def PrepareCalculations(all_points, kmax):
  print("Assuming you have your INCAR, POTCAR, POSCAR and all important CHGCAR in this directory...")
  total = math.ceil(len(all_points)/kmax)
  print("Making",total,"subdirectories for all of these calculations.")
  for i in range(total):
    subprocess.call(["mkdir", "{0:0={nrange}}".format(i, nrange=math.floor(math.log(total)))],
    stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
  try:
    f = open("POSCAR", 'r')
    print("Attempting to copy POSCAR to each subdirectory...")
    for i in range(total):
      subprocess.call(["cp", "POSCAR", "{0:0={nrange}}".format(i, nrange=math.floor(math.log(total)))],
      stdout=subprocess.DEVNULL)
    f.close()
  except:
    try:
      f = open("CONTCAR", 'r')
      print("Attempting to copy CONTCAR to POSCAR in each subdirectory...")
      for i in range(total):
        subprocess.call(["cp", "CONTCAR", "{0:0={nrange}}/POSCAR".format(i, nrange=math.floor(math.log(total)))], stdout=subprocess.DEVNULL)
      f.close()
    except:
      print("ERROR: No POSCAR or CONTCAR available...")
      return
  try:
    print("Attempting to copy INCAR to each subdirectory...")
    incar = GetBandINCARCopy()
    print(incar)
    for i in range(total):
      with open("{0:0={nrange}}/INCAR".format(i, nrange=math.floor(math.log(total))), 'w') as out:
        out.write(incar)
      out.close()
#      subprocess.call(["cp", "INCAR", "{0:0={nrange}}".format(i, nrange=math.floor(math.log(total)))],
#      stdout=subprocess.DEVNULL)
  except:
    print("ERROR: No INCAR available...")
    return
  try:
    f = open("POTCAR", 'r')
    print("Attempting to copy POTCAR to each subdirectory...")
    for i in range(total):
      subprocess.call(["cp", "POTCAR", "{0:0={nrange}}".format(i, nrange=math.floor(math.log(total)))],
      stdout=subprocess.DEVNULL)
    f.close()
  except:
    print("ERROR: No POTCAR available...")
    return
  try:
    f = open("CHGCAR", 'r')
    print("Attempting to copy CHGCAR to each subdirectory...")
    for i in range(total):
      subprocess.call(["cp", "CHGCAR", "{0:0={nrange}}".format(i, nrange=math.floor(math.log(total)))],
      stdout=subprocess.DEVNULL)
    f.close()
  except:
    print("ERROR: No CHGCAR available...")
    return
  print("Now the KPOINT file generation.")
  p = 0 # Current kpoint
  for i in range(total):
      print('File',i+1,'of',total)
      filename = "{0:0={nrange}}/KPOINTS".format(i, nrange=math.floor(math.log(total)))
      try:
        f = open(filename, 'w')
      except:
        print("Can't open file",filename," for writing.")
        return
      f.write("Explicit k-points\n")
      numberOfPoints = min(len(all_points) - p, kmax)
      k = 0 # Kpoints in this calculation
      f.write(str(int(numberOfPoints))+"\n") # Number of points in this run
      f.write("Reciprocal") # Assuming points were given in reciprocal space
      while k < numberOfPoints:
        print('Point', k+1,'of',numberOfPoints)
        f.write('\n')
        # Assuming weights are all equal
        f.write("{0: <10.7f} {1: <10.7f} {2: <10.7f} {3: <10.6f}".format(
        all_points[p].x, all_points[p].y, all_points[p].z, 1.0))
        p += 1
        k += 1
  f.close()
  print("All done.")
  WriteRunScript(total)
  
def OutputExample():
  print("A '-' on a line will start a new search with a discontinuity.")
  print("The first 3 numbers on each line is read as a point.")
  print("Everything else in the file will be ignored.")
  with open("Example-kpoints.txt", 'w') as f:
    f.write("0.0 0.5 0.0\n")
    f.write("0.0 0.0 0.0\n")
    f.write("0.25 0.25 0.25\n")
    f.write("0.25 0.5 0.0\n")
    f.write("0.0 0.0 0.0\n")
    f.write("- Break in line search\n")
    f.write("0.25 0.5  0.0  W\n")
    f.write("0.0  0.5  0.0  X\n")
    f.write("0.25 0.25 0.25 L\n")
  f.close()
  print("File written: Example-kpoints.txt")

def Setup(args):
  print("We're going to set up an extended VASP band structure calculation...")
  print("Expecting commandline arguments:\n\
  # of k-points per calculation, # of k-points per line search,\n\
  k-points along which to calculate the structure.\n\
  Option: -r - specify a file with desired k-points for line search\n\
  Option: -o - output an example k-points file")
  print("Example: runbands.py 10 10 -r kpoints.txt")
  print("Example: runbands.py 10 10 0.0 0.0 0.0 0.5 0.0 0.0 0.25 0.25 0.25")
  if (len(args) < 5):
    for i in args:
      if '-o' in i:
        OutputExample()
        return
    print("Not enough command line arguments")
    return
  max_kpoints = int(args[1])
  line_kpoints = int(args[2])
  points = []
  if '-r' in args[3]:
    #Read input from file
    infile = args[4]
    with open(infile, 'r') as f:
      while True:
        l = f.readline()
        if not l:
          break
        a = l.split()
        if len(a[0]) == 1 and a[0][0] == '-':
          points.append("-")
        else:
          points.append(vec3(l.split()[:3]))
  else:
    i = 3
    while i < (len(args)-2):
      points.append(vec3(args[i:i+3]))
      print(points[-1])
      i += 3
  all_points = []
  all_points.append(points[0])
  i = 0 ; start = 1 ; restart = False
  while i < (len(points)-1):
    p0 = points[i]; p1 = points[i+1]
    if p0 == "-" or p1 == "-":
      i += 1
      start = 0
      restart = True
      continue
    else:
      if restart:
        restart = False
      else:
        start = 1
    scale = (p1 - p0)/float(line_kpoints)
    for j in range(start,line_kpoints+1):
      all_points.append(p0+scale*float(j))
    i += 1
  print("All k-points along desired path.")
  for i in all_points:
    print(i)
  PrepareCalculations(all_points, max_kpoints)
  print("Good Day!")

if __name__ == "__main__":
  import sys
  Setup(sys.argv)

