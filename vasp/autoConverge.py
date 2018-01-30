#!/bin/python3

import math, os, subprocess

# TODO
# Add option like check to keep running a
#   calculation until there is only one geometry step before moving
#   to handle variable volume calculations
# Add Tolerance and dP to the program commandline options
# Add option to select shell environment or change default to bash
# Maybe add options for non-cubic MP grids
# Add more convergence options
# Mimic the program to find minimum energy volume through SPE calculations
# **Implement restarting from part way through

# Math helpers
def Dot(a,b):
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
def Cross(a,b):
  # | a0 a1 a2 |
  # | b0 b1 b2 |
  # | i  j  k  |
  return [(a[1]*b[2]-a[2]*b[1]), -(a[0]*b[2]-a[2]*b[0]), (a[0]*b[1]-a[1]*b[0])]

def RunCalculation():
  f = open('vasp.out', 'w')
  # Stackexchange method to use aliases
  subprocess.call(["/bin/tcsh","-i","-c", "vasp"], env = os.environ,
  stdout=f, stderr=subprocess.STDOUT)
  # stdout=subprocess.DEVNULL) Hide output
  # stderr=subprocess.DEVNULL) Hide output of errors
  f.close()

# Call mkdir with subprocess
def MKDIR(name):
  subprocess.call(["mkdir", name], env = os.environ,
  stdout=subprocess.DEVNULL)
  
def CopyVASPIn(directory, origin=".", useCONTCAR=False):
  o = origin + "/"
  if subprocess.check_call(["cp", o+"INCAR", directory]):
    return False
  elif subprocess.check_call(["cp", o+"POTCAR", directory]):
    return False
  elif subprocess.check_call(["cp", o+"KPOINTS", directory]):
    return False
  if useCONTCAR:
    if subprocess.check_call(["cp", o+"CONTCAR", directory+"/"+"POSCAR"]):
      return False
  else:
    if subprocess.check_call(["cp", o+"POSCAR", directory]):
      return False
  return True

# Read for ENCUT
def GetENCUT():
  try:
    f = open("INCAR", 'r')
  except:
    print("ERROR: Failed to open INCAR")
    f.close()
    return 0
  try:
    while True:
      l = f.readline()
      if not l: # EOF
        break
      if "ENCUT" in l:
        l = l.split()
        if l[0][-1] == '=':
          ENCUT = float(l[1])
        else:
          ENCUT = float(l[2])
    f.close()
    return ENCUT
  except:
    print("ERROR: INCAR reading error")
    f.close()
    return 0
# sed old ENCUT
def SetENCUT(ENCUT):
  subprocess.call(["sed", "-i", "s/ENCUT.*/ENCUT = {0}/g".format(float(ENCUT)), "INCAR"], env = os.environ,
  stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
# Read for Monkhorst Pack grid
def GetKPT():
  try:
    f = open("KPOINTS", 'r')
  except:
    print("ERROR: Failed to open KPOINTS")
    f.close()
    return 1
  try:
    l = f.readline() # Automatic Mesh
    l = f.readline() # 0
    l = f.readline()
    if "ex" in l.lower():
      # Explicit kpoints
      f.close()
      return 1
    else:
      l = f.readline()
      f.close()
      return int(l.split()[0])
  except:
    print("ERROR: KPOINTS reading error")
    f.close()
    return 1
# Writes a new KPOINTS file
def SetKPT(KPT):
  try:
    f = open("KPOINTS", 'w')
    f.write("Automatic mesh\n0\nMonkhorst Pack\n")
    f.write("{0:>3}".format(int(KPT))*3 + "\n" + " 0."*3+"\n")
    f.close()
  except:
    print("ERROR: KPOINTS writing error")
    f.close()
    

#  return (EwE, eF, elapsed, NKPT,Ng,Nscf,Nelectrons) or (False)
def GetOUTCARData():
  try:
    f = open("OUTCAR", 'r')
  except:
    print("ERROR: Failed to open OUTCAR")
    f.close()
    return (False)
  try:
    EwE, eF, elapsed, NKPT, Ng, Nscf, Nelectrons = 0.,0,0.,0,0,0,0.
    while True:
      l = f.readline()
      if not l: # EOF
        break
      if "Elapsed time" in l:
        elapsed = float(l.split()[3])
      elif "energy  without entropy" in l:
        l = l.split()
        EwE = float(l[3])
        Ng += 1
      elif "energy without entropy" in l:
        Nscf += 1
      elif "NKPTS =" in l:
        NKPT = int(l.split()[3])
      elif "E-fermi" in l:
        l = l.split()
        eFerm = float(l[2])
        alphaBeta = 0.
        if len(l[6]) > 1:
          alphaBeta = float(l[6][1:])
        else:
          alphaBeta = float(l[7])
        eF = eFerm + alphaBeta
      elif "NELECT =" in l:
        Nelectrons = float(l.split()[2])
    # Wait for last read value
    f.close()
    return (EwE, eF, elapsed, NKPT,Ng,Nscf,Nelectrons)
  except IOError as err:
    print("ERROR: Failed to open OUTCAR or error reading", err)
    f.close()
    return 0.
  except ValueError as err:
    print("ERROR: Failed processing OUTCAR data")
    f.close()
    return 0.
  except:
    print("Unaccounted for Error reading OUTCAR")
    f.close()
    return 0.

def GetEnergyT0():
  try:
    f = open("OUTCAR", 'r')
    EwE = 0.
    while True:
      l = f.readline()
      if not l: # EOF
        break
      if "energy  without entropy" in l:
        l = l.split()
        EwE = float(l[3])
    # Wait for last read value
    f.close()
    return EwE
  except IOError as err:
    print("ERROR: Failed to open OUTCAR or error reading", err)
    f.close()
    return 0.
  except ValueError as err:
    print("ERROR: Failed processing OUTCAR data")
    f.close()
    return 0.
  except:
    print("Unaccounted for Error reading OUTCAR")
    f.close()
    return 0.

def GetFermiLevel():
  try:
    f = open("OUTCAR", 'r')
    eF = 0.
    alphaBeta = 0.
    while True:
      l = f.readline()
      if not l: # EOF
        break
      if "E-fermi" in l:
        l = l.split()
        eF = float(l[2])
        if len(l[6]) > 1:
          alphaBeta = float(l[6][1:])
        else:
          alphaBeta = float(l[7])
    # Wait for last read value
    f.close()
    return eF + alphaBeta
  except IOError as err:
    print("ERROR: Failed to open OUTCAR or error reading", err)
    f.close()
    return 0.
  except ValueError as err:
    print("ERROR: Failed processing OUTCAR data")
    f.close()
    return 0.
  except:
    print("Unaccounted for Error reading OUTCAR")
    f.close()
    return 0.

def GetVolume():
  try:
    f = open("CONTCAR", 'r')
    l = f.readline() # title
    scale = float(f.readline().split()[0]) # scale
    a = [float(x) for x in f.readline().split()]
    b = [float(x) for x in f.readline().split()]
    c = [float(x) for x in f.readline().split()]
    f.close()
    return Dot(Cross(a,b), c)
  except IOError as err:
    print("ERROR: Failed to open CONTCAR or error reading", err)
    f.close()
    return 0.
  except ValueError as err:
    print("ERROR: Failed processing CONTCAR data")
    f.close()
    return 0.
  except:
    print("Unaccounted for Error reading CONTCAR")
    f.close()
    return 0.

# Target is value to converge
# Vary is parameter to vary until convergence
# Check is a flag to double last calculation
# Tol is tolerance to consider converged
# dP is magnitude of parameter variation
#   0 will later set defaults 2 for MP and 25 for ENCUT
def Converge(target, vary, check, convRepeat=False, tol = 1e-2, dP = 0):
  # target
  #  1 = 'energy  without entropy'
  #  2 = 'E-fermi'
  #  3 = 'Volume'
  # vary
  #  1 = ENCUT
  #  2 = MP KPOINTS Grid
  R = None # Calculation Return
  if target == 1:
    R = GetEnergyT0
  elif target == 2:
    R = GetFermiLevel
  elif target == 3:
    R = GetVolume
  else:
    print("Bad input to converge target", target)
    return
  G = None # Calculation Get Parameter
  S = None # Calculation Set Parameter
  prefix = ""
  if vary == 1:
    G = GetENCUT
    S = SetENCUT
    prefix = "ENCUT-"
    if dP <= 0:
      dP = float(25)
    dP = float(dP)
  elif vary == 2:
    G = GetKPT
    S = SetKPT
    prefix = "MP-"
    if dP <= 0:
      dP = int(2)
    dP = int(dP)
  else:
    print("Bad input to converge parameter", vary)
    return
  # Long loop of calculations that may run for days
  parameter = G()
  lastVal = 0.
  ldir = "." # Last directory
  # First Calculation
  ndir = prefix+"{0}".format(parameter)
  MKDIR(ndir)
  CopyVASPIn(ndir, useCONTCAR=False)
  os.chdir(ndir)
  
  plabels = ["", "ENCUT", "KPT Grid"]
  label = ["", "Energy without entropy", "E-Fermi", "Volume"]
  unit = ["", "eV", "eV", "Ang^3"]
  
  print("Running first calculation...")
  print(plabels[vary] + ": " + str(parameter)) # e.g. "KPT Grid: 3"
  RunCalculation()
  lastVal = R()
  print(plabels[vary]+"\t\t"+label[target])                 # e.g. "KPT Grid         Volume"
  print("{0:<10}\t\t{1:>20.5f}".format(parameter, lastVal)) # e.g. "1                 100.00000"
  #TODO Clean this all up with recursive functions to loop calculations
  #  Inner loop for convRepeat
  #  Outer loop for value
  checkCount = 1
  checked = False
  conv = True
  if convRepeat:
    conv = False
  first = True
  while True:
    ## Move to original directory
    os.chdir("..")
    ldir = ndir
    ## Setup new calc.
     # New directory
    ndir = prefix+"{0}".format(parameter) 
    if checkCount > 0:
      ndir = ndir + "-{0}".format(checkCount)
     # mkdir $ndir
    MKDIR(ndir)
     # Copy files to new directory
    CopyVASPIn(ndir, origin = ldir, useCONTCAR=True)
     # cd $ndir
    os.chdir(ndir)
     # Set new parameter
    S(parameter)
    ## Run
    RunCalculation()
    ## Repeat calc if not a single step
    if convRepeat:
      results = GetOUTCARData()
      #  return (EwE, eF, elapsed, NKPT,Ng,Nscf,Nelectrons) or (False)
      if len(results) > 1:
        if results[4] > 1:
          checkCount += 1
          conv = False
          print("\t\tN Geometry Steps:", results[4])
          continue
        else:
          print("\tConverged value (1 geometry step)...")
          conv = True
      else:
        print("Calculation Failed? Check directory:",ndir)
        return
    ## Check Result
    Val = R()
    dV = lastVal - Val
    lastVal = Val
    print("{0:<10}\t\t{1:>20.5f}".format(parameter, Val))
    if abs(dV) <= tol and not first:
      if not conv:
        continue
      if check and not checked:
        checkCount += 1
        checked = True
      else:
        print(label[target],"converged to",Val,unit[target],"with delta:", dV)
        return
    else:
      checked = False
      first = False
      parameter += dP
      if checkCount > 0:
        checkCount = 0

def FileExists(filename):
  try:
    f = open(filename, 'r')
    f.seek(0,2)
    if f.tell() == 0:
      print("Empty file:",filename)
      f.close()
      return False
  except IOError as e:
    print("Can't find",filename)
    f.close()
    return False
  except:
    print("Error reading",filename)
    f.close()
    return False
  f.close()
  return True
  
def ValidateDirectory():
  return FileExists("INCAR") and FileExists("POSCAR") and FileExists("KPOINTS") and FileExists("POTCAR")
  
def PrintHelp():
  print("\
autoConverge.py v0.1 for VASP 4.3?-5.3?\nIntended to run a series of vasp calculations\n\
until a desired parameter is converged.\n\
Should be run in a new directory containing vasp input files.\n\
You should also have a 'vasp' alias setup to run vasp.\n\
NOTE: The vasp alias should be in ~/.cshrc, it is run in /bin/tcsh.")
  print("\
  Options:\n\
    -EwE   - Converge Energy without entropy value\n\
    -EF    - Converge E-fermi value (E-fermi + alpha+bet)\n\
    -Vol   - Converge axb.c (Triple Scalar Product of lattice parameters)\n\
    -ENCUT - Vary ENCUT to Converge\n\
    -KPT   - Vary MP kpoint grid size.\n\
    -C     - Confirm convergence by repeating last calculation (Volume)\n\
    -C2    - Repeats all calculations until a single geometry step is taken\n\
             Important for Volume relaxation with constant Planewaves\n\
    -t #   - Sets tolerance to consider converged\n\
             e.g. -t 0.001 with -EwE will converge Enthalpy to within 1 meV.\n\
             Default: 0.001\n\
    -dp #  - Sets value to vary parameter on convergence search.\n\
             Default: 25 eV (ENCUT) or 2 KPTs (Odd grid search 1,3,5,...)\n")

def Run(args):
  # Check arguments
  if len(args)<3 or ' -h' in args:
    PrintHelp()
  target = 0
  vary = 0
  convRepeat = False
  check = False
  tolerance = 0.001
  dP = 0.
  for i,a in enumerate(args):
    if a == '-t' and len(args)>i+1:
      tolerance = float(args[i+1])
    elif a == '-EwE':
      target = 1
    elif a == '-EF':
      target = 2
    elif a == '-Vol':
      target = 3
    elif a == '-ENCUT':
      vary = 1
    elif a == '-KPT':
      vary = 2
    elif a == '-C2':
      convRepeat = True
    elif a == '-C':
      check = True
    elif a == '-dp' and len(args)>i+1:
      dP = float(args[i+1])
  if target == 0 or vary == 0:
    print("Invalid Convergence Selection")
    PrintHelp()
    return
  # Now run program
  if ValidateDirectory():
    Converge(target, vary, check, convRepeat=convRepeat, tol=tolerance, dP=dP)
  print("Good Day!")

if __name__ == "__main__":
  import sys
  Run(sys.argv)
