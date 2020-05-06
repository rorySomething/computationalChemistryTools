#!/bin/python3

"""
autoConverge.py

Rory Vander Valk 2016 - 2020

Automatically runs VASP calculations in a series.
Relies heavily on linux sed command.
Need the user to set an alias in ~/.cshrc to run VASP.
  e.g. alias vasp 'mpirun -np 8 vasp_std'
Or use sed on the script for now for your desired shell
  sed -i "s|bin/tcsh|bin/bash|g" autoConverge.py
It will make a series of subdirectories for each calculation.
Given valid INCAR, KPOINTS, POSCAR, POTCAR it will vary chosen
setting ENCUT, KPOINTS grid, vacuum length (c axis).
These setting are varied until a chosen property converges
-Energy  without entropy from the OUTCAR
-E-Fermi (E-fermi + alpha+bet) from the OUTCAR
-Volume of CONTCAR lattice axb.c

Plenty of room for clean up.
"""

__all__ = [Converge]
__version__ = "0.9999999"  # Almost...
__author__ = "Rory Vander Valk"

import math
import os
import subprocess
import argparse

# TODO
# Run through a linter?
# Incorporate checks
#   Check that the user has alias vasp in their shell
#   Make sure KPOINTS grid is X X 1 for surface calculations
#   Make sure ISIF=2 for -Vacuum calculations
#   If ISIF=2 or other constant volume settings, ignore -C2 and -C only need to optimize the 1 time
# Add option to select shell environment or change default to bash
# Add options for non-cubic MP grids
#   -surface flag to ensure -kpt keeps grid as X X 1 as it increases X
# Add more convergence options
# Mimic the program to find minimum energy volume through SPE calculations
# **Implement restarting from part way through
# Sort of handled by -min, -dP options

# Math helpers
def Dot(a, b):
  """Dot product of 2 lists containing 3 floats each."""
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
def Cross(a, b):
  """Cross product of 2 lists containing 3 floats each."""
  # | a0 a1 a2 |
  # | b0 b1 b2 |
  # | i  j  k  |
  return [(a[1]*b[2]-a[2]*b[1]), -(a[0]*b[2]-a[2]*b[0]), (a[0]*b[1]-a[1]*b[0])]


def RunCalculation():
  """Call 'vasp' alias in a tcsh shell.
  
  Assumes we are in a directory containing valid VASP input files.
  Pipes output to 'vasp.out' text file.
  -i makes the shell interactive
    I think this was necessary to load alias from resource files
  -c starts the command
  """
  f = open('vasp.out', 'w')
  # Stackexchange method to use aliases
  subprocess.call(["/bin/tcsh","-i","-c", "vasp"], env = os.environ,
  stdout=f, stderr=subprocess.STDOUT)
  # stdout=subprocess.DEVNULL) Hide output
  # stderr=subprocess.DEVNULL) Hide output of errors
  f.close()


def MKDIR(name):
  """Calls mkdir via subprocess."""
  subprocess.call(["mkdir", name], env = os.environ,
  stdout=subprocess.DEVNULL)


def CopyVASPIn(directory, origin=".", useCONTCAR=False):
  """Copy VASP input files from one directory to another.
  
  Returns False if any files fail to copy."""
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



# Could probably make and reuse a function that just returns
# a dict of all set tags in an INCAR file for many VASP related scripts
def GetENCUT():
  """Read INCAR file for ENCUT tag setting."""
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
    # Or float casting error...
    f.close()
    return 0


def SetENCUT(ENCUT):
  """Calls sed to update INCAR file ENCUT tag."""
  subprocess.call(["sed", "-i", "s/ENCUT.*/ENCUT = {0}/g".format(float(ENCUT)), "INCAR"], env = os.environ,
  stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)


def GetISIF():
  """Read INCAR file for ISIF tag."""
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
      if "ISIF" in l:
        # ISIF=2 or ISIF = 2 # Comment
        l = l.split('=')
        isif = l[1]
        if '#' in isif:
          isif = l.split('#')[0]
        break
    f.close()
    return int(isif.strip())
  except:
    print("ERROR: INCAR reading error")
    f.close()
    return 0


def SetISIF(constantVolume=True):
  """Conditionally call sed to update ISIF in an INCAR file.
  
  Set ISIF to an equivalent setting based on the parameter constantVolume.
  If the calculation is set ISIF=3, but the converging property should be
  constant volume it will set the INCAR to ISIF=2.
  If we are changing KPOINTS to converge Volume,
  the calculation should not be constant volume (ISIF=2),
  so we'd set it to 3.
  https://www.vasp.at/wiki/index.php/ISIF
  """
  isif = GetISIF()
  newISIF = isif
  if constantVolume:
    # Want constant volume, ISIF = 0,1,2,4,5
    if isif in [6,7]:
      newISIF = 5
    elif isif == 3:
      newISIF = 2 # Or 4 to allow cell shape to change?
  else:
    # Want variable volume, ISIF = 3,6,7
    if isif == 5:
      newISIF = 6
    if isif in [0,1,2,4]:
      newISIF = 3
  # Only edit file if we need to update it
  if isif != newISIF:
    subprocess.call(["sed", "-i", "s/ISIF.*/ISIF = {0}/g".format(float(newISIF)), "INCAR"], env = os.environ,
    stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)


def GetKPT():
  """Read KPOINTS file for grid settings.
  
  Only uses the first value for the grid size.
  So a 'weird' setting of 1 2 2 will be read as 1.
  Ignores whether it is set as Gamma or Monkhorst-Pack.
  """
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


def SetKPT(KPT, gamma=True, surface=False):
  """Writes/overwrites a new KPOINTS file.
  
  KPT - integer, grid size
  gamma - bool, is not yet called outside of this function so it's always true.
    gamma sets the grid type to Gamma or Monkhorst-Pack
  surface - bool, determines if the last grid point should stay 1, default: False
    Generally assuming a surface calculation keeps the vacuum along the c lattice vector.
  """
  try:
    f = open("KPOINTS", 'w')
    if gamma:
      f.write("Automatic mesh\n0\nGamma\n")
    else:
      f.write("Automatic mesh\n0\nMonkhorst Pack\n")
    if surface:
      f.write("{0:>3}".format(int(KPT))*2 \
      + "  1\n 0. 0. 0.\n")
    else:
      f.write("{0:>3}".format(int(KPT))*3 + "\n" + " 0."*3+"\n")
    f.close()
  except:
    print("ERROR: KPOINTS writing error")
    f.close()


# TODO
def getENMIN(debug=False):
  """Return max(ENMIN for each species) from POTCAR"""
  enmin = 0.
  try:
    f = open("POTCAR", 'r')
    while True:
      l = f.readline()
      if "ENMIN" in l:
        enmin = max(enmin, float(l.split()[1]))
        f.seek(500) # Figure out how to skip most of the file
  except:
    print("Unaccounted for Error reading POTCAR")
    f.close()
  return enmin


def GetScaleValue(debug=False, file="CONTCAR"):
  """Return lattice scale from a VASP POSCAR/CONTCAR file"""
  try:
    f = open(file, 'r')
    l = f.readline() # title
    scale = float(f.readline().split()[0]) # scale
    f.close()
    return scale
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


def GetCVectorString(debug=False, file="CONTCAR"):
  """Returns the 5th line of a VASP POSCAR/CONTCAR file.
  
  Presently a quick hack job from GetLengthC function.
  """
  try:
    f = open(file, 'r')
    l = f.readline() # title
    scale = float(f.readline().split()[0]) # scale
    a = [float(x) for x in f.readline().split()]
    b = [float(x) for x in f.readline().split()]
    c = f.readline()
    f.close()
    return c
  except IOError as err:
    print("ERROR: Failed to open CONTCAR or error reading", err)
    f.close()
    return ""
  except ValueError as err:
    print("ERROR: Failed processing CONTCAR data")
    f.close()
    return ""
  except:
    print("Unaccounted for Error reading CONTCAR")
    f.close()
    return ""


def SetLengthC(length, debug=False):
  """Set the c axis in a POSCAR file to our chosen length via sed."""
  # String to replace
  poscarC = GetCVectorString(debug, file = "POSCAR")
  # Get existing length
  cLength = GetLengthC(debug, file = "POSCAR")
  # length parameter is target length, scale will change the lattice vector in the file,
  # but also changes the vector length reported from GetLengthC
  # So we don't need to consider the scale value
  # Need to multiply length by a normalized c vector
  # Normalize by dividing each component by the length
  # Then multiplying each component by the target length to get final vector
  C = [float(x)/cLength*length for x in poscarC.split()]
  # TODO Work out proper float length formatting in POSCAR
  newString = "{:>12.6f}{:>12.6f}{:>12.6f}".format(*C)
  # Set new length by editing the POSCAR in place with sed
  subprocess.call(
    ["sed", "-i", "s/{}/{}/g".format(poscarC, newString), "INCAR"], env = os.environ,
  stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)


#  return (EwE, eF, elapsed, NKPT,Ng,Nscf,Nelectrons) or (False,)
def GetOUTCARData():
  """Returns a tuple of data from VASP OUTCAR file, or (False,).
  
  Reads through the VASP OUTCAR file looking for keyword's values
  and counting SCF and geometry steps.
  Only fails from IO errors, doesn't currently care if the calculation failed or data is missing.
  On success returns a tuple of floats and ints:
    (Energy  without Entropy, E-Fermi, Elasped Time,
    NKPT, N geometry steps, N SCF steps, N electrons)
  E-Fermi is taken as 'E-fermi'+'alpha+bet'
  On failure returns (False,)
  """
  try:
    f = open("OUTCAR", 'r')
  except:
    print("ERROR: Failed to open OUTCAR")
    f.close()
    return (False,)
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
    return (False,)
  except ValueError as err:
    print("ERROR: Failed processing OUTCAR data")
    f.close()
    return (False,)
  except:
    print("Unaccounted for Error reading OUTCAR")
    f.close()
    return (False,)


def GetEnergyT0():
  """Reads VASP OUTCAR and returns value of final 'energy  without entropy'."""
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
  """Reads VASP OUTCAR and returns value of final 'E-fermi'+'alpha+bet'."""
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


def GetVolume(filename="CONTCAR"):
  """Return volume of VASP CONTCAR/POSCAR file. Lattice axb.c."""
  try:
    f = open(filename, 'r')
    l = f.readline() # title
    scale = float(f.readline().split()[0]) # scale
    a = [float(x)*scale for x in f.readline().split()]
    b = [float(x)*scale for x in f.readline().split()]
    c = [float(x)*scale for x in f.readline().split()]
    f.close()
    return Dot(Cross(a,b), c)
  except IOError as err:
    print("ERROR: Failed to open {} or error reading".format(filename), str(err))
    f.close()
    return 0.
  except ValueError as err:
    print("ERROR: Failed processing {} data".format(filename), str(err))
    f.close()
    return 0.
  except:
    print("Unaccounted for Error reading {}".format(filename))
    f.close()
    return 0.


def GetLengthC(debug=False, file="CONTCAR"):
  """Return length of the c axis from a VASP POSCAR/CONTCAR file."""
  try:
    f = open(file, 'r')
    l = f.readline() # title
    scale = float(f.readline().split()[0]) # scale
    a = [float(x) for x in f.readline().split()]
    b = [float(x) for x in f.readline().split()]
    c = [float(x) for x in f.readline().split()]
    f.close()
    return math.sqrt(Dot(c, c)) * scale
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



inf = float('inf') # For setting default maximum, it's overkill
def Converge(target, vary, dP = 0, minimum=0, maximum=inf, tolerance=1e-2,
    check=False, convRepeat=False, surface=False, debug=False):
  """Run VASP calculations in a loop until a chosen value is converged.
  
  Note: Use named parameters to keep compatible function calls if the
  parameter order is rearranged and/or more options are added in the future...
  As if this function will be used in any complex way...
  
  The workhorse function, creates named subdirectories for each varied
  calculation parameter and starts the VASP calculation in those
  subdirectories.
  Along the way it will print simple calculation results to the terminal for
  each step that it has completed.
  TODO Add timing, and a timelimit parameter
  It will keep going until the difference in the results of 2 calculations are
  within a tolerance value, the parameters reach a prescribed maximum value, you cancel the command, or the your
  computer fails (power outage, network outage, Segmentation Fault, etc...).
  You can ideally restart the calculations by choosing the minimum value for
  the next incompleted step.  It is not yet tested but it should overwrite
  an old incomplete calculation without any hiccups.
  Be sure to specify surface=True for a slab calculation to ensure the
  KPOINTS are kept to 1 along the c axis and the calculation is
  performed with constant volume.
  NOTE: Some combinations of parameters are invalid, like varying vacuum
  length and converging cell volume...
  or doing any surface=True calculation while trying to converge volume.
  That would be the work of a different program to vary slab a and b vectors
  while keeping vacuum space present.
  Optimizing a and b vectors based on the stress tensors in the cell.
  
  
  Parameters
  
  target - integer, enumerates system property we want to converge.
    0 - uninitialized, unaccepted
    1 - 'Energy  without entropy' from the OUTCAR
        This is the energy at 0K.
    2 - 'E-fermi'+'alpha+bet' from the OUTCAR
        This is an attempt to get vacuum referenced fermi levels.
        I'm not sure if adding alpha+bet is doing that correctly
        but it seems okay in my silicon calculations.
    3 - Volume from CONTCAR
        Triple scalar product of the lattice vectors.
        (A cross B) dot C
        Keep varying parameters until the system geometry is completely stable (Up to your chosen EDIFF/EDIFFG).
  
  vary - integer, enumerates input parameter to vary.
      Change the parameter to see its effect on the calculated ground state.
    0 - uninitialized, unaccepted
    1 - ENCUT, INCAR tag
          As the ENCUT increases more planewaves are included in your given
          cell volume.  The calculation then takes longer but the wavefunction
          gains more precision and possibly accuracy.
    2 - KPOINT grid, KPOINTS file
          The grid sampling the wavefunction in the reciprocal lattice.
          More points means more sampling locations giving greater accuracy.
          More points increases calculation time, not quite linearly because
          by default they are symmetry reduced so 2 2 2 might only give
          2 points instead of 8 in a high symmetry cell.
          It has large effects on the cell if there are no points at 0. 0. 0.
          0. 0. 0. corresponds to the infinite repeating bulk in real space.
          Gamma centered or odd numbered Monkhorst-Pack should be preferred.
          I haven't encountered cases where even MP grids worked out well.
          You can easily test with this program... Once I include the
          gamma flag in this function to call SetKPT.
    3 - Lattice vector c length
          Surface calculations involve looking at a slab of material.
          So 2 axis remain periodic and the 3rd contains vacuum space.
          Here we assume the c axis is the vacuum direction.
          We extend the cell gradually along the c axis to find the smallest
          amount of vacuum that leaves the wavefunction or ground state
          energy unchanged.  So if we put a molecule on the surface we want to
          minimize interactions wrapping around to the bottom of the slab.
          The larger the cell and vacuum the more planewaves you can fit
          with a given energy below ENCUT which makes the calculation take
          longer.  We want to find the smallest cell we can use for long
          term band structure, adsorption, transition state calculations, etc.
  
  dP - int of float depending on choice of vary
    This is the amount to increment your parameter between steps.
    Defaults
      ENCUT - 25, 25 eV, arbitrarily chosen
      KPT - 2, gives odd MP grid 1, 3, 5, ...
      Vacuum - 10% of the initial c length, arbitrarily chosen
        10% seems like a nice round number if no other is given.
        Could be better as 1-3 Angstroms.
  
  minimum - int or float depending on choice of vary
    This is the lowest starting point of the calculations so it can be set
    high to use as a restart from previous series of calculations.
    If you stopped with 1 3 5 7 you can start up again with minimum = 9
    Defaults
      ENCUT - ENMIN from POTCAR
      KPT - 1
      Vacuum - POSCAR c vector length (The initial structure)
  
  maximum - int or float depending on choice of vary
    Set this to stop the calculations at some point even if they don't
    converge.  Like if it is a large cell and KPT 5 would be infeasible,
    set maximum to 3 or 4.  Likewise an ENCUT of 2000 likely doesn't make
    sense because it's way outside the limits of precision of the data
    used to creation the pseudopotentials.
    Default is inf
  
  tolerance - float, tolerance before considering two points equal,
    giving a negligible difference between them.
    So when the change in a value between two steps is less than tolerance the
    calculation is considered complete.
    Default is 1e-2
      Corresponding to 10 meV of precision or .01 cubic Angstroms.
  
  check - bool, do the final calculation twice to confirm it is the final
          state.  It only really makes sense to do it if you are monitoring
          the change in volume.
          As the volume changes, the basis set changes, because the planewaves
          that fit in the new volume change.  They do eventually converge
          but sometimes it takes a while for a cell to relax fully.
          Especially when you change functionals.
          PW91 may give a=b=c=4.1 while PBE gives a=b=c=4.9
          Each VASP optimization starts with a fixed set of planewaves
          as the calculation continues the volume hits an optimization limit.
          The next calculation starts with a new set of planewaves and changes
          the volume further, so it may change over numerous calculations.
          e.g. 4.1 -> 4.5, 4.5 -> 4.8, 4.8 -> 4.9
          Default is False
  
  convRepeat - bool, keep doing repeated calculations at each parameter step
    until only 1 geometry step is taken in the calculation.  This will ensure the value you get for each setting is the absolute precise ground structure.  This is again related to the change in volume.
    Default is False
  
  surface - bool, Flag to tell if this is a surface/slab calculation.
    If it is a surface calculation the k-point grid is kept as N N 1
    So we don't waste time sampling vacuum space.
    It is unnecessary to have k-points along the vacuum axis because of the
    large number of planewaves that fit in a large vacuum axis length.
  
  debug - bool, print out a larger wall of text
      Not very well or always implemented.
    Default is False
  """
  R = None # Calculation Return
  targetFormatString = ""
  prec = GetPrintPrecision(tolerance)
  if target == 1:
    R = GetEnergyT0
    targetFormatString = "{0}:"+">{width}.{prec}f".format(width=prec+5, prec=prec)
    # {0}:7.2f -> '-182.36 eV'
  elif target == 2:
    R = GetFermiLevel
    targetFormatString = "{0}:"+">{width}.{prec}f".format(width=prec+4, prec=prec)
    # {0}:6.2f -> ' -6.06 eV'
  elif target == 3:
    R = GetVolume
    targetFormatString = "{0}:"+">{width}.{prec}f".format(width=prec+6, prec=prec)
    # {0}:8.2f -> '10597.06 Å³'
  else:
    print("Bad input to converge target", target)
    return
  G = None # Calculation Get Parameter
  S = None # Calculation Set Parameter
  prefix = ""
  parameterFormatString = ""
  if vary == 1:
    G = GetENCUT
    S = SetENCUT
    prefix = "ENCUT-"
    parameterFormatString = "{0:.1f}"
    if dP <= 0:
      dP = float(25)
    dP = float(dP)
    if minimum <= 0:
      minimum = getENMIN()
    minimum = float(minimum)
  elif vary == 2:
    G = GetKPT
    S = SetKPT
    prefix = "MP-"
    parameterFormatString = "{0:d}"
    if dP <= 0:
      dP = int(2)
    dP = int(dP)
    if minimum <= 0:
      minimum = 1
    minimum = int(minimum)
  elif vary == 3:
    G = GetLengthC
    S = SetLengthC
    prefix = "Vacuum-"
    parameterFormatString = "{0:.3f}"
    c0 = GetLengthC(debug, file="POSCAR")
    if dP == 0:
      dP = float(c0*0.1) # 10% c0
    dP = float(dP)
    if minimum == 0:
      minimum = c0
    minimum = float(minimum)
  else:
    print("Bad input to converge parameter", vary)
    return
  # Long loop of calculations that may run for days
  # Start from minimum value or starting value
  parameter = max(G(), minimum)
  lastVal = 0.
  ldir = "." # Last directory
  # First Calculation
  ndir = prefix+parameterFormatString.format(parameter)
  MKDIR(ndir)
  CopyVASPIn(ndir, useCONTCAR=False)
  os.chdir(ndir)
  
  plabels = ["", "ENCUT", "KPT Grid", "Vacuum"]
  punit =   ["", "eV",    "",         "Å³"]
  label = ["", "Energy without entropy", "E-Fermi", "Volume"]
  unit =  ["", "eV",                     "eV",      "Å³"]
  
  print("Running first calculation...")
  print(plabels[vary] + ": " + parameterFormatString.format(parameter)) # e.g. "KPT Grid: 3"
  if surface or vary == 3:
    SetISIF(True) # Make sure the calculation is constant volume
  if surface and vary == 2:
    S(parameter, gamma=True, surface=surface)
  else:
    S(parameter)
  RunCalculation()
  lastVal = R()
  print(plabels[vary]+"\t"+label[target])
  # e.g. "KPT Grid         Volume"
  print(parameterFormatString.format(parameter)+" "+punit[vary]+"\t"+targetFormatString.format(lastVal)+" "+unit[target]) # e.g. "1                 100.00000 Å³"
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
    ndir = prefix+parameterFormatString.format(parameter) 
    if checkCount > 0:
      ndir = ndir + "-{0}".format(checkCount)
     # mkdir $ndir
    MKDIR(ndir)
     # Copy files to new directory
    CopyVASPIn(ndir, origin = ldir, useCONTCAR=True)
     # cd $ndir
    os.chdir(ndir)
     # Set new parameter
    if surface or vary == 3:
      SetISIF(True) # Make sure the calculation is constant volume
    if surface and vary == 2:
      S(parameter, gamma=True, surface=surface)
    else:
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
    print(parameterFormatString.format(parameter)+" "+punit[vary]+"\t"+targetFormatString.format(lastVal)+" "+unit[target])
    # 700.0 eV      12345.12 Å³
    if abs(dV) <= tolerance and not first:
      if not conv:
        continue
      if check and not checked:
        checkCount += 1
        checked = True
      else:
        print(label[target],"converged to",targetFormatString.format(Val),unit[target],"with delta:", dV)
        return
    else:
      checked = False
      first = False
      parameter += dP
      if parameter > maximum:
        print(plabels[vary],"passed maximum value.")
        print(label[target],"failed to converge before then.")
        print("Last delta: "+targetFormatString.format(dV)+" {}".format(unit[target]))
        return
      if checkCount > 0:
        checkCount = 0


def FileExists(filename, debug=False):
  """Tries to open a file and confirm it isn't empty.
  
  Returns false (and warns with debug=True) if it is empty or inaccessible
  """
  try:
    f = open(filename, 'r')
    f.seek(0,2)
    if f.tell() == 0:
      if debug:
        print("Empty file:",filename)
      f.close()
      return False
  except IOError as e:
    print("ERROR: Can't find",filename)
    f.close()
    return False
  except:
    print("ERROR: Error reading",filename)
    f.close()
    return False
  f.close()
  return True


def ValidateDirectory(debug=False):
  """Confirm necessary VASP inputs are present."""
  return FileExists("INCAR", debug)
    and FileExists("POSCAR", debug)
    and FileExists("KPOINTS", debug)
    and FileExists("POTCAR", debug)


def GetPrintPrecision(tolerance):
  """Returns int, number of decimals to display in a format string.
  
  tolerance - float target difference between 2 values
  {0:10.{precision}f}.format(10.0, precision=printPrecision(tolerance))
  """
  if tolerance >= 0:
    return 0
  # Log10 will be negative because tolerance with be below 0
  width = math.log10(tolerance)
  return math.ceil(-width)


def Run():
  """Run a series of VASP calculations to converge a calculated property while varying calculation parameters.
  
  This function just processes command line input to feed to Converge()
  Ground State Energy, Fermi level, geometry and volume can all change based on calculation parameters.
  You can vary parameters such as ENCUT, KPOINTS grid, PREC
  This program doesn't set or change PREC
  """
  parser = argparse.ArgumentParser(
    description=
    """Runs a series of VASP calculations to converge calculated properties.
    
    In a folder with valid VASP INCAR, KPOINTS, POTCAR, and POSCAR files.
    Creates new sub directories for the calculated series.
    Ground State Energy, Fermi level, geometry, and volume can all change based on calculation parameters.  You can vary parameters such as ENCUT, KPOINTS grid, PREC, etc.
    This program only modifies ENCUT, c lattice constant (vacuum space), or KPOINTS grid.
    TODO Add varying a,b,c or a=b=c lattice constants
    IMPORTANT: An alias 'vasp' is expected to be present in your .cshrc files set to start a VASP calculation. i.e. alias vasp 'mpirun -np 8 vasp_std'""")
  parser.add_argument("calculationDirectory", nargs='?',
  default='.', help="Directory containing VASP input files, defaults to current directory.")
  parser.add_argument("-d", "--debug", action="store_true",
  help="Print verbose debug information")
  parser.add_argument("-EwE", "-ewe", action="store_true",
  dest="EwE",
  help="Converge calculations with respect to ground state Energy without entropy value.")
  parser.add_argument("-Ef", "-ef", "-fermi", action="store_true",
  dest="Efermi",
  help="Converge calculations with respect to ground state Fermi level (E-fermi + alpha+bet.")
  parser.add_argument("-Vol", "-volume", action="store_true",
  dest="Vol",
  help="Converge calculations with respect to ground state lattice volume axb.c (Triple scalar product).")
  parser.add_argument("-ENCUT", "-encut", action="store_true",
  dest="ENCUT",
  help="Converge calculations with varying ENCUT (Changes number of planewaves / ~Wavefunction size).")
  parser.add_argument("-KPT", "-kpoints", "-kpt", action="store_true",
  dest="KPT",
  help="Converge calculations with varying Gamma centered KPOINTS grid (1x1x1, 2x2x2, ...).")
  parser.add_argument("-Vacuum", "-vacuum", action="store_true",
  dest="Vacuum",
  help="Converge calculations with varying c (Vacuum space).")
  parser.add_argument("-min", nargs='?', default=0, type=float,
  dest="minimum",
  help="Lowest starting value, defaults ENCUT: ENMIN from POTCAR, KPT: 1, Vacuum: length(c).")
  parser.add_argument("-max", nargs='?', default=0, type=float,
  dest="maximum",
  help="Highest value to try to converge to, default is unlimited.  Until you run out of memory (SEG FAULT), can't wait any longer (Ctrl-C), or the process is killed by technology failure (power outage, terminal exit HUP, etc.)")
  parser.add_argument("-dp", "-delta", nargs='?', default=0, type=float,
  dest="dp",
  help="Amount to vary parameter, defaults ENCUT: 25 eV, KPT: 2 (Odd grid 1, 3, 5, ...), Vacuum: 10% initial length(c)")
  parser.add_argument("-t", "-tolerance", nargs='?', default=0.001, type=float,
  dest="tolerance",
  help="Tolerance of convergence, default 0.001, converges when energy changes by <= 1 meV or volume changes by <= 0.001 cubic Angstrom, or Fermi level changes by <= 1 meV.")
  parser.add_argument("-c", "-check", action="store_true",
  dest="check",
  help="Confirm convergence by repeating the last calculation (Make sure it only takes one geometry step starting from the converged CONTCAR), cheaper than -C2")
  parser.add_argument("-c2", "-geometry", action="store_true",
  dest="convRepeat",
  help="Repeats all calculations until a single geometry step is taken.\n\
  This is important for volume relaxation (ISIF=3,?) because new Planewave sets will be used as the volume changes altering convergence results.\n\
  e.g. If 1x1x1 changes the volume 10 - 20 a 2nd run might continue it 20-22\n\
       Then 2x2x2 changes volume 10 - 18, a 2nd run might continue to 18-20\n\
       Further down it may falsely converge the parameters early or late resulting in lower accuracy/precision of results or extra unnecessary calculation time.  Though you may have already wasted a lot running this program...")
  parser.add_argument("-s", "-surface", "-surf", action="store_true",
  dest="surface",
  help="Tip the program off if this is intended to be a surface, not a bulk calculation.  This will ensure ISIF does not optimize volume and KPOINTS grid stays X X 1 to avoid sampling along the vacuum.")
  args = parser.parse_args()
  # Check arguments
  debug = args.debug
  directory = args.calculationDirectory
  convRepeat = args.convRepeat
  check = args.check
  tolerance = args.tolerance
  dP = args.dp
  surface = args.surface
  target = 0
  vary = 0
  if args.EwE:
    target = 1
  elif args.Efermi:
    target = 2
  elif args.Vol:
    target = 3
  if args.ENCUT:
    vary = 1
  elif args.KPT:
    vary = 2
  elif args.Vacuum:
    vary = 3
  if target == 0 or vary == 0:
    print("Invalid convergence selection")
    parser.print_help()
    return
  if target == 3 and vary == 3:
    print("Invalid convergence selection")
    print("""Changing vacuum length should be done with constant volume... I don't know how to selectively optimize the a and b lattice vectors while keeping vector c constant.
    That might be a job for a different program running single point calculations while varying a and b and monitoring the cell stress output to minimize only those 2 vectors.""")
    parser.print_help()
    return
  if debug and sum([args.ENCUT, args.KPT, args.Vacuum]) > 1:
    print("WARNING: This program only modifies 1 parameter at a time.")
    print("It is set up so Vacuum > KPT > ENCUT, so ENCUT or KPT will be ignored if Vacuum is also selected.")
  if debug and sum([args.EwE, args.Efermi, args.Vol]) > 1:
    print("WARNING: This program currently only tracks on convergence value. TODO wait for all selected parameters to converge")
    print("Priority is Volume > E-fermi > Energy without entropy.  So it will favor just focusing on volume or E-fermi changes if multiple are selected.")
  if not os.exists(directory):
    print("Can't access directory {}".format(directory))
    return
  # Now run program
  if directory != '.':
    os.chdir(directory)
  if directory != os.getcwd():
    print("Failed to change directory?, program error?")
    print("{} != {}".format(directory, os.getcwd()))
    return
  if ValidateDirectory(debug):
    Converge(target, vary, dP=dP, minimum=args.minimum, maximum=args.maximum, tolerance=tolerance, check=check, convRepeat=convRepeat, surface=surface, debug=debug)
  print("Good Day!")

if __name__ == "__main__":
  Run()
