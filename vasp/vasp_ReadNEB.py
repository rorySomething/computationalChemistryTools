#!/bin/python3 -u

import os
import subprocess

def CopyFile(File, Dest):
  subprocess.call(["cp", File, Dest], env = os.environ)

def Run(args):
  d = "."
  if not os.path.exists("00"):
    print("No 00 directory... Not an NEB calc?")
    return
  if not os.path.exists("00/CONTCAR"):
    print("No CONTCAR in 00")
  if not os.path.exists("00/POSCAR"):
    print("No POSCAR data either")
  else:
    print("Using POSCAR")
  if not os.path.exists("frames"):
    print("Creating Frame Directory: frames")
    env = os.environ
    subprocess.call(["mkdir", "frames"], env = env)
    if not os.path.exists("frames"):
      print("Error making directory: frames")
      return
  cpto = os.path.abspath("frames")
  path = os.path.abspath(d)
  finished = False
  index = 0
  # First get frame 0
  p = "00"
  if os.path.exists(p):
    if os.path.exists(p+"/CONTCAR"):
      print("CONTCAR in",p)
      CopyFile(p+"/CONTCAR", cpto+"/POSCAR"+str(index))
      index += 1
    elif os.path.exists(p+"/POSCAR"):
      print("No CONTCAR, using POSCAR in",p)
      CopyFile(p+"/POSCAR", cpto+"/POSCAR"+str(index))
      index += 1
    else:
      print("No POSCAR or CONTCAR data...")
      return
  depth = 0
  up01 = False; up12 = False
  b01 = False; b12 = False
  second = False
  while not finished:
    if depth < 0:
      break
    # Check if 0-1
    if not up01 and not up12 and os.path.exists("0-1"):
      os.chdir("0-1")#subprocess.call(["cd","0-1"], env = os.environ)
      depth += 1
      b01 = False; b12 = False
    elif not b01 and not up12 and os.path.exists("01"):
      print("depth:",depth)
      p = "01"
      if os.path.exists(p+"/CONTCAR"):
        print("CONTCAR in",p)
        CopyFile(p+"/CONTCAR", cpto+"/POSCAR"+str(index))
        index += 1
      elif os.path.exists(p+"/POSCAR"):
        print("No CONTCAR, using POSCAR in",p)
        CopyFile(p+"/POSCAR", cpto+"/POSCAR"+str(index))
        index += 1
      else:
        print("No POSCAR or CONTCAR data...")
        return
      b01 = True
    elif not second and not up12 and os.path.exists("1-2"):
      os.chdir("1-2")#subprocess.call(["cd","1-2"], env = os.environ)
      if depth == 0 and up01:
        print("depth:",depth)
        second = True
      depth += 1
      up01 = False; up12 = False; b01 = False; b12 = True
    else:
      if depth == 0 and not os.path.exists("1-2") and up01:
        finished = True
        continue
      if "1-2" in os.path.basename(os.path.abspath(".")):
        up12 = True; up01 = False
      else:
        up01 = True; up12 = False
      if second and depth == 1: # End of the line
        finished = True
      os.chdir("..")#subprocess.call(["cd",".."], env = os.environ)
      depth -= 1
      b01 = False; b12 = False
  # Get last frame
  p = "02"
  if os.path.exists(p):
    if os.path.exists(p+"/CONTCAR"):
      print("CONTCAR in",p)
      CopyFile(p+"/CONTCAR", cpto+"/POSCAR"+str(index))
      index += 1
    elif os.path.exists(p+"/POSCAR"):
      print("No CONTCAR, using POSCAR in",p)
      CopyFile(p+"/POSCAR", cpto+"/POSCAR"+str(index))
      index += 1
    else:
      print("No POSCAR or CONTCAR data...")
      return
  with open("frames/frames.vmd", 'w') as f:
    f.write("#!/usr/local/bin/vmd\n")
    f.write("# VMD version: 1.9.1\n")
    f.write("display projection   Orthographic\n")
    f.write("mol new {0} type POSCAR first 0 last -1 filebonds 1 autobonds 0 waitfor all\n".
    format(os.path.abspath("frames/POSCAR0")))
    for i in range(1,index):
      f.write("mol addfile {0} type POSCAR first 0 last -1 filebonds 1 autobonds 0 waitfor all\n".
      format(os.path.abspath("frames/POSCAR"+str(i))))
    f.write("mol delrep 0 top\n")
    f.write("mol representation CPK 1.000000 0.300000 10.000000 10.000000\n")
    f.write("mol color Name\n")
    f.write("mol selection {all}\n")
    f.write("mol material Opaque\n")
    f.write("mol addrep top\n")
  f.close()
  print("Good Day.")

if __name__ == "__main__":
  import sys
  Run(sys.argv)
