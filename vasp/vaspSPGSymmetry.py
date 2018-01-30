#!/bin/python3

import spglib

def SPGReadCONTCAR(filename = "CONTCAR"):
  # spg data
  # cell = (lattice, positions, numbers)
  # lattice = 3x3 matrix a, b, c in rows
  # = [ [ax,ay,az], [bx,by,bz], [cx,cy,cz] ]
  # positions = Nx3 matrix fractional
  #  = [ [x,y,z],[x,y,z],...]
  # numbers = N list of atom species
  #  = [1,1,1,1,1,2,2,2,1,1,1,1,1,2,2,2,...]
  lattice = []
  positions = []
  numbers = []
  with open(filename, 'r') as f:
    f.readline() # Title
    scale = float(f.readline())
    lattice = [
      [float(x) for x in f.readline().split()],
      [float(x) for x in f.readline().split()],
      [float(x) for x in f.readline().split()] ]
    l = f.readline().split()
    labelled = False; names = []
    counts = []
    try:
      n = int(l[0])
    except:
      labelled = True
    if labelled:
      names = l[:]
      counts = [int(x) for x in f.readline().split()]
    else:
      counts = [int(x) for x in f.readline().split()]
    i = 0
    for c in counts:
      i += 1
      for x in range(c):
        numbers.append(i)
    direct = (f.readline().lower() == "direct")
    for i in counts:
      for j in range(i):
        l = f.readline().split()
        if not direct:
          l = l[:3] # Remove rear flags for selective dynamics
        positions.append([float(x) for x in l])
    print("{0} atoms".format(len(positions)))
    return (lattice, positions, numbers)

if __name__ == "__main__":
  import sys
  args = sys.argv
  cell = []
  if len(args) > 1:
    filename = args[1]
    cell = SPGReadCONTCAR(filename)
  else:
    cell = SPGReadCONTCAR()
  print(spglib.get_spacegroup(cell))
  print(spglib.get_symmetry_dataset(cell)['choice'])
  print("Good Day!")
  
  
