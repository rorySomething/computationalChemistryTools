#!/bin/python3 -u

# Reference !!!
# http://www4.ncsu.edu/~franzen/public_html/CH795N/lecture/XX/tcf.html


import math

def Average(x):
  total = 0.
  for i in x:
    total += i
  return total/float(len(x))

elapsed = 0.
def ReadNAMDLog(filename, data, readVolume = False):
  try:
    f = open(filename)
  except:
    print("Can't open file:", filename)
    return
  if data == None:
    data = []
  timestep = None
  running = False
  minimize = False
  gbis = False
  slic = len("ENERGY: ")
  for l in f:
    if "TIMESTEP" in l:
      a = l.split()
      if len(a) == 3:
        timestep = float(a[2])
    elif "TCL:" in l:
      if "Running" in l:
        running = True
        minimize = False
      elif "Minimizing" in l:
        running = False
        minimize = True
    if "GBIS GENERALIZED BORN IMPLICIT SOLVENT ACTIVE" in l:
      gbis = True
    if running:
      if "ENERGY: " == l[:slic]:
        #print("Found info")
        a = l.split()
        #print("Step", a[1])
        time = int(a[1])
        total = float(a[11])
        volume = 0.0
        if not gbis:
          volume = float(a[18])
        if readVolume:
          data.append([time, volume])
        else:
          data.append([time, total])
  f.close()
  global elapsed
  if timestep:
    time = float(data[-1][0] - data[0][0])*timestep/1000.0
    elapsed += time
    print("Elapsed Run Time:",elapsed,"ps\r")
  else:
    time = data[-1][0] - data[0][0]
    elapsed += time
    print("Elapsed Run Time:",elapsed,"steps\r")
  return data

def ReadIn(filename, data):
  if data == None:
    data = []
  with open(filename, 'r') as f:
    while True:
      l = f.readline()
      if not l:
        break
      l = l.split()
      if len(l) == 2:
        data.append([int(l[0]), float(l[1])])
  f.close()
  return data

def CorrelationTimeFunction(data, out, mAverage = False, optimized = True):
  if mAverage: # Alter values to be relative to the average value +/- around 0 being the average
    total = 0.
    for i in data:
      total += i[1]
    average = total/float(len(data))
    for i in range(len(data)):
      data[i][1] -= average
  print("Solving correlation time function")
  timeAverages = []
  totalSteps = int(float(len(data))/2.0 * (len(data)-1))  # N/2 * (N-1), Guessed?
  sCount = 0 # Step counter out of total steps above to give user status update
  shortEnd = False
  endCount = 0
  for i in range(len(data)):
    print("{0: 4.1f}%".format(float(sCount)/totalSteps*100.),end='\r')
    total = 0.
    count = 0
    for j in range(len(data)):
      if j+i >= len(data):
        break
      total += (data[j][1] * data[j+i][1])  # A(0)*A(0+0), then A(0)*A(0+1), then A(0)*A(0+2)
      count += 1
      sCount += 1
    timeAverages.append(total / float(count))
    endCount += 1
    if optimized and timeAverages[-1]/timeAverages[0] < 0.001 and not shortEnd:
      # Begin early cutoff
      shortEnd = True
      endCount = 0
    if shortEnd and float(endCount)/len(data) > .05: # 5% more after the correlation cutoff is reached
      break
  ctime = False # Only report first find
  with open(out, 'w') as f:
    f.write("# Time\tC(t)\tNormalized C(t)\n")
    for t,i in enumerate(timeAverages):
      norm = i/timeAverages[0]
      f.write("{0: <8}\t{1: <8.3f}\t{2: <8.3f}\n".format(t, i, norm))
      if norm < 0.001 and not ctime:
        print("Correlation Time:", t," {x Timestep} (Normalized value less than 0.001)")
        totalSteps = len(data)
        print("Correlation Correction: {0:< 8.3e}".format(math.sqrt(2*t/totalSteps)))
        ctime = True # Found first instance
        if shortEnd:
          f.write("# Short end to correlation function evaluation.\n")
  f.close()

def Run(args):
  if len(args) < 2:
    print("Input text file list of '(Timestep) (Energy or value)'")
    print("Option: -a to shift values relative to the average value (v - average)")
    print("\t-n namd logfile as input")
    print("\t-v With -n read volume instead of energy")
    print("autocorrelation.py input.txt -o output.txt")
    print("Assumes decay time is when normalized correlation value reaches 0.001")
    return
  out = "auto-correllation-out.txt"
  i = 1
  mAverage = False
  data = None
  volume = False
  namd = False
  frame = None
  fileName = []
  while i < len(args):
    if args[i][0] == '-':
      if args[i][1] == 'o':
        out = args[i+1]
        i += 1
      elif args[i][1] == 'a':
        mAverage = True
      elif args[i][1] == 'n':
        namd = True
        i += 1
        while i < len(args):
          if args[i][0] != '-':
            fileName.append(args[i])
            i += 1
          else:
            i -= 1 # To counter end of loop i+=1
            break
      elif args[i][1] == 'v':
        volume = True
      elif args[i][1] == 'f':
        frame = int(args[i+1])
        i += 1
    else:
      fileName.append(args[i])
      data = ReadIn(fileName[-1], data)
    i += 1
  if namd:
    for i in fileName:
      data = ReadNAMDLog(i, data, volume)
  if not data:
    print("Trouble reading input", fileName)
  if frame:
    last = 0
    for i in range(0, len(data), frame):
      end = min(last + frame, len(data)-1)
      if end == last:
        break
      out1 = out[:-4] + "{0} to {1}".format(last, end) + out[-4:]
      CorrelationTimeFunction(data[last:end], out1, mAverage)
      last = end
  else:
    CorrelationTimeFunction(data, out, mAverage)
  print("Good day!")
  return

if __name__ == "__main__":
  import sys
  Run(sys.argv)
