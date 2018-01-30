#! /bin/python3

if __name__ == "__main__":
  # Collect data
  f = open('OUTCAR', 'r')
  forces = []
  for l in f:
    if "FORCES:" in l:
      forces.append([float(x) for x in l.split()[-2:]])
  f.close()
  # Early exit check to avoid errors from empty arrays
  if(len(forces) == 0):
    print("No 'FORCES:' entries found")
  else:
    print(forces[-1])
    window = 5
    ave = []
    # Average for first set of numbers
    # Store sums initially and loop to
    # divide into an average at the end
    ave.append([
      sum([x[0] for x in forces[:window]]),
      sum([x[1] for x in forces[:window]])
    ])
    avMax = ave[0][0]
    avRMS = ave[0][1]
    for i in range(window, len(forces)):
      avMax = avMax - forces[i - window][0] + forces[i][0]
      avRMS = avRMS - forces[i - window][1] + forces[i][1]
      ave.append([avMax, avRMS])
    # Print results
    print("Running Average (window size "+str(window)+")")
    print("Max, RMS")
    for x,y in ave:  # zip(*ave):
      print("{0:>8.4f}{1:>8.4f}".format(x/window,y/window))
