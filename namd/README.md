# computationalChemistryTools.namd

Most highly recommended
  autocorrelation.py
  
Some of these may need minor cleaning up... (from 2015-2017)
They mostly involved a quick copy and paste from a previous program.
So there is left over cruft that I'm cleaning out before uploading.

autocorrelation.py - Calculates the correlation time of a series of values over a trajectory  
                     A usual MD run will have regular energy fluctuations from the thermostat
                     The correlation time is a measure of those oscillations in the system
                     This can be used to reduce the uncertainty in your average value by
                       accounting for the periodic deviations from the average (autocorrelation correction)
                     As is it can be used for any type of time varying data with proper input file format