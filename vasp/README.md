# computationalChemistryTools.vasp

Most highly recommended
  vasp_rms
  autoConverge
  nearest_neighbor_data
  supercell & reflect
  vasp_runbands & vaspAssembleBands
  vasp_basic_analysis
  

Some of these may need minor cleaning up... (from 2015-2017)
They mostly involved a quick copy and paste from a previous program.
So there is left over cruft that I'm cleaning out before uploading.

autoConverge.py - Run calculations repeatedly with increased ENCUT or MP Grid
                  to attempt to converge total energy, volume, or fermi level
                  Currently counts on linux, csh only needs a quick edit to use bash instead

centerAndTrim.py - Cut down the atoms in a cell
                   Was for post processing by another program
                   Not practical for maintaining proper lattice
                   REQUIRES NUMPY was experimenting with it
                     Or update with code from utility.vec3

nearest_neighbor_data.py - Gives a list of nearest neighbors for each atom
                           Helpful when looking for structural changes quantitatively
                           Could be improved slightly with atom naming included with index in the output

print_bands.py - Writes band energies (EIGENVAL) to a text file for easier reading

read_eigen.py - Prints details about the band structure to the terminal
                Band gap, band edges, etc. from EIGENVAL and OUTCAR
                For quick analysis of many band structure calculations

read_jmol-anime.py - Prints vibration energies from output of Phonopy
                     anime.xyz_jmol output from Phonopy of crystal vibrations
                     Easier than parsing the OUTCAR, and usually want to see the vibrations anyway
                     Jmol used to view the output
                     REQUIRES compChemTools.utility.vec3

reflect.py - Moves atom coordinates around lattice edges
             Could probably use a better name
             e.g. if an atom is at 0,0,.99 it will be moved to 0,0,-0.01
             Useful with the supercell.py script to connect cells
             without having atoms overlap
             
rotate.py - Rotates atoms and lattice about a given axis and angle
            WARNING: Newest, least tested, quickly thrown together
                     Not personally used yet
            TODO Add option to only rotate the atoms
                 Test

supercell.py - Create a supercell from a given input cell
               Option to surround a cell with a different cell
               Needs the same lattice
               REQUIRES compChemTools.utility.vec3

vasp_basic_analysis.py - Reads a directory and prints general information about
                         the VASP calculation...
                         iterations, final energy, run time, cell volume, atoms

vasp_directory_peek.py - Reads through a hierarchy of directories, looks for
                         vasp calculations, and makes a very bare html page to
                         display some data and results of each calculation.
                         Index of all calculations on the left of the page

vasp_ReadNEB.py - Reads through a directory of a VASP NEB calculation
                  Creates a frames directory and .vmd file to load them
                  all into VMD to easily watch the trajectory
                  WARNING: May need some testing for complicated deep directories
                  
vasp_rms.py - Reports the RMSD of atom positions between two crystals CONTCAR/POSCAR
              Uses cartesian coordinates, so lattice mismatch should be watched for
                  
vasp_runbands.py - Creates directories and run script to run band structure
                   calculations along your defined band lines.
                   Then just run the output script assuming your aliases are
                   setup and wait for your results

vaspAssembleBands.py - Creates a gnuplot script from band structure calculation output
                       Run gnuplot for a quick pretty plot
                       WARNING: gnuplot labeling settings specific to Silicon
                       Intended to be used with vasp_runbands.py
                        
vaspForcesRunningAverage.py - Glorified shell script
                              Reads an OUTCAR and displays running average
                              of force during optimization of the system
                              May not work without vasp_tst to add the FORCES: output

vaspSPGSymmetry.py - Reads CONTCAR/POSCAR and reports crystal symmetry
                     REQUIRES spglib for python