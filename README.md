# computationalChemistryTools
Scripts that help with computational chemistry tasks

I assume a good amount of this is available elsewhere.
These were done for personal familiarity, programming experience, and learning.

Keeping them all in a folder structure
  data & utility - common scripts used in many others (keep in same relative location)
  namd - scripts specifically to work with NAMD, VMD, PSF, PBD, CHARMM files
  vasp - scripts for VASP I/O
More specifics in those folders

Linux setup for beginners
  1. Copy scripts wherever you like (e.g. ~/scripts/python or /usr/local/scripts/python)
       So it should look like
         ~/scripts/python/compChemTools/data
                                       /utility
                                       /vasp
  2. Update environmental variables
     BASH
     vi ~/.bashrc
       PATH = $PATH:~/scripts/python/compChemTools/vasp   <- so you can directly type the script name from any folder
       PYTHONPATH = $PYTHONPATH:~/scripts/python

Slowly converting everything I've done to python for easy distribution.

-Rory