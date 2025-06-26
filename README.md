# Digital-Appendix
Digital Appendix to the thesis: Systematic creation of complex geometries and wave functions for transition metals of groups 3-7

All Programs are written in Python 3.9

## 01 Complex Generation
### Database for Ligands and Metals
Ligands and Metals Database (Ligans.db and Metals.db respectively) can be reviewed via Database.py with the Database folder in the same location.

Required Packages: tkinter, sqlite3, matplotlib and numpy

### Program for combining Ligands and Metals (Startup.py)
Constructs the input data for the gfn2-xTB geometry optimization. The folder Database needs to be in the same location as the script Startup.py.

Required Packages:

## 02 XTB Calculations
### Input Files and Output Files
Contains the generated output files from the geometry optimization, as well as the Program to convert these into the needed GOAT Input files. The input files are needed as well for the conversion as they store the charge and multiplicity.
### Scripts
XTB_to_GOAT.py: Script which transforms the output from the XTB calculations to the GOAT input files. Requires the path to the XTB calculations folder containing both the input and output files as well as the desired path for the GOAT input files.

Required Packages:

## 03 GOAT
### Input Files and Output Files
Input and Output Files from the GOAT conformer search.
### Scripts
SpinHandling.py: Sorts complexes with both spin states according to their energies and only keeps the favorable spin state.

Required Packages:

Extract_Data.py: Extracts the bond lengths, angles and torsions from the calculated data and saves them to three .json files respectively

Required Packages:

## 04 Comparison
Contains three folders with generated .json files from the GOAT-data, the CCDC and the OMol25 Database. The data can be compared via the modified Program Interactive.py from Lisa Mattiussi (https://github.com/chaosliza/Interactive-Graphs)

