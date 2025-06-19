# Digital-Appendix
Digital Appendix to the thesis: Systematic creation of complex geometries and wave functions for transition metals of groups 3-7

## 01 Complex Generation
### Database for Ligands and Metals
Ligands and Metals Database (Ligans.db and Metals.db respectively) can be reviewed via Database.py with both Ligands.db and Metals.db in the same folder.

Requires:

### Program for combining Ligands and Metals
Can be used by downloading the whole folder "Generation" and using the Program StartUp.py in this folder.

Requires:

## 02 XTB Calculations
Contains both the generated input files and the resulting output files from the geometry optimization.
## 03 GOAT
### Input Files and Output Files
Input and Output Files from the GOAT conformer search.
### Scripts
XTB_to_GOAT.py: Script which transforms the output from the XTB calculations to the GOAT input files. Requires the path to the XTB calculations folder containing both the input and output files as well as the desired path for the GOAT input files.

SpinHandling.py: Sorts complexes with both spin states according to their energies and only keeps the favorable spin state.

Extract_Data.py: Extracts the bond lengths, angles and torsions from the calculated data and saves them to three .json files respectively
 ## 04 Comparison
 Contains three folders with generated .json files from the GOAT-data, the CCDC and the OMol25 Database. The data can be compared via the modified Program Interactive.py from Lisa Matttiussi (https://github.com/chaosliza/Interactive-Graphs)

