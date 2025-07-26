# Digital-Appendix
Digital Appendix to the thesis: Systematic creation of complex geometries and wave functions for transition metals of groups 3-7

All Programs are written in Python 3.9

## 01 Complex Generation
### Database for Ligands and Metals
Ligands and Metals Database (Ligans.db and Metals.db respectively) can be reviewed via Database.py with the Database folder in the same location.

### Program for combining Ligands and Metals (Startup.py)
Constructs the input data for the gfn2-xTB geometry optimization. The folder Database needs to be in the same location as the script Startup.py.

## 02 XTB Calculations

### Script
XTB_to_GOAT.py: Script which transforms the output from the XTB calculations to the GOAT input files. Requires the path to the XTB calculations folder containing both the input and output files as well as the desired path for the GOAT input files.

## 03 GOAT
### Scripts
Rename.py: Groups complexes in subfolders and renames them.

sort_multiplicities.py: Sorts complexes with multiple multiplicities according to their energies and only keeps the favorable spin state.

## 04 Comparison
### Data
All .json files are in the zip folder named Data.
Each folder contains the generated .json files from the GOAT-data, the CCDC and the OMol25 Database, respectively. The data can be compared via the modified Program Interactive.py from Lisa Mattiussi (https://github.com/chaosliza/Interactive-Graphs)

### Scripts
tools.py: Contains the function for simplifying the smiles codes which was implemented in the program Interactive.py from Lisa Mattiussi (https://github.com/chaosliza/Interactive-Graphs)

xyzTomol.py: Function to transfer the xyz data to RDKit Mol objects for extracting bond lengths etc. The script was implemented in a program by Lisa Matiussi for extracting the bond lengths, angles and torsions into .json files.

### Modified_Scripts
Contains the modified scripts by Lisa Matiussi (https://github.com/chaosliza). All of these are Jupyter Notebooks utilized in an Ubuntu environment.

Interactive_Clean.ipynb: Generates the .json files with the simplified SMILES codes with the help of tools.py
Interactive.ipynb: Uses theses cleaned .json files to generate the desired plots.
own_data.ipynb: Generates the original .json files for the generated Data with the help of xyz_to_mol.py