# Bond_Analysis_Script
An Automated Python Script for Bond Analysis; it measures the bond length (Intramolecular) and Both bond length and the corresponding angle (Intermolecular).

Two scripts will be found here: 1) count_atoms11 and 2) Bond_Analysis.py

1) count_atoms11: it is a bash script that extracts Atomic Positions and converts it from fractional to Cartesian. Moreover, It extracts Cell Parameters and writes them in two different forms which are Lattice (A B C) and Cell Matrix (3x3).
   
2) Bond_Analysis.py: The script asks the user to enter the file name which it will be the start of the iteration. After that, It will show intra and inter-molecular length and angle as indices to highlight the donor and acceptor atoms along with their bond length and angle to be tracked all over the iterations. Lastly, the chosen indices will be appended in a txt file titled with the corresponding bond and atom type.

   Note: It is advised to view the desired bond using any visualization package before the first run or ignition. 
