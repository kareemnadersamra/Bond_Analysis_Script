# Bond_Analysis_Script
An Automated Python Script for Bond Analysis; it measures the bond length (Intramolecular) and Both bond length and the corresponding angle (Intermolecular).

Four scripts will be found here: 
1) count_atoms11 and 2) Bond_Analysis.py 3)Intra_NH.py 4)Intra_OHO.py

1) count_atoms11: it is a bash script that extracts Atomic Positions and converts it from fractional to Cartesian considering Periodic Boundry Condition (PBC). Moreover, It extracts Cell Parameters and writes them in two different forms which are Lattice (A B C) and Cell Matrix (3x3).
   
2) Bond_Analysis.py: The script asks the user to enter the file name which it will be the start of the iteration. After that, It will show intra-molecular length as indices. Lastly, the chosen indices will be appended in a txt file titled with the corresponding bond and atom type.
3) Please ignore the inter-molecular data shown from Bond_Analysis as it's not accurate and use Intra NH & OHO instead.



   Notes:
            1. It is advised to view the desired bond using any visualization package before the first run or ignition.
            2. count_Atoms11 is developed based on Quantum Espresso inputs and outputs.
            3. Bond_Analysis.py starts working on files as XYZ prepared from count_atoms11.
            4. Intra_NH & OHO.py: you need to put the directory path that contains XYZ files generated from count_Atoms11, then, it will play randomly any xyz file bonds indices; however, it will track all the chosen indices for all the files regardless of the random start.
