## Package name: pdbKit
Rename the folder name to "pdbKit" once you start the development.
Examine the packages present in this folder. Create a python package with python version 3.10 and above which will integrate these into a single package with added functionality as mentioned below. Make sure that the tool functions equally well with .pdb and .mmcif three dimensional biomolecular structure files. Add the following features: 

1. - merge the packages with their functionality into a single command line python package which can be called from terminal or from other python code or used in ipython notebooks. 
2. - add options to fetch files three dimansional structure from PDB and alphafold.
3. - split chains and write to separate files (option to include only c-alpha, all-atom, c-alpha with all/specific hetatoms, all-atom with all/specific hetatoms).
4. - select chains based on user provided chain-id and provide functionality as per #3. 
5. - calculate center of mass (entire structure or specific chain) and set the center of mass to zero. Output in a new file. 
6. - include a het-pdb-navigator parser- provided in the "hetpdbnavi_parser.py" file.
7. - get domin and family information- provided in the "mmcif_domain_family.py" file.
8. - integrate the code for "hetpadnavi_parser.py" and "mmcf_domain_family.py" in the python package that
you develop. 
9. - Any resulting pdb or mmcif file produced from this package should be compatible with other macromolecular viewers like pymol, mol* etc.
10. - Provide clear installation and usage instructtions in the form of a user manual with examples of each scenario as an obsidian-compatible a markdown file.
11. Perform testing and add the testing results in the user manual as a separate chapter.

For the development purposes, use the "biopyc" conda environment.
