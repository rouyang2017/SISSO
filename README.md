SISSO: A compressed-sensing method for identifying the best low-dimensional descriptor for material properties.
This code is licensed under the Apache License, Version 2.0

If you are using this software, please cite
R. Ouyang, S. Curtarolo, E. Ahmetcik, M. Scheffler, and L. M. Ghiringhelli, Phys. Rev. Mater. 2, 083802 (2018).


Installation:
-------------
Go to folder 'src' and do: mpiifort var_global.f90 libsisso.f90 DI.f90 FC.f90 SISSO.f90 -o ~/bin/your_code_name
(Currently, only Intel compilers are acceptable)

Modules:
- var_global.f90 module for declaring global variables
- libsisso.f90   module for mathematical functions
- DI.f90         module for descriptor identification
- FC.f90         module for feature construction


Usage:
-------------
Input Files: "SISSO.in" and "train.dat"
(see folder "input_template" for templates; as a toy example, just run SISSO with the templates without any changes)

Output:
- File SISSO.out: all the information regarding parameter setting, feature space, and the best descriptors/models
- Folder models: the top ranked candidate descriptors/models
- Folder feature_space: SIS-selected subspaces (feature data and names)
- Folder desc_dat: the data for the best descriptors/models
- Folder residual: residual data generated at each iteration


About:
-------------
Created and maintained by Runhai Ouyang (Fritz Haber Institute of the Max Planck Society, Berlin, Germany).
Feel free to contact Runhai Ouyang in case of any problems/comments/suggestions in using the code.
(ouyang@fhi-berlin.mpg.de)
