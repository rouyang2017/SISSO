
Version SISSO.3.0, June, 2019.
This code is licensed under the Apache License, Version 2.0

References:
R. Ouyang, S. Curtarolo, E. Ahmetcik, M. Scheffler, and L. M. Ghiringhelli, Phys. Rev. Mater. 2, 083802 (2018).  
R. Ouyang, E. Ahmetcik, C. Carbogno, M. Scheffler, and L. M. Ghiringhelli, J. Phys.: Mater. 2, 024002 (2019).

Please use corresponding input templates when switching to a new version.  
See the wiki page for the list of published works using SISSO for materials discovery.


Installation
-------------
This program is written in Fortran 90. A MPI Fortran compiler is needed for the installation.  
E.g.: go to the folder "src" and do:  
"mpiifort -O2 var_global.f90 libsisso.f90 DI.f90 FC.f90 SISSO.f90 -o ~/bin/your_code_name"  
or  
"mpigfortran -O2 var_global.f90 libsisso.f90 DI.f90 FC.f90 SISSO.f90 -o ~/bin/your_code_name"  
(The code compiled using mpiifort was found ~1.5X faster than that using mpigfortran from my tests)

Note: if 'mpi' related errors present during the compilation, please open the file 'var_global.f90' and replace
the line " use mpi " with " include 'mpif.h' ". However, " use mpi " is strongly encouraged  
(https://www.mpi-forum.org/docs/mpi-3.1/mpi31-report/node411.htm).

Modules in the code:
- var_global.f90     global variables
- libsisso.f90       library of subroutines and functions 
- DI.f90             for model sparsification
- FC.f90             for feature construction
- SISSO.f90


Running SISSO
-------------
Input Files: "SISSO.in" and "train.dat". The input templates can be found in the folder "input_template". 
To run SISSO, put in your job-submission script e.g.: 'mpirun -np xxx SISSO >log ' or 'srun SISSO >log'

Output:
- File SISSO.out: all the information regarding parameter setting, feature space, and the best descriptors/models
- Folder models: the top ranked candidate descriptors/models
- Folder feature_space: SIS-selected subspaces (feature data and names)
- Folder desc_dat: the data for the best descriptors/models
- Folder residual: residual data generated at each iteration
- Files convex2d_hull (convex3d_hull): the vertices of the 2D (3D) convex hulls for classification


About
-------------
Created and maintained by Runhai Ouyang. Please feel free to open issues in the Github or contact Ouyang  
(rouyang@shu.edu.cn) in case of any problems/comments/suggestions in using the code. 


