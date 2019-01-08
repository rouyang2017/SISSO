
Version SISSO.2.4, January 8, 2019.

SISSO: A compressed-sensing method for identifying the best low-dimensional descriptor for material properties.
This code is licensed under the Apache License, Version 2.0

If you are using this software, please cite:  
R. Ouyang, S. Curtarolo, E. Ahmetcik, M. Scheffler, and L. M. Ghiringhelli, Phys. Rev. Mater. 2, 083802 (2018).   
R. Ouyang, E. Ahmetcik, C. Carbogno, M. Scheffler, and L. M. Ghiringhelli, "Simultaneous Learning of Several Materials 
Properties from Incomplete Databases with Multi-Task SISSO", submitted to J. Phys.: Mater. (2019) (Preprint: arXiv:1901.00948)


Please read the file Changelog for the update from previous to this new version. Please use the new input templates.


Installation
-------------
This program is written in Fortran 90. A MPI Fortran compiler (from Intel or GNU) is needed for the installation.
E.g.: go to the folder "src" and do:
"mpiifort -O2 var_global.f90 libsisso.f90 DI.f90 FC.f90 SISSO.f90 -o ~/bin/your_code_name"
or
"mpigfortran -O2 var_global.f90 libsisso.f90 DI.f90 FC.f90 SISSO.f90 -o ~/bin/your_code_name"
(The code compiled using mpiifort was found ~1.5X faster than that using mpigfortran from my tests)

Note: if 'mpi' related errors present during the compilation, please open the file 'var_global.f90' and comment out 
the line " use mpi " and add a new line " include 'mpif.h' " right after the line " implicit none ".
However, " use mpi " is strongly encouraged (https://www.mpi-forum.org/docs/mpi-3.1/mpi31-report/node411.htm).

The modules:
- var_global.f90     for declaring global variables
- libsisso.f90       for mathematical functions
- DI.f90             for sparsifying operators
- FC.f90             for feature construction
- SISSO.f90


Running SISSO
-------------
Input Files: "SISSO.in" and "train.dat". The input templates can be found in the folder "input_template". 
To run SISSO, just put in your submission script something like 'mpirun -np xxx SISSO' or 'srun SISSO' or else, depending on your batch system.

Output:
- File SISSO.out: all the information regarding parameter setting, feature space, and the best descriptors/models
- Folder models: the top ranked candidate descriptors/models
- Folder feature_space: SIS-selected subspaces (feature data and names)
- Folder desc_dat: the data for the best descriptors/models
- Folder residual: residual data generated at each iteration
- Files convex2d_hull (convex3d_hull): the vertices of the 2D (3D) convex hulls for classification


About
-------------
Created and maintained by Runhai Ouyang (Fritz Haber Institute of the Max Planck Society, Berlin, Germany).
Please feel free to contact Runhai Ouyang (ouyang@fhi-berlin.mpg.de) or open issues in the Github in case of any problems/comments/suggestions in using the code.


