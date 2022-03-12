Version SISSO.3.1, February, 2022.   
This code is licensed under the Apache License, Version 2.0  

If you are using this code, please cite:   
R. Ouyang, S. Curtarolo, E. Ahmetcik, M. Scheffler, and L. M. Ghiringhelli, Phys. Rev. Mater. 2, 083802 (2018).  

Features:   
- Regression & Classification    
  Ref.: [R. Ouyang et al., Phys. Rev. Mater. 2, 083802 (2018)]   
- Multi-Task Learning    
  Ref.: [R. Ouyang et al., J. Phys.: Mater. 2, 024002 (2019)]   
- Variables Selection for SISSO (see the program VarSelect.py in 'utilities')   
  Ref.: [Z. Guo et al., submitted]   
- Sign-Constrained Multi-Task Learning    
  Ref.: [J. Wang et al., in preparation]  

(Please refer to the Refs. and the SISSO.in template for more details in using these features)  


Installation
-------------
A Fortran mpi compiler is required to compile the SISSO parallel program. Below are two options for compiling the program using an IntelMPI compiler (other compilers may work as well). In the folder 'src', do:    
(1)  mpiifort -fp-model precise var_global.f90 libsisso.f90 DI.f90 FC.f90 SISSO.f90 -o ~/bin/your_code_name  
(2)  mpiifort -O2 var_global.f90 libsisso.f90 DI.f90 FC.f90 SISSO.f90 -o ~/bin/your_code_name  
  
Note:
- option (1) enables better accuracy and run-to-run reproducibility of floating-point calculations; (2) is ~ 2X faster 
  than (1) but tiny run-to-run variations may happen between processors of different types, e.g. Intel and AMD.   
- if 'mpi' related errors present during the compilation, try opening the file 'var_global.f90' and replace
  the line "use mpi" with "include 'mpif.h'". However, " use mpi " is strongly encouraged 
  (see https://www.mpi-forum.org/docs/mpi-3.1/mpi31-report/node411.htm).

Modules in the code:  
- var_global.f90     ! declaring the global variables
- libsisso.f90       ! needed subroutines and functions 
- DI.f90             ! model sparsification (descriptor identification)
- FC.f90             ! feature construction
- SISSO.f90


Running SISSO
-------------
Input Files: "SISSO.in" and "train.dat", whose templates can be found in 'input_templates'.  
Command-line usage:   
 SISSO > log  ! You may need to remove resource limit first by running the command 'ulimit -s unlimited'  
Running on computer clusters, for example, using this command in your submission script:   
 mpirun -np number_of_cores SISSO >log    

Output: 
- File "SISSO.out": overall information from feature construction to model building
- Folder "models": the top ranked descriptors/models
- Folder "feature_space": SIS-selected subspaces (feature data and formulas)
- Folder "desc_dat": the data for the best descriptor/model
- Folder "residual": residual data generated at each iteration
- File "convexnd_hull": the vertices of the nD convex hulls in classification


About
------
Created and maintained by Runhai Ouyang. Please feel free to open issues in the Github or contact Ouyang  
(rouyang@shu.edu.cn) in case of any problems/comments/suggestions in using the code. 


Other SISSO-related codes
-------------------------
MATLAB: https://github.com/NREL/SISSORegressor_MATLAB  
Python interface: https://github.com/Matgenix/pysisso  


