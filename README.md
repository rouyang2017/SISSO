Version SISSO.3.3, July, 2023.   
This code is licensed under the Apache License, Version 2.0  

If you are using this code, please cite:   
R. Ouyang, S. Curtarolo, E. Ahmetcik, M. Scheffler, and L. M. Ghiringhelli, Phys. Rev. Mater. 2, 083802 (2018).  

Features   
--------
- Regression & Classification    
  Ref.: [R. Ouyang et al., Phys. Rev. Mater. 2, 083802 (2018)]   
- Multi-Task Learning (MT-SISSO)    
  Ref.: [R. Ouyang et al., J. Phys.: Mater. 2, 024002 (2019)]   
- Variables Selection assisted Symbolic Regression (VS-SISSO, see the VarSelect.py in 'utilities')   
  Ref.: [Z. Guo et al., J. Chem. Theory Comput. 18, 4945 (2022)]
- Sign-Constrained Multi-Task Learning (SCMT-SISSO)    
  Ref.: [J. Wang et al., J. Am. Chem. Soc. 145, 11457 (2023)]  

(Please refer to the Refs. and the SISSO_guide.pdf for more details in using the code)  


Installation
------------
A Fortran mpi compiler is required to compile the SISSO parallel program. Below are two options for compiling the program using an IntelMPI compiler (other compilers may work as well). In the folder 'src', do:    
(1)  mpiifort -fp-model precise var_global.f90 libsisso.f90 DI.f90 FC.f90 SISSO.f90 -o ~/bin/SISSO    
or (2)  mpiifort -O2 var_global.f90 libsisso.f90 DI.f90 FC.f90 SISSO.f90 -o ~/bin/SISSO
  
Note:
- option (1) enables better accuracy and run-to-run reproducibility of floating-point calculations; (2) is ~ 2X faster 
  than (1) but tiny run-to-run variations may happen between processors of different types, e.g. Intel and AMD.   
- if 'mpi' related errors present during the compilation, try opening the file 'var_global.f90' and replace
  the line "use mpi" with "include 'mpif.h'". However, " use mpi " is strongly encouraged 
  (see https://www.mpi-forum.org/docs/mpi-3.1/mpi31-report/node411.htm).

Modules of the program:  
- var_global.f90     ! declaring the global variables
- libsisso.f90       ! subroutines and functions for mathematical operations
- DI.f90             ! model sparsification (descriptor identification)
- FC.f90             ! feature construction
- SISSO.f90          ! the main program


Running SISSO
-------------
Input Files: SISSO.in and train.dat, whose templates can be found in the folder input_templates.  
Note that the input templates and the tools in the folder utilities may be modified accordingly when a new version of the code is released. Thus, users are recommended to always use the updated files, in particular the SISSO.in.

Command-line usage:   
 SISSO > log  ! You may need to remove resource limit first by running the command 'ulimit -s unlimited'  
Running on computer clusters, for example, using this command in your submission script:   
 mpirun -np number_of_cores SISSO >log    

Primary Output Files: 
- File "SISSO.out": overall information from feature construction to model building
- Folder "models": the top ranked descriptors/models
- Folder "SIS_subspaces": SIS-selected subspaces (feature data and expressions)
- Folder "desc_dat": the data for the best descriptors/models
- File "convexnd_hull": the vertices of the nD convex hulls in classification
- File "VS_results": the results from the VS-SISSO run.


User guide
----------
More details on using this code can be found in the SISSO_guide.pdf


About
------
Created and maintained by Runhai Ouyang. Please feel free to open issues in the Github or contact Ouyang  
(rouyang@shu.edu.cn) in case of any problems/comments/suggestions in using the code. 


Other SISSO-related codes
-------------------------
SISSO++: https://gitlab.com/sissopp_developers/sissopp
MATLAB: https://github.com/NREL/SISSORegressor_MATLAB  
Python interface: https://github.com/Matgenix/pysisso  


