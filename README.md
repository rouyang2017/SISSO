Version SISSO.3.5, August, 2024.   
This code is licensed under the Apache License, Version 2.0  

If you are using this code, please cite:   
R. Ouyang, S. Curtarolo, E. Ahmetcik, M. Scheffler, and L. M. Ghiringhelli, Phys. Rev. Mater. 2, 083802 (2018).  

Key update in v3.5: Implementation of representing features in memory by S-expression tree. 
In all previous versions, features are stored in memory by data. Now, user can choose which scheme to use by specifying fstore =1 (features stored by data) or fstore=2 ( features stored by expression tree) in SISSO.in. 'fstore=1' is fast but high
memory demand; 'fstore=2' is low memory demand but could be several times slower than the former. Therefore, if you meet memory bottleneck because of large dataset (e.g. >5K), then set 'fstore=2'; otherwise, use 'fstore=1'. To give a feel on this, a performance comparison is shown in the SISSO_Guide.pdf at the introduction of the keyword 'fstore'.


Features   
--------
- Regression & Classification    
  Ref.: [R. Ouyang et al., Phys. Rev. Mater. 2, 083802 (2018)]   
- Multi-Task Learning (MT-SISSO)    
  Ref.: [R. Ouyang et al., J. Phys.: Mater. 2, 024002 (2019)] 
        [J. Wang et al., J. Am. Chem. Soc. 145, 11457 (2023)]  
- Variables Selection assisted Symbolic Regression (VS-SISSO, see the VarSelect.py in 'utilities')   
  Ref.: [Z. Guo et al., J. Chem. Theory Comput. 18, 4945 (2022)]

(See the Refs. and the SISSO_guide.pdf for more details in using the code)  


Installation
------------
A Fortran mpi compiler is required to compile the SISSO parallel program. Below are two options for compiling the program using an IntelMPI compiler (other compilers may work as well). In the folder 'src', do:    
(1)  mpiifort -fp-model precise *f90 -o ~/bin/SISSO    
(2)  mpiifort -O2 *f90 -o ~/bin/SISSO    
  
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
- FC.f90             ! feature construction with features stored in memory as numerical data
- FCse.f90           ! feature construction with features stored in memory as expression tree
- SISSO.f90          ! the main program


Running SISSO
-------------
Input Files: SISSO.in and train.dat, whose templates can be found in the folder input_templates.  
Note that the input templates may be different from version to version. Thus, it is always recommended to take the corresponding input templates from the package where the code is in use.

Command-line usage:   
 SISSO > log  ! You may need to remove resource limit first by running the command 'ulimit -s unlimited'  
Running on computer clusters, for example, using this command in your submission script:   
 mpirun -np number_of_cores SISSO >log    

Primary Output Files: 
- File "SISSO.out": overall information from feature construction to model building
- Folder "Models": list of the top ranked models, and the data for the top-1 model
- Folder "SIS_subspaces": SIS-selected subspaces (feature data and expressions)


User guide
----------
More details on using this code can be found in the SISSO_guide.pdf


About
------
Created and maintained by Runhai Ouyang. Please feel free to open issues in the Github or contact Ouyang  
(rouyang@shu.edu.cn) in case of any problems/comments/suggestions in using the code. 


Other SISSO-related codes
-------------------------
MATLAB: https://github.com/NREL/SISSORegressor_MATLAB  
Python interface: https://github.com/Matgenix/pysisso  
SISSO++: https://gitlab.com/sissopp_developers/sissopp


