! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.


module var_global
! global variabls
use mpi

implicit none
real*8    stime_FCDI,etime_FCDI,stime_FC,etime_FC,stime_DI,etime_DI,decorr_theta,decorr_delta,decorr_alpha
integer   nsf,nvf,ntask,vfsize,ndimtype,fcomplexity,rung,maxrung,lname,maxcomb,fs_size_DI,fs_size_L0,ptype,&
          L1_max_iter,L1_nlambda,L1_dens,desc_dim,nmodels,CV_fold,CV_repeat,iFCDI,fileunit,task_weighting,&
          nreaction,npoints,restart
parameter (lname=150,maxrung=20,maxcomb=10)
character vf2sf*10,ops(maxrung)*200,method_so*10,metric*10
real*8    fmax_min,fmax_max,bwidth,L1_tole,L1_minrmse,L1_elastic,PI
integer*8 nf_sis(10000),nsis(10000)
logical L1_warm_start,L1_weighted,fit_intercept,ffdecorr,scmt
integer,allocatable:: nsample(:),ngroup(:,:),isconvex(:,:),react_speciesID(:,:)
real*8,allocatable:: prop_y(:),psfeat(:,:),res(:),feature_units(:,:),pvfeat(:,:,:),react_coeff(:,:)
character(len=lname),allocatable:: pfname(:)
parameter (PI=3.14159265358979d0)
integer   mpierr,mpirank,mpisize,status(MPI_STATUS_SIZE)
end module
