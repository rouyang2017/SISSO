! Copyright 2016-2018 The NOMAD Developers Group
!
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
real*8    stime_FCDI,etime_FCDI,stime_FC,etime_FC,stime_DI,etime_DI
integer   nsf,nvf,ntask,vfsize,ndimtype,maxcomplexity,rung,maxrung,lname,maxcomb,fs_size_DI,fs_size_L0,&
          L1_max_iter,L1_nlambda,L1_dens,desc_dim,nm_output,CV_fold,CV_repeat,iFCDI,npf_must,funit,task_weighting
parameter (lname=150,maxrung=20,maxcomb=10)
character vf2sf*10,opset(maxrung)*200,ptype*10,method*10,metric*10,calc*4
real*8    maxfval_lb,maxfval_ub,width,L1_tole,L1_minrmse,L1_elastic,PI
integer*8 subs_sis
logical   L1_warm_start,L1_weighted,restart
integer,allocatable:: nsample(:),ngroup(:,:)
real*8,allocatable::  prop_y(:),psfeat(:,:),res(:),pfdim(:,:),pvfeat(:,:,:)
character(len=lname),allocatable:: pfname(:)
parameter (PI=3.14159265358979d0)
integer   mpierr,mpirank,mpisize,status(MPI_STATUS_SIZE)
end module
