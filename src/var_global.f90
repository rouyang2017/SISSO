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
! global variabls used by multiple modules
use mpi

implicit none

type run_time ! parameters to show the program execution time
  real*8 sFCDI,eFCDI,sFC,eFC,sDI,eDI
end type

type LASSOpara  ! parameters for LASSO
  integer max_iter,nlambda,dens,nl1l0
  real*8  tole,minrmse,elastic
  logical warm_start,weighted
end type LASSOpara

type(run_time) mytime
type(LASSOpara) L1para

integer   nsf,ntask,nunit,fcomplexity,rung,str_len,nf_DI,nf_L0,ptype,Smaxlen,&
          desc_dim,nmodel,iFCDI,fileunit,task_weighting,npoint,restart,fstore
real*8    fmax_min,fmax_max,bwidth,PI
parameter (str_len=150,PI=3.14159265d0,Smaxlen=60)
character ops(20)*200,method_so*10,metric*10
integer*8 nf_sis(10000),nf_sis_avai(10000)
logical   fit_intercept,scmt
integer,allocatable:: nsample(:),ngroup(:,:),isconvex(:,:)
real*8,allocatable:: target_y(:),pfdata(:,:),res(:),feature_units(:,:),ypred(:)
character(len=30),allocatable:: pfname(:)
integer   mpierr,mpirank,mpisize,status(MPI_STATUS_SIZE)

end module 

