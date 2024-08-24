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


program SISSO

use var_global
use libsisso
use FCse
use FC
use DI
use ifport
!-------------------
implicit none

integer i,j,k,l,icontinue,iostatus,line_len
character tcontinue*2,nsample_line*500,isconvex_line*500,funit_line*500,ops_line*500 !,sysdate*8,systime*10
logical fexist

call mpi_init(mpierr)
call mpi_comm_size(mpi_comm_world,mpisize,mpierr)
call mpi_comm_rank(mpi_comm_world,mpirank,mpierr)
!---

call random_seed()
call initialization  ! parameters initialization
call read_para_a     ! read from SISSO.in


fileunit=100

if(mpirank==0) then
 mytime.sFCDI=mpi_wtime()

! inquire(file='SISSO.out',exist=fexist)
! call date_and_time(date=sysdate,time=systime)
! if(fexist) iostatus=rename('SISSO.out','SISSO.out'//sysdate//systime(:6))
 
 open(9,file='SISSO.out',status='replace')
 write(9,'(a)') '****************************************************************'
 write(9,'(a)') '  Sure Independence Screening and Sparsifying Operator (SISSO)  '
 write(9,'(a)') '             Version SISSO.3.5, August, 2024.                '
 write(9,'(a/)')'****************************************************************'
end if

allocate(nsample(ntask))  ! regression: number of samples for each task
allocate(ngroup(ntask,1000))  ! classification: maximally 1000 groups inside each task
allocate(isconvex(ntask,1000)) ! restricted to be convex domain for each data group?
allocate(feature_units(nunit,nsf)) ! the units of features
call read_para_b 
npoint=sum(nsample)

allocate(target_y(npoint))   ! target property
allocate(pfdata(npoint,nsf))  ! primary scalar features                                                    
allocate(pfname(nsf)) ! primary-feature name
allocate(res(npoint))   ! residual errors
allocate(ypred(npoint))  ! fitted values at each dimension

line_len=1000+nsf*20 ! maximal line length in train.dat
call read_data  ! from train.dat 
if(mpirank==0) call output_para 


!---------------------
! FC and DI starts ...
!---------------------

if(restart==1) then
  open(fileunit,file='CONTINUE',status='old')
     read(fileunit,*) icontinue
     read(fileunit,'(a)') tcontinue
     if (tcontinue=='FC') then
        read(fileunit,*) nf_sis_avai(:icontinue-1)
     else if (tcontinue=='DI') then
        read(fileunit,*) nf_sis_avai(:icontinue)
     end if
  close(fileunit)
else
  icontinue=1
  tcontinue='FC'
  if(mpirank==0) then
     iostatus=delfilesqq('Models/data_top1/*')
     iostatus=delfilesqq('SIS_subspaces/*')
     iostatus=delfilesqq('Models/*')

     iostatus=deldirqq('Models/data_top1')
     iostatus=deldirqq('SIS_subspaces')
     iostatus=deldirqq('Models')

     iostatus=makedirqq('Models')
     iostatus=makedirqq('Models/data_top1')
     iostatus=makedirqq('SIS_subspaces')
!     iostatus=makedirqq('residual')
  end if
end if

! iterations
do iFCDI=icontinue,desc_dim
   if(mpirank==0) then
      write(*,'(/a,i3)') 'Dimension: ',iFCDI
      write(*,'(a)') '-------------------'
      write(9,'(/a,i3)') 'Dimension: ',iFCDI
      write(9,'(a)') '-------------------'
   end if

   if(iFCDI>icontinue .or. (iFCDI==icontinue .and. tcontinue=='FC') ) then
     ! run FC
     if(mpirank==0) then
         write(*,'(a)') 'Feature Construction (FC) starts ...'
         write(9,'(a)') 'Feature Construction (FC) starts ...'
     end if

     if(mpirank==0) call writeCONTINUE('FC')
     if(iFCDI==1) then
       res=target_y
     else
       if (ptype==1) res=target_y-ypred
       if (ptype==2) res=ypred
     end if
     call mpi_barrier(mpi_comm_world,mpierr)
     if(fstore==1) then
         call feature_construction   ! Storing features by data (fast, high-memory demand)
     else if(fstore==2) then        
         call feature_construction_se   ! Storing features by S-expression (low-momery demand, slower)
     end if
   end if

   ! run DI
   if(mpirank==0) then
        write(*,'(/a)') 'Descriptor Identification (DI) starts ...'
        write(9,'(/a)') 'Descriptor Identification (DI) starts ...'
   end if

   nf_DI=sum(nf_sis_avai(:iFCDI))  ! SIS-selected features for DI
   if(mpirank==0) &
       write(9,'(a,i10)') 'Total number of SIS-selected features for this DI: ',sum(nf_sis_avai(:iFCDI))
   if(trim(adjustl(method_so))=='L0') then ! number of features for L0
      nf_L0=nf_DI
   else if(trim(adjustl(method_so))=='L1L0') then  ! number of features for the L0 from L1
      nf_L0=L1para.nl1l0
   end if

   if(mpirank==0) call writeCONTINUE('DI')
   call mpi_barrier(mpi_comm_world,mpierr)
   call descriptor_identification

   call flush(9)   !Flushes Fortran unit(s)
   call flush(6)
end do

if(mpirank==0) then
 write(*,'(/a/)') 'SISSO done successfully!'
 open(1,file='CONTINUE',iostat=iostatus,status='old')
 if(iostatus==0) close(1,status='delete')
end if

deallocate(nsample)
deallocate(ngroup)
deallocate(isconvex)
deallocate(feature_units)
deallocate(target_y)
deallocate(pfdata)
deallocate(pfname)
deallocate(res)
deallocate(ypred)

call mpi_barrier(mpi_comm_world,mpierr)

if(mpirank==0) then
   mytime.eFCDI=mpi_wtime()
   write(9,'(a,f15.2)') 'Total time (second): ',mytime.eFCDI-mytime.sFCDI
   write(9,'(a/)') 'Have a nice day !    '
   close(9)
end if

call mpi_finalize(mpierr)

contains

subroutine prepare4FC
if(iFCDI==1) then
  res=target_y
else
  ! get the residual
  if (ptype==1) res=target_y-ypred
  if (ptype==2) res=ypred   
end if

end subroutine


subroutine writeCONTINUE(AA)
character AA*2
2005 format(*(i10))
   open(1234,file='CONTINUE',status='replace')
   write(1234,'(i2.2)') iFCDI
   write(1234,'(a)') AA
   write(1234,2005) nf_sis_avai(:iFCDI)
   close(1234)
end subroutine

!read in parameters from SISSO.in
subroutine read_para_a
integer i,j,k,l,ioerr
character line_short*500        

open(fileunit,file='SISSO.in',status='old')
do while(.true.)
   read(fileunit,'(a)',iostat=ioerr) line_short
   if(ioerr<0) exit
   if(index(line_short,'!')/=0) line_short(index(line_short,'!'):)=''
   i=index(line_short,'=')
   if(i>0) then
   select case (trim(adjustl(line_short(1:i-1))))
   case('restart')
   read(line_short(i+1:),*,err=1001) restart
   case('nsf')
   read(line_short(i+1:),*,err=1001) nsf
   case('ntask')
   read(line_short(i+1:),*,err=1001) ntask
   case('task_weighting')
   read(line_short(i+1:),*,err=1001) task_weighting
   case('scmt')
   read(line_short(i+1:),*,err=1001) scmt
   case('nsample') 
   read(line_short(i+1:),'(a)',err=1001) nsample_line
   case('isconvex')
   read(line_short(i+1:),'(a)',err=1001) isconvex_line
   case('desc_dim')
   read(line_short(i+1:),*,err=1001) desc_dim
   case('funit')
     read(line_short(i+1:),'(a)',err=1001) funit_line
     if(index(funit_line(k+1:),'(')>0) nunit=0
     k=0
     do while (index(funit_line(k+1:),'(')>0)  ! calculate the nunit
       k=index(funit_line(k+1:),'(')+k
       nunit=nunit+1
     end do
   case('fmax_min')
   read(line_short(i+1:),*,err=1001) fmax_min
   case('fmax_max')
   read(line_short(i+1:),*,err=1001) fmax_max
   case('fcomplexity')
   read(line_short(i+1:),*,err=1001) fcomplexity
   case('ops')
   read(line_short(i+1:),'(a)',err=1001) ops_line
   case('nf_sis')
   if(index(line_short(i+1:),',')/=0) then
       read(line_short(i+1:),*,err=1001) nf_sis(:desc_dim) ! set the subspace size for each dimension
   else 
       read(line_short(i+1:),*,err=1001) nf_sis(1) ! set the same subspace size for all dimensions
       nf_sis=nf_sis(1)  
   end if
   case('ptype')
   read(line_short(i+1:),*,err=1001) ptype
   case('fstore')
   read(line_short(i+1:),*,err=1001) fstore
   case('bwidth')
   read(line_short(i+1:),*,err=1001) bwidth
   case('nmodel')
   read(line_short(i+1:),*,err=1001) nmodel
   case('metric')
   read(line_short(i+1:),*,err=1001) metric
   case('method_so')
   read(line_short(i+1:),*,err=1001) method_so
   case('fit_intercept')
   read(line_short(i+1:),*,err=1001) fit_intercept
   case('L1para.max_iter')
   read(line_short(i+1:),*,err=1001) L1para.max_iter
   case('L1para.tole')
   read(line_short(i+1:),*,err=1001) L1para.tole
   case('L1para.nlambda')
   read(line_short(i+1:),*,err=1001) L1para.nlambda
   case('L1para.dens')
   read(line_short(i+1:),*,err=1001) L1para.dens
   case('L1para.minrmse')
   read(line_short(i+1:),*,err=1001) L1para.minrmse
   case('L1para.warm_start')
   read(line_short(i+1:),*,err=1001) L1para.warm_start
   case('L1para.elastic')
   read(line_short(i+1:),*,err=1001) L1para.elastic
   case('L1para.weighted')
   read(line_short(i+1:),*,err=1001) L1para.weighted
   end select
   end if
end do
close(fileunit)

if(fcomplexity==0) then
rung=0
elseif(fcomplexity==1) then
rung=1
elseif(fcomplexity>1 .and. fcomplexity <=3) then
rung=2
elseif(fcomplexity>3 .and. fcomplexity<=7 ) then
rung=3
elseif(fcomplexity>7 .and. fcomplexity<=15) then
rung=4
end if

return
1001 print *, 'Error: Cannot read the parameter: '//trim(line_short); stop

end subroutine

! read in the data for line nsample, funit, ops, and isconvex
subroutine read_para_b
integer*8 i,j,k,kk,l,ll

! nsample
nsample=0
ngroup=0
isconvex=1

if(ptype==1) then
  read(nsample_line,*) nsample
else 
  do ll=1,ntask
    i=index(nsample_line,'(')
    j=index(nsample_line,')')
    l=0
    do k=i,j
       if(nsample_line(k:k)==',') l=l+1
    end do
    if((i+1)>(j-1)) goto 1002
    read(nsample_line(i+1:j-1),*,err=1002) ngroup(ll,1:l+1)
    ngroup(ll,1000)=l+1   ! number of groups
    nsample(ll)=sum(ngroup(ll,1:l+1))
    nsample_line(:j)=''
  end do

  if(index(isconvex_line,'(')/=0) then
     do ll=1,ntask
       i=index(isconvex_line,'(')
       j=index(isconvex_line,')')
       l=0
       do k=i,j
          if(isconvex_line(k:k)==',') l=l+1
       end do
       if((i+1)>(j-1)) goto 1003
       read(isconvex_line(i+1:j-1),*,err=1003) isconvex(ll,1:l+1)
       isconvex(ll,1000)=ngroup(ll,1000) 
       isconvex_line(:j)=''
     end do
  end if
end if

! ops
ops=''    
if(index(ops_line,'(')/=0) then
   i=index(ops_line,'=')
   if(i>0) ops_line(:i)=''
   if(index(ops_line,',')/=0) then
       read(ops_line,*,err=10031) ops(:rung)  ! operators for each rung
   else
       read(ops_line,*,err=10031) ops(1)  ! same operators for all the rungs
       ops=ops(1)
   end if
end if

! funit
feature_units=0.d0  ! default
if(index(funit_line,'(')/=0) then 
  do ll=1,nunit
    i=index(funit_line,'(')
    j=index(funit_line,':')
    kk=index(funit_line,')')
    if((i+1)>(j-1)) goto 1004
    if((j+1)>(kk-1)) goto 1004
    if(i>0 .and. j>0) then
      read(funit_line(i+1:j-1),*,err=1004) k
      read(funit_line(j+1:kk-1),*,err=1004) l
      feature_units(ll,k:l)=1.d0
      funit_line(:kk)=''
    end if
  end do
end if

inquire(file='feature_units',exist=fexist) ! detect if the file 'feature_units' exists
if(fexist) then
  open(1,file='feature_units',status='old')
  do i=1,nsf
     read(1,*) feature_units(:,i)   ! unit of each feature represented by a vector (row)
  end do  
  close(1)
end if

return
1002 print *, 'Error: Cannot read the parameter "nsample"'; stop
1003 print *, 'Error: Cannot read the parameter "isconvex"'; stop
10031 print *, 'Error: Cannot read the parameter "ops"'; stop
1004 print *, 'Error: Cannot read the parameter "funit"'; stop

end subroutine



!read in the data from train.dat 
subroutine read_data
integer*8 i,j,k,l,mm1,mm2
character(len=str_len) string_tmp(2+nsf),reactionline(100)*10000,samplename(sum(nsample)),line_verylong*line_len
real*8 SD(ntask)

 if(mpirank==0) write(9,'(a)') 'Read in data from train.dat'

 open(fileunit,file='train.dat',status='old')
   read(fileunit,'(a)') line_verylong  ! feature names in the first line
   call sepchange(line_verylong)
   call string_split(line_verylong,string_tmp,' ')

   if(ptype==1 ) then
      pfname=string_tmp(3:2+nsf)
   else
      pfname=string_tmp(2:1+nsf)
   end if
   do i=1,sum(nsample)
       read(fileunit,'(a)') line_verylong   ! data
       call sepchange(line_verylong)
       line_verylong=adjustl(line_verylong)
       j=index(line_verylong,' ')
       samplename(i)=line_verylong(:j)
       line_verylong=line_verylong(j:)
       if(ptype==1 ) then
         read(line_verylong,*,err=1005) target_y(i),pfdata(i,:)
       else  ! ptype==2 
         read(line_verylong,*,err=1005) pfdata(i,:)   ! no y values in the train.dat file for classification
         target_y(i)=0.d0    ! 0 denote unclassified
       end if
   end do
 close(fileunit)

 1999 format(a,i3.3,a,*(f10.5))
 if(mpirank==0 .and. ptype==1) then
   do j=1,ntask
    mm1=sum(nsample(:j-1))+1
    mm2=sum(nsample(:j))
    SD(j)=sqrt(sum((target_y(mm1:mm2)-sum(target_y(mm1:mm2))/(mm2-mm1+1))**2)/(mm2-mm1+1))   ! population SD
    if(ntask==1) then
       write(9,'(a,f10.5)') 'Standard Deviation (sqrt[sum(y-y_mean)^2/n]) of the target property:',SD(j)
    else
       write(9,1999) 'Standard Deviation (sqrt[sum(y-y_mean)^2/n]) of the target property of task ',j,': ',SD(j)
    end if
   end do
 end if

return

1005 print *, 'Error while reading the file "train.dat", at: '//trim(line_verylong); stop
1007 print *, 'Error while reading the file "reaction.dat", at: '//trim(line_verylong); stop

end subroutine

! converting the separator in train.dat from Tab, comma, to space
subroutine sepchange(line)
character(len=*) line
do while (index(line,char(9))/=0)   ! separator TAB to space
 line(index(line,char(9)):index(line,char(9)))=' '
end do
do while (index(line,',')/=0)    ! separator comma to space
 line(index(line,','):index(line,','))=' '
end do
end subroutine


subroutine initialization
ptype=1               
ntask=1               
scmt=.false.          
desc_dim=2           
nsample=5             
nsf= 5                
fcomplexity=3         
nunit=1
nf_sis=50000
fstore=1
fmax_min=1e-3         
fmax_max=1e5          
restart=0            
method_so='L0'       
fit_intercept=.true.  
metric='RMSE'         
nmodel=100          
bwidth=0.001         
nf_sis_avai=0 
task_weighting=1      
!---------------------
L1para.max_iter=1e6         ! Max iteration for LASSO (given a lambda) to stop
L1para.tole=1e-6            ! Convergence criteria for LASSO to stop
L1para.dens=120             ! Density of lambda grid = number of points in [0.001*max,max]
L1para.nlambda=1e3          ! Max number of lambda points
L1para.minrmse=1e-3         ! Min RMSE for the LASSO to stop
L1para.warm_start=.true.    ! Using previous solution for the next step
L1para.nl1l0= 30            ! Number of features for L0 from L1
L1para.weighted=.false.       ! Weighted learning for L1 (provide file lasso.weight if yes)
!----------------------

end subroutine


subroutine output_para
!----------------------
! output the parameters
!----------------------
2001   format(a,*(i6))
2002   format(*(f8.2))
2003   format(a,i3.3,a,*(i5))
2004   format(*(a))
   write(9,'(a)') 'Read in data from SISSO.in'
   write(9,'(a,i3)') 'Property type:   ',ptype
   write(9,'(a,i8)') 'Number of tasks: ',ntask
   write(9,'(a,i8)') 'Descriptor dimension: ',desc_dim
   write(9,2001)  'Number of samples for the task(s): ',nsample
   write(9,'(a,i3)') 'Restarts :',restart

   if(ptype==1) then
      if(ntask>1) then
        write(9,'(a,i8)') 'Task_weighting: ',task_weighting
        write(9,'(a,l6)') 'Sign-constrained multi-task learning: ',scmt
      end if
   end if

   if(ptype==2) then
     do i=1,ntask
      write(9,2003) 'Number of samples in each category of the task ',i,': ',ngroup(i,:ngroup(i,1000))
      write(9,2003) 'Convexity of the domains of the task ',i,': ',isconvex(i,:isconvex(i,1000))
     end do
     write(9,'(a,f10.6)') 'Domain-boundary tolerance: ',bwidth
   end if

   write(9,'(a,i8)') 'Scheme for feature storing in memory: ',fstore
   write(9,'(a,i8)') 'Number of scalar features: ',nsf
   write(9,'(a,i8)')  'Tier of the feature space (max depth of expression tree): ',rung
   write(9,'(a,i8)')  'Maximal feature complexity (number of operators in a feature): ',fcomplexity
   write(9,'(a)') 'Unit of input primary feature, each represented by a row vector: '
   do i=1,nsf
     write(9,2002) feature_units(:,i)
   end do
   write(9,'(a,e15.5)') 'The feature will be discarded if the minimum of the maximal abs. value in it <',fmax_min
   write(9,'(a,e15.5)') 'The faature will be discarded if the maximum of the maximal abs. value in it > ',fmax_max
   write(9,2001) 'Size of the SIS-selected (single) subspace : ',nf_sis(:desc_dim)
   write(9,2004)  'Operators for feature construction: ',(trim(ops(j)),' ',j=1,rung)

   write(9,'(a,a)') 'Method for sparse regression:  ',method_so
   if(ptype==1) then
     write(9,'(a,l6)') 'Fitting intercept: ',fit_intercept
     write(9,'(a,a)')  'Metric for model selection: ',trim(metric)
     if(trim(adjustl(method_so))=='L1L0') then
       write(9,'(a,i10)') 'Max iterations for LASSO (with given lambda) to stop: ',L1para.max_iter
       write(9,'(a,e20.10)') 'Convergence criterion for LASSO: ',L1para.tole
       write(9,'(a,i8)') 'Number of lambda trial: ',L1para.nlambda
       write(9,'(a,i8)') 'Density of lambda points: ',L1para.dens
       write(9,'(a,e20.10)') 'Minimal RMSE for LASSO to stop: ',L1para.minrmse
       write(9,'(a,l6)') 'Warm start?  ',L1para.warm_start
       write(9,'(a,e20.10)') 'Elastic net: ',L1para.elastic
       write(9,'(a,l6)') 'Weighted LASSO (file lasso.weight required)? ',L1para.weighted
     end if
     if(trim(adjustl(method_so))=='L1L0') then
         write(9,'(a,i8)') 'Number of selected features by L1 for L0 in L1L0:', L1para.nl1l0
     end if
   end if

   write(9,'(a,i8)') 'Number of the top-ranked models to output: ',nmodel
   write(9,'(a)') '--------------------------------------------------------------------------------'

end subroutine

end program


