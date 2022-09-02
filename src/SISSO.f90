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
use FC
use DI
use ifport
!-------------------

integer i,j,k,l,icontinue,iostatus,datalen
character tcontinue*2,nsample_line*500,isconvex_line*500,funit*500,sysdate*8,systime*10
logical fexist

! mpi initialization
call mpi_init(mpierr)
call mpi_comm_size(mpi_comm_world,mpisize,mpierr)
call mpi_comm_rank(mpi_comm_world,mpirank,mpierr)
!---

call random_seed()
call initialization  ! parameters initialization
call read_para_a  ! from SISSO.in
if(nreaction>0 .and. (ptype==2 .or. ntask>1)) stop "Error: Reaction ML works only with ptype=1 and ntask=1 !"
fileunit=100

if(mpirank==0) then
 stime_FCDI=mpi_wtime()

 inquire(file='SISSO.out',exist=fexist)
 call date_and_time(date=sysdate,time=systime)
 if(fexist) iostatus=rename('SISSO.out','SISSO.out'//sysdate//systime(:6))
 
 open(9,file='SISSO.out',status='replace')
 write(9,'(a)') '****************************************************************'
 write(9,'(a)') '  Sure Independence Screening and Sparsifying Operator (SISSO)  '
 write(9,'(a)') '             Version SISSO.3.2, September, 2022.                '
 write(9,'(a/)')'****************************************************************'
end if

! array allocation
allocate(nsample(ntask))  ! for regression, # samples for each task
allocate(ngroup(ntask,1000))  ! for classification, # samples in each group of a task
allocate(isconvex(ntask,1000)) 
allocate(feature_units(ndimtype,nsf+nvf))
call read_para_b
if(nreaction>0) then
  npoints=nreaction
else
  npoints=sum(nsample)
end if

allocate(prop_y(npoints))
allocate(psfeat(sum(nsample),nsf))  ! primary scalar features                                                    
allocate(pfname(nsf+nvf)) ! primary-feature name
allocate(pvfeat(sum(nsample),vfsize,nvf))  ! primary vector features             
allocate(res(npoints))

datalen=1000+max(nsf,nvf)*20
call read_data  ! from train.dat (and train_vf.dat, reaction.dat if available)
if(mpirank==0) call output_para ! output the parameters setting to SISSO.out


!---------------------
! FC and DI starts ...
!---------------------

if(restart==1) then
  open(fileunit,file='CONTINUE',status='old')
     read(fileunit,*) icontinue
     read(fileunit,'(a)') tcontinue
     if (tcontinue=='FC') then
        read(fileunit,*) nsis(:icontinue-1)
     else if (tcontinue=='DI') then
        read(fileunit,*) nsis(:icontinue)
     end if
  close(fileunit)
else
  icontinue=1
  tcontinue='FC'
  if(mpirank==0) then
     iostatus=makedirqq('desc_dat')
     iostatus=makedirqq('SIS_subspaces')
     iostatus=makedirqq('models')
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
     if(mpirank==0) call prepare4FC
     call mpi_bcast(res,npoints,mpi_double_precision,0,mpi_comm_world,mpierr)
     call mpi_barrier(mpi_comm_world,mpierr)
     call feature_construction
   end if

   ! run DI
   if(mpirank==0) then
        write(*,'(/a)') 'Descriptor Identification (DI) starts ...'
        write(9,'(/a)') 'Descriptor Identification (DI) starts ...'
   end if

   fs_size_DI=sum(nsis(:iFCDI))  ! SIS-selected features for DI
   if(mpirank==0) write(9,'(a,i10)') 'Total number of SIS-selected features from all iterations: ',sum(nsis(:iFCDI))
   if(trim(adjustl(method_so))=='L0') then  ! features for L0
      fs_size_L0=fs_size_DI
   else if(trim(adjustl(method_so))=='L1L0') then
      fs_size_L0=nl1l0
   end if

   if(mpirank==0) call writeCONTINUE('DI')
   if(mpirank==0) call prepare4DI
   call mpi_barrier(mpi_comm_world,mpierr)
   call descriptor_identification

   call flush(9)   !Flushes Fortran unit(s) currently open for output
   call flush(6)
end do

if(mpirank==0) then
 write(*,'(/a/)') 'SISSO done successfully!'
 open(1,file='CONTINUE',iostat=iostatus,status='old')
 if(iostatus==0) close(1,status='delete')
end if

3001  continue

deallocate(nsample)
deallocate(ngroup)
deallocate(isconvex)
deallocate(feature_units)
deallocate(prop_y)
deallocate(psfeat)
deallocate(pfname)
deallocate(pvfeat)
deallocate(res)
if(nreaction>0) deallocate(react_speciesID)
if(nreaction>0) deallocate(react_coeff)

call mpi_barrier(mpi_comm_world,mpierr)

if(mpirank==0) then
   etime_FCDI=mpi_wtime()
   write(9,'(/a,f15.2)') 'Total time (second): ',etime_FCDI-stime_FCDI
   write(9,'(/a/)') 'Have a nice day !    '
   close(9)
end if

call mpi_finalize(mpierr)

contains

subroutine prepare4FC
integer i,j,k,l,ntmp
real*8 rtmp,fit(npoints)
character line_short1*500,line_short2*500,line_long*10000

if(iFCDI==1) then
  res=prop_y
else
  ! get the residual
  i=0
  do k=1,ntask
     i=i+1 
     write(line_short1,'(a,i3.3,a,i3.3,a)') 'desc_',iFCDI-1,'d_p',i,'.dat'
     write(line_short2,'(a,i3.3,a,i3.3,a)') 'res_',iFCDI-1,'d_p',i,'.dat'
     open(fileunit,file='desc_dat/'//trim(adjustl(line_short1)),status='old')
!     open(1,file='residual/'//trim(adjustl(line_short2)),status='replace')
     read(fileunit,*)
     if(nreaction==0) then
         do l=1,nsample(k)
            read(fileunit,'(a)') line_long
            if (ptype==1) then
                read(line_long,*) ntmp,rtmp,fit(sum(nsample(:k-1))+l)
!                write(1,'(e20.10)') prop_y(sum(nsample(:k-1))+l)-fit(sum(nsample(:k-1))+l)
            else if (ptype==2) then
                if(index(line_long,'YES')/=0) line_long(index(line_long,'YES'):index(line_long,'YES')+2)='  1'
                if(index(line_long,'NO')/=0) line_long(index(line_long,'NO'):index(line_long,'NO')+1)=' 0'
                read(line_long,*) ntmp,res(sum(nsample(:k-1))+l)
!                write(1,'(e20.10)') res(sum(nsample(:k-1))+l)
            end if
         end do
     elseif(nreaction>0) then   ! for ptype=1 only
         do l=1,nreaction
            read(fileunit,'(a)') line_long
            read(line_long,*) ntmp,rtmp,fit(l)
!            write(1,'(e20.10)') prop_y(l)-fit(l)
         end do
     end if
     close(fileunit)
!     close(1)
  end do
  if (ptype==1) res=prop_y-fit
end if

end subroutine

subroutine prepare4DI
integer i,j,k,l
character line_short*500
real*8 feat(npoints,sum(nf_sis(:desc_dim))),rtmp

IF(nsis(iFCDI)==0) return

! Uspace - feature expressions
if(iFCDI==1) then
  open(1,file='SIS_subspaces/Uspace.expressions',status='replace')
else
  open(1,file='SIS_subspaces/Uspace.expressions',position='append',status='old')
end if
write(line_short,'(a,i3.3,a)') 'space_',iFCDI,'d.expressions'
open(2,file='SIS_subspaces/'//trim(adjustl(line_short)),status='old')
do j=1,nsis(iFCDI)
  read(2,'(a)') line_short
  write(1,'(a)') trim(adjustl(line_short))
end do
close(1)
close(2)

! Uspace - feature data
i=0
do k=1,ntask
    i=i+1
    ! previous subspaces
    if(iFCDI>1) then
        write(line_short,'(a,i3.3,a)') 'Uspace_p',i,'.dat'
        open(fileunit,file='SIS_subspaces/'//trim(adjustl(line_short)),status='old')
        if(nreaction==0) then
             if(ptype==1) then
                do l=1,nsample(k)
                  read(fileunit,*) rtmp,feat(sum(nsample(:k-1))+l,:sum(nsis(:iFCDI-1)))
                end do
             else
                do l=1,nsample(k)
                  read(fileunit,*) feat(sum(nsample(:k-1))+l,:sum(nsis(:iFCDI-1)))
                end do
             end if
        elseif(nreaction>0) then
             do l=1,nreaction
               read(fileunit,*) rtmp,feat(l,:sum(nsis(:iFCDI-1)))
             end do
        end if
        close(fileunit)
    end if

    ! this subspace
    write(line_short,'(a,i3.3,a,i3.3,a)') 'space_',iFCDI,'d_p',i,'.dat'
    open(fileunit,file='SIS_subspaces/'//trim(adjustl(line_short)),status='old')
    if(nreaction==0) then
        if(ptype==1) then
           do l=1,nsample(k)
             read(fileunit,*) rtmp,feat(sum(nsample(:k-1))+l,sum(nsis(:iFCDI-1))+1:sum(nsis(:iFCDI)))
           end do
        else
           do l=1,nsample(k)
             read(fileunit,*) feat(sum(nsample(:k-1))+l,sum(nsis(:iFCDI-1))+1:sum(nsis(:iFCDI)))
           end do
        end if
    elseif(nreaction>0) then
        do l=1,nreaction
          read(fileunit,*) rtmp,feat(l,sum(nsis(:iFCDI-1))+1:sum(nsis(:iFCDI)))
        end do
    end if
    close(fileunit)

    write(line_short,'(a,i3.3,a)') 'Uspace_p',i,'.dat'
    ! update Uspace

2000 format(*(e20.10))
     open(fileunit,file='SIS_subspaces/'//trim(adjustl(line_short)),status='replace')
     if(nreaction==0) then
       if(ptype==1) then
          do l=1,nsample(k)
           write(fileunit,2000) prop_y(sum(nsample(:k-1))+l),feat(sum(nsample(:k-1))+l,:sum(nsis(:iFCDI)))
          end do
       else
          do l=1,nsample(k)
           write(fileunit,2000) feat(sum(nsample(:k-1))+l,:sum(nsis(:iFCDI)))
          end do
       end if
     elseif(nreaction>0) then
       do l=1,nreaction
        write(fileunit,2000) prop_y(l),feat(l,:sum(nsis(:iFCDI)))
       end do
     end if
     close(fileunit)
end do
end subroutine


subroutine writeCONTINUE(AA)
character AA*2
2005 format(*(i10))
   open(1234,file='CONTINUE',status='replace')
   write(1234,'(i2.2)') iFCDI
   write(1234,'(a)') AA
   write(1234,2005) nsis(:iFCDI)
   close(1234)
end subroutine

subroutine read_para_a
integer i,j,k,l,ioerr
character line_short*500        

!read parameters from SISSO.in
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
   case('nvf')
   read(line_short(i+1:),*,err=1001) nvf
   case('vfsize')
   read(line_short(i+1:),*,err=1001) vfsize
   case('vf2sf')
   read(line_short(i+1:),*,err=1001) vf2sf
   case('desc_dim')
   read(line_short(i+1:),*,err=1001) desc_dim
   case('funit')
     read(line_short(i+1:),'(a)',err=1001) funit
     ndimtype=0
     k=0
     do while (index(funit(k+1:),'(')>0)
       k=index(funit(k+1:),'(')+k
       ndimtype=ndimtype+1
     end do
   case('fmax_min')
   read(line_short(i+1:),*,err=1001) fmax_min
   case('fmax_max')
   read(line_short(i+1:),*,err=1001) fmax_max
   case('fcomplexity')
   read(line_short(i+1:),*,err=1001) fcomplexity
   case('ffdecorr')
   read(line_short(i+1:),*,err=1001) ffdecorr
   case('decorr_theta')
   read(line_short(i+1:),*,err=1001) decorr_theta
   case('decorr_delta')
   read(line_short(i+1:),*,err=1001) decorr_delta
   case('decorr_alpha')
   read(line_short(i+1:),*,err=1001) decorr_alpha
   case('ops')
   if(index(line_short(i+1:),',')/=0) then
       read(line_short(i+1:),*,err=1001) ops(:rung)  ! multiple values
   else
       read(line_short(i+1:),*,err=1001) ops(1)  ! one value
       ops=ops(1) ! same operators for all
   end if
   case('nf_sis')
   if(index(line_short(i+1:),',')/=0) then
       read(line_short(i+1:),*,err=1001) nf_sis(:desc_dim) ! multiple values
   else 
       read(line_short(i+1:),*,err=1001) nf_sis(1) ! one value
       nf_sis=nf_sis(1)  ! same size for all
   end if
   case('ptype')
   read(line_short(i+1:),*,err=1001) ptype
   case('bwidth')
   read(line_short(i+1:),*,err=1001) bwidth
   case('nmodels')
   read(line_short(i+1:),*,err=1001) nmodels
   case('metric')
   read(line_short(i+1:),*,err=1001) metric
   case('method_so')
   read(line_short(i+1:),*,err=1001) method_so
   case('fit_intercept')
   read(line_short(i+1:),*,err=1001) fit_intercept
   case('nl1l0')
   read(line_short(i+1:),*,err=1001) nl1l0
   case('nreaction')
   read(line_short(i+1:),*,err=1001) nreaction
!--
   case('L1_max_iter')
   read(line_short(i+1:),*,err=1001) L1_max_iter
   case('L1_tole')
   read(line_short(i+1:),*,err=1001) L1_tole
   case('L1_nlambda')
   read(line_short(i+1:),*,err=1001) L1_nlambda
   case('L1_dens')
   read(line_short(i+1:),*,err=1001) L1_dens
   case('L1_minrmse')
   read(line_short(i+1:),*,err=1001) L1_minrmse
   case('L1_warm_start')
   read(line_short(i+1:),*,err=1001) L1_warm_start
   case('L1_weighted')
   read(line_short(i+1:),*,err=1001) L1_weighted
   case('L1_elastic')
   read(line_short(i+1:),*,err=1001) L1_elastic

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
1001 stop 'Error while reading file "SISSO.in" !!!'

end subroutine


subroutine read_para_b
! get the input for nsample and funit
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
    read(nsample_line(i+1:j-1),*,err=1001) ngroup(ll,1:l+1)
    ngroup(ll,1000)=l+1   ! number of groups
    nsample(ll)=sum(ngroup(ll,1:l+1))
    nsample_line(:j)=''
  end do

  do ll=1,ntask
    i=index(isconvex_line,'(')
    j=index(isconvex_line,')')
    l=0
    do k=i,j
       if(isconvex_line(k:k)==',') l=l+1
    end do
    read(isconvex_line(i+1:j-1),*,err=1001) isconvex(ll,1:l+1)
    isconvex(ll,1000)=l+1 
  end do

end if

! funit
feature_units=0.d0   ! dimensionless for default
do ll=1,ndimtype
  i=index(funit,'(')
  j=index(funit,':')
  kk=index(funit,')')
  if(i>0 .and. j>0) then
    read(funit(i+1:j-1),*,err=1001) k
    read(funit(j+1:kk-1),*,err=1001) l
    feature_units(ll,k:l)=1.d0
    funit(:kk)=''
  end if
end do

inquire(file='feature_units',exist=fexist) ! detect if the file 'feature_units' exists
if(fexist) then
  open(1,file='feature_units',status='old')
  do i=1,nsf+nvf
     read(1,*) feature_units(:,i)   ! one rwo, one feature
  end do  
  close(1)
end if

return
1001 stop 'Error while reading file "SISSO.in" !!!'

end subroutine



subroutine read_data
integer*8 i,j,k,l,ll,nentry
character(len=lname) string_tmp(2+nsf+nvf),reactionline(100)*10000,samplename(sum(nsample)),line_verylong*datalen

 if(mpirank==0) write(9,'(a)') 'Read in data from train.dat.'

 !read train.dat 
 open(fileunit,file='train.dat',status='old')
   read(fileunit,'(a)',err=1002) line_verylong  ! feature expressions
   call sepchange(line_verylong)
   call string_split(line_verylong,string_tmp,' ')

   if(ptype==1 .and. nreaction==0) then
      pfname=string_tmp(3:2+nsf)
   else
      pfname=string_tmp(2:1+nsf)
   end if
   do i=1,sum(nsample)
       read(fileunit,'(a)',err=1002) line_verylong   ! data
       call sepchange(line_verylong)
       line_verylong=adjustl(line_verylong)
       j=index(line_verylong,' ')
       samplename(i)=line_verylong(:j)
       line_verylong=line_verylong(j:)
       if(ptype==1 .and. nreaction==0) then
         read(line_verylong,*,err=1002) prop_y(i),psfeat(i,:)
       else  ! ptype==2 or nreaction>0
         read(line_verylong,*,err=1002) psfeat(i,:)   ! no y value in the train.dat file for classification
         if(nreaction==0) prop_y(i)=0.d0    ! 0 denote unclassified
       end if
   end do
 close(fileunit)


 ! read train_vf.dat
 if(nvf>0) then
   if(mpirank==0) write(9,'(a)') 'Read in data from train_vf.dat.'
   open(fileunit,file='train_vf.dat',status='old')
   read(fileunit,'(a)',err=1003) line_verylong  

   call sepchange(line_verylong)
   call string_split(line_verylong,pfname(nsf+1:),' ')  ! save vector-feature expressions

   do i=1,sum(nsample)
      read(fileunit,*,err=1003) (pvfeat(i,:,j),j=1,nvf)
   end do
   close(fileunit)
 end if


! read reaction.dat
if(nreaction>0)then
   if(mpirank==0) write(9,'(a)') 'Read in data from reaction.dat.'
   open(fileunit,file='reaction.dat',status='old')
   read(fileunit,'(a)',err=1004)    ! comment line
   prop_y=0.d0
   do i=1,nreaction
      read(fileunit,'(a)',err=1004) line_verylong 

      call sepchange(line_verylong)
      line_verylong(index(line_verylong,'!'):)=''
      reactionline=''
      call string_split(line_verylong,reactionline,' ') 
 
      nentry=0
      do j=1,100
        if(reactionline(j)/='') nentry=nentry+1
      end do
      if(i==1) then
        allocate(react_speciesID(nreaction,100))
        allocate(react_coeff(nreaction,100))
        react_speciesID=0
        react_coeff=0.d0
      end if
      do j=1,(nentry-1)/2
        do k=1,sum(nsample)
           if(trim(adjustl(reactionline(j)))==trim(adjustl(samplename(k)))) exit
        end do
        react_speciesID(i,j)=k
        read(reactionline(j+(nentry-1)/2),*,err=1004) react_coeff(i,j)
      end do
      read(reactionline(nentry),*,err=1004) prop_y(i)
   end do
   close(fileunit)
 end if

return

1002 stop 'Error while reading file "train.dat" !!!'
1003 stop 'Error while reading file "train_vf.dat" !!!'
1004 stop 'Error while reading file "reaction.dat" !!!'

end subroutine


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
task_weighting=1      
scmt=.false.          
desc_dim=1           
nsample=1             
restart=0            

nsf= 1                
ops='(+)(-)(*)(/)'    
fcomplexity=3         
funit='(1:1)'         
fmax_min=1e-3         
fmax_max=1e5          
nf_sis=20            

method_so='L0'       
nl1l0= 1       
fit_intercept=.true.  
metric='RMSE'         
nmodels=100          
isconvex=(1,1)        
bwidth=0.001         

!---------------------
nreaction=0             ! number of reactions in reaction machine learning.
ffdecorr=.false.        ! feature-feature decorrelation
decorr_theta=1.0        ! Threshold (<=1) for feature decorrelation (one of the highly correlated feature will be removed).
decorr_delta=0.0        ! Size of score-window for evaluating correlation of features with similar SIS-scores (>=0).
decorr_alpha=1.0        ! Preselecting alpha*nf_sis features to ensure the final size nf_sis after decorrelation.
vfsize= 0               ! size of each vector feature (all vectors have the same size)
vf2sf= 'sum'            ! transforming vector to scalar features: sum,norm,min,max
L1_max_iter=1e6         ! max iteration for LASSO (given a lambda) to stop
L1_tole=1e-6            ! convergence criteria for LASSO to stop
L1_dens=120             ! density of lambda grid = number of points in [0.001*max,max]
L1_nlambda=1e3          ! max number of lambda points
L1_minrmse=1e-3         ! Min RMSE for the LASSO to stop
L1_warm_start=.true.    ! using previous solution for the next step
L1_weighted=.false.     ! weighted observations? (provide file prop.weight if yes)
nsis=0
!----------------------

end subroutine


subroutine output_para
!----------------------
! output the parameters
!----------------------
2001   format(a,*(i8))
2002   format(*(f6.2))
2003   format(a,i3,a,*(i5))
2004   format(*(a))
   write(9,'(a)') 'Read in parameters from SISSO.in'
   write(9,'(a,i3)') 'Property type:   ',ptype
   write(9,'(a,i8)') 'Number of tasks: ',ntask
   if(ntask>1) then
     write(9,'(a,i8)') 'Task_weighting: ',task_weighting
     write(9,'(a,l6)') 'Sign-constrained multi-task learning: ',scmt
   end if
   write(9,'(a,i8)') 'Descriptor dimension: ',desc_dim
   write(9,2001)  'Number of samples for each task: ',nsample
   if(ptype==2) then
     do i=1,ntask
      write(9,2003) 'Number of samples in each category of task ',i,': ',ngroup(i,:ngroup(i,1000))
      write(9,2003) 'Convexity of the domains of task ',i,': ',isconvex(i,:isconvex(i,1000))
     end do
     write(9,'(a,f10.6)') 'Domain-boundary tolerance: ',bwidth
   end if
   write(9,'(a,i3)') 'Restarts :',restart
   if(nreaction>0) write(9,'(a,i8)') 'Number of training reactions: ',nreaction


   write(9,'(a,i8)') 'Number of scalar features: ',nsf
   if(nvf>0) then
     write(9,'(a,i8)') 'Number of vector features: ',nvf
     write(9,'(a,i8)') 'Size of the vector features: ',vfsize
     write(9,'(a,a)')  'Method for transforming vectors to scalars? ',trim(vf2sf)
   end if
   write(9,'(a,i8)')  'Tier of the feature space: ',rung
   write(9,'(a,i8)')  'Maximal feature complexity (number of operators in a feature): ',fcomplexity
   write(9,'(a)') 'Units of the input primary features (each represented by a vector): '
   do i=1,nsf+nvf
     write(9,2002) feature_units(:,i)
   end do
   write(9,'(a,e15.5)') 'The feature will be discarded if the minimum of the maximal abs. value in it <',fmax_min
   write(9,'(a,e15.5)') 'The faature will be discarded if the maximum of the maximal abs. value in it > ',fmax_max
   if(ffdecorr) then
      write(9,'(a,f15.5)') 'Correlation threshold for feature decorrelation: ',decorr_theta
      write(9,'(a,f15.5)') 'Size of the score-window for feature decorrelation: ',decorr_delta
      write(9,'(a,f15.5)') 'Selecting alpha*nf_sis features to ensure the final size being &
          nf_sis after decorrelation: ',decorr_alpha
   end if
   write(9,2001) 'Size of the SIS-selected (single) subspace : ',nf_sis(:desc_dim)
   write(9,2004)  'Operators for feature construction: ',(trim(ops(j)),' ',j=1,rung)


   if(ptype==1) then
     write(9,'(a,l6)') 'Fitting intercept: ',fit_intercept
     write(9,'(a,a)')  'Metric for model selection: ',trim(metric)
     if(trim(adjustl(method_so))=='L1L0') then
       write(9,'(a,i10)') 'Max iterations for LASSO (with given lambda) to stop: ',L1_max_iter
       write(9,'(a,e20.10)') 'Convergence criterion for LASSO: ',L1_tole
       write(9,'(a,i8)') 'Number of lambda trial: ',L1_nlambda
       write(9,'(a,i8)') 'Density of lambda points: ',L1_dens
       write(9,'(a,e20.10)') 'Minimal RMSE for LASSO to stop: ',L1_minrmse
       write(9,'(a,l6)') 'Weighted observations (if yes, provide file prop.weight)? ',L1_weighted
       write(9,'(a,l6)') 'Warm start?  ',L1_warm_start
       write(9,'(a,e20.10)') 'Elastic net: ',L1_elastic
     end if
   end if
   write(9,'(a,a)') 'Method for sparse regression:  ',method_so
   if(trim(adjustl(method_so))=='L1L0') then
       write(9,'(a,i8)') 'Number of selected features by L1 for L0 in L1L0:', nl1l0
   end if
   write(9,'(a,i8)') 'Number of the top-ranked models to output: ',nmodels
   write(9,'(a)') '--------------------------------------------------------------------------------'
end subroutine

end program


