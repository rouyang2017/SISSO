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


module FCse
! feature construction
!----------------------------------------------------------------------

use var_global
use libsisso
implicit none
! variables used by subroutines and functions in this module

type Sexpression
  integer list_id(Smaxlen),list_len,list_pointer(2,Smaxlen)  ! A list for each S-expression
  character list_var(Smaxlen)*30,list_op(Smaxlen)*10
end type

type feature
  type(Sexpression),allocatable:: Sexpr(:)
  character(len=str_len),allocatable:: feat_name(:),lastop(:)*10
  integer*8,allocatable:: feat_comp(:)
  real*8,allocatable:: feat_unit(:,:)
end type

type feature_selected
  integer*8 nselect
  type(Sexpression),allocatable:: Sexpr(:)
  real*8,allocatable:: feat_score(:,:)
  character(len=str_len),allocatable:: feat_name(:)
  integer*8,allocatable:: feat_comp(:)
end type

type feature_sis
  type(Sexpression),allocatable:: Sexpr(:)
  real*8,allocatable:: feat_score(:,:)
  character(len=str_len),allocatable:: feat_name(:)
end type

type(feature) gen
type(feature_selected) sel
type (feature_sis) sis

integer*8 ntot,nthis,nreject,nbasic_select,nextra_select,icomb
real*8 score_threshold
real*8,allocatable:: trainy(:),trainy_c(:)
character(len=str_len),allocatable::  reject(:)

contains


subroutine feature_construction_se
implicit none
integer   loc(1),ioerr
character line*500,phiname*5,reject_file_name*100,superline*(20*(1+sum(nf_sis(:desc_dim))))
real*8    bbb,ccc,aaa(mpisize),tag(npoint)
integer*8 i,j,k,l,ll,mm1,mm2,nf(20),mpii,mpij,mpik,total_comb,nfpcore_this(mpisize),mpin(mpisize),mpin2(mpisize)
real*8,allocatable:: fID(:),tmp_data(:,:)
integer*8,allocatable:: order(:)
character mpisync*1
logical,allocatable:: available(:)
type(feature) inp

mpisync='Y'

! Stop if the whole feature space had been selected.
IF (iFCDI>1 .and. nf_sis_avai(iFCDI-1)<nf_sis(iFCDI-1)) THEN
   if(mpirank==0) then
      write(*,'(a)') 'The whole feature-space has been selected! No more FC will be performed.' 
      write(9,'(a)') 'The whole feature-space has been selected! No more FC will be performed.'
   end if
   return
END IF

! running time by FC
if(mpirank==0) mytime.sFC=mpi_wtime()

!------------------------------------------------------
allocate(trainy(npoint))
allocate(trainy_c(npoint))
i=max(1000,nsf)
allocate(inp.Sexpr(i))   
allocate(inp.feat_name(i))
allocate(inp.lastop(i))
allocate(inp.feat_unit(nunit,i))
allocate(inp.feat_comp(i))

allocate(gen.Sexpr(i))
allocate(gen.feat_name(i))
allocate(gen.lastop(i))
allocate(gen.feat_comp(i))
allocate(gen.feat_unit(nunit,i))

nbasic_select=nf_sis(iFCDI)  ! parallel
if(mpisize==1) nbasic_select=2*nf_sis(iFCDI)  ! serial
nextra_select=min(nf_sis(iFCDI),50000)
j=nbasic_select+nextra_select

allocate(sel.Sexpr(j))
allocate(sel.feat_score(2,j))
allocate(sel.feat_name(j))
allocate(sel.feat_comp(j))

! initialization
!------------------
trainy=res
score_threshold=-1.0
inp.feat_unit(:,:nsf)=feature_units
inp.feat_name(:nsf)=pfname(:nsf)
inp.lastop=''     ! the last operation of the generated feature
inp.feat_comp=0

do i=1,nsf    
  inp.Sexpr(i).list_id(1)=1   ! first entry of the list of expression i
  inp.Sexpr(i).list_len=1   ! current list length
  inp.Sexpr(i).list_op(1)='var'  ! is variable, no operation
  inp.Sexpr(i).list_pointer(:,1)=0   ! the 2 operands specified by the pointer. No operand for variable
  inp.Sexpr(i).list_var(1)=trim(adjustl(pfname(i)))  ! the variable name
end do


! The 'tag' to be used by 'fID' during redundant check
do j=1,ntask
do i=1,nsample(j)    
  tag(sum(nsample(:j-1))+i)=1.0d0+0.001d0*i
end do
end do

! avoid duplication in the selected space
if(iFCDI>1) reject_file_name='SIS_subspaces/Uspace.expressions'

nreject=0
if(iFCDI >1 .and. mpirank==0) then
   open(fileunit,file=trim(adjustl(reject_file_name)),status='old')

   ! find the number of expressions in the file
   do while(.true.)
      read(fileunit,*,iostat=ioerr)
      if(ioerr<0) exit
      nreject=nreject+1
   end do

   if(nreject>0) allocate(reject(nreject))
   rewind(fileunit)
   
   do j=1,nreject
      read(fileunit,'(a)') line
      call string_split(line,reject(j:j),' ')
      reject(j)=adjustl(reject(j))
   end do
   close(fileunit)
   !ordering
   do i=1,nreject-1
       loc(1)=i
       do j=i+1,nreject
         if(reject(j)>reject(loc(1))) loc(1)=j
       end do
       if(loc(1)>i) then
         line=trim(reject(i))
         reject(i)=reject(loc(1))
         reject(loc(1))=trim(line)
       end if
   end do
end if

call mpi_bcast(nreject,1,mpi_integer8,0,mpi_comm_world,mpierr)
if(mpirank/=0 .and. nreject>0) allocate(reject(nreject))
if(nreject>0) call mpi_bcast(reject,nreject*str_len,mpi_character,0,mpi_comm_world,mpierr)
!---------------------------------
! Population Standard Deviation
!---------------------------------

if(ptype==1) then
   do j=1,ntask
    mm1=sum(nsample(:j-1))+1
    mm2=sum(nsample(:j))
    trainy_c(mm1:mm2)= trainy(mm1:mm2)-sum(trainy(mm1:mm2))/(mm2-mm1+1)   ! centered
   end do
end if

!---------------------------------------------------
! feature construction start ...
! Phi0 is the initial space with primary features
!---------------------------------------------------
sel.nselect=0 ! number of selected features
nf=0      ! number of generated features from each combination
i=nsf ! total number of primary features
j=0

! no combination, just primary features
total_comb=0
call combine_se(inp,1,i,0,0,'NO',j,total_comb) 
! ntot is forced to equal total number of pf as all are already stored in inp.feat, yet sel.nselect is not
! necessarily to be ntot if some pf is not good (e.g. constant).
! The purpose here is store good primary features in sel.feat

ntot=nsf    ! total number of features
nthis=ntot  ! number of features generated from this combine_se()
write(phiname,'(a,i2.2)') 'Phi',0
if(mpirank==0)  call writeout_se(phiname,sel.nselect,ntot)  

if(rung==0) then
  deallocate(gen.Sexpr)
  deallocate(gen.feat_name)
  deallocate(gen.lastop)
  deallocate(gen.feat_comp)
  deallocate(gen.feat_unit)
end if

!-------------------------
!creating PHI1, PHI2, ...
!-------------------------
do icomb=1,rung
    if(mpirank==0) write(*,'(/a,i2.2,a)') 'Generating Phi',icomb,' ...'

    call mpi_barrier(mpi_comm_world,mpierr)

    if(icomb==rung) then
       deallocate(gen.Sexpr)
       deallocate(gen.feat_name)
       deallocate(gen.lastop)
       deallocate(gen.feat_comp)
       deallocate(gen.feat_unit)
    end if

    ! allocating equivalent workfload for each core
    ! Total jobs: T=N*M + N(N-1)/2, where N=nthis, M = ntot-nthis
    ! Divide N into n parts: X1, X2, ..., Xi, ..., where n=mpisize
    ! X1=X1(X1-1)/2+X1(N-X1+M) = T/n --> X1
    ! X2=X2(X2-1)/2+X2(N-X1-X2+M) = T/n --> X2
    ! Xi=Xi(Xi-1)/2+Xi(N-X1-X2-...-Xi+M) = T/n -->Xi
    ! Let bbb=N+M-1/2-sum(X1+X2+...+X{i-1}),ccc=-(N*M+N(N-1)/2)/n
    ! Xi= bbb-sqrt(bbb^2+2ccc)
    ! mpin: store the Xi, mpin2: starting position of the Xi part in nthis

      mpin=0
      mpin2=0
      ccc=-(nthis*(ntot-nthis)+nthis*(nthis-1)/2.d0)/float(mpisize)
      do mpii=1,mpisize
         bbb=ntot-0.5-sum(aaa(:mpii-1))
         if((bbb**2+2*ccc)>0) then
            aaa(mpii)=bbb-sqrt(bbb**2+2*ccc)
            mpin(mpii)=ceiling(aaa(mpii))
            if (sum(mpin(:mpii))>nthis) then
                  mpin(mpii)=nthis-sum(mpin(:mpii-1))
                  exit
            end if
         end if
      end do

      if(sum(mpin)<nthis) then  ! it can happen when ntot==nthis
         mpin(mpisize)=mpin(mpisize)+nthis-sum(mpin)
      end if

      do i=1,mpisize
        if(mpin(i)>0) mpin2(i)=sum(mpin(1:i-1))+1
      end do

     ! check if array size need to be increased
      i=ntot-nthis+(nthis-mpin2(mpirank+1)+1)   ! total number of features to be stored in this core
      if(ubound(inp.Sexpr,1)< i)  call addm_inp_se(i-ubound(inp.Sexpr,1),inp)

     ! broadcast ntot-nthis
       if(ntot>nthis) then
        i=ntot-nthis
        do ll=1,i
          call mpi_bcast(inp.Sexpr(ll).list_id(:Smaxlen),Smaxlen,mpi_integer,0,mpi_comm_world,mpierr)
          call mpi_bcast(inp.Sexpr(ll).list_len,1,mpi_integer,0,mpi_comm_world,mpierr)
          call mpi_bcast(inp.Sexpr(ll).list_pointer(:2,:Smaxlen),Smaxlen*2,mpi_integer,0,mpi_comm_world,mpierr)
          call mpi_bcast(inp.Sexpr(ll).list_var(:Smaxlen),Smaxlen*30,mpi_character,0,mpi_comm_world,mpierr)
          call mpi_bcast(inp.Sexpr(ll).list_op(:Smaxlen),Smaxlen*10,mpi_character,0,mpi_comm_world,mpierr)
        end do
        call mpi_bcast(inp.feat_name(:i),i*str_len,mpi_character,0,mpi_comm_world,mpierr)
        call mpi_bcast(inp.lastop(:i),i*10,mpi_character,0,mpi_comm_world,mpierr)
        call mpi_bcast(inp.feat_unit(:,:i),i*nunit,mpi_double_precision,0,mpi_comm_world,mpierr)
        call mpi_bcast(inp.feat_comp(:i),i,mpi_integer8,0,mpi_comm_world,mpierr)
       end if

      if(mpirank==0) then  ! send nthis 
        do i=2,mpisize
          j=ntot-nthis+mpin2(i)  ! starting position
          k=nthis-mpin2(i)+1   ! size
          if(mpin(i)>0) then
          do ll=j,ntot
           call mpi_send(inp.Sexpr(ll).list_id(:Smaxlen),Smaxlen,mpi_integer,i-1,5,mpi_comm_world,status,mpierr)
           call mpi_send(inp.Sexpr(ll).list_len,1,mpi_integer,i-1,6,mpi_comm_world,status,mpierr)
           call mpi_send(inp.Sexpr(ll).list_pointer(:2,:Smaxlen),Smaxlen*2,mpi_integer,i-1,7,mpi_comm_world,status,mpierr)
           call mpi_send(inp.Sexpr(ll).list_var(:Smaxlen),Smaxlen*30,mpi_character,i-1,8,mpi_comm_world,status,mpierr)
           call mpi_send(inp.Sexpr(ll).list_op(:Smaxlen),Smaxlen*10,mpi_character,i-1,9,mpi_comm_world,status,mpierr)
           call mpi_recv(mpisync,1,mpi_character,i-1,10,mpi_comm_world,status,mpierr)
          end do
          call mpi_send(inp.feat_name(j:ntot),k*str_len,mpi_character,i-1,1,mpi_comm_world,status,mpierr)
          call mpi_send(inp.lastop(j:ntot),k*10,mpi_character,i-1,2,mpi_comm_world,status,mpierr)
          call mpi_send(inp.feat_comp(j:ntot),k,mpi_integer8,i-1,3,mpi_comm_world,status,mpierr)
          call mpi_send(inp.feat_unit(:,j:ntot),k*nunit,mpi_double_precision,i-1,4,mpi_comm_world,status,mpierr)
          end if
        end do
      else  ! receive
        j=ntot-nthis+1  ! start
        l=ntot-nthis+(nthis-mpin2(mpirank+1)+1)  ! end
        k=nthis-mpin2(mpirank+1)+1  ! size
        if(mpin(mpirank+1)>0) then
          do ll=j,l
           call mpi_recv(inp.Sexpr(ll).list_id(:Smaxlen),Smaxlen,mpi_integer,0,5,mpi_comm_world,status,mpierr)
           call mpi_recv(inp.Sexpr(ll).list_len,1,mpi_integer,0,6,mpi_comm_world,status,mpierr)
           call mpi_recv(inp.Sexpr(ll).list_pointer(:2,:Smaxlen),Smaxlen*2,mpi_integer,0,7,mpi_comm_world,status,mpierr)
           call mpi_recv(inp.Sexpr(ll).list_var(:Smaxlen),Smaxlen*30,mpi_character,0,8,mpi_comm_world,status,mpierr)
           call mpi_recv(inp.Sexpr(ll).list_op(:Smaxlen),Smaxlen*10,mpi_character,0,9,mpi_comm_world,status,mpierr)
           call mpi_send(mpisync,1,mpi_character,0,10,mpi_comm_world,status,mpierr)
          end do
          call mpi_recv(inp.feat_name(j:l),k*str_len,mpi_character,0,1,mpi_comm_world,status,mpierr)
          call mpi_recv(inp.lastop(j:l),k*10,mpi_character,0,2,mpi_comm_world,status,mpierr)
          call mpi_recv(inp.feat_comp(j:l),k,mpi_integer8,0,3,mpi_comm_world,status,mpierr)
          call mpi_recv(inp.feat_unit(:,j:l),k*nunit,mpi_double_precision,0,4,mpi_comm_world,status,mpierr)
        end if
      end if

      ! feature combination
      i=ntot-nthis+1  ! start1
      j=ntot-nthis+mpin(mpirank+1)  ! end1
      k=ntot-nthis+(nthis-mpin2(mpirank+1)+1) ! end2
      total_comb=mpin(mpirank+1)*(ntot-nthis)+(mpin(mpirank+1)-1)*mpin(mpirank+1)/2.0+&
                 mpin(mpirank+1)*(nthis-mpin2(mpirank+1)+1-mpin(mpirank+1))+mpin(mpirank+1) ! binary+unary

      if(mpin(mpirank+1)>0) &
      call combine_se(inp,i,j,1,k,trim(adjustl(ops(icomb))),nf(icomb),total_comb)
      call mpi_barrier(mpi_comm_world,mpierr)

  IF (icomb < rung) THEN

       !---------------------------
       ! redundant check
       !---------------------------
       ! create nfpcore_this: number of new features generated in each CPU core
       if(mpirank/=0) then
          call mpi_send(nf(icomb),1,mpi_integer8,0,1,mpi_comm_world,status,mpierr)
       else
          nfpcore_this(1)=nf(icomb)
          do mpii=1,mpisize-1
               call mpi_recv(nfpcore_this(mpii+1),1,mpi_integer8,mpii,1,mpi_comm_world,status,mpierr)
          end do
       end if
       call mpi_bcast(nfpcore_this,mpisize,mpi_integer8,0,mpi_comm_world,mpierr)

       ! creating a unique scalar number fID for each feature
       allocate(fID(nf(icomb))) 
       do i=1,nf(icomb)
          fID(i)=sqrt(sum((tag+evaluator_se(gen.Sexpr(i)))**2))
       end do

       ! create the "available" to store the repetition information.
        allocate(available(nf(icomb)))
        available=.true.
       ! create the "order" to store the ordering information for bisection method.
        allocate(order(nf(icomb)+1))
 
       !----------------------------------------------------------------------------------------
       ! redundant check inside each core: creating the "order" and "available".
       if(mpirank==0) then
            write(*,'(a,i15)') 'Number of newly generated features: ',sum(nfpcore_this) 
            write(*,'(a)') 'Redundant check on the newly generated features ...'
       end if
        
       ! inside each core
       call dup_scheck_se(nfpcore_this(mpirank+1),fID,gen.feat_name,gen.feat_comp,gen.Sexpr,order,available)
       ! between cores
       if(mpisize>1) call dup_pcheck_se(nfpcore_this,fID,gen.feat_name,gen.feat_comp,gen.Sexpr,order,available)
       !-----------------------------------------------------------------------------------------


       !---------------------------
       ! remove duplicated features
       !---------------------------
       j=0
       do i=1,nf(icomb)
          if(available(i)) then
             j=j+1
             gen.Sexpr(j)=gen.Sexpr(i)
             gen.feat_name(j)=gen.feat_name(i)
             gen.lastop(j)=gen.lastop(i)
             gen.feat_comp(j)=gen.feat_comp(i)
             gen.feat_unit(:,j)=gen.feat_unit(:,i)
             fID(j)=fID(i)
          end if
       end do
       nf(icomb)=j

       !renew nfpcore_this
       if(mpirank/=0) then
          call mpi_send(nf(icomb),1,mpi_integer8,0,1,mpi_comm_world,status,mpierr)
       else
          nfpcore_this(1)=nf(icomb)
          do mpii=1,mpisize-1
               call mpi_recv(nfpcore_this(mpii+1),1,mpi_integer8,mpii,1,mpi_comm_world,status,mpierr)
          end do
       end if
       call mpi_bcast(nfpcore_this,mpisize,mpi_integer8,0,mpi_comm_world,mpierr)

       nthis=sum(nfpcore_this)
       ntot=ntot+nthis
       if(mpirank==0) then
           write(*,'(a,i15)') 'Number of newly generated features after redundant check: ',nthis
       end if
      !--------------------------------------------------------------------
      ! cp results of gen.feat in all cores to the inp.feat in core0
      !--------------------------------------------------------------------
       i=nf(icomb)
       if(mpirank>0 .and. i>0) then
          do ll=1,i
           call mpi_send(gen.Sexpr(ll).list_id(:Smaxlen),Smaxlen,mpi_integer,0,5,mpi_comm_world,status,mpierr)
           call mpi_send(gen.Sexpr(ll).list_len,1,mpi_integer,0,6,mpi_comm_world,status,mpierr)
           call mpi_send(gen.Sexpr(ll).list_pointer(:2,:Smaxlen),Smaxlen*2,mpi_integer,0,7,mpi_comm_world,status,mpierr)
           call mpi_send(gen.Sexpr(ll).list_var(:Smaxlen),Smaxlen*30,mpi_character,0,8,mpi_comm_world,status,mpierr)
           call mpi_send(gen.Sexpr(ll).list_op(:Smaxlen),Smaxlen*10,mpi_character,0,9,mpi_comm_world,status,mpierr)
           call mpi_recv(mpisync,1,mpi_character,0,10,mpi_comm_world,status,mpierr)
          end do
          call mpi_send(gen.feat_name(:i),i*str_len,mpi_character,0,1,mpi_comm_world,status,mpierr)
          call mpi_send(gen.lastop(:i),i*10,mpi_character,0,2,mpi_comm_world,status,mpierr)
          call mpi_send(gen.feat_comp(:i),i,mpi_integer8,0,3,mpi_comm_world,status,mpierr)
          call mpi_send(gen.feat_unit(:,:i),nunit*i,mpi_double_precision,0,4,mpi_comm_world,status,mpierr)
       else if(mpirank==0) then
          if(ntot>ubound(inp.Sexpr,1)) call addm_inp_se(ntot-ubound(inp.Sexpr,1),inp)
          ! from core0 to core0
          inp.Sexpr(ntot-nthis+1:ntot-nthis+i)=gen.Sexpr(:i)
          inp.lastop(ntot-nthis+1:ntot-nthis+i)=gen.lastop(:i)
          inp.feat_comp(ntot-nthis+1:ntot-nthis+i)=gen.feat_comp(:i)
          inp.feat_name(ntot-nthis+1:ntot-nthis+i)=gen.feat_name(:i)
          inp.feat_unit(:,ntot-nthis+1:ntot-nthis+i)=gen.feat_unit(:,:i)
          ! from all other cores to core0
          do mpii=1,mpisize-1
              j=ntot-nthis+sum(nfpcore_this(:mpii))+1    ! start
              k=ntot-nthis+sum(nfpcore_this(:mpii+1))    ! end
              l=nfpcore_this(mpii+1)     ! size
              if(l>0) then
              do ll=j,k
         call mpi_recv(inp.Sexpr(ll).list_id(:Smaxlen),Smaxlen,mpi_integer,mpii,5,mpi_comm_world,status,mpierr)
         call mpi_recv(inp.Sexpr(ll).list_len,1,mpi_integer,mpii,6,mpi_comm_world,status,mpierr)
         call mpi_recv(inp.Sexpr(ll).list_pointer(:2,:Smaxlen),Smaxlen*2,mpi_integer,mpii,7,mpi_comm_world,status,mpierr)
         call mpi_recv(inp.Sexpr(ll).list_var(:Smaxlen),Smaxlen*30,mpi_character,mpii,8,mpi_comm_world,status,mpierr)
         call mpi_recv(inp.Sexpr(ll).list_op(:Smaxlen),Smaxlen*10,mpi_character,mpii,9,mpi_comm_world,status,mpierr)
         call mpi_send(mpisync,1,mpi_character,mpii,10,mpi_comm_world,status,mpierr)
              end do
              call mpi_recv(inp.feat_name(j:k),l*str_len,mpi_character,mpii,1,mpi_comm_world,status,mpierr)
              call mpi_recv(inp.lastop(j:k),l*10,mpi_character,mpii,2,mpi_comm_world,status,mpierr)
              call mpi_recv(inp.feat_comp(j:k),l,mpi_integer8,mpii,3,mpi_comm_world,status,mpierr)
              call mpi_recv(inp.feat_unit(:,j:k),nunit*l,mpi_double_precision,mpii,4,mpi_comm_world,status,mpierr)
              end if
         end do
        end if

       !----------------------------------------------------------------
       ! collect the newly selected features from each core to mpirank0
       !----------------------------------------------------------------
       if(mpirank/=0) then
         call mpi_send(sel.nselect,1,mpi_integer8,0,5,mpi_comm_world,status,mpierr)
       else
         mpik=sel.nselect
         do mpii=1,mpisize-1
           call mpi_recv(mpij,1,mpi_integer8,mpii,5,mpi_comm_world,status,mpierr)
           mpik=mpik+mpij  ! count the total number of selected features
         end do
       end if
       !--------------------------------------
       ! print the space information
       !--------------------------------------
        if(mpirank==0) then
          write(phiname,'(a,i2.2)') 'Phi',icomb
          call writeout_se(phiname,mpik,ntot)
        end if
       ! delete useless space
       deallocate(order)
       deallocate(available)
       deallocate(fID)
  END IF

end do
! -------- end of feature combination ------

! release the spaces
deallocate(inp.Sexpr)
deallocate(inp.feat_name)
deallocate(inp.lastop)
deallocate(inp.feat_comp)
deallocate(inp.feat_unit)

!---------------------------------------------------------------
!  collect the information of selected features from all cores
!---------------------------------------------------------------

if(rung==0) then
 nfpcore_this=sel.nselect
else
 if(mpirank/=0) then
     call mpi_send(nf(rung),1,mpi_integer8,0,1,mpi_comm_world,status,mpierr)
     call mpi_send(sel.nselect,1,mpi_integer8,0,2,mpi_comm_world,status,mpierr)
 else
      nfpcore_this(1)=sel.nselect
      do mpii=1,mpisize-1
           call mpi_recv(mpij,1,mpi_integer8,mpii,1,mpi_comm_world,status,mpierr)
           nf(rung)=nf(rung)+mpij
           call mpi_recv(nfpcore_this(mpii+1),1,mpi_integer8,mpii,2,mpi_comm_world,status,mpierr)
      end do
     write(phiname,'(a,i2.2)') 'Phi',rung
     write(*,'(a,i15)') 'Total number of newly generated features: ',nf(rung)
     call writeout_se(phiname,sum(nfpcore_this),ntot+nf(rung))
 end if
end if

call mpi_bcast(nfpcore_this,mpisize,mpi_integer8,0,mpi_comm_world,mpierr)
call mpi_barrier(mpi_comm_world,mpierr)
!-----------------------------------------------------------------------------

! redundant check for selected features
!---------------------------------------
allocate(available(nfpcore_this(mpirank+1)))
available=.true.
allocate(order(nfpcore_this(mpirank+1)+1))

if(mpirank==0) write(*,'(/a)') 'Redundant check on selected features ...'

! serial redundant check
call dup_scheck_se(nfpcore_this(mpirank+1),sel.feat_score(1,:),sel.feat_name,sel.feat_comp,sel.Sexpr,order,available)

! parallel redundant check
if(mpisize>1) call dup_pcheck_se(nfpcore_this,sel.feat_score(1,:),sel.feat_name,sel.feat_comp,sel.Sexpr,order,available)

!---------------------------------------
! sure independence screening
!---------------------------------------
if(mpirank==0) then
    allocate(sis.Sexpr(nf_sis(iFCDI)))
    allocate(sis.feat_name(nf_sis(iFCDI)))
    allocate(sis.feat_score(2,nf_sis(iFCDI)))
end if

! selecting the best features from sel.XXX to sis.XXX
call sure_indep_screening_se(nfpcore_this,available)

! output the selected feature space
!--------------------------------------
if(mpirank==0) then
   ! output the space.expressions
   if(iFCDI==1) then
     open(fileunit,file='SIS_subspaces/Uspace.expressions',status='replace')
   else
     open(fileunit,file='SIS_subspaces/Uspace.expressions',position='append',status='old')
   end if
   if(ptype==1) then
      do i=1,nf_sis_avai(iFCDI)
        write(fileunit,'(2a,f12.4)') trim(sis.feat_name(i)),'     SIS_score =',sis.feat_score(1,i)
      end do
   else if(ptype==2) then
      do i=1,nf_sis_avai(iFCDI)
        sis.feat_score(1,i)=1.d0/sis.feat_score(1,i)-1.d0   ! score 1: overlap_n, score 2: normalized overlap_length
        if(abs(sis.feat_score(2,i))>1d9) then
           write(fileunit,'(2a,i6,a,f12.4)') trim(sis.feat_name(i)),'    N_overlap =',nint(sis.feat_score(1,i))        
        else
            if( nint(sis.feat_score(1,i))/=0 ) then  ! overlapped
               sis.feat_score(2,i)=1.d0/sis.feat_score(2,i)-1.d0  ! overlapped length
            else   
               sis.feat_score(2,i)=-sis.feat_score(2,i)  ! separation distance
            end if
           write(fileunit,'(2a,i6,a,f12.4)') trim(sis.feat_name(i)),'    N_overlap =',nint(sis.feat_score(1,i)),&
                                                                   '    S_overlap =',sis.feat_score(2,i)
        end if
      end do
   end if
   close(fileunit)

     allocate(tmp_data(npoint,nf_sis_avai(iFCDI)))
     do j=1,nf_sis_avai(iFCDI)
        tmp_data(:,j)=evaluator_se(sis.Sexpr(j))
     end do
   ! output the data
2001  format(*(e20.10))
     do j=1,ntask
        mm1=sum(nsample(:j-1))+1
        mm2=sum(nsample(:j))
        if(ntask>1) then
            write(line,'(a,i3.3,a)') 'Uspace_t',j,'.dat'
        else
            write(line,'(a)') 'Uspace.dat'
        end if
        
        if(iFCDI>1) then
            call rename('SIS_subspaces/'//trim(adjustl(line)),'SIS_subspaces/'//trim(adjustl(line))//'_tmp')
            open(fileunit+1,file='SIS_subspaces/'//trim(adjustl(line))//'_tmp',status='old')
        end if
        open(fileunit,file='SIS_subspaces/'//trim(adjustl(line)),status='replace')

        do ll=mm1,mm2
           superline=''
           if(iFCDI==1) then
             if(ptype==1) then
                  write(fileunit,2001) trainy(ll),tmp_data(ll,:nf_sis_avai(iFCDI))
             elseif(ptype==2) then
                  write(fileunit,2001) tmp_data(ll,:nf_sis_avai(iFCDI))
             endif
           elseif(iFCDI>1) then
             read(fileunit+1,'(a)') superline
             if(ptype==1) then
               write(superline(20*(sum(nf_sis_avai(:iFCDI-1))+1)+1:),2001) tmp_data(ll,:nf_sis_avai(iFCDI)) 
             elseif(ptype==2) then
               write(superline(20*(sum(nf_sis_avai(:iFCDI-1)))+1:),2001) tmp_data(ll,:nf_sis_avai(iFCDI))
             endif
             write(fileunit,'(a)') trim(superline)
          end if
        end do
        close(fileunit)
        if(iFCDI>1) close(fileunit+1,status='delete')
     end do
     deallocate(tmp_data)

   write(*,'(3a,i10)') 'Size of the SIS-selected subspace from ',phiname,': ',nf_sis_avai(iFCDI)
   write(9,'(3a,i10)') 'Size of the SIS-selected subspace from ',phiname,': ',nf_sis_avai(iFCDI)
end if  

! release all the spaces
deallocate(trainy)
deallocate(trainy_c)
deallocate(sel.Sexpr)
deallocate(sel.feat_comp)
deallocate(sel.feat_score)
deallocate(sel.feat_name)
deallocate(available)
deallocate(order)

if(mpirank==0) then
  deallocate(sis.Sexpr)
  deallocate(sis.feat_name)
  deallocate(sis.feat_score)
end if
if(nreject>0) deallocate(reject)

call mpi_barrier(mpi_comm_world,mpierr)
if(mpirank==0) then
   mytime.eFC=mpi_wtime()
   write(9,'(a,f15.2)') 'Time (second) used for this FC: ',mytime.eFC-mytime.sFC
end if

end subroutine


subroutine combine_se(myinp,s1,e1,s2,e2,myops,nf,total_comb) 
implicit none
type(feature) myinp
type(Sexpression) Sexprtmp
real progress
real*8 unit_tmp(ubound(myinp.feat_unit,1))
integer*8 i,j,k,kk,kkk,l,nf,comp_tmp,counter,total_comb,s1,e1,s2,e2
character(len=*) myops
character(len=str_len) name_tmp,lastop_tmp*10
logical skip

counter=0
progress=0.2

do i=s1,e1
  
! no operation
      IF(trim(adjustl(myops))=='NO') THEN
          lastop_tmp=''
          comp_tmp=myinp.feat_comp(i)
          name_tmp='('//trim(adjustl(myinp.feat_name(i)))//')'
          unit_tmp=myinp.feat_unit(:,i)
          Sexprtmp=myinp.Sexpr(i)
          call isgoodf_se(Sexprtmp,name_tmp,lastop_tmp,comp_tmp,unit_tmp,nf)
          cycle
      END IF

! unary operators
      counter=counter+1
      comp_tmp=myinp.feat_comp(i)+1
      if(comp_tmp>fcomplexity) goto 599

      ! exp
      IF(index(myops,'(exp)')/=0 ) then
      if(  index(myinp.lastop(i),'(exp')==0 .and. index(myinp.lastop(i),'(log)')==0 ) then ! avoid exp(exp( and exp(log(
        lastop_tmp='(exp)'
        name_tmp='exp('//trim(adjustl(myinp.feat_name(i)))//')'
        unit_tmp=dimcomb_se(myinp.feat_unit(:,i),myinp.feat_unit(:,i),'(exp)')
        call Sexpr_ucomb_se(myinp.Sexpr(i),Sexprtmp,'(exp)')
        call isgoodf_se(Sexprtmp,name_tmp,lastop_tmp,comp_tmp,unit_tmp,nf)
      end if
      END IF

      ! exp-
      IF(index(myops,'(exp-)')/=0 ) then
      if( index(myinp.lastop(i),'(exp')==0  .and. index(myinp.lastop(i),'(log)')==0 ) then ! avoid exp(exp( and exp(log(
        lastop_tmp='(exp-)'
        name_tmp='exp(-'//trim(adjustl(myinp.feat_name(i)))//')'
        unit_tmp=dimcomb_se(myinp.feat_unit(:,i),myinp.feat_unit(:,i),'(exp-)')
        call Sexpr_ucomb_se(myinp.Sexpr(i),Sexprtmp,'(exp-)')
        call isgoodf_se(Sexprtmp,name_tmp,lastop_tmp,comp_tmp,unit_tmp,nf)
      end if
      END IF

      ! ^-1
      IF(index(myops,'(^-1)')/=0) then
        if(minval(abs(evaluator_se(myinp.Sexpr(i))))>1d-50 ) then  ! avoid divided by zero
           lastop_tmp='(^-1)'
           name_tmp='('//trim(adjustl(myinp.feat_name(i)))//')^-1'
           unit_tmp=dimcomb_se(myinp.feat_unit(:,i),myinp.feat_unit(:,i),'(^-1)')
           call Sexpr_ucomb_se(myinp.Sexpr(i),Sexprtmp,'(^-1)')
           call isgoodf_se(Sexprtmp,name_tmp,lastop_tmp,comp_tmp,unit_tmp,nf)
        end if
       END IF

      ! scd: Standard Cauchy Distribution
      IF(index(myops,'(scd)')/=0) then
           lastop_tmp='(scd)'
           name_tmp='scd('//trim(adjustl(myinp.feat_name(i)))//')'
           unit_tmp=dimcomb_se(myinp.feat_unit(:,i),myinp.feat_unit(:,i),'(scd)')
           call Sexpr_ucomb_se(myinp.Sexpr(i),Sexprtmp,'(scd)')
           call isgoodf_se(Sexprtmp,name_tmp,lastop_tmp,comp_tmp,unit_tmp,nf)
       END IF

      ! ^2
      IF(index(myops,'(^2)')/=0) then
        if(index(myinp.lastop(i),'(sqrt)')==0 ) then ! avoid (sqrt())^2
           lastop_tmp='(^2)'
           name_tmp='('//trim(adjustl(myinp.feat_name(i)))//')^2'
           unit_tmp=dimcomb_se(myinp.feat_unit(:,i),myinp.feat_unit(:,i),'(^2)')
           call Sexpr_ucomb_se(myinp.Sexpr(i),Sexprtmp,'(^2)')
           call isgoodf_se(Sexprtmp,name_tmp,lastop_tmp,comp_tmp,unit_tmp,nf)
        end if
       END IF

      ! ^3
      IF(index(myops,'(^3)')/=0) then
        if(index(myinp.lastop(i),'(cbrt)')==0 ) then ! avoid (cbrt())^3
         lastop_tmp='(^3)'
         name_tmp='('//trim(adjustl(myinp.feat_name(i)))//')^3'
         unit_tmp=dimcomb_se(myinp.feat_unit(:,i),myinp.feat_unit(:,i),'(^3)')
         call Sexpr_ucomb_se(myinp.Sexpr(i),Sexprtmp,'(^3)')
         call isgoodf_se(Sexprtmp,name_tmp,lastop_tmp,comp_tmp,unit_tmp,nf)
       end if
     END IF

      ! ^6
      IF(index(myops,'(^6)')/=0) then
         lastop_tmp='(^6)'
         name_tmp='('//trim(adjustl(myinp.feat_name(i)))//')^6'
         unit_tmp=dimcomb_se(myinp.feat_unit(:,i),myinp.feat_unit(:,i),'(^6)')
         call Sexpr_ucomb_se(myinp.Sexpr(i),Sexprtmp,'(^6)')
         call isgoodf_se(Sexprtmp,name_tmp,lastop_tmp,comp_tmp,unit_tmp,nf)
       END IF

      ! sqrt
      IF(index(myops,'(sqrt)')/=0) then
        if(index(myinp.lastop(i),'(^2)')==0 ) then  ! avoid sqrt((^2))
          if( minval(evaluator_se(myinp.Sexpr(i)))>0  ) then
              lastop_tmp='(sqrt)'     
              name_tmp='sqrt('//trim(adjustl(myinp.feat_name(i)))//')'
              unit_tmp=dimcomb_se(myinp.feat_unit(:,i),myinp.feat_unit(:,i),'(sqrt)')
              call Sexpr_ucomb_se(myinp.Sexpr(i),Sexprtmp,'(sqrt)')
              call isgoodf_se(Sexprtmp,name_tmp,lastop_tmp,comp_tmp,unit_tmp,nf)
          end if
         end if
       END IF

      ! cbrt: cube root
      IF(index(myops,'(cbrt)')/=0) then
        if(index(myinp.lastop(i),'(^3)')==0 ) then  ! avoid cbrt((^3))
              lastop_tmp='(cbrt)'
              name_tmp='cbrt('//trim(adjustl(myinp.feat_name(i)))//')'
              unit_tmp=dimcomb_se(myinp.feat_unit(:,i),myinp.feat_unit(:,i),'(cbrt)')
              call Sexpr_ucomb_se(myinp.Sexpr(i),Sexprtmp,'(cbrt)')
              call isgoodf_se(Sexprtmp,name_tmp,lastop_tmp,comp_tmp,unit_tmp,nf)
         end if
       END IF

      ! log
      IF(index(myops,'(log)')/=0) then
        if( index(myinp.lastop(i),'(exp')==0 .and. index(myinp.lastop(i),'(log)')==0  ) then ! avoid log(exp( and log(log(
          if( minval(evaluator_se(myinp.Sexpr(i)))>0 ) then
              lastop_tmp='(log)'
              name_tmp='log('//trim(adjustl(myinp.feat_name(i)))//')'
              unit_tmp=dimcomb_se(myinp.feat_unit(:,i),myinp.feat_unit(:,i),'(log)')
              call Sexpr_ucomb_se(myinp.Sexpr(i),Sexprtmp,'(log)')
              call isgoodf_se(Sexprtmp,name_tmp,lastop_tmp,comp_tmp,unit_tmp,nf)
          end if
         end if
       END IF

      ! sin
      IF(index(myops,'(sin)')/=0) then
         lastop_tmp='(sin)'
         name_tmp='sin('//trim(adjustl(myinp.feat_name(i)))//')'
         unit_tmp=dimcomb_se(myinp.feat_unit(:,i),myinp.feat_unit(:,i),'(sin)')
         call Sexpr_ucomb_se(myinp.Sexpr(i),Sexprtmp,'(sin)')
         call isgoodf_se(Sexprtmp,name_tmp,lastop_tmp,comp_tmp,unit_tmp,nf)
     END IF

      ! cos
      IF(index(myops,'(cos)')/=0) then
         lastop_tmp='(cos)'
         name_tmp='cos('//trim(adjustl(myinp.feat_name(i)))//')'
         unit_tmp=dimcomb_se(myinp.feat_unit(:,i),myinp.feat_unit(:,i),'(cos)')
         call Sexpr_ucomb_se(myinp.Sexpr(i),Sexprtmp,'(cos)')
         call isgoodf_se(Sexprtmp,name_tmp,lastop_tmp,comp_tmp,unit_tmp,nf)
     END IF

599  continue

  ! binary operators
  do j=s2,e2

      if( j>=s1 .and. j<=i ) cycle

      counter=counter+1
      comp_tmp=myinp.feat_comp(i)+myinp.feat_comp(j)+1
      if(comp_tmp>fcomplexity) goto 602

      ! sum and subtract 
       IF(index(myops,'(+)')/=0 .or. index(myops,'(-)')/=0  .or. index(myops,'(|-|)')/=0   ) THEN

               ! different units
               if ( maxval(abs(myinp.feat_unit(:,i)-myinp.feat_unit(:,j)))>1d-8 ) goto 600

               !---
               IF(index(myops,'(+)')/=0) then
                   lastop_tmp='(+)'
                   name_tmp='('//trim(adjustl(myinp.feat_name(i)))//'+'//trim(adjustl(myinp.feat_name(j)))//')'
                   unit_tmp=dimcomb_se(myinp.feat_unit(:,i),myinp.feat_unit(:,j),'(+)')  
                   call Sexpr_bcomb_se(myinp.Sexpr(i),myinp.Sexpr(j),Sexprtmp,'(+)')
                   call isgoodf_se(Sexprtmp,name_tmp,lastop_tmp,comp_tmp,unit_tmp,nf)
                END IF
                !---

                IF(index(myops,'(-)')/=0) then  ! A-B
                   lastop_tmp='(-)'
                   name_tmp='('//trim(adjustl(myinp.feat_name(i)))//'-'//trim(adjustl(myinp.feat_name(j)))//')'
                   if(index(myops,'(+)')==0) unit_tmp=dimcomb_se(myinp.feat_unit(:,i),myinp.feat_unit(:,j),'(-)')   
                   call Sexpr_bcomb_se(myinp.Sexpr(i),myinp.Sexpr(j),Sexprtmp,'(-)')
                   call isgoodf_se(Sexprtmp,name_tmp,lastop_tmp,comp_tmp,unit_tmp,nf)

                  if(icomb<rung) then  ! B-A (no need if icomb=rung because of the coeff.)
                    ! lastop_tmp='(-)'   ! same with above
                     name_tmp='('//trim(adjustl(myinp.feat_name(j)))//'-'//trim(adjustl(myinp.feat_name(i)))//')'
                    ! unit_tmp=dimcomb_se(inp.feat_unit(:,i),dim_in2(:,j),'(-)')   ! same with above
                     call Sexpr_bcomb_se(myinp.Sexpr(j),myinp.Sexpr(i),Sexprtmp,'(-)')
                     call isgoodf_se(Sexprtmp,name_tmp,lastop_tmp,comp_tmp,unit_tmp,nf)
                  end if
                 END IF

                !-----
                IF(index(myops,'(|-|)')/=0) then  ! absolute difference
                  lastop_tmp='(|-|)'
                  name_tmp='abs('//trim(adjustl(myinp.feat_name(i)))//'-'//trim(adjustl(myinp.feat_name(j)))//')'
                  if(index(myops,'(+)')==0) unit_tmp=dimcomb_se(myinp.feat_unit(:,i),myinp.feat_unit(:,j),'(|-|)')  
                  call Sexpr_bcomb_se(myinp.Sexpr(i),myinp.Sexpr(j),Sexprtmp,'(|-|)')
                  call isgoodf_se(Sexprtmp,name_tmp,lastop_tmp,comp_tmp,unit_tmp,nf)
                 END IF
                 !----
         END IF

600   continue

      ! multiplication and division

           !multiplication
            IF(index(myops,'(*)')/=0) then
                lastop_tmp='(*)'
                name_tmp='('//trim(adjustl(myinp.feat_name(i)))//'*'//trim(adjustl(myinp.feat_name(j)))//')'
                unit_tmp=dimcomb_se(myinp.feat_unit(:,i),myinp.feat_unit(:,j),'(*)')
                call Sexpr_bcomb_se(myinp.Sexpr(i),myinp.Sexpr(j),Sexprtmp,'(*)')
                call isgoodf_se(Sexprtmp,name_tmp,lastop_tmp,comp_tmp,unit_tmp,nf)
            END IF

          ! division
           IF(index(myops,'(/)')/=0) then

              ! avoid A/(B/C)=(A*C)/B and (A/B)/C=A/(B*C) which already exist via previous (*)
              if(index(myops,'(*)')/=0 .and. &
                (index(myinp.lastop(i),'(/)')/=0 .or. index(myinp.lastop(j),'(/)')/=0) .and. &
                (icomb < rung)  )  goto 602   ! do nothing

              skip=.false.
              ! avoid A*B/B and A*B/A
              if(index(myinp.lastop(i),'(*)')/=0) then
                   k=index(myinp.feat_name(i),trim(adjustl(myinp.feat_name(j))))
                 if(k>0) then
                   kk=len_trim(myinp.feat_name(i)(k:))
                   kkk=len_trim(adjustl(myinp.feat_name(j)))
                  if( (trim(adjustl(myinp.feat_name(i)(:k-1)))=='(' &   ! j being the left part of i
                      .and. myinp.feat_name(i)(k+kkk:k+kkk)=='*') .or. &
                      (kk==kkk+1 .and. myinp.feat_name(i)(k-1:k-1)=='*')) then ! j being the right part of i
                      if(index(myops,'(^-1)')/=0) then  ! i/j is no need. what about j/i?
                          goto 602   ! no j/i
                      else
                          goto 601   ! do j/i
                      end if
                   end if
                 end if
              end if
                 
              ! avoid A/(A*B) and B/(A*B)
             if(index(myinp.lastop(j),'(*)')/=0) then
                 k=index(myinp.feat_name(j),trim(adjustl(myinp.feat_name(i))))
                if(k>0) then
                  kk=len_trim(myinp.feat_name(j)(k:))
                  kkk=len_trim(adjustl(myinp.feat_name(i)))
                  if((trim(adjustl(myinp.feat_name(j)(:k-1)))=='(' &
                      .and. myinp.feat_name(j)(k+kkk:k+kkk)=='*') .or. &   ! i being the left part of j
                     (kk==kkk+1 .and. myinp.feat_name(j)(k-1:k-1)=='*')) then  ! i beingn the right part of j
                      if(index(myops,'(^-1)')/=0) then  ! j/i is no need. what about i/j?
                          goto 602  ! no i/j
                      else
                          skip=.true.  ! do i/j
                      end if
                  end if
                end if
             end if

              !i/j
              !---------
              lastop_tmp='(/)'
              name_tmp='('//trim(adjustl(myinp.feat_name(i)))//'/'//trim(adjustl(myinp.feat_name(j)))//')'
              unit_tmp=dimcomb_se(myinp.feat_unit(:,i),myinp.feat_unit(:,j),'(/)')
              if(minval(abs(evaluator_se(myinp.Sexpr(j))))>1d-50 ) then 
              call Sexpr_bcomb_se(myinp.Sexpr(i),myinp.Sexpr(j),Sexprtmp,'(/)')
              call isgoodf_se(Sexprtmp,name_tmp,lastop_tmp,comp_tmp,unit_tmp,nf)
              end if
              !------

     601      continue
              if(skip) goto 602

            !j/i
            !---------
                 lastop_tmp='(/)'
                 name_tmp='('//trim(adjustl(myinp.feat_name(j)))//'/'//trim(adjustl(myinp.feat_name(i)))//')'
                 unit_tmp=dimcomb_se(myinp.feat_unit(:,j),myinp.feat_unit(:,i),'(/)')
                 if(minval(abs(evaluator_se(myinp.Sexpr(i))))>1d-50) then
                 call Sexpr_bcomb_se(myinp.Sexpr(j),myinp.Sexpr(i),Sexprtmp,'(/)')
                 call isgoodf_se(Sexprtmp,name_tmp,lastop_tmp,comp_tmp,unit_tmp,nf)
                 end if
           END if

   602   continue
  end do

  if(float(counter)/float(total_comb)>=progress) then
     do while( (float(counter)/float(total_comb)) >= (progress +0.2) )
           progress = progress+0.2
     end do
     write(*,'(a,i4,2(a,i15),a,f6.1,a)') &
     'mpirank = ',mpirank,'  #generated =',nf,'  #selected =',sel.nselect,'  progress =',progress*100,'%'
     progress=progress+0.2
  end if

end do

end subroutine


function goodf_se(mySexpr,feat_name,feat_unit,feat_comp)
type(Sexpression) mySexpr
integer*8 i,j,k,l,ll,mm1,mm2,feat_comp
real*8 feat_data(npoint),feat_unit(:),scoretmp(2),maxabs
character(len=*) feat_name
logical goodf_se,lsame

goodf_se=.true.
   feat_data=evaluator_se(mySexpr)

   mm1=1
   mm2=sum(nsample(:ntask))
  if(maxval(abs(feat_data(mm1:mm2)-feat_data(mm1)))<=1d-8) then
    goodf_se=.false. ! constant feature.
    return
  end if

  maxabs=maxval(abs(feat_data(mm1:mm2)))
  if(maxabs>1d50 .or. maxabs<=1d-50) then
    goodf_se=.false. ! intinity or zero
    return
  end if

! not selected but can be used for further combination
if(maxabs>fmax_max .or. maxabs<fmax_min) return

! feature score 
if(ptype==1) scoretmp=sis_score_se(mySexpr,trainy_c)
if(ptype==2) scoretmp=sis_score_se(mySexpr,trainy) !classification

if (scoretmp(1)<score_threshold) return

! reject the previously selected features
if(nreject>0) then
  feat_name=adjustl(feat_name)
  lsame=.false.
  i=0; j=nreject; k=i+ceiling((j-i)/2.0)
  do while(i/=j)
     if(trim(feat_name)==trim(reject(k))) then
       lsame=.true.
       i=j
     else if(feat_name<reject(k)) then
       i=k;
     else if(feat_name>reject(k)) then
       j=k;
     end if

     if(k==i+ceiling((j-i)/2.0)) then
         i=j
     else
         k=i+ceiling((j-i)/2.0)
     end if
  end do
  if(lsame) return
end if

!--------------------------
! selected
!--------------------------
sel.nselect=sel.nselect+1
sel.Sexpr(sel.nselect)=mySexpr
sel.feat_comp(sel.nselect)=feat_comp
sel.feat_name(sel.nselect)=feat_name
sel.feat_score(:,sel.nselect)=scoretmp
if( sel.nselect== nbasic_select+nextra_select ) call update_select_se

end function


subroutine addm_gen_se(n,gen)
! increase array size
real*8,allocatable:: real2d(:,:)
integer*8,allocatable:: integer1d(:)
character(len=str_len),allocatable:: char1d(:),char1d2(:)*10
integer*8 i,j,k,n
type(feature) gen
type(Sexpression),allocatable:: Sexprtmp(:)
i=npoint
j=ubound(gen.Sexpr,1)

! gen.Sexpr
allocate(Sexprtmp(j))
Sexprtmp=gen.Sexpr
deallocate(gen.Sexpr)
allocate(gen.Sexpr(j+n))
gen.Sexpr(:j)=Sexprtmp
deallocate(Sexprtmp)
!-----
!gen.feat_name
allocate(char1d(j))
char1d=gen.feat_name
deallocate(gen.feat_name)
allocate(gen.feat_name(j+n))
gen.feat_name(:j)=char1d
deallocate(char1d)
!----
!gen.lastop
allocate(char1d2(j))
char1d2=gen.lastop
deallocate(gen.lastop)
allocate(gen.lastop(j+n))
gen.lastop(:j)=char1d2
deallocate(char1d2)
!----
!gen.feat_comp
allocate(integer1d(j))
integer1d=gen.feat_comp
deallocate(gen.feat_comp)
allocate(gen.feat_comp(j+n))
gen.feat_comp(:j)=integer1d
deallocate(integer1d)
!gen.feat_unit
i=ubound(gen.feat_unit,1)
allocate(real2d(i,j))
real2d=gen.feat_unit
deallocate(gen.feat_unit)
allocate(gen.feat_unit(i,j+n))
gen.feat_unit(:,:j)=real2d
deallocate(real2d)
!--
end subroutine


subroutine addm_inp_se(n,inp)
! increase array size
real*8,allocatable:: real2d(:,:)
type(feature) inp
type(Sexpression),allocatable:: Sexprtmp(:)
integer*8,allocatable:: integer1d(:)
character(len=str_len),allocatable:: char1d(:),char1d2(:)*10
integer*8 i,j,k,n
i=npoint
j=ubound(inp.Sexpr,1)

! inp.Sexpr
allocate(Sexprtmp(j))
Sexprtmp=inp.Sexpr
deallocate(inp.Sexpr)
allocate(inp.Sexpr(j+n))
inp.Sexpr(:j)=Sexprtmp
deallocate(Sexprtmp)
!-----
!inp.feat_name
allocate(char1d(j))
char1d=inp.feat_name
deallocate(inp.feat_name)
allocate(inp.feat_name(j+n))
inp.feat_name(:j)=char1d
deallocate(char1d)
!----
!inp.lastop
allocate(char1d2(j))
char1d2=inp.lastop
deallocate(inp.lastop)
allocate(inp.lastop(j+n))
inp.lastop(:j)=char1d2
deallocate(char1d2)
!---
!inp.feat_comp
allocate(integer1d(j))
integer1d=inp.feat_comp
deallocate(inp.feat_comp)
allocate(inp.feat_comp(j+n))
inp.feat_comp(:j)=integer1d
deallocate(integer1d)
!inp.feat_unit
i=ubound(inp.feat_unit,1)
allocate(real2d(i,j))
real2d=inp.feat_unit
deallocate(inp.feat_unit)
allocate(inp.feat_unit(i,j+n))
inp.feat_unit(:,:j)=real2d
deallocate(real2d)
!--
end subroutine


function dimcomb_se(dim1,dim2,op)
! calculate the units for new features
! unary operator: set dim1 and dim2 the same
real*8 dim1(:),dim2(:),dimcomb_se(ubound(dim1,1))
character(len=*) op
integer i,j,k

if(trim(adjustl(op))=='(+)' .or. trim(adjustl(op))=='(-)' .or. trim(adjustl(op))=='(|-|)'  ) then
   do i=1,ubound(dim1,1)
      if(abs(dim1(i)-dim2(i))>1d-8) stop 'Error: different units in linear combination'
   end do
   dimcomb_se=dim1
else if(trim(adjustl(op))=='(*)') then
   dimcomb_se=dim1+dim2
else if(trim(adjustl(op))=='(/)') then
   dimcomb_se=dim1-dim2
else if(trim(adjustl(op))=='(exp)') then   
   dimcomb_se=0.d0
else if(trim(adjustl(op))=='(exp-)') then
   dimcomb_se=0.d0
else if(trim(adjustl(op))=='(log)') then
   dimcomb_se=0.d0
else if(trim(adjustl(op))=='(scd)') then 
   dimcomb_se=0.d0
else if(trim(adjustl(op))=='(sin)') then
   dimcomb_se=0.d0
else if(trim(adjustl(op))=='(cos)') then
   dimcomb_se=0.d0
else if(trim(adjustl(op))=='(^-1)') then
   dimcomb_se=dim1*(-1)
else if(trim(adjustl(op))=='(^2)') then
   dimcomb_se=dim1*2
else if(trim(adjustl(op))=='(^3)') then
   dimcomb_se=dim1*3
else if(trim(adjustl(op))=='(^6)') then
   dimcomb_se=dim1*6
else if(trim(adjustl(op))=='(sqrt)') then
   dimcomb_se=dim1/2.d0
else if(trim(adjustl(op))=='(cbrt)') then
   dimcomb_se=dim1/3.d0
end if
end function


subroutine writeout_se(phiname,i,k)
integer*8 i,j,k
character(len=*) phiname
write(*,'(3a,i15)') 'Total number of features in the space ',trim(phiname),':',k
write(9,'(3a,i15)') 'Total number of features in the space ',trim(phiname),':',k
end subroutine


function simpler_se(complexity1,complexity2,name1,name2)
integer*8 simpler_se,complexity1,complexity2
character(len=*) name1,name2

if(complexity1<complexity2) then
   simpler_se=1
else if(complexity1>complexity2) then
   simpler_se=2
else if(complexity1==complexity2) then
   if(len_trim(name1)<len_trim(name2) .or. &
     (len_trim(name1)==len_trim(name2) .and. trim(name1)<=trim(name2) ) ) then
      simpler_se=1
   else
      simpler_se=2
   end if
end if
end function


subroutine update_availability_se(available,pavailable)
logical available(:),pavailable(:)
if( any(available) ) then
   pavailable(mpirank+1)=.true.
else
   pavailable(mpirank+1)=.false.
end if
end subroutine


subroutine dup_pcheck_se(nfpcore_this,fID,fname,complexity,mySexpr,order,available)
! output array 'available'
real*8 fID(:)
integer*8 i,j,k,l,ll,mpii,mpij,nfpcore_this(:),mpim(mpisize),mpiloc(mpisize),order(:),loc(1),&
          complexity(:),simpler_result
character(len=*) fname(:)
logical available(:)
real*8,allocatable:: compidentity(:)
character(len=str_len),allocatable:: compname(:)
integer*8,allocatable:: compcomplexity(:)
type(Sexpression),allocatable:: Sexpr4comp(:)
type(Sexpression) mySexpr(:)

IF(mpisize>1) THEN
   mpim=nfpcore_this  ! number of features in each core
   do i=1,mpisize 
      loc(1:1)=maxloc(mpim)
      mpiloc(i)=loc(1)  ! sort the cores by the number of features from high to low
      mpim(loc(1))=-1
   end do

   do k=1,mpisize-1

     if(nfpcore_this(mpiloc(k))>0) then
         allocate(compidentity(nfpcore_this(mpiloc(k))))
         allocate(compname(nfpcore_this(mpiloc(k))))
         allocate(compcomplexity(nfpcore_this(mpiloc(k))))
         allocate(Sexpr4comp(nfpcore_this(mpiloc(k))))

         if(mpirank==mpiloc(k)-1) then   ! every core knows the same mpiloc
           compidentity=fID(:nfpcore_this(mpiloc(k)))
           compname=fname(:nfpcore_this(mpiloc(k)))
           compcomplexity=complexity(:nfpcore_this(mpiloc(k)))
           Sexpr4comp=mySexpr(:nfpcore_this(mpiloc(k)))
         end if

         ! broadcast compidentity from core mpiloc(k)-1 to all other cores
          call mpi_bcast(compidentity,nfpcore_this(mpiloc(k)),mpi_double_precision,mpiloc(k)-1,mpi_comm_world,mpierr)
         ! broadcast the feature expressions
          call mpi_bcast(compname,nfpcore_this(mpiloc(k))*str_len,mpi_character,mpiloc(k)-1,mpi_comm_world,mpierr)
         ! broadcast the feature complexities
          call mpi_bcast(compcomplexity,nfpcore_this(mpiloc(k)),mpi_integer8,mpiloc(k)-1,mpi_comm_world,mpierr)
         ! broadcast the features
          do ll=1,nfpcore_this(mpiloc(k))
           call mpi_bcast(Sexpr4comp(ll).list_id(:Smaxlen),Smaxlen,mpi_integer,mpiloc(k)-1,mpi_comm_world,mpierr)
           call mpi_bcast(Sexpr4comp(ll).list_len,1,mpi_integer,mpiloc(k)-1,mpi_comm_world,mpierr)
           call mpi_bcast(Sexpr4comp(ll).list_pointer(:2,:Smaxlen),Smaxlen*2,mpi_integer,mpiloc(k)-1,mpi_comm_world,mpierr)
           call mpi_bcast(Sexpr4comp(ll).list_var(:Smaxlen),Smaxlen*30,mpi_character,mpiloc(k)-1,mpi_comm_world,mpierr)
           call mpi_bcast(Sexpr4comp(ll).list_op(:Smaxlen),Smaxlen*10,mpi_character,mpiloc(k)-1,mpi_comm_world,mpierr)
          end do
         ! do the comparision. All the cores of mpiloc(k:mpisize) will be compared with
         ! that on the core mpiloc(k)-1. 
         if( all( mpirank/=(mpiloc(:k)-1) ) .and. nfpcore_this(mpirank+1)>0 ) then ! this core vs core mpiloc(k)-1
            do mpij=1,nfpcore_this(mpiloc(k))  ! core mpiloc(k)-1
              l=0; ll=order(nfpcore_this(mpirank+1)+1)  ! order stores features in descending order 
              i=l+ceiling(float(ll-l)/2.0)   ! bisection method

              124 continue
              simpler_result=simpler_se(compcomplexity(mpij),complexity(order(i)),compname(mpij),fname(order(i)))

                 ! check for duplication
                 if(equivalent_se(compidentity(mpij),fID(order(i)),Sexpr4comp(mpij),mySexpr(order(i)) ) ) then 
                       if(simpler_result==1) then
                          available(order(i))=.false.
                       else
                          compidentity(mpij)=0.d0
                       end if
                    cycle
                 end if
                 ! if not equal
                 if(compidentity(mpij)>fID(order(i)) .or.  (compidentity(mpij)==fID(order(i)) .and. &
                    simpler_result==1 ) ) then   ! mpij is better than i
                    ll=i   ! replace the right end with i to find a new middle point
                 else
                     l=i  ! replace the left end with i to findn a new middle point
                 end if
                 if(i==l+ceiling(float(ll-l)/2.0)) cycle  ! i is already the left or right end; no more
                 i=l+ceiling(float(ll-l)/2.0)   ! the new middle point
                 goto 124
            end do
            call mpi_send(compidentity,nfpcore_this(mpiloc(k)),mpi_double_precision,mpiloc(k)-1,33,&
                          mpi_comm_world,status,mpierr)

         else if(mpirank==mpiloc(k)-1) then
             do l=1,mpisize
                 if( any( l==mpiloc(:k) ) .or. nfpcore_this(l)==0 ) cycle
            call mpi_recv(compidentity,nfpcore_this(mpiloc(k)),mpi_double_precision,mpi_any_source,33,&
                      mpi_comm_world,status,mpierr)
               do i=1,nfpcore_this(mpiloc(k))
                if(abs(compidentity(i))<1d-8) available(i)=.false.
               end do
            end do
         end if

         deallocate(compidentity)
         deallocate(compname)
         deallocate(compcomplexity)
         deallocate(Sexpr4comp)
     end if

  end do

  j=0
  do i=1,nfpcore_this(mpirank+1)
   if(available(i)) j=j+1
  end do

  if(mpirank/=0)  then
    call mpi_send(j,1,mpi_integer8,0,1,mpi_comm_world,status,mpierr)
  else
    do k=1,mpisize-1
      call mpi_recv(i,1,mpi_integer8,k,1,mpi_comm_world,status,mpierr)
      j=j+i
    end do
  end if

END IF

end subroutine

subroutine sure_indep_screening_se(nfpcore_this,available)
! sure independence screening
real*8 tmp,pscore(2,mpisize)
integer*8 i,j,k,l,ll,mpii,mpij,nfpcore_this(:),loc(2),simpler_result,pcomplexity(mpisize)
character(len=str_len) pname(mpisize)
logical available(:),pavailable(mpisize)


if(mpirank==0) write(*,'(/a)') 'Sure Independence Screening on selected features ...'

pavailable=.false.
if(nfpcore_this(mpirank+1)>0) call update_availability_se(available,pavailable)
do l=1,mpisize
 call mpi_bcast(pavailable(l),1,mpi_logical,l-1,mpi_comm_world,mpierr) 
end do

! selection starts ...
i=0   ! count of selected features
do while( any(pavailable) .and. i<nf_sis(iFCDI) )

      i=i+1

      ! find the max score on each core
      if(nfpcore_this(mpirank+1)>0) then
        loc(2:2)=maxloc(sel.feat_score(1,:sel.nselect)) ! score_2 is less important or not used
        tmp=sel.feat_score(1,loc(2))  ! (loc(1) to be used for other purpose
      end if
      do l=1,nfpcore_this(mpirank+1)  ! equal scores
        if(  abs(tmp-sel.feat_score(1,l))<=1d-8 ) then
           if( sel.feat_score(2,l)-sel.feat_score(2,loc(2))>1d-8 .or. &
               (abs(sel.feat_score(2,l)-sel.feat_score(2,loc(2)))<=1d-8 .and. &
            simpler_se(sel.feat_comp(l),sel.feat_comp(loc(2)),sel.feat_name(l),sel.feat_name(loc(2)))==1 ) )  loc(2)=l
        end if
      end do
      if(nfpcore_this(mpirank+1)>0) then
        pscore(:,mpirank+1)=sel.feat_score(:,loc(2))  ! location of max score of this core
        pcomplexity(mpirank+1)=sel.feat_comp(loc(2))  ! corresponding feature complexity
        pname(mpirank+1)=sel.feat_name(loc(2))  ! corresponding feature name
      else
        pscore(:,mpirank+1)=-1   ! simply a negative value to indicate empty
      end if

      do l=1,mpisize  ! inform other cores about local best features
        call mpi_bcast(pscore(:,l),2,mpi_double_precision,l-1,mpi_comm_world,mpierr)
        call mpi_bcast(pcomplexity(l),1,mpi_integer8,l-1,mpi_comm_world,mpierr)
        call mpi_bcast(pname(l),str_len,mpi_character,l-1,mpi_comm_world,mpierr)
      end do

      !find the max score between cores
      loc(1:1)=maxloc(pscore(1,:))
      tmp=maxval(pscore(1,:))
      do l=1,mpisize  ! if equal score between features
        if(  abs(tmp-pscore(1,l))<=1d-8 ) then
           if( pscore(2,l)-pscore(2,loc(1))>1d-8 .or.  ( abs(pscore(2,l)-pscore(2,loc(1)))<=1d-8 .and. &
               simpler_se(pcomplexity(l),pcomplexity(loc(1)),pname(l),pname(loc(1)))==1 ) ) loc(1)=l
        end if
      end do

      ! save the highest-scored feature
      if((loc(1)-1)==mpirank) then
         if(mpirank==0) then
             sis.Sexpr(i)=sel.Sexpr(loc(2))
             sis.feat_name(i)=sel.feat_name(loc(2))
             sis.feat_score(:,i)=sel.feat_score(:,loc(2))
         else
        call mpi_send(sel.Sexpr(loc(2)).list_id(:Smaxlen),Smaxlen,mpi_integer,0,96,mpi_comm_world,status,mpierr)
        call mpi_send(sel.Sexpr(loc(2)).list_len,1,mpi_integer,0,97,mpi_comm_world,status,mpierr)
        call mpi_send(sel.Sexpr(loc(2)).list_pointer(:2,:Smaxlen),Smaxlen*2,mpi_integer,0,98,mpi_comm_world,status,mpierr)
        call mpi_send(sel.Sexpr(loc(2)).list_var(:Smaxlen),Smaxlen*30,mpi_character,0,99,mpi_comm_world,status,mpierr)
        call mpi_send(sel.Sexpr(loc(2)).list_op(:Smaxlen),Smaxlen*10,mpi_character,0,100,mpi_comm_world,status,mpierr)
        call mpi_send(sel.feat_name(loc(2)),str_len,mpi_character,0,101,mpi_comm_world,status,mpierr)
        call mpi_send(sel.feat_score(:,loc(2)),2,mpi_double_precision,0,103,mpi_comm_world,status,mpierr)
         end if
         available(loc(2))=.false.   ! avoid to be selected again
         sel.feat_score(:,loc(2))=-1  

      end if
      if(mpirank==0 .and. mpirank/=(loc(1)-1) ) then
     call mpi_recv(sis.Sexpr(i).list_id(:Smaxlen),Smaxlen,mpi_integer,loc(1)-1,96,mpi_comm_world,status,mpierr)
     call mpi_recv(sis.Sexpr(i).list_len,1,mpi_integer,loc(1)-1,97,mpi_comm_world,status,mpierr)
     call mpi_recv(sis.Sexpr(i).list_pointer(:2,:Smaxlen),Smaxlen*2,mpi_integer,loc(1)-1,98,mpi_comm_world,status,mpierr)
     call mpi_recv(sis.Sexpr(i).list_var(:Smaxlen),Smaxlen*30,mpi_character,loc(1)-1,99,mpi_comm_world,status,mpierr)
     call mpi_recv(sis.Sexpr(i).list_op(:Smaxlen),Smaxlen*10,mpi_character,loc(1)-1,100,mpi_comm_world,status,mpierr)
     call mpi_recv(sis.feat_name(i),str_len,mpi_character,loc(1)-1,101,mpi_comm_world,status,mpierr)
     call mpi_recv(sis.feat_score(:,i),2,mpi_double_precision,loc(1)-1,103,mpi_comm_world,status,mpierr)
      end if
      !---

      if(nfpcore_this(mpirank+1)>0) call update_availability_se(available,pavailable)
      do l=1,mpisize
       call mpi_bcast(pavailable(l),1,mpi_logical,l-1,mpi_comm_world,mpierr)
      end do

end do
!--------------------


nf_sis_avai(iFCDI)=i    ! the actual number of selected features

end subroutine


subroutine dup_scheck_se(num,fID,fname,complexity,mySexpr,order,available)
! duplication check within each core
! output the arrays "order" and "available"
integer*8 i,j,l,ll,order(:),n,num,complexity(:),simpler_result
real*8 fID(:)
character(len=*) fname(:)
logical available(:)
type(Sexpression) mySexpr(:)

IF(num==0) THEN
  n=0
ELSE

! creating the descending ordering via bisection method 
order(1)=1  ! initial first entry
n=1
do i=2,num

   l=0; ll=n;   ! left and right end
   j=l+ceiling(float(ll-l)/2.0)  ! the middle one

   123 continue
   simpler_result=simpler_se(complexity(i),complexity(order(j)),fname(i),fname(order(j)))

   ! check for duplication
      if(equivalent_se(fID(i),fID(order(j)),mySexpr(i),mySexpr(order(j))) ) then
           if( simpler_result==1 ) then
              available(order(j))=.false.
              order(j)=i
           else
              available(i)=.false.
           end if
           cycle
      end if
     ! if not equal
     if(fID(i)>fID(order(j)) .or. (fID(i)==fID(order(j)) .and. simpler_result==1 ) ) then  ! i is better
        ll=j    ! replace the ll with the j to find a new middle point
        if(j==l+ceiling(float(ll-l)/2.0)) then   ! if j is already the left end
          order(j+1:n+1)=order(j:n)  
          order(j)=i   ! insert i to that before the original j
          n=n+1
          cycle
        end if
      else
         l=j   ! j is better; replace the left end with j to find a new middle point
         if(j==l+ceiling(float(ll-l)/2.0)) then  ! if j is already the right end
           if(n>j) order(j+2:n+1)=order(j+1:n)
           order(j+1)=i   ! insert i to that after j
           n=n+1
           cycle
         end if
      end if

      j=l+ceiling(float(ll-l)/2.0)  ! the new middle point
      goto 123
end do

END IF

order(num+1)=n  ! store the number of features after check
if(mpirank/=0)  then
  call mpi_send(n,1,mpi_integer8,0,1,mpi_comm_world,status,mpierr)
else
  do l=1,mpisize-1
    call mpi_recv(i,1,mpi_integer8,l,1,mpi_comm_world,status,mpierr)
    n=n+i
 end do
end if

end subroutine


function sis_score_se(mySexpr,yyy)   
! correlation between a feature 'feat' and the target 'yyy'
! sis_score_se returns a vector with 2 elements
integer i,j,mm1,mm2,mm3,mm4,k,kk,l,overlap_n,nf1,nf2,itask,nconvexpair
real*8 feat(npoint),sdfeat(npoint),tmp(ntask),sis_score_se(2),yyy(:),xnorm(ntask),xmean(ntask),&
       overlap_length,length_tmp,feat_tmp1(npoint),feat_tmp2(npoint),mindist,minlen
logical isoverlap
type(Sexpression) mySexpr

feat=evaluator_se(mySexpr)

if(ptype==1) then  
     do j=1,ntask
         mm1=sum(nsample(:j-1))+1
         mm2=sum(nsample(:j))
         xmean(j)=sum(feat(mm1:mm2))/nsample(j)
         sdfeat(mm1:mm2)=feat(mm1:mm2)-xmean(j)  
         xnorm(j)=sqrt(sum((sdfeat(mm1:mm2))**2))
         if(xnorm(j)>1d-50) then
            sdfeat(mm1:mm2)=sdfeat(mm1:mm2)/xnorm(j)  ! standardization to the features
            tmp(j)=abs(sum(sdfeat(mm1:mm2)*yyy(mm1:mm2))) ! |xy|. Donot use normalized yyy here!
         else
            tmp(j)=0.d0
         end if
     end do
   ! score ranges from 0 to 1
   sis_score_se(1)=sqrt(sum(tmp**2)/ntask)  ! quadratic mean of the scores of different tasks
   sis_score_se(1)=sis_score_se(1)/sqrt(sum(yyy**2)/ntask)  ! normalization
   sis_score_se(2)=1.d0   ! not used

else if(ptype==2) then  
   mindist=-1d10
   overlap_n=0
   overlap_length=0.d0
   isoverlap=.false.
   nconvexpair=0
   do itask=1,ntask

   ! calculate overlap between domains of property i
        do i=1,ngroup(itask,1000)-1   ! ngroup(itask,1000) store the number of groups in this task
            if(itask==1) then
              mm1=sum(ngroup(itask,:i-1))+1
              mm2=sum(ngroup(itask,:i))
            else            
              mm1=sum(nsample(:itask-1))+sum(ngroup(itask,:i-1))+1
              mm2=sum(nsample(:itask-1))+sum(ngroup(itask,:i))
            end if
            nf1=0
            do k=mm1,mm2
               if(yyy(k)<0.5) then  ! y=1: classified, y=0: unclassified
                   nf1=nf1+1
                   feat_tmp1(nf1)=feat(k) 
               end if
            end do       
          do j=i+1,ngroup(itask,1000) 
               if(itask==1) then
                 mm3=sum(ngroup(itask,:j-1))+1
                 mm4=sum(ngroup(itask,:j))
               else
                 mm3=sum(nsample(:itask-1))+sum(ngroup(itask,:j-1))+1
                 mm4=sum(nsample(:itask-1))+sum(ngroup(itask,:j))
               end if
               nf2=0
               do k=mm3,mm4
                   if(yyy(k)<0.5) then
                       nf2=nf2+1
                       feat_tmp2(nf2)=feat(k) 
                   end if
               end do
               if(isconvex(itask,i)==0 .and. isconvex(itask,j)==0) cycle ! both are not convex domains

               if(isconvex(itask,i)==1 .and. isconvex(itask,j)==1) then
                  call convex1d_overlap(feat_tmp1(:nf1),feat_tmp2(:nf2),bwidth,k,length_tmp)
                  overlap_n=overlap_n+k
                  nconvexpair=nconvexpair+1
                  if(length_tmp>=0.d0) isoverlap=.true.
                   ! which feature is shorter
                   minlen=min(maxval(feat_tmp1(:nf1))-minval(feat_tmp1(:nf1)),maxval(feat_tmp2(:nf2))&
                          -minval(feat_tmp2(:nf2)))
                   if(length_tmp<0.d0) then  ! if separated
                      if(length_tmp>mindist) mindist=length_tmp  ! renew the worst separation
                   else if(length_tmp>=0.d0 .and. minlen==0.d0) then  ! if overlapped and one feature is 0D
                      overlap_length=overlap_length+1.d0  ! totally overlapped
                   else if(length_tmp>=0.d0 .and. minlen>0.d0) then  ! if not separated and no 0D feature
                      overlap_length=overlap_length+length_tmp/minlen  ! calculate total overlap
                   end if
               else if (isconvex(itask,i)==0 .and. isconvex(itask,j)==1) then  !count the number of i-data in j-domain
                  overlap_n=overlap_n+convex1d_in(feat(mm3:mm4),feat_tmp1(:nf1),bwidth)
               else if (isconvex(itask,i)==1 .and. isconvex(itask,j)==0) then !count the number of j-data in i-domain
                  overlap_n=overlap_n+convex1d_in(feat(mm1:mm2),feat_tmp2(:nf2),bwidth)
               end if

          end do  ! j
        end do  ! i
   end do  ! itask

   sis_score_se(1)=float(overlap_n)  ! >=0, higher sis_score --> worse feature
   sis_score_se(1)=1.d0/(sis_score_se(1)+1.d0)  ! transform to <=1, higher sis_score --> better feature
   if(isoverlap) then  ! there are domains overlapped
     sis_score_se(2)=overlap_length/float(nconvexpair)     ! second metric, <=1, higher --> worse
     sis_score_se(2)=1.d0/(sis_score_se(2)+1.0)  ! transform to <=1, higher --> better
   else ! separated length
     sis_score_se(2)=-mindist   ! separated length between domains,positive
   end if

end if
end function


subroutine isgoodf_se(mySexpr,feat_name,lastop,feat_comp,feat_unit,nf)
real*8 feat_unit(:)
character(len=*) feat_name,lastop
integer*8 nf,feat_comp
type(Sexpression) mySexpr

if(goodf_se(mySexpr,feat_name,feat_unit,feat_comp)) then
  nf=nf+1
  if(icomb < rung) then
     if(nf>ubound(gen.Sexpr,1)) call addm_gen_se(int8(ceiling(100000.0/mpisize)),gen)
     gen.Sexpr(nf)=mySexpr
     gen.feat_name(nf)=feat_name
     gen.lastop(nf)=lastop
     gen.feat_comp(nf)=feat_comp
     gen.feat_unit(:,nf)=feat_unit
  end if
end if
end subroutine


subroutine update_select_se
! update the selected space
! bisection method for descending order
integer*8 i,j,k,l,ll,order(sel.nselect),n,tmp,tmpcomplexity(nbasic_select),simpler_result
real*8 tmpscore(2,nbasic_select)
character(len=str_len) tmpname(nbasic_select)
type(Sexpression) Sexprtmp(nbasic_select)

! creating an ordering from large to small
order(1)=1   ! assuming the first feature being the best (highest score)
n=1  ! count of features saved in 'order'

do i=2,sel.nselect   ! compare feaure i with the j located at middle of order(1:n)

  l=0; ll=n;
  j=l+ceiling(float(ll-l)/2.0)   ! j is at middle between l and ll.
 
  125 continue
  simpler_result=simpler_se(sel.feat_comp(i),sel.feat_comp(order(j)),sel.feat_name(i),sel.feat_name(order(j)))

  ! get rid of duplicated features
  if(equivalent_se(sel.feat_score(1,i),sel.feat_score(1,order(j)),sel.Sexpr(i),sel.Sexpr(order(j))) ) then
        if( simpler_result==1) order(j)=i  ! update the order
        cycle
  end if

  ! bisection method for descending order
  if( (sel.feat_score(1,i)>sel.feat_score(1,order(j))) .or.  ((sel.feat_score(1,i)==sel.feat_score(1,order(j))) &
       .and. sel.feat_score(2,i)>sel.feat_score(2,order(j))) .or.  ((sel.feat_score(1,i)==sel.feat_score(1,order(j))) &
       .and. sel.feat_score(2,i)==sel.feat_score(2,order(j)) .and. simpler_result==1) )then 
      ll=j   ! i is better. Replace the right end ll with j, and find a new middle point
      if(j==l+ceiling(float(ll-l)/2.0)) then  ! if j is already the left end
        order(j+1:n+1)=order(j:n)   ! move all the original j:n features by 1 step
        order(j)=i ! insert the i right before the j, because score_i > score_j
        n=n+1  ! size increased by 1
        cycle
      end if
   else   ! j is better. Replace the left end l with j, and find a new middle point
       l=j
       if(j==l+ceiling(float(ll-l)/2.0)) then
         if(n>j) order(j+2:n+1)=order(j+1:n)
         order(j+1)=i
         n=n+1
         cycle
       end if
   end if

   j=l+ceiling(float(ll-l)/2.0)  ! the new middle point
   goto 125
end do

sel.nselect=min(n,nbasic_select)   ! reduce the space by removing features with low ranking

! reordering
do i=1,sel.nselect
Sexprtmp(i)=sel.Sexpr(order(i))
tmpcomplexity(i)=sel.feat_comp(order(i))
tmpname(i)=sel.feat_name(order(i))
tmpscore(:,i)=sel.feat_score(:,order(i))
end do

! update the selected features
score_threshold=tmpscore(1,sel.nselect)
sel.Sexpr(:sel.nselect)=Sexprtmp(:sel.nselect)
sel.feat_comp(:sel.nselect)=tmpcomplexity(:sel.nselect)
sel.feat_name(:sel.nselect)=tmpname(:sel.nselect)
sel.feat_score(:,:sel.nselect)=tmpscore(:,:sel.nselect)

end subroutine


function equivalent_se(score1,score2,mySexpr1,mySexpr2)
! check if two features are the same or highly correlated.
type(Sexpression) mySexpr1,mySexpr2
real*8 score1,score2,diff,feat1(npoint),feat2(npoint),mean1,mean2,sd1,sd2,ffcorr
logical equivalent_se

equivalent_se=.false.
diff=score1-score2
if(abs(diff)>1d-8) return

feat1=evaluator_se(mySexpr1)
feat2=evaluator_se(mySexpr2)

mean1=sum(feat1)/float(npoint)
mean2=sum(feat2)/float(npoint)
sd1=sqrt((sum((feat1-mean1)**2))/float(npoint))
sd2=sqrt((sum((feat2-mean2)**2))/float(npoint))
ffcorr=(sum((feat1-mean1)*(feat2-mean2)))/float(npoint)/sd1/sd2

if(ffcorr>1.d0) ffcorr=1.d0    ! avoid numerical noise
if(ffcorr<-1.d0) ffcorr=-1.d0

if( abs(ffcorr)> 0.99 ) equivalent_se=.true.

end function


subroutine Sexpr_bcomb_se(Sexpr1,Sexpr2,Sexpr12,op)
type(Sexpression) Sexpr1,Sexpr2,Sexpr12
integer i,j,k1,k2
character(len=*) op

     k1=Sexpr1.list_len
     k2=Sexpr2.list_len
     Sexpr12=Sexpr1
     Sexpr12.list_id(k1+1:k1+k2)=Sexpr2.list_id(1:k2)+k1
     Sexpr12.list_op(k1+1:k1+k2)=Sexpr2.list_op(1:k2)
     Sexpr12.list_var(k1+1:k1+k2)=Sexpr2.list_var(1:k2)
     do i=1,k2
       do j=1,2
        if(Sexpr2.list_pointer(j,i)/=0) then
           Sexpr12.list_pointer(j,k1+i)=Sexpr2.list_pointer(j,i)+k1
        else
           Sexpr12.list_pointer(j,k1+i)=0
        end if
      end do
     end do

     Sexpr12.list_len=k1+k2+1
     Sexpr12.list_id(Sexpr12.list_len)=Sexpr12.list_len
     Sexpr12.list_var(Sexpr12.list_len)=' '   
     Sexpr12.list_pointer(:,Sexpr12.list_len)=(/k1,k1+k2/)
     select case(trim(adjustl(op)))
        case('(+)')
          Sexpr12.list_op(Sexpr12.list_len)='(+)'
        case('(-)')
          Sexpr12.list_op(Sexpr12.list_len)='(-)'
        case('(*)')
          Sexpr12.list_op(Sexpr12.list_len)='(*)'
        case('(/)')
          Sexpr12.list_op(Sexpr12.list_len)='(/)'
        case('(|-|)')
          Sexpr12.list_op(Sexpr12.list_len)='(|-|)'
      end select
end subroutine


subroutine Sexpr_ucomb_se(Sexpr1,Sexpr2,op)
type(Sexpression) Sexpr1,Sexpr2
integer i,j,k1
character(len=*) op

     Sexpr2=Sexpr1
     Sexpr2.list_len=Sexpr2.list_len+1
     Sexpr2.list_id(Sexpr2.list_len)=Sexpr2.list_len
     Sexpr2.list_var(Sexpr2.list_len)=' '
     Sexpr2.list_pointer(:,Sexpr2.list_len)=(/Sexpr2.list_len-1,0/)
     select case(trim(adjustl(op)))
        case('(exp)')
          Sexpr2.list_op(Sexpr2.list_len)='(exp)'
        case('(exp-)')
          Sexpr2.list_op(Sexpr2.list_len)='(exp-)'
        case('(^-1)')
          Sexpr2.list_op(Sexpr2.list_len)='(^-1)'
        case('(^2)')
          Sexpr2.list_op(Sexpr2.list_len)='(^2)'
        case('(^3)')
          Sexpr2.list_op(Sexpr2.list_len)='(^3)'
        case('(sqrt)')
          Sexpr2.list_op(Sexpr2.list_len)='(sqrt)'
        case('(cbrt)')
          Sexpr2.list_op(Sexpr2.list_len)='(cbrt)'
        case('(log)')
          Sexpr2.list_op(Sexpr2.list_len)='(log)'
        case('(scd)')
          Sexpr2.list_op(Sexpr2.list_len)='(scd)'
        case('(^6)')
          Sexpr2.list_op(Sexpr2.list_len)='(^6)'
        case('(sin)')
          Sexpr2.list_op(Sexpr2.list_len)='(sin)'
        case('(cos)')
          Sexpr2.list_op(Sexpr2.list_len)='(cos)'
     end select

end subroutine

function evaluator_se(mySexpr)
! evaluator_se
integer i,j,k,l,length
type(Sexpression) mySexpr
real*8 evaluator_se(npoint)
real*8 val(npoint,mySexpr.list_len)

length=mySexpr.list_len
val=0.0
do i=1,length
  if(index(mySexpr.list_op(i),'var')/=0) then
     do j=1,nsf
        if(trim(adjustl(mySexpr.list_var(i)))==trim(adjustl(pfname(j)))) then
            val(:,i)=pfdata(:,j)
            exit
        end if
     end do
  else if (index(mySexpr.list_op(i),'(+)')/=0) then
     val(:,i)=val(:,mySexpr.list_pointer(1,i))+val(:,mySexpr.list_pointer(2,i))
  else if (index(mySexpr.list_op(i),'(-)')/=0) then
     val(:,i)=val(:,mySexpr.list_pointer(1,i))-val(:,mySexpr.list_pointer(2,i))
  else if (index(mySexpr.list_op(i),'(*)')/=0) then
     val(:,i)=val(:,mySexpr.list_pointer(1,i))*val(:,mySexpr.list_pointer(2,i))
  else if (index(mySexpr.list_op(i),'(/)')/=0) then
     val(:,i)=val(:,mySexpr.list_pointer(1,i))/val(:,mySexpr.list_pointer(2,i))
  else if (index(mySexpr.list_op(i),'(|-|)')/=0) then
     val(:,i)=abs(val(:,mySexpr.list_pointer(1,i))-val(:,mySexpr.list_pointer(2,i)))
  else if (index(mySexpr.list_op(i),'(exp)')/=0) then
     val(:,i)=exp(val(:,mySexpr.list_pointer(1,i)))
  else if (index(mySexpr.list_op(i),'(exp-)')/=0) then
     val(:,i)=exp(-val(:,mySexpr.list_pointer(1,i)))
  else if (index(mySexpr.list_op(i),'(^2)')/=0) then
     val(:,i)=(val(:,mySexpr.list_pointer(1,i)))**2
  else if (index(mySexpr.list_op(i),'(^-1)')/=0) then
     val(:,i)=(val(:,mySexpr.list_pointer(1,i)))**(-1)
  else if (index(mySexpr.list_op(i),'(^3)')/=0) then
     val(:,i)=(val(:,mySexpr.list_pointer(1,i)))**3
  else if (index(mySexpr.list_op(i),'(sqrt)')/=0) then
     val(:,i)=sqrt(val(:,mySexpr.list_pointer(1,i)))
  else if (index(mySexpr.list_op(i),'(cbrt)')/=0) then
     val(:,i)=(val(:,mySexpr.list_pointer(1,i)))**(1.d0/3.d0)
  else if (index(mySexpr.list_op(i),'(log)')/=0) then
     val(:,i)=log(val(:,mySexpr.list_pointer(1,i)))
  else if (index(mySexpr.list_op(i),'(scd)')/=0) then
     val(:,i)=1.0d0/(PI*(1.0d0+(val(:,mySexpr.list_pointer(1,i)))**2)) 
  else if (index(mySexpr.list_op(i),'(^6)')/=0) then
     val(:,i)=(val(:,mySexpr.list_pointer(1,i)))**6
  else if (index(mySexpr.list_op(i),'(sin)')/=0) then
     val(:,i)=sin(val(:,mySexpr.list_pointer(1,i)))
  else if (index(mySexpr.list_op(i),'(cos)')/=0) then
     val(:,i)=cos(val(:,mySexpr.list_pointer(1,i)))
  end if
end do
evaluator_se(:)=val(:,length)

end function

subroutine printSlist_se(mySexpr)
type(Sexpression) mySexpr
integer i,j,k
k=mySexpr.list_len
do i=1,k
write(*,'(a,i3.3,a,a,2i5,2a)') '(',mySexpr.list_id(i),')  ',&
      trim(mySexpr.list_op(i)),mySexpr.list_pointer(:,i),'  ',trim(mySexpr.list_var(i))
end do
end subroutine


end module

