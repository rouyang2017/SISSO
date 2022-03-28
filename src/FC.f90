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


module FC
! feature construction
!----------------------------------------------------------------------

use var_global
use libsisso
! variables for this module
integer*8 ntot,nthis,nreject,nvf_new,nfinal_select,nbasic_select,nextra_select,nselect,icomb
real*8,allocatable:: trainy(:),trainy_c(:),fout(:,:),dim_out(:,:),vfeat(:,:,:),f_select(:,:),score_select(:,:),tag(:)
character(len=lname),allocatable::  name_out(:),lastop_out(:)*10,name_select(:),reject(:)
integer*8,allocatable:: complexity_out(:),mpin(:),mpin2(:),complexity_select(:),usedup_warning(:)
logical threshold_select,foutsave

contains


subroutine feature_construction
implicit none
integer   loc(1),nvf_old,ioerr
logical   lreject
character line*500,phiname*5,rejectname*100
real*8    tmp,bbb,tmp_score(2)
integer*8 aaa,i,j,k,l,ll,mm1,mm2,nf(maxcomb),mpii,mpij,mpik
real*8,allocatable:: ftag(:,:),SD(:),sisfeat(:,:),feat_in1(:,:),dim_in1(:,:)
integer*8,allocatable:: nfpcore_this(:),nfpcore_tot(:),order(:),nvfpcore(:),complexity_in1(:)
character(len=lname),allocatable:: name_sisfeat(:),name_in1(:),lastop_in1(:)*10
logical,allocatable:: available(:)


! Stop if the whole feature space had been selected.
IF (iFCDI>1 .and. nsis(iFCDI-1)<subs_sis(iFCDI-1)) THEN
   if(mpirank==0) then
      write(*,'(a)') 'Warning: Last selected subspace < subs_sis! No more FC will be performed.' 
      write(9,'(a)') 'Wrning: Last selected subspace < subs_sis! No more FC will be performed.'
   end if
   return
end if

! running time for FC
if(mpirank==0) stime_FC=mpi_wtime()

!------------------------------------------------------
! array allocation
!------------------------------------------------------
allocate(trainy(npoints))
allocate(trainy_c(npoints))
i=max(1000,nsf+nvf)
allocate(feat_in1(sum(nsample),i))   
allocate(name_in1(i))
allocate(lastop_in1(i))
allocate(dim_in1(ndimtype,i))
allocate(complexity_in1(i))

allocate(fout(sum(nsample),i))
allocate(name_out(i))
allocate(lastop_out(i))
allocate(complexity_out(i))
allocate(dim_out(ndimtype,i))
if(nvf>0) allocate(vfeat(sum(nsample),vfsize,i))

if(ffdecorr) then
  if(mpisize>1) then
     nbasic_select=ceiling(decorr_alpha*subs_sis(iFCDI))
  else if(mpisize==1) then
     nbasic_select=2*ceiling(decorr_alpha*subs_sis(iFCDI))
  end if
else
  if(mpisize>1) then
     nbasic_select=subs_sis(iFCDI)
  else if(mpisize==1) then
     nbasic_select=2*subs_sis(iFCDI)
  end if
end if

nextra_select=min(subs_sis(iFCDI),50000)
j=nbasic_select+nextra_select

allocate(f_select(npoints,j))
allocate(score_select(j,2))
allocate(name_select(j))
allocate(complexity_select(j))

allocate(mpin(mpisize))
allocate(mpin2(mpisize))
allocate(nfpcore_this(mpisize))
allocate(nfpcore_tot(mpisize))
allocate(SD(ntask))
allocate(nvfpcore(mpisize))
allocate(tag(sum(nsample))) 
allocate(usedup_warning(mpisize))

! initialization
!------------------
dim_in1(:,:nsf+nvf)=pfdim
threshold_select=.false.
nvf_new=nvf
feat_in1(:,1:nsf)=psfeat(:,1:nsf)
name_in1(:nsf+nvf)=pfname(:nsf+nvf)
trainy=res

! A vector to asist redundant check
do j=1,ntask
do i=1,nsample(j)    
  tag(sum(nsample(:j-1))+i)=1.0d0+0.001d0*i
end do
end do

if(nvf>0) then  ! there are vectors features
  vfeat(:,:,:nvf)=pvfeat(:,:,:nvf)
  do i=1,nvf  
    feat_in1(1,nsf+i)=i    ! vf ID to the 1st element
    feat_in1(2,nsf+i)=sqrt(sum((tag+sum(vfeat(:,:,i),2))**2))  ! norm of the new vector
  end do
end if

! avoid repetition of features from a existing space
inquire(file='feature_space/reject.name',exist=lreject)
if(lreject) then
   rejectname='feature_space/reject.name'
   if(mpirank==0) write(9,'(2a)') 'File containing the features to be rejected: ',trim(rejectname)
else
 if(iFCDI>1) inquire(file='feature_space/Uspace.name',exist=lreject)
 if(lreject) then
  rejectname='feature_space/Uspace.name'
!   if(mpirank==0) write(9,'(2a)') 'File containing the features to be rejected: ',trim(rejectname)
 end if
end if

nreject=0
if(lreject .and. mpirank==0) then
   open(funit,file=trim(adjustl(rejectname)),status='old')

   ! check file length
   do while(.true.)
      read(funit,*,iostat=ioerr)
      if(ioerr<0) exit
      nreject=nreject+1
   end do

   allocate(reject(nreject))
   rewind(funit)
   
   do j=1,nreject
      read(funit,*) line
      call string_split(line,reject(j:j),' ')
      reject(j)=adjustl(reject(j))
   end do
   close(funit)
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
if(nreject>0) call mpi_bcast(reject,nreject*lname,mpi_character,0,mpi_comm_world,mpierr)

!------------------------------------------------------------
! Population Standard Deviation
! For rough guess of fitting quality by comparing with RMSE
!------------------------------------------------------------
2000 format(a,i3.3,a,*(f10.5))

if(nreaction==0) then
   do j=1,ntask
    mm1=sum(nsample(:j-1))+1
    mm2=sum(nsample(:j))
    trainy_c(mm1:mm2)=trainy(mm1:mm2)-sum(trainy(mm1:mm2))/(mm2-mm1+1)  ! centered
    SD(j)=sqrt(sum((trainy_c(mm1:mm2))**2)/(mm2-mm1+1))  ! population standard deviation
    if(mpirank==0 .and. ptype==1 .and. iFCDI==1) &
        write(9,2000) 'Population Standard Deviation (SD) of task ',j,': ',SD(j)
   end do
elseif(nreaction>0) then
    trainy_c(1:nreaction)=trainy(1:nreaction)-sum(trainy(1:nreaction))/nreaction  ! centered
    SD(j)=sqrt(sum((trainy_c(1:nreaction))**2)/nreaction)  ! population standard deviation
    if(mpirank==0 .and. ptype==1 .and. iFCDI==1) &
        write(9,2000) 'Population Standard Deviation (SD) of task ',j,': ',SD(j)
end if

!-------------------------------------------------------------------------
! feature combinations start ...
! \PHI0 is the initial primary feature space
!-------------------------------------------------------------------------
nselect=0 ! # of selected features
nf=0      ! # generated features from each combination
foutsave=.true.   ! save the generated features for this combination ?
nvf_old=nvf_new   ! update # vfeat
lastop_in1=''     ! the last operation of the generated feature
complexity_in1=0

i=nsf+nvf ! total number of primary features
j=0

! no combination, but to check if primary features are good to be included in the 'selected' set.
call combine(feat_in1(:,1:i),name_in1(1:i),lastop_in1(1:i),complexity_in1(1:i),dim_in1(:,1:i),&  
             feat_in1(:,1:i),name_in1(1:i),lastop_in1(1:i),complexity_in1(1:i),dim_in1(:,1:i),'NO',j) 

nvf_new=nvf_old  
ntot=nsf+nvf_new
nthis=ntot
write(phiname,'(a,i2.2)') 'phi',0
if(mpirank==0)  call writeout(phiname,nselect,nvf_new,ntot)

if(rung==0) then
  deallocate(fout)
  deallocate(name_out)
  deallocate(lastop_out)
  deallocate(complexity_out)
  deallocate(dim_out)
end if

!-------------------------
!creating PHI1, PHI2, ...
!-------------------------
do icomb=1,rung
    if(mpirank==0) write(*,'(/a,i2.2,a)') 'Generating phi',icomb,' ...'

    call mpi_barrier(mpi_comm_world,mpierr)

    if(icomb==rung) then
       foutsave=.false.
       deallocate(fout)
       deallocate(name_out)
       deallocate(lastop_out)
       deallocate(complexity_out)
       deallocate(dim_out)
    end if

    nvf_old=nvf_new

    !allocating jobs for each core
    ! divide the combination into n parts with equal load
    ! X1, X2, ..., Xi, ...
    ! #X1=X1(X1-1)/2+X1(N-X1) = N(N-1)/2/n -> X1
    ! #X2=X2(X2-1)/2+X2(N-X1-X2)
    ! #Xi=Xi(Xi-1)/2+Xi(N-X1-X2-...-Xi)
    ! let a=sum(#X1+#X2+...+#X{i-1}),b=#PHI*(#PHI-1)/2/n
    ! #Xi = #PHI-1/2-a-sqrt((PHI-1/2-a)^2-2*b)
    ! mpin: size for each core, mpin2: starting id
      mpin=0
      mpin2=0
      bbb=nthis*(nthis-1)/2.d0/float(mpisize)
      do mpii=1,mpisize-1
         aaa=sum(mpin(:mpii-1))
         if((nthis-aaa-0.5)**2-2*bbb>0) mpin(mpii)=ceiling(nthis-aaa-0.5-sqrt((nthis-aaa-0.5)**2-2*bbb))
      end do
      mpin(mpisize)=max(nthis-sum(mpin),0)
      do i=1,mpisize
        if(mpin(i)>0) mpin2(i)=sum(mpin(:i-1))+1
      end do

     ! balance memory_requirment for each node
     ! exchange large and small ones
      do i=2,mpisize,2
         j=mpisize-i+2
         if(j<=i) exit
         mpii=mpin(i)
         mpij=mpin2(i)
         mpin(i)=mpin(j)
         mpin2(i)=mpin2(j)
         mpin(j)=mpii
         mpin2(j)=mpij
      end do

     ! check if need more memory 
      i=ntot-nthis+(nthis-mpin2(mpirank+1)+1)
      if(ubound(feat_in1,2)< i)  call addm_in1(i-ubound(feat_in1,2),feat_in1,name_in1,lastop_in1,complexity_in1,dim_in1)

     ! broadcast ntot-nthis
       if(ntot>nthis) then
        i=ntot-nthis
        call mpi_bcast(feat_in1(:,:i),i*sum(nsample),mpi_double_precision,0,mpi_comm_world,mpierr)
        call mpi_bcast(name_in1(:i),i*lname,mpi_character,0,mpi_comm_world,mpierr)
        call mpi_bcast(lastop_in1(:i),i*10,mpi_character,0,mpi_comm_world,mpierr)
        call mpi_bcast(dim_in1(:,:i),i*ndimtype,mpi_double_precision,0,mpi_comm_world,mpierr)
        call mpi_bcast(complexity_in1(:i),i,mpi_integer8,0,mpi_comm_world,mpierr)
       end if

      if(mpirank==0) then  ! send nthis 
        do i=2,mpisize
          j=ntot-nthis+mpin2(i)  ! start
          k=nthis-mpin2(i)+1   ! size
          if(mpin(i)>0) then
          call mpi_send(feat_in1(:,j:ntot),k*sum(nsample),mpi_double_precision,i-1,31,mpi_comm_world,status,mpierr)
          call mpi_send(name_in1(j:ntot),k*lname,mpi_character,i-1,32,mpi_comm_world,status,mpierr)
          call mpi_send(lastop_in1(j:ntot),k*10,mpi_character,i-1,33,mpi_comm_world,status,mpierr)
          call mpi_send(complexity_in1(j:ntot),k,mpi_integer8,i-1,34,mpi_comm_world,status,mpierr)
          call mpi_send(dim_in1(:,j:ntot),k*ndimtype,mpi_double_precision,i-1,35,mpi_comm_world,status,mpierr)
          end if
        end do
      else  ! receive
        j=ntot-nthis+1  ! start
        l=ntot-nthis+(nthis-mpin2(mpirank+1)+1)  ! end
        k=nthis-mpin2(mpirank+1)+1  ! size
        if(mpin(mpirank+1)>0) then
        call mpi_recv(feat_in1(:,j:l),k*sum(nsample),mpi_double_precision,0,31,mpi_comm_world,status,mpierr)
        call mpi_recv(name_in1(j:l),k*lname,mpi_character,0,32,mpi_comm_world,status,mpierr)
        call mpi_recv(lastop_in1(j:l),k*10,mpi_character,0,33,mpi_comm_world,status,mpierr)
        call mpi_recv(complexity_in1(j:l),k,mpi_integer8,0,34,mpi_comm_world,status,mpierr)
        call mpi_recv(dim_in1(:,j:l),k*ndimtype,mpi_double_precision,0,35,mpi_comm_world,status,mpierr)
        end if
      end if

      ! combination
      i=ntot-nthis+1  ! start1
      j=ntot-nthis+mpin(mpirank+1)  ! end1
      k=ntot-nthis+(nthis-mpin2(mpirank+1)+1) ! end2
      if(mpin(mpirank+1)>0) &
      call combine(feat_in1(:,i:j),name_in1(i:j),lastop_in1(i:j),complexity_in1(i:j),dim_in1(:,i:j),&
                   feat_in1(:,1:k),name_in1(1:k),lastop_in1(1:k),complexity_in1(1:k),dim_in1(:,1:k),& 
                   trim(adjustl(opset(icomb))),nf(icomb))

      call mpi_barrier(mpi_comm_world,mpierr)

  IF (foutsave) THEN

       !---------------------------
       ! redundant check
       !---------------------------
       ! create nfpcore_this: # new features generated on each core
       if(mpirank/=0) then
          call mpi_send(nf(icomb),1,mpi_integer8,0,1,mpi_comm_world,status,mpierr)
       else
          nfpcore_this(1)=nf(icomb)
          do mpii=1,mpisize-1
               call mpi_recv(nfpcore_this(mpii+1),1,mpi_integer8,mpii,1,mpi_comm_world,status,mpierr)
          end do
       end if
       call mpi_bcast(nfpcore_this,mpisize,mpi_integer8,0,mpi_comm_world,mpierr)

       ! create identity for the new features
       allocate(ftag(nf(icomb),2)) 
       do i=1,nf(icomb)
         if( isscalar(name_out(i)) ) then
          ftag(i,1)=sqrt(sum((tag+fout(:,i))**2))
         else
          ftag(i,1)=fout(2,i)
         end if
       end do

       ! create the "available" to store the repetition information.
       ! available(i)=.false.: duplicate
        allocate(available(nf(icomb)))
        available=.true.
       ! create the "order" to store the ordering information for bisection method.
        allocate(order(nf(icomb)+1))
 
       !----------------------------------------------------------------------------------------
       ! redundant check on each core: output order and available.
       if(mpirank==0) then
            write(*,'(a,i15)') 'Number of newly generated features: ',sum(nfpcore_this) 
            write(*,'(a)') 'Redundant check on the newly generated features ...'
       end if
       ! serial duplication check
       call dup_scheck(nfpcore_this(mpirank+1),ftag,name_out,complexity_out,order,available,'all')

       ! parallel duplication check
       ! redundant check between cores,output available
       if(mpisize>1) call dup_pcheck(nfpcore_this,ftag,name_out,complexity_out,order,available,'all')
       !-----------------------------------------------------------------------------------------

       ! renew the nvf_new and nvfpcore on each core
       if(mpirank==0) then
          nvf_new=nvf_old
       else
          nvf_new=0
       end if
       do i=1,nf(icomb)
          if(available(i)) then
             if( .not. isscalar(name_out(i)) ) nvf_new=nvf_new+1
          end if
       end do

       if(mpirank/=0) then
          call mpi_send(nvf_new,1,mpi_integer8,0,1,mpi_comm_world,status,mpierr)
       else
          nvfpcore(1)=nvf_new
          do mpii=1,mpisize-1
               call mpi_recv(nvfpcore(mpii+1),1,mpi_integer8,mpii,1,mpi_comm_world,status,mpierr)
          end do
       end if
       call mpi_bcast(nvfpcore,mpisize,mpi_integer8,0,mpi_comm_world,mpierr)

       !---------------------------
       ! remove duplicate features
       !---------------------------
       j=0
       if(mpirank==0) then
          nvf_new=nvf_old
       else
          nvf_new=0
       end if
       do i=1,nf(icomb)
          if(available(i)) then
             j=j+1
             fout(:,j)=fout(:,i)
             name_out(j)=name_out(i)
             lastop_out(j)=lastop_out(i)
             complexity_out(j)=complexity_out(i)
             dim_out(:,j)=dim_out(:,i)
             ftag(j,1)=ftag(i,1)
             if(nvf>0) then
             if( .not. isscalar(name_out(j)) ) then
               nvf_new=nvf_new+1
               vfeat(:,:,nvf_new)=vfeat(:,:,nint(fout(1,j)))
               if(mpirank==0) fout(1,j)=nvf_new
               if(mpirank/=0) fout(1,j)=nvf_new+sum(nvfpcore(:mpirank))
             end if
             end if
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
       nfpcore_tot=nfpcore_tot+nfpcore_this
       if(mpirank==0) then
           write(*,'(a,i15)') 'Number of newly generated features after redundant check: ',nthis
       end if
      !--------------------------------------------------------------------
      ! cp results from _out to _in1 of mpirank0
      !--------------------------------------------------------------------
       i=nf(icomb)
       if(mpirank>0 .and. i>0) then
          call mpi_send(fout(:,:i),sum(nsample)*i,mpi_double_precision,0,2,mpi_comm_world,status,mpierr)
          call mpi_send(name_out(:i),i*lname,mpi_character,0,3,mpi_comm_world,status,mpierr)
          call mpi_send(lastop_out(:i),i*10,mpi_character,0,4,mpi_comm_world,status,mpierr)
          call mpi_send(complexity_out(:i),i,mpi_integer8,0,5,mpi_comm_world,status,mpierr)
          call mpi_send(dim_out(:,:i),ndimtype*i,mpi_double_precision,0,6,mpi_comm_world,status,mpierr)
       else if(mpirank==0) then
         if(ntot>ubound(feat_in1,2)) &
                 call addm_in1(ntot-ubound(feat_in1,2),feat_in1,name_in1,lastop_in1,complexity_in1,dim_in1)
         feat_in1(:,ntot-nthis+1:ntot-nthis+i)=fout(:,:i)
         lastop_in1(ntot-nthis+1:ntot-nthis+i)=lastop_out(:i)
         complexity_in1(ntot-nthis+1:ntot-nthis+i)=complexity_out(:i)
         name_in1(ntot-nthis+1:ntot-nthis+i)=name_out(:i)
         dim_in1(:,ntot-nthis+1:ntot-nthis+i)=dim_out(:,:i)
         do mpii=1,mpisize-1
              j=ntot-nthis+sum(nfpcore_this(:mpii))+1    ! start
              k=ntot-nthis+sum(nfpcore_this(:mpii+1))    ! end
              l=nfpcore_this(mpii+1)     ! size
              if(l>0) then
              call mpi_recv(feat_in1(:,j:k),sum(nsample)*l,mpi_double_precision,mpii,2,mpi_comm_world,status,mpierr)
              call mpi_recv(name_in1(j:k),l*lname,mpi_character,mpii,3,mpi_comm_world,status,mpierr)
              call mpi_recv(lastop_in1(j:k),l*10,mpi_character,mpii,4,mpi_comm_world,status,mpierr)
              call mpi_recv(complexity_in1(j:k),l,mpi_integer8,mpii,5,mpi_comm_world,status,mpierr)
              call mpi_recv(dim_in1(:,j:k),ndimtype*l,mpi_double_precision,mpii,6,mpi_comm_world,status,mpierr)
              end if
         end do
        end if

        !------------------------------------------------------
        ! broadcast vector features 
        !------------------------------------------------------
         if(nvf>0) then
            if(mpirank/=0) then
               if(nvf_new>0) call mpi_send(vfeat(:,:,:nvf_new),sum(nsample)*vfsize*nvf_new,&
                        mpi_double_precision,0,7,mpi_comm_world,status,mpierr)
            else
              if(sum(nvfpcore)>ubound(vfeat,3)) call addm_vf(sum(nvfpcore)-ubound(vfeat,3)) ! vector feature
              do mpii=1,mpisize-1
                   i=ntot-nthis+sum(nfpcore_this(:mpii))+1
                   j=ntot-nthis+sum(nfpcore_this(:mpii+1))
                   if(nvfpcore(mpii+1)>0) call mpi_recv(vfeat(:,:,nvf_new+1:nvf_new+nvfpcore(mpii+1)),&
                     sum(nsample)*vfsize*nvfpcore(mpii+1),mpi_double_precision,mpii,7,mpi_comm_world,status,mpierr)
                   nvf_new=nvf_new+nvfpcore(mpii+1)
              end do
             end if

             nvf_new=sum(nvfpcore)
             if(nvf_new>0 ) then
                if(nvf_new>ubound(vfeat,3)) call addm_vf(nvf_new-ubound(vfeat,3))
                k=100000000.d0/sum(nsample)/vfsize
                i=nvf_new/k
                j=mod(nvf_new,k)
                do mpii=1,i
                   call mpi_bcast(vfeat(:,:,(mpii-1)*k+1:mpii*k),sum(nsample)*vfsize*k,mpi_double_precision,0,&
                                   mpi_comm_world,mpierr)
                end do
                call mpi_bcast(vfeat(:,:,nvf_new-j+1:nvf_new),sum(nsample)*vfsize*j,mpi_double_precision,0,&
                     mpi_comm_world,mpierr)
             end if
          end if

       !----------------------------------------------------------------
       ! collect the newly selected features from each core to mpirank0
       !----------------------------------------------------------------
       if(mpirank/=0) then
         call mpi_send(nselect,1,mpi_integer8,0,5,mpi_comm_world,status,mpierr)
       else
         mpik=nselect
         do mpii=1,mpisize-1
           call mpi_recv(mpij,1,mpi_integer8,mpii,5,mpi_comm_world,status,mpierr)
           mpik=mpik+mpij  ! count the total number of selected features
         end do
       end if
       !--------------------------------------
       ! print the space information
       !--------------------------------------
        if(mpirank==0) then
          write(phiname,'(a,i2.2)') 'phi',icomb
          call writeout(phiname,mpik,nvf_new,ntot)
        end if
       ! delete useless space
       deallocate(order)
       deallocate(available)
       deallocate(ftag)
  END IF

end do
! -------- end of feature combination ------

! release spaces
deallocate(feat_in1)
deallocate(name_in1)
deallocate(lastop_in1)
deallocate(complexity_in1)
deallocate(dim_in1)

if(nvf>0) deallocate(vfeat)

!---------------------------------------------------------------
!  collect the information of selected features from all cores
!---------------------------------------------------------------

if(rung==0) then
 nfpcore_this=nselect
else
 if(mpirank/=0) then
     call mpi_send(nf(rung),1,mpi_integer8,0,1,mpi_comm_world,status,mpierr)
     call mpi_send(nselect,1,mpi_integer8,0,2,mpi_comm_world,status,mpierr)
     call mpi_send(nvf_new-nvf_old,1,mpi_integer8,0,3,mpi_comm_world,status,mpierr)
 else
      nfpcore_this(1)=nselect
      nvfpcore(1)=nvf_new
      do mpii=1,mpisize-1
           call mpi_recv(mpij,1,mpi_integer8,mpii,1,mpi_comm_world,status,mpierr)
           nf(rung)=nf(rung)+mpij
           call mpi_recv(nfpcore_this(mpii+1),1,mpi_integer8,mpii,2,mpi_comm_world,status,mpierr)
           call mpi_recv(nvfpcore(mpii+1),1,mpi_integer8,mpii,3,mpi_comm_world,status,mpierr)
      end do
     write(phiname,'(a,i2.2)') 'phi',rung
     write(*,'(a,i15)') 'Total number of newly generated features: ',nf(rung)
     call writeout(phiname,sum(nfpcore_this),sum(nvfpcore),ntot+nf(rung))
 end if
end if

call mpi_bcast(nfpcore_this,mpisize,mpi_integer8,0,mpi_comm_world,mpierr)
call mpi_barrier(mpi_comm_world,mpierr)
!-----------------------------------------------------------------------------
! saved information for selected features on each core ->
! nselect: number of features selected on each core
! f_select: selected features data on each core
! complexity_select: complexity of the selected features on each core
! name_select: name of the selected features on each core
! score_select: score of the selected features on each core

! redundant check for selected features
!---------------------------------------
allocate(available(nfpcore_this(mpirank+1)))
available=.true.
allocate(order(nfpcore_this(mpirank+1)+1))

! serial check
call dup_scheck(nfpcore_this(mpirank+1),score_select,name_select,complexity_select,order,available,'selected')

! parallel redundant check
if(mpisize>1) call dup_pcheck(nfpcore_this,score_select,name_select,complexity_select,order,available,'selected')

!---------------------------------------
! sure independence screening
!---------------------------------------
if(mpirank==0) then
    allocate(sisfeat(npoints,subs_sis(iFCDI)))
    allocate(name_sisfeat(subs_sis(iFCDI)))
end if

call sure_indep_screening(nfpcore_this,available,complexity_select,name_select,f_select,sisfeat,name_sisfeat)

!--------------------------------------
!--------------------------------------
! output of the selected feature space
!--------------------------------------
!--------------------------------------
if(mpirank==0) then
   ! output space.name
   write(line,'(a,i3.3,a)') 'space_',iFCDI,'d.name'
   open(funit,file='feature_space/'//trim(adjustl(line)),status='replace')
   if(ptype==1) then
      do i=1,nsis(iFCDI)
        tmp_score=sis_score(sisfeat(:,i),trainy_c)
        write(funit,'(2a,f12.4)') trim(name_sisfeat(i)),'  corr=',tmp_score(1)
      end do
   else if(ptype==2) then
      do i=1,nsis(iFCDI)
        tmp_score=sis_score(sisfeat(:,i),trainy) ! for classification
        tmp_score(1)=1.d0/tmp_score(1)-1.d0   ! score 1: overlap_n, score 2: normalized overlap_length
        if(abs(tmp_score(2))>1d9) then
           write(funit,'(2a,i6,a,f12.4)') trim(name_sisfeat(i)),'  N=',nint(tmp_score(1))        
        else
            if( nint(tmp_score(1))/=0 ) then  ! overlapped
               tmp_score(2)=1.d0/tmp_score(2)-1.d0  ! overlapped length
            else   
               tmp_score(2)=-tmp_score(2)  ! separation distance
            end if
           write(funit,'(2a,i6,a,f12.4)') trim(name_sisfeat(i)),'  N=',nint(tmp_score(1)),'  S=',tmp_score(2)
        end if
      end do
   end if
   close(funit)

   ! output space data
2001  format(*(e20.10))
   if(nreaction==0) then
     do j=1,ntask
        mm1=sum(nsample(:j-1))+1
        mm2=sum(nsample(:j))
        write(line,'(a,i3.3,a,i3.3,a)') 'space_',iFCDI,'d_p',j,'.dat'
        open(funit,file='feature_space/'//trim(adjustl(line)),status='replace')
        if(ptype==1) then
          do ll=mm1,mm2
            write(funit,2001) trainy(ll),sisfeat(ll,:nsis(iFCDI))
          end do
        else
          do ll=mm1,mm2
            write(funit,2001) sisfeat(ll,:nsis(iFCDI))  ! no y value for categorical properties
          end do
        end if
        close(funit)
     end do
   elseif(nreaction>0) then
        write(line,'(a,i3.3,a,i3.3,a)') 'space_',iFCDI,'d_p',1,'.dat'
        open(funit,file='feature_space/'//trim(adjustl(line)),status='replace')
        do ll=1,nreaction
          write(funit,2001) trainy(ll),sisfeat(ll,:nsis(iFCDI))
        end do
        close(funit)
   end if

   write(*,'(3a,i10)') 'Size of the SIS-selected subspace from ',phiname,': ',nsis(iFCDI)
   if(nsis(iFCDI)<subs_sis(iFCDI)) then
   write(*,'(a)') 'Warning: The selected subspace is smaller than subs_sis (SISSO.in) !!!'
   end if
   write(9,'(3a,i10)') 'Size of the SIS-selected subspace from ',phiname,': ',nsis(iFCDI)
   if(any(usedup_warning>0) .and. ffdecorr ) &
     write(9,'(/a,i5,a/)') 'Warning: local subspaces used up by SIS. Please increase decorr_alpha and rerun SISSO'
end if  

! release all the spaces
deallocate(trainy)
deallocate(trainy_c)
deallocate(f_select)
deallocate(complexity_select)
deallocate(score_select)
deallocate(name_select)
deallocate(mpin)
deallocate(mpin2)
deallocate(nfpcore_this)
deallocate(nfpcore_tot)
deallocate(available)
deallocate(order)
deallocate(SD)
deallocate(usedup_warning)

if(mpirank==0) then
  deallocate(sisfeat)
  deallocate(name_sisfeat)
end if
deallocate(tag)
deallocate(nvfpcore)
if(lreject) deallocate(reject)

call mpi_barrier(mpi_comm_world,mpierr)
if(mpirank==0) then
   etime_FC=mpi_wtime()
   write(9,'(a,f15.2)') 'Time (second) used for this FC: ',etime_FC-stime_FC
end if

end subroutine


subroutine combine(fin1,name_in1,lastop_in1,compl_in1,dim_in1,fin2,name_in2,lastop_in2,compl_in2,dim_in2,op,nf) 
implicit none
real progress
real*8 fin1(:,:),fin2(:,:),tmp(ubound(fin1,1)),dim_in1(:,:),dim_in2(:,:),&
       dimtmp(ubound(dim_in1,1)),tmp_vf(ubound(fin1,1),vfsize)
integer*8 nfin1,nfin2,i,j,k,kk,kkk,l,nf,compl_in1(:),compl_in2(:),compl_tmp,counter,total_comb
character(len=*) name_in1(:),name_in2(:),op,lastop_in1(:),lastop_in2(:)
character(len=lname) name_tmp,lastop_tmp*10
logical first,skip

nfin1=ubound(fin1,2)  ! from the subspace
nfin2=ubound(fin2,2)  ! from the total
counter=0
progress=0.2
total_comb=mpin(mpirank+1)*(ntot-nthis)+(mpin(mpirank+1)-1)*mpin(mpirank+1)/2.0+&
           mpin(mpirank+1)*(nthis-mpin2(mpirank+1)+1-mpin(mpirank+1))+mpin(mpirank+1)

do i=1,nfin1
  
! no operation
      IF(trim(adjustl(op))=='NO') THEN
          lastop_tmp=''
          compl_tmp=compl_in1(i)
          name_tmp='('//trim(adjustl(name_in1(i)))//')'
          dimtmp=dim_in1(:,i)
        if( isscalar(name_in1(i)) ) then  
             tmp=fin1(:,i)
             call isgoodf(tmp,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
        else
             tmp_vf=vfeat(:,:,nint(fin1(1,i)))
             call isgoodvf(tmp_vf,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
        end if
        cycle
      END IF

! unary operators
! element-wise operation on vector features
      counter=counter+1
      compl_tmp=compl_in1(i)+1
      if(compl_tmp>fcomplexity) goto 599

      ! exp
      IF(index(op,'(exp)')/=0 ) then
      if(  index(lastop_in1(i),'(exp')==0 .and. index(lastop_in1(i),'(log)')==0 ) then ! avoid exp(exp( and exp(log(
        lastop_tmp='(exp)'
        name_tmp='exp('//trim(adjustl(name_in1(i)))//')'
        dimtmp=dimcomb(dim_in1(:,i),dim_in1(:,i),'(exp)')
        if( isscalar(name_in1(i)) ) then  ! vector feature or not
             tmp=exp(fin1(:,i))
             call isgoodf(tmp,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
        else
             tmp_vf=exp(vfeat(:,:,nint(fin1(1,i))))
             call isgoodvf(tmp_vf,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
        end if
      end if
      END IF

      ! exp-
      IF(index(op,'(exp-)')/=0 ) then
      if( index(lastop_in1(i),'(exp')==0  .and. index(lastop_in1(i),'(log)')==0 ) then ! avoid exp(exp( and exp(log(
        lastop_tmp='(exp-)'
        name_tmp='exp(-'//trim(adjustl(name_in1(i)))//')'
        dimtmp=dimcomb(dim_in1(:,i),dim_in1(:,i),'(exp-)')
        if( isscalar(name_in1(i)) ) then  ! vector feature or not
             tmp=exp(-fin1(:,i))
             call isgoodf(tmp,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
        else
             tmp_vf=exp(-vfeat(:,:,nint(fin1(1,i))))
             call isgoodvf(tmp_vf,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
        end if
      end if
      END IF

      ! ^-1
      IF(index(op,'(^-1)')/=0) then
        if(minval(abs(fin1(:,i)))>1d-50 ) then
           lastop_tmp='(^-1)'
           name_tmp='('//trim(adjustl(name_in1(i)))//')^-1'
           dimtmp=dimcomb(dim_in1(:,i),dim_in1(:,i),'(^-1)')
           if( isscalar(name_in1(i)) ) then  ! vector feature or not
                tmp=(fin1(:,i))**(-1)
                call isgoodf(tmp,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
           else
                tmp_vf=(vfeat(:,:,nint(fin1(1,i))))**(-1)
                call isgoodvf(tmp_vf,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
           end if
        end if
       END IF

      ! scd: Standard Cauchy Distribution
      IF(index(op,'(scd)')/=0) then
           lastop_tmp='(scd)'
           name_tmp='scd('//trim(adjustl(name_in1(i)))//')'
           dimtmp=dimcomb(dim_in1(:,i),dim_in1(:,i),'(scd)')
           if( isscalar(name_in1(i)) ) then  ! vector feature or not
                tmp=1.0d0/(PI*(1.0d0+(fin1(:,i))**2))
                call isgoodf(tmp,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
           else
                tmp_vf=1.0d0/(PI*(1.0d0+(vfeat(:,:,nint(fin1(1,i))))**2))
                call isgoodvf(tmp_vf,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
           end if
       END IF

      ! ^2
      IF(index(op,'(^2)')/=0) then
        if(index(lastop_in1(i),'(sqrt)')==0 ) then ! avoid (sqrt())^2
           lastop_tmp='(^2)'
           name_tmp='('//trim(adjustl(name_in1(i)))//')^2'
           dimtmp=dimcomb(dim_in1(:,i),dim_in1(:,i),'(^2)')
           if( isscalar(name_in1(i)) ) then  ! vector feature or not
                tmp=(fin1(:,i))**2
                call isgoodf(tmp,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
           else
                tmp_vf=(vfeat(:,:,nint(fin1(1,i))))**2
                call isgoodvf(tmp_vf,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
           end if
        end if
       END IF

      ! ^3
      IF(index(op,'(^3)')/=0) then
        if(index(lastop_in1(i),'(cbrt)')==0 ) then ! avoid (cbrt())^3
         lastop_tmp='(^3)'
         name_tmp='('//trim(adjustl(name_in1(i)))//')^3'
         dimtmp=dimcomb(dim_in1(:,i),dim_in1(:,i),'(^3)')
         if(  isscalar(name_in1(i)) ) then  ! vector feature or not
              tmp=(fin1(:,i))**3
              call isgoodf(tmp,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
         else
              tmp_vf=(vfeat(:,:,nint(fin1(1,i))))**3
              call isgoodvf(tmp_vf,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
         end if
       end if
     END IF

      ! ^6
      IF(index(op,'(^6)')/=0) then
         lastop_tmp='(^6)'
         name_tmp='('//trim(adjustl(name_in1(i)))//')^6'
         dimtmp=dimcomb(dim_in1(:,i),dim_in1(:,i),'(^6)')
         if( isscalar(name_in1(i)) ) then  ! vector feature or not
              tmp=(fin1(:,i))**6
              call isgoodf(tmp,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
         else
              tmp_vf=(vfeat(:,:,nint(fin1(1,i))))**6
              call isgoodvf(tmp_vf,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
         end if
       END IF

      ! sqrt
      IF(index(op,'(sqrt)')/=0) then
        if(index(lastop_in1(i),'(^2)')==0 ) then  ! avoid sqrt((^2))
          if( minval(fin1(:,i))>0 .and. (isscalar(name_in1(i)))  ) then
              lastop_tmp='(sqrt)'     
              name_tmp='sqrt('//trim(adjustl(name_in1(i)))//')'
              dimtmp=dimcomb(dim_in1(:,i),dim_in1(:,i),'(sqrt)')
              tmp=sqrt(fin1(:,i))
              call isgoodf(tmp,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
          else if ( .not. isscalar(name_in1(i)) )  then
              if(minval(vfeat(:,:,nint(fin1(1,i))))>0) then
              lastop_tmp='(sqrt)'
              name_tmp='sqrt('//trim(adjustl(name_in1(i)))//')'
              dimtmp=dimcomb(dim_in1(:,i),dim_in1(:,i),'(sqrt)')
              tmp_vf=sqrt(vfeat(:,:,nint(fin1(1,i))))
              call isgoodvf(tmp_vf,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
              end if
          end if
         end if
       END IF

      ! cbrt: cube root
      IF(index(op,'(cbrt)')/=0) then
        if(index(lastop_in1(i),'(^3)')==0 ) then  ! avoid cbrt((^3))
          if( isscalar(name_in1(i)) ) then
              lastop_tmp='(cbrt)'
              name_tmp='cbrt('//trim(adjustl(name_in1(i)))//')'
              dimtmp=dimcomb(dim_in1(:,i),dim_in1(:,i),'(cbrt)')
              tmp=(fin1(:,i))**(1.d0/3.d0)
              call isgoodf(tmp,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
          else if ( .not. isscalar(name_in1(i)) )  then
              lastop_tmp='(cbrt)'
              name_tmp='cbrt('//trim(adjustl(name_in1(i)))//')'
              dimtmp=dimcomb(dim_in1(:,i),dim_in1(:,i),'(cbrt)')
              tmp_vf=(vfeat(:,:,nint(fin1(1,i))))**(1.d0/3.d0)
              call isgoodvf(tmp_vf,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
          end if
         end if
       END IF

      ! log
      IF(index(op,'(log)')/=0) then
        if( index(lastop_in1(i),'(exp')==0 .and. index(lastop_in1(i),'(log)')==0  ) then ! avoid log(exp( and log(log(
          if( minval(fin1(:,i))>0 .and. (isscalar(name_in1(i)))  ) then
              lastop_tmp='(log)'
              name_tmp='log('//trim(adjustl(name_in1(i)))//')'
              dimtmp=dimcomb(dim_in1(:,i),dim_in1(:,i),'(log)')
              tmp=log(fin1(:,i))
              call isgoodf(tmp,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
          else if ( .not. isscalar(name_in1(i)) )  then
              if(minval(vfeat(:,:,nint(fin1(1,i))))>0) then
              lastop_tmp='(log)'
              name_tmp='log('//trim(adjustl(name_in1(i)))//')'
              dimtmp=dimcomb(dim_in1(:,i),dim_in1(:,i),'(log)')
              tmp_vf=log(vfeat(:,:,nint(fin1(1,i))))
              call isgoodvf(tmp_vf,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
              end if
          end if
         end if
       END IF

      ! sin
      IF(index(op,'(sin)')/=0) then
         lastop_tmp='(sin)'
         name_tmp='sin('//trim(adjustl(name_in1(i)))//')'
         dimtmp=dimcomb(dim_in1(:,i),dim_in1(:,i),'(sin)')
         if(  isscalar(name_in1(i)) ) then  ! vector feature or not
              tmp=sin(fin1(:,i))
              call isgoodf(tmp,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
         else
              tmp_vf=sin(vfeat(:,:,nint(fin1(1,i))))
              call isgoodvf(tmp_vf,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
         end if
     END IF

      ! cos
      IF(index(op,'(cos)')/=0) then
         lastop_tmp='(cos)'
         name_tmp='cos('//trim(adjustl(name_in1(i)))//')'
         dimtmp=dimcomb(dim_in1(:,i),dim_in1(:,i),'(cos)')
         if(  isscalar(name_in1(i)) ) then  ! vector feature or not
              tmp=cos(fin1(:,i))
              call isgoodf(tmp,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
         else
              tmp_vf=cos(vfeat(:,:,nint(fin1(1,i))))
              call isgoodvf(tmp_vf,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
         end if
     END IF

599  continue

  ! binary operators
  ! element-wise operations between vector features
  do j=1,nfin2

      if( j-(ntot-nthis)>0 .and. j-(ntot-nthis)<=i ) cycle

      counter=counter+1
      compl_tmp=compl_in1(i)+compl_in2(j)+1
      if(compl_tmp>fcomplexity) goto 602

      if(simpler(compl_in1(i),compl_in2(j),name_in1(i),name_in2(j))==1) then
         first=.true.
      else
         first=.false.
      end if

      ! sum and subtract 
       IF(index(op,'(+)')/=0 .or. index(op,'(-)')/=0  .or. index(op,'(|-|)')/=0   ) THEN
               ! linear combination for exp() and log()
               if( (index(lastop_in1(i),'(exp')==0 .and. index(lastop_in2(j),'(exp')/=0)    .or. &
                    (index(lastop_in1(i),'(exp')/=0 .and. index(lastop_in2(j),'(exp')==0)    .or. &
                    (index(lastop_in1(i),'(log)')==0 .and. index(lastop_in2(j),'(log)')/=0)  .or. &
                    (index(lastop_in1(i),'(log)')/=0 .and. index(lastop_in2(j),'(log)')==0) ) goto 600
               ! homogeneous combination
               if ( maxval(abs(dim_in1(:,i)-dim_in2(:,j)))>1d-8 ) goto 600

               !---
               IF(index(op,'(+)')/=0) then
                   lastop_tmp='(+)'
                   if(first) then
                      name_tmp='('//trim(adjustl(name_in1(i)))//'+'//trim(adjustl(name_in2(j)))//')'
                   else
                      name_tmp='('//trim(adjustl(name_in2(j)))//'+'//trim(adjustl(name_in1(i)))//')'
                   end if
                   dimtmp=dimcomb(dim_in1(:,i),dim_in2(:,j),'(+)')  
                  if( isscalar(name_in1(i)) .and. isscalar(name_in2(j))  ) then 
                       tmp=fin1(:,i)+fin2(:,j)
                       call isgoodf(tmp,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
                  else 
                   ! combination between a vector and a scalar feature
                   !  if(  .not. isscalar(name_in1(i))  .and. isscalar(name_in2(j))  )then 
                   !    tmp=sum(vfeat(:,:,nint(fin1(1,i))),2)+fin2(:,j)
                   !    call isgoodf(tmp,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
                   !  else if( isscalar(name_in1(i)) .and. (.not. isscalar(name_in2(j))) )then
                   !    tmp=fin1(:,i)+sum(vfeat(:,:,nint(fin2(1,j))),2)
                   !    call isgoodf(tmp,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
                     if (  (.not. isscalar(name_in1(i)))  .and. (.not. isscalar(name_in2(j))) ) then
                        tmp_vf=vfeat(:,:,nint(fin1(1,i)))+vfeat(:,:,nint(fin2(1,j)))
                        call isgoodvf(tmp_vf,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
                      end if
                  end if
                END IF
                !---

                IF(index(op,'(-)')/=0) then  ! A-B
                  lastop_tmp='(-)'
                  name_tmp='('//trim(adjustl(name_in1(i)))//'-'//trim(adjustl(name_in2(j)))//')'
                  if(index(op,'(+)')==0) dimtmp=dimcomb(dim_in1(:,i),dim_in2(:,j),'(-)')   
                  if( isscalar(name_in1(i)) .and. isscalar(name_in2(j))  ) then 
                       tmp=fin1(:,i)-fin2(:,j)
                       call isgoodf(tmp,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
                  else 
                     ! if( (.not. isscalar(name_in1(i)))  .and. isscalar(name_in2(j)) ) then
                     !     tmp=sum(vfeat(:,:,nint(fin1(1,i))),2)-fin2(:,j)
                     !     call isgoodf(tmp,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
                     ! else if( isscalar(name_in1(i)) .and. (.not. isscalar(name_in2(j))) )then
                     !     tmp=fin1(:,i)-sum(vfeat(:,:,nint(fin2(1,j))),2)
                     !     call isgoodf(tmp,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
                      if (  (.not. isscalar(name_in1(i)))  .and. (.not. isscalar(name_in2(j)))  ) then
                           tmp_vf=vfeat(:,:,nint(fin1(1,i)))-vfeat(:,:,nint(fin2(1,j)))
                          call isgoodvf(tmp_vf,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
                      end if
                   end if

                  if(icomb<rung) then  ! B-A (no need for the last combination)
                    ! lastop_tmp='(-)'   ! same with above
                     name_tmp='('//trim(adjustl(name_in2(j)))//'-'//trim(adjustl(name_in1(i)))//')'
                    ! dimtmp=dimcomb(dim_in1(:,i),dim_in2(:,j),'(-)')   ! same with above
                     if( isscalar(name_in2(j)) .and. isscalar(name_in1(i))  ) then
                          tmp=-tmp
                          call isgoodf(tmp,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
                     else
                      !if( (.not. isscalar(name_in1(i)))  .and. isscalar(name_in2(j)) ) then
                      !    tmp=-tmp
                      !    call isgoodf(tmp,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
                      !else if( isscalar(name_in1(i)) .and. (.not. isscalar(name_in2(j))) )then
                      !    tmp=-tmp
                      !    call isgoodf(tmp,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
                      if (  (.not. isscalar(name_in1(i)))  .and. (.not. isscalar(name_in2(j)))  ) then
                          tmp_vf=-tmp_vf
                          call isgoodvf(tmp_vf,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
                      end if
                   end if
                  end if
                 END IF

                !-----
                IF(index(op,'(|-|)')/=0) then  ! absolute difference
                  lastop_tmp='(|-|)'
                  if(first) then
                    name_tmp='abs('//trim(adjustl(name_in1(i)))//'-'//trim(adjustl(name_in2(j)))//')'
                  else
                    name_tmp='abs('//trim(adjustl(name_in2(j)))//'-'//trim(adjustl(name_in1(i)))//')'
                  end if
                  if(index(op,'(+)')==0) dimtmp=dimcomb(dim_in1(:,i),dim_in2(:,j),'(|-|)')  
                  if( isscalar(name_in1(i)) .and. isscalar(name_in2(j))  ) then
                       if(index(op,'(-)')/=0) then
                         tmp=abs(tmp)   ! if (-) exists, save time
                       else
                         tmp=abs(fin1(:,i)-fin2(:,j))
                       end if
                       call isgoodf(tmp,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
                  else
                     ! if( (.not. isscalar(name_in1(i)))  .and. isscalar(name_in2(j)) ) then
                     !    if(index(op,'(-)')/=0) then
                     !      tmp=abs(tmp)   ! if (-) exists, save time
                     !    else
                     !     tmp=abs(sum(vfeat(:,:,nint(fin1(1,i))),2)-fin2(:,j))
                     !   end if
                     !     call isgoodf(tmp,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
                     ! else if( isscalar(name_in1(i)) .and. (.not. isscalar(name_in2(j))) )then
                     !    if(index(op,'(-)')/=0) then
                     !      tmp=abs(tmp)   ! if (-) exists, save time
                     !    else
                     !     tmp=abs(fin1(:,i)-sum(vfeat(:,:,nint(fin2(1,j))),2))
                     !    end if
                     !     call isgoodf(tmp,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
                       if (  (.not. isscalar(name_in1(i)))  .and. (.not. isscalar(name_in2(j)))  ) then
                         if(index(op,'(-)')/=0) then
                           tmp_vf=abs(tmp_vf)   ! if (-) exists, save time
                         else
                           tmp_vf=abs(vfeat(:,:,nint(fin1(1,i)))-vfeat(:,:,nint(fin2(1,j))))
                         end if
                          call isgoodvf(tmp_vf,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
                       end if
                  end if
                 END IF
                 !----
         END IF

600   continue

      ! multiplication and division
           ! exp(A)*exp(B)=exp(A+B)
           if( (index(lastop_in1(i),'(exp)')/=0 .and. index(lastop_in2(j),'(exp)')/=0) .or. &
               (index(lastop_in1(i),'(exp-)')/=0 .and. index(lastop_in2(j),'(exp-)')/=0) .or. &
               (index(lastop_in1(i),'(log)')/=0 .and. index(lastop_in2(j),'(log)')/=0) ) goto 602

           !multiplication
            IF(index(op,'(*)')/=0) then
                lastop_tmp='(*)'
               if(first) then
                  name_tmp='('//trim(adjustl(name_in1(i)))//'*'//trim(adjustl(name_in2(j)))//')'
               else
                  name_tmp='('//trim(adjustl(name_in2(j)))//'*'//trim(adjustl(name_in1(i)))//')'
               end if            
                dimtmp=dimcomb(dim_in1(:,i),dim_in2(:,j),'(*)')
               if( isscalar(name_in1(i)) .and. isscalar(name_in2(j)) ) then 
                    tmp=fin1(:,i)*fin2(:,j)
                    call isgoodf(tmp,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
               else 
                  !  if( (.not. isscalar(name_in1(i)))  .and. isscalar(name_in2(j))  ) then
                  !     do k=1,vfsize
                  !        tmp_vf(:,k)=vfeat(:,k,nint(fin1(1,i)))*fin2(:,j)
                  !     end do
                  !  else if(isscalar(name_in1(i)) .and. (.not. isscalar(name_in2(j))) ) then
                  !     do k=1,vfsize
                  !        tmp_vf(:,k)=fin1(:,i)*vfeat(:,k,nint(fin2(1,j)))
                  !     end do
                    if (  (.not. isscalar(name_in1(i)))  .and. (.not. isscalar(name_in2(j))) ) then
                       tmp_vf=vfeat(:,:,nint(fin1(1,i)))*vfeat(:,:,nint(fin2(1,j)))
                    end if
                    call isgoodvf(tmp_vf,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
               end if
            END IF

          ! division
           IF(index(op,'(/)')/=0) then

              ! avoid repeatition: A/(B/C)=(A*C)/B, (A/B)/C=A/(B*C), A*B/B=A, A/(A*B)=B^-1
              if(index(op,'(*)')/=0 .and. (index(lastop_in1(i),'(/)')/=0 .or. index(lastop_in2(j),'(/)')/=0) .and. &
                 icomb<rung ) goto 602

              skip=.false.
              if(index(lastop_in1(i),'(*)')/=0) then
                 l=index(adjustl(name_in1(i)),trim(adjustl(name_in2(j))))
                 if(l>0) then
                   k=index(name_in1(i),trim(adjustl(name_in2(j))))
                   kk=len_trim(name_in1(i)(k:))
                   kkk=len_trim(adjustl(name_in2(j)))
                   if( (l==2 .and. name_in1(i)(k+kkk:k+kkk)=='*') .or. (kk==kkk+1 .and. name_in1(i)(k-1:k-1)=='*')) then
                      if(index(op,'(^-1)')/=0) then
                          goto 602
                      else
                          goto 601
                      end if
                   end if
                 end if
              end if
                 
             if(index(lastop_in1(i),'(*)')/=0) then
                l=index(adjustl(name_in2(j)),trim(adjustl(name_in1(i))))
                if(l>0) then
                  k=index(name_in2(j),trim(adjustl(name_in1(i))))
                  kk=len_trim(name_in2(j)(k:))
                  kkk=len_trim(adjustl(name_in1(i)))
                  if((l==2 .and. name_in2(j)(k+kkk:k+kkk)=='*') .or. (kk==kkk+1 .and. name_in2(j)(k-1:k-1)=='*')) then
                      if(index(op,'(^-1)')/=0) then
                          goto 602
                      else
                          skip=.true.
                      end if
                  end if
                end if
             end if

              !i/j
              !---------
              lastop_tmp='(/)'
              name_tmp='('//trim(adjustl(name_in1(i)))//'/'//trim(adjustl(name_in2(j)))//')'
              dimtmp=dimcomb(dim_in1(:,i),dim_in2(:,j),'(/)')
               if( isscalar(name_in1(i)) .and. isscalar(name_in2(j))  ) then
                    if(minval(abs(fin2(:,j)))>1d-50 ) then 
                    tmp=fin1(:,i)/fin2(:,j)
                    call isgoodf(tmp,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
                    end if
               else 
                 !   if( (.not. isscalar(name_in1(i)))  .and. isscalar(name_in2(j)) ) then
                 !      if(minval(abs(fin2(:,j)))>1d-50) then
                 !        do k=1,vfsize
                 !           tmp_vf(:,k)=vfeat(:,k,nint(fin1(1,i)))/fin2(:,j)
                 !        end do
                 !      end if
                 ! else if(isscalar(name_in1(i)) .and. (.not. isscalar(name_in2(j))) )then
                 !      if(minval(abs(vfeat(:,:,nint(fin2(1,j)))))>1d-50) then
                 !        do k=1,vfsize
                 !          tmp_vf(:,k)=fin1(:,i)/vfeat(:,k,nint(fin2(1,j)))
                 !        end do
                 !      end if
                  if (  (.not. isscalar(name_in1(i)))  .and. (.not. isscalar(name_in2(j)))  ) then
                       if(minval(abs(vfeat(:,:,nint(fin2(1,j)))))>1d-50 ) then
                       tmp_vf=vfeat(:,:,nint(fin1(1,i)))/vfeat(:,:,nint(fin2(1,j)))
                       end if
                  end if
                  call isgoodvf(tmp_vf,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
              end if
              !------

     601      continue
              if(skip) goto 602

            !j/i
            !---------
                 lastop_tmp='(/)'
                 name_tmp='('//trim(adjustl(name_in2(j)))//'/'//trim(adjustl(name_in1(i)))//')'
                 dimtmp=dimcomb(dim_in2(:,j),dim_in1(:,i),'(/)')
                if( isscalar(name_in1(i)) .and.  isscalar(name_in2(j)) ) then
                     if(minval(abs(fin1(:,i)))>1d-50) then
                     tmp=fin2(:,j)/fin1(:,i)
                     call isgoodf(tmp,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
                     end if
                else
                  !   if( (.not. isscalar(name_in1(i)))  .and.  isscalar(name_in2(j)) ) then
                  !      if(minval(abs(vfeat(:,:,nint(fin1(1,i)))))>1d-50) then
                  !        do k=1,vfsize
                  !          tmp_vf(:,k)=fin2(:,j)/vfeat(:,k,nint(fin1(1,i)))
                  !        end do
                  !      end if
                  !   elseif( isscalar(name_in1(i)) .and. (.not. isscalar(name_in2(j))) )then
                  !      if(minval(abs(fin1(:,i)))>1d-50) then
                  !         do k=1,vfsize
                  !           tmp_vf(:,k)=vfeat(:,k,nint(fin2(1,j)))/fin1(:,i)
                  !         end do
                  !      end if
                     if ( (.not. isscalar(name_in1(i)))  .and. (.not. isscalar(name_in2(j)))  ) then
                        if(minval(abs(vfeat(:,:,nint(fin1(1,i)))))>1d-50 ) then
                          tmp_vf=vfeat(:,:,nint(fin2(1,j)))/vfeat(:,:,nint(fin1(1,i)))
                        end if
                     end if
                     call isgoodvf(tmp_vf,name_tmp,lastop_tmp,compl_tmp,dimtmp,nf)
               end if
           END if

   602   continue

  end do

  if(float(counter)/float(total_comb)>=progress) then
  write(*,'(a,i4,2(a,i15),a,f6.1,a)') &
  'mpirank = ',mpirank,'  #generated =',nf,'  #selected =',nselect,'  progress =',progress*100,'%'
  progress=progress+0.2
  end if

end do

end subroutine


function goodf(feat,name_feat,dimens,compl)
integer*8 i,j,k,nf,l,ll,mm1,mm2,compl
real*8 feat(:),dimens(:),scoretmp(2),maxabs,react_feat(npoints)
character(len=*) name_feat
logical goodf,lsame

goodf=.true.

! delete constant/zero/infinity features 
!do ll=1,ntask
!  mm1=sum(nsample(:ll-1))+1
!  mm2=sum(nsample(:ll))
   mm1=1
   mm2=sum(nsample(:ntask))
  if(maxval(abs(feat(mm1:mm2)-feat(mm1)))<=1d-8) then
    goodf=.false. ! constant feature.
    return
  end if
  maxabs=maxval(abs(feat(mm1:mm2)))
  if(maxabs>1d50 .or. maxabs<=1d-50) then
    goodf=.false. ! intinity or zero
    return
  end if
!end do

! not to be selected but can be used for further transformation
if(maxabs>maxfval_ub .or. maxabs<maxfval_lb) return

if(nreaction==0) then
  ! calculate the score for each feature
  if(ptype==1) scoretmp=sis_score(feat,trainy_c)
  if(ptype==2) scoretmp=sis_score(feat,trainy) !classification
else if(nreaction>0) then
  ! create reaction-feature
  do i=1,nreaction
     react_feat(i)=sum(react_coeff(i,:)*feat([react_speciesID(i,:)]))
  end do
  scoretmp=sis_score(react_feat,trainy_c)
end if

if(threshold_select) then
   if (scoretmp(1)<score_select(nfinal_select,1)) return
end if

!----
if(nreject>0) then
  name_feat=adjustl(name_feat)
  lsame=.false.
  i=0; j=nreject; k=i+ceiling((j-i)/2.0)
  do while(i/=j)
     if(trim(name_feat)==trim(reject(k))) then
       lsame=.true.
       i=j
     else if(name_feat<reject(k)) then
       i=k;
     else if(name_feat>reject(k)) then
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
nselect=nselect+1
if(nreaction==0) then
  f_select(:,nselect)=feat
elseif(nreaction>0) then
  f_select(:,nselect)=react_feat
endif
complexity_select(nselect)=compl
name_select(nselect)=name_feat
score_select(nselect,:)=scoretmp
if( nselect== nbasic_select+nextra_select )then
  call update_select
  threshold_select=.true.
end if

end function


subroutine addm_out(n)
! increase array size
real*8,allocatable:: real2d(:,:)
integer*8,allocatable:: integer1d(:)
character(len=lname),allocatable:: char1d(:),char1d2(:)*10
integer*8 i,j,k,n
i=ubound(fout,1)
j=ubound(fout,2)

! fout
allocate(real2d(i,j))
real2d=fout
deallocate(fout)
allocate(fout(i,j+n))
fout(:,:j)=real2d
deallocate(real2d)
!-----
!name_out
allocate(char1d(j))
char1d=name_out
deallocate(name_out)
allocate(name_out(j+n))
name_out(:j)=char1d
deallocate(char1d)
!----
!lastop_out
allocate(char1d2(j))
char1d2=lastop_out
deallocate(lastop_out)
allocate(lastop_out(j+n))
lastop_out(:j)=char1d2
deallocate(char1d2)
!----
!complexity_out
allocate(integer1d(j))
integer1d=complexity_out
deallocate(complexity_out)
allocate(complexity_out(j+n))
complexity_out(:j)=integer1d
deallocate(integer1d)
!dim_out
i=ubound(dim_out,1)
allocate(real2d(i,j))
real2d=dim_out
deallocate(dim_out)
allocate(dim_out(i,j+n))
dim_out(:,:j)=real2d
deallocate(real2d)
!--
end subroutine


subroutine addm_in1(n,fin1,name_in1,lastop_in1,complexity_in1,dim_in1)
! increase array size
real*8,allocatable:: real2d(:,:),fin1(:,:),dim_in1(:,:)
integer*8,allocatable:: integer1d(:),complexity_in1(:)
character(len=lname),allocatable:: char1d(:),char1d2(:)*10
character(len=*),allocatable:: name_in1(:),lastop_in1(:)
integer*8 i,j,k,n
i=ubound(fin1,1)
j=ubound(fin1,2)

! fin1
allocate(real2d(i,j))
real2d=fin1
deallocate(fin1)
allocate(fin1(i,j+n))
fin1(:,:j)=real2d
deallocate(real2d)
!-----
!name_in1
allocate(char1d(j))
char1d=name_in1
deallocate(name_in1)
allocate(name_in1(j+n))
name_in1(:j)=char1d
deallocate(char1d)
!----
!lastop_in1
allocate(char1d2(j))
char1d2=lastop_in1
deallocate(lastop_in1)
allocate(lastop_in1(j+n))
lastop_in1(:j)=char1d2
deallocate(char1d2)
!---
!complexity_in1
allocate(integer1d(j))
integer1d=complexity_in1
deallocate(complexity_in1)
allocate(complexity_in1(j+n))
complexity_in1(:j)=integer1d
deallocate(integer1d)
!dim_in1
i=ubound(dim_in1,1)
allocate(real2d(i,j))
real2d=dim_in1
deallocate(dim_in1)
allocate(dim_in1(i,j+n))
dim_in1(:,:j)=real2d
deallocate(real2d)
!--
end subroutine

subroutine addm_vf(n)
! increase array size for vector feature
real*8,allocatable:: real3d(:,:,:)
integer*8 i,j,k,n
i=ubound(vfeat,1)
j=ubound(vfeat,2)
k=ubound(vfeat,3)
allocate(real3d(i,j,k))
real3d=vfeat
deallocate(vfeat)
allocate(vfeat(i,j,k+n))
vfeat(:,:,:k)=real3d
deallocate(real3d)
end subroutine


function dimcomb(dim1,dim2,op)
! calculate the unit for the new feature
! unary operator: set dim1 and dim2 the same
real*8 dim1(:),dim2(:),dimcomb(ubound(dim1,1))
character(len=*) op
integer i,j,k

if(trim(adjustl(op))=='(+)' .or. trim(adjustl(op))=='(-)' .or. trim(adjustl(op))=='(|-|)'  ) then
   do i=1,ubound(dim1,1)
      if(abs(dim1(i)-dim2(i))>1d-8) stop 'Error: dimension not matched'
   end do
   dimcomb=dim1
else if(trim(adjustl(op))=='(*)') then
   dimcomb=dim1+dim2
else if(trim(adjustl(op))=='(/)') then
   dimcomb=dim1-dim2
else if(trim(adjustl(op))=='(exp)') then   
   dimcomb=0.d0
else if(trim(adjustl(op))=='(exp-)') then
   dimcomb=0.d0
else if(trim(adjustl(op))=='(log)') then
   dimcomb=0.d0
else if(trim(adjustl(op))=='(scd)') then 
   dimcomb=0.d0
else if(trim(adjustl(op))=='(sin)') then
   dimcomb=0.d0
else if(trim(adjustl(op))=='(cos)') then
   dimcomb=0.d0
else if(trim(adjustl(op))=='(^-1)') then
   dimcomb=dim1*(-1)
else if(trim(adjustl(op))=='(^2)') then
   dimcomb=dim1*2
else if(trim(adjustl(op))=='(^3)') then
   dimcomb=dim1*3
else if(trim(adjustl(op))=='(^6)') then
   dimcomb=dim1*6
else if(trim(adjustl(op))=='(sqrt)') then
   dimcomb=dim1/2.d0
else if(trim(adjustl(op))=='(cbrt)') then
   dimcomb=dim1/3.d0
end if
end function


subroutine writeout(phiname,i,j,k)
integer*8 i,j,k
character(len=*) phiname
if(nvf>0) write(*,'(a,i15)') 'Total Number vector features: ',j
!write(*,'(a,i15)') 'Total mumber of selected features by SIS: ',i
write(*,'(3a,i15)') 'Total number of features in the space ',trim(phiname),':',k
write(9,'(3a,i15)') 'Total number of features in the space ',trim(phiname),':',k
end subroutine


function simpler(complexity1,complexity2,name1,name2)
integer*8 simpler,complexity1,complexity2
character(len=*) name1,name2

if(complexity1<complexity2) then
   simpler=1
else if(complexity1>complexity2) then
   simpler=2
else if(complexity1==complexity2) then
! the names make no sense, but to make sure a unique return when two features are exactly the same except names. 
   if(len_trim(name1)<len_trim(name2) .or. &
     (len_trim(name1)==len_trim(name2) .and. trim(name1)<=trim(name2) ) ) then
      simpler=1
   else
      simpler=2
   end if
end if
end function


subroutine update_availability(available,pavailable)
logical available(:),pavailable(:)
if( any(available) ) then
   pavailable(mpirank+1)=.true.
else
   pavailable(mpirank+1)=.false.
   usedup_warning(mpirank+1)=usedup_warning(mpirank+1)+1
   if(usedup_warning(mpirank+1)==1 .and. ffdecorr ) then
     write(*,'(/a,i5,a)') 'Warning: the local subspace at mpirank ',mpirank,' is used up by SIS. &
                    Please increase decorr_alpha and rerun SISSO'
   end if
end if
end subroutine


subroutine dup_pcheck(nfpcore_this,fidentity,fname,complexity,order,available,ftype)
real*8 fidentity(:,:),time_start,time_end
integer*8 i,j,k,l,ll,mpii,mpij,nfpcore_this(:),mpim(mpisize),mpiloc(mpisize),order(:),loc(1),&
          complexity(:),simpler_result
character(len=*) ftype,fname(:)
logical available(:)
real*8,allocatable:: compidentity(:,:),compfeat(:,:)
character(len=lname),allocatable:: compname(:)
integer*8,allocatable:: compcomplexity(:)

IF(mpisize>1) THEN
   mpim=nfpcore_this  ! number of features on each core
   do i=1,mpisize 
      loc(1:1)=maxloc(mpim)
      mpiloc(i)=loc(1)  ! ordering the cores with the number of features from high to low
      mpim(loc(1))=-1
   end do

   do k=1,mpisize-1

     ! time recording
!     if(mpirank==0) time_start=mpi_wtime()
     if(nfpcore_this(mpiloc(k))>0) then
         ! allocate arrays
         allocate(compidentity(nfpcore_this(mpiloc(k)),2))
         allocate(compname(nfpcore_this(mpiloc(k))))
         allocate(compcomplexity(nfpcore_this(mpiloc(k))))
         allocate(compfeat(npoints,nfpcore_this(mpiloc(k))))

         if(mpirank==mpiloc(k)-1) then
           compidentity=fidentity(:nfpcore_this(mpiloc(k)),:)
           compname=fname(:nfpcore_this(mpiloc(k)))
           compcomplexity=complexity(:nfpcore_this(mpiloc(k)))
           if(index(ftype,'selected')/=0) compfeat=f_select(:,:nfpcore_this(mpiloc(k)))
         end if

         !broadcast compidentity from core mpiloc(k)-1 to all other cores
          call mpi_bcast(compidentity,2*nfpcore_this(mpiloc(k)),mpi_double_precision,mpiloc(k)-1,mpi_comm_world,mpierr)
         ! broadcast the feature names
          call mpi_bcast(compname,nfpcore_this(mpiloc(k))*lname,mpi_character,mpiloc(k)-1,mpi_comm_world,mpierr)
         ! broadcast the feature complexities
          call mpi_bcast(compcomplexity,nfpcore_this(mpiloc(k)),mpi_integer8,mpiloc(k)-1,mpi_comm_world,mpierr)
         ! broadcast the selected features
         if (index(ftype,'selected')/=0) call mpi_bcast(compfeat,npoints*nfpcore_this(mpiloc(k)),&
                                                        mpi_double_precision,mpiloc(k)-1,mpi_comm_world,mpierr)

         ! do the comparision. All cores of mpiloc(k:mpisize) will make comparison with
         ! that on core mpiloc(k)-1
         if( all( mpirank/=(mpiloc(:k)-1) ) .and. nfpcore_this(mpirank+1)>0 ) then 
            do mpij=1,nfpcore_this(mpiloc(k))
              l=0; ll=order(nfpcore_this(mpirank+1)+1);
              i=l+ceiling(float(ll-l)/2.0)   ! bisection method

              124 continue
              simpler_result=simpler(compcomplexity(mpij),complexity(order(i)),compname(mpij),fname(order(i)))

              ! the check for 'selected' is more strict than for 'all'.
              if(index(ftype,'selected')/=0) then
                 ! get rid of duplication and correlation for selected features
                 if(equivalent(compidentity(mpij,:),fidentity(order(i),:),compfeat(:npoints,mpij),& ! mpij from mpiloc(k)-1
                               f_select(:npoints,order(i))) ) then ! order(i) mpirank
                       if(simpler_result==1) then
                          available(order(i))=.false.
                       else
                          compidentity(mpij,1)=0.d0
                       end if
                    cycle
                 end if
                 ! if not equal
                 if(compidentity(mpij,1)>fidentity(order(i),1) .or. (compidentity(mpij,1)==fidentity(order(i),1) &
                    .and. compidentity(mpij,2)>fidentity(order(i),2)) .or. &
                   (compidentity(mpij,1)==fidentity(order(i),1) .and. &
                    compidentity(mpij,2)==fidentity(order(i),2) .and. simpler_result==1 ) ) then
                    ll=i
                 else
                     l=i
                 end if
                 if(i==l+ceiling(float(ll-l)/2.0)) cycle
                 i=l+ceiling(float(ll-l)/2.0)
                 goto 124
              else if(index(ftype,'all')/=0) then
                 if(abs(compidentity(mpij,1)-fidentity(order(i),1))<=1d-8 ) then
                 ! if equal
                    if(simpler_result==1) then
                       available(order(i))=.false.
                    else
                       compidentity(mpij,1)=0.d0
                    end if
                    cycle
                 else   ! if not equal
                   if(compidentity(mpij,1)>fidentity(order(i),1)) then
                      ll=i
                   else
                       l=i
                   end if
                   if(i==l+ceiling(float(ll-l)/2.0)) cycle
                    i=l+ceiling(float(ll-l)/2.0)
                    goto 124
                 end if
              end if
            end do
        call mpi_send(compidentity,2*nfpcore_this(mpiloc(k)),mpi_double_precision,mpiloc(k)-1,33,&
                      mpi_comm_world,status,mpierr)
          else if(mpirank==mpiloc(k)-1) then
            do l=1,mpisize
             if( any( l==mpiloc(:k) ) .or. nfpcore_this(l)==0 ) cycle
        call mpi_recv(compidentity,2*nfpcore_this(mpiloc(k)),mpi_double_precision,mpi_any_source,33,&
                      mpi_comm_world,status,mpierr)
               do i=1,nfpcore_this(mpiloc(k))
                if(abs(compidentity(i,1))<1d-8) available(i)=.false.
               end do
            end do
         end if

         ! deallocate
         deallocate(compidentity)
         deallocate(compname)
         deallocate(compcomplexity)
         deallocate(compfeat)
     end if

     ! if(mpirank==0) then
     ! time_end=mpi_wtime()
     ! write(*,'(a,i5,a,i4.4,a,f15.2)') &
     !  'Time (s) spent on block ',k,'/',mpisize-1,'   is:',(time_end-time_start)
     ! end if
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
!   write(*,'(a,i10)') 'Number of features after the redundant check is: ',j
  end if

END IF

end subroutine


subroutine sure_indep_screening(nfpcore_this,available,complexity,fname,feat,sisfeat,name_sisfeat)
! sure independence screening
real*8 tmp,sisfeat(:,:),pscore(mpisize,2),feat(:,:),dd1,dd2
real*8,allocatable:: score(:,:),score_sisfeat(:,:)
integer*8,allocatable:: complexity_sisfeat(:)
integer*8 i,j,k,l,ll,mpii,mpij,nfpcore_this(:),loc(2),simpler_result,complexity(:),pcomplexity(mpisize)
character(len=lname) fname(:),name_sisfeat(:),pfname(mpisize)
logical available(:),pavailable(mpisize),correlated

allocate(score(nfpcore_this(mpirank+1),2))
if(mpirank==0) then
  allocate(score_sisfeat(subs_sis(iFCDI),2))
  allocate(complexity_sisfeat(subs_sis(iFCDI)))
end if

usedup_warning=0
if(mpirank==0) write(*,'(/a)') 'Sure Independence Screening ...'

pavailable=.false.
if(nfpcore_this(mpirank+1)>0) call update_availability(available,pavailable)
do l=1,mpisize
 call mpi_bcast(pavailable(l),1,mpi_logical,l-1,mpi_comm_world,mpierr)
end do
if(nfpcore_this(mpirank+1)>0) score=-1   ! initial score

! get the scores
if(ptype==1) then
   do k=1,nfpcore_this(mpirank+1)
     if(available(k)) score(k,:)=sis_score(feat(:,k),trainy_c)
   end do
else if(ptype==2) then
   do k=1,nfpcore_this(mpirank+1)
     if(available(k)) score(k,:)=sis_score(feat(:,k),trainy) ! for classification
   end do
end if

! selection starts ...
i=0   ! count of selected features
do while( any(pavailable) .and. i<subs_sis(iFCDI) )

      i=i+1

      ! find the max score on each core
      if(nfpcore_this(mpirank+1)>0) then
        loc(2:2)=maxloc(score(:,1))
        tmp=score(loc(2),1)
      end if
      do l=1,nfpcore_this(mpirank+1)  ! equal scores
        if(  abs(tmp-score(l,1))<=1d-8 ) then
           if( score(l,2)-score(loc(2),2)>1d-8 .or. (abs(score(l,2)-score(loc(2),2))<=1d-8 .and. &
               simpler(complexity(l),complexity(loc(2)),fname(l),fname(loc(2)))==1 ) )  loc(2)=l
        end if
      end do
      if(nfpcore_this(mpirank+1)>0) then
        pscore(mpirank+1,:)=score(loc(2),:)
        pcomplexity(mpirank+1)=complexity(loc(2))
        pfname(mpirank+1)=fname(loc(2))
      else
        pscore(mpirank+1,:)=-1   ! simply a negative value to indicate empty
      end if

      !broadcast to send the information to all other cores
      do l=1,mpisize
        call mpi_bcast(pscore(l,:2),2,mpi_double_precision,l-1,mpi_comm_world,mpierr)
        call mpi_bcast(pcomplexity(l),1,mpi_integer8,l-1,mpi_comm_world,mpierr)
        call mpi_bcast(pfname(l),lname,mpi_character,l-1,mpi_comm_world,mpierr)
      end do

      !find the max score between cores
      loc(1:1)=maxloc(pscore(:,1))
      tmp=maxval(pscore(:,1))
      do l=1,mpisize
        if(  abs(tmp-pscore(l,1))<=1d-8 ) then
           if( pscore(l,2)-pscore(loc(1),2)>1d-8 .or.  ( abs(pscore(l,2)-pscore(loc(1),2))<=1d-8 .and. &
               simpler(pcomplexity(l),pcomplexity(loc(1)),pfname(l),pfname(loc(1)))==1 ) ) loc(1)=l
        end if
      end do

      ! save the highest-scored feature
      if((loc(1)-1)==mpirank) then
         if(mpirank==0) then
             sisfeat(:,i)=feat(:,loc(2))
             name_sisfeat(i)=fname(loc(2))
             score_sisfeat(i,:)=score(loc(2),:)
             complexity_sisfeat(i)=complexity(loc(2))
         else
            call mpi_send(feat(:,loc(2)),npoints,mpi_double_precision,0,100,mpi_comm_world,status,mpierr)
            call mpi_send(fname(loc(2)),lname,mpi_character,0,101,mpi_comm_world,status,mpierr)
            call mpi_send(complexity(loc(2)),1,mpi_integer8,0,102,mpi_comm_world,status,mpierr)
            call mpi_send(score(loc(2),:),2,mpi_double_precision,0,103,mpi_comm_world,status,mpierr)
         end if
         available(loc(2))=.false.   ! avoid to be selected again
         score(loc(2),:)=-1  
      end if

      if(mpirank==0 .and. mpirank/=(loc(1)-1) ) then
           call mpi_recv(sisfeat(:,i),npoints,mpi_double_precision,loc(1)-1,100,mpi_comm_world,status,mpierr)
           call mpi_recv(name_sisfeat(i),lname,mpi_character,loc(1)-1,101,mpi_comm_world,status,mpierr)
           call mpi_recv(complexity_sisfeat(i),1,mpi_integer8,loc(1)-1,102,mpi_comm_world,status,mpierr)
           call mpi_recv(score_sisfeat(i,:),2,mpi_double_precision,loc(1)-1,103,mpi_comm_world,status,mpierr)
      end if
      !---

      if(nfpcore_this(mpirank+1)>0) call update_availability(available,pavailable)
      do l=1,mpisize
       call mpi_bcast(pavailable(l),1,mpi_logical,l-1,mpi_comm_world,mpierr)
      end do

      if(ffdecorr) then
         correlated=.false.
         dd1=max(1d-8,decorr_delta)
         dd2=min(decorr_theta,1.d0-1.d-8)
         if(mpirank==0 .and. i>1) then
            do j=i-1,1,-1
            if(abs(score_sisfeat(j,1)-score_sisfeat(i,1))>dd1 .or. &
              (abs(score_sisfeat(j,1)-score_sisfeat(i,1))<=dd1 .and. &
               abs(score_sisfeat(j,2)-score_sisfeat(i,2))>dd1) )  exit ! if not similar scores

            if(abs(ffcorr(sisfeat(:,j),sisfeat(:,i)))>= dd2  ) then  ! decorrelation
                correlated=.true.
               if(score_sisfeat(j,1)<score_sisfeat(i,1) .or. &
                  (score_sisfeat(j,1)==score_sisfeat(i,1) .and. score_sisfeat(j,2)<score_sisfeat(i,2)) .or. &
                  (score_sisfeat(j,1)==score_sisfeat(i,1) .and. score_sisfeat(j,2)==score_sisfeat(i,2) .and. &
                  simpler(complexity_sisfeat(i),complexity_sisfeat(j),name_sisfeat(i),name_sisfeat(j))==1)) then
                  sisfeat(:,j)=sisfeat(:,i)
                  name_sisfeat(j)=name_sisfeat(i)
                  complexity_sisfeat(j)=complexity_sisfeat(i)
                  score_sisfeat(j,:)=score_sisfeat(i,:)
               end if
               exit
            end if
            end do
         end if
         call mpi_barrier(mpi_comm_world,mpierr)
         call mpi_bcast(correlated,1,mpi_logical,0,mpi_comm_world,mpierr)
         if(correlated) i=i-1
      end if

end do
!--------------------

do l=1,mpisize
  call mpi_bcast(usedup_warning(l),1,mpi_integer8,l-1,mpi_comm_world,mpierr)
end do

nsis(iFCDI)=i    ! the actual number of selected features
deallocate(score)
if(mpirank==0) then
  deallocate(score_sisfeat)
  deallocate(complexity_sisfeat)
end if

end subroutine


subroutine dup_scheck(nf,fidentity,fname,complexity,order,available,ftype)
! duplication check
! output order and available
integer*8 i,j,l,ll,order(:),n,nf,complexity(:),simpler_result
real*8 fidentity(:,:)
character(len=*) fname(:),ftype
logical available(:)

IF(nf==0) THEN
  n=0
ELSE

  ! ordering features fidentity from large to small
  order(1)=1
  n=1
  do i=2,nf

    l=0; ll=n;
    j=l+ceiling(float(ll-l)/2.0)

    123 continue
    simpler_result=simpler(complexity(i),complexity(order(j)),fname(i),fname(order(j)))

    ! the check for 'selected' is more strict than for 'all'.
    if(index(ftype,'selected')/=0) then
       ! get rid of duplication and correlation for selected features
       if(equivalent(fidentity(i,1:2),fidentity(order(j),1:2),f_select(:npoints,i),&
                     f_select(:npoints,order(j))) ) then
            if( simpler_result==1 ) then
               available(order(j))=.false.
               order(j)=i
            else
               available(i)=.false.
            end if
            cycle
       end if
      ! if not equal
      if(fidentity(i,1)>fidentity(order(j),1) .or. (fidentity(i,1)==fidentity(order(j),1) .and. &
         fidentity(i,2)>fidentity(order(j),2)) .or. (fidentity(i,1)==fidentity(order(j),1) .and. & 
         fidentity(i,2)==fidentity(order(j),2) .and. simpler_result==1 ) ) then
         ll=j
         if(j==l+ceiling(float(ll-l)/2.0)) then
           order(j+1:n+1)=order(j:n)
           order(j)=i
           n=n+1
           cycle
         end if
       else
          l=j
          if(j==l+ceiling(float(ll-l)/2.0)) then
            if(n>j) order(j+2:n+1)=order(j+1:n)
            order(j+1)=i
            n=n+1
            cycle
          end if
       end if
       j=l+ceiling(float(ll-l)/2.0)
       goto 123
    else if(index(ftype,'all')/=0) then
       if(abs(fidentity(i,1)-fidentity(order(j),1))<=1d-8 ) then  ! if equal
            if( simpler_result==1 ) then
               available(order(j))=.false.
               order(j)=i
            else
               available(i)=.false.
            end if
            cycle
        else  ! not equal  
            if(fidentity(i,1)>fidentity(order(j),1)) then
               ll=j
               if(j==l+ceiling(float(ll-l)/2.0)) then
                 order(j+1:n+1)=order(j:n)
                 order(j)=i
                 n=n+1
                 cycle
               end if
             else
                l=j
                if(j==l+ceiling(float(ll-l)/2.0)) then
                  if(n>j) order(j+2:n+1)=order(j+1:n)
                  order(j+1)=i
                  n=n+1
                  cycle
                end if
             end if

             j=l+ceiling(float(ll-l)/2.0)
             goto 123
        end if
    end if
 end do

END IF

order(nf+1)=n  ! store the number of features after check
if(mpirank/=0)  then
  call mpi_send(n,1,mpi_integer8,0,1,mpi_comm_world,status,mpierr)
else
  do l=1,mpisize-1
    call mpi_recv(i,1,mpi_integer8,l,1,mpi_comm_world,status,mpierr)
    n=n+i
 end do
end if
end subroutine


function sis_score(feat,yyy)   
! correlation between a feature feat and the property yyy
! sis_score returns a vector with 2 elements
integer i,j,mm1,mm2,mm3,mm4,k,kk,l,overlap_n,nf1,nf2,itask,nconvexpair
real*8 feat(:),sdfeat(ubound(feat,1)),tmp(ntask),sis_score(2),yyy(:),xnorm(ntask),xmean(ntask),&
       overlap_length,length_tmp,feat_tmp1(ubound(feat,1)),feat_tmp2(ubound(feat,1)),mindist,minlen
logical isoverlap


if(ptype==1) then  
   if(nreaction==0) then
     do j=1,ntask
         mm1=sum(nsample(:j-1))+1
         mm2=sum(nsample(:j))
         xmean(j)=sum(feat(mm1:mm2))/nsample(j)
         sdfeat(mm1:mm2)=feat(mm1:mm2)-xmean(j)  
         xnorm(j)=sqrt(sum((sdfeat(mm1:mm2))**2))
         if(xnorm(j)>1d-50) sdfeat(mm1:mm2)=sdfeat(mm1:mm2)/xnorm(j)  ! standardization to the features
         tmp(j)=abs(sum(sdfeat(mm1:mm2)*yyy(mm1:mm2)))   ! |xy|
     end do
   elseif(nreaction>0) then
     xmean(1)=sum(feat(:nreaction))/nreaction
     sdfeat(:nreaction)=feat(:nreaction)-xmean(1)
     xnorm(1)=sqrt(sum((sdfeat(:nreaction))**2))
     if(xnorm(1)>1d-50) sdfeat(:nreaction)=sdfeat(:nreaction)/xnorm(1)
     tmp(1)=abs(sum(sdfeat(:nreaction)*yyy))   ! |xy|
   end if
   ! score ranges from 0 to 1
   sis_score(1)=sqrt(sum(tmp**2)/ntask)  ! quadratic mean of the scores of different tasks
   sis_score(1)=sis_score(1)/sqrt(sum(yyy**2)/ntask)  ! normalization
   sis_score(2)=1.d0   ! not used

else if(ptype==2) then  
   mindist=-1d10
   overlap_n=0
   overlap_length=0.d0
   isoverlap=.false.
   nconvexpair=0
   do itask=1,ntask

   ! calculate overlap between domains of property i
        do i=1,ngroup(itask,1000)-1   ! ngroup(itask,1000) record the number of groups in this task
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
                  call convex1d_overlap(feat_tmp1(:nf1),feat_tmp2(:nf2),width,k,length_tmp)
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
                  overlap_n=overlap_n+convex1d_in(feat(mm3:mm4),feat_tmp1(:nf1),width)
               else if (isconvex(itask,i)==1 .and. isconvex(itask,j)==0) then !count the number of j-data in i-domain
                  overlap_n=overlap_n+convex1d_in(feat(mm1:mm2),feat_tmp2(:nf2),width)
               end if

          end do  ! j
        end do  ! i
   end do  ! itask

   sis_score(1)=float(overlap_n)  ! >=0,larger sis_score,worse feature
   sis_score(1)=1.d0/(sis_score(1)+1.d0)  ! transform to <=1, larger sis_score,better feature
   if(isoverlap) then  ! there are domains overlapped
     sis_score(2)=overlap_length/float(nconvexpair)     ! second metric, <=1, larger, worse
     sis_score(2)=1.d0/(sis_score(2)+1.0)  ! transform to <=1, larger, better
   else ! separated length
     sis_score(2)=-mindist   ! separated length between domains,positive
   end if

end if
end function


subroutine isgoodf(feat,name_feat,lastop,compl,dimens,nf)
real*8 feat(:),dimens(:)
character(len=*) name_feat,lastop
integer*8 nf,compl

if(goodf(feat,name_feat,dimens,compl)) then
  nf=nf+1
  if(foutsave) then
  if(nf>ubound(fout,2)) call addm_out(int8(ceiling(100000.0/mpisize)))
  fout(:,nf)=feat
  name_out(nf)=name_feat
  lastop_out(nf)=lastop
  complexity_out(nf)=compl
  dim_out(:,nf)=dimens
  end if
end if
end subroutine

subroutine isgoodvf(feat,name_feat,lastop,compl,dimens,nf)
real*8 feat(:,:),dimens(:),feat2(ubound(feat,1))
character(len=*) name_feat,lastop
integer*8 nf,d1,i,compl

d1=ubound(feat,1)

! transforming vector to scalar
do i=1,d1
  if(trim(vf2sf)=='sum') then
     feat2(i)=sum(feat(i,:))  ! sum of elements
  else if (trim(vf2sf)=='norm') then
     feat2(i)=sqrt(sum((feat(i,:))**2))
  else if (trim(vf2sf)=='min') then
     feat2(i)=minval(feat(i,:))
  else if(trim(vf2sf)=='max') then
     feat2(i)=maxval(feat(i,:))
  end if
end do

if(goodf(feat2,name_feat,dimens,compl)) then
  nf=nf+1
  nvf_new=nvf_new+1
  if(foutsave) then
  if(nf>ubound(fout,2)) call addm_out(int8(ceiling(100000.0/mpisize)))
  if(nvf_new>ubound(vfeat,3)) call addm_vf(int8(ceiling(100000.0/mpisize)))
  fout(1,nf)=nvf_new
  fout(2,nf)=sqrt(sum((tag+feat2)**2))
  vfeat(:,:,nvf_new)=feat
  name_out(nf)=name_feat
  lastop_out(nf)=lastop
  complexity_out(nf)=compl
  dim_out(:,nf)=dimens
  end if
end if
end subroutine


function isscalar(fname)
logical isscalar
character(len=*) fname
isscalar=.false.  
if(nvf==0) then
isscalar=.true.
return
else if(index(fname,'v_')==0) then
isscalar=.true.
end if
end function

subroutine update_select
! update the selected space
! bisection method for descending order
integer*8 i,j,k,l,ll,order(nselect),n,tmp,tmpcomplexity(nbasic_select),simpler_result
real*8 tmpf(ubound(f_select,1),nbasic_select),tmpscore(nbasic_select,2)
character(len=lname) tmpname(nbasic_select)

! ordering from large to small
order(1)=1   ! assuming the first feature being the best (highest score)
n=1  ! count of features saved in 'order'

do i=2,nselect   ! compare feaure i with the j located at middle of order(1:n)

  l=0; ll=n;
  j=l+ceiling(float(ll-l)/2.0)   ! j is at middle between l and ll.
 
  125 continue
  simpler_result=simpler(complexity_select(i),complexity_select(order(j)),name_select(i),name_select(order(j)))

  ! get rid of duplicate features
  if(equivalent(score_select(i,1:2),score_select(order(j),1:2),f_select(:npoints,i),&
                f_select(:npoints,order(j))) ) then
        if( simpler_result==1) order(j)=i  ! retain the simpler feature, and remove the other
        cycle
  end if

  ! bisection method for descending order
  if( (score_select(i,1)>score_select(order(j),1)) .or.  ((score_select(i,1)==score_select(order(j),1)) &
       .and. score_select(i,2)>score_select(order(j),2) ) .or.  ((score_select(i,1)==score_select(order(j),1)) &
       .and. score_select(i,2)==score_select(order(j),2) .and. simpler_result==1) )then
      ll=j   ! i is better, let the middle at j to be the new right end
      if(j==l+ceiling(float(ll-l)/2.0)) then  ! if the middle equal the last, insert i to position j
        order(j+1:n+1)=order(j:n)   ! shift the original j:n features backward by 1 step
        order(j)=i ! ordering over for this cycle
        n=n+1  ! size of 'order' increased by 1
        cycle
      end if
   else   ! j is better, let the middle at j to be the new left end
       l=j
       if(j==l+ceiling(float(ll-l)/2.0)) then
         if(n>j) order(j+2:n+1)=order(j+1:n)
         order(j+1)=i
         n=n+1
         cycle
       end if
   end if

   j=l+ceiling(float(ll-l)/2.0)  ! the new middle
   goto 125
end do

nfinal_select=min(n,nbasic_select)

! reordering
do i=1,nfinal_select
tmpf(:,i)=f_select(:,order(i))
tmpcomplexity(i)=complexity_select(order(i))
tmpname(i)=name_select(order(i))
tmpscore(i,:)=score_select(order(i),:)
end do

nselect=nfinal_select
f_select(:,:nselect)=tmpf(:,:nselect)
complexity_select(:nselect)=tmpcomplexity(:nselect)
name_select(:nselect)=tmpname(:nselect)
score_select(:nselect,:)=tmpscore(:nselect,:)
end subroutine


function ffcorr(feat1,feat2)
! calculate the correlation coefficient between two features
real*8 ffcorr,feat1(:),feat2(:),mean1,mean2,sd1,sd2
integer ns

ns=sum(nsample)
mean1=sum(feat1)/float(ns)
mean2=sum(feat2)/float(ns)
sd1=sqrt((sum((feat1-mean1)**2))/float(ns))
sd2=sqrt((sum((feat2-mean2)**2))/float(ns))
ffcorr=(sum((feat1-mean1)*(feat2-mean2)))/float(ns)/sd1/sd2
if(ffcorr>1.d0) ffcorr=1.d0
if(ffcorr<-1.d0) ffcorr=-1.d0
end function


function equivalent(score1,score2,feat1,feat2)
! calculate if two features are the same or highly correlated.
real*8 score1(:),score2(:),diff1,diff2,diff3,feat1(:),feat2(:)
logical equivalent

equivalent=.false.

diff1=score1(1)-score2(1)
diff2=score1(2)-score2(2)
diff3=1.d0-1d-8
if( abs(diff1)<=1d-8 .and. abs(diff2)<=1d-8 .and. abs(ffcorr(feat1,feat2))>=diff3 ) equivalent=.true.

end function


end module

