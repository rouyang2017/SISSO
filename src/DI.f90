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


module DI
! descriptor identification
!-------------------------------------------------------------

use var_global
use libsisso
! variables for this module
integer dim_DI

contains




subroutine descriptor_identification
real*8,allocatable:: xinput(:,:,:),yinput(:,:),beta(:,:),beta_init(:,:),coeff(:,:),xprime(:,:,:),&
        xdprime(:,:,:),xmean(:,:),norm_xprime(:,:),yprime(:,:),prod_xty(:,:),prod_xtx(:,:,:),lassormse(:),&
        ymean(:),intercept(:),weight(:,:)
real*8  lambda,lambda_max,alpha,L1_minrmse,totalrmse
integer run_iter,i,j,k,l,ll,nf,maxns,ntry,nactive
integer,allocatable:: idf(:),activeset(:),ncol(:)
character line*500
character,allocatable:: fname(:)*200
logical isnew,dat_readerr
real*8 stime_DI,etime_DI


if(mpirank==0) stime_DI=mpi_wtime()
if(nreaction==0) then
   maxns=maxval(nsample)
elseif(nreaction>0) then
   maxns=nreaction
endif
dim_DI=iFCDI

!--------------------------
allocate(xinput(maxns,fs_size_DI,ntask))
allocate(yinput(maxns,ntask))
allocate(fname(fs_size_DI))
allocate(weight(maxns,ntask))
allocate(activeset(fs_size_L0))

!------------------------
! job allocation
!------------------------
allocate(ncol(mpisize))
mpii=fs_size_DI/mpisize
mpij=mod(fs_size_DI,mpisize)
ncol(1)=mpii
do mpin=1,mpisize-1
if(mpin<=mpij) then
ncol(1+mpin)=mpii+1
else
ncol(1+mpin)=mpii
end if
end do

!-----------------------
! data read in
!-----------------------

dat_readerr=.true.

yinput=0
xinput=0
do i=1,ntask
  write(line,'(a,i3.3,a)') 'Uspace_p',i,'.dat'
  if(mpirank==0) then
  open(fileunit,file='SIS_subspaces/'//trim(adjustl(line)),status='old')
  if(nreaction==0) then
     if(ptype==1) then
        do j=1,nsample(i)
           read(fileunit,*,err=334) yinput(j,i),xinput(j,:,i)
        end do
     else
        do j=1,nsample(i)
           read(fileunit,*,err=334) xinput(j,:,i)
           yinput(j,i)=0.d0  
        end do
     end if
   elseif(nreaction>0) then ! ptype=1 only for reactionML
        do j=1,nreaction
           read(fileunit,*,err=334) yinput(j,i),xinput(j,:,i)
        end do
   end if    
   close(fileunit)
   end if
   call mpi_bcast(yinput,maxns*ntask,mpi_double_precision,0,mpi_comm_world,mpierr)
   call mpi_bcast(xinput,maxns*fs_size_DI*ntask,mpi_double_precision,0,mpi_comm_world,mpierr)
end do

dat_readerr=.false.
334 if(dat_readerr) stop 'Error when reading the space data files !!!'


if(mpirank==0) then
  open(fileunit,file='SIS_subspaces/Uspace.expressions',status='old')
  
  do i=1,fs_size_DI
  read(fileunit,'(a)') line
  line=trim(adjustl(line))
  fname(i)=line(1:index(line,' '))
  end do
  close(fileunit)
end if

call mpi_bcast(fname,fs_size_DI*200,mpi_character,0,mpi_comm_world,mpierr)

weight=1.0
if(task_weighting==1) then ! tasks treated equally
  weight=1.0/float(ntask)
else if(task_weighting==2) then ! tasks weighted by #sample_task_i/total_sample
  do i=1,ntask
     if(nreaction==0) then
        weight(:,i)=float(nsample(i))/float(sum(nsample))
     elseif(nreaction>0) then
        weight(:,i)=1.d0
     endif
  end do
end if

if(L1_weighted) then
open(fileunit,file='lasso.weight',status='old')
read(fileunit,*) ! title line
do i=1,ntask
  if(nreaction==0) then
    do j=1,nsample(i)
      read(fileunit,*) weight(j,i)
    end do
  elseif(nreaction>0) then
    do j=1,nreaction
      read(fileunit,*) weight(j,i)
    end do
  endif   
end do
close(fileunit)
end if


if(fs_size_L0>fs_size_DI) stop "Error: fs_size_L0 must not larger than fs_size_DI ! "

!-----------------------------------------------------
! MODEL SELECTION (L0 or L1L0)
!-----------------------------------------------------

if(trim(adjustl(method_so))=='L1L0' .and. fs_size_DI/=fs_size_L0 .and. ptype==1) then
    !----------------
    ! L1 of the L1L0
    ! see J. Friedman, T. Hastie, and R. Tibshirani 2010,"Regularization Paths for Generalized Linear Models via CD"
    ! "J. Friedman, T. Hastie, and R. Tibshirani, 2010 "A note on the group lasso and a sparse group lasso"
    !----------------
    allocate(xprime(maxns,fs_size_DI,ntask))
    allocate(xdprime(maxns,fs_size_DI,ntask))
    allocate(xmean(fs_size_DI,ntask))
    allocate(norm_xprime(fs_size_DI,ntask))
    allocate(prod_xtx(fs_size_DI,ncol(mpirank+1),ntask))
    allocate(intercept(ntask))
    allocate(beta(fs_size_DI,ntask)) 
    allocate(beta_init(fs_size_DI,ntask))
    allocate(yprime(maxns,ntask))
    allocate(prod_xty(fs_size_DI,ntask))
    allocate(lassormse(ntask))
    allocate(ymean(ntask))
    
    ! initialization for L1
    xprime=0
    xdprime=0
    xmean=0
    norm_xprime=0
    prod_xtx=0
    beta=0
    beta_init=0.0
    yprime=0
    prod_xty=0
    lassormse=0
    ymean=0
    
    !-------------------------------------------------------------------------------------------------------
    ! standardize and precompute the data for lasso
    ! transform y=a+xb into y=xb; y'=y-ymean; x'(:,j)=x(:,j)-xmean(j); 
    ! x(:,j)''=x(:,j)'/nomr(x'(j)); a=ymean-sum(xmean*b); b'=b*norm(x')
    ! after solving b', do b'->b and calculate a to get y=a+bx
    ! (observation) weighting is applied to square error, not the training property.
    ! see "J. Friedman, T. Hastie, and R. Tibshirani, 2010 "A note on the group lasso and a sparse group lasso"
    !--------------------------------------------------------------------------------------------------------
    if(mpirank==0) write(9,'(a)') 'data standardizing ... '
    !yprime
    do i=1,ntask
      if(nreaction==0) then
        ymean(i)=sum(yinput(:nsample(i),i))/nsample(i)
        yprime(:nsample(i),i)=yinput(:nsample(i),i)-ymean(i)
      elseif(nreaction>0) then
        ymean(i)=sum(yinput(:nreaction,i))/nreaction
        yprime(:nreaction,i)=yinput(:nreaction,i)-ymean(i)
      endif
    end do
    
    do i=1,ntask
       !xprime and xdprime
      if(nreaction==0) then
         do j=1,fs_size_DI
           xmean(j,i)=sum(xinput(:nsample(i),j,i))/nsample(i)
           xprime(:nsample(i),j,i)=xinput(:nsample(i),j,i)-xmean(j,i)
           norm_xprime(j,i)=sqrt(sum(xprime(:nsample(i),j,i)**2))
           if(abs(norm_xprime(j,i))<1d-10) then
           if(mpirank==0) print *, 'Error: norm of feature ',j,' of task ',i,' is zero'
           stop
           end if
           xdprime(:nsample(i),j,i)=xprime(:nsample(i),j,i)/norm_xprime(j,i)
         end do
         !xty
         prod_xty(:,i)=matmul(transpose(xdprime(:nsample(i),:,i)),weight(:nsample(i),i)*yprime(:nsample(i),i))
         !xtx=(xtx)t !prod_xtx(mpii,:)=prod_xtx(:,mpii)
         do mpii=1,ncol(mpirank+1)
         prod_xtx(:,mpii,i)=matmul(xdprime(:nsample(i),mpii+sum(ncol(1:mpirank)),i)*weight(:nsample(i),i),&
                                            xdprime(:nsample(i),:,i))
         end do
      elseif(nreaction>0) then
         do j=1,fs_size_DI
           xmean(j,i)=sum(xinput(:nreaction,j,i))/nreaction
           xprime(:nreaction,j,i)=xinput(:nreaction,j,i)-xmean(j,i)
           norm_xprime(j,i)=sqrt(sum(xprime(:nreaction,j,i)**2))
           if(abs(norm_xprime(j,i))<1d-10) then
           if(mpirank==0) print *, 'Error: norm of feature ',j,' of task ',i,' is zero'
           stop
           end if
           xdprime(:nreaction,j,i)=xprime(:nreaction,j,i)/norm_xprime(j,i)
         end do
         !xty
         prod_xty(:,i)=matmul(transpose(xdprime(:nreaction,:,i)),weight(:nreaction,i)*yprime(:nreaction,i))
         !xtx=(xtx)t !prod_xtx(mpii,:)=prod_xtx(:,mpii)
         do mpii=1,ncol(mpirank+1)
         prod_xtx(:,mpii,i)=matmul(xdprime(:nreaction,mpii+sum(ncol(1:mpirank)),i)*weight(:nreaction,i),&
                                            xdprime(:nreaction,:,i))
         end do
      endif
    end do
    
    !--------------------------------------------------------------------------------------------------------------
    ! get lambda_max
    !---------------------------------------------------------------------------------------------------------------
    lambda_max=0
    do j=1,fs_size_DI
    lambda=sqrt(sum(prod_xty(j,:)**2))
    if(lambda>lambda_max) lambda_max=lambda
    end do
    
    if(mpirank==0) then
    write(9,'(a)') 'feature selection by multi-task lasso starts ... '
    end if
    
    !----------------------------------------
    ! (Elastic net) H. Zou, and T. Hastie, J. R. Statist. Soc. B 67, 301 (2005).
    ! L1_elastic=1 is back to lasso (Eq. 14)
    !----------------------------------------
    !if(L1_elastic<1.0) then
    !prod_xtx=L1_elastic*prod_xtx
    !do j=1,fs_size_DI
    !prod_xtx(j,j,:)=prod_xtx(j,j,:)+(1-L1_elastic)
    !end do
    !end if
    
    !---------------------------------------------------------------------------
    ! Calculating lasso path with decreasing lambdas
    !---------------------------------------------------------------------------
    
    !counting the total size of the active set
    nactive=0
    ! iteration starts here
    do ntry=1,L1_nlambda
        ! create lambda sequence
        ! max and min in log space: max=log10(lambda_max),min=log10(0.001*lambda_max)
        ! density, i.e.: 100 points, interval=(max-min)/100=0.03
        ! lambda_i=10^(max-0.03*(i-1))
        lambda=10**(log10(lambda_max)-3.0/L1_dens*(ntry-1))
    
        ! call mtlasso
        call mtlasso_mpi(prod_xty,prod_xtx,lambda,L1_max_iter,L1_tole,beta_init,run_iter,beta,nf,ncol)
        
        ! lasso rmse
        do i=1,ntask
         if(nreaction==0) then
           lassormse(i)=sqrt(sum(weight(:nsample(i),i)*(yprime(:nsample(i),i)-matmul(xdprime(:nsample(i),:,i),&
                       beta(:,i)))**2)/nsample(i))  ! RMSE
         elseif(nreaction>0) then
           lassormse(i)=sqrt(sum(weight(:nreaction,i)*(yprime(:nreaction,i)-matmul(xdprime(:nreaction,:,i),&
                       beta(:,i)))**2)/nreaction)  ! RMSE
         endif
        end do
        totalrmse=sqrt(sum(lassormse**2))  ! the 'weight' already considered the averaging.
    
        ! L1_warm_start
        if(L1_warm_start) beta_init=beta
        
        ! intercept and beta
        do i=1,ntask
          intercept(i)=ymean(i)-sum(xmean(:,i)*beta(:,i)/norm_xprime(:,i))
        end do   
                       
        ! store the coeff. and id of features of this active set
        allocate(coeff(nf,ntask))
        allocate(idf(nf))
        coeff=0
        idf=0
        k=0
        do j=1,fs_size_DI
          if(maxval(abs(beta(j,:)))>1d-10) then
           k=k+1
           idf(k)=j
           coeff(k,:)=beta(j,:)
          end if
        end do
        if(k/=nf .and. mpirank==0 ) stop 'ID of selected features not correctly identified !'
        
        ! add to the total active set
        isnew=.true.  
        do k=1,nf
          if(any(idf(k)==activeset(:nactive))) isnew=.false.
          if(isnew) then
             nactive=nactive+1
             activeset(nactive)=idf(k)
          end if
          isnew=.true.
          if(nactive==fs_size_L0) exit
        end do
        
        ! output information of this lambda
        if(mpirank==0) then
        write(9,'(/)') 
        write(9,'(a,i5)') 'Try: ',ntry
        write(9,'(a)') '----------------'
        write(9,'(a30,e20.10)')  'Lambda: ', lambda
        write(9,'(a30,i10)')   'Number of iteration: ',run_iter
2004    format(a30,*(f10.6))
        write(9,2004)  ' Total LASSO RMSE: ', totalrmse
        write(9,2004)  ' LASSO RMSE for each property: ', lassormse
        write(9,'(a30,i6)')   'Size of this active set: ',nf
2005    format(a30,*(i8))
        if(nf/=0) write(9,2005)   'This active set: ',idf
        if(nf/=0) then
2006    format(a25,i3.3,a2,*(e20.10))
        do i=1,ntask
        write(9,2006) 'LASSO Coeff._',i,': ', coeff(:,i)
        end do
        end if
        write(9,'(a30,i6)') 'Size of the tot act. set: ',nactive
2007    format(a30,*(i8))
        if(nactive/=0) write(9,2007) 'Total active set: ',(activeset(j),j=1,nactive)
        end if
        deallocate(coeff)
        deallocate(idf)
        
        ! exit
        if(nactive>=fs_size_L0) then
        if(mpirank==0) write(9,'(/a,i3)') 'Total number of selected features >= ',fs_size_L0
        exit
        else if (nactive==fs_size_DI) then
        if(mpirank==0) write(9,'(/a)') 'The whole feature space is already selected!'
        exit
        else if (ntry==L1_nlambda) then
        if(mpirank==0) write(9,'(/a)') 'Number of lambda trials hits the max !'
        exit
        else if (maxval(lassormse)<L1_minrmse) then
        if(mpirank==0) write(9,'(/a)') 'RMSE is already smaller than the criterion!'
        exit
        end if
    
    end do
    deallocate(prod_xty)
    deallocate(prod_xtx)
    deallocate(ncol)
    deallocate(beta)
    deallocate(beta_init)
    deallocate(xprime)
    deallocate(xdprime)
    deallocate(xmean)
    deallocate(yprime)
    deallocate(intercept)
    deallocate(norm_xprime)
    deallocate(ymean)
    deallocate(lassormse)
    if(mpirank==0) write(9,'(a)') ' L1 finished ! '
    !---------------
    ! L0 of the L1L0
    !---------------
     call model(xinput,yinput,fname,nactive,activeset)
 
else if (trim(adjustl(method_so))=='L0' .or. fs_size_DI==fs_size_L0) then
    !--------
    ! L0
    !--------
    nactive=fs_size_L0
    do i=1,fs_size_L0
       activeset(i)=i
    end do
    if(ptype==1) then
     call model(xinput,yinput,fname,nactive,activeset)
    else if (ptype==2) then
     call model2(xinput,fname,nactive,activeset)
    end if
end if


! release spaces 
!--------------------
deallocate(xinput)
deallocate(yinput)
deallocate(fname)
deallocate(activeset)
deallocate(weight)

call mpi_barrier(mpi_comm_world,mpierr)
if(mpirank==0) then
  etime_DI=mpi_wtime()
  write(9,'(a,f15.2)') 'Time (second) used for this DI: ',etime_DI-stime_DI
end if

end subroutine


subroutine model(x,y,fname,nactive,activeset)
integer nactive,activeset(:),i,j,k,l,loc(1),isc,ii(dim_DI),idimen,select_model(max(nmodels,1),dim_DI),&
        mID(max(nmodels,1),dim_DI),sc_loc(1)
integer*8 totalm,bigm,nall,nrecord,nupper,nlower,njob(mpisize)
real*8 x(:,:,:),y(:,:),tmp,tmp2,yfit,wrmse(ntask),rmse(ntask),sc_rmse(2**dim_DI,ntask),maxae(ntask),&
sc_maxae(2**dim_DI,ntask),beta(dim_DI,ntask),intercept(ntask),select_rmse(max(nmodels,1)),&
select_score(max(nmodels,1),2,1+ntask),select_maxae(max(nmodels,1)),&
select_coeff(max(nmodels,1),dim_DI+1,ntask),mscore(max(nmodels,1),2,1+ntask),&
select_metric(max(nmodels,1)),mcoeff(max(nmodels,1),dim_DI+1,ntask),sc_beta(dim_DI,2**dim_DI,ntask),&
sc_intercept(2**dim_DI,ntask),sc_tmp(2**dim_DI),sc_tmp2(2**dim_DI),mpicollect(mpisize)
character fname(:)*200,line_name*100
logical isgood
real   progress

if(nmodels<1) nmodels=1

do idimen=1,dim_DI

   progress=0.2

   if(idimen<dim_DI) cycle ! .and. trim(adjustl(calc))=='FCDI') cycle

   ! job assignment for each CPU
   nall=1
   do mpii=1,idimen
      nall=nall*(nactive-mpii+1)/mpii
   end do
   if(nactive==idimen) nall=1
   njob=nall/mpisize

   mpij=mod(nall,int8(mpisize))
   do mpii=1,mpisize-1
     if(mpii<=mpij) njob(mpii+1)=njob(mpii+1)+1   
   end do

   ! initialization   
   select_rmse=10*maxval(abs(y))
   select_maxae=10*maxval(abs(y))
   totalm=0

! get the best nmodels models on each CPU core
   nrecord=0
   do k=1,nactive   ! features
      ii(1)=k   ! d1

      do j=2,idimen   ! d2, ..., d_idimen
        ii(j)=ii(j-1)+1   ! d1=k, d2=k+1,d3=k+2,...
      end do      

      nlower=sum(njob(:mpirank))+1
      nupper=sum(njob(:mpirank+1))
      123 continue 
      nrecord=nrecord+1
      if( nrecord< nlower) goto 124
      if( nrecord> nupper) exit

      if(float(nrecord-nlower+1)/float(njob(mpirank+1))>=progress) then
       write(*,'(a,i4,a,i4,a,f6.1,a)') 'dimension =',idimen,'   mpirank = ',mpirank,'   progress =',progress*100,'%'
        progress=progress+0.2
      end if

      if ( trim(adjustl(metric))=='RMSE' .or. trim(adjustl(metric))=='MaxAE') then
         do i=1,ntask
         if(nreaction==0) then
            if(.not. scmt) then
                 if(fit_intercept) then
                    call orth_de(x(:nsample(i),[activeset(ii(:idimen))],i),y(:nsample(i),i),&
                                 intercept(i),beta(:idimen,i),rmse(i))
                 else
                    call orth_de_nointercept(x(:nsample(i),[activeset(ii(:idimen))],i),y(:nsample(i),i),&
                                             beta(:idimen,i),rmse(i))
                    intercept(i)=0.d0
                 end if
                 if(isnan(rmse(i))) then
                    rmse(i)=sqrt(sum((y(:nsample(i),i)-sum(y(:nsample(i),i))/nsample(i))**2)/nsample(i))
                    intercept(i)=sum(y(:nsample(i),i))/nsample(i)
                    beta(:idimen,i)=0.d0
                 end if
                 maxae(i)= maxval(abs(y(:nsample(i),i)-(intercept(i)+matmul(x(:nsample(i),&
                                  [activeset(ii(:idimen))],i),beta(:idimen,i)))))
            else  ! sign-constrained multi-task learning
                 if(fit_intercept) then
                    call sc_coord_descent(x(:nsample(i),[activeset(ii(:idimen))],i),y(:nsample(i),i),10000,1d-5,&
                         idimen,sc_intercept(:2**idimen,i),sc_beta(:idimen,:2**idimen,i),sc_rmse(:2**idimen,i))
                 else
                    call sc_coord_descent_nointercept(x(:nsample(i),[activeset(ii(:idimen))],i),y(:nsample(i),i),&
                         10000,1d-5,idimen,sc_beta(:idimen,:2**idimen,i),sc_rmse(:2**idimen,i))
                    sc_intercept(:2**idimen,i)=0.d0
                 end if
                 if(isnan(rmse(i))) then
                    sc_rmse(:,i)=sqrt(sum((y(:nsample(i),i)-sum(y(:nsample(i),i))/nsample(i))**2)/nsample(i))
                    sc_intercept(:,i)=sum(y(:nsample(i),i))/nsample(i)
                    sc_beta(:idimen,:,i)=0.d0
                 end if
                 do isc=1,2**idimen
                    sc_maxae(isc,i)=maxval(abs(y(:nsample(i),i)-(sc_intercept(isc,i)+matmul(x(:nsample(i),&
                                  [activeset(ii(:idimen))],i),sc_beta(:idimen,isc,i)))))
                 end do
            end if
         elseif(nreaction>0) then
            if(fit_intercept) then
                call orth_de(x(:nreaction,[activeset(ii(:idimen))],i),y(:nreaction,i),intercept(i),&
                             beta(:idimen,i),rmse(i))
            else
                call orth_de_nointercept(x(:nreaction,[activeset(ii(:idimen))],i),y(:nreaction,i),&
                                         beta(:idimen,i),rmse(i))
                intercept(i)=0.d0
            end if
            if(isnan(rmse(i))) then
               rmse(i)=sqrt(sum((y(:nreaction,i)-sum(y(:nreaction,i))/nreaction)**2)/nreaction)
               intercept(i)=sum(y(:nreaction,i))/nreaction
               beta(:idimen,i)=0.d0
            end if
            maxae(i)=maxval(abs(y(:nreaction,i)-(intercept(i)+matmul(x(:nreaction,[activeset(ii(:idimen))],i),&
                        beta(:idimen,i)))))
         endif
         end do

         if( scmt ) then  ! best choice of coeff. signs
             do isc=1,2**idimen
                 if(task_weighting==1) then
                     sc_tmp(isc)=sqrt(sum(sc_rmse(isc,:)**2)/ntask)
                 else if(task_weighting==2) then
                     sc_tmp(isc)=sqrt(sum((sc_rmse(isc,:)**2)*nsample)/sum(nsample))
                 end if
                 sc_tmp2(isc)=maxval(sc_maxae(isc,:))
             end do
             if( trim(adjustl(metric))=='RMSE' ) then
                tmp=minval(sc_tmp)
                sc_loc=minloc(sc_tmp)
                tmp2=sc_tmp2(sc_loc(1))
             elseif( trim(adjustl(metric))=='MaxAE' ) then
                tmp2=minval(sc_tmp2)
                sc_loc=minloc(sc_tmp2)
                tmp=sc_tmp(sc_loc(1))
             endif
         else
              if(task_weighting==1) then
                 tmp=sqrt(sum(rmse**2)/ntask)  !overall error
              else if(task_weighting==2) then
                 if(nreaction==0) then
                    tmp=sqrt(sum((rmse**2)*nsample)/sum(nsample))
                 elseif(nreaction>0) then
                    tmp=rmse(1)
                 endif
              end if
              tmp2=maxval(maxae)        ! max(maxAE)
         end if
      end if
    
      isgood=.false.
      if ( trim(adjustl(metric))=='RMSE' .and. any(tmp <select_rmse)) then
         isgood = .true.
         loc=maxloc(select_rmse)
      else if ( trim(adjustl(metric))=='MaxAE' .and. any(tmp2 <select_maxae)) then
         isgood = .true.
         loc=maxloc(select_maxae)
      end if

      ! accept this model
      if (isgood) then
         totalm=totalm+1
         select_rmse(loc(1))=tmp
         select_maxae(loc(1))=tmp2
         select_model(loc(1),:idimen)=activeset(ii(:idimen))
         select_score(loc(1),1,1)=tmp
         select_score(loc(1),2,1)=tmp2
         do i=1,ntask
             if( scmt ) then
               select_coeff(loc(1),1,i)=sc_intercept(sc_loc(1),i)
               select_coeff(loc(1),2:idimen+1,i)=sc_beta(:idimen,sc_loc(1),i)
               select_score(loc(1),1,1+i)=sc_rmse(sc_loc(1),i)
               select_score(loc(1),2,1+i)=sc_maxae(sc_loc(1),i)
             else
               select_coeff(loc(1),1,i)=intercept(i)
               select_coeff(loc(1),2:idimen+1,i)=beta(:idimen,i)
               select_score(loc(1),1,1+i)=rmse(i)
               select_score(loc(1),2,1+i)=maxae(i)
             end if
         end do
      end if

      124 continue
      ! update the models(123,124,125,134,135,145,234,...)
      if(idimen==1) cycle
      ii(idimen)=ii(idimen)+1  ! add 1 for the highest dimension

      do j=idimen,3,-1
        if(ii(j)> (nactive-(idimen-j)) ) ii(j-1)=ii(j-1)+1 
      end do
      do j=3,idimen
        if(ii(j)> (nactive-idimen+j) ) ii(j)=ii(j-1)+1
      end do

      if(ii(2)>(nactive-(idimen-2)))  cycle

      goto 123
   end do

! collecting the best models
   if(mpirank>0) then
     call mpi_send(totalm,1,mpi_integer8,0,1,mpi_comm_world,status,mpierr)
   else
     do i=1,mpisize-1
       call mpi_recv(bigm,1,mpi_integer8,i,1,mpi_comm_world,status,mpierr)
       totalm=totalm+bigm
     end do
   end if
   call mpi_bcast(totalm,1,mpi_integer8,0,mpi_comm_world,mpierr)

   if(trim(adjustl(metric))=='RMSE') then
     select_metric=select_rmse
   else if (trim(adjustl(metric))=='MaxAE') then
     select_metric=select_maxae     
   end if

   do j=1,min(int8(nmodels),totalm)
      loc=minloc(select_metric)
      k=loc(1)
      if(mpirank>0) then
        call mpi_send(select_metric(loc(1)),1,mpi_double_precision,0,1,mpi_comm_world,status,mpierr)
      else
        mpicollect(1)=select_metric(loc(1))
        do i=1,mpisize-1
          call mpi_recv(mpicollect(i+1),1,mpi_double_precision,i,1,mpi_comm_world,status,mpierr)
        end do
        loc=minloc(mpicollect)
      end if
      call mpi_bcast(loc(1),1,mpi_integer,0,mpi_comm_world,mpierr)

      if(mpirank==loc(1)-1) then
        mID(j,:idimen)=select_model(k,:idimen)
        mcoeff(j,:idimen+1,:)=select_coeff(k,:idimen+1,:)
        mscore(j,:2,:)=select_score(k,:2,:)
        select_metric(k)=10*maxval(abs(y))
      end if
      call mpi_bcast(mID(j,:idimen),idimen,mpi_integer,loc(1)-1,mpi_comm_world,mpierr)
      call mpi_bcast(mscore(j,:2,:),2*(1+ntask),mpi_double_precision,loc(1)-1,mpi_comm_world,mpierr)
      call mpi_bcast(mcoeff(j,:idimen+1,:),(idimen+1)*ntask,mpi_double_precision,loc(1)-1,mpi_comm_world,mpierr)
   end do
 

   if(mpirank==0) then
      call writeout(idimen,mID(1,:idimen),mscore(1,1,2:1+ntask),mscore(1,2,2:1+ntask),&
                    mcoeff(1,2:,:),mcoeff(1,1,:),x,y,fname)
      write(line_name,'(a,i4.4,a,i3.3,a)') 'top',min(int8(nmodels),totalm),'_',idimen,'d'
      open(fileunit,file='models/'//trim(adjustl(line_name)),status='replace')
      write(fileunit,'(4a12)') 'Rank','RMSE','MaxAE','Feature_ID'
2008  format(i12,2f12.6,a,*(i8))
      do i=1,min(int8(nmodels),totalm)
        write(fileunit,2008,advance='no') i,mscore(i,:,1),'  (',mID(i,:idimen)
        write(fileunit,'(a)') ')'
      end do
      close(fileunit)

      write(line_name,'(a,i4.4,a,i3.3,a)') 'top',min(int8(nmodels),totalm),'_',idimen,'d_coeff'
      open(fileunit,file='models/'//trim(adjustl(line_name)),status='replace')
      write(fileunit,'(a)') 'Model_ID, [(c_i,i=0,n)_j,j=1,ntask]'
20081 format(i12,*(e20.10))
      do i=1,min(int8(nmodels),totalm)
        write(fileunit,20081) i,(mcoeff(i,:idimen+1,j),j=1,ntask)
      end do
      close(fileunit)

   end if

end do

call mpi_barrier(mpi_comm_world,mpierr)
end subroutine


subroutine writeout(idimen,id,rmse,maxae,betamin,interceptmin,x,y,fname)
integer idimen,i,j,id(:),k
real*8 rmse(:),maxae(:),betamin(:,:),interceptmin(:),x(:,:,:),y(:,:),yfit
character fname(:)*200,line_name*100

if(idimen==desc_dim) then ! .and. trim(adjustl(calc))=='FCDI') then
  write(9,'(/a)')  'Final model/descriptor !'
  write(9,'(a)')'================================================================================'
!else if(idimen<desc_dim) then ! .and. trim(adjustl(calc))=='FCDI') then
!  write(9,'(/a)')  'Model/descriptor for generating residual:'
end if
write(9,'(i3,a)') idimen,'D descriptor (model): '

if(task_weighting==1) then  ! tasks are treated equally
   write(9,'(a,2f10.6)')  'RMSE and MaxAE: ',sqrt(sum(rmse**2)/ntask),maxval(maxae)
else if(task_weighting==2) then  ! weighted by sample size
   if(nreaction==0) then
     write(9,'(a,2f10.6)')  'RMSE and MaxAE: ',sqrt(sum((rmse**2)*nsample)/sum(nsample)),maxval(maxae)
   elseif(nreaction>0) then
     write(9,'(a,2f10.6)')  'RMSE and MaxAE: ',rmse,maxval(maxae)
   endif
end if

write(9,'(a)')  '@@@descriptor: '
do i=1,idimen
write(9,'(i23,3a)') id(i),':[',trim(adjustl(fname(id(i)))),']'
end do

2009 format(a20,i3.3,a2,*(e20.10))
2010 format(a20,a3,a2,*(e20.10))
2011 format(a10,2a20,*(a19,i1.1))
2012 format(i10,*(e20.10))

do i=1,ntask
  write(9,2009) ' coefficients_',i,': ', betamin(:idimen,i)
  write(9,'(a20,i3.3,a2,e20.10)') '  Intercept_',i,': ', interceptmin(i)
  write(9,'(a20,i3.3,a2,2e20.10)') ' RMSE,MaxAE_',i,': ', rmse(i),maxae(i)
  write(line_name,'(a,i3.3,a,i3.3,a)') 'desc_',idimen,'d_p',i,'.dat'
  open(fileunit,file='desc_dat/'//trim(adjustl(line_name)),status='replace')
  write(fileunit,2011) 'Index','y_true','y_pred',('descriptor_',j,j=1,idimen)
  if(nreaction==0) then
    do j=1,nsample(i)
       yfit=interceptmin(i)+sum(x(j,[id(:idimen)],i)*betamin(:idimen,i))
       write(fileunit,2012) j,y(j,i),yfit,x(j,[id(:idimen)],i)
    end do
  elseif(nreaction>0) then
    do j=1,nreaction
       yfit=interceptmin(i)+sum(x(j,[id(:idimen)],i)*betamin(:idimen,i))
       write(fileunit,2012) j,y(j,i),yfit,x(j,[id(:idimen)],i)
    end do
  endif
  close(fileunit)
end do

if(idimen==desc_dim) then 
  if(ntask>1) then
     if(task_weighting==1) then
     write(9,'(a,f10.6)') ' RMSE(task_weighting=1) = sqrt(sum(RMSE_i^2)/ntask) = ',sqrt(sum(rmse**2)/ntask)
     else if(task_weighting==2) then
     write(9,'(a,f10.6)') ' RMSE(task_weighting=2) = sqrt(sum(RMSE_i^2*n_i)/sum(n_i)) = ', &
                           sqrt(sum((rmse**2)*nsample)/sum(nsample))
     end if
     write(9,'(a,f10.6)') ' MaxAE = max(MaxAE_i) = ',maxval(maxae)
  end if
  write(9,'(a)')'================================================================================'
else
  write(9,'(a)')'--------------------------------------------------------------------------------'
end if

end subroutine


! descriptor for classification
subroutine model2(x,fname,nactive,activeset)
integer nactive,activeset(:),i,j,k,l,loc(1),ii(dim_DI),itask,mm1,mm2,mm3,mm4,ns,&
idimen,select_model(max(nmodels,1),dim_DI),mID(max(nmodels,1),dim_DI),overlap_n,overlap_n_tmp,&
select_overlap_n(max(nmodels,1)),nh,ntri,nconvexpair
integer*8 totalm,bigm,nall,nrecord,nlower,nupper
real*8  x(:,:,:),mscore(max(nmodels,1),2),overlap_area,overlap_area_tmp,select_overlap_area(max(nmodels,1)),&
hull(ubound(x,1),2),area(maxval(ngroup(:,1000))),mindist,xtmp1(ubound(x,1),3),xtmp2(ubound(x,1),3)
character fname(:)*200,line_name*100
integer*8 njob(mpisize)
real*8 mpicollect(mpisize,2)
logical isoverlap
real progress
integer,allocatable:: triangles(:,:)

if(nmodels<1) nmodels=1

do idimen=1,dim_DI

  progress=0.2

  if(idimen<dim_DI) cycle ! .and. trim(adjustl(calc))=='FCDI') cycle
  if(idimen>3) exit  ! up to 3D

   ! job assignment for each CPU core
   nall=1
   do mpii=1,idimen
      nall=nall*(nactive-mpii+1)/mpii
   end do
   if(nactive==idimen) nall=1
   njob=nall/mpisize

   mpij=mod(nall,int8(mpisize))
   do mpii=1,mpisize-1
     if(mpii<=mpij) njob(mpii+1)=njob(mpii+1)+1   
   end do

   ! initialization   
   select_overlap_n=sum(nsample)*sum(ngroup(:,1000))  ! set to a large value
   totalm=0

! get the best nmodels models on each CPU core
   nrecord=0
   do k=1,nactive
      ii(1)=k

      do j=2,idimen
        ii(j)=ii(j-1)+1   ! starting number
      end do      

      nlower=sum(njob(:mpirank))+1
      nupper=sum(njob(:mpirank+1))
      123 continue 
      nrecord=nrecord+1

      if( nrecord< nlower) goto 124
      if( nrecord> nupper) exit

      if(float(nrecord-nlower+1)/float(njob(mpirank+1))>=progress) then
        write(*,'(a,i4,a,i4,a,f6.1,a)') 'dimension =',idimen,'   mpirank =',mpirank,'   progress =',progress*100,'%'
        progress=progress+0.2
      end if

      ! initialization
      overlap_n=0
      overlap_area=0.d0
      mindist=-1d10
      isoverlap=.false.
      nconvexpair=0

      !------ get the score for this model ----------
      do itask=1,ntask
        ! area
           if(idimen==1) then
             do i=1,ngroup(itask,1000) ! number of groups in the task $itask
              if(isconvex(itask,i)==0) cycle
              mm1=sum(ngroup(itask,:i-1))+1 
              mm2=sum(ngroup(itask,:i))
              area(i)=maxval(x(mm1:mm2,[activeset(ii(:idimen))],itask))-minval(x(mm1:mm2,[activeset(ii(:idimen))],itask))
             end do
           else if(idimen==2) then
             do i=1,ngroup(itask,1000)
              if(isconvex(itask,i)==0) cycle
              mm1=sum(ngroup(itask,:i-1))+1 
              mm2=sum(ngroup(itask,:i))
              area(i)=convex2d_area(x(mm1:mm2,[activeset(ii(:idimen))],itask))
             end do
           else if(idimen==3) then
             area(:ngroup(itask,1000))=1.d0  ! not yet implemented for 3D volumes
           end if

        !  model scores
           do i=1,ngroup(itask,1000)-1
             mm1=sum(ngroup(itask,:i-1))+1
             mm2=sum(ngroup(itask,:i))
             do j=i+1,ngroup(itask,1000)
               mm3=sum(ngroup(itask,:j-1))+1
               mm4=sum(ngroup(itask,:j))

               IF(isconvex(itask,i)==1 .and. isconvex(itask,j)==1) then
                 nconvexpair=nconvexpair+1
                 if(idimen==1) then
                     xtmp1(mm1:mm2,1)=x(mm1:mm2,activeset(ii(idimen)),itask)
                     xtmp2(mm3:mm4,1)=x(mm3:mm4,activeset(ii(idimen)),itask)
                     call convex1d_overlap(xtmp1(mm1:mm2,1),xtmp2(mm3:mm4,1),bwidth,overlap_n_tmp,overlap_area_tmp)
                 else if(idimen==2) then
                     xtmp1(mm1:mm2,:2)=x(mm1:mm2,activeset(ii(:2)),itask)
                     xtmp2(mm3:mm4,:2)=x(mm3:mm4,activeset(ii(:2)),itask)
                     call convex2d_overlap(xtmp1(mm1:mm2,:2),xtmp2(mm3:mm4,:2),bwidth,overlap_n_tmp,overlap_area_tmp)
                 else if(idimen==3) then
                     xtmp1(mm1:mm2,:3)=x(mm1:mm2,activeset(ii(:3)),itask)
                     xtmp2(mm3:mm4,:3)=x(mm3:mm4,activeset(ii(:3)),itask)
                     call convex3d_overlap(xtmp1(mm1:mm2,:3),xtmp2(mm3:mm4,:3),bwidth,overlap_n_tmp)
                     overlap_area_tmp=0.d0  ! Calculating overlap_area is not implemented
                 end if
                 overlap_n=overlap_n+overlap_n_tmp
                 if(overlap_area_tmp>=0) isoverlap=.true.
 
                 if(.not. isoverlap) then  ! if separated
                     if(mindist<overlap_area_tmp) then
                         mindist=overlap_area_tmp ! renew the worst separation
                         overlap_area=mindist  ! smaller, better
                     end if
                 else 
                     if(overlap_area_tmp>=0.d0 .and. min(area(i),area(j))==0.d0) then !  0D feature
                         overlap_area=max(0.d0,overlap_area)+1.d0   ! totally overlapped
                    else if (overlap_area_tmp>=0.d0 .and. min(area(i),area(j))>0.d0) then ! not 0D feature
                         overlap_area=max(0.d0,overlap_area)+overlap_area_tmp/(min(area(i),area(j)))  ! total overlap
                    end if
                 end if

              ELSE IF(isconvex(itask,i)==1 .and. isconvex(itask,j)==0) then
                 if(idimen==1) then
                     xtmp1(mm1:mm2,1)=x(mm1:mm2,activeset(ii(idimen)),itask)
                     xtmp2(mm3:mm4,1)=x(mm3:mm4,activeset(ii(idimen)),itask)
                     overlap_n=overlap_n+convex1d_in(xtmp1(mm1:mm2,1),xtmp2(mm3:mm4,1),bwidth)
                 else if(idimen==2) then
                     xtmp1(mm1:mm2,:2)=x(mm1:mm2,activeset(ii(:2)),itask)
                     xtmp2(mm3:mm4,:2)=x(mm3:mm4,activeset(ii(:2)),itask)
                     overlap_n=overlap_n+convex2d_in(xtmp1(mm1:mm2,:2),xtmp2(mm3:mm4,:2),bwidth)
                 else if(idimen==3) then
                     xtmp1(mm1:mm2,:3)=x(mm1:mm2,activeset(ii(:3)),itask)
                     xtmp2(mm3:mm4,:3)=x(mm3:mm4,activeset(ii(:3)),itask)
                     call convex3d_hull(xtmp1(mm1:mm2,:3),ntri,triangles)
                     overlap_n=overlap_n+convex3d_in(xtmp1(mm1:mm2,:3),xtmp2(mm3:mm4,:3),bwidth,ntri,triangles)
                     deallocate(triangles)
                 end if
              ELSE IF(isconvex(itask,i)==0 .and. isconvex(itask,j)==1) then
                 if(idimen==1) then
                     xtmp1(mm1:mm2,1)=x(mm1:mm2,activeset(ii(idimen)),itask)
                     xtmp2(mm3:mm4,1)=x(mm3:mm4,activeset(ii(idimen)),itask)
                     overlap_n=overlap_n+convex1d_in(xtmp2(mm3:mm4,1),xtmp1(mm1:mm2,1),bwidth)
                 else if(idimen==2) then
                     xtmp1(mm1:mm2,:2)=x(mm1:mm2,activeset(ii(:2)),itask)
                     xtmp2(mm3:mm4,:2)=x(mm3:mm4,activeset(ii(:2)),itask)
                     overlap_n=overlap_n+convex2d_in(xtmp2(mm3:mm4,:2),xtmp1(mm1:mm2,:2),bwidth)
                 else if(idimen==3) then
                     xtmp1(mm1:mm2,:3)=x(mm1:mm2,activeset(ii(:3)),itask)
                     xtmp2(mm3:mm4,:3)=x(mm3:mm4,activeset(ii(:3)),itask)
                     call convex3d_hull(xtmp2(mm3:mm4,:3),ntri,triangles)
                     overlap_n=overlap_n+convex3d_in(xtmp2(mm3:mm4,:3),xtmp1(mm1:mm2,:3),bwidth,ntri,triangles)
                     deallocate(triangles)
                 end if
              END IF

             end do ! j
           end do ! i
         end do ! itask
         
         if(isoverlap)  overlap_area=overlap_area/float(nconvexpair) ! smaller, better
         !----------------------------------

      ! store good models 
      if (any(overlap_n<select_overlap_n)) then
         totalm=totalm+1
         loc=maxloc(select_overlap_n)
         select_overlap_n(loc(1))=overlap_n
         select_overlap_area(loc(1))=overlap_area
         select_model(loc(1),:idimen)=activeset(ii(:idimen))
      else if (overlap_n==maxval(select_overlap_n) .and. any(overlap_area<select_overlap_area) )  then
         totalm=totalm+1
         loc=maxloc(select_overlap_area)
         select_overlap_n(loc(1))=overlap_n
         select_overlap_area(loc(1))=overlap_area
         select_model(loc(1),:idimen)=activeset(ii(:idimen))         
      end if

      ! update models(123,124,125,134,135,145,234,...)
      124 continue
      if(idimen==1) cycle
      ii(idimen)=ii(idimen)+1  ! add 1 for the highest dimension

      do j=idimen,3,-1
        if(ii(j)> (nactive-(idimen-j)) ) ii(j-1)=ii(j-1)+1 
      end do
      do j=3,idimen
        if(ii(j)> (nactive-idimen+j) ) ii(j)=ii(j-1)+1
      end do
      if(ii(2)>(nactive-(idimen-2)))  cycle
      goto 123
   end do

! collecting the best models from all CPU cores
   if(mpirank>0) then
     call mpi_send(totalm,1,mpi_integer8,0,1,mpi_comm_world,status,mpierr)
   else
     do i=1,mpisize-1
       call mpi_recv(bigm,1,mpi_integer8,i,1,mpi_comm_world,status,mpierr)
       totalm=totalm+bigm
     end do
   end if
   call mpi_bcast(totalm,1,mpi_integer8,0,mpi_comm_world,mpierr)

   do j=1,min(int8(nmodels),totalm)
      loc=minloc(select_overlap_n)
      do k=1,min(int8(nmodels),totalm)
        if(select_overlap_n(k)==select_overlap_n(loc(1)) .and. &
           select_overlap_area(k)<select_overlap_area(loc(1))) loc(1)=k
      end do
      k=loc(1)
      mscore(j,:)=(/dble(select_overlap_n(k)),select_overlap_area(k)/)

      if(mpirank>0) then
        call mpi_send(mscore(j,:),2,mpi_double_precision,0,1,mpi_comm_world,status,mpierr)
      else
        mpicollect(1,:)=mscore(j,:)
        do i=1,mpisize-1
          call mpi_recv(mpicollect(i+1,:),2,mpi_double_precision,i,1,mpi_comm_world,status,mpierr)
        end do
        loc=minloc(mpicollect(:,1))
        do l=1,mpisize
          if(mpicollect(l,1)==mpicollect(loc(1),1) .and. &
             mpicollect(l,2)<mpicollect(loc(1),2) ) loc(1)=l
        end do
        l=loc(1)
      end if
      call mpi_bcast(l,1,mpi_integer,0,mpi_comm_world,mpierr)

      if(mpirank==l-1) then
        mID(j,:idimen)=select_model(k,:idimen)
        select_overlap_n(k)=sum(nsample)*sum(ngroup(:,1000))   ! set to large to avoid to be selected again
      end if
      call mpi_bcast(mID(j,:idimen),idimen,mpi_integer,l-1,mpi_comm_world,mpierr)
      call mpi_bcast(mscore(j,:),2,mpi_double_precision,l-1,mpi_comm_world,mpierr)
   end do

   if(mpirank==0) then
      call writeout2(x,idimen,mID(1,:idimen),fname,mscore)

2013  format(i10,i20,f20.5,a,*(i8))
2014  format(i10,i20,a,*(i8))
      if(nmodels>1) then
        write(line_name,'(a,i4.4,a,i3.3,a)') 'top',min(int8(nmodels),totalm),'_',idimen,'d'
        open(fileunit,file='models/'//trim(adjustl(line_name)),status='replace')
        if(idimen<=2) then
            write(fileunit,'(a10,3a20)') 'Rank','#data_overlap','size_overlap','Feature_ID'
        else 
            write(fileunit,'(a10,2a20)') 'Rank','#data_overlap','Feature_ID'
        end if
        do i=1,min(int8(nmodels),totalm)
           if(idimen<=2) then
               write(fileunit,2013,advance='no') i,nint(mscore(i,1)),mscore(i,2),'   (',mID(i,:idimen)
           else
               write(fileunit,2014,advance='no') i,nint(mscore(i,1)),'   (',mID(i,:idimen)
           end if
           write(fileunit,'(a)') ')'
        end do
        close(fileunit)
      end if

      if(idimen>=2) then
        if(idimen==2) then
           open(fileunit,file='convex2d_hull',status='replace')
           write(fileunit,'(a)') 'The coordinates (x,y) of 2D-vertices arranged in clockwise direction'
        end if
        if(idimen==3) then
           open(fileunit,file='convex3d_hull',status='replace')
           write(fileunit,'(a)') 'Facet-triangles with each indicated by three indexes of the data of this class'
           write(fileunit,'(a)') 'E.g.: ploting a 3D convex hull by using the MATLAB function: trisurf(Tri,X,Y,Z)'
        end if

         do itask=1,ntask
         write(fileunit,'(a,i5)') 'Task :',itask
         do i=1,ngroup(itask,1000)
           mm1=sum(ngroup(itask,:i-1))+1
           mm2=sum(ngroup(itask,:i))
           if(isconvex(itask,i)==0) cycle
           write(fileunit,'(a,i5,g20.10)') 'Hull',i
           if(idimen==2) then
              call convex2d_hull(x(mm1:mm2,[mID(1,:idimen)],itask),j,hull)
              do k=1,j
                write(fileunit,'(2e20.10)') hull(k,:)
              end do
           elseif(idimen==3) then
              call convex3d_hull(x(mm1:mm2,[mID(1,:idimen)],itask),j,triangles)
              do k=1,j
                write(fileunit,'(3i5)') triangles(k,:)
              end do
              deallocate(triangles) ! allocated in the convex3d_hull subroutine
           end if
         end do
         end do
         close(fileunit)
      end if
   end if
end do

call mpi_barrier(mpi_comm_world,mpierr)
end subroutine


subroutine writeout2(x,idimen,id,fname,mscore)
integer idimen,i,j,id(:),k,l,itask,mm1,mm2,mm3,mm4,ndata_ol,ntri
integer,allocatable:: triangles(:,:)
real*8 mscore(:,:),tmp,tmp2,x(:,:,:)
character fname(:)*200,line_name*100
logical inside,convexpair

2015 format(2a20,*(a19,i1.1))
2016 format(i20,a20,*(e20.10))

ndata_ol=0
do itask=1,ntask
  write(line_name,'(a,i3.3,a,i3.3,a)') 'desc_',idimen,'d_p',itask,'.dat'
  open(fileunit,file='desc_dat/'//trim(adjustl(line_name)),status='replace')
  write(fileunit,2015) 'index','classified?',('descriptor_',j,j=1,idimen)
  do i=1,ngroup(itask,1000)  ! the groups in this task $itask
  do j=sum(ngroup(itask,:i-1))+1,sum(ngroup(itask,:i))  ! samples in this group
     inside=.false.
     do k=1,ngroup(itask,1000)  ! check if sample j of group i is in the domain of group k
       if(i==k) cycle
       if(isconvex(itask,k)==0) cycle ! group k is not convex domain
       mm1=sum(ngroup(itask,:k-1))+1
       mm2=sum(ngroup(itask,:k))
       if(idimen==1) then
            if(convex1d_in(x(mm1:mm2,id(idimen),itask),x(j:j,id(idimen),itask),bwidth)==1) then
              inside=.true.
              exit
            end if
       elseif(idimen==2) then
           if(convex2d_in(x(mm1:mm2,[id(:idimen)],itask),x(j:j,[id(:idimen)],itask),bwidth)==1) then
             inside=.true.
             exit
           end if
       elseif(idimen==3) then
           call convex3d_hull(x(mm1:mm2,[id(:idimen)],itask),ntri,triangles)
           if(convex3d_in(x(mm1:mm2,[id(:idimen)],itask),x(j:j,[id(:idimen)],itask),bwidth,ntri,triangles)==1) then
             inside=.true.
             deallocate(triangles)
             exit
           end if
           deallocate(triangles)
      end if
    end do

    if(inside) then
        write(fileunit,2016) j,'NO',x(j,[id(:idimen)],itask)  ! unclassified
        ndata_ol=ndata_ol+1
    else
        write(fileunit,2016) j,'YES',x(j,[id(:idimen)],itask) ! classified
    end if
  end do
  end do
  close(fileunit)
end do

if(idimen==desc_dim) then
  write(9,'(/a)')  'Final model/descriptor !'
  write(9,'(a)')'================================================================================'
!else if(idimen<desc_dim) then 
!  write(9,'(/a)')  'Model/descriptor for generating residual:'
end if

write(9,'(i3,a)') idimen,'D descriptor (model): '
write(9,'(a,i10)') 'Number of data in all overlap regions:',int(mscore(1,1))
convexpair=.false.
do itask=1,ntask
  do i=1,ngroup(itask,1000)
    do j=i+1,ngroup(itask,1000)
       if(isconvex(itask,i)==1 .and. isconvex(itask,j)==1) convexpair=.true.
    end do
  end do
end do
if(idimen<=2 .and. convexpair) then
  write(9,'(a,i10)') 'Actual number (without double counting) of data in all overlap regions: ', ndata_ol
  write(9,'(a,f15.5)') 'Size of the overlap:',mscore(1,2)
  write(9,'(a)') 'Note: positive value indicates the overlap length/area/volumn; negative means no overlap,'
  write(9,'(a)') '      with the magnitude being the smallest distance between all domains.'
end if
write(9,'(a)')  '@@@descriptor: '
do i=1,idimen
write(9,'(i23,3a)') id(i),':[',trim(adjustl(fname(id(i)))),']'
end do

if(idimen==desc_dim) then
  write(9,'(a)')'================================================================================'
else
  write(9,'(a)')'--------------------------------------------------------------------------------'
end if

end subroutine

end module

