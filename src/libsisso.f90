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


module libsisso
!***********************************************************************
! Subroutines:
! ----------
! det: calculates the determinant of a given matrix
! inverse: calculates the inverse of a given matrix
! orth_de: linear least square fitting by orthogonal decomposition
! worth_de: observation weighted orth_de
! lls: linear least fitting by conventional method (not robust)
! coord_descent: coordinate descent
! crosspro: crossproduct between two vectors
! normp: norm vector for a plane 
! LASSO: least absolute shrinkage and selection operator
! mtlasso_mpi Multi-task LASSO
! corr:  pearson's correlation
! string_split: break a string into sub-strings.
! kfoldCV: k-fold CV for linear model
! convex2d_xxx: construction of convex hull
!************************************************************************

use var_global

contains





subroutine string_split(instr,outstr,sp)
! break a string into sub-strings
! input: instr, input string; sp, separator
! output: outstr, output sub-strings
character(len=*) instr,outstr(:),sp
integer n
logical isend

isend=.false.
n=0
outstr=''

i=index(instr,sp)
if(i/=0) then
  if(i/=1) then
    n=n+1
    outstr(n)=instr(1:i-1)
  end if
  do while ((.not. isend) .and. n<ubound(outstr,1))
    if ((i+len(sp))<len(instr)) then
        j=index(instr(i+len(sp):),sp)
        if(j==0) then
          isend=.true.
          n=n+1
          outstr(n)=instr(i+len(sp):)
        else if(j>1) then
          n=n+1
          outstr(n)=instr(i+len(sp):i+len(sp)-1+j-1)
        end if
        i=i+len(sp)+j-1
     else
       isend=.true.
     end if
  end do
end if

end subroutine


function corr(x,y)
! pearson's correlation
! input: vector x and y
! output: the correlation coefficient
real*8 x(:),y(:),meanx,meany,sigmax,sigmay,m,corr
m=size(x)

meanx=sum(x)/m
meany=sum(y)/m
sigmax=sqrt(sum((x-meanx)**2)/m)
sigmay=sqrt(sum((y-meany)**2)/m)

corr=sum((x-meanx)*(y-meany))/(m*sigmax*sigmay)
end function

function det(mat)
! LU docomposition for sovling determinant of a matrix
! input: matrix mat
! output: the determinant.

integer i,j,k,n
real*8 mat(:,:),um(ubound(mat,1),ubound(mat,1)),s,det,temp(ubound(mat,1))
s=0
um=mat
n=ubound(mat,1)

do i=1,n
do j=i+1,n

if (um(i,i)==0.d0) then
  do k=i+1,n
    if(um(k,i) /=0.d0) then
      temp=um(i,:)
      um(i,:)=um(k,:)
      um(k,:)=-1*temp
      exit
    end if
  end do
end if
if (um(i,i)==0.d0) then
det=0;return
end if

um(j,:)=um(j,:)-um(j,i)/um(i,i)*um(i,:)
end do
end do

do i=1,n
if(i==1) then
s=um(i,i)
else
s=s*um(i,i)
endif
end do

det=s
end function

function inverse(mat)
! calculate the inverse of a given matrix
! input: the matrix
! output: the inversed matrix

real*8  mat(:,:),um(ubound(mat,1),ubound(mat,1)),lm(ubound(mat,1),ubound(mat,1))
real*8  x(ubound(mat,1),ubound(mat,1)),y(ubound(mat,1),ubound(mat,1)),inverse(ubound(mat,1),ubound(mat,1))
integer i,j,k,n
um=mat
n=ubound(mat,1)
lm=0;x=0;y=0
if (abs(det(mat))==0) stop 'Error: value of the determinant is 0, the matrix is singular'
if (abs(det(mat))<1d-30) write(*,'(a,E15.6E3,a)') 'Warning: value of the determinant is:',&
det(mat),' matrix might be singular'

!!!!!! construct up triangle matrix
do i=1,n
do j=i+1,n
um(j,:)=um(j,:)-um(j,i)/um(i,i)*um(i,:)
end do
end do

! construct low triangle matrix
lm(1,1)=mat(1,1)/um(1,1)
do i=2,n
do j=1,i
s=sum(lm(i,1:(i-1))*um(1:(i-1),j))
lm(i,j)=(mat(i,j)-s)/um(j,j)
end do
end do

! construct y matrix
y(1,1)=1/lm(1,1)
do i=2,n
do j=1,i
s=sum(lm(i,1:(i-1))*y(1:(i-1),j))
if(j/=i) y(i,j)=-s/lm(i,i)
if(j==i) y(i,j)=(1-s)/lm(i,i)
end do
end do

! construct x matrix, which is also the inverse matrix
do i=1,n
x(n,i)=y(n,i)/um(n,n)
end do
do i=n-1,1,-1
do j=1,n,1
s=sum(um(i,i+1:n)*x(i+1:n,j))
x(i,j)=(y(i,j)-s)/um(i,i)
end do
end do

inverse=x

end function

subroutine qr_de(a,q,r)
! QR decomposition of input matrix a
! https://en.wikipedia.org/wiki/QR_decomposition
integer i,j,k,m,n
real*8 a(:,:),q(ubound(a,1),ubound(a,2)),r(ubound(a,2),ubound(a,2))
real*8 u(ubound(a,1),ubound(a,2)),e(ubound(a,1),ubound(a,2))
m=ubound(a,1)
n=ubound(a,2)
u=0
e=0

do i=1,n
if(i==1) then
u(:,1)=a(:,1)
else
u(:,i)=a(:,i)
do j=1,i-1
u(:,i)=u(:,i)-dot_product(u(:,j),a(:,i))/dot_product(u(:,j),u(:,j))*u(:,j)
end do
end if
e(:,i)=u(:,i)/sqrt(dot_product(u(:,i),u(:,i)))
end do
q=e

do j=1,n
do i=1,j
r(i,j)=sum(e(:,i)*a(:,j))
end do

do k=j+1,n
r(k,j)=0.0
end do
end do
end subroutine


subroutine orth_de(x,y,intercept,beta,rmse)
! linear least square fit to y=a+x*b by Orthogonal decomposition methods method
! input: matrix x,vector y; 
! output intercept,beta,rmse
! https://en.wikipedia.org/wiki/Linear_least_squares_(mathematics)

real*8  x(:,:),y(:),beta(ubound(x,2)),xprime(ubound(x,1),ubound(x,2)),yprime(ubound(y,1)),intercept,rmse
real*8  q(ubound(x,1),ubound(x,2)),r(ubound(x,2),ubound(x,2)),qty(ubound(x,2)),xmean(ubound(x,2)),ymean
integer i,j,k,m,n

m=ubound(x,1)
n=ubound(x,2)
do i=1,n
xmean(i)=sum(x(:,i))/m
xprime(:,i)=x(:,i)-xmean(i)
end do
ymean=sum(y)/m
yprime=y-ymean

call qr_de(xprime,q,r) ! results stored in q and r
qty=matmul(transpose(q),yprime)
beta(n)=qty(n)/r(n,n)

do i=n,1,-1
beta(i)=qty(i)
do j=i+1,n
beta(i)=beta(i)-r(i,j)*beta(j)
end do
beta(i)=beta(i)/r(i,i)
end do

rmse=sqrt(sum((yprime-matmul(xprime,beta))**2)/m)
intercept=ymean-sum(xmean*beta)

end subroutine

subroutine worth_de(x,y,weight,intercept,beta,rmse,wrmse)
! weighted linear least sqaure
implicit none
real*8  x(:,:),y(:),beta(ubound(x,2)),xprime(ubound(x,1),ubound(x,2)),yprime(ubound(y,1)),intercept,wrmse,rmse,&
q(ubound(x,1),ubound(x,2)),r(ubound(x,2),ubound(x,2)),qty(ubound(x,2)),xmean(ubound(x,2)),ymean,weight(:),&
wx(ubound(x,1),ubound(x,2)),wy(ubound(y,1))
integer i,j,k,m,n

m=ubound(x,1)
n=ubound(x,2)

do i=1,n
xmean(i)=sum(x(:,i))/m
xprime(:,i)=x(:,i)-xmean(i)
wx(:,i)=xprime(:,i)*sqrt(weight)
end do
ymean=sum(y)/m
yprime=y-ymean
wy=yprime*sqrt(weight)


call qr_de(wx,q,r) ! results stored in q and r
qty=matmul(transpose(q),wy)
beta(n)=qty(n)/r(n,n)

do i=n,1,-1
beta(i)=qty(i)
do j=i+1,n
beta(i)=beta(i)-r(i,j)*beta(j)
end do
beta(i)=beta(i)/r(i,i)
end do

wrmse=sqrt(sum(weight*(yprime-matmul(xprime,beta))**2)/m)
rmse=sqrt(sum((yprime-matmul(xprime,beta))**2)/m)
intercept=ymean-sum(xmean*beta)

end subroutine


subroutine lls(x,y,intercept,beta,rmse)
! linear least square fit by the standard method
! input: matrix x, vector y
! output: beta,rmse, intercept
! https://en.wikipedia.org/wiki/Linear_least_squares_(mathematics)
implicit none
real*8  x(:,:),y(ubound(x,1)),beta0(ubound(x,2)+1),rmse,x0(ubound(x,1),ubound(x,2)+1),&
        beta(ubound(x,2)),intercept
integer m,n
m=ubound(x,1)
n=ubound(x,2)

x0(:,1)=1
x0(:,2:n+1)=x

beta0=matmul(inverse(matmul(transpose(x0),x0)),matmul(transpose(x0),y))
rmse=sqrt(sum((y-matmul(x0,beta0))**2)/m)
intercept=beta0(1)
beta=beta0(2:n)

end subroutine


function crosspro(a,b)
! calculate the cross product: a x b
! input: vector a and b
! output: a vector

real*8  a(3),b(3),crosspro(3)
crosspro(1)=a(2)*b(3)-a(3)*b(2);
crosspro(2)=a(3)*b(1)-a(1)*b(3);
crosspro(3)=a(1)*b(2)-a(2)*b(1);
end function


subroutine lasso(prod_xty,prod_xtx,lambda,max_iter,tole,beta_init,run_iter,beta,nf)
! min: f(beta)=1/2*||y-x*beta||**2 + lambda*||beta||_l1
! prod_xty=XtY; prod_xtx=XtX; max_iter:allowed max cycle; 
! tole: convergence criteria for cd to stop
! beta_init: initial coefficient
! run_iter, actual run cycles; nf: number of selected features; 
! input prod_xty,prod_xtx,lambda,max_iter,tole,beta_init; 
! output run_iter,beta, nf
! (LASSO) J. Friedman, T. Hastie, R. Tibshirani, J. Stat. Softw. 33, 1(2010).
! (LASSO) J. Friedman, T. Hastie, H. Hofling, R. Tibshirani, Ann. Appl. Stat. 1, 302 (2007).
implicit none
real*8 prod_xty(:),prod_xtx(:,:),lambda,beta(:),beta_init(:),&
beta_old(ubound(beta,1)),z,tole,rand(ubound(beta,1)),beta_tmp,dbeta,xxbeta(ubound(beta,1))
integer ntotf,i,i1,i2(1),j,k,max_iter,run_iter,nf,ac1(ubound(beta,1)),ac2(ubound(beta,1))
logical active_change

ntotf=ubound(prod_xtx,1)

beta_old=beta_init
beta=beta_init
xxbeta=matmul(prod_xtx,beta)

!-----------------------------------
! entering the main updation cycles
!-----------------------------------

do i=1,max_iter

   if(i==1) then
      active_change=.true.
      ac1=0; ac2=0
   end if

  do j=1,ntotf
   if((.not. active_change) .and. abs(beta(j))<1d-10) then
     cycle
   else
       beta_tmp=beta(j)
       z=prod_xty(j)-xxbeta(j)+beta(j)
      if(z>0 .and. lambda<abs(z)) then
       beta(j)=z-lambda
      else if(z<0 .and. lambda<abs(z)) then
       beta(j)=z+lambda
      else if(lambda>=abs(z)) then
       beta(j)=0
      end if
      xxbeta=xxbeta-prod_xtx(:,j)*beta_tmp
      xxbeta=xxbeta+prod_xtx(:,j)*beta(j)
   end if
  end do   
  active_change=.false. 

  if(i==max_iter) print*, 'Warning: lasso run hits the max_iter.'
  if(maxval(abs(beta-beta_old))<tole .or. i==max_iter) then
    ! check if the active set converge
     do k=1,ntotf
        if(abs(beta(k))>1d-10) then
          ac2(k)=1
        else
          ac2(k)=0
        end if
     end do

     do k=1,ntotf
       if(ac1(k)/=ac2(k)) then
         ac1=ac2
         active_change=.true.
         exit
       end if
     end do
     if((.not. active_change)) exit
  else
    beta_old=beta
  end if

end do

nf=sum(ac2)
run_iter=i

end subroutine

subroutine mtlasso_mpi(prod_xty,prod_xtx,lambda,max_iter,tole,beta_init,run_iter,beta,nf,ncol)
! prod_xty=XtY; prod_xtx=XtX; max_iter:allowed max cycle; 
! tole: convergence criteria for cd to stop
! beta_init: initial coefficient
! run_iter, actual run cycles; nf: number of selected features; 
! ncol: mpi jobs assignment
! output run_iter,beta, nf
! (LASSO) J. Friedman, T. Hastie, R. Tibshirani, J. Stat. Softw. 33, 1(2010).
! (LASSO) J. Friedman, T. Hastie, H. Hofling, R. Tibshirani, Ann. Appl. Stat. 1, 302 (2007).
! (MTLASSO) G. Obozinski, B. Taskar, M. Jordan, 2006 "Multi-task feature selection"
! multi-task lasso return back to lasso when the number of task is one
! algorithm for sovling multi-task lasso can be found from that of group lasso
! (GLASSO) M. Yuan and Y. Lin, J. R. Statist. Soc. B 68,49(2006)
! (GLASSO) J. Friedman, T. Hastie, and R. Tibshirani, 2010 "A note on the group lasso and a sparse group lasso"
! (Elastic net) H. Zou, and T. Hastie, J. R. Statist. Soc. B 67, 301 (2005).

real*8 prod_xty(:,:),prod_xtx(:,:,:),lambda,beta(:,:),beta_init(:,:),beta_old(ubound(beta,1),ubound(beta,2)),tole,&
beta_tmp(ubound(beta,2)),beta_tmp2(ubound(beta,2)),dbeta,xxbeta(ubound(beta,1),ubound(beta,2)),&
Sj(ubound(beta,2)),norm
integer ntotf,i,j,k,max_iter,run_iter,nf,ac1(ubound(beta,1)),ac2(ubound(beta,1)),ntask,mpii,mpij,mpik,mpin,ncol(:)
logical active_change
!----

ntask=ubound(prod_xty,2)
ntotf=ubound(prod_xty,1)
beta_old=beta_init
beta=beta_init

! calculate xtx*beta
do mpij=1,ncol(mpirank+1)
do mpik=1,ntask
  xxbeta(mpij,mpik)=sum(prod_xtx(:,mpij,mpik)*beta(:,mpik))
end do
end do

! rank0 collect data and broadcast
mpin=ncol(mpirank+1)
if (mpirank /=0) then
call mpi_send(xxbeta(1:mpin,:),mpin*ntask,mpi_double_precision,0,1,mpi_comm_world,mpierr)
else
do mpii=1,mpisize-1
mpij=sum(ncol(1:mpii))
mpik=ncol(mpii+1)
call mpi_recv(xxbeta(mpij+1:mpij+mpik,:),mpik*ntask,mpi_double_precision,mpii,1,mpi_comm_world,status,mpierr)
end do
end if

call mpi_barrier(mpi_comm_world,mpierr)
call mpi_bcast(xxbeta,ntotf*ntask,mpi_double_precision,0,mpi_comm_world,mpierr)

!-----------------------------------
! entering the main updation cycles
!-----------------------------------
do i=1,max_iter

   if(i==1) then
      active_change=.true.
      ac1=0; ac2=0
   end if

  if(mpirank/=0) then  ! do the job processor after processor,rank0 first
   call mpi_recv(xxbeta,ntotf*ntask,mpi_double_precision,mpirank-1,1,mpi_comm_world,status,mpierr)
  end if

  do j=1+sum(ncol(1:mpirank)),sum(ncol(1:mpirank+1))

    if((.not. active_change) .and. maxval(abs(beta(j,:)))<1d-10) then
     cycle
    else
       ! Multi-task lasso and group lasso share the similar solution
       ! see paper "A note on the group lasso and a sparse group lasso" by J. Friedman, et al., 2010
       ! and also note that features are the same for all responses
       ! if ntask==1, back to lasso automatically.
       Sj=prod_xty(j,:)-xxbeta(j,:)+prod_xtx(j,j-sum(ncol(1:mpirank)),:)*beta(j,:) 
       norm=sqrt(sum(Sj**2))
       beta_tmp=beta(j,:)
       if(norm<=lambda) then
         beta(j,:)=0.0
       else
         ! only when features are orthogonal
         !beta(j,:)=Sj*(1-lambda/norm)

         ! do it iteratively when features are not orthogonal
           beta_tmp2=beta(j,:)+1
           norm=sqrt(sum((beta(j,:))**2))
           if(norm<1d-10) norm=1.0
           do while(maxval(abs(beta(j,:)-beta_tmp2))>1d-10)
             beta_tmp2=beta(j,:)
             beta(j,:)=Sj/(prod_xtx(j,j-sum(ncol(1:mpirank)),:)+lambda/norm)
             norm=sqrt(sum((beta(j,:))**2))
           end do
         !--
       end if
       do k=1,ntask
       xxbeta(:,k)=xxbeta(:,k)-prod_xtx(:,j-sum(ncol(1:mpirank)),k)*beta_tmp(k)
       xxbeta(:,k)=xxbeta(:,k)+prod_xtx(:,j-sum(ncol(1:mpirank)),k)*beta(j,k)
       end do
   end if

 end do


 if(mpirank/=mpisize-1) then
 call mpi_send(xxbeta,ntotf*ntask,mpi_double_precision,mpirank+1,1,mpi_comm_world,mpierr)
 end if
 call mpi_barrier(mpi_comm_world,mpierr)
 call mpi_bcast(xxbeta,ntotf*ntask,mpi_double_precision,mpisize-1,mpi_comm_world,mpierr)

  mpin=ncol(mpirank+1)
  if (mpirank /=0) then
  call mpi_send(beta(1+sum(ncol(1:mpirank)):sum(ncol(1:mpirank+1)),:),ncol(mpirank+1)*ntask,&
                     mpi_double_precision,0,1,mpi_comm_world,mpierr)
  else
  do mpii=1,mpisize-1
  mpij=sum(ncol(1:mpii))
  mpik=ncol(mpii+1)
  call mpi_recv(beta(mpij+1:mpij+mpik,:),mpik*ntask,mpi_double_precision,mpii,1,mpi_comm_world,status,mpierr)
  end do
  end if
  call mpi_bcast(beta,ntotf*ntask,mpi_double_precision,0,mpi_comm_world,mpierr)
  call mpi_barrier(mpi_comm_world,mpierr)

  active_change=.false.

  if(i==max_iter) print*, 'Warning: lasso run hits the max_iter.'
  if(maxval(abs(beta-beta_old))<tole .or. i==max_iter) then
    ! check if the active set converge
     do k=1,ntotf
        if(maxval(abs(beta(k,:)))>1d-10) then
          ac2(k)=1
        else
          ac2(k)=0
        end if
     end do

     do k=1,ntotf
       if(ac1(k)/=ac2(k)) then
         ac1=ac2
         active_change=.true.
         exit
       end if
     end do
     if((.not. active_change)) exit
  else
    beta_old=beta
  end if

end do

nf=sum(ac2)
run_iter=i

end subroutine


subroutine coord_descent(x,y,max_iter,tole,beta,rmse)
! using coordinate descent method to solve the linear equation y=xb
! input: x(m,n),y(m),max_iter,tole; output beta(n),rmse
! x and y are already standardized, and y=xb
! output: rmse,beta
implicit none
real*8 x(:,:),y(ubound(x,1)),beta(ubound(x,2)),beta_old(ubound(x,2)),&
prod_xty(ubound(x,2)),prod_xtx(ubound(x,2),ubound(x,2)),rmse,tole,norm_beta
integer m,n,i,j,k,max_iter

m=ubound(x,1)
n=ubound(x,2)

beta_old=0.0
beta=0.0

! definition to make fast
prod_xty=matmul(transpose(x),y)
prod_xtx=matmul(transpose(x),x)

! entering the main updation cycles
do i=1,max_iter
 do j=1,n
       beta(j)=prod_xty(j)-sum(prod_xtx(j,:)*beta)+beta(j)
 end do
   norm_beta=sqrt(sum(beta**2))
  if(norm_beta==0 .or. (sqrt(sum((beta-beta_old)**2))/norm_beta)<tole .or. i==max_iter) then
    exit
  else
    beta_old=beta
  end if
end do

rmse=sqrt(sum((y-matmul(x,beta))**2)/m)

end subroutine

subroutine kfoldCV(x,y,random,fold,noise,CVrmse,CVmax)
! k-fold cross validation
! input: x, y, random, fold, noise
! output: CVrmse, CVmax,

integer ns,fold,random(:),mm1,mm2,mm3,mm4,i,j,k,kk,l
real*8 x(:,:),y(:),beta(ubound(x,2)),intercept,rmse,CVse(fold),CVrmse,CVmax,pred(ubound(y,1)),noise(:)

if(fold<2) stop 'Error: the fold of CV must be >=2 !'
ns=ubound(y,1)

CVmax=0
CVrmse=0
k=int(ns/fold)
kk=mod(ns,fold)
do l=1,fold
   mm1=1  ! sample start
   mm2=ns      ! sample end
   mm3=(l-1)*k+min((l-1),kk)+1 ! test start
   mm4=mm3+k-1+int(min(l,kk)/l)  ! test end
   if(mm1==mm3) then
      call orth_de(x([random(mm4+1:mm2)],:),y([random(mm4+1:mm2)])+noise([random(mm4+1:mm2)]),intercept,beta(:),rmse)
   else if(mm4==mm2) then
      call orth_de(x([random(mm1:mm3-1)],:),y([random(mm1:mm3-1)])+noise([random(mm1:mm3-1)]),intercept,beta(:),rmse)
   else
      call orth_de(x([random(mm1:mm3-1),random(mm4+1:mm2)],:),&
                   y([random(mm1:mm3-1),random(mm4+1:mm2)])+noise([random(mm1:mm3-1),random(mm4+1:mm2)]),&
                   intercept,beta(:),rmse)
   end if
   pred(mm3:mm4)=(intercept+matmul(x([random(mm3:mm4)],:),beta(:)))
end do
CVmax=maxval(abs(y([random])-pred))
CVrmse=sqrt(sum((y([random])-pred)**2)/ns)  ! quadratic mean over samples
end subroutine


subroutine convex2d_hull(set,numb,hull)
! calculate the convex hull for a given data set
! input: set, the data set {(x1,y1),(x2,y2),...}
! output: numb (number of vertices); hull, the vertices stored in clockwise direction
real*8 set(:,:),hull(:,:),tmp,vjj(2),vij(2),vkj(2),normij,normkj,set_line(ubound(set,1),2)
integer numb,i,j,k,ntot,loc(1),nrecord
logical used(ubound(set,1)),distinct

ntot=ubound(set,1)
used=.false.
numb=0


! find the the most left and lowest point
loc=minloc(set(:,1)) 
j=loc(1)
do i=1,ntot  
  if(i==j) cycle
  if(set(loc(1),1)==set(i,1) .and. set(loc(1),2)>set(i,2)) loc(1)=i  
end do

j=loc(1)
numb=1
hull(1,:)=set(j,:)

! at the same positions with j
do i=1,ntot
   if(i==j) cycle
   if(abs(set(j,1)-set(i,1))<=1d-10 .and. abs(set(j,2)-set(i,2))<=1d-10) then
   numb=numb+1
   hull(numb,:)=set(i,:)
   used(i)=.true.
  end if
end do


if(numb==ntot) return

! start searching vertices

vjj=(/0.d0,-1.d0/)  ! initial vector must point downward

nrecord=0
do while(any(used .eqv. .false.) )
  nrecord=nrecord+1
  distinct=.false.
  vkj=vjj
  normkj=-1.d0

  ! update j, vjj, numb, hull, used
  do i=1,ntot

    if(used(i) .or. i==j ) cycle

    if(abs(set(j,1)-set(i,1))<=1d-10 .and. abs(set(j,2)-set(i,2))<=1d-10  ) then   ! i and j at the same position
       numb=numb+1
       hull(numb,:)=set(i,:)
       used(i)=.true.
       cycle
    end if

    vij=set(i,:)-set(j,:)     
    normij=sqrt(sum(vij**2))
    vij=vij/normij

    tmp=vkj(1)*vij(2)-vij(1)*vkj(2)  ! cross product between vector vij and vkj

    if(i==loc(1)) normij=1d10   

    if(tmp>1d-10   .or. & 
       (normkj>0.d0 .and. tmp>0.d0 .and. tmp<=1d-10 .and. sum(vij*vkj)>0.d0 ) .or. &
       (normkj<0.d0 .and. tmp>0.d0 .and. tmp<=1d-10 .and. sum(vij*vjj)<0.d0 ) ) then  ! vij is on right size of vkj/vjj
       k=i
       vkj=vij
       normkj=normij
       distinct=.true.
    else if ( (normkj>0.d0 .and. tmp==0.d0  .and. sum(vij*vkj)>0.d0 ) .or. &
              (normkj<0.d0 .and. tmp==0.d0  .and. sum(vij*vjj)<0.d0 ) ) then 
       if(normkj<0.d0 .or. normij<normkj) then   ! i is closer than k to i or i=loc(1)
         k=i
         vkj=vij
         normkj=normij
         distinct=.true.
       end if
    end if
  end do

  if(nrecord>ntot) then
     k=loc(1)
     distinct=.true.
  end if

  if(distinct) then
    if(k==loc(1)) exit
     numb=numb+1
     hull(numb,:)=set(k,:)
     used(k)=.true.
     vjj=-vkj
     j=k
  end if


end do

end subroutine


subroutine convex2d_overlap(set1,set2,width,numb,area)
! counting the number of data and area in the overlapped region between two convex hulls
! input: data set 1, data set 2,width(boundary tolerance)
! output: number of data in the overlapped region, area of the overlapped region
! segment intersection: (http://dec3.jlu.edu.cn/webcourse/t000096/graphics/chapter5/01_1.html)
! if two convex domain has more than two intersection points, the two domains are considered totally overlap !

integer i,j,i2,j2,k,numb,nh1,nh2,ns1,ns2,n_intersect
real*8 set1(:,:),set2(:,:),area,hull1(ubound(set1,1),2),hull2(ubound(set2,1),2),width,norm1,norm2,vunit1(2),vunit2(2),&
set3(10*(ubound(set1,1)+ubound(set2,1)),2),tmp,segp(4,2),delta,lambda,mu,xa,xb,ya,yb,xc,xd,yc,yd
logical inside,polygon1,polygon2

call convex2d_hull(set1,nh1,hull1)  ! convex hull of data set 1 
call convex2d_hull(set2,nh2,hull2)  ! convex hull of data set 2
ns1=ubound(set1,1)
ns2=ubound(set2,1)
numb=0
area=0.d0
polygon1=.true.
polygon2=.true.

if(convex2d_area(hull1(:nh1,:))<1d-10) polygon1=.false.
if(convex2d_area(hull2(:nh2,:))<1d-10) polygon2=.false.

! both the sets are not polygons
if((ns1>2 .and. (.not. polygon1)) .or. (ns2>2 .and. (.not. polygon2)) ) then   ! two lines
  numb=ns1+ns2      ! set to max since this is not wanted.
  area=0.d0
else
  ! set1 is polygon, count the number of set2 data in hull1
  if(polygon1) then
  do i=1,ns2
    ! check if data from set2 is inside hull1
    inside=.true.
    do j=1,nh1
       k=j+1
       if(j==nh1) k=1

       norm1=sqrt(sum((hull1(k,:)-hull1(j,:))**2))
       if(norm1<1d-10) cycle
       vunit1=(hull1(k,:)-hull1(j,:))/norm1
       norm2=sqrt(sum((set2(i,:)-hull1(j,:))**2))
       if(norm2<1d-10) cycle
       vunit2=(set2(i,:)-hull1(j,:))/norm2
       tmp=vunit1(1)*vunit2(2)-vunit2(1)*vunit1(2)  ! cross product

       if(tmp>0.d0) then  
          if(width>0.d0) then
            if(convex2d_dist(set2(i:i,:2),hull1(:nh1,:2))>width) then  ! boundary tolerance
             inside=.false.
             exit  
            end if
          else
             inside=.false.
             exit
          end if
      end if
    end do
    if(inside) then
       numb=numb+1
       set3(numb,:)=set2(i,:)
    end if
  end do
  end if

  ! set1 data in hull2
  if(polygon2) then
  do i=1,ns1
    ! check if set1 data is inside hull2
    inside=.true.
    do j=1,nh2
       k=j+1
       if(j==nh2) k=1

       norm1=sqrt(sum((hull2(k,:)-hull2(j,:))**2))
       if(norm1<1d-10) cycle
       vunit1=(hull2(k,:)-hull2(j,:))/norm1
       norm2=sqrt(sum((set1(i,:)-hull2(j,:))**2))
       if(norm2<1d-10) cycle
       vunit2=(set1(i,:)-hull2(j,:))/norm2
       tmp=vunit1(1)*vunit2(2)-vunit2(1)*vunit1(2)  ! cross product


       if(tmp>0.d0) then  ! point on the left side of the edge 
          if(width>0.d0) then
             if(convex2d_dist(set1(i:i,:2),hull2(:nh2,:2))>width) then  ! boundary tolerance
                inside=.false.
                exit
             end if
          else
                inside=.false.
                exit
          end if              
       end if
    end do
    if(inside) then
       numb=numb+1
       set3(numb,:)=set1(i,:)
    end if
  end do
  end if

  ! if has common region, solve edge intersection
  k=numb
  n_intersect=0
   do i=1,nh1
        i2=i+1
        if(i==nh1) i2=1
        do j=1,nh2
          j2=j+1
          if(j==nh2) j2=1
          xa=hull1(i,1); xb=hull1(i2,1); ya=hull1(i,2); yb=hull1(i2,2)
          xc=hull2(j,1); xd=hull2(j2,1); yc=hull2(j,2); yd=hull2(j2,2)
          delta=-(xb-xa)*(yd-yc)+(yb-ya)*(xd-xc)
          if(delta/=0.d0) then   ! intersect
             lambda=(-(xc-xa)*(yd-yc)+(yc-ya)*(xd-xc))/delta
             mu=((xb-xa)*(yc-ya)-(yb-ya)*(xc-xa))/delta
             if(lambda>=0.d0 .and. lambda<=1.d0 .and. mu>=0.d0 .and. mu<=1.d0) then
                k=k+1
                set3(k,:)=(/xa+lambda*(xb-xa),ya+lambda*(yb-ya)/)
             end if
             if(lambda>0 .and. lambda<1.d0 .and. mu>0 .and. mu<1.d0) n_intersect=n_intersect+1
          end if
        end do
   end do
   if(k>=3) then
      area=convex2d_area(set3(:k,:))   ! there must be intersection if there is overlap
   else if (k==0) then
      area= -convex2d_dist(hull1(:nh1,:),hull2(:nh2,:)) ! set a negative value
   end if
   if(n_intersect>2) numb=ns1+ns2  ! presumed totally overlap if there are more than 2 intersect points

end if
end subroutine


function convex2d_area(set)
! calculate the area of a 2d convex hull
! input: set, a matrix
! output: area
real*8 va(2),vb(2),area,convex2d_area,set(:,:),hull(ubound(set,1),ubound(set,2))
integer i,j,k,nh

call convex2d_hull(set,nh,hull)

area=0.d0
do i=3,nh   ! number of vertices must be >=3
  j=i-1
  va=(/hull(j,1)-hull(1,1),hull(j,2)-hull(1,2)/)
  vb=(/hull(i,1)-hull(1,1),hull(i,2)-hull(1,2)/)
  area=area+0.5d0*abs(va(1)*vb(2)-va(2)*vb(1))
end do
convex2d_area=area
if(area<1d-10) area=0.d0

end function


function convex2d_in(set,point,width)
! check if a point is inside the hull
! input: the data set forming the hull; the point; and boundary tolerence
! output: .false. or .true.

integer i,j,k,nh
real*8 point(2),set(:,:),hull(ubound(set,1),ubound(set,2)),tmp,point2(1,2),width,norm1,norm2,vunit1(2),vunit2(2)
logical convex2d_in

call convex2d_hull(set,nh,hull)
convex2d_in=.true.
point2(1,:)=point

if(nh==1) then
  if(sqrt(sum((point-set(1,:))**2))>1d-10) then
      convex2d_in=.false.
  end if
  return
end if

do j=1,nh
   k=j+1
   if(j==nh) k=1

    norm1=sqrt(sum((hull(k,:)-hull(j,:))**2))
    if(norm1<1d-10) cycle
    vunit1=(hull(k,:)-hull(j,:))/norm1
    norm2=sqrt(sum((point-hull(j,:))**2))
    if(norm2<1d-10) cycle
    vunit2=(point-hull(j,:))/norm2
    tmp=vunit1(1)*vunit2(2)-vunit2(1)*vunit1(2)  ! cross product

   if(tmp>0.d0) then 
        if(width>0.d0) then
           if(convex2d_dist(point2(1:1,:2),hull(:nh,:2))>width) then  ! boundary tolerance
              convex2d_in=.false.
              exit
           end if
        else
              convex2d_in=.false.
              exit             
        end if  
   end if
end do

end function

function convex2d_dist(set1,set2)
! calculate the distance between two data sets (convex hulls)

real*8 set1(:,:),set2(:,:),convex2d_dist,hull1(ubound(set1,1),2),hull2(ubound(set2,1),2),&
p1(2),p2(2),p3(2),va(2),vb(2),vv(2),dot,len_sq,param,dist
integer nh1,nh2,i,j,k

if(ubound(set1,1)==1) then
  nh1=1
  hull1=set1
else
  call convex2d_hull(set1,nh1,hull1)
end if

if(ubound(set2,1)==1) then
  nh2=1
  hull2=set2
else
  call convex2d_hull(set2,nh2,hull2)
end if

convex2d_dist=1d10

do i=1,nh1
   if(nh1==1) exit
   j=i+1
   if(i==nh1) j=1
   p1=hull1(i,:)
   p2=hull1(j,:)
   va=p2-p1
   len_sq=sum(va**2)
   do k=1,nh2     ! distance of data from set2 to hull1
      p3=hull2(k,:)
      vb=p3-p1
      dot=sum(va*vb)
      if(len_sq/=0.d0) then
         param=dot/len_sq
      else
         param=-1
      end if

      if(param<0.d0) then
         vv=p1
      else if (param>1.d0) then
         vv=p2
      else
        vv(1)=p1(1)+param*va(1)
        vv(2)=p1(2)+param*va(2)
      end if
      dist=sqrt(sum((p3-vv)**2))
      if(dist<convex2d_dist) convex2d_dist=dist
   end do
end do

do i=1,nh2
   if(nh2==1) exit
   j=i+1
   if(i==nh2) j=1
   p1=hull2(i,:)
   p2=hull2(j,:)
   va=p2-p1
   len_sq=sum(va**2)
   do k=1,nh1
      p3=hull1(k,:)
      vb=p3-p1
      dot=sum(va*vb)
      if(len_sq/=0.d0) then
         param=dot/len_sq
      else
         param=-1
      end if

      if(param<0.d0) then
         vv=p1
      else if (param>1.d0) then
         vv=p2
      else
        vv(1)=p1(1)+param*va(1)
        vv(2)=p1(2)+param*va(2)
      end if
      dist=sqrt(sum((p3-vv)**2))
      if(dist<convex2d_dist) convex2d_dist=dist
   end do
end do
end function


subroutine convex1d_overlap(set1,set2,width,numb,length)
! input: set1, set2, width
! output: number of data, and the length, in the overlapped region
real*8 set1(:),set2(:),length,mini,maxi,width
integer numb,i,ns1,ns2
ns1=ubound(set1,1)
ns2=ubound(set2,1)

numb=0
length=0.d0

mini=minval(set2); maxi=maxval(set2)
do i=1,ns1
 if(set1(i)>=mini-width .and. set1(i)<=maxi+width) numb=numb+1
end do

mini=minval(set1); maxi=maxval(set1)
do i=1,ns2
 if(set2(i)>=mini-width .and. set2(i)<=maxi+width) numb=numb+1
end do

length=(min(maxval(set1),maxval(set2))-max(minval(set1),minval(set2)))
! if overlapped when length>0, and separated distance when length<0.

end subroutine


end module

