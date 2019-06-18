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
! Subroutines/Functions:
! ----------
! det: calculates the determinant of a given matrix
! inverse: calculates the inverse of a given matrix
! orth_de: linear least square fitting by orthogonal decomposition
! worth_de: observation weighted orth_de
! lls: linear least fitting by conventional method (not robust)
! coord_descent: coordinate descent
! crosspro: crossproduct between two vectors
! LASSO: least absolute shrinkage and selection operator
! mtlasso_mpi Multi-task LASSO
! corr:  pearson's correlation
! string_split: break a string into sub-strings.
! kfoldCV: k-fold CV for linear model
! convex2d_xxx: construction of 2D convex hull
! intriangle: check if a 3D point is inside a triangle
! convex3d_xxx: construction of 3D convex hull
! interlp: intersection of a line and a plane
! dispp: distance from a point to a plane
!************************************************************************

use var_global

contains




function dispp(p,a,b,c)
! Distance of a point from a plane
! input: a,b,c are the three points determing the plane, p is the point outside of the plane.
! output: the distance

real*8 dispp,a(3),b(3),c(3),p(3),normal(3)
normal=crosspro((a-b),(c-b))
dispp=abs(sum(normal*(p-a)))/sqrt(sum(normal**2))

end function

function interlp(p1,p2,p3,p4,p5)
! intersection between a line and a plane
! input: p1,p2,p3 for determing the plane, p4 and p5 determing the line
! output: the interception point

real*8 interlp(3),p1(3),p2(3),p3(3),p4(3),p5(3),normal(3),a,b,c,d,&
       m,n,p,t
! for plane
normal=crosspro((p1-p2),(p3-p2))
a=normal(1)
b=normal(2)
c=normal(3)
d=-(a*p1(1)+b*p1(2)+c*p1(3))  ! aX+bY+cZ+d=0 (X,Y,Z) is the normal
! for line
m=p5(1)-p4(1)
n=p5(2)-p4(2)
p=p5(3)-p4(3)  !line equation: (x-x0)/m=(y-y0)/n=(z-z0)/p=t
! intersection
t=-(a*p4(1)+b*p4(2)+c*p4(3)+d)/(m*a+n*b+p*c) ! X=mt+x0; Y=nt+y0; Z=pt+z0
interlp(1)=m*t+p4(1)
interlp(2)=n*t+p4(2)
interlp(3)=p*t+p4(3)
end function


function intriangle(p,a,b,c)
! p,a,b,c are in the same plane
! check if point p is inside the triangle formed by points a,b,c
real*8 a(3),b(3),c(3),p(3),v1(3),v2(3),normal_0(3),normal_1(3),normal_2(3),normal_3(3)
logical intriangle
intriangle=.false.
v1=a-b
v2=c-b
normal_0=crosspro_abnormalized(v1,v2) ! triangle normal
v1=p-b
v2=c-b
normal_1=crosspro_abnormalized(v1,v2)
v1=p-c
v2=a-c
normal_2=crosspro_abnormalized(v1,v2)
v1=p-a
v2=b-a
normal_3=crosspro_abnormalized(v1,v2)
if(sum(normal_0*normal_1)>=0.d0 .and. sum(normal_0*normal_2)>=0.d0 .and. sum(normal_0*normal_3)>=0.d0) intriangle=.true.
end function


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
! linear least square fit to y=a+x*b by orthogonal decomposition
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

subroutine orth_de_nointercept(x,y,beta,rmse)
! linear least square fit to y=x*b by orthogonal decomposition
! input: matrix x,vector y; 
! output beta,rmse

real*8 x(:,:),y(:),beta(ubound(x,2)),rmse
real*8 q(ubound(x,1),ubound(x,2)),r(ubound(x,2),ubound(x,2)),qty(ubound(x,2))
integer i,j,k,m,n

m=ubound(x,1)
n=ubound(x,2)

call qr_de(x,q,r) ! results stored in q and r
qty=matmul(transpose(q),y)
beta(n)=qty(n)/r(n,n)

do i=n,1,-1
beta(i)=qty(i)
do j=i+1,n
beta(i)=beta(i)-r(i,j)*beta(j)
end do
beta(i)=beta(i)/r(i,i)
end do

rmse=sqrt(sum((y-matmul(x,beta))**2)/m)

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

function crosspro_abnormalized(a,b)
! a and b are normalized before calculating their cross product
real*8  a(3),b(3),aa(3),bb(3),crosspro_abnormalized(3),norma,normb
aa=a
bb=b
norma=sqrt(sum(aa**2))
normb=sqrt(sum(bb**2))
if(norma/=0.d0 .and. normb/=0.d0) then
  aa=aa/norma
  bb=bb/normb
else
  crosspro_abnormalized=0.d0
end if

crosspro_abnormalized(1)=aa(2)*bb(3)-aa(3)*bb(2);
crosspro_abnormalized(2)=aa(3)*bb(1)-aa(1)*bb(3);
crosspro_abnormalized(3)=aa(1)*bb(2)-aa(2)*bb(1);

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
! input: set, a matrix N x 2
! output: numb (number of vertices); hull, the vertices stored in clockwise direction
real*8 set(:,:),hull(:,:),tmp,vjj(2),vij(2),vkj(2),normij,normkj
integer numb,i,j,k,ntot,loc(1),nrecord
logical used(ubound(set,1)),isvertex

ntot=ubound(set,1)
used=.false.
numb=0

! find the initial point with min(x) and min(y)
! j can't be set 'used' as it is waiting to be found to close the hull
loc=minloc(set(:,1)) 
j=loc(1)
do i=1,ntot  
  if(i==j) cycle
  if(abs(set(loc(1),1)-set(i,1))<=1d-10 .and. set(loc(1),2)>set(i,2)) loc(1)=i  
end do
j=loc(1)
numb=1
hull(1,:)=set(j,:)

! find points at the same position with j
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
! from known vector vjj, find the next vertex k by searching all i

vjj=(/0.d0,-1.d0/)  ! initial vector pointing downward
nrecord=0

do while(any(used .eqv. .false.) )
  nrecord=nrecord+1
  if(nrecord>ntot) exit ! in case of endless loop
  isvertex=.false.
  vkj=vjj
  normkj=-1.d0  ! normkj<0 if vkj==vjj; normkj>0 if vkj\=vjj

  ! find the next vertex k to update j, vjj, numb, hull, used
  do i=1,ntot

    if(used(i) .or. i==j) cycle   ! loc(1) will be the last point to find

    if(abs(set(j,1)-set(i,1))<=1d-10 .and. abs(set(j,2)-set(i,2))<=1d-10  ) then  ! all the points at the same position 
       numb=numb+1
       hull(numb,:)=set(i,:)
       used(i)=.true.
       cycle
    end if

    vij=set(i,:)-set(j,:)     
    normij=sqrt(sum(vij**2))
    vij=vij/normij
    if(i==loc(1)) normij=1d10  ! ensure loc(1) is the last vertex 

    tmp=vkj(1)*vij(2)-vij(1)*vkj(2)  ! cross product between vector vij and vkj

    if (tmp>1d-10   .or. &  ! on the left side of vkj(vkj==vjj) or vij
       (normkj>0.d0 .and. tmp>0.d0 .and. tmp<=1d-10 .and. sum(vij*vkj)>0.d0 ) .or. &  ! vkj/=vjj
       (normkj<0.d0 .and. tmp>=0.d0 .and. tmp<=1d-10 .and. sum(vij*vjj)<0.d0 ) .or. &  ! vkj==vjj
       (normkj>0.d0 .and. tmp==0.d0  .and. sum(vij*vkj)>0.d0 .and. normij<normkj ) ) then 
          k=i
          vkj=vij
          normkj=normij  ! positive
          isvertex=.true.
    end if
  end do

  if(isvertex) then  ! assign the new vertex k to j
    if(k==loc(1)) exit  ! loc(1) is already in the hull(1,:)
     numb=numb+1
     hull(numb,:)=set(k,:)
     used(k)=.true.
     vjj=-vkj   ! reverse the direction of vkj and assign to vjj
     j=k   ! new j
  end if

end do

end subroutine


subroutine convex3d_hull(set,ntri,triangles)
! find the convex hull by finding all the triangle facets
! input: set, a matrix N x 3
! output:  triangles, the ids of vertices forming traingles; ntri, number of triangles
! NOTE: triangles is an allocatable array
real*8  set(:,:),tmp,vij(2),vkj(2),normal_a(3),normal_b(3),vtmp(3,3),vedge(3),aa,bb,cc
integer ntri,i,j,k,l,m,ntot,loc(2),nedge,iedge,leftright,comp(3),size_tri,size_edge,ntri_coplane,ntri_saved,nvert
integer,allocatable:: triangles(:,:),edges(:,:),change(:,:),recorder(:),coplane(:,:),ver_index(:)
logical overlap(ubound(set,1)),used,planar,goodedge

ntot=ubound(set,1)
overlap=.false.  ! points at the same position
ntri=0
nedge=0

size_tri=10000
size_edge=10000
allocate(triangles(size_tri,3))  ! size to be adjusted
allocate(edges(size_edge,2))

! find the initial edge (the most left one projected on the x-y plane)

!-------------------------------------
! first point: min(x),min(y),min(z)
loc(1:1)=minloc(set(:,1)) 
k=loc(1)
do i=1,ntot  
  if(i==k) cycle
  if(abs(set(loc(1),1)-set(i,1))<=1d-10 .and. set(loc(1),2)>set(i,2)) loc(1)=i  
end do
k=loc(1)
do i=1,ntot
 if(i==k) cycle
 if(abs(set(loc(1),1)-set(i,1))<=1d-10 .and. abs(set(loc(1),2)-set(i,2))<=1d-10 .and. set(loc(1),3)>set(i,3)) loc(1)=i
end do

! second point for the initial edge
vkj=(/0.d0,-1.d0/)
vtmp(3,:)=(/0.d0,0.d0,-1.d0/)
loc(2)=0
do i=1,ntot   ! find a i so that vector i<-j is to the left side of k<-j
  if(i==loc(1)) cycle 
  vij=set(i,:2)-set(loc(1),:2)     
  tmp=vkj(1)*vij(2)-vij(1)*vkj(2)  ! cross product vkj x vij

  if (tmp>1d-10 ) then  ! on the left side of vkj or vij
     loc(2)=i
     vkj=vij
     vtmp(3,:)=set(i,:3)-set(loc(1),:3) ! 3D form of vkj
  else if(abs(tmp)<=1d-10) then  ! considering 3D
    vtmp(1,:)=(/0.d0,0.d0,1.d0/)
    vtmp(2,:)=set(i,:3)-set(loc(1),:3)
    aa=sum(crosspro_abnormalized(vtmp(2,:),vtmp(1,:))*crosspro_abnormalized(vtmp(2,:),vtmp(3,:)))
    bb=sum(crosspro_abnormalized(vtmp(2,:),vtmp(1,:))**2)  ! between vi and z
    cc=sum(crosspro_abnormalized(vtmp(2,:),vtmp(3,:))**2)  ! between vi and vk
    if(aa<0 .or. (bb==0.d0 .and. cc/=0.d0) .or. (cc==0.d0 .and. sum(vtmp(2,:)**2)>sum(vtmp(3,:)**2)) ) then
        loc(2)=i
        vkj=vij
        vtmp(3,:)=set(i,:3)-set(loc(1),:3) ! 3D form of vkj
    end if
  end if
end do
if(loc(2)==0) return  ! all points are at the same one point
!-------------------------------------------------------

nedge=1
ntri=0
if(loc(1)<loc(2)) then ! let edges(i,1) always < edges(i,2)
  edges(1,:)=(/loc(1),loc(2)/)
else
  edges(1,:)=(/loc(2),loc(1)/)
end if

! find all points at the same position with the initial edge points
do i=1,ntot
   if(i==edges(1,1) .or. i==edges(1,2)) cycle
   if( sqrt(sum((set(edges(1,1),:)-set(i,:))**2))<1d-10 .or. sqrt(sum((set(edges(1,2),:)-set(i,:))**2))<1d-10) &
      overlap(i)=.true.
end do

! With the initial edge and plane, find all the triangles
leftright=1
iedge=1
do while(iedge<=nedge)

  ! set a initial plane
   vedge=set(edges(iedge,2),:)-set(edges(iedge,1),:)
   goodedge=.true.

   do i=1,ntot
     if(i==edges(iedge,1) .or. i==edges(iedge,2) .or. overlap(i) ) cycle
     vtmp(1,:)=set(i,:)-set(edges(iedge,1),:)
     vtmp(2,:)=crosspro_abnormalized(vedge,vtmp(1,:)) ! normal
     if(sqrt(sum(vtmp(2,:)**2))>1d-10) then
       normal_a=vtmp(2,:)      ! the normal vector of the initial plane
       k=i
       exit
     end if
     if(i==ntot) return  ! all points are on a line
   end do

 ! check if all the points are in one plane
   if(iedge==1) then
     planar=.true.
     do i=1,ntot
         if( sum((set(i,:)-set(edges(iedge,1),:))*normal_a) /=0.d0) then
           planar=.false.
           exit
         end if 
     end do
     if(planar) return
   end if

 !-------------------------------------
 ! For each edge, find the leftmost and rightmost facet!
  do i=1,ntot
    if(overlap(i) .or. k==i .or. edges(iedge,1)==i .or. edges(iedge,2)==i) cycle 
    vtmp(1,:)=set(i,:)-set(edges(iedge,1),:) 
    vtmp(2,:)=crosspro_abnormalized(vedge,vtmp(1,:)) ! normal of the considered plane
    normal_b=vtmp(2,:)
    if(sqrt(sum(normal_b**2))<1d-10) cycle  ! point i is on the edge-line
    vtmp(3,:)=crosspro_abnormalized(normal_a,normal_b) ! cross product between the normal of two plane
    tmp=sum(vtmp(3,:)*vedge) ! to the left or right

    if(abs(tmp)>1d-10) then  !/=0
       if((leftright>0 .and. tmp>0.d0) .or. (leftright<0 .and. tmp<0.d0)) then
          k=i
          normal_a=normal_b
       end if
    else if(.not. intriangle(set(i,:),set(edges(iedge,1),:),set(edges(iedge,2),:),set(k,:)) ) then
       if(leftright>0) then
         vtmp(1,:)=set(i,:)-set(edges(iedge,1),:)
         vtmp(2,:)=set(k,:)-set(edges(iedge,1),:)
         aa=sum(crosspro_abnormalized(vedge,vtmp(2,:))*crosspro_abnormalized(vtmp(2,:),vtmp(1,:))) !edge x vk,vk x vi
         bb=sum(crosspro_abnormalized(vedge,vtmp(2,:))*crosspro_abnormalized(vedge,vtmp(1,:))) 
       else ! reverse the edge direction
         vtmp(1,:)=set(i,:)-set(edges(iedge,2),:)
         vtmp(2,:)=set(k,:)-set(edges(iedge,2),:)
         aa=sum(crosspro_abnormalized(-vedge,vtmp(2,:))*crosspro_abnormalized(vtmp(2,:),vtmp(1,:)))!edge x vk,vk x vi
         bb=sum(crosspro_abnormalized(-vedge,vtmp(2,:))*crosspro_abnormalized(-vedge,vtmp(1,:)))
       end if

       if( bb <-1d-10) then   ! the edge is not really edge
          goodedge=.false.
          exit
       end if
  
       if( aa >=0.d0) then ! the edge repel the 'i' to the other edge
          k=i
          normal_a=normal_b
       end if

    end if
  end do

 !---------------------------------------

  if(goodedge) then
      ! add new triangle
       comp(1)=min(edges(iedge,1),edges(iedge,2),k)
       comp(3)=max(edges(iedge,1),edges(iedge,2),k)
       if(k/=comp(1) .and. k/=(comp(3))) then
           comp(2)=k
       elseif(edges(iedge,1)/=comp(1) .and. edges(iedge,1)/=comp(3)) then
           comp(2)=edges(iedge,1)
       else
           comp(2)=edges(iedge,2)
       end if  

       used=.false.
       do i=1,ntri
         if(all(triangles(i,:) == comp(:))) then
            used=.true.
            exit
         end if
       end do

       if(.not. used ) then
          ntri=ntri+1
          triangles(ntri,:3)=comp(:3)
       end if
       if(ntri==size_tri) then
         allocate(change(size_tri,3))
         change=triangles
         deallocate(triangles)
         size_tri=size_tri*2
         allocate(triangles(size_tri,3))
         triangles(:ntri,:)=change
         deallocate(change)
       end if

       ! find overlap points
       do i=1,ntot
          if(i==k .or. overlap(i) ) cycle
          if( (abs(set(k,1)-set(i,1))<=1d-10 .and. abs(set(k,2)-set(i,2))<=1d-10 .and. &
              abs(set(k,3)-set(i,3))<=1d-10 )) overlap(i)=.true.
       end do

       do j=1,2
           comp(1)=min(edges(iedge,j),k)
           comp(2)=max(edges(iedge,j),k)
           used=.false.
           do i=1,nedge
              if(all(edges(i,:)==comp(:2))) then
                  used=.true.
                  exit
              end if
           end do
           if(.not. used) then
              nedge=nedge+1
              edges(nedge,:)=comp(:2)
           end if

          if(nedge==size_edge) then
            allocate(change(size_edge,2))
            change=edges
            deallocate(edges)
            size_edge=size_edge*2
            allocate(edges(size_edge,2))
            edges(:nedge,:)=change
            deallocate(change)
          end if

       end do
   end if

   leftright=-leftright   ! each edge is to be used twice.
   if(leftright>0) iedge=iedge+1

end do


! rearrangement of the triangles in the same plane
allocate(change(5*ntri,3))
allocate(coplane(5*ntri,3))
allocate(recorder(ntri))
allocate(ver_index(ntri*3))
recorder=1  ! 1: unaccessed; 0: accessed
ntri_saved=0  ! number of rearranged triangles
ntri_coplane=0  ! number of triangles in one plane

do i=1,ntri
  if(recorder(i)==0) cycle
  recorder(i)=0
  ntri_coplane=1
  coplane(ntri_coplane,:)=triangles(i,:)
  vtmp(1,:)=set(triangles(i,2),:)-set(triangles(i,1),:)
  vtmp(2,:)=set(triangles(i,3),:)-set(triangles(i,1),:)
  normal_a=crosspro_abnormalized(vtmp(1,:),vtmp(2,:)) ! the reference normal
  do j=i+1,ntri  ! find the coplane triangles if there are
     if(recorder(j)==0) cycle
     vtmp(1,:)=set(triangles(j,2),:)-set(triangles(j,1),:)
     vtmp(2,:)=set(triangles(j,3),:)-set(triangles(j,1),:)
     normal_b=crosspro_abnormalized(vtmp(1,:),vtmp(2,:)) !the normal of candidate triangle
     if(sqrt(sum(crosspro_abnormalized(normal_a,normal_b)**2))<=1d-10 .and. &
        dispp(set(triangles(j,1),:),set(triangles(i,1),:),set(triangles(i,2),:),set(triangles(i,3),:))<1d-6) then 
        ntri_coplane=ntri_coplane+1
        coplane(ntri_coplane,:)=triangles(j,:)
        recorder(j)=0
     end if
  end do
 
  ver_index=0
  ver_index(1:3)=coplane(1,:) 
  nvert=3 
  do j=2,ntri_coplane  ! find the vertice of all the coplane triangles
    if(all(coplane(j,1)/=ver_index(:nvert))) then
       nvert=nvert+1; ver_index(nvert)=coplane(j,1)
    end if
    if(all(coplane(j,2)/=ver_index(:nvert))) then
       nvert=nvert+1; ver_index(nvert)=coplane(j,2)
    end if
    if(all(coplane(j,3)/=ver_index(:nvert))) then
       nvert=nvert+1; ver_index(nvert)=coplane(j,3)
    end if
  end do

  if(nvert>3) then
  ! direction
    vtmp(1,:)=set(ver_index(2),:)-set(ver_index(1),:)
    vtmp(2,:)=set(ver_index(3),:)-set(ver_index(1),:)
    normal_a=crosspro_abnormalized(vtmp(1,:),vtmp(2,:))  

  ! find the starting point
    k=1
    do j=3,nvert
       vtmp(1,:)=set(ver_index(2),:)-set(ver_index(k),:)
       vtmp(2,:)=set(ver_index(j),:)-set(ver_index(k),:)
       vtmp(3,:)=crosspro_abnormalized(vtmp(1,:),vtmp(2,:))
       if(sum(vtmp(3,:)**2)==0.d0 .and. sum(vtmp(1,:)*vtmp(2,:))<0.d0) k=j
    end do
    ! switch points 1 and k
    j=ver_index(1)
    ver_index(1)=k
    ver_index(k)=j

   ! find the second point
    k=2
    do j=3,nvert  !find the first most left edge: 1-k
       vtmp(1,:)=set(ver_index(k),:)-set(ver_index(1),:)
       vtmp(2,:)=set(ver_index(j),:)-set(ver_index(1),:)
       normal_b=crosspro_abnormalized(vtmp(1,:),vtmp(2,:))
       aa=sum(normal_b*normal_a)
       if( aa>0.d0 .or. (aa==0.d0 .and.  sum(vtmp(2,:)**2)>sum(vtmp(1,:)**2) ) ) k=j
    end do
    
    ! find all the triangles
    m=2  ! point 1 and k
    vtmp(1,:)=set(ver_index(1),:)-set(ver_index(k),:) 
    do while(m<nvert)
       l=0
       do j=2,nvert  ! find triangle 1-l-k
         if(j==k .or. ver_index(j)==0) cycle ! skip the evaluated points
         vtmp(2,:)=set(ver_index(j),:)-set(ver_index(k),:)
         normal_b=crosspro_abnormalized(vtmp(1,:),vtmp(2,:))
         aa=sum(normal_b*normal_a)
         if( aa> 0.d0) then
           l=j
           vtmp(1,:)=set(ver_index(l),:)-set(ver_index(k),:)
         elseif(aa==0.d0) then   ! linear
           if(sum(vtmp(2,:)**2)>sum(vtmp(1,:)**2)) then
             l=j
             vtmp(1,:)=set(ver_index(l),:)-set(ver_index(k),:)
           else
             ver_index(j)=0    ! remove points inside a line segment
           end if
         end if
       end do

       if(l/=0) then
         vtmp(2,:)=crosspro_abnormalized(vtmp(1,:),set(ver_index(l),:)-set(ver_index(1),:))
         if(sqrt(sum(vtmp(2,:)**2))>1d-10) then  ! the three points are not on the same line
            ntri_saved=ntri_saved+1 ! every loop there must be a new point found
            change(ntri_saved,:)=(/ver_index(1),ver_index(k),ver_index(l)/)
         end if
          ver_index(k)=0
          vtmp(1,:)=-vtmp(1,:)  ! reverse the direction for the next loop
          k=l
       end if

       m=m+1  ! every loop find one vertex,totally nvert-2 points
    end do

  else
     ntri_saved=ntri_saved+1
     change(ntri_saved,:)=ver_index(1:3)
  end if    

end do

ntri=ntri_saved
triangles(:ntri_saved,:)=change(:ntri_saved,:)

deallocate(change)
deallocate(recorder)
deallocate(coplane)
deallocate(ver_index)
deallocate(edges)

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


function convex2d_in(set1,set2,width)
! check how many points of set2 are inside the hull of set1
! input: the data set forming the hull; the point; and boundary tolerence
! output: .false. or .true.

integer i,j,k,nh,np,convex2d_in
real*8 set1(:,:),set2(:,:),hull(ubound(set1,1),ubound(set1,2)),tmp,width,norm1,norm2,vunit1(2),vunit2(2)

call convex2d_hull(set1,nh,hull)
convex2d_in=0

if(nh==1) return
np=ubound(set2,1)


do i=1,np
    convex2d_in=convex2d_in+1
    do j=1,nh
       k=j+1
       if(j==nh) k=1
    
        norm1=sqrt(sum((hull(k,:)-hull(j,:))**2))
        if(norm1<1d-10) cycle
        vunit1=(hull(k,:)-hull(j,:))/norm1
        norm2=sqrt(sum((set2(i,:)-hull(j,:))**2))
        if(norm2<1d-10) cycle
        vunit2=(set2(i,:)-hull(j,:))/norm2
        tmp=vunit1(1)*vunit2(2)-vunit2(1)*vunit1(2)  ! cross product
    
       if(tmp>0.d0) then 
            if(width>0.d0) then
               if(convex2d_dist(set2(i:i,:2),hull(:nh,:2))>width) then  ! boundary tolerance
                  convex2d_in=convex2d_in-1
                  exit
               end if
            else
                  convex2d_in=convex2d_in-1
                  exit             
            end if  
       end if
    end do
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

mini=minval(set2)-width
maxi=maxval(set2)+width
do i=1,ns1
 if(set1(i)>=mini .and. set1(i)<=maxi) numb=numb+1
end do

mini=minval(set1)-width
maxi=maxval(set1)+width
do i=1,ns2
 if(set2(i)>=mini .and. set2(i)<=maxi) numb=numb+1
end do

length=(min(maxval(set1),maxval(set2))-max(minval(set1),minval(set2)))
! if overlapped when length>0, and separated distance when length<0.

end subroutine


function convex1d_in(set1,set2,width)
! check how many data points of set2 are inside the segment by set1
real*8 set1(:),set2(:),mini,maxi,width
integer i,convex1d_in
ns1=ubound(set1,1)
ns2=ubound(set2,1)

convex1d_in=0

mini=minval(set1)-width
maxi=maxval(set1)+width
do i=1,ubound(set2,1)
 if(set2(i)>=mini .and. set2(i)<=maxi) convex1d_in=convex1d_in+1
end do

end function


function convex3d_in(set1,set2,width,ntri,triangles)
! check how many points of set2 are inside a 3d convex hull of set1
! input: the full data set1,set2, the boundary tolerence, and the hull-triangles

integer iii,i,j,k,ntri,ninter,triangles(:,:),convex3d_in
real*8 pref(3),set1(:,:),set2(:,:),width,pproj(3),pinter(3),v1(3),v2(3),v3(3),&
       normal(3),ddd,ttt,area,dist(3),vtmp(3,3)
real*8,allocatable:: inter_all(:,:)
logical inside

convex3d_in=0
allocate(inter_all(ntri,3))

do iii=1,ubound(set2,1)
     ninter=0
     pref=(/set2(iii,1)+1.0,set2(iii,2)+1.0,set2(iii,3)+1.0/) ! point-pref(arbitrary) is the line
     inside=.false.
     do i=1,ntri
       v1=set1(triangles(i,2),:)-set1(triangles(i,1),:)
       v2=set1(triangles(i,3),:)-set1(triangles(i,1),:)
       v3=set1(triangles(i,3),:)-set1(triangles(i,2),:)
       normal=crosspro(v1,v2)  ! normal of the triangle plane
       ddd=-sum(normal*set1(triangles(i,1),:)) ! the d of plane: ax+by+cz+d=0
       ttt=-(sum(normal*set2(iii,:))+ddd)/sum(normal**2)  ! the t of line equation
       pproj=set2(iii,:)+normal*ttt  ! projection of the point on the triangle plane
     
       ! measure the distance of the point from the triangle
       if(sqrt(sum((set2(iii,:)-pproj)**2))<=width) then
          if(intriangle(pproj,set1(triangles(i,1),:),set1(triangles(i,2),:),set1(triangles(i,3),:))) then !projection inside
             convex3d_in=convex3d_in+1
             inside=.true.
          else
             vtmp(1,:)=set2(iii,:)-set1(triangles(i,1),:)
             vtmp(2,:)=set2(iii,:)-set1(triangles(i,2),:)
             vtmp(3,:)=set2(iii,:)-set1(triangles(i,3),:)
             if(sum(vtmp(1,:)*v1) <=0.d0 .or. &
                sum(vtmp(2,:)*(-v1)) <=0.d0 ) then 
                 dist(1)=min(sqrt(sum((vtmp(1,:))**2)),sqrt(sum((vtmp(2,:))**2)))
             else ! only acute angles in the triangle
                 area=sqrt(sum((crosspro(vtmp(1,:),v1))**2)) 
                 dist(1)=area/sqrt(sum(v1**2))
             end if
             if(sum(vtmp(1,:)*v2) <=0.d0 .or. &
                sum(vtmp(3,:)*(-v2)) <=0.d0 ) then
                 dist(2)=min(sqrt(sum((vtmp(1,:))**2)),sqrt(sum((vtmp(3,:))**2)))
             else ! only acute angles in the triangle
                 area=sqrt(sum((crosspro(vtmp(1,:),v2))**2))
                 dist(2)=area/sqrt(sum(v2**2))
             end if
             if(sum(vtmp(3,:)*(-v3)) <=0.d0 .or. &
                sum(vtmp(2,:)*v3) <=0.d0 ) then
                 dist(3)=min(sqrt(sum((vtmp(3,:))**2)),sqrt(sum((vtmp(2,:))**2)))
             else ! only acute angles in the triangle
                 area=sqrt(sum((crosspro(vtmp(3,:),(-v3)))**2))
                 dist(3)=area/sqrt(sum(v3**2))
             end if
             if(min(dist(1),dist(2),dist(3))<=width) then
                 convex3d_in=convex3d_in+1
                 inside=.true.
              end if
          end if
       end if
     
       if(inside) exit

       ! if not in the triangle, check intersections
       if(sum((set2(iii,:)-pref)*normal)/=0.d0 ) then ! line-plane is not parallel
         pinter=interlp(set1(triangles(i,1),:),set1(triangles(i,2),:),set1(triangles(i,3),:),set2(iii,:),pref)
         if(intriangle(pinter,set1(triangles(i,1),:),set1(triangles(i,2),:),set1(triangles(i,3),:))) then
             ninter=ninter+1
             inter_all(ninter,:)=pinter
         end if
       end if
     end do
     
     if(inside) cycle

     do i=1,ninter-1
       do j=i+1,ninter
         if(sqrt(sum((inter_all(i,:)-inter_all(j,:))**2))<1d-10) cycle ! the two points are at the same position
         v1=inter_all(i,:)-set2(iii,:)
         v2=inter_all(j,:)-set2(iii,:)
         if(sum(v1*v2)<0.d0) then  ! if v1 and v2 point to different direction, the point is inside
            convex3d_in=convex3d_in+1     ! otherwise the point is outside
            inside=.true.
         end if
         if(inside) exit
       end do
       if(inside) exit
     end do

end do

deallocate(inter_all)

end function

subroutine convex3d_overlap(set1,set2,width,numb)
! counting the number of overlap data
! input: data set 1, data set 2,width(boundary tolerance)
! output: number of data in the overlapped region

integer no1,no2,i,j,k,nset1,nset2,ntri1,ntri2
integer,allocatable:: triangles(:,:)
real*8 set1(:,:),set2(:,:),width
logical inside

nset1=ubound(set1,1)
nset2=ubound(set2,1)

call convex3d_hull(set1,ntri1,triangles)
if(ntri1==0) then  ! <3D are not prefered
   numb=nset1+nset2
   return
end if

! counts data of set2 in set1
no1=0
no1=no1+convex3d_in(set1,set2,width,ntri1,triangles)
deallocate(triangles)  ! allocated in the convex3d_hull


call convex3d_hull(set2,ntri2,triangles)
if(ntri2==0) then  ! <3D are not prefered
   numb=nset1+nset2
   return
end if

! counts data of set1 in set2
no2=0
no2=no2+convex3d_in(set2,set1,width,ntri2,triangles)
deallocate(triangles)

numb=no1+no2

end subroutine


end module


