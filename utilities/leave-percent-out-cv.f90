program kfoldcv
! creating train.dat for leave-percent-out cross validation 

integer i,j,k,l,niter,nsample
character jobname*6,line*100000
real rand

! USER input
parameter(niter=50,nsample=192,percent=0.1)
logical selected(nsample)

call random_seed()

!----------------------

do i=1,niter

  k=0
  selected=.false.
  l=nsample*percent

  do while(k<l) 
    call random_number(rand)
    j=ceiling(rand*nsample)
    if(selected(j)==.false.) then
       selected(j)=.true.
       k=k+1
    end if
  end do

  jobname(1:3)='sub'
  write(jobname(4:6),'(i3.3)'),i
  call system('mkdir '//trim(jobname)//'')
  
  open(1,file='train.dat',status='old')
  open(2,file=jobname//'/train.dat',status='replace')
  open(3,file=jobname//'/predict.dat',status='replace')
  open(4,file=jobname//'/rand.dat',status='replace')
  read(1,'(a)'),line
  write(2,'(a)'),trim(line)
  write(3,'(a)'),trim(line)

  do j=1,nsample
  read(1,'(a)'),line
  if(selected(j)) then
   write(2,'(a)'),trim(line)  ! train.dat
   write(4,'(i5)'),j
  else  
   write(3,'(a)'),trim(line)  ! predict.dat
  end if
  end do

  close(1)
  close(2)
  close(3)
  close(4)

end do


end program

