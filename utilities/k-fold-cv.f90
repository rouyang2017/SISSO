program kfoldcv
! creating the train.dat and SISSO.in for the subsets of samples for  k-fold cross validation 
! applicable to both classification and regression

integer i,j,k,l,iii,kfold,ptype,ngroup
integer,allocatable:: nsample(:),msample(:)
character jobname*7,line*100000,nsample_line*100000
real rand
logical,allocatable:: selected(:,:),selected_this(:,:)

! User input
 parameter(kfold=10)    


! read SISSO.in for ptype, nsample
open(1,file='SISSO.in',status='old')
do while (.not. eof(1))
   read(1,'(a)') line
   i=index(line,'!')
   if(i/=0) line(i:)=''
   if(index(line,'ptype')/=0) then
      j=index(line,'=')
      read(line(j+1:),*) ptype
   else if(index(line,'nsample')/=0) then
      j=index(line,'=')
      read(line(j+1:),'(a)') nsample_line
   end if
end do
close(1)

if(ptype==1) then
   allocate(nsample(1))
   read(nsample_line,*) nsample(1)
   ngroup=1
else if(ptype==2) then
   i=index(nsample_line,'(')
   j=index(nsample_line,')')
   l=0
   do k=i,j
      if(nsample_line(k:k)==',') l=l+1
   end do
   allocate(nsample(l+1))
   read(nsample_line(i+1:j-1),*) nsample(1:l+1)
   ngroup=l+1
end if

allocate(msample(ngroup))
allocate(selected(ngroup,maxval(nsample)))
allocate(selected_this(ngroup,maxval(nsample)))

call random_seed()

!----------------------

selected=.false.

do i=1,kfold

  selected_this=.false.
  do iii=1,ngroup  ! randomly selecting training data for this fold
     k=0
     msample(iii)=nsample(iii)/kfold  
     if(i<=mod(nsample(iii),kfold)) msample(iii)=msample(iii)+1
     do while(k<msample(iii)) 
       call random_number(rand)
       j=ceiling(rand*nsample(iii))
       if(selected(iii,j)==.false.) then
          selected(iii,j)=.true.
          selected_this(iii,j)=.true.
          k=k+1
       end if
     end do
  end do


  jobname(1:4)='iter'
  write(jobname(5:7),'(i3.3)'),i
  call system('mkdir '//trim(jobname)//'')
  
  open(1,file='train.dat',status='old')
  open(2,file=jobname//'/train.dat',status='replace')
  open(3,file=jobname//'/predict.dat',status='replace')
  open(4,file=jobname//'/rand.dat',status='replace')
  read(1,'(a)'),line
  write(2,'(a)'),trim(line)
  write(3,'(a)'),trim(line)

  do iii=1,ngroup
     do j=1,nsample(iii)
       read(1,'(a)'),line
       if(.not. selected_this(iii,j)) then
        write(2,'(a)'),trim(line)  ! train.dat
       else  
        write(3,'(a)'),trim(line)  ! predict.dat
        if (ptype==1) then
             write(4,'(i5)'),j
        else if (ptype==2) then
             write(4,'(i5,a,i5)'),j,'  in group ',iii
        end if
       end if
     end do
  end do

  close(1)
  close(2)
  close(3)
  close(4)

  open(1,file='SISSO.in',status='old')
  open(2,file=jobname//'/SISSO.in',status='replace')
  do while(.not. eof(1))
     read(1,'(a)') line
     j=index(line,'!')
     if(j/=0) line(j:)=''
     if(line=='') cycle
     if(index(line,'nsample')==0) then
       write(2,'(a)') trim(line)
     else
       if(ptype==1) then
         write(line,'(a,i5)') 'nsample = ',nsample(1)-msample(1)
       else if(ptype==2) then
         write(line,'(a,<ngroup-1>(i5,a),i5,a)') &
              'nsample = (',((nsample(k)-msample(k),','),k=1,ngroup-1),nsample(ngroup)-msample(ngroup),')'
       end if
       write(2,'(a)') trim(line)
     end if
  end do

  close(1)
  close(2)

end do


deallocate(nsample)
deallocate(msample)
deallocate(selected)
deallocate(selected_this)


end program

