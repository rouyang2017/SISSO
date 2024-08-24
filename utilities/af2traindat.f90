program af2traindat
! Creating the train.dat based on the files atom_features (list of atomic features) and
! samplelist (list of the training samples). 

real*8,allocatable:: feat(:,:)
character line*1000,atomline*1000,strings(500)*50,letter(14)*2
character,allocatable:: featname(:,:)*50,element(:)*10
integer i,j,k,k1,k2,k3,k4,k5,k6,natom,naf,nele,npair,length

! mark for the constituent elements of a compound
letter=(/'_A','_B','_C','_D','_E','_F','_G','_H','_I','_J','_K','_L','_M','_N'/)

! files atom_features and samplelist must be ready
write(*,'(a)') 'Number of atomic features to input from the file "atom_features"'
read(*,*) naf
write(*,'(a)'),'Number of elements in each sample in the file "samplelist" '
read(*,*),nele

allocate(feat(nele,naf))
allocate(featname(nele,naf))
allocate(element(nele))

open(1,file='train.dat',status='replace')
open(2,file='samplelist',status='old')
open(3,file='atom_features',status='old')

! read the title line from atom_features for the feature names
read(3,'(a)'),line
call string_split(line,strings,' ')
do i=1,naf
do j=1,nele
featname(j,i)=trim(adjustl(strings(i+1)))//letter(j) 
end do
end do
rewind(3)

! get the element-name length
read(2,*)  ! skip the title line
read(2,'(a)'),line
line=adjustl(line)
length=(index(line,' ')-1)/nele
rewind(2)
! write the title line to the train.dat
write(1,'(a<max(8,length*nele)>,<nele*naf>a15)'),'compound',(((trim(adjustl(featname(k,i)))),k=1,nele),i=1,naf)

! creating the file train.dat
read(2,*)  ! skip the title line
do while(.not. eof(2))
read(2,'(a)'),line
   ! get the compound name
   call string_split(line,strings,' ')
   strings(1)=trim(adjustl(strings(1)))
   do i=1,nele
   element(i)=strings(1)((i-1)*length+1:i*length)
   end do

   ! get the feature data from atom_features
   do i=1,nele
     rewind(3)
     100 read(3,'(a)',end=101),atomline
     k=index(atomline,trim(adjustl(element(i))))
     if(k/=0) then
        read(atomline(k+length:),*),feat(i,:naf)
        cycle
     else 
        goto 100
     end if
     101 continue
     print *,'atom not found: ', trim(adjustl(element(i)))
     stop
   end do     


   ! write the feature data
   write(1,'(<max(0,8-nele*length)>X,<nele>a,<nele*naf>f15.5)'),(trim(adjustl(element(i))),i=1,nele),&
         ((feat(j,i),j=1,nele),i=1,naf)
end do

close(1)
close(2)
close(3)

deallocate(feat)
deallocate(featname)
deallocate(element)

contains


subroutine string_split(instr,outstr,sp)
! break a string into sub-strings
! input: instr, input string; sp, separator
! output: outstr, output sub-strings
character(len=*) instr,outstr(:),sp
integer n,i,j
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


end


