program relax
implicit none

real(kind=8),parameter :: kb=1.380648813d-23 ! J/K
real(kind=8),parameter :: hbar=1.05457172647d-22 ! J/THz
integer:: nm,im
integer:: nq,iq
real(4),allocatable ::freq(:,:)
integer,allocatable ::weight(:)
real(4):: Tem,M
real(4),allocatable ::gru(:,:),velo(:,:,:),debye2(:),debye(:)
real(4),allocatable ::scat(:,:),tao(:,:)
character(10)::str1,str2,str3,str4,str5
real(4)::Factor1,Factor2
real(4)::factor01,factor02,factor03
!弛豫时间的单位是s

open(unit=103,file="mesh.dat")
open(unit=104,file="gru.dat")
open(unit=105,file="velo.dat")
open(unit=108,file="debye.dat")
open(unit=109,file="relax.dat")
open(unit=106,file="tcond.in")


read(106,*)str1,Tem
read(106,*)
read(106,*)str2,M
read(106,*)str3,nq
read(106,*)str4,nm

!write(*,*)Tem
!write(*,*)M
!write(*,*)nq
!write(*,*)nm
M=M/(6.02d23)    !原子密度改成g/atom
write(*,*)M


allocate(freq(nq,nm))
allocate(weight(nq))
allocate(gru(nq,nm))
allocate(velo(nq,nm,3))
allocate(debye2(nm))
allocate(debye(nm))
allocate(scat(nq,nm))
allocate(tao(nq,nm))



do im=1,nm
   read(108,*)str5,debye2(im) 
   debye(im)=abs(sqrt(debye2(im)))
   !write(*,*)debye(im)
   do iq=1,nq
      read(104,*)gru(iq,im)
      read(105,*)velo(iq,im,:)
      read(103,*)weight(iq),freq(iq,im)
      !write(*,*)freq(iq,im)  
      write(*,*)gru(iq,im)
      
      factor01=3*Tem
      factor02=exp(-debye(im)/factor01)
      !write(*,*)factor02
      factor03=velo(iq,im,1)**2+velo(iq,im,2)**2+velo(iq,im,3)**2 
      
      Factor1=hbar*gru(iq,im)**2*freq(iq,im)**2*Tem*factor02
      Factor2=M*debye(im)*factor03
      
      scat(iq,im)=Factor1/Factor2
      tao(iq,im)=1/scat(iq,im)
      tao(iq,im)=tao(iq,im)*1d-11
      if (freq(iq,im)>0)then
          write(109,'(e14.6)')tao(iq,im)
          else
          write(109,'(e14.6)')0
      endif
   end do
      read(104,*)
      read(105,*)
      read(103,*)
      write(109,*)
end do

   
end
