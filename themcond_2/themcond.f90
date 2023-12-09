program themcond

implicit none

integer:: nm,im
integer:: nq,iq
integer,allocatable ::weight(:)
real(4),allocatable ::velo(:,:,:),tao(:,:),C(:,:),freq(:,:)
real(4)::Tem,klx,kly,scal
character(10)::str1,str2,str3

!晶格热导率的单位是W m-1k-1


open(unit=103,file="mesh.dat")
open(unit=105,file="velo.dat")
open(unit=107,file="cv.dat")
open(unit=109,file="relax.dat")
open(unit=110,file="kappL.dat")
open(unit=106,file="tcond.in")


read(106,*)str1,Tem
read(106,*)
read(106,*)
read(106,*)str2,nq
read(106,*)str3,nm






allocate(velo(nq,nm,3))
allocate(C(nq,nm))
allocate(tao(nq,nm))
allocate(freq(nq,nm))
allocate(weight(nq))



klx=0
kly=0
do im=1,nm
   do iq=1,nq
      read(103,*)weight(iq),freq(iq,im)
      read(105,*)velo(iq,im,:)
      read(107,*)C(iq,im)
      read(109,*)tao(iq,im)
      klx=klx+velo(iq,im,1)**2*C(iq,im)*tao(iq,im)*weight(iq)
      kly=kly+velo(iq,im,2)**2*C(iq,im)*tao(iq,im)*weight(iq)
      !write(110,'(3e15.8)')klx,kly
   end do
      !write(110,'(3e15.8)')klx,kly
      !write(110,'(3e15.8)')
end do
scal=sum(weight(:))
write(110,*)scal
!write(110,'(3e15.8)')Tem,klx*1d4,kly*1d4
write(110,'(3e15.8)')Tem,klx*1d4/scal,kly*1d4/scal
end
