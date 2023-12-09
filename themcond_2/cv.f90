program cv

implicit none

real(kind=8),parameter :: kb=1.380648813d-23 ! J/K
real(kind=8),parameter :: hbar=1.05457172647d-22 ! J/THz
integer:: nm,im                      !模式数,im遍历模式
integer:: nq,iq                      !q点个数,iq遍历q点
real(4),allocatable ::freq(:,:)      !声子频率
integer,allocatable ::weight(:)      !q点权重
real(4):: x,Tem,volum  !,Cvei        !hbar*freq/kb*T,温度T,体积
real(8),allocatable ::C(:,:)         !热容
character(10)::str1,str2,str3,str4   !不需要的字符串
real(4)::cj

!!!!!!!!!最终热容的单位是J K-1 m-3


open(unit=103,file="mesh.dat")
open(unit=106,file="tcond.in")
open(unit=107,file="cv.dat")


read(106,*)str1,Tem
read(106,*)str2,volum
read(106,*)
read(106,*)str3,nq
read(106,*)str4,nm



allocate(freq(nq,nm))
allocate(weight(nq))
allocate(C(nq,nm))


volum=volum*1d-30

do im=1,nm
   do iq=1,nq
      read(103,*)weight(iq),freq(iq,im)
      if (freq(iq,im)>0)then
         x=hbar*freq(iq,im)/(kb*Tem)
         !Cvei=kb*x**2*exp(x)/(volum*(exp(x)-1)**2)
         C(iq,im)=kb*x**2*exp(x)/(volum*(exp(x)-1)**2) 
         !write(107,'(f15.8)')Cvei
         write(107,'(f15.8)')C(iq,im)
      else
         write(107,'(f15.8)')0
      end if
   end do
   read(103,*)
   write(107,*)
end do

cj=0
do im=1,nm
do iq=1,nq
cj=cj+C(iq,im)*weight(iq)
end do
end do
write(*,*)cj
end
