program deal_dos

implicit none

REAL(8), PARAMETER :: pi     = 3.141592653589793238462643383279502884197169399375105820974944d0
integer:: iskip0,iskip1,nskip1
integer:: natoms,nm,im
integer:: nq,iq
real(4),allocatable ::freq(:,:)
integer,allocatable ::weight(:)
character(10)::str1,str2,str3,str4,str5


open(unit=100,file="BTE.omega")                 !打开要处理的文件
open(unit=111,file="BTE.qpoints")
open(unit=103,file="mesh.dat")                 !处理之后的数据放在这里 

nq=3661
nm=36

allocate(freq(nq,nm))
allocate(weight(nq))

do iq=1,nq
read(100,*)freq(iq,:)
read(111,*)str1,str2,weight(iq),str3,str4,str5
write(*,'(36f15.8)')freq(iq,:)
end do


do im=1,nm
   do iq=1,nq
      write(103,'(I2,f15.8)')weight(iq),freq(iq,im)/(2*pi)
   end do
   write(103,*)
end do


end
