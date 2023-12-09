program deal_velo

implicit none

integer:: iskip0,iskip1,nskip1,iskip2,nskip2
integer:: natoms,nm,im
integer:: nq,iq
real(4),allocatable ::velo(:,:,:)
integer,allocatable ::weight(:)

character(10)::str1,str2,str3,str4,str5,str6
real(4)::vlo

!!!!!!!速度的单位是A.THZ

open(unit=102,file="BTE.v")                 !打开要处理的文件
open(unit=105,file="velo.dat")                 !处理之后的数据放在这里 

nq=3661
nm=36


allocate(velo(nq,nm,3))
!allocate(weight(nq))

do im=1,nm
   do iq=1,nq
      read(102,*) velo(iq,im,:)
      write(105,'(3f15.8)')velo(iq,im,1)*10,velo(iq,im,2)*10,velo(iq,im,3)*10
   end do
      write(105,'(3f15.8)')
end do


end
