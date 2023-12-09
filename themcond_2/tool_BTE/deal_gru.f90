program deal_gru

implicit none

integer:: iskip0,iskip1,nskip1
integer:: natoms,nm,im
integer:: nq,iq
real(4),allocatable ::gru(:,:)
integer,allocatable ::weight(:)
character(10)::str1,str2,str3,str4


open(unit=101,file="BTE.gruneisen")                 !打开要处理的文件
open(unit=104,file="gru.dat")                 !处理之后的数据放在这里 

nq=3661
nm=36

allocate(gru(nq,nm))


do iq=1,nq
   read(101,*)gru(iq,:)
end do



do im=1,nm
   do iq=1,nq
      write(104,'(f15.8)')gru(iq,im)
   end do
   write(104,*)
end do

end
