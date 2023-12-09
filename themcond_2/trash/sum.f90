program san0

implicit none




real(4):: x(5)=(/1, 2, 3, 4, 5/)
real(4):: y(5)=(/0.1, 0.2, 0.3, 0.4, 0.5/)
real(4),allocatable:: z(:),k(:)
real(4)::w,p
real(4) :: j
integer::i


allocate(z(10))




z(1)=sum(x(:)*y(:))

write(*,*)z(1)


w=1d-30

p=123456789

write(*,*)w
write(*,'(2e12.5)')p


!j=0
!do w=1,2
!   do i=1,5
!      j=i**2
!      write(*,*)j
!   end do
!write(*,*)j

!#end do

j=5*2**2*3
write(*,*)j

write(*,*)sqrt(j)

end
