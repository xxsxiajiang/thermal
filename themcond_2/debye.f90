program debye

implicit none

REAL(8), PARAMETER :: pi     =3.141592653589793238462643383279502884197169399375105820974944d0
real(kind=8),parameter :: kb=1.380648813d-23 ! J/K
real(kind=8),parameter :: hbar=1.05457172647d-22 ! J/THz
integer :: nw,iw             !频率个数,第几个频率
integer :: nm,im             !模式数,第几个模式
integer :: nq,iq             !q点数,第几个q 
real(4) :: minw,maxw,dew     !最小频率,最大频率,频率间隔
real(4) :: smear,dos      !q点权重
real(4) :: factor1,factor2,debye2   
real(4),allocatable :: freq(:,:),w(:)
integer,allocatable ::weight(:)
character(10)::str1,str2,str3,str4,str5,str6,str7

!!!!!!!debye温度平方的单位是K2


open(unit=103,file="mesh.dat")
open(unit=108,file="debye.dat")
open(unit=106,file="tcond.in")


read(106,*)
read(106,*)
read(106,*)
read(106,*)str1,nq
read(106,*)str2,nm
read(106,*)
read(106,*)str3,nw
read(106,*)str4,minw
read(106,*)str5,maxw
read(106,*)str6,dew
read(106,*)
read(106,*)str7,smear


!!!!!!!!设置这些参数
allocate(w(nw))
do iw=1,nw                  !iw表示第iw个w
   w(iw)=minw+dew*(iw-1)         !第iw个频率的频率大小
end do
!!!!!!!!!!每一个模式求这些能量的态密度



allocate(freq(nq,nm))
allocate(weight(nq))


do im=1,nm
   do iq=1,nq
      read(103,*)weight(iq),freq(iq,im)
   end do
   read(103,*)
end do


do im=1,nm                    
   factor1=0                         
   factor2=0                                           
   do iw=0,nw
      dos=0
      do iq=1,nq     
         dos=dos+weight(iq)*exp(-(freq(iq,im)-w(iw))**2/(2*smear**2))/(smear*sqrt(2*pi)) 
      end do   
      !dos=sum(weight(1:nq)*exp(-(freq(1:nq,im)-w(iw))**2/(2*smear**2)))                       
      factor1= factor1+w(iw)**2*dos
      factor2= factor2+dos
      debye2=(5*hbar**2/(3*kb**2))*(factor1/factor2)
      !write(108,'(I2,f15.8)')im,debye2
   end do
   !write(108,'(I2,f15.8)')
   !debye2=(5*hbar**2/(3*kb**2))*(factor1/factor2)                                           
   write(108,'(I5,f15.8)')im,debye2
   !write(108,'(I2,f15.8)')
end do



end


