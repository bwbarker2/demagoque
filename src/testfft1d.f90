program testfft1d
  !! testfft1d - tests fft with simple N=4 function. See derivation on paper notes BWB 2010-08-30
  use lib_fftw
  implicit none

  integer, parameter :: Long=8

  integer, parameter :: N=4 !size of array

  complex (Long) :: array(0:N-1)
  complex*16     :: arrayend(0:N-1)
  
  complex*16, parameter :: imagi=cmplx(0d0,1d0,8)
  real*8, parameter :: pi=3.1415926535897932d0

  integer :: i

  do i=0,N-1
   array(i)=cos(-pi/2d0+i*pi/4d0)
   write(*,*)i,array(i)
  enddo

  !shift from (-L to L) transform to (0,2L) transform
  write(*,*)'shift rho_x'
  do i=0,N-1
   array(i)=array(i)*exp(imagi*pi*i)
   write(*,*)array(i)
  enddo

  call ft_z2z_1d_naive(array,arrayend,N)

  !shift from (-L to L) transform to (0,2L) transform
  write(*,*)'shift rho_k'
  do i=0,N-1
   arrayend(i)=arrayend(i)*exp(-imagi*(-4d0+2d0*i)*(-pi/2d0))*2d0
   write(*,*)arrayend(i)
  enddo

!  write(*,*)arrayend

!  array=cmplx(arrayr,arrayi,8)

!  write(*,*)array

end program testfft1d

subroutine ft_z2z_1d_naive(arrayin,arrayout,num)
 implicit none

 integer, intent(in) :: num
 complex*16, intent(in) :: arrayin(0:num-1)
 complex*16, intent(out) :: arrayout(0:num-1)

 complex*16, parameter :: imagi=cmplx(0d0,1d0,8)
 real*8, parameter :: pi=3.1415926535897932d0

 integer :: i,j

 do i=0,num-1

  arrayout(i)=cmplx(0d0,0d0,8)

  do j=0,num-1

   arrayout(i)=arrayout(i)+arrayin(j)*exp(-imagi*2d0*pi*i*j/num)

  enddo

   arrayout(i)=arrayout(i)/sqrt(2d0*pi)
   write(*,*)i,arrayout(i)

 enddo 

end subroutine ft_z2z_1d_naive
