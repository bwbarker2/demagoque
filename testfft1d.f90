program testfft1d
  implicit none

  integer, parameter :: Long=8

  integer, parameter :: N=6 !size of array

  real (Long) :: arrayr(N)
  real (Long) :: arrayi(N)=0.0
  complex (Long) :: array(N)
  integer :: ier,in
  integer,parameter :: lensav=18, lenwrk=12
  real (Long) :: wsave(lensav),work(lenwrk)
  real (Long) :: cvec(N)

  arrayr=(/ 0.0, -0.86603,-0.86603, 0.0, 0.86603, 0.86603 /)

  array=CMPLX(arrayr,arrayi)

  ier=0
  call DCFFT1I(N,WSAVE,LENSAV,IER)
!  write(*,*)WSAVE
!  write(*,*)
  call DCFFT1F(N,1,array,N,WSAVE,LENSAV,WORK,LENWRK,IER)
  write(*,*)real(array)
  write(*,*)aimag(array)
  write(*,*)

!  call CFFT1B(N,1,array,N,WSAVE,LENSAV,WORK,LENWRK,IER)
!  write(*,*)real(array)
!  write(*,*)aimag(array)
!
!  do in=1,100
!     call CFFT1F(N,1,array,N,WSAVE,LENSAV,WORK,LENWRK,IER)
!     call CFFT1B(N,1,array,N,WSAVE,LENSAV,WORK,LENWRK,IER)
!  enddo
!
!  write(*,*)
!  write(*,*)real(array)
!  write(*,*)aimag(array)


!  call FT(L,M,arrayr,arrayi)
!  arrayr=TRANSPOSE(arrayr)
!  arrayi=TRANSPOSE(arrayi)
!  call FT(L,M,arrayr,arrayi)
!  arrayr=TRANSPOSE(arrayr)
!  arrayi=TRANSPOSE(arrayi)


!  call IFT(L,M,arrayr,arrayi)


!  array=CMPLX(arrayr,arrayi)


  !columns:
!  allocate(cvec(M))
  !initialize work arrays
!  write(*,*)'starting CFF'
!  !transform columns
!  do il=1,L
!     cvec=array(il,:)
!     write(*,*)cvec
!     write(*,*)'calling CFFT1F'
!     array(il,:)=cvec
!     write(*,*)'iterating'
!  enddo

!  write(*,*)array

end program testfft1d
