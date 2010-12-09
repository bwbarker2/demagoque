program testfft
  implicit none

  integer, parameter :: L=5,M=5 !size of array
  integer, parameter :: LM=L*M

  real :: farrayr(LM)=(/1,2,3,4,5 &
                    ,6,7,8,9,10 &
                    ,11,12,13,14,15 &
                    ,16,17,18,19,20 &
                    ,21,22,23,24,25/)
  real :: farrayi(LM)=0.0
  real :: arrayr(L,M),arrayi(L,M)
  complex :: array(L,M)
  integer :: ier,il
  integer,parameter :: lensav=50,lenwrk=50
  real :: wsave(lensav),work(lenwrk)
  real :: cvec(M)

  arrayr = reshape(farrayr,(/5,5/))
  arrayi = reshape(farrayi,(/5,5/))

  array=CMPLX(arrayr,arrayi)

  

  call FT(L,M,arrayr,arrayi)
  arrayr=TRANSPOSE(arrayr)
  arrayi=TRANSPOSE(arrayi)
  call FT(L,M,arrayr,arrayi)
  arrayr=TRANSPOSE(arrayr)
  arrayi=TRANSPOSE(arrayi)


!  call IFT(L,M,arrayr,arrayi)


  array=CMPLX(arrayr,arrayi)


  !columns:
!  allocate(cvec(M))
  !initialize work arrays
!  write(*,*)'starting CFF'
!  call CFFT1I(M,WSAVE,LENSAV,IER)
!  !transform columns
!  do il=1,L
!     cvec=array(il,:)
!     write(*,*)cvec
!     write(*,*)'calling CFFT1F'
!     call CFFT1F(5,1,cvec,5,WSAVE,LENSAV,WORK,LENWRK,IER)
!     array(il,:)=cvec
!     write(*,*)'iterating'
!  enddo

  write(*,*)array

end program testfft
