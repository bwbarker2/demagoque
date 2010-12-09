!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     FOURIER TRANSFORM THE DENSITY MATRIX WITH THE HELP OF  c
!c     SOME FFT NAG ROUTINES. THE DIFFICULTY COMES FROM THE   c
!c     FACT THAT RIGHT AND LEFT COMPLEX EXPONENTIALS HAVE     c
!c     DIFFERENT SIGNS                                        c
!c     MAXIMUM N+1=1000!!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
SUBROUTINE FT(Nx,Ny,xre,xim)
  use prec_def
  implicit none

  ! matrices in the routine
  integer :: Nx,Ny
  real (Long) :: xre(Nx,Ny),xim(Nx,Ny)

  ! parameters of NAG routine for FFT
  integer,parameter :: MMAX=1000, NMAX=1000
  integer :: i,ifail,j,m,n
  real (Long) :: trig(2*NMAX),work(2*MMAX*NMAX)
  external c06frf,c06gcf

  ! parameters of NAG routine for matrix transpose
  integer :: mmat,nmat, mnmat
  integer,parameter :: LMOVE=(MMAX+NMAX)/2
  integer :: move(LMOVE)
  external f01crf

  mmat=Nx
  nmat=Ny
  mnmat=mmat*nmat

  ifail=0

  !fourier transform the columns
  call C06FRF(MMAT,NMAT,xre,xim,'I',TRIG,WORK,IFAIL)

  !transpose the arrays
  call F01CRF(xre,MMAT,NMAT,MNMAT,MOVE,LMOVE,IFAIL)
  call F01CRF(xim,MMAT,NMAT,MNMAT,MOVE,LMOVE,IFAIL)
  
  ! THE MINUS SIGN IS ACHIEVED BY DOUBLE COMPLEX CONJUGATION
  call C06GCF(xim,MNMAT,IFAIL)
  ! FOURIER TRANSFORM THE ROWS
  call C06FRF(NMAT,MMAT,xre,xim,'I',TRIG,WORK,IFAIL)
  ! THE MINUS SIGN IS ACHIEVED BY DOUBLE COMPLEX CONJUGATION
  call C06GCF(xim,MNMAT,IFAIL)
  
  ! TRANSPOSE BACK
  call F01CRF(xre,NMAT,MMAT,MNMAT,MOVE,LMOVE,IFAIL)
  call F01CRF(xim,NMAT,MMAT,MNMAT,MOVE,LMOVE,IFAIL)
  
end SUBROUTINE FT


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     INVERSE FOURIER TRANSFORM THE DENSITY MATRIX WITH THE  c
!c     HELP OF SOME FFT NAG ROUTINES. THE DIFFICULTY COMES    c
!c     FROM THE FACT THAT RIGHT AND LEFT COMPLEX EXPONENTIALS c
!c     HAVE DIFFERENT SIGNS                                   c
!c     MAXIMUM N+1=1000!!                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
SUBROUTINE IFT(Nx,Ny,xre,xim)
  use prec_def
  implicit none
  
  ! matrices in the routine
  integer :: Nx,Ny
  real (Long) :: xre(Nx,Ny),xim(Nx,Ny)

  ! parameters of NAG routine for FFT
  integer,parameter :: MMAX=1000, NMAX=1000
  integer :: i,ifail,j,m,n
  real (Long) :: trig(2*NMAX),work(2*MMAX*NMAX)
  external c06frf,c06gcf

  ! parameters of NAG routine for matrix transpose
  integer :: mmat,nmat, mnmat
  integer,parameter :: LMOVE=(MMAX+NMAX)/2
  integer :: move(LMOVE)
  external f01crf

  mmat=Nx
  nmat=Ny
  mnmat=mmat*nmat

  ifail=0

  ! THE MINUS SIGN IS ACHIEVED BY DOUBLE COMPLEX CONJUGATION
  call C06GCF(xim,MNMAT,IFAIL)
  ! FOURIER TRANSFORM THE COLUMNS
  call C06FRF(MMAT,NMAT,xre,xim,'I',TRIG,WORK,IFAIL)
  ! THE MINUS SIGN IS ACHIEVED BY DOUBLE COMPLEX CONJUGATION
  call C06GCF(xim,MNMAT,IFAIL)

  !transpose the arrays
  call F01CRF(xre,MMAT,NMAT,MNMAT,MOVE,LMOVE,IFAIL)
  call F01CRF(xim,MMAT,NMAT,MNMAT,MOVE,LMOVE,IFAIL)

  ! FOURIER TRANSFORM THE ROWS
  call C06FRF(NMAT,MMAT,xre,xim,'I',TRIG,WORK,IFAIL)

  ! TRANSPOSE BACK
  call F01CRF(xre,NMAT,MMAT,MNMAT,MOVE,LMOVE,IFAIL)
  call F01CRF(xim,NMAT,MMAT,MNMAT,MOVE,LMOVE,IFAIL)

end SUBROUTINE IFT
