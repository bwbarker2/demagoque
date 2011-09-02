!> \brief   Produces amplitudes of harmonic oscillator wave functions
!! \details Wavefunction is defined as sum of N=0 to N=n hamonic oscillator
!!          eigenvectors.
!!
!! \author  Arnau Rios
!!
!! \returns amplitude of sum of harmonic oscillator wavefunctions at position x
!!
real (Long) FUNCTION wfnho(x,n,whm)
  use phys_cons
  use prec_def
  implicit none

  real (Long), intent(in) :: x   !< position at which to calculate wavefunction
  integer,     intent(in) :: n   !< order of maximum harmonic oscillator used
  real (Long), intent(in) :: whm !< angular wavenumber squared

  integer :: nfac
  integer :: ifac
  real (long) :: xnorm &  !< normalization of wavefunction
                 ,xx,Hn

  nfac=1

  if(n .gt. 1) then
     do ifac=2,n
        nfac=nfac*ifac
     enddo
  endif
  
!  write(*,*)nfac

  xnorm=1/sqrt((2.d0**n)*dble(nfac))*(whm/pi)**0.25d0
  xx=x*sqrt(whm)
  
  wfnho=xnorm*dexp(-xx**2/2.d0)*Hn(xx,n)

!  write(*,*)wfnho
  
end FUNCTION wfnho


!> \brief Computes Hermite polynomial of order n at position x
!!
!! \author Arnau Rios
!!
real (Long) FUNCTION Hn(x,n)
  use prec_def
  implicit none

  real (long), intent(in) :: x  !< position at which to compute polynomial
  integer,     intent(in) :: n  !< order of polynomial

  real (long) :: hm,hm1,hm2
  integer :: m
  
  if(n.eq.0) then
     Hn=1d0
  else if(n.eq.1) then
     Hn=2.d0*x
  else if(n.gt.1) then
     hm2=1d0
     hm1=2d0*x
     do m=2,n
        hm=2d0*x*hm1 - 2d0*(m-1)*hm2
        hm2=hm1
        hm1=hm
     enddo
     Hn=hm
   else
     write(*,*)
     write(*,*)'*****'
     write(*,*)'Hn: cannot have negative order of Hermite polynomial, n=,',n
     write(*,*)'Hn: setting n=0 -> Hn=1'
     write(*,*)'*****'
     write(*,*)
     Hn=1d0
  endif
  
!  write(*,*)x,Hn
  
end FUNCTION Hn
