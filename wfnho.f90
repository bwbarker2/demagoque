!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c ... HARMONIC OSCILLATOR EIGENVECTORS
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
real*8 FUNCTION wfnho(x,n,whm)
  use phys_cons
  use prec_def
  implicit none

  integer :: nfac
  integer :: n,ifac
  real (long) :: x,whm,xnorm,xx,Hn

  nfac=1

  if(n .gt. 1) then
     do ifac=2,n
        nfac=nfac*ifac
     enddo
  endif
  
!  write(*,*)nfac

  xnorm=1/dsqrt((2.d0**n)*dble(nfac))*(whm/pi)**0.25d0
  xx=x*dsqrt(whm)
  
  wfnho=xnorm*dexp(-xx**2/2.d0)*Hn(xx,n)

!  write(*,*)wfnho
  
end FUNCTION wfnho


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c ... HERMITE POLYNOMIAL OF ORDER N
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
real*8 FUNCTION Hn(x,n)
  use prec_def
  implicit none

  real (long) :: x,hm,hm1,hm2
  integer :: n,m
  
  if(n.eq.0) then
     Hn=1
  else if(n.eq.1) then
     Hn=2.d0*x
  else if(n.gt.1) then
     hm2=1
     hm1=2.0*x
     do m=2,n
        hm=2.0*x*hm1 - 2.0*(dble(m-1))*hm2
        hm2=hm1
        hm1=hm
     enddo
     Hn=hm
  endif
  
!  write(*,*)x,Hn
  
end FUNCTION Hn
