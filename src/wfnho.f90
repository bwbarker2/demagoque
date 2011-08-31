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

  xnorm=1/sqrt((2.d0**n)*dble(nfac))*(whm/pi)**0.25d0
  xx=x*sqrt(whm)
  
  wfnho=xnorm*dexp(-xx**2/2.d0)*Hn(xx,n)

!  write(*,*)wfnho
  
end FUNCTION wfnho


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c ... HERMITE POLYNOMIAL OF ORDER N
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
real*8 FUNCTION Hn(x,n)
  use prec_def
  implicit none

  real (long), intent(in) :: x
  integer, intent(in) :: n

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
