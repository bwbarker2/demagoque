subroutine outAnalHarmonic
 !! outAnalHarmonic - outputs analytic time evolution of gaussian wavepacket
  ! in free space using same frequency as initialState_gaussianNuclear
 use input_parameters
 use mesh
 use time
 implicit none

  INTEGER :: ixa
  REAL*8 :: ddre
  real*8 :: calcHarmonicEv

  WRITE(45,*)'# time=',t,'fm/c'
  WRITE(45,*)'# x [fm], real, imaginary amplitudes'

  DO ixa=-Nxa2,Nxa2-1
     !add wave packets from surrounding cells, important for long times
     ddre=calcHarmonicEv(xa(ixa),t) &
          +calcHarmonicEv(2d0*xLa+xa(ixa),t) &
          +calcHarmonicEv(2d0*xLa-xa(ixa),t) &
          +calcHarmonicEv(4d0*xLa+xa(ixa),t) &
          +calcHarmonicEv(4d0*xLa-xa(ixa),t)
         
     WRITE(45,93) xa(ixa),ddre,0.d0
  ENDDO

  WRITE(45,*)
  WRITE(45,*)

93 format(3e16.8)

end subroutine outAnalHarmonic

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real*8 function calcHarmonicEv(xx,tt)
 !! computes the time evolution of the square modulus of a normalized
  ! Gaussian wave packet. For derivation, see written notes BWB 2010-12-15.
  ! Form of Gaussian:
  ! rho(x,t) = sqrt(1/(pi*sigma^2))*exp(-x^2/sigma^2)
  ! , where
  ! sigma(t)^2 = sigma0^2 (1 + (hbar*t/(msigma0^2)^2))
 use input_parameters
 use phys_cons
 implicit none

 real*8, intent(in) :: xx, tt
 real*8 :: sigma2   !sigma squared
 real*8 :: sigma02  !sigma-naught squared

 sigma02=1d0/whm

 sigma2=sigma02*(1d0+(hbar*tt/(m0*sigma02))**2)

 calcHarmonicEv=sqrt(1/(pi*sigma2))*exp(-xx**2/sigma2)

!do we need this check anymore?
! if(calcHarmonicEv<1.d-30)calcHarmonicEv=0.d0

end function calcHarmonicEv
