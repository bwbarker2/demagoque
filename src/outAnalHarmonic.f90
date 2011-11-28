subroutine outAnalHarmonic
 !! outAnalHarmonic - outputs analytic time evolution of gaussian wavepacket
  ! in free space using same frequency as initialState_gaussianNuclear
 use input_parameters
 use mesh
 use time
 implicit none

  INTEGER :: ixa
  REAL (Long) :: ddre
  real (Long) :: calcHarmonicEv

  WRITE(45,*)'# time=',t,'fm/c'
  WRITE(45,*)'# x [fm], real, imaginary amplitudes'

  DO ixa=-Nxa2,Nxa2-1
     !add wave packets from surrounding cells, important for long times
     ddre=calcHarmonicEv(xa(ixa),t) &
          +calcHarmonicEv(2e0_Long*xLa+xa(ixa),t) &
          +calcHarmonicEv(2e0_Long*xLa-xa(ixa),t) &
          +calcHarmonicEv(4e0_Long*xLa+xa(ixa),t) &
          +calcHarmonicEv(4e0_Long*xLa-xa(ixa),t)
         
     WRITE(45,93) xa(ixa),ddre,0e0_Long
  ENDDO

  WRITE(45,*)
  WRITE(45,*)

93 format(3e16.8)

end subroutine outAnalHarmonic

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function calcHarmonicEv(xx,tt) result(self)
 !! computes the time evolution of the square modulus of a normalized
  ! Gaussian wave packet. For derivation, see written notes BWB 2010-12-15.
  ! Form of Gaussian:
  ! rho(x,t) = sqrt(1/(pi*sigma^2))*exp(-x^2/sigma^2)
  ! , where
  ! sigma(t)^2 = sigma0^2 (1 + (hbar*t/(msigma0^2)^2))
 use input_parameters
 use phys_cons
 implicit none

 real (Long) :: self

 real (Long), intent(in) :: xx, tt
 real (Long)             :: sigma2   !sigma squared
 real (Long)             :: sigma02  !sigma-naught squared

 sigma02=1e0_Long/whm

 sigma2=sigma02*(1e0_Long+(hbar*tt/(m0*sigma02))**2)

 self=sqrt(1e0_Long/(pi*sigma2))*exp(-xx**2/sigma2)

!do we need this check anymore?
! if(calcHarmonicEv<1.d-30)calcHarmonicEv=0.d0

end function calcHarmonicEv
