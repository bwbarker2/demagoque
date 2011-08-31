MODULE phys_cons
  ! physical constants
  USE prec_def
  IMPLICIT NONE

  complex*16,   parameter :: imagi=cmplx(0.d0,1.d0,8)  !sqrt(-1)

!  REAL (long) , PARAMETER :: pi=3.1415926535897932d0
  REAL (long) , PARAMETER :: pi=4d0*atan(1d0)
  real (long) , parameter :: invpi=1d0/pi  !1/pi
  real (long) , parameter :: invsqrt2pi=1d0/sqrt(2d0*pi)

  REAL (long) , PARAMETER :: rho0=0.16d0       !nuclear saturation density [fm^-3]
  REAL (long) , PARAMETER :: hbc=197.326963d0  !h-bar*c [MeV fm]
  REAL (long) , PARAMETER :: hbc2=hbc*hbc
  REAL (long) , PARAMETER :: mp=938.272013d0   !proton mass [MeV/c^2]
  REAL (long) , PARAMETER :: mn=939.565560d0   !neutron mass
  REAL (long) , PARAMETER :: m0=(mp+mn)*0.5d0
  REAL (Long) , PARAMETER :: a0=931.494028d0   !atomic mass unit
  REAL (long) , PARAMETER :: hm=hbc*hbc/(2.d0*m0)   !hbar^2/2m (useful for kinetic 
                                             ! time evolution)
  REAL (long) , PARAMETER :: deg=4.d0 ! degeneracy

END MODULE phys_cons


