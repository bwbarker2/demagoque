!> \brief Stores mathematical and physical constants, both in SI units and user-defined units.
!!
!! \todo give citations for all units
!! \todo code them in SI, provide conversion routines, see
!!       class_PhysicalConstants for current work
MODULE phys_cons
  USE prec_def
  IMPLICIT NONE

  !math constants (unitless)
  complex*16,   parameter :: imagi=cmplx(0.d0,1.d0,8)  !sqrt(-1)
  REAL (long) , PARAMETER :: pi=4d0*atan(1d0)
  real (long) , parameter :: invpi=1d0/pi  !1/pi
  real (long) , parameter :: invsqrt2pi=1d0/sqrt(2d0*pi)

  !physical constants in SI units, with "current units" as well

  !> relative atomic mass of 87 Rubidium, in u, from Audi 2003
  real (Long), parameter :: &
       MASS_REL_RUBIDIUM_87 = 86.909180527_Long, &
       MASS_REL_RUBIDIUM_87_D = 0.000000013_Long

  !> Mass of the neutron, in MeV/c<sup>^2</sup>, from 2010 CODATA
  real (Long), parameter :: &
       MASS_NEUTRON_MEV_C2 = 939.565379_Long, &
       MASS_NEUTRON_MEV_C2_D = 0.000021_Long  !< uncertainty

  !> Mass of the proton, in MeV/c<sup>^2</sup>, from 2010 CODATA
  real (Long), parameter :: MASS_PROTON_MEV_C2 = 938.272046_Long
  !> uncertainty
  real (Long), parameter :: MASS_PROTON_MEV_C2_D = 0.000021_Long

  !> Planck's constant divided by 2π, in MeV fm/c, from 2010 CODATA
  real (Long), parameter :: &
       PLANCKS_CONSTANT_HBAR_MEV_FM_C = 197.3269718_Long, &
       PLANCKS_CONSTANT_HBAR_MEV_FM_C_D = 0.0000044_Long  !< uncertainty

  !> speed of light in vacuum, c, in m s<sup>-1</sup>, from 2010 CODATA
  !! recommended values
  real (Long), parameter :: SI_SPEED_OF_LIGHT   = 299792458d0
  real (Long), parameter :: SI_SPEED_OF_LIGHT_D = 0d0 !< uncertainty
  real (Long)            :: SPEED_OF_LIGHT            !< in current units

  !> Planck's constant divided by 2π, ℏ, in J s, from 2010 CODATA
  real (Long), parameter :: SI_PLANCKS_CONSTANT_HBAR   = 1.054571726e-34_Long
  !> uncertainty
  real (Long), parameter :: SI_PLANCKS_CONSTANT_HBAR_D = 0.000000047e-34_Long
  !> in current units
  real (Long)            :: PLANCKS_CONSTANT_HBAR

  !> atomic mass constant, u, in kg, from 2010 CODATA
  real (Long), parameter :: &
       SI_ATOMIC_MASS_CONSTANT = 1.660538921e-27_Long, &
       SI_ATOMIC_MASS_CONSTANT_D = 0.000000073e-27_Long  !< uncertainty

  !> mass of the neutron, m<sub>n</sub>, in kg, from 2010 CODATA
  real (Long), parameter :: SI_MASS_NEUTRON   = 1.674927351e-27_Long
  !> uncertainty
  real (Long), parameter :: SI_MASS_NEUTRON_D = 0.000000074e-27_Long
  !> in current units
  real (Long)            :: MASS_NEUTRON

  !> mass of the proton, m<sub>p</sub>, in kg.
  real (Long), parameter :: SI_MASS_PROTON   = 1.672621777e-27_Long
  !> uncertainty
  real (Long), parameter :: SI_MASS_PROTON_D = 0.000000074e-27_Long
  !> in current units
  real (Long)            :: MASS_PROTON

  REAL (long) , PARAMETER :: rho0=0.16d0       !nuclear saturation density [fm^-3]
  REAL (long) :: hbar  !< Planck's constant in current unit system
  REAL (long) :: m0   !< mass of particle being evolved in time
  REAL (long) , PARAMETER :: deg=4.d0 ! degeneracy

contains

 !> \brief Initializes physical constants, using nuclear system:
 !!
 !! -# length - femtometer - fm
 !! -# mass   - mega-electron volt per speed of light squared - MeV/c^2
 !! -# time   - femtometer / speed of light - fm/c
 !!
 subroutine phys_cons_initializeNuclear
  implicit none

  hbar = PLANCKS_CONSTANT_HBAR_MEV_FM_C  !MeV fm
  m0 = 0.5_Long*(MASS_NEUTRON_MEV_C2+MASS_PROTON_MEV_C2)
  

 end subroutine phys_cons_initializeNuclear

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11

 !> \brief Initializes physical constants for BEC evolution (87-Rb)
 !!
 !! -# length = micrometer - um
 !! -# mass   = atomic mass unit - u
 !! -# time   = millisecond - ms
 subroutine phys_cons_initializeBEC
  implicit none

  hbar = SI_PLANCKS_CONSTANT_HBAR &
         *1_Long    / SI_ATOMIC_MASS_CONSTANT  &  ! u / kg
         * 1_Long / MASS_REL_RUBIDIUM_87       &  ! 1 mass of 87Rb / u
         *(1e6_Long / 1_Long)**2               &  ! um / m twice
         *1_Long    / 1e3_Long                    ! s / ms
 
  write(*,*)'hbar=',hbar
 
!  m0 = MASS_REL_RUBIDIUM_87
  m0 = 1_Long

 end subroutine phys_cons_initializeBEC
 

END MODULE phys_cons


