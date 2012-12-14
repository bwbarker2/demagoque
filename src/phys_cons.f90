!> \brief Stores mathematical and physical constants, both in SI units and user-defined units.
!!
!! \todo give citations for all units
!! \todo code them in SI, provide conversion routines, see
!!       class_PhysicalConstants for current work
MODULE phys_cons
  use iso_varying_string
  USE prec_def
  IMPLICIT NONE

  logical :: phys_cons_isInitialized

  !math constants (unitless)
  complex (Long) :: czero
  complex (Long) :: imagi
  REAL (long) :: pi
  real (Long) :: inv2pi
  real (long) :: invpi
  real (long) :: invsqrt2pi

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

  !> Planck's constant divided by 2*pi, in MeV fm/c, from 2010 CODATA
  real (Long), parameter :: &
       PLANCKS_CONSTANT_HBAR_MEV_FM_C = 197.3269718_Long, &
       PLANCKS_CONSTANT_HBAR_MEV_FM_C_D = 0.0000044_Long  !< uncertainty

  !> elementary charge, e, in C= A s, from 2010 CODATA
  real (Long), parameter :: SI_ELEMENTARY_CHARGE = 1.602176565e-19_Long &
                         ,SI_ELEMENTARY_CHARGE_D = 0.000000035e-19_Long

  !> speed of light in vacuum, c, in m s<sup>-1</sup>, from 2010 CODATA
  !! recommended values
  real (Long), parameter :: SI_SPEED_OF_LIGHT   = 299792458e0_Long
  real (Long), parameter :: SI_SPEED_OF_LIGHT_D = 0e0_Long !< uncertainty
  real (Long)            :: SPEED_OF_LIGHT            !< in current units

  !> Planck's constant divided by 2 pi, hbar, in J s, from 2010 CODATA
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

  REAL (long) , PARAMETER :: rho0=0.16_Long       !nuclear saturation density [fm^-3]
  REAL (long) :: hbar  !< Planck's constant in current unit system
  REAL (long) :: m0   !< mass of particle being evolved in time
  REAL (long) , PARAMETER :: deg=4e0_Long ! degeneracy

  real (Long) :: kilogram !< magnitude of kilogram in current unit system
  real (Long) :: joule !< magnitude of joule in current unit system
  real (Long) :: meter !< magnitude of meter in current unit system
  real (Long) :: second

  !> Abbreviation of unit of time in current unit system
  type(varying_string) :: phys_cons_unit_abbrev_time

contains

 subroutine phys_cons_init
  implicit none

  !math constants (unitless)
  czero=cmplx(0e0_Long,0e0_Long,Long)
  imagi=cmplx(0e0_Long,1e0_Long,Long)  !sqrt(-1)
  pi=4e0_Long*atan(1e0_Long)
  inv2pi=1e0_Long/(2e0_Long*pi)
  invpi=1e0_Long/pi  !1/pi
  invsqrt2pi=sqrt(inv2pi)

 end subroutine phys_cons_init


 !> \brief Initializes physical constants, using nuclear system:
 !!
 !! -# length - femtometer - fm
 !! -# mass   - mega-electron volt per speed of light squared - MeV/c^2
 !! -# time   - femtometer / speed of light - fm/c
 !!
 subroutine phys_cons_initializeNuclear
  implicit none

  call phys_cons_init
!  fm=1e0_Long
!  fm_c=1e0_Long
!  Mev=1e0_Long
!  Mev_c2=1e0_Long
!
!  amp=1e0_Long
!  meter=1e15_Long*fm
!  second=SI_SPEED_OF_LIGHT*meter
!  joule=
!  kilogram=1e-6_Long*amp*second**3*SI_SPEED_OF_LIGHT**2 &
!           /(SI_meter**2)
!
  hbar = PLANCKS_CONSTANT_HBAR_MEV_FM_C  !MeV fm/c
  m0 = 0.5_Long*(MASS_NEUTRON_MEV_C2+MASS_PROTON_MEV_C2)

  phys_cons_unit_abbrev_time='fm/c'

  phys_cons_isInitialized=.true.

 end subroutine phys_cons_initializeNuclear

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11

 !> \brief Initializes physical constants for BEC evolution (87-Rb)
 !!
 !! -# length = 10 micrometer - 10 um
 !! -# mass   = atomic mass unit - u
 !! -# time   = millisecond - ms
 !!
 !! -# energy \f$\approx\f$ 0.104 millielectronvolt - meV

 subroutine phys_cons_initializeBEC
  implicit none

  call phys_cons_init

  kilogram = 1e0_Long / SI_ATOMIC_MASS_CONSTANT
  meter = 1e5_Long
  second = 1e3_Long

  joule=kilogram*(meter/second)**2

  hbar = SI_PLANCKS_CONSTANT_HBAR*joule*second

!  write(*,*)'hbar=',hbar
 
  m0 = MASS_REL_RUBIDIUM_87

!  m0 = 1

  phys_cons_unit_abbrev_time='ms'

  phys_cons_isInitialized=.true.

 end subroutine phys_cons_initializeBEC
 

END MODULE phys_cons


