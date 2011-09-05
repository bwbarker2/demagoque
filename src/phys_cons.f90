!> \brief Stores mathematical and physical constants, both in SI units and user-defined units.
!!
!! \todo give citations for all units
!! \todo code them in SI, provide conversion routines
MODULE phys_cons
  ! physical constants
  USE prec_def
  IMPLICIT NONE

  !math constants (unitless)
  complex*16,   parameter :: imagi=cmplx(0.d0,1.d0,8)  !sqrt(-1)
  REAL (long) , PARAMETER :: pi=4d0*atan(1d0)
  real (long) , parameter :: invpi=1d0/pi  !1/pi
  real (long) , parameter :: invsqrt2pi=1d0/sqrt(2d0*pi)

  !physical constants in SI units, with "current units" as well

  ! conversion factors from base units. Multiply SI unit by this to get
  ! current unit
  real (Long), private :: convertLength
  real (Long), private :: convertMass
  real (Long), private :: convertTime
  real (Long), private :: convertElectricCurrent
  real (Long), private :: convertThermodynamicTemperature
  real (Long), private :: convertAmountOfSubstance
  real (Long), private :: convertLuminousIntensity

  ! names, symbols of units in current system of units
  !> \todo change to using a string class that has variable length
  character(10) :: unitLengthName
  character(10) :: unitMassName
  character(10) :: unitTimeName
  character(10) :: unitElectricCurrentName

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
  REAL (long) , PARAMETER :: hbc=197.326963d0  !h-bar*c [MeV fm]
  REAL (long) , PARAMETER :: hbc2=hbc*hbc
  REAL (long) , PARAMETER :: mp=938.272013d0   !proton mass [MeV/c^2]
  REAL (long) , PARAMETER :: mn=939.565560d0   !neutron mass
  REAL (long) , PARAMETER :: m0=(mp+mn)*0.5d0
  REAL (Long) , PARAMETER :: a0=931.494028d0   !atomic mass unit
  REAL (long) , PARAMETER :: hm=hbc*hbc/(2.d0*m0)   !hbar^2/2m (useful for kinetic 
                                             ! time evolution)
  REAL (long) , PARAMETER :: deg=4.d0 ! degeneracy

contains

 !> \brief Initializes module, assumes current system of units is SI.
 subroutine physCons_initialize()
  implicit none


 end subroutine physCons_initialize

END MODULE phys_cons


