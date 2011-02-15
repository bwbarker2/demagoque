! defines precision of types
MODULE prec_def
  IMPLICIT NONE
  
  INTEGER,PARAMETER :: long=8  !just use double precision for now
!  integer,parameter :: long=selected_real_kind(15,307)
  
END MODULE prec_def

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module skyrme_params
 implicit none

 !mean field interaction parameters
 ! Vautherin, Treiner, Veneroni PB191,6(1987)
 real*8, parameter :: t0  = -2164.1d0
 real*8, parameter :: t3  = 15079.0d0
 real*8, parameter :: sig = 0.25d0
 
end module skyrme_params

module fftw_constants
 use prec_def
 implicit none

 include '/usr/include/fftw3.f'

 integer*8 :: planf, planb !plan forward, plan reverse 2D full transforms
 integer*8 :: planwcos, planwsin !plan for cos,sin transform
 integer*8 :: planwx  ! plan for wigner to space transform

 real*8, allocatable :: arraycos(:), arraysin(:)
 complex*16, allocatable :: arraycnum(:)

! complex (Long), dimension(:,:), allocatable :: tempdenmat

end module


module cons_laws
  use prec_def
  implicit none

  real (Long) :: ekin   ! total energy calculated in k-space
  real (Long) :: ekerr  ! uncertainty of above
  real (Long) :: ek0    ! initial energy calc'd in k-space
  real (Long) :: ek0err ! error in initial

  real (Long), allocatable :: potx(:) ! potential on diagonal
  real (Long) :: epot   ! total potential energy in x-space
  real (Long) :: eperr  ! uncertainty of above
  real (Long) :: ep0    ! initial pot energy
  real (Long) :: ep0err ! error in initial

  real (Long) :: nnum   ! number of particles/(Nmax+1)

end module cons_laws

module formatting
 implicit none

 character(len=20), parameter :: fr5 = "(5E17.9)"

end module formatting

MODULE fft_vars
!! fft_vars - variables for FFT
  USE prec_def
  IMPLICIT NONE

  ! lensav(1) is for L, lensav(2) is for M
  INTEGER :: lensav(2), lenwrk(2)


  ! for the following, the first index: 1 - L, 2 - M

  REAL (Long), ALLOCATABLE :: work(:,:)

  !arrays for cosine transforms
  REAL (Long), ALLOCATABLE :: wsavec(:,:)

  !arrays for sine transforms
  REAL (Long), ALLOCATABLE :: wsaves(:,:)

END MODULE fft_vars


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! physical constants
MODULE phys_cons
  USE prec_def
  IMPLICIT NONE

  complex*16,   parameter :: imagi=cmplx(0.d0,1.d0,8)  !sqrt(-1)
  REAL (long) , PARAMETER :: pi=3.1415926535897932d0
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

!  integer :: usePotHOext ! use external HO potential in evolution?
!  integer :: usePotHOmf  ! use meanfield HO potential in evolution?
!  integer :: useAdiabaticSwitching 
!  integer :: useSkyrmeMF

END MODULE phys_cons


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! harmonic oscillator parameters
MODULE osc_pars
  USE prec_def
  IMPLICIT NONE

  REAL (long) :: w
  REAL (long) :: whm
  INTEGER     :: Nmax   ! maximum oscillator shell
  
END MODULE osc_pars


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! time parameters
MODULE time
  USE prec_def

  INTEGER     :: it     ! current iteration of time
  REAL*8      :: t      ! current time during time evolution [fm/c]
  REAL*8      :: delt   ! timestep [fm/c]
  integer     :: Nevt   ! number of timesteps for evolution
  INTEGER     :: Nt     ! number of timesteps in current mode (adiabatic or
                        !  time)

  logical     :: useImEvol   ! use imaginary evolution?
  integer     :: Nimev  ! number of timesteps for imaginary evolution

  real*8      :: ea     ! energy per particle [MeV]

  integer     :: iadib  ! 1 = run adiabatic, 0 = read adiabatic from file
  integer     :: Nad    !adiabatic switching parameters
  real*8      :: tad,wtad  
END MODULE time


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! output variables
MODULE out
  USE prec_def

  INTEGER :: ntime   ! write data every ntime timesteps
END MODULE out
