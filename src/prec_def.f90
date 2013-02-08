!> defines precision of numbers to use.
MODULE prec_def
 use ISO_FORTRAN_ENV
  IMPLICIT NONE

 !> Maximum positive value of a double precision number. From
 !! https://en.wikipedia.org/w/index.php?title=Double-precision_floating-point_format&oldid=489929288
 real(kind=8), parameter :: maxdouble=(1.d0+(1.d0-2.d0**(-52)))*2.d0**1023

  !> KIND of numbers. REAL64 is the KIND number of a 64-bit real number
  !! in this architecture (provided by ISO_FORTRAN_ENV)
!  INTEGER,PARAMETER :: long=REAL64
  integer, parameter :: long=selected_real_kind(15,50)

  !> Machine epsilon for KIND parameter given above. Use this if you are
  !! determining whether two floating-point numbers are equal. For example,
  !! if(abs(a-b)<epzero) write(*,*)'a and be are equal'
  real(Long) :: epzero = (2._Long)**(-52)

  !> double precision hexadecimal representation of a signalling-not-a-number. Great for suppressing warnings about a variable maybe not being initialized.
  real(kind=8) :: snan

  real (Long), parameter :: maxlong=maxdouble

  integer,parameter :: stderr=102  !file unit of standard error

contains

 subroutine init_prec_def
  implicit none

  snan = real(Z'FFFFFFFF',8)

 end subroutine init_prec_def
END MODULE prec_def
