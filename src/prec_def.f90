!> defines precision of numbers to use.
MODULE prec_def
 use ISO_FORTRAN_ENV
  IMPLICIT NONE

  !> KIND of numbers. REAL64 is the KIND number of a 64-bit real number
  !! in this architecture (provided by ISO_FORTRAN_ENV)
  INTEGER,PARAMETER :: long=REAL64

  integer,parameter :: stderr=102  !file unit of standard error
  
END MODULE prec_def
