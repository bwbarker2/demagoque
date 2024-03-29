module functions
 implicit none

contains

 pure real(Long) function fsin(xx)
  use prec_def

  real(Long), intent(in) :: xx

  fsin=sin(xx)
 end function fsin

end module functions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_bmath_LDiffRichardson
 use bmath
 use functions
 use phys_cons
 implicit none

 real(Long) :: diff, err

 call phys_cons_init

 diff = bmath_LDiffRichardson(fsin,acos(1._Long/3._Long),pi*0.1_Long,err)

end program test_bmath_LDiffRichardson
