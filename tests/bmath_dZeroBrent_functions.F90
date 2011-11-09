module bmath_dZeroBrent_functions
 use prec_def
 implicit none

contains

 real(Long) function flin(x)
  implicit none

  real (Long), intent(in) :: x

  flin=(x+3_Long)*(x-1_Long)**2
!  flin=10e0_Long/(1_Long+tan(x)**2)-x**2

 end function flin

end module bmath_dZeroBrent_functions


