program bmath_dZeroBrent_run
 use bexception
 use bmath
 use bmath_dZeroBrent_functions
 use iso_fortran_env
 use prec_def
 implicit none

 real (Long) :: answer
 real (Long) :: eps

! write(OUTPUT_UNIT,*)flin(2e0_Long)

 eps=1e-10_Long
 answer=bmath_dZeroBrent(-4e0_Long,4e0_Long/3e0_Long,flin,err=eps)
 if(abs(answer+3_Long)>eps) call throwException('answer not within epsilon',BEXCEPTION_FATAL)

 

! stop 1

end program bmath_dZeroBrent_run


