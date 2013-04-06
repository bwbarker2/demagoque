module fcns
 use prec_def
 use bexception

contains

 real (Long) function fcn_plus(ab)
  implicit none
  
  real (Long), dimension(:), intent(in) :: ab !< parameters a and b

  if(size(ab)/=2) call throwException('fcns: fcn_plus: needs 2 arguments',BEXCEPTION_FATAL)

  fcn_plus = ab(1) + ab(2)

 end function fcn_plus

end module fcns

program bmath_propagateError_1
 use bmath
 use fcns
 use phys_cons

 real(Long), dimension(2) :: params2 = (/10._Long,0.5_Long/)
 real(Long), dimension(2) :: errs2 = (/1._Long,0.01_Long/)

 real(Long) :: newerr,pdiff

 integer :: i, n, clock
 INTEGER, DIMENSION(:), ALLOCATABLE :: seed

 call phys_cons_init
 call random_seed(size=n)
  
ALLOCATE(seed(n))
          
            CALL SYSTEM_CLOCK(COUNT=clock)
          
            seed = clock + 37 * (/ (i - 1, i = 1, n) /)
            CALL RANDOM_SEED(PUT = seed)
          
            DEALLOCATE(seed)

! write(*,*)fcn_plus(params2)

 

 newerr = bmath_propagateError(fcn_plus,params2,errs2)
 pdiff = sqrt(errs2(1)**2 + errs2(2)**2)
 pdiff = (newerr-pdiff)/pdiff*100._Long

 

! write(*,*)'pdiff=',newerr,pdiff

end program bmath_propagateError_1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


