PROGRAM testint
  !! testint: tests my integrator
  implicit none

  integer , parameter :: Long = 8

  integer , parameter :: n = 20  ! (number-1)/2 of array elements
  real (Long) :: f(-n:n)         ! array
  real (Long) :: x, dx              ! x, delta-x
  real (Long) :: sum, err        ! integral, uncertainty in integral

  integer :: i

  dx = 0.1

  do i=-n,n
     x = i*dx
     f(i)=dexp(-x*x)
     write(*,*)x, f(i)
  enddo

  call dint_simp1(2*n+1, f, dx, sum, err)

  WRITE(*,*)'sum,err=',sum,err

END PROGRAM testint
