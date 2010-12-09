program interp_test
 use prec_def
 implicit none

 ! test find_points
 integer, parameter :: nx=10
 integer, parameter :: ny=10
 real (Long) :: xx(nx) = (/1.0,3.0,5.0,7.0,9.0,11.0,13.0,15.0,17.0,19.0/)
 real (Long) :: zz(nx) = (/2.0,6.0,10.0,14.0,18.0,22.0,26.0,30.0,34.0,38.0/)
 integer :: x1,x2
 real (Long) :: z

 call LIN_INT_1D(xx,zz,nx,23.0_Long,z)
! call find_points(xx,10,16.0_Long,x1,x2)

 write(*,*)z

end program interp_test

