subroutine LIN_INT_1D(xx,zz,nx,x,z)
  !! LIN_INT_1D - linear interpolator for 1D array of values
  use prec_def
  implicit none

  integer, intent(in)      :: nx        ! number of points in x direction
  real (Long), intent(in)  :: xx(nx)    ! coordinates of interpolation points
  real (Long), intent(in)  :: zz(nx)    ! values of interpolation points
  real (Long), intent(in)  :: x         ! coordinates to interpolate to
  real (Long), intent(out) :: z         ! end value

  integer :: x1,x2                   ! endpoints for interpolation
!  real (Long) :: sw,nw,se,ne  ! distance from point to each endpoint
                              ! (think compass directions, sw=southwest, etc)
!  real (Long) :: xx1,xx2,yy1,yy2     ! values of endpoints

  call find_points(xx,nx,x,x1,x2)

  z=zz(x1)*(xx(x2)-x)+zz(x2)*(x-xx(x1))
  z=z/(xx(x2)-xx(x1))

end subroutine LIN_INT_1D
  


subroutine LIN_INT_2D(xx,yy,zz,nx,ny,x,y,z)
  !! LIN_INT_2D - linear interpolater for 2D array of values
  use prec_def
  implicit none

  integer, intent(in) :: nx,ny             ! number of points in x,y direction
  real (Long), intent(in) :: xx(nx),yy(ny) ! coordinates of interpolation points
  real (Long), intent(in) :: zz(nx,ny)     ! values of interpolation points
  real (Long), intent(in) :: x,y           ! coordinates to interpolate to
  real (Long), intent(out) :: z            ! end value

  integer :: x1,x2,y1,y2                   ! endpoints for interpolation
  real (Long) :: sw,nw,se,ne  ! distance from point to each endpoint
                              ! (think compass directions, sw=southwest, etc)
  real (Long) :: xx1,xx2,yy1,yy2     ! values of endpoints

  ! find points to interpolate between
  call find_points(xx,nx,x,x1,x2)
  call find_points(yy,ny,y,y1,y2)

  xx1=xx(x1)
  xx2=xx(x2)
  yy1=yy(y1)
  yy2=yy(y2)

 ! calculate weighted values
 sw=(xx2-x)*(yy2-y)
 nw=(xx2-x)*(y-yy1)
 se=(x-xx1)*(yy2-y)
 ne=(x-xx1)*(y-yy1)
 
 ! calculate weights (distances from point to each endpoint)
 !sw=(x-xx(x1))**2+(y-yy(y1))**2
 !sw=dsqrt(sw)
 !nw=(x-xx(x1))**2+(y-yy(y2))**2
 !nw=dsqrt(nw)
 !se=(x-xx(x2))**2+(y-yy(y1))**2
 !se=dsqrt(se)
 !ne=(x-xx(x2))**2+(y-yy(y2))**2
 !ne=dsqrt(ne)

 ! interpolate! (average weighted with distance to endpoints)
 z=(sw*zz(x1,y1)+nw*zz(x1,y2)+se*zz(x2,y1)+ne*zz(x2,x2))/((xx2-xx1)*(yy2-yy1))

end subroutine LIN_INT_2D

subroutine find_points(xx,n,x,x1,x2)
  !! find_points - Given an array xx of length n, sorted from lowest to highest
  !  values, find the two array elements, x1 and x2, that surround the value x.
  use prec_def
  implicit none

  integer, intent(in) :: n  ! length of array x
  real(Long), intent(in) :: xx(n) ! array of coordinates
  real(Long), intent(in) :: x     ! coordinate desired
  integer, intent(out) :: x1,x2   ! endpoints to use for lin_int

  integer :: xpiv                 ! pivot point

  ! check x is within xx
  if(x.lt.xx(1).or.x.gt.xx(n))then
   write(*,*)'LIN_INT: point is out of range of array. Cannot interpolate.'
   write(*,*)'xx_min,xx_max,x_desired:',xx(1),xx(n),x
   return
  endif

  ! set initial endpoints
  x1=1
  x2=n

  do while((x2-x1).gt.1)
    !move to new pivot
    xpiv=(x1+x2)/2
    !write(*,*)'xpiv:',xpiv
    if (x.lt.xx(xpiv)) then
      x2=xpiv
    elseif (xx(xpiv).lt.x) then
      x1=xpiv
    else
      x2=x1+1
!      write(*,*)'find_points: point is equal to grid point'
    endif
  enddo

end subroutine find_points
