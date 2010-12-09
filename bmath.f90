complex*16 function zdet2d(cmat,n)
 !! zdet2d - nonoptimized calculation of determinant of 2D complex double matrix - order N^2
 implicit none

 integer,                    intent(in) :: n    !dimension of square matrix cmat
 complex*16, dimension(n,n), intent(in) :: cmat   !matrix

 integer, dimension(n) :: l  !row ordering after Gaussian elimination
 integer, dimension(n) :: lo !work array for ordering

 integer :: i !loop variables
 integer :: j !dummy variable for exchange
 real*8 :: rsign  !row interchange sign factor for determinant

 !diagonalize matrix
 call zGauss(n,cmat,l)

 lo=l

 !determine number of row interchanges to get actual diagonal matrix - this
 !determines sign of determinant.
 rsign=1
! write(*,*)lo
 do i=1,n-1
  if(lo(i).ne.i) then
   j=lo(lo(i))
   lo(lo(i))=lo(i)
   lo(i)=j
   rsign=rsign*(-1)
  endif
! write(*,*)lo
 enddo

! write(*,*)'rsign:',rsign

 !multiply diagonal terms together to get determinant
 zdet2d=dsign(1.d0,rsign)
 do i=1,n
  zdet2d=zdet2d*cmat(l(i),i)
 enddo

end function zdet2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

logical function isOdd(num)
 implicit none

 integer, intent(in) :: num

 if(MOD(num,2)==1) then
  isOdd=.true.
 else
  isOdd=.false.
 endif

end function isOdd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine zGauss(n,a,l)
 !! zGauss - performs Gaussian Elimination with scaled partial pivoting as
 !           described in W.~Cheney, D.~Kincaid, Numerical Mathematics and
 !           Computing, 5th ed., Thomson (2004). Extended for complex values.
 !
 ! output: The output matrix 'a' is a triangular matrix with the rows out of
 ! order. The row order is specified in the array 'l'. Scaling factors used
 ! in the process of the Gaussian elimination are stored in 'a' in the elements
 ! that are eliminated. For example, if L=[3,1,2], then row 3 is the 'first'
 ! row, and all other rows have '0' in the first column. But instead of a '0',
 ! the value stored is the multiplier. In the above case, the diagonal elements
 ! are L(3,1), L(1,2), and L(2,3], and the row elements before the diagonal
 ! elements have been eliminated.
 !
 ! This format reduces the amount of redundant copying and memory needed for
 ! storing scaling factors needed for solving a system of linear equations.
 implicit none

 integer,                    intent(in)    :: n !dimension of square matrix
 complex*16, dimension(n,n), intent(inout) :: a !matrix
 integer,    dimension(n),   intent(out)   :: l !index array

 integer                  :: i,j,k  !loop variables
 complex*16               :: r,rmax,smax,xmult
 complex*16, dimension(n) :: s
 
 do i=1,n
   l(i)=i
   smax=0.d0
   do j=1,n
     if(abs(a(i,j))>abs(smax)) smax=a(i,j)
   enddo
   s(i)=smax
 enddo

 do k=1,n-1
   rmax=0.d0
   do i=k,n
     r=a(l(i),k)/s(l(i))
     if(abs(r)>abs(rmax)) then
       rmax=r
       j=i
     endif
   enddo
   i=l(j)
   l(j)=l(k)
   l(k)=i
   do i=k+1,n
     xmult=a(l(i),k)/a(l(k),k)
     a(l(i),k)=xmult  !where zeroes would be, put multiplier instead
     do j=k+1,n
       a(l(i),j)=a(l(i),j)-xmult*a(l(k),j)
     enddo
   enddo
 enddo

end subroutine zGauss
