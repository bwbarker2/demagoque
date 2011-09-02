module bmath
! Copyright (C) 2011  Brent W. Barker
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program (gpl-3.0.txt).  If not, see
!    <http://www.gnu.org/licenses/>.
!
!    Author: Brent W. Barker
!            barker at nscl dot msu dot edu
!            National Superconducting Cyclotron Laboratory
!            Michigan State University
!            1 Cyclotron, East Lansing, MI 48824-1321

contains

!logical function isPrime(num)
! implicit none
!
! integer, intent(in) :: num
!
! integer, parameter :: isPrime_u=1
! logical :: datexists
! integer :: numPrimes
! integer, dimension(:), allocatable :: primes
!
! integer :: i
!
! inquire(file='isPrime.dat', exist=datexists)
!
! if(datexists) then
!  open(unit=isPrime_u, file='isPrime.dat',status='old')
!  read(isPrime_u,*)numPrimes
!  allocate(primes(numPrimes))
!  do i=1,numPrimes
!   read(
! endif
! isPrime=.true.

! do i=2,num/2

!end function isPrime


!integer function findLargestPrimeFactor(num)
! implicit none
!
! integer, intent(in) :: num
!
! integer :: i,factor,remainder
!
! do i=2,num
!  remainder=mod(num,i)
!  if(remainder==0) then
!   factor=num/i
!   if(isPrime(factor))exit
!  endif
! enddo
!
! findLargestPrimeFactor=factor
!
!end function findLargestPrimeFactor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

complex*16 function zdet2d(cmat,n)
 !! zdet2d - nonoptimized calculation of determinant of 2D complex double matrix - order N^2
 implicit none

 integer,                    intent(in) :: n    !dimension of square matrix cmat
 complex*16, dimension(n,n), intent(in) :: cmat   !matrix

 integer, dimension(n) :: l  !row ordering after Gaussian elimination
 integer, dimension(n) :: lo !work array for ordering
 complex*16, dimension(n,n) :: wmat  !writable matrix for using in zGauss

 integer :: i !loop variables
 integer :: j !dummy variable for exchange
 real*8 :: rsign  !row interchange sign factor for determinant

 wmat=cmat

 !diagonalize matrix
 call zGauss(n,wmat,l)

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
  zdet2d=zdet2d*wmat(l(i),i)
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
   j=k
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine zlin_int(xa,ya,n,x,y,ki)
 !! zlin_int - linear interpolation for complex double arrays
 implicit none

 integer,                  intent(in) :: n  !size of xa,ya arrays
 real*8,     dimension(n), intent(in) :: xa !independent values
 complex*16, dimension(n), intent(in) :: ya !dependent values
 real*8,                   intent(in) :: x  !x value desired
 
 complex*16, intent(out)   :: y  !function value at x
 integer,    intent(inout) :: ki !index to interpolate around, or index of lower value used for interpolation (good for incremented calls)

 integer :: k,khi,klo
 real*8  :: h

 klo=1
 khi=n

 if(ki.ne.1.and.xa(ki).lt.x.and.xa(ki+1).gt.x) then
  klo=ki
  khi=klo+1
 elseif(ki.ne.1.and.xa(ki+1).lt.x .and. xa(ki+2).gt.x) then
  klo=ki+1
  khi=klo+1

 else
  do while (khi-klo.gt.1)
   k=(khi+klo)/2
   if(xa(k).gt.x) then
    khi=k
   else
    klo=k
   endif !xa.gt.x
  enddo !while khi-klo.gt.1
 endif !entire if

 !store lower index for future calls to this subroutine
 ki=klo

 h=xa(khi)-xa(klo)
 if(h.eq.0.d0) then
  write(*,*) 'bad xa input'
  read(*,*)
 endif

 y=ya(klo)+(ya(khi)-ya(klo))/(xa(khi)-xa(klo))*(x-xa(klo))

end subroutine zlin_int


end module bmath