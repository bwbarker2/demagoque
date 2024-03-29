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
 use specfun
 use prec_def
 implicit none

! private :: calck

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

!> Calculates modified Bessel function of the second kind (K) for any integer order and real argument.
!! Uses specfun.f90 to get besk0 and besk1, then uses recurrence to get other orders.
real (Long) function bmath_besk(order,arg)
 implicit none

 integer, intent(in) :: order !< order of function (subscript of K)
 real (Long), intent(in) :: arg !< argument of function
 
 real (Long) :: bk0, bk1 !< calculation of besselk0 and k1

 if(order==0) then
  bmath_besk = besk0(arg)
  return
 endif

 if(order==1) then
  bmath_besk = besk1(arg)
  return
 endif

 bk0 = besk0(arg)
 bk1 = besk1(arg)

 ! if not 0 or 1, then use recursion relation
 bmath_besk=calcbesk(order,arg,bk0,bk1)

end function bmath_besk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Calculates first derivative of modified Bessel function of the second kind (K) for any integer order,
!! using recursion relation.
real (Long) function bmath_dbeskdx(order,arg) result (fres)
 implicit none

 integer, intent(in) :: order !< order of function (subscript of K)
 real (Long), intent(in) :: arg !< argument of function

 fres=-0.5_Long*(bmath_besk(order-1,arg)+bmath_besk(order+1,arg))

end function bmath_dbeskdx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Calculates modified Bessel function of the second kind (K) for any integer order, given the values of
!! orders 1 and 2.
real (Long) recursive function calcbesk(order,arg,bk0,bk1) result(bes)
 implicit none

 integer, intent(in) :: order !< order of function
 real (Long), intent(in) :: arg !< argument of function
 real (Long), intent(in) :: bk0 !< value of K_0(arg)
 real (Long), intent(in) :: bk1 !< value of K_1(arg)

 if(order>1) then
  bes =  2._Long*(order-1)/arg*calcbesk(order-1,arg,bk0,bk1) + calcbesk(order-2,arg,bk0,bk1)
 elseif(order<0) then
  bes = -2._Long*(order+1)/arg*calcbesk(order+1,arg,bk0,bk1) + calcbesk(order+2,arg,bk0,bk1)
 elseif(order==0) then
  bes=bk0
 else !if(order==1) then
  bes=bk1
 endif

end function calcbesk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Swaps two real(Long) values
subroutine dswap(a,b)
 implicit none

 real(Long), intent(inout) :: a,b !< variables to swap

 real(Long) :: c

 c=a
 a=b
 b=c

end subroutine dswap

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Finds if a is between b and c
logical function bmath_dIsBetween(a,bin,cin) result(ans)
 implicit none

 real(Long), intent(in) :: a !< value to test
 real(Long), intent(in) :: bin,cin !< bounds

 real(Long) :: b,c

 b=bin
 c=cin

 if(b>c)then
  call dswap(b,c)
 endif

 if ( (b<a) .and. (a<c) ) then
  ans=.true.
 else
  ans=.false.
 endif

end function bmath_dIsBetween


!> Finds the zero of a function f between a and b.
!! f(a) and f(b) must have opposite signs.
!! Official reference: R.P. Brent (1973). Algorithms for Minimization without Derivatives, Chapter 4. Prentice-Hall, Englewood Cliffs, NJ. ISBN 0-13-022335-2.
!! My reference: http://en.wikipedia.org/w/index.php?title=Brent%27s_method&oldid=457696075#Algorithm .
function bmath_dZeroBrent(ain,bin,f,iflag,err,maxiterin) result(zero)
 use bexception
 implicit none

 real (Long) :: zero ! return variable
 real (Long), intent(in) :: ain !< one bound of search interval
 real (Long), intent(in) :: bin !< the other bound
 real (Long), external :: f   !< function to find roots of
 integer, optional, intent(out) :: iflag !< error flag, 0 = no errors
 real (Long), optional, intent(in) :: err !< convergence tolerance, default=1e-10
 integer, optional, intent(in) :: maxiterin !< maximum iterations, default=1000

 real (Long) :: fa,fb,fc,fs !< function evaluated at a,b,c,s

 real (Long) :: a,b,c,d,s !< trial positions

 logical :: mflag

 real (Long) :: erruse  !convergence interval to use
 integer :: maxiters !maximum iterations to use
 integer :: iters  !number of iterations used

 character(len=130) :: errmsg

 zero=snan

 if(present(iflag))iflag=0

 if(present(err)) then
  erruse=err
 else
!  erruse=10.0**(-Long)
  erruse=epzero*2._Long**3
 endif

 if(present(maxiterin)) then
  maxiters=maxiterin
 else
  maxiters=1000
 endif

 a=ain
 b=bin

 fa=f(a)
 fb=f(b)
 fs=fb

 if (fa*fb>=0) then
  write(ERROR_UNIT,*)'bmath_dZeroBrent: a,f(a),b,f(b)=',a,fa,b,fb
  write(errmsg,*)'bmath_dZeroBrent: root is not bracketed'
  call throwException(errmsg,BEXCEPTION_FATAL)
  if(present(iflag))iflag=1
  return
 endif

 if (abs(fa) < abs(fb)) then
  call dswap(a,b)
  call dswap(fa,fb)
 endif

 c=a
 fc=fa

 mflag=.true.

 iters=0

 d=0e0_Long !set this to suppress compiler warning

 do while ((fb /= 0) .and. (fs /= 0) .and. abs(b-a)>erruse)
  iters=iters+1
  if(iters>maxiters) call throwException('bmath_dZeroBrent: iterations exceeded maximum',BEXCEPTION_FATAL)

  if ( (fa/=fc) .and. (fb/=fc) ) then
   s= a*fb*fc/((fa-fb)*(fa-fc)) &
     +b*fc*fa/((fb-fa)*(fb-fc)) &
     +c*fa*fb/((fc-fa)*(fc-fb))
!   write(ERROR_UNIT,*)'tried quad'
  else
   s=b-fb*(b-a)/(fb-fa)
!   write(ERROR_UNIT,*)'tried secant'
  endif

  if(     .not.bmath_dIsBetween(s,(3_Long*a+b)*0.25_Long,b) &
     .or. (     mflag .and. abs(s-b)>=abs(b-c)*0.5_Long) &
     .or. (.not.mflag .and. abs(s-b)>=abs(c-d)*0.5_Long) &
     .or. (     mflag .and. abs(b-c)<abs(erruse)) &
     .or. (.not.mflag .and. abs(c-d)<abs(erruse)) ) then
   s = (a+b)*0.5_Long
   mflag=.true.
!   write(ERROR_UNIT,*)'...but used bisection'
  else
   mflag=.false.
!   write(ERROR_UNIT,*)'...and used it'
  endif

  fs=f(s)
  d=c
  c=b
  fc=fb

  if(fa*fs<0) then
   b=s
   fb=fs
  else
   a=s
   fa=fs
  endif

  if (abs(fa) < abs(fb)) then
   call dswap(a,b)
   call dswap(fa,fb)
  endif
!  write(ERROR_UNIT,*)'s:',s
 end do

 zero=s

end function bmath_dZeroBrent

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Returns the inverse of a matrix calculated by finding the LU
!! decomposition.  Depends on LAPACK.
!!
!! Original code public domain, by Jason Blevins, found at
!! http://fortranwiki.org/fortran/show/inv
!!
function bmath_inv(A) result(Ainv)
  real(long), dimension(:,:), intent(in) :: A
  real(long), dimension(size(A,1),size(A,2)) :: Ainv

  real(long), dimension(size(A,1)) :: work  ! work array for LAPACK
  integer, dimension(size(A,1)) :: ipiv   ! pivot indices
  integer :: n, info

  ! External procedures defined in LAPACK
!  external DGETRF
!  external DGETRI

  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  Ainv = A
  n = size(A,1)

  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  call DGETRF(n, n, Ainv, n, ipiv, info)

  if (info /= 0) then
     stop 'Matrix is numerically singular!'
  end if

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  call DGETRI(n, Ainv, n, ipiv, work, n, info)

  if (info /= 0) then
     stop 'Matrix inversion failed!'
  end if
end function bmath_inv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Finds the first derivative of a function of a single variable.
!!
!! Uses the 3rd order finite different formula with Richardson's extrapolation
!! to gain nth order convergence
!! Source: W.~Cheney, D.~Kincaid, Numerical Mathematics and
!! Computing, 5th ed., Thomson (2004).
!!
!! Can be called with optional arguments to set minimum error tolerance
!! and maximum number of function calls.
real(Long) function bmath_LDiffRichardson( &
                    ff,xx,hin,derr,nmaxin,errin,istat) result(dfdx)
 use bexception
 implicit none

 interface
  real(Long) function ff(xxx)
   use prec_def
   implicit none
   real(Long), intent(in) :: xxx
  end function ff
 end interface

! real (Long), external         :: ff     !< function to differentiate
 real (Long)      ,intent(in)  :: xx     !< argument of function
 real (Long)      ,intent(out) :: derr   !< uncertainty of result
 real (Long)      ,intent(in)  :: hin    !< initial step, set to ~1/10 of 
                                         !! length scale
 integer ,optional,intent(in)  :: nmaxin !< max steps, default 25
 real (Long)  ,optional,intent(in)  :: errin  !< tolerance, default epzero
 integer ,optional,intent(inout) :: &
  istat  !< error flag. 0=success, 1=error above tolerance and algorithm is
         !! diverging, 2=hit max iterations without reaching
         !! desired error and is still converging

 real (Long), allocatable, dimension(:,:) :: dd !< where the different terms are stored
 integer :: ii,jj,nmax
 real (Long) :: hh,err
 real (Long) :: try

 if(present(errin)) then
  err=errin
 else
  err=epzero
 endif

 ! make sure stepsize doesn't go below epzero
 if(present(nmaxin)) then
  nmax=min(nmaxin,nint(log(hin)-log(epzero)-1._Long))
 else
  nmax=nint(log(hin)/log(2._Long)-log(epzero)/log(2._Long)-1._Long)
 endif

 if(present(istat)) istat=0

 allocate(dd(0:nmax,0:nmax))

 derr=maxdouble
 try=derr
 hh=hin
 dfdx=snan

 dd(0,0)=(ff(xx+hh)-ff(xx-hh))/(2._Long*hh)
 hh=hh*0.5_Long

 do ii=1,nmax

  dd(ii,0)=(ff(xx+hh)-ff(xx-hh))/(2._Long*hh)

  do jj=0,ii-1

   dd(ii,jj+1)=dd(ii,jj)+(dd(ii,jj)-dd(ii-1,jj))/real(4**(jj+1)-1,Long)
   try=max(abs(dd(ii,jj+1)-dd(ii,jj)),abs(dd(ii,jj+1)-dd(ii-1,jj)))

   if(try<=derr) then
    derr=try
    dfdx=dd(ii,jj+1)
    if(derr<=err) return   ! if desired tolerance is reached, exit early
   endif

  end do

  hh = hh*0.5_Long

  !if relative error starts to get worse instead of better
  if(abs(dd(ii,ii)-dd(ii-1,ii-1))>=2._Long*derr) then
   !if desired error is specified
   if(present(errin)) then
    !and our error is greater than that
    if(derr>errin) then
     !then either set istat or throw an exception
     if(present(istat)) then
      istat=1
     else
      call throwException('bmath_LDiffRichardson: error tolerance not reached'&
                           //' and algorithm no longer converging' &
                          ,BEXCEPTION_FATAL)
     endif !present istat
    endif !derr>errin
   else !if no error tolerance given
    return  ! just exit early
   endif !present errin
  endif !relative error
 end do !ii

 if(present(istat)) then
  istat=2
 else
  call throwException('bmath_LDiffRichardson: error above desired tolerance' &
                      //' and max iterations reached '&
                      //' (suggest increasing nmaxin).',BEXCEPTION_FATAL)
 endif

end function bmath_LDiffRichardson

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

complex(Long) function zdet2d(cmat,n)
 !! zdet2d - nonoptimized calculation of determinant of 2D complex double matrix - order N^2
 implicit none

 integer,                    intent(in) :: n    !dimension of square matrix cmat
 complex (Long), dimension(n,n), intent(in) :: cmat   !matrix

 integer, dimension(n) :: l  !row ordering after Gaussian elimination
 integer, dimension(n) :: lo !work array for ordering
 complex (Long), dimension(n,n) :: wmat  !writable matrix for using in zGauss

 integer :: i !loop variables
 integer :: j !dummy variable for exchange
 real (Long) :: rsign  !row interchange sign factor for determinant

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
 zdet2d=sign(1.0_Long,rsign)
 do i=1,n
  zdet2d=zdet2d*wmat(l(i),i)
 enddo

end function zdet2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

logical function isEven(num)
 implicit none

 integer, intent(in) :: num

 isEven=.not.isOdd(num)

end function isEven

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
 complex (Long), dimension(n,n), intent(inout) :: a !matrix
 integer,    dimension(n),   intent(out)   :: l !index array

 integer                      :: i,j,k  !loop variables
 complex (Long)               :: r,rmax,smax,xmult
 complex (Long), dimension(n) :: s
 
 do i=1,n
   l(i)=i
   smax=0.0_Long
   do j=1,n
     if(abs(a(i,j))>abs(smax)) smax=a(i,j)
   enddo
   s(i)=smax
 enddo

 do k=1,n-1
   rmax=0.0_Long
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

 integer,                      intent(in) :: n  !size of xa,ya arrays
 real (Long),    dimension(n), intent(in) :: xa !independent values
 complex (Long), dimension(n), intent(in) :: ya !dependent values
 real (Long),                  intent(in) :: x  !x value desired
 
 complex (Long), intent(out)   :: y  !function value at x
 integer,        intent(inout) :: ki !index to interpolate around, or index of lower value used for interpolation (good for incremented calls)

 integer     :: k,khi,klo
 real (Long) :: h

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
 if(h.eq.0.0_Long) then
  write(*,*) 'bad xa input'
  read(*,*)
 endif

 y=ya(klo)+(ya(khi)-ya(klo))/(xa(khi)-xa(klo))*(x-xa(klo))

end subroutine zlin_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Calculates uncertainty of function, given uncertainties of parameters.
!! Assumes that errors are random and normally distributed, and that the
!! error given represents one standard deviation. Result is one standard
!! deviation, assuming values make a normal distribution.
!!
!! \TODO sometimes this takes forever or hangs. Debug
real(Long) function bmath_propagateError(fcn,params,errs) result(propErr)
 use bexception

 real (Long), dimension(:), intent(in) :: params !< list of parameters
 real (Long), dimension(:), intent(in) :: errs   !< list of uncertainties of parameters

 real (Long), dimension(:), allocatable :: params1
 real (Long) :: fcn1,ave,var2,varvar,avevar,var,jvar2,jave

 integer :: ii,jj,kk

 interface
  real(Long) function fcn(params)
   use prec_def
   implicit none
   real(Long), dimension(:), intent(in) :: params
  end function fcn
 end interface

 if(size(params)/=size(errs)) then
  call throwException('bmath_propagateError: number of parameters must equal' &
   //' number of errors',BEXCEPTION_FATAL)
 endif

 allocate(params1(size(params)))
 propErr = 0._Long

 jvar2=0._Long
 jave=0._Long
 do kk=1,1000!00
  ave = 0._Long
  var2 = 0._Long
  varvar=0._Long
  avevar=0._Long
  do jj=1,100000
   do ii=1,size(params)
    params1(ii)=bmath_rand_normal(params(ii),errs(ii))
   enddo
 
   fcn1 = fcn(params1)
 
 
   ! calculate running sig^2 and ave to prevent large-number errors. (from https://en.wikipedia.org/w/index.php?title=Standard_deviation&oldid=541302131)
   var2 = var2 + (jj-1._Long)/jj*(fcn1-ave)**2
   ave = ave + (fcn1 - ave)/jj

   var = sqrt(var2/jj)

   varvar=varvar + (jj-1._Long)/jj*(var-avevar)**2
   avevar=avevar + (var - avevar) / jj
 
 !  write(error_unit,*)jj,fcn1,var,sqrt(varvar/jj)/var
 
   if(jj>=10) then
    ! if variance of the variances is low enough, then call it converged
    if(sqrt(varvar/jj)/var<=0.01) then
!     propErr=var
 !    write(error_unit,*)'jj=',jj
     exit
    endif
   endif !jj>1
!   prevvar=var
  enddo !jj
  jvar2 = jvar2 + (kk - 1._Long) / kk * (var - jave)**2
  jave = jave + (var - jave) / kk

  ! if relative standard deviation of the mean is small enough
  if(sqrt(jvar2)/(kk*jave)<0.0000001) then
   propErr=var
   return
  endif
 enddo !kk

 call throwException('bmath_propagateError: max iterations reached: error not certain',BEXCEPTION_WARNING)
 propErr=var
end function bmath_propagateError

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Returns random number from the normal distribution with PDF
!!   phi(x) = exp(-(x-mean)^2/(2*sig**2))
!! using Marsaglia polar method and PRNG from Fortran 95 intrinsic.
!! (implementation converted from Java here: https://en.wikipedia.org/w/index.php?title=Marsaglia_polar_method&oldid=539981976 )
!! 
!! \TODO replace call to random_number with custom PRNG for more repeatability, since
!!       it looks like implementation is compiler-specific.
!!
!! \TODO for faster algorithm, implement Ziggurat Algorithm
real(Long) function bmath_rand_normal(mean,sig) result(rand)
 use phys_cons

 real(Long), intent(in) :: mean !< mean of PDF of normal distribution
 real(Long), intent(in) :: sig !< width of PDF of normal distribution

 real(Long), save :: spare
 logical, save :: isSpareReady = .false.

 real (Long) :: u,v,s,mul

 if(isSpareReady) then
  isSpareReady=.false.
  rand = spare * sig + mean
 else
  s = 2._Long
  do while ((s>=1._Long).or.(s<epzero))
   call random_number(u)
   call random_number(v)
   u = u * 2._Long - 1._Long
   v = v * 2._Long - 1._Long
   s = u*u + v*v
  enddo
  mul = sqrt(-2._Long*log(s)/s)
  spare = v * mul
  isSpareReady = .true.
  rand = mean + sig * u * mul
 endif

end function bmath_rand_normal
 
end module bmath
