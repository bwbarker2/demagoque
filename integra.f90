SUBROUTINE dint_simp1(n, f, h, sum, err)

  !***************************************************************************80
  !! int_simp1 - integrates with Composite Simpson's Rule, double precision, 1D
  !
  ! License:
  !
  !    Copyright (C) 2009  Brent W. Barker
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
  !    along with this program.  If not, see <http://www.gnu.org/licenses/>.
  !
  !    The error given is the error used to guide an adaptive Simpson's Rule,
  !    and I'm not sure what relation it has to a standard deviation.
  !
  !  Modified:
  !
  !    8 September 2009
  !
  !  Author:
  !
  !    Brent W. Barker
  !
  !  Reference:
  !
  !    W. Cheney, D. Kincaid, Numerical Mathematics and Computing, 5th ed.,
  !    Brooks/Cole-Thomson, pages 237-41, 2004.
  !
  !  Parameters:
  !
  !    Input, integer n, length of array to be integrated. Must satisfy
  !    MOD(n-1,4)=0
  !  
  !    Input, real (Long) f, array to be integrated.
  !
  !    Input, real (Long) h, interval between 2 elements of array.
  !
  !    Output, real (Long) sum, result of integration
  !
  !    Output, real (Long) err, uncertainty in result
  IMPLICIT NONE

  INTEGER , PARAMETER :: Long = 8

  INTEGER, INTENT(in) :: n
  REAL (Long), INTENT(in) :: f(n)
  REAL (Long), INTENT(in) :: h
  REAL (Long), INTENT(out) :: sum
  REAL (Long), INTENT(out) :: err

  REAL (Long) :: sum2  ! integral with half the points
  INTEGER :: i

  sum = f(1) + f(n)
  sum2 = sum
  err = 0

  ! calculate sum
  DO i=3,n-2,2
     sum = sum + 2*f(i)
  ENDDO

  DO i=2,n-1,2
     sum = sum + 4*f(i)
  ENDDO

  sum = sum * h/3


  ! calculate sum2
  DO i=5,n-4,4
     sum2 = sum2 + 2*f(i)
  ENDDO

  DO i=3,n-2,4
     sum2 = sum2 + 4*f(i)
  ENDDO

  sum2 = sum2 * 2*h/3

  err = ( sum - sum2 ) / 15

  sum = ( 16 * sum - sum2 ) / 15

  err = dabs(err)

END SUBROUTINE dint_simp1
