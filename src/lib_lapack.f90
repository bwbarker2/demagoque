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

module lib_lapack

 contains

  subroutine getEigenSq(mat,num,evals,evecs)

   complex*16, dimension(0:num-1,0:num-1), intent(inout) :: mat ! input matrix
   integer,                                intent(in) :: num  ! dimension and rank of matrix
   complex*16, dimension(0:num-1),         intent(out) :: evals !eigenvalues
   complex*16, dimension(0:num-1,0:num-1), intent(inout) :: evecs !eigenvectors, ith column corresponds to ith eigenvalue in evals array

   complex*16, dimension(0:num-1,0:num-1) :: mat2 ! copy of input matrix
   integer                           :: lwork  !length of work array
   integer                           :: info   !return statement from zgeev
   complex*16, dimension(1)          :: vl  !left eigenvalues, not calculated
   complex*16, allocatable, dimension(:) :: work  !working array
   complex*16, dimension(0:2*num-1)  :: rwork      !other working array

   !initialize values for first ZGEEV call
   mat2=mat
   lwork=-1
   info=-99  !arbitrary value
   allocate(work(0:0))

   ! calculates correct dimension of working array 'work'
   call zgeev('N','V',num,mat2,num,evals,vl,1,evecs,num,work,lwork,rwork,info)
   if(info.ne.0)write(*,*)'error with zgeev, info=',info 

!   write(*,*)'info=',info
!   write(*,*)'lwork should be',work(0)
!   write(*,*)

   lwork=nint(dble(work(0)))
   deallocate(work)
   allocate(work(0:lwork-1))

   ! calculates eigenvalues/vectors
   call zgeev('N','V',num,mat2,num,evals,vl,1,evecs,num,work,lwork,rwork,info)
   if(info.ne.0)write(*,*)'error with zgeev, info=',info

!   write(*,*)'evals='
!   write(*,*)evals
!   write(*,*)
!   write(*,*)'evecs='
!   write(*,*)evecs
!   write(*,*)

  end subroutine getEigenSq

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine getInvMat(mat,num,matinv)
   implicit none

   complex*16, dimension(0:num-1,0:num-1), intent(in) :: mat !matrix to invert
   integer, intent(in) :: num
   complex*16, dimension(0:num-1,0:num-1), intent(out) :: matinv !inverted matrix

   integer, dimension(0:num-1) :: ipiv ! pivot array for LU factorization
   integer :: info ! return status of zgetrf
   complex*16, allocatable, dimension(:) :: work !work array for zgetri
   integer :: lwork !size of work array for zgetri

   matinv=mat
   info=-99

   !factors matrix into LU
   call zgetrf(num,num,matinv,num,ipiv,info)
   if(info.ne.0)write(*,*)'error with zgetrf, info=',info

!   write(*,*)'LU='
!   write(*,*)matinv
!   write(*,*)

   allocate(work(0:0))

   lwork=-1

   !calculates inverse using LU factorizations
   call zgetri(num,matinv,num,ipiv,work,lwork,info)
   if(info.ne.0)write(*,*)'error with zgetri, info=',info

!   write(*,*)'lwork should be',work(0)
 
   lwork=nint(dble(work(0)))

   deallocate(work)
   allocate(work(0:lwork-1))

   call zgetri(num,matinv,num,ipiv,work,lwork,info)
   if(info.ne.0)write(*,*)'error with zgetri, info=',info

!   write(*,*)'matinv:'
!   write(*,*)matinv
!   write(*,*)

  end subroutine getInvMat

end module lib_lapack
