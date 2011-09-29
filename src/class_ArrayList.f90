!> partial implementation of Java's "ArrayList", an array that has whatever 
!! size is necessary. One is able to append a value to the end of the array,
!! and it will grow if needed.
module class_ArrayList
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
implicit none

 private :: INITIAL_LENGTH, reallocate

 integer, parameter :: INITIAL_LENGTH=10

 type dArrayList
  private
  integer :: capacity
  integer :: sizeOfIt  !"OfIt" is so it is different from size below
  real*8, dimension(:), allocatable  :: values
 contains
  procedure :: add => dArrayList_add
  procedure :: ensureCapacity => dArrayList_ensureCapacity
  procedure :: get => dArrayList_get
  procedure :: set => dArrayList_set
  procedure :: size => dArrayList_size
 end type dArrayList
 
contains

 function make_dArrayList(initLength) result (self)
  implicit none

  integer, optional, intent(in) :: initLength
  type (dArrayList) :: self
  real*8, dimension(:), allocatable :: list_construct

  if(present(initLength)) then
   allocate(list_construct(initLength))
   self = dArrayList(initLength,0,list_construct)
  else
   allocate(list_construct(INITIAL_LENGTH))
   self = dArrayList(INITIAL_LENGTH,0,list_construct)
  endif

 end function make_dArrayList

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine dArrayList_add(self, newValue)
  implicit none

  class (dArrayList), intent(inout) :: self
  real*8,            intent(in)    :: newValue

  if(self%size()==self%capacity) then
   call reallocate(self,self%capacity*2)
  endif

  self%values(self%size()+1)=newValue
  self%sizeOfIt=self%size()+1

 end subroutine dArrayList_add

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine dArrayList_ensureCapacity(self, newCapacity)
  implicit none

  class (dArrayList), intent(inout) :: self
  integer,           intent(in)    :: newCapacity

  if(self%capacity<newCapacity) then
   call reallocate(self, newCapacity)
  endif

 end subroutine dArrayList_ensureCapacity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 function dArrayList_get(self, indix) result(val)
  implicit none

  class(dArrayList), intent(in) :: self
  integer,         intent(in) :: indix
  real*8                      :: val

  if(indix<=self%size()) then
   val=self%values(indix)
  else
   write(102,*)'dArrayList_get: ArrayOutOfBoundsError, index,size=', &
                indix,self%size()
   write(102,*)'dArrayList_get: Returning value -1'
   val=-1
  endif

 end function dArrayList_get

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine dArrayList_set(self, indix, val)
  implicit none

  class(dArrayList), intent(inout) :: self
  integer,          intent(in)    :: indix
  real*8,           intent(in)    :: val

  if(indix<=self%size()) then
   self%values(indix)=val
  else
   write(102,*)'dArrayList_set: ArrayOutOfBoundsError, indix,size=' &
               ,indix,self%size()
   write(102,*)'dArrayList_set: Doing nothing...'
  endif

 end subroutine dArrayList_set

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 integer function dArrayList_size(self)
  implicit none

  class(dArrayList), intent(in) :: self

  dArrayList_size=self%sizeOfIt

 end function dArrayList_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine reallocate(self, newCapacity)
  implicit none

  class (dArrayList), intent(inout) :: self
  integer,           intent(in)    :: newCapacity

  real*8, dimension(self%capacity) :: newList

  if(newCapacity==self%capacity) return

  newList=self%values
  deallocate(self%values)
  allocate(self%values(newCapacity))
  if(newCapacity<self%capacity) then
   self%values=newList(1:newCapacity)
  else
   self%values(1:self%capacity)=newList
  endif
  self%capacity=newCapacity

 end subroutine reallocate


end module class_ArrayList
