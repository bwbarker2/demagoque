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
use bexception
use prec_def
implicit none

 private

 public :: dArrayList, iArrayList, new_dArrayList, new_iArrayList

 integer, parameter :: INITIAL_LENGTH=10

 type dArrayList
  private
  integer :: capacity
  integer :: sizeOfIt  !"OfIt" is so it is different from size below
  real(Long), dimension(:), allocatable  :: values
 contains
  procedure :: add => dArrayList_add
  procedure :: ensureCapacity => dArrayList_ensureCapacity
  procedure :: get => dArrayList_get
  procedure :: set => dArrayList_set
  procedure :: size => dArrayList_size
 end type dArrayList

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 type iArrayList
  private
  integer :: capacity
  integer :: sizeOfIt  !"OfIt" is so it is different from size below
  integer, dimension(:), allocatable  :: values
 contains
  procedure :: add => iArrayList_add
  procedure :: ensureCapacity => iArrayList_ensureCapacity
  procedure :: get => iArrayList_get
  procedure :: set => iArrayList_set
  procedure :: size => iArrayList_size
 end type iArrayList
 
contains

 !> constructs new dArrayList with optional initial length
 function new_dArrayList(initLength) result (self)
  implicit none

  integer, optional, intent(in) :: initLength !< initial length of list
  type (dArrayList) :: self
  real(Long), dimension(:), allocatable :: list_construct

  if(present(initLength)) then
   allocate(list_construct(initLength))
   self = dArrayList(initLength,0,list_construct)
  else
   allocate(list_construct(INITIAL_LENGTH))
   self = dArrayList(INITIAL_LENGTH,0,list_construct)
  endif

 end function new_dArrayList

 !> constructs new iArrayList with optional initial length
 function new_iArrayList(initLength) result (self)
  implicit none

  integer, optional, intent(in) :: initLength !< initial length of list
  type (iArrayList) :: self
  integer, dimension(:), allocatable :: list_construct

  if(present(initLength)) then
   allocate(list_construct(initLength))
   self = iArrayList(initLength,0,list_construct)
  else
   allocate(list_construct(INITIAL_LENGTH))
   self = iArrayList(INITIAL_LENGTH,0,list_construct)
  endif

 end function new_iArrayList

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !> Appends item to end of list
 subroutine dArrayList_add(self, newValue)
  implicit none

  class (dArrayList), intent(inout) :: self
  real(Long),            intent(in)    :: newValue ! value to add to list

  if(self%size()==self%capacity) then
   call dReallocate(self,self%capacity*2)
  endif

  self%values(self%size()+1)=newValue
  self%sizeOfIt=self%size()+1

 end subroutine dArrayList_add

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !> Appends item to end of list
 subroutine iArrayList_add(self, newValue)
  implicit none

  class (iArrayList), intent(inout) :: self
  integer,            intent(in)    :: newValue ! value to add to list

  if(self%size()==self%capacity) then
   call iReallocate(self,self%capacity*2)
  endif

  self%values(self%size()+1)=newValue
  self%sizeOfIt=self%size()+1

 end subroutine iArrayList_add

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !> Ensures that list can store 'newCapacity' items. Does not reduce
 !! capacity.
 subroutine dArrayList_ensureCapacity(self, newCapacity)
  implicit none

  class (dArrayList), intent(inout) :: self
  integer,           intent(in)    :: newCapacity !< capacity to ensure

  if(self%capacity<newCapacity) then
   call dReallocate(self, newCapacity)
  endif

 end subroutine dArrayList_ensureCapacity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !> Ensures that list can store 'newCapacity' items. Does not reduce
 !! capacity.
 subroutine iArrayList_ensureCapacity(self, newCapacity)
  implicit none

  class (iArrayList), intent(inout) :: self
  integer,           intent(in)    :: newCapacity !< capacity to ensure

  if(self%capacity<newCapacity) then
   call iReallocate(self, newCapacity)
  endif

 end subroutine iArrayList_ensureCapacity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !> Returns value at a given index.
 function dArrayList_get(self, indix) result(val)
  implicit none

  class(dArrayList), intent(in) :: self
  integer,         intent(in) :: indix  !< index of value to return
  real(Long)                      :: val

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

 !> Returns value at a given index.
 function iArrayList_get(self, indix) result(val)
  implicit none

  class(iArrayList), intent(in) :: self
  integer,           intent(in) :: indix  !< index of value to return
  integer                      :: val

  if(indix<=self%size()) then
   val=self%values(indix)
  else
   write(ERROR_UNIT,*)'iArrayList_get: index,size=',indix,self%size()
   call throwException('iArrayList_get: ArrayOutOfBoundsError, index,size=' &
                       ,BEXCEPTION_FATAL)
  endif

 end function iArrayList_get

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !> Sets a element to a value. That element must already have a value
 !! in it (size of the list must be >= index requested)
 subroutine dArrayList_set(self, indix, val)
  implicit none

  class(dArrayList), intent(inout) :: self 
  integer,          intent(in)    :: indix !< index of element to set
  real(Long),           intent(in)    :: val   !< new value for the element

  if(indix<=self%size()) then
   self%values(indix)=val
  else
   write(102,*)'dArrayList_set: ArrayOutOfBoundsError, indix,size=' &
               ,indix,self%size()
   write(102,*)'dArrayList_set: Doing nothing...'
  endif

 end subroutine dArrayList_set

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !> Sets a element to a value. That element must already have a value
 !! in it (size of the list must be >= index requested)
 subroutine iArrayList_set(self, indix, val)
  implicit none

  class(iArrayList), intent(inout) :: self 
  integer,           intent(in)    :: indix !< index of element to set
  integer,           intent(in)    :: val   !< new value for the element

  if(indix<=self%size()) then
   self%values(indix)=val
  else
   write(ERROR_UNIT,*)'iArrayList_set: indix,size=',indix,self%size()
   call throwException('iArrayList_set: ArrayOutOfBoundsError',BEXCEPTION_FATAL)
  endif

 end subroutine iArrayList_set

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !> Returns logical size of list.
 integer function dArrayList_size(self)
  implicit none

  class(dArrayList), intent(in) :: self

  dArrayList_size=self%sizeOfIt

 end function dArrayList_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !> Returns logical size of list.
 integer function iArrayList_size(self)
  implicit none

  class(iArrayList), intent(in) :: self

  iArrayList_size=self%sizeOfIt

 end function iArrayList_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !> Changes capacity of storage array. If newCapacity < current capacity,
 !! then it truncates, potentially losing information and reducing size.
 subroutine dReallocate(self, newCapacity)
  implicit none

  class (dArrayList), intent(inout) :: self
  integer,           intent(in)    :: newCapacity

  real(Long), dimension(self%capacity) :: newList

  if(newCapacity==self%capacity) return

  newList=self%values
  deallocate(self%values)
  allocate(self%values(newCapacity))
  if(newCapacity<self%capacity) then
   self%values=newList(1:newCapacity)
   if(newCapacity<self%sizeOfIt) self%sizeOfIt=newCapacity
  else
   self%values(1:self%capacity)=newList
  endif
  self%capacity=newCapacity

 end subroutine dReallocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !> Changes capacity of storage array. If newCapacity < current capacity,
 !! then it truncates, potentially losing information and reducing size.
 subroutine iReallocate(self, newCapacity)
  implicit none

  class (iArrayList), intent(inout) :: self
  integer,           intent(in)    :: newCapacity

  integer, dimension(self%capacity) :: newList

  if(newCapacity==self%capacity) return

  newList=self%values
  deallocate(self%values)
  allocate(self%values(newCapacity))
  if(newCapacity<self%capacity) then
   self%values=newList(1:newCapacity)
   if(newCapacity<self%sizeOfIt) self%sizeOfIt=newCapacity
  else
   self%values(1:self%capacity)=newList
  endif
  self%capacity=newCapacity

 end subroutine iReallocate

end module class_ArrayList
