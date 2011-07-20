module class_ArrayList
 implicit none

 private :: INITIAL_LENGTH, reallocate

 integer, parameter :: INITIAL_LENGTH=10

 type dArrayList
  private
  integer :: capacity
  integer :: size
  real*8, dimension(:), allocatable  :: values
 end type dArrayList
 
contains

 function make_dArrayList(initLength) result (name)
  implicit none

  integer, optional, intent(in) :: initLength
  type (dArrayList) :: name
  real*8, dimension(:), allocatable :: list_construct

  if(present(initLength)) then
   allocate(list_construct(initLength))
   name = dArrayList(initLength,0,list_construct)
  else
   allocate(list_construct(INITIAL_LENGTH))
   name = dArrayList(INITIAL_LENGTH,0,list_construct)
  endif

 end function make_dArrayList

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine dArrayList_add(list, newValue)
  implicit none

  type (dArrayList), intent(inout) :: list
  real*8,            intent(in)    :: newValue

  if(list%size==list%capacity) then
   call reallocate(list,list%capacity*2)
  endif

  list%values(list%size+1)=newValue
  list%size=list%size+1

 end subroutine dArrayList_add

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine dArrayList_ensureCapacity(list, newCapacity)
  implicit none

  type (dArrayList), intent(inout) :: list
  integer,           intent(in)    :: newCapacity

  if(list%capacity<newCapacity) then
   call reallocate(list, newCapacity)
  endif

 end subroutine dArrayList_ensureCapacity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 function dArrayList_get(list, indix) result(val)
  implicit none

  type(dArrayList), intent(in) :: list
  integer,         intent(in) :: indix
  real*8                      :: val

  if(indix<=list%size) then
   val=list%values(indix)
  else
   write(102,*)'dArrayList_get: ArrayOutOfBoundsError, index,size=', &
                indix,list%size
   write(102,*)'dArrayList_get: Returning value -1'
   val=-1
  endif

 end function dArrayList_get

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine dArrayList_set(list, indix, val)
  implicit none

  type(dArrayList), intent(inout) :: list
  integer,          intent(in)    :: indix
  real*8,           intent(in)    :: val

  if(indix<=list%size) then
   list%values(indix)=val
  else
   write(102,*)'dArrayList_set: ArrayOutOfBoundsError, indix,size=' &
               ,indix,list%size
   write(102,*)'dArrayList_set: Doing nothing...'
  endif

 end subroutine dArrayList_set

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 integer function dArrayList_size(list)
  implicit none

  type(dArrayList), intent(in) :: list

  dArrayList_size=list%size

 end function dArrayList_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine reallocate(list, newCapacity)
  implicit none

  type (dArrayList), intent(inout) :: list
  integer,           intent(in)    :: newCapacity

  real*8, dimension(list%capacity) :: newList

  if(newCapacity==list%capacity) return

  newList=list%values
  deallocate(list%values)
  allocate(list%values(newCapacity))
  if(newCapacity<list%capacity) then
   list%values=newList(1:newCapacity)
  else
   list%values(1:list%capacity)=newList
  endif
  list%capacity=newCapacity

 end subroutine reallocate


end module class_ArrayList
