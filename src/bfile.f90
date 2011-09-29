!> Common file operations
module bfile
 implicit none

! I want to put a linked list here of free file units. Maybe an ArrayList is fine?

contains

 !> closes FUnit and frees up unit number
 subroutine closeFUnit(funit)
  implicit none


 end subroutine closeFUnit

 !> returns next un-used file unit
 integer function getNextFUnit()
  implicit none
  
 end function getNextFUnit

end module bfile
