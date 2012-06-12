!> Stores a superposition of potentials of type Potential and is 
!! itself a Potential
module class_PotentialList
 use class_Potential
 use prec_def
 implicit none

 private

 public :: new_PotentialList, PotentialList

 !> Superposition of potentials
 type, extends(Potential) :: PotentialList
  private
  class(Link_Potential), pointer :: firstLink => null()  !< first link in linked list
  class(Link_Potential), pointer :: lastLink => null()   !< last link
 contains
  procedure,public :: add => PotentialList_add  !< add a potential to the list
  procedure,public :: potV_denx => potV_denx !< get sum of potentials
  procedure,public :: potV_x => potV_x !< get sum of potentials
 end type PotentialList

 !> Individual link in linked list
 type :: Link_Potential
  private
  class(Potential), pointer :: myPot => null()     !< ptr to potential
  class(Link_Potential), pointer :: next => null() !< ptr to next link
 contains
  private
  procedure :: getPot => PotentialList_getPot      !< get potential
  procedure :: nextLink => PotentialList_nextLink  !< returns ptr to next link
  procedure :: setNextLink => PotentialList_setNextLink !< sets 'next'
 end type Link_Potential

contains

 function PotentialList_getPot(this) result(getPot)
  class(Potential), pointer :: getPot
  class(Link_Potential), intent(in) :: this

  getPot => this%myPot
 end function PotentialList_getPot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 function PotentialList_nextLink(this) result(nextLink)
  class(Link_Potential), pointer :: nextLink

  class(Link_Potential), intent(in) :: this

  nextLink => this%next
 end function PotentialList_nextLink

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 function new_Link_Potential(newPot,next)
  class(Link_Potential), pointer :: new_Link_Potential

  class(Potential), intent(in) :: newPot
  class(Link_Potential), pointer, intent(in) :: next

  allocate(new_Link_Potential)
  new_Link_Potential%next => next
  allocate(new_Link_Potential%myPot, source=newPot)

 end function new_Link_Potential

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine PotentialList_setNextLink(this,next)

  class(Link_Potential), intent(inout) :: this
  class(Link_Potential), pointer, intent(in) :: next

  this%next => next
 end subroutine PotentialList_setNextLink

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 pure function new_PotentialList() result(new)
  implicit none

  class(PotentialList), pointer :: new

  allocate(new)

 end function new_PotentialList

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine PotentialList_add(this,pottoadd)
  implicit none

  class(PotentialList), intent(inout) :: this
  class(Potential), intent(in) :: pottoadd

  class(Link_Potential), pointer :: newLink

!write(ERROR_UNIT,*)'Entering PotentialList_add'

  if (.not. associated(this%firstLink)) then
   this%firstLink => new_Link_Potential(pottoadd,this%firstLink)
   this%lastLink => this%firstLink
  else
   newLink => new_Link_Potential(pottoadd, this%lastLink%nextLink())
   call this%lastLink%setNextLink(newLink)
   this%lastLink => newLink
  endif

!write(ERROR_UNIT,*)'Leaving PotentialList_add'

 end subroutine PotentialList_add

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 real(Long) function potV_denx(this,rho) result(pot)
  implicit none

  class(PotentialList), intent(in) :: this
  real (Long), intent(in) :: rho

  class(Link_Potential), pointer :: curr

  pot=0._Long

  curr => this%firstLink

  do while(associated(curr))
   pot=pot+curr%myPot%potV_denx(rho)
   curr => curr%next
  enddo

 end function potV_denx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 real(Long) function potV_x(this,xx) result(pot)
  implicit none

  class(PotentialList), intent(in) :: this
  real (Long), intent(in) :: xx

  class(Link_Potential), pointer :: curr

  pot=0._Long

  curr => this%firstLink

  do while(associated(curr))
   pot=pot+curr%myPot%potV_x(xx)
   curr => curr%next
  enddo

 end function potV_x

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

pure function PotentialList_isEmpty(this) result(isEmpty)
    logical isEmpty
    class(PotentialList), intent(in) :: this

    if (associated(this%firstLink)) then
       isEmpty = .false.
    else
       isEmpty = .true.
    endif
  end function PotentialList_isEmpty

end module class_PotentialList

