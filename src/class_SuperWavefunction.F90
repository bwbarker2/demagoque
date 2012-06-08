!> Stores a superposition of wavefunctions of type Wavefunction and is 
!! itself a Wavefunction
module class_SuperWavefunction
 use class_Wavefunction
 use prec_def
 implicit none

 private

 public :: new_SuperWavefunction, SuperWavefunction

 integer, parameter :: INITIAL_CAPACITY = 10

 !> Superposition of wavefunctions
 type, extends(Wavefunction) :: SuperWavefunction
  private
  class(Link_Wavefunction), pointer :: firstLink => null()  !< first link in linked list
  class(Link_Wavefunction), pointer :: lastLink => null()   !< last link
 contains
  procedure,public :: add => SuperWavefunction_add  !< add a wavefunction to the list
  procedure,public :: getWavefn => SuperWavefunction_getWavefn !< get sum of wavefunctions
 end type SuperWavefunction

 !> Individual link in linked list
 type :: Link_Wavefunction
  private
  class(Wavefunction), pointer :: myWf => null()     !< ptr to wavefunction
  class(Link_Wavefunction), pointer :: next => null() !< ptr to next link
 contains
  private
  procedure :: getWf => SuperWavefunction_getWf      !< get wavefunction
  procedure :: nextLink => SuperWavefunction_nextLink  !< returns ptr to next link
  procedure :: setNextLink => SuperWavefunction_setNextLink !< sets 'next'
 end type Link_Wavefunction

contains

 function SuperWavefunction_getWf(this) result(getWf)
  class(Wavefunction), pointer :: getWf
  class(Link_Wavefunction), intent(in) :: this

  getWf => this%myWf
 end function SuperWavefunction_getWf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 function SuperWavefunction_nextLink(this) result(nextLink)
  class(Link_Wavefunction), pointer :: nextLink

  class(Link_Wavefunction), intent(in) :: this

  nextLink => this%next
 end function SuperWavefunction_nextLink

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 function new_Link_Wavefunction(newWf,next)
  class(Link_Wavefunction), pointer :: new_Link_Wavefunction

  class(Wavefunction), intent(in) :: newWf
  class(Link_Wavefunction), pointer, intent(in) :: next

  allocate(new_Link_Wavefunction)
  new_Link_Wavefunction%next => next
  allocate(new_Link_Wavefunction%myWf, source=newWf)

 end function new_Link_Wavefunction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine SuperWavefunction_setNextLink(this,next)

  class(Link_Wavefunction), intent(inout) :: this
  class(Link_Wavefunction), pointer, intent(in) :: next

  this%next => next
 end subroutine SuperWavefunction_setNextLink

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 pure function new_SuperWavefunction() result(new)
  implicit none

  class(SuperWavefunction), pointer :: new

  allocate(new)

 end function new_SuperWavefunction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine SuperWavefunction_add(this,wftoadd)
  implicit none

  class(SuperWavefunction), intent(inout) :: this
  class(Wavefunction), intent(in) :: wftoadd

  class(Link_Wavefunction), pointer :: newLink

!write(ERROR_UNIT,*)'Entering SuperWavefunction_add'

  if (.not. associated(this%firstLink)) then
   this%firstLink => new_Link_Wavefunction(wftoadd,this%firstLink)
   this%lastLink => this%firstLink
  else
   newLink => new_Link_Wavefunction(wftoadd, this%lastLink%nextLink())
   call this%lastLink%setNextLink(newLink)
   this%lastLink => newLink
  endif

!write(ERROR_UNIT,*)'Leaving SuperWavefunction_add'

 end subroutine SuperWavefunction_add

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 complex(Long) function SuperWavefunction_getWavefn(this,xx,tt) result(wf)
  use phys_cons
  implicit none

  class(SuperWavefunction), intent(in) :: this
  real (Long), intent(in) :: xx
  real (Long), optional, intent(in) :: tt

  class(Link_Wavefunction), pointer :: curr

  wf=czero

  curr => this%firstLink

  do while(associated(curr))
   wf=wf+curr%myWf%getWavefn(xx,tt)
   curr => curr%next
  enddo

 end function SuperWavefunction_getWavefn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

pure function SuperWavefunction_isEmpty(this) result(isEmpty)
    logical isEmpty
    class(SuperWavefunction), intent(in) :: this

    if (associated(this%firstLink)) then
       isEmpty = .false.
    else
       isEmpty = .true.
    endif
  end function SuperWavefunction_isEmpty

end module class_SuperWavefunction

