module class_SuperWavefunction
 use class_OneWavefunction
 use class_Wavefunction
 use class_WfKronigPenney
 use prec_def
 implicit none

 private

 public :: SuperWavefunction, new_SuperWavefunction

 integer, parameter :: INITIAL_CAPACITY = 1

 type, extends(Wavefunction) :: SuperWavefunction
  type(OneWavefunction), allocatable, dimension(:) :: wavefunctions
!  type(WfKronigPenney), allocatable, dimension(:) :: WfKronigPenneys
  integer :: numWf  !< number of wf's
 contains
!  procedure,public :: add => SuperWavefunction_add
  procedure,public :: getWavefn => SuperWavefunction_getWavefn
 end type SuperWavefunction

contains

 function new_SuperWavefunction() result(new)
  implicit none

  type(SuperWavefunction) :: new

  type(OneWavefunction),allocatable, dimension(:) :: wfs

  allocate(wfs(INITIAL_CAPACITY))

  new = SuperWavefunction(wfs,0)

 end function new_SuperWavefunction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine SuperWavefunction_add(this,wftoadd)
  implicit none

  class(SuperWavefunction), intent(inout) :: this
  class(Wavefunction), intent(in) :: wftoadd

  if(size(this%wavefunctions)<=this%numWf) then
   call SuperWavefunction_grow(this%wavefunctions)
  endif
  this%numWf=this%numWf+1
  this%wavefunctions(this%numWf)=new_OneWavefunction(wftoadd)

 end subroutine SuperWavefunction_add

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 complex(Long) function SuperWavefunction_getWavefn(this,xx,tt) result(wf)
  use phys_cons
  implicit none

  class(SuperWavefunction), intent(in) :: this
  real (Long), intent(in) :: xx
  real (Long), optional, intent(in) :: tt

  integer :: ii

  wf=czero

  do ii=1,this%numWf
   wf=wf+this%wavefunctions(ii)%getWavefn(xx,tt)
  enddo

 end function SuperWavefunction_getWavefn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine SuperWavefunction_grow(wfarray)
  type(OneWavefunction), allocatable, dimension(:), intent(inout) :: wfarray

  type(OneWavefunction), dimension(size(wfarray)) :: wftemp
  integer :: newSize


  wftemp=wfarray

  deallocate(wfarray)

  if(size(wftemp)>1) then
   newSize=nint(size(wftemp)*1.5)
  else if(size(wftemp)==1) then
   newSize=2
  else
   newSize=1
  endif

  allocate(wfarray(newSize))

  wfarray(1:size(wftemp))=wftemp

 end subroutine SuperWavefunction_grow


end module class_SuperWavefunction

