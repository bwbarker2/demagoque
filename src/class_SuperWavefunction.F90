module class_SuperWavefunction
 use class_Wavefunction
 use prec_def
 implicit none

 private

 public :: SuperWavefunction, new_SuperWavefunction

 integer, parameter :: INITIAL_CAPACITY = 10

 type, extends(Wavefunction) :: SuperWavefunction
  type(Skin_Wavefunction), allocatable, dimension(:) :: wavefunctions
  integer :: numWf  !< number of wf's
 contains
  procedure,public :: add => SuperWavefunction_add
  procedure,public :: getWavefn => SuperWavefunction_getWavefn
  procedure :: grow => SuperWavefunction_grow
 end type SuperWavefunction

contains

 function new_SuperWavefunction() result(new)
  implicit none

  type(SuperWavefunction) :: new

  type(Skin_Wavefunction),allocatable, dimension(:) :: wfs

  allocate(wfs(INITIAL_CAPACITY))

  new = SuperWavefunction(wfs,0)

 end function new_SuperWavefunction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine SuperWavefunction_add(this,wftoadd)
  implicit none

  class(SuperWavefunction), intent(inout) :: this
  class(Wavefunction), intent(in) :: wftoadd

  if(size(this%wavefunctions)<=this%numWf) then
   call this%grow()
  endif
  this%numWf=this%numWf+1
  this%wavefunctions(this%numWf)=new_Skin_Wavefunction(wftoadd)

!  write(*,*)'SuperWavefunction_add: wftoaddwf=',this%wavefunctions(this%numWf)%getWavefn(-30)
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
   write(ERROR_UNIT,*)'SupWav_getWav=',wf
  enddo

 end function SuperWavefunction_getWavefn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !> Doesn't work yet :(
 subroutine SuperWavefunction_grow(this)
  class(SuperWavefunction), intent(inout) :: this
!  type(Skin_Wavefunction), allocatable, dimension(:), intent(inout) :: wfarray

  type(Skin_Wavefunction), dimension(size(this%wavefunctions)) :: wftemp
  integer :: newSize


  wftemp=this%wavefunctions

  if(size(wftemp)>1) then
   newSize=nint(size(wftemp)*1.5)
  else if(size(wftemp)==1) then
   newSize=2
  else
   newSize=1
  endif

write(ERROR_UNIT,*)'SuperWavefunction_grow: dim=',size(wftemp)
  deallocate(this%wavefunctions)
!deallocate(this%wavefunctions(1)%oneWf)

write(ERROR_UNIT,*)'SuperWavefunction_grow: deallocated wfarray'
  allocate(this%wavefunctions(newSize))

  this%wavefunctions(1:size(wftemp))=wftemp

 end subroutine SuperWavefunction_grow


end module class_SuperWavefunction

