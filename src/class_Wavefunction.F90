!> parent class of wavefunctions, so they can be put in an array together.
module class_Wavefunction
 use prec_def
 implicit none

 private

 public :: Wavefunction, Skin_Wavefunction, new_Skin_Wavefunction

 type, abstract :: Wavefunction 
  contains
   procedure(getWavefn), deferred :: getWavefn

 end type Wavefunction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !> Wrapper (or 'skin') for class Wavefunction, so we can have an array of
 !! different children of Wavefunction
 type, extends(Wavefunction) :: Skin_Wavefunction
  class(Wavefunction), pointer :: oneWf
 contains
  procedure :: getWavefn => Skin_Wavefunction_getWavefn

 end type Skin_Wavefunction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 interface

  complex(Long) function getWavefn(this,xx,tt)
   use prec_def
   import Wavefunction
 
   class(Wavefunction),  intent(in) :: this
   real(Long),           intent(in) :: xx
   real(Long), optional, intent(in) :: tt

  end function getWavefn

 end interface

contains

 !> Constructs new Wavefunction skin out of a Wavefunction
 function new_Skin_Wavefunction(wavefin) result(this)
  implicit none

  type(Skin_Wavefunction),pointer :: this

  class(Wavefunction), target, intent(in) :: wavefin

!  class(Wavefunction), allocatable :: wf

!  if(allocated(wf)) then
!   deallocate(wf)
!  endif
  allocate(this)
  this%oneWf => wavefin
!  allocate(wf,source=wavefin)
!  this = Skin_Wavefunction(wavefin)
!  deallocate(this%oneWf)
!  allocate(this%oneWf, source=wavefin)

 end function new_Skin_Wavefunction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 complex(Long) function Skin_Wavefunction_getWavefn(this,xx,tt) result(wf)
  implicit none

  class(Skin_Wavefunction), intent(in) :: this
  real (Long), intent(in) :: xx
  real (Long), optional, intent(in) :: tt

  wf=this%oneWf%getWavefn(xx,tt)

 end function Skin_Wavefunction_getWavefn

end module class_Wavefunction

