!> parent class of wavefunctions, so they can be put in an array together.
module class_Wavefunction
 use prec_def
 implicit none

 private

 public :: Wavefunction

 type, abstract :: Wavefunction 
  contains
   procedure(getWavefn), deferred :: getWavefn

 end type Wavefunction

 interface

  complex(Long) function getWavefn(this,xx,tt)
   use prec_def
   import Wavefunction
 
   class(Wavefunction),  intent(in) :: this
   real(Long),           intent(in) :: xx
   real(Long), optional, intent(in) :: tt

  end function getWavefn

 end interface

end module class_Wavefunction

