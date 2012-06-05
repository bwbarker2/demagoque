module class_OneWavefunction
 use class_Wavefunction
 use prec_def
 implicit none

 private

 public :: OneWavefunction, new_OneWavefunction

 type, extends(Wavefunction) :: OneWavefunction
  class(Wavefunction), allocatable :: oneWf
 contains
  procedure :: getWavefn => OneWavefunction_getWavefn

 end type OneWavefunction

contains

 function new_OneWavefunction(wavefin) result(this)
  implicit none

  type(OneWaveFunction) :: this

  class(Wavefunction), intent(in) :: wavefin

  this = OneWavefunction(wavefin)

 end function new_OneWaveFunction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 complex(Long) function OneWavefunction_getWavefn(this,xx,tt) result(wf)
  implicit none

  class(OneWavefunction), intent(in) :: this
  real (Long), intent(in) :: xx
  real (Long), optional, intent(in) :: tt

  wf=this%oneWf%getWavefn(xx,tt)

 end function OneWavefunction_getWavefn

end module class_OneWavefunction

