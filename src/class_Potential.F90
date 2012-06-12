module class_Potential
 use prec_def
 implicit none

 private

 public :: Potential

 type :: Potential
 contains
  procedure :: potV_x => potV_x
  procedure :: potV_denx => potV_denx
 end type Potential

contains

 !> returns potential as a function of the density
 real(Long) function potV_denx(this,rho)
  class(Potential), intent(in) :: this
  real(Long), intent(in) :: rho  !< probability density

  potV_denx=0._Long

 end function potV_denx

 !> returns potential as a function of position
 real(Long) function potV_x(this, xx)
  class(Potential), intent(in) :: this
  real(Long), intent(in) :: xx !< position

  potV_x=0._Long

 end function potV_x

end module class_Potential
