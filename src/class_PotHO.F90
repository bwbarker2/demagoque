module class_PotHO
 use class_Potential
 use prec_def
 implicit none

 private

 public :: PotHO, new_PotHO

 type, extends(Potential) :: PotHO
  private
  real(Long) :: mass  !< mass of particle in oscillator
  real(Long) :: freq !< frequency (strength) of oscillator field
 contains
  procedure :: potV_x => potV_x
 end type PotHO

contains

 function new_PotHO (newMass,newFreq)
  class(PotHO), pointer :: new_PotHO

  real(Long), intent(in) :: newMass
  real(Long), intent(in) :: newFreq

  allocate(new_PotHO)

  new_PotHO%mass=newMass
  new_PotHO%freq=newFreq

!write(*,*)'new_PotHO:',new_PotHO%mass,new_PotHO%freq

 end function new_PotHO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 real(Long) function potV_x(this,xx)
  class(PotHO), intent(in) :: this
  real (Long), intent(in) :: xx

  potV_x=0.5_Long*this%mass*(this%freq*xx)**2

!write(*,*)'potV_x:',this%mass,this%freq

 end function


end module class_PotHO

