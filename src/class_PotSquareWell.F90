module class_PotSquareWell
 use class_Potential
 use prec_def
 implicit none

 private

 public :: PotSquareWell, new_PotSquareWell

 type, extends(Potential) :: PotSquareWell
  private
  real(Long) :: hwidth  !< half width of well (-hwidth to +hwidth)
  real(Long) :: v0 !< height of barrier
 contains
  procedure :: potV_x => potV_x
 end type PotSquareWell

contains

 function new_PotSquareWell (hwidthn,v0n)
  class(PotSquareWell), pointer :: new_PotSquareWell

  real(Long), intent(in) :: hwidthn
  real(Long), intent(in) :: v0n

  allocate(new_PotSquareWell)

  new_PotSquareWell%hwidth=hwidthn
  new_PotSquareWell%v0=v0n

 end function new_PotSquareWell

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 real(Long) function potV_x(this,xx)
  class(PotSquareWell), intent(in) :: this
  real (Long), intent(in) :: xx

  if(abs(xx)>=this%hwidth) then
   potV_x=this%v0
  else
   potV_x=0._Long
  endif

!write(*,*)'PotHO: potV_x = ',xx,potV_x

 end function

end module class_PotSquareWell

