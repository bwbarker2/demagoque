module class_StateKronigPenney
 use bmath
 use phys_cons
 use prec_def
 implicit none

 private

 public :: StateKronigPenney, new_StateKronigPenney

 type StateKronigPenney
  private
  real(Long)  :: mass  !< mass of particle
  integer     :: level !< excitation level, 0=ground state
  real (Long) :: v0    !< height of rectangular barrier
  real (Long) :: bwidth    !< width of barrier
  real (Long) :: period    !< length of period
  real (Long) :: energy !< energy of state
  real (Long) :: a1,a2,b1,b2 !< coefficients of wavefunctions
  
 contains
  procedure,public :: getWavefn => stateKronigPenney_getWavefn
  procedure,private :: energyRoot => stateKronigPenney_energyRoot

 end type StateKronigPenney

contains
 
 pure function new_StateKronigPenney(mass,level,v0,bwidth,period) result(this)
  type (StateKronigPenney) :: this

  real (Long) ,intent(in) :: mass
  integer     ,intent(in) :: level
  real (Long) ,intent(in) :: v0
  real (Long) ,intent(in) :: bwidth
  real (Long) ,intent(in) :: period

  this=stateKronigPenney(mass,level,v0,bwidth,period &
                         ,0e0_Long,0e0_Long,0e0_Long,0e0_Long,0e0_Long)

 end function new_StateKronigPenney

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 pure real(Long) function stateKronigPenney_getWaveFn(this,xx) result(wf)
  implicit none

  class(StateKronigPenney), intent(in) :: this
  real (Long)            , intent(in) :: xx

  wf=0e0_Long
 end function stateKronigPenney_getWaveFn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 pure real (Long) function stateKronigPenney_dRdE(this,eek) result(dRdE)
  implicit none

  class(stateKronigPenney), intent(in) :: this
  real (Long), intent(in) :: eek !< sqrt of energy of state

  

 end function stateKronigPenney_dRdE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 pure real (Long) function stateKronigPenney_energyRoot(this,eek) result(en)
  implicit none

  class(stateKronigPenney), intent(in) :: this
!  real(Long), intent(in) :: mass & !< mass of particle
!                           ,aa & !< width of periodic box
!                           ,bb & !< width of barrier
!                           ,v0 & !< height of barrier
   real (Long), intent(in) :: eek  !< sqrt of energy of state

  real (Long) :: alpha,beta

  alpha=sqrt(2e0_Long*this%mass)*eek/hbar
  beta=sqrt(2e0_Long*this%mass*(eek**2-this%v0))/hbar

  en =  cos(alpha*(this%period-this%bwidth))*cos(beta*this%bwidth) &
      - (alpha**2+beta**2)/(2e0_Long*alpha*beta) &
        * sin(alpha*(this%period-this%bwidth))*sin(beta*this%bwidth) &
      - 1e0_Long

 end function stateKronigPenney_energyRoot

end module class_StateKronigPenney

