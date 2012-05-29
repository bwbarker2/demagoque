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
  real (Long) :: alpha,beta  !< phase factors in exponents
  
 contains
  procedure,public :: getWavefn => stateKronigPenney_getWavefn
  procedure,private :: energyRoot => stateKronigPenney_energyRoot
  procedure,private :: dRdE => stateKronigPenney_dRdE

 end type StateKronigPenney

 !> Need this to call differentiation function which takes a function with only
 !! one argument.
 type(StateKronigPenney) :: currState

contains
 
 !> Constructor of StateKronigPenney object
 function new_StateKronigPenney(mass,level,v0,bwidth,period) result(this)
  type (StateKronigPenney) :: this

  real (Long) ,intent(in) :: mass
  integer     ,intent(in) :: level
  real (Long) ,intent(in) :: v0
  real (Long) ,intent(in) :: bwidth
  real (Long) ,intent(in) :: period

  write(*,*)'Entering function new_StateKronigPenney'

  this=stateKronigPenney(mass,level,v0,bwidth,period &
       ,0._Long,0._Long,0._Long,0._Long,0._Long,0._Long,0._Long)

  currState=this

  call stateKronigPenney_calcEnergy(this)

  write(*,*)'Leaving function new_StateKronigPenney'

 end function new_StateKronigPenney

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !> Calculates the energy of this level, specified by this%level.
 !!
 !! This involves finding the zeroes of the equation given by the energyRoot
 !! function. Since the zeroes occur in pairs, very close to each other, it is
 !! difficult to make initial guesses for these zeroes. Since they are on either
 !! side of an extremum, the initial guesses can be taken as the zeroes of the 1st
 !! derivative, which are much easier to guess. This is what is done here. dEdR is
 !! the derivative of energyRoot.
 subroutine stateKronigPenney_calcEnergy(this)
  use iso_fortran_env
  implicit none

  class(StateKronigPenney), intent(inout) :: this

  real (Long) :: hinitial !< initial guess for dRdE

  integer :: levelsofar !< energy level found so far
  real (Long) :: elast,ecurr !< energy of level found so far
  real (Long) :: erootlast,erootcurr !< value of level

!  write(*,*)'Entering subroutine stateKronigPenney_calcEnergy'

  levelsofar=-1

  ! 1/10 of a coefficient of eek in the energyRoot equation
  hinitial = 0.1_Long/(sqrt(2._Long*this%mass)/hbar*(this%period-this%bwidth))

  elast=2._Long*hinitial
  erootlast=this%dRdE(elast)

  !work up the energy levels until we find the one we want
  do while(levelsofar<this%level)
   ecurr=elast+hinitial
   erootcurr=this%dRdE(ecurr)

   !find initial guesses for dRdE zero
   do while(erootlast*erootcurr>=0._Long)
    ecurr=ecurr+hinitial
    erootcurr=this%dRdE(ecurr)
   enddo

   !using initial guesses, find zero of dRdE (thus extremum of energyRoot)
   ecurr=bmath_dZeroBrent(ecurr-hinitial,ecurr,getCurrDRDE)

   !find zero of energyRoot, using elast and extremum of energyRoot found above
   this%energy=bmath_dZeroBrent(elast,ecurr,getCurrEnergyRoot)

   levelsofar=levelsofar+1

   !now we know a better scale for the initial guesses for zeroes of drde
   hinitial=0.333*(ecurr-elast)

   !set up for next energy
   elast=ecurr+epzero*2._Long**10  ! start just after the last extremum
                                   ! (otherwise, it finds this root again and
                                   !  gives a root-not-bracketed error)
   erootlast=this%DRDE(elast)

!   write(ERROR_UNIT,*)'level,eek,energyRoot(eek)=' &
!                     ,levelsofar,this%energy,this%energyRoot(this%energy)

 enddo !while levelsofar


!  write(*,*)'Leaving subroutine stateKronigPenney_calcEnergy'

 end subroutine stateKronigPenney_calcEnergy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 pure real(Long) function stateKronigPenney_getWaveFn(this,xx) result(wf)
  implicit none

  class(StateKronigPenney), intent(in) :: this
  real (Long)            , intent(in) :: xx

  wf=1._Long
 end function stateKronigPenney_getWaveFn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 real (Long) function stateKronigPenney_dRdE(this,eek) result(dRdE)
  implicit none

  class(stateKronigPenney), intent(in) :: this
  real (Long), intent(in) :: eek !< sqrt of energy of state

  real (Long) :: hinitial  ! initial step for derivative - length scale

  real (Long) :: error  !output tolerance of derivative

!  write(*,*)'Entering function stateKronigPenney_dRdE'

  ! 1/10 of a coefficient of eek in the energyRoot equation
  hinitial = 0.1_Long/(sqrt(2._Long*this%mass)/hbar*(this%period-this%bwidth))

  select type(this)
   type is (stateKronigPenney)
    currState=this
  end select
  dRdE=bmath_LDiffRichardson(getCurrEnergyRoot,eek,hinitial,error)

 end function stateKronigPenney_dRdE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 real (Long) function stateKronigPenney_energyRoot(this,eek) result(en)
  implicit none

  class(stateKronigPenney), intent(in) :: this
  real (Long), intent(in) :: eek  !< sqrt of energy of state

  real (Long) :: alpha  ! trial alpha and beta
  complex (Long) :: beta

  alpha=sqrt(2._Long*this%mass)*eek/hbar
  beta=sqrt(cmplx(2._Long*this%mass*(eek**2-this%v0),KIND=Long))/hbar

  en = real(cos(alpha*(this%period-this%bwidth))*cos(beta*this%bwidth) &
            - (alpha**2+beta**2)/(2._Long*alpha*beta) &
              * sin(alpha*(this%period-this%bwidth))*sin(beta*this%bwidth) &
            - 1._Long &
           )

 end function stateKronigPenney_energyRoot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 real (Long) function getCurrEnergyRoot(eek) result(root)
  implicit none

  real (Long), intent(in) :: eek

  root=currState%energyRoot(eek)

 end function getCurrEnergyRoot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 real (Long) function getCurrDRDE(eek) result(drde)
  implicit none

  real (Long), intent(in) :: eek

  drde=currState%dRdE(eek)

 end function getCurrDRDE

end module class_StateKronigPenney

