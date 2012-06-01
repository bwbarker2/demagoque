module class_StateKronigPenney
 use bmath
 use phys_cons
 use prec_def
 implicit none

 private

 public :: StateKronigPenney, new_StateKronigPenney

 type StateKronigPenney
  private
  real(Long)     :: mass        !< mass of particle
  integer        :: level       !< excitation level, 0=ground state
  real (Long)    :: v0          !< height of rectangular barrier
  real (Long)    :: bwidth      !< width of barrier
  real (Long)    :: period      !< length of period
  real (Long)    :: energy      !< energy of state
  complex (Long) :: a1,a2,b1,b2 !< coefficients of wavefunctions
  real (Long)    :: alpha       !< phase factor in exponent
  complex (Long) :: beta        !< another phase factor
  
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

!  write(*,*)'Entering function new_StateKronigPenney'

  this=stateKronigPenney(mass,level,v0,bwidth,period &
       ,0._Long,czero,czero,czero,czero,0._Long,czero)

  currState=this

  call stateKronigPenney_calcEnergy(this)

  call stateKronigPenney_calcCoeffs(this)

  write(*,*)'this=',this

!  write(*,*)'Leaving function new_StateKronigPenney'

 end function new_StateKronigPenney

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine stateKronigPenney_calcCoeffs(this)
  implicit none

  class(StateKronigPenney), intent(inout) :: this

  complex(Long), dimension(4,4) :: coeffeqn

  ! needed for zcgesv below
  integer, dimension(4) :: ipivs ! pivot indices, row i was interchanged with
                                 ! row ipivs(i)
  complex(Long), dimension(4) :: rhs,sols,rworks  !RHS, solutions of linear equation to solve
  complex(Long), dimension(4) :: works
  complex, dimension(20) :: sworks  !single-precision work array
  integer :: iters !result of iterative refinement
  integer :: infos !result of solver

  complex(Long) :: norm !normalization factor

!  integer :: ii,jj
 

  coeffeqn(1,1) =  1._Long
  coeffeqn(2,1) =  1._Long
  coeffeqn(3,1) = -1._Long
  coeffeqn(4,1) = -1._Long
  coeffeqn(1,2) =  this%alpha
  coeffeqn(2,2) = -this%alpha
  coeffeqn(3,2) = -this%beta
  coeffeqn(4,2) =  this%beta
  coeffeqn(1,3) =             exp( imagi*this%alpha*(this%period-this%bwidth))
  coeffeqn(2,3) =             exp(-imagi*this%alpha*(this%period-this%bwidth))
  coeffeqn(3,3) =            -exp(-imagi*this%beta *(this%bwidth)            )
  coeffeqn(4,3) =            -exp( imagi*this%beta *(this%bwidth)            )
  coeffeqn(1,4) =  this%alpha*exp( imagi*this%alpha*(this%period-this%bwidth))
  coeffeqn(2,4) = -this%alpha*exp(-imagi*this%alpha*(this%period-this%bwidth))
  coeffeqn(3,4) = -this%beta *exp(-imagi*this%beta *(this%bwidth)            )
  coeffeqn(4,4) =  this%beta *exp( imagi*this%beta *(this%bwidth)            )

  sols=czero

!evidently zcgesv (from LAPACK) wants a row-major array? Who knew?
coeffeqn=transpose(coeffeqn)

  call zcgesv(4,1,coeffeqn,4,ipivs,rhs,4,sols,4,works,sworks,rworks,iters,infos)

!  write(*,*)'sols=',sols

  this%a1=sols(1)
  this%a2=sols(2)
  this%b1=sols(3)
  this%b2=sols(4)

  !normalize to 1 over the box length, integral from -bwidth to period-bwidth
  !this is the completed integral, from notes BWB 2012-05-08p2

  norm =   (this%a1*conjg(this%a1)+this%a2*conjg(this%a2)) &
           * (this%period-this%bwidth) &
     + 2._Long*real(this%a1*conjg(this%a2)/(2._Long*imagi*this%alpha) &
                    * exp(2._Long*imagi*this%alpha*(this%period-this%bwidth))) &
     + (this%b1*conjg(this%b1)+this%b2*conjg(this%b2))*this%bwidth &
     + 2._Long*real(conjg(this%b1)*this%b2/(2._Long*imagi*this%beta) &
                    * exp(2._Long*imagi*this%beta*this%bwidth))

  norm=1._Long/sqrt(norm)

!write(*,*)'norm=',norm

  this%a1=this%a1*norm
  this%a2=this%a2*norm
  this%b1=this%b1*norm
  this%b2=this%b2*norm

!  do ii=1,4
!   write(*,'(4(2ES12.3,4X))')(real(coeffeqn(jj,ii)),aimag(coeffeqn(jj,ii)),jj=1,4)
!  enddo

!  write(*,*)'infos=',infos
!  write(*,*)'iters=',iters
!  write(*,*)'sols=',sols
!  write(*,*)'ipivs=',ipivs

 end subroutine stateKronigPenney_calcCoeffs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

 ! we were solving for sqrt(ee). Square to get the actual energy of the state
 this%energy=this%energy**2

 this%alpha=sqrt(2._Long*this%mass*this%energy)/hbar
 this%beta=sqrt(cmplx(2._Long*this%mass*(this%energy-this%v0),KIND=Long))/hbar

!  write(*,*)'Leaving subroutine stateKronigPenney_calcEnergy'

 end subroutine stateKronigPenney_calcEnergy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !> Returns the wavefunction at position xx.
 !!
 !! We did a coordinate shift in the calculation of the wavefunction, so,
 !! assuming that the periodic barriers are symmetric about the origin, and
 !! that at the origin there is no barrier, then this does the proper
 !! coordinate shift.
 complex(Long) function stateKronigPenney_getWaveFn(this,xxin) result(wf)
  implicit none

  class(StateKronigPenney), intent(in) :: this
  real (Long)            , intent(in) :: xxin

  real (Long) :: xx !position in the coordinate frame of the calculation
  complex(Long) :: wfhere

!  select type(this)
!   type is (stateKronigPenney)
!    write(*,*)this
!  end select

  xx=xxin+(this%period-this%bwidth)*0.5_Long

  do while ( xx < (-this%bwidth) )
   xx=xx+this%period
  enddo

  do while(xx>=(this%period-this%bwidth))
   xx=xx-this%period
  enddo

  if (xx >= 0._Long ) then
   wf= this%a1*exp( imagi*this%alpha*xx) &
      +this%a2*exp(-imagi*this%alpha*xx)
  else
   wf= this%b1*exp( imagi*this%beta*xx) &
      +this%b2*exp(-imagi*this%beta*xx)
  endif

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

