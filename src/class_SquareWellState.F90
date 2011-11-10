!> calculates properties of the square well potential quantum system
module class_SquareWellState
 use bexception
 use bmath
 use iso_fortran_env
 use phys_cons
 use prec_def
 implicit none

  private :: currentState, squareWellState_calcEnergy, &
             squareWellState_calcNorm, squareWellState_energyRoot

  public
  type squareWellState
   real (Long) :: mass     !< mass of particle
   integer     :: level    !< excitation level, 0=ground state, 1=first excited
   real (Long) :: v0       !< height of barrier
   real (Long) :: d        !< half width of well
   real (Long) :: energy   !< energy of state
   real (Long) :: norm     !< normalization constant
   real (Long) :: normTail !< normalization of wavefunction within barrier

   contains
    procedure :: getWavefn => squareWellState_getWavefn
    procedure :: getPotential => squareWellState_getPotential
  end type squareWellState

  type(squareWellState) :: currentState !< for root-finding function, it can only have one argument, so set this as current state before running root finder

 contains

  function make_squareWellState (mass,level,v0,d) result(self)
   implicit none

   type (squareWellState)  :: self
   real (Long), intent(in) :: mass
   integer, intent(in)     :: level
   real (Long), intent(in) :: v0
   real (Long), intent(in) :: d

   real (Long) :: energy, norm, normTail

   ! initial values
   energy=0
   norm=0
   normTail=0

   self=squareWellState(mass,level,v0,d,energy,norm,normTail)

   self%energy=squareWellState_calcEnergy(self)
   call squareWellState_calcNorm(self)

  end function make_squareWellState

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculates energy of this level by matching first derivations of the wf
  !! at the barrier.
  !!
  !! This means numerically solving the equality:
  !! \f[
  !!  U^2=\frac{2 m V_0 d^2 / \hbar^2}{1+\tan^2 U}, \quad U=\frac{2 m E}{\hbar} d
  !! \f]
  !!
  !! For more information, see BWB Notes 2011-11-02p1 and 2011-11-08p1.
  !! 
  !! Correction: \f$ \tan U \f$ must be positive to satisfy the original
  !! equality, so only look for solutions in 1st and 3rd quadrants of unit
  !! circle.
  function squareWellState_calcEnergy(state) result(energy)
   implicit none

   type(squareWellState), intent(in) :: state
   real (Long) :: energy

   real (Long) :: lowbound,upbound !< bounds for zero-finding

   lowbound = state%level*pi
   upbound = lowbound+0.5_Long*pi

   currentState=state

   energy=bmath_dZeroBrent(lowbound,upbound,squareWellState_energyRoot)
!   write(ERROR_UNIT,*)'level,root:',state%level,energy
   energy=energy**2*hbar**2/(2e0_Long*state%mass*state%d**2)
!   write(ERROR_UNIT,*)'energy:',energy

  end function squareWellState_calcEnergy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculates normalization constant
  subroutine squareWellState_calcNorm(state)
   implicit none

   type(squareWellState), intent(inout) :: state

   real (Long) :: sqrt2mEhbar

   sqrt2mEhbar=sqrt(2e0_Long*state%mass*state%energy)/hbar

   state%norm=( hbar*cos(sqrt2mEhbar*state%d)**2 &
                /sqrt(2e0_Long*state%mass*(state%v0-state%energy)) &
               +state%d &
               +0.5_Long/sqrt2mEhbar*sin(2e0_Long*sqrt2mEhbar*state%d) &
              )**(-0.5_Long)

   state%normTail=state%norm &
                  *cos(sqrt2mEhbar*state%d) &
                  *exp(sqrt(2e0_Long*state%mass*(state%v0-state%energy)) &
                       /hbar*state%d &
                      )

  end subroutine squareWellState_calcNorm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculates equation to find zero of, in order to find energy level
  function squareWellState_energyRoot(uu) result(fn)
   implicit none

   real(Long) :: fn
   real(Long) :: uu !< dimensionless energy parameter

   fn=2e0_Long*currentState%mass*currentState%v0*currentState%d**2 &
      /(hbar**2*(1e0_Long+tan(uu)**2))-uu**2

  end function squareWellState_energyRoot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns potential of state
  function squareWellState_getPotential(self,xx) result(pot)
   implicit none

   class (squareWellState), intent(in) :: self !< state
   real (Long) :: xx  !< position at which to get potential
   real (Long) :: pot

   if(abs(xx)<=self%d) then
    pot=self%v0
   else
    pot=0e0_Long
   endif

  end function squareWellState_getPotential

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function squareWellState_getWavefn(self,xx) result(wf)
   implicit none

   real (Long) :: wf
   class (squareWellState), intent(inout) :: self
   real (Long), intent(in) :: xx

   if(xx<-self%d) then
    wf=self%normTail*exp(sqrt(2e0_Long*self%mass*(self%v0-self%energy))/hbar*xx)
   elseif(xx>self%d) then
    wf=self%normTail*exp(-sqrt(2e0_Long*self%mass*(self%v0-self%energy))/hbar*xx)
   else
    wf=self%norm*cos(sqrt(2e0_Long*self%mass*self%energy)/hbar*xx)
   endif

  end function squareWellState_getWavefn

end module class_SquareWellState

