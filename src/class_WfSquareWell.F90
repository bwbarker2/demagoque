!> calculates properties of the square well potential quantum system
module class_WfSquareWell
 use bexception
 use bmath
 use class_Wavefunction
 use iso_fortran_env
 use phys_cons
 use prec_def
 implicit none

  private :: currentState, WfSquareWell_calcEnergy, &
             WfSquareWell_calcNorm, WfSquareWell_energyRoot

  public
  type, extends(Wavefunction) :: WfSquareWell
   real (Long) :: mass     !< mass of particle
   integer     :: level    !< excitation level, 0=ground state, 1=first excited
   real (Long) :: v0       !< height of barrier
   real (Long) :: d        !< half width of well
   real (Long) :: energy   !< energy of state
   real (Long) :: norm     !< normalization constant
   real (Long) :: normTail !< normalization of wavefunction within barrier

   contains
    procedure :: getWavefn => WfSquareWell_getWavefn
    procedure :: getPotential => WfSquareWell_getPotential
  end type WfSquareWell

  type(WfSquareWell) :: currentState !< for root-finding function, it can only have one argument, so set this as current state before running root finder

 contains

  function make_WfSquareWell (mass,level,v0,d) result(self)
   implicit none

   type (WfSquareWell)  :: self
   real (Long), intent(in) :: mass
   integer, intent(in)     :: level
   real (Long), intent(in) :: v0
   real (Long), intent(in) :: d

   real (Long) :: energy, norm, normTail

   ! initial values
   energy=0e0_Long
   norm=0e0_Long
   normTail=0e0_Long

   self=WfSquareWell(mass,level,v0,d,energy,norm,normTail)

!   write(*,*)self
!   write(*,*)'hbar=',hbar

   self%energy=WfSquareWell_calcEnergy(self)
   call WfSquareWell_calcNorm(self)

  end function make_WfSquareWell

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
  function WfSquareWell_calcEnergy(state) result(energy)
   implicit none

   type(WfSquareWell), intent(in) :: state
   real (Long) :: energy

   real (Long) :: lowbound,upbound !< bounds for zero-finding

   lowbound = state%level*pi
   upbound = lowbound+0.5_Long*pi

   currentState=state

!   write(ERROR_UNIT,*)WfSquareWell_energyRoot(lowbound)

   energy=bmath_dZeroBrent(lowbound,upbound,WfSquareWell_energyRoot)
!   write(ERROR_UNIT,*)'level,root:',state%level,energy
   energy=energy**2*hbar**2/(2e0_Long*state%mass*state%d**2)
!   write(ERROR_UNIT,*)'energy:',energy

  end function WfSquareWell_calcEnergy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculates normalization constant
  subroutine WfSquareWell_calcNorm(state)
   implicit none

   type(WfSquareWell), intent(inout) :: state

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

  end subroutine WfSquareWell_calcNorm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculates equation to find zero of, in order to find energy level
  function WfSquareWell_energyRoot(uu) result(fn)
   implicit none

   real(Long) :: fn
   real(Long) :: uu !< dimensionless energy parameter

   fn=2e0_Long*currentState%mass*currentState%v0*currentState%d**2 &
      /(hbar**2*(1e0_Long+tan(uu)**2))-uu**2

  end function WfSquareWell_energyRoot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns potential of state
  function WfSquareWell_getPotential(self,xx) result(pot)
   implicit none

   class (WfSquareWell), intent(in) :: self !< state
   real (Long) :: xx  !< position at which to get potential
   real (Long) :: pot

   if(abs(xx)<=self%d) then
    pot=0e0_Long
   else
    pot=self%v0
   endif

  end function WfSquareWell_getPotential

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function WfSquareWell_getWavefn(this,xx,tt) result(wf)
   implicit none

   complex (Long) :: wf
   class (WfSquareWell), intent(in) :: this
   real (Long), intent(in) :: xx
   real (Long), optional, intent(in) :: tt

   if(xx<-this%d) then
    wf=this%normTail*exp(sqrt(2e0_Long*this%mass*(this%v0-this%energy))/hbar*xx)
   elseif(xx>this%d) then
    wf=this%normTail*exp(-sqrt(2e0_Long*this%mass*(this%v0-this%energy))/hbar*xx)
   else
    wf=this%norm*cos(sqrt(2e0_Long*this%mass*this%energy)/hbar*xx)
   endif

  end function WfSquareWell_getWavefn

end module class_WfSquareWell

