module class_WfKronigPenney
 use bmath
 use class_Wavefunction
 use phys_cons
 use prec_def
 implicit none

 private

 public :: WfKronigPenney, new_WfKronigPenney

 type, extends(Wavefunction) :: WfKronigPenney
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
  procedure,public :: getWavefn => WfKronigPenney_getWavefn
  procedure,private :: energyRoot => WfKronigPenney_energyRoot
  procedure,private :: dRdE => WfKronigPenney_dRdE

 end type WfKronigPenney

 !> Need this to call differentiation function which takes a function with only
 !! one argument.
 class(WfKronigPenney), pointer :: currState

contains
 
 !> Constructor of WfKronigPenney object
 function new_WfKronigPenney(mass,level,v0,bwidth,period) result(this)
  class (WfKronigPenney), pointer :: this

  real (Long) ,intent(in) :: mass
  integer     ,intent(in) :: level
  real (Long) ,intent(in) :: v0
  real (Long) ,intent(in) :: bwidth
  real (Long) ,intent(in) :: period

  write(*,*)'Entering function new_WfKronigPenney'

  allocate(this)

  this%mass=mass
  this%level=level
  this%v0=v0
  this%bwidth=bwidth
  this%period=period
!  this=WfKronigPenney(mass,level,v0,bwidth,period &
!       ,0._Long,czero,czero,czero,czero,0._Long,czero)

  currState=>this

  call WfKronigPenney_calcEnergy(this)

  call WfKronigPenney_calcCoeffs(this)

  write(ERROR_UNIT,*)'new_WfKronigPenney: mass=',this%mass &
  ,this%level, this%v0,this%bwidth, this%period, this%energy, this%a1 &
  ,this%a2, this%b1, this%b2, this%alpha, this%beta

  write(*,*)'Leaving function new_WfKronigPenney'

 end function new_WfKronigPenney

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine WfKronigPenney_calcCoeffs(this)
  implicit none

  class(WfKronigPenney), intent(inout) :: this

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

  integer :: ii,jj
 

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

!  !eliminate extra row by gaussian elimination
!  call zGauss(4,coeffeqn,ipivs)
!
!  !erase lower triangle that actually has weights in it
!  do ii=2,4
!   do jj=1,ii-1
!    coeffeqn(ii,jj)=czero
!   enddo
!  enddo
!
!  !set the last row (which is all zero right now) to 0 0 0 1 = 1, which will
!  ! set b2 to 1. Don't worry, normalization will happen soon.
!  coeffeqn(4,4)=cmplx(1._Long,0._Long)
!  rhs(4)=1._Long
!
!  write(*,*)'coeffeqn:'
!  do ii=1,4
!   write(*,'(8E14.5)')(coeffeqn(ii,jj),jj=1,4)
!  enddo
!  write(*,*)'ipivs=',ipivs

  call zcgesv(4,1,coeffeqn,4,ipivs,rhs,4,sols,4,works,sworks,rworks,iters,infos)

  !if I comment this out, then a1,a2,b1,b2 are all NaN... I don't know why :(
  write(*,*)'sols=',sols

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

 end subroutine WfKronigPenney_calcCoeffs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !> Calculates the energy of this level, specified by this%level.
 !!
 !! This involves finding the zeroes of the equation given by the energyRoot
 !! function. Since the zeroes occur in pairs, very close to each other, it is
 !! difficult to make initial guesses for these zeroes. Since they are on either
 !! side of an extremum, the initial guesses can be taken as the zeroes of the 1st
 !! derivative, which are much easier to guess. This is what is done here. dEdR is
 !! the derivative of energyRoot.
 subroutine WfKronigPenney_calcEnergy(this)
  use iso_fortran_env
  implicit none

  class(WfKronigPenney), intent(inout) :: this

  real (Long) :: hinitial !< initial guess for dRdE

  integer :: levelsofar !< energy level found so far
  real (Long) :: elast,ecurr !< energy of level found so far
  real (Long) :: erootlast,erootcurr !< value of level

!  write(*,*)'Entering subroutine WfKronigPenney_calcEnergy'

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

!  write(*,*)'Leaving subroutine WfKronigPenney_calcEnergy'

 end subroutine WfKronigPenney_calcEnergy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !> Returns the wavefunction at position xx.
 !!
 !! We did a coordinate shift in the calculation of the wavefunction, so,
 !! assuming that the periodic barriers are symmetric about the origin, and
 !! that at the origin there is no barrier, then this does the proper
 !! coordinate shift.
 complex(Long) function WfKronigPenney_getWaveFn(this,xx,tt) result(wf)
  implicit none

  class(WfKronigPenney), intent(in) :: this
  real (Long)            , intent(in) :: xx
  real (Long), optional, intent(in) :: tt   !< time to evolve forward

  real (Long) :: xc !position in the coordinate frame of the calculation
  complex(Long) :: tfac !< factor for time evolution

!  select type(this)
!   type is (WfKronigPenney)
!    write(*,*)this
!  end select
!write(ERROR_UNIT,*)'WfKronigPenney_getWavefn: this=',this%mass,this%level,this%v0 &
!          ,this%bwidth,this%period, this%energy,this%a1,this%a2,this%b1 &
!          ,this%b2,this%alpha,this%beta

  if(present(tt)) then
   tfac= -imagi*tt/hbar
  else
   tfac=czero
  endif

  xc=xx+(this%period-this%bwidth)*0.5_Long

  do while ( xc < (-this%bwidth) )
!   write(*,*)'xc=',xc,this%period
   xc=xc+this%period
  enddo

  do while(xc>=(this%period-this%bwidth))
   xc=xc-this%period
  enddo

  ! if in the section with no potential barrier
  if (xc >= 0._Long ) then
   wf= this%a1*exp( imagi*this%alpha*xc &
                   + tfac*hbar**2/(2._Long*this%mass)*this%alpha**2) &
      +this%a2*exp(-imagi*this%alpha*xc &
                   + tfac*hbar**2/(2._Long*this%mass)*this%alpha**2)
  else
   wf= this%b1*exp( imagi*this%beta*xc + tfac*this%beta**2 &
                   + tfac*(hbar**2/(2._Long*this%mass)*this%beta**2 + this%v0)) &
      +this%b2*exp(-imagi*this%beta*xc + tfac*this%beta**2 &
                   + tfac*(hbar**2/(2._Long*this%mass)*this%beta**2 + this%v0))
  endif

!  write(ERROR_UNIT,*)'WfKroPen_getWav: a1,a2,wf=',this%a1,this%a2,wf

 end function WfKronigPenney_getWaveFn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 real (Long) function WfKronigPenney_dRdE(this,eek) result(dRdE)
  implicit none

  class(WfKronigPenney), intent(in) :: this
  real (Long), intent(in) :: eek !< sqrt of energy of state

  real (Long) :: hinitial  ! initial step for derivative - length scale

  real (Long) :: error  !output tolerance of derivative

!  write(*,*)'Entering function WfKronigPenney_dRdE'

  ! 1/10 of a coefficient of eek in the energyRoot equation
  hinitial = 0.1_Long/(sqrt(2._Long*this%mass)/hbar*(this%period-this%bwidth))

!  select type(this)
!   type is (WfKronigPenney)
!    currState=>this
!  end select
  dRdE=bmath_LDiffRichardson(getCurrEnergyRoot,eek,hinitial,error)

 end function WfKronigPenney_dRdE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 real (Long) function WfKronigPenney_energyRoot(this,eek) result(en)
  implicit none

  class(WfKronigPenney), intent(in) :: this
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

 end function WfKronigPenney_energyRoot

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

end module class_WfKronigPenney

