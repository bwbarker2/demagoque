! Copyright (C) 2011  Brent W. Barker
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program (gpl-3.0.txt).  If not, see
!    <http://www.gnu.org/licenses/>.
!
!    Author: Brent W. Barker
!            barker at nscl dot msu dot edu
!            National Superconducting Cyclotron Laboratory
!            Michigan State University
!            1 Cyclotron, East Lansing, MI 48824-1321

!> \brief Evolves the density matrix forward in time
SUBROUTINE time_evolution
  use input_parameters
  use mesh
  use time
  IMPLICIT NONE
  
  REAL (Long) :: dt2  ! half delt
  real (Long) :: soms5  !split operator method s_5

  integer :: ii

 firstOutput=.true.

  it=0
  t=0d0

  call calcPotDiag
  do ii=1,4
   call output
  enddo

!return

 dt2=delt*0.5d0

 soms5=1d0/3d0*(2d0+2d0**(-1d0/3d0)+2d0**(1d0/3d0))

!call mesh_setReflectedLR(.true.)

  DO it=1,Nt  !Nt  changed for debugging

!     write(*,*)'running step',it    ! moved to output.f90

     ! update current time
     t=it*delt

!   call mesh_setReflectedLR(.true.)


!call mesh_setReflectedLR(.true.)

   if(splitOperatorMethod==3)then
!    if(denState==SPACE)then
!call mesh_setReflectedLR(.true.)
!     CALL evol_x(dt2)
!call mesh_setReflectedLR(.true.)
!     call output
!     call output
!     call output
!     CALL evol_k(delt)
!call mesh_setReflectedLR(.true.)
!     call output
!     CALL evol_x(dt2)
!call mesh_setReflectedLR(.false.)
!    else
!     call setState(MOMENTUM)
     CALL evol_k(dt2)
!     call output
     CALL evol_x(delt)
!     call output
     CALL evol_k(dt2)
!    endif
   elseif(splitOperatorMethod==5)then
    call evol_x(soms5*dt2)
    call evol_k(soms5*delt)
    call evol_x((1d0-soms5)*dt2)
    call evol_k((1d0-2d0*soms5)*delt)
    call evol_x((1d0-soms5)*dt2)
    call evol_k(soms5*delt)
    call evol_x(soms5*dt2)
   endif

!   call mesh_setReflectedLR(.false.)

     call output
!   call mesh_setReflectedLR(.false.)

! NOTE: calling this subroutine during normal evolution causes divergences if ntime.ne.1
   if(useImEvol)call renormalizeDM

  ENDDO


!call mesh_setReflectedLR(.false.)

END SUBROUTINE time_evolution

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> evolves density matrix according to the semiclassical kinetic energy
SUBROUTINE evol_k(dtim)
  use input_parameters
  use mesh
  use phys_cons
  use prec_def
  IMPLICIT NONE

  real (Long), intent(in) :: dtim  !< timestep to use

  INTEGER ika, ikr !loop variables
!  real*8 :: cos2k, sin2k, xre, xim,xre2,xim2   ! cos,sin part of exp, exponent itself, den_re, den_im
  real*8 :: edt,k1,k2

!call mesh_setReflectedLR(.true.)

  call setState(MOMENTUM)

!call mesh_setReflectedLR(.true.)
!call mesh_setReflectedLR(.false.)
!  call makeMomentumHermitian()

!  write(*,*)'dtim=',dtim

!  write(*,*)'starting evol_k loop'
  !loop over all grid points

  DO ikr=-Nkr2,Nkr2-1
!  DO ikr=Nkr2-1,-Nkr2,-1  !reverse direction of indices
     
     DO ika=-Nka,Nka-1
!     DO ika=Nka-1,-Nka,-1  !reverse direction of indices
!        call getDenPtsK(ikr,ika,iikr,iika)
        
!call mesh_setReflectedLR(.true.)

        call getK12(ika,ikr,k1,k2)

       if(potFinal==3)then
        edt=sin(w*delt)/(w*delt)
       else
        edt=1d0
       endif

       if(useImEvol)then
        edt=edt*(-hbar/m0*0.5d0*(k1*k1+k2*k2)*dtim)
        call setDenK(ikr,ika,exp(edt)*getDenK(ikr,ika))
       else
        !time evolution operator = exp(-i(E-E')t/h)
        !                        = exp(-ih/2m(k^2-k'^2))
        edt=edt*(-hbar/m0*0.5d0*(k1*k1-k2*k2)*dtim)
!        edt=edt*(-hbar/m0*ka(ika)*kr(ikr)*dtim)
        if(ika==-Nka) then
         call setDenK(ikr,ika,cos(edt)*getDenK(ikr,ika))
        else
         call setDenK(ikr,ika,exp(imagi*edt)*getDenK(ikr,ika))
        endif
!        cos2k=cos(edt)
!        sin2k=sin(edt)

!if(abs(sin2k)>0.001)then
! write(*,*)'**** sin2k too big. it,ikr,ika,sin2k:',it,ikr,ika,sin2k
!endif

!        xre=DBLE(getDenK(ikr,ika))
!        xim=DIMAG(getDenK(ikr,ika))
        
        ! exp(i*edt) = cos2k + i*sin2k
!        xre2=xre*cos2k - xim*sin2k
!        xim2=xre*sin2k + xim*cos2k

!if(xim2>0.0001)then
! write(*,*)'*** xim2 > 0.0001. it,ikr,ika,xim2',it,ikr,ika,xim2
!endif

!        call setDenK(ikr,ika,cmplx(xre2,xim2,8))
!        denmat(ikr,ika)=cmplx(xre2,xim2,8)
       endif

!call mesh_setReflectedLR(.false.)
     ENDDO
     
  ENDDO

!  write(*,*)'ending evol_k loop'

!call mesh_setReflectedLR(.false.)

END SUBROUTINE evol_k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Enforces Hermiticity by setting cells equal to the average of left and right.
subroutine makeMomentumHermitian()
 use mesh
 use phys_cons
 implicit none

 integer :: ikr

 do ikr=1,Nkr2-1
  denmat(ikr,:)=(denmat(ikr,:)+conjg(denmat(-ikr,:)))*0.5d0
 enddo

 do ikr=-Nkr2+1,-1
  denmat(ikr,:)=conjg(denmat(-ikr,:))
 enddo

end subroutine makeMomentumHermitian

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!> Evolves density matrix in position space
SUBROUTINE evol_x(dtim)
 use cons_laws
 use input_parameters
 use mesh
 use phys_cons
 use prec_def
 IMPLICIT NONE

 real*8, intent(in) :: dtim !< timestep

 INTEGER :: ixr, ixa !loop variables
 real*8 :: cos2k, sin2k, udt !cos,sin part of exp, exponent itself
 real*8 :: xim,xre,xre2,xim2 !x density matrix, imaginary, real
 real*8 :: x1,x2,ux1,ux2 !position in x,x' basis, potential at x,x'
 real*8 :: cutfac   !factor for imaginary off-diagonal cutoff

 integer :: ki !used by LIN_INT

 logical :: debugxall
 real*8 :: debugudt

!call mesh_setReflectedLR(.true.)

! real*8, dimension(-Nxr2:Nxr2):: tpots !debugging vars

 ! initialize potential (esp. if no potential)
 ux1=0.0d0
 ux2=0.0d0

! tpots=0.d0

 call setState(SPACE)

!call mesh_setReflectedLR(.true.)
 !call makeSpaceHermitian()

!  write(*,*)'debug: dtim2=',dtim
! write(101,*)'# time=',t

!  write(*,*)'starting evol_k loop'

 ! calculates and stores the potential at each grid point in x
 call calcPotDiag()
! write(*,*)'potDiag:',potDiag

!call mesh_setReflectedLR(.true.)
! write(*,*)'finished calcing pot diag'

 !loop over all grid points

 DO ixr=-Nxr2,Nxr2-1

  !get imaginary cutoff factor if needed (not for adiabatic evolution)
!  if(useImCutoff.and..not.useAdiabatic)then
  if(useImCutoff)then
   call getImCutoff(cutfac, ixr,dtim)
  else
   cutfac=1d0
  endif

!  write(*,*)'timestep,ixr,cutfac=',it,ixr,cutfac

!  DO ixa=Nxa2-1,-Nxa2,-1  !reverse index direction
  DO ixa=-Nxa2,Nxa2-1

!  if(ixa.eq.-24.and.ixr.eq.-49)then
!   debugxall=.true.
!  else
   debugxall=.false.
!  endif

  if(debugxall)then
   write(*,*)denmat(ixa,ixr)-denmat(-ixa,ixr)
  endif

   call getX12(ixa,ixr,x1,x2)

   ki=1
   call LIN_INT(xa,potDiag,Nxa+1,x1,ux1,ki)
   call LIN_INT(xa,potDiag,Nxa+1,x2,ux2,ki)

! debugging test - straight-calc HO:
!     if(ixr==0) then
!      write(101,*)x1,0.5*m0*(w*x1)**2,ux1
!     endif

   if(useImEvol) then
    udt=-(ux1+ux2)*dtim/hbar
    call setDenX(ixa,ixr,cutfac*exp(udt)*getDenX(ixa,ixr))
   else
     !time evolution operator = exp(-i(U(x)-U(x'))t/h)
    udt=-(ux1-ux2)*dtim/hbar
   if(debugxall)debugudt=udt
 !   tpots(ixr)=udt  !debugging
    cos2k=dcos(udt)
    sin2k=dsin(udt)
 !     if(ixr==0)write(*,*)'ixa,ixr,cos2k,sin2k=',ixa,ixr,cos2k,sin2k
 
 !     if(ixr==0)write(*,*)'ixa,ixr,den_im-pree=',ixa,ixr,den_im(iixa,iixr)     
    xre=DBLE(getDenX(ixa,ixr))
    xim=DIMAG(getDenX(ixa,ixr))

!    if(ixa==2)then
!     if(ixr==1.or.ixr==-1)then
!      write(*,'(A,O24)')'denmat=',dble(denmat(ixa,ixr))
!     endif
!    endif
 
    ! exp(i*edt) = cos2k + i*sin2k
    xre2=xre*cos2k - xim*sin2k
    xim2=xre*sin2k + xim*cos2k
 !   if(ixr==0)write(*,*)'ixa,ixr,den_im-post=',ixa,ixr,den_im(iixa,iixr)


 
    call setDenX(ixa,ixr,cutfac*cmplx(xre2,xim2,8))

   if(debugxall.and.ixa==24.and.ixr==-49)then
    write(*,*)'debugudt,udt',debugudt,udt
   endif

   if(debugxall.and.ixa>0)then
    if(dble(denmat(ixa,ixr)).ne.dble(denmat(-ixa,ixr)))then
     write(*,*)'den ne!',ixa,ixr,denmat(ixa,ixr)-denmat(-ixa,ixr)
    endif
   endif
!    if(ixa>0)then
!     if(cutfac*xre2-dble(denmat(-ixa,ixr))>1e-16)then
!      write(*,*)ixa,ixr,x1,x2,cutfac*xre2-dble(denmat(-ixa,ixr))
!     endif !>1e-16
!    endif !ixa>0

!    if(ixa==2)then
!     if(ixr==1.or.ixr==-1)then
!      write(*,*)ixa,ixr,x1,x2,ux1,ux2,xre2
!      write(*,'(A,O30)')'udt=',udt
!      write(*,'(A,O24,O24)')'cos2k,sin2k=',cos2k,sin2k
!     endif
!    endif

   endif !not useImEvol

   !find maximum value of imaginary component, BWB 2011-03-11
   if(DIMAG(getDenX(ixa,ixr))>maxxim)maxxim=DIMAG(getDenX(ixa,ixr))
!   if(ixr>0) then
!    write(*,*)'ixr,tpots-diff:',ixr,tpots(ixr)+tpots(-ixr)
!   endif

!call mesh_setReflectedLR(.false.)

   ENDDO
 ENDDO

 !copy cells to extra redundant parts
 call copyExtra

! write(*,*)
! write(*,*)

END SUBROUTINE evol_x

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Enforces Hermiticity by setting cells equal to the average of top and bottom.
subroutine makeSpaceHermitian()
 use mesh
 use phys_cons
 implicit none

 integer :: ixr

 do ixr=1,Nxr2-1
  denmat(:,ixr)=(denmat(:,ixr)+conjg(denmat(:,-ixr)))*0.5d0
 enddo

 do ixr=-Nxr2+1,-1
  denmat(:,ixr)=conjg(denmat(:,-ixr))
 enddo

 call copyExtra

end subroutine makeSpaceHermitian

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Calculates and stores the potential at each grid point in x
subroutine calcPotDiag()
 use input_parameters
 use mesh
 use prec_def
 implicit none

 integer :: ixa
 real (Long) :: weight !< weighting for adiabatic switching
 real (Long) :: potI, potF !< potentials
 real (Long) :: getWeight !< functions

!call mesh_setReflectedLR(.true.)

! write(*,*)'debug: starting calcPotDiag'
 if (useAdiabatic)then
  weight=getWeight()
 else
  weight=1d0
 endif


! write(*,*)'time,weight:',t,weight,1.0-weight

 do ixa=-Nxa2,Nxa2-1
!  write(*,*)'debug: ixa=',ixa
  if (useAdiabatic) then
!call mesh_setReflectedLR(.true.)
   call getPotX(potI,potInitial,ixa)
   call getPotX(potF,potFinal,ixa)
   potDiag(ixa)=weight*potI + (1.d0-weight)*potF
!  write(*,*)'sofar:',ixa,potDiag(ixa)
!   write(*,*)potDiag(ixa),getPotX(potInitial,ixa)
!call mesh_setReflectedLR(.false.)
  else
!  write(*,*)'debug: ixa,potDiag initial=',ixa,potDiag(ixa)
   call getPotX(potDiag(ixa),potFinal,ixa)
!   write(*,*)'debug: ixa,potDiag:',ixa,potDiag(ixa)
  endif
 enddo
 potDiag(Nxa2)=potDiag(-Nxa2)

!call mesh_setReflectedLR(.false.)

end subroutine calcPotDiag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11

subroutine getImCutoff(cutfac, ixr,dtim)
 use input_parameters
 use mesh
 use phys_cons
 implicit none

 integer, intent(in) :: ixr
 real*8, intent(out) :: cutfac
 real*8, intent(in) :: dtim !timestep

 real*8 :: xxr

 xxr=abs(xr(ixr))

 if(xxr<=cutoff_x0) then
  cutfac=0d0
 elseif(xxr<=cutoff_x0+cutoff_d0/2d0) then
  cutfac=2d0*(xxr-cutoff_x0)**2/cutoff_d0**2
 elseif(xxr<=cutoff_x0+cutoff_d0) then
  cutfac=1d0-2d0*(xxr-(cutoff_x0+cutoff_d0))**2/cutoff_d0**2
 else
  cutfac=1d0
 endif
! write(*,*)cutfac
 cutfac=exp(-2d0*cutoff_w0*cutfac*dtim/hbar)
! write(*,*)ixr,xxr,cutfac

end subroutine getImCutoff


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111

real*8 function getWeight()
 use input_parameters
 use time
 implicit none

! real (Long), intent(out) :: weight

 getWeight=1.d0/(1.d0+exp((t-tad)/wtad))

end function getWeight

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine getPotX(potX,potType,ix)
 use bexception
 use input_parameters
 use mesh
 use prec_def
 implicit none

 real (Long), intent(out) :: potX
 integer,     intent(in)  :: potType
 integer,     intent(in)  :: ix
! real (Long) :: potSkyrme  !functions

 select case (potType)
  case (-1)
   potX=0d0
  case (0)
   call potHO(potX,ix)
  case (1)
   call potHOmf(potX,ix)
  case (2)
   call potSkyrme(potX,ix)
  case (3)
   call potHOexact(potX,ix)
  case (4)
   call potBEC_1D_HO_Mateo2011(potX,ix,ho_mateo_wz,ho_mateo_wt,ho_mateo_scat,ho_mateo_Npart)
  case default
   call throwException('getPotX: improper potential type', BEXCEPTION_FATAL)
   potX=0.d0
 end select

end subroutine getPotX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> calculates potential due to a 3D harmonic trap, but with trapping potential
!! much greater in z direction, such that it can be described with modified
!! 1D GPE.
subroutine potBEC_1D_HO_Mateo2011(potX,ix,wz,wt,scat,Npart)
 use mesh
 use phys_cons
 use prec_def
 implicit none

 real (Long),  intent(out) :: potX
 integer,      intent(in)  :: ix

 !> angular frequency of harmonic trap in z (elongated) direction
 real (Long),  intent(in)  :: wz

 !> angular frequency of harmonic trap in transverse direction
 real (Long),  intent(in)  :: wt

 !> s-wave scattering length of particle
 real (Long), intent(in) :: scat

 !> total number of particles in condensate
 real (Long), intent(in) :: Npart

!potX=1
 potX=0.5_Long*m0*wz**2*xa(ix)**2 &
      +hbar*wt*sqrt(1_Long+4_Long*scat*Npart*getDenX(ix,0))
! write(*,*)'m0,wz,xa=',m0,wz,xa(ix)

end subroutine potBEC_1D_HO_Mateo2011

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine potHO(potX,ix)
!! pot_ho - Harmonic oscillator potential in x,x' basis, centered at x=0
  use input_parameters
  use phys_cons
  use prec_def
  use mesh
  implicit none

  real*8,  intent(out) :: potX
  integer, intent(in)  :: ix

!  write(*,*)'debug: m0,w,ix,delxa',m0,w,ix,delxa

  potX=0.5d0*m0*(w*xa(ix))**2
!  write(*,*)potX

end subroutine potHO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine potHOexact(potX,ix)
 !! potHOexact - computes exact evolution of external HO, from Chin, Krotsheck, Phys. Rev. E 72 (2005) 036705.
 use input_parameters
 use mesh
 use phys_cons
 implicit none

 real*8,  intent(out) :: potX
 integer, intent(in)  :: ix

 potX=0.5d0*m0*(w*xa(ix))**2*2*(1-cos(w*delt))/(w*delt*sin(w*delt))

end subroutine potHOexact

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine potHOmf(potX,ixa1)
  !! calcPotMF - calculates meanfield potential at all diagonal points, stores
  ! in potMF
  use input_parameters
  use mesh
  use phys_cons
  use prec_def
  implicit none

  real (Long), intent(out) :: potX
  integer, intent(in) :: ixa1
  integer :: ixa2,id,itry   ! loop variables
!  real (Long) :: testtot

!  testtot=0

!  iixa1=iNxa2(ixa1)
  potX=0

  do id=-Nxa2,Nxa2      ! id is the difference between ixa1 and ixa2
   ! to obey periodic boundary conditions, it is necessary to construct the
   ! meanfield with this periodicity, and treat the distances properly.
   itry=ixa1+id
   if(itry.LT.-Nxa2)then
    ixa2=itry+Nxa+1
   elseif(itry.GT.Nxa2)then
    ixa2=itry-Nxa-1
   else
    ixa2=itry
   endif
      

    !write(*,*)ixa1,id,ixa2
!      potMF(ixa1)=potMF(ixa1)+den_re(iixa1,iixr0)*den_re(iixa2,iixr0) &
   potX=potX+DBLE(getDenX(ixa2,0)) &
                           *(xa(id))**2
!      testtot=testtot+den_re(iixa1,iixr0)*den_re(iixa2,iixr0) &
!                      *2*xa(ixa1)*xa(ixa2)
  enddo !id
!    read(*,*)itry
  potX=potX*0.25d0*m0*w*w/(Nmax+1)
!    write(*,*)'ixa,potMF=',ixa1,potMF(ixa1)


!  write(*,*)'testtot=',testtot

end subroutine potHOmf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111111111

subroutine potSkyrme(potX,ix)
 use mesh
 use prec_def
 use skyrme_params
 implicit none

 real*8,  intent(out) :: potX
 integer, intent(in)  :: ix
 real*8               :: xxr ! density corrected for dimensionality
 real*8               :: skyContact  !contact skyrme potential



 xxr=DBLE(getDen(ix,0))*facd
! write(*,*)ix,facd,den_re(iNxa2(ix),iNxr2(0)),xxr
! write(*,*)xxr,abs(xxr)


! !prescription for negative density 1
! if(xxr <= 0) then
!  potX=0.0d0
! else
!  potX=skyContact(xxr)
! endif

 !prescription for negative density 2
 potX=skyContact(xxr)
 if(xxr<0.d0)potX=potX*(-1.d0)

! !prescription 3 - if surrounding cells are opposite sign, then pot=0
! if(ix.lt.Nxa2.and.ix.gt.-Nxa2) then
!  if(xxr*den_re(iNxa2(ix+1),iNxr2(0))<0.d0 &
!     .and.xxr*den_re(iNxa2(ix-1),iNxr2(0))<0.d0) potX=0.d0
! elseif(ix.eq.-Nxa2) then
!  if(xxr*den_re(iNxa2(ix+1),iNxr2(0))<0.d0) potX=0.d0
! elseif(ix.eq.Nxa2) then
!  if(xxr*den_re(iNxa2(ix-1),iNxr2(0))<0.d0) potX=0.d0
! else
!  write(*,*)'counted something wrong in potSkryme'
! endif
  

end subroutine potSkyrme

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real*8 function skyContact(rho)
 use skyrme_params

 real*8, intent(in) :: rho

 skyContact=0.75d0*t0*dabs(rho) &
       + (2.0d0+sig)/16.0d0*t3*dabs(rho)**(1.0d0+sig)
end function skyContact

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!subroutine enforceHermiticityX
! !! enforceHermiticityX - The real part of the spatial density matrix should be symmetric around the xr=0 axis, and the imaginary part should be antisymmetric. This subroutine enforces this by averaging the 2 values that should be equal in magnitude and assigning that value to both of them.
! use mesh
! implicit none
!
! integer :: ixa,ixr,iixa,iixr,iixrn
! real*8  :: ave
!
! do ixa=-Nxa2,Nxa2-1
!  iixa=iNxa2(ixa)
!  do ixr=1,Nxr2
!   iixr=iNxr2(ixr)
!   iixrn=iNxr2(-ixr)
!   ave=0.5d0*(den_re(iixa,iixr)+den_re(iixa,iixrn))
!   den_re(iixa,iixr)=ave
!   den_re(iixa,iixrn)=ave
!   ave=0.5d0*(den_im(iixa,iixr)-den_im(iixa,iixrn))
!   den_im(iixa,iixr)=ave
!   den_im(iixa,iixrn)=-ave
!  enddo
! enddo
!end subroutine enforceHermiticityX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!subroutine enforceHermiticityK
! !! enforceHermiticityK - The real part of the spectral density matrix should be symmetric around the kr=0 axis, and the imaginary part should be antisymmetric. This subroutine enforces this by averaging the 2 values that should be equal in magnitude and assigning that value to both of them.
! use mesh
! implicit none
!
! integer :: ika,ikr,iika,iikr,iikrn
! real*8  :: ave
!
! do ika=-Nka2,Nka2-1
!  iika=iNka2(ika)
!  do ikr=1,Nkr2
!   iikr=iNkr2(ikr)
!   iikrn=iNkr2(-ikr)
!   ave=0.5d0*(den_re(iikr,iika)+den_re(iikrn,iika))
!   den_re(iikr,iika)=ave
!   den_re(iikrn,iika)=ave
!   ave=0.5d0*(den_im(iikr,iika)-den_im(iikrn,iika))
!   den_im(iikr,iika)=ave
!   den_im(iikrn,iika)=-ave
!  enddo
! enddo
!end subroutine enforceHermiticityK
