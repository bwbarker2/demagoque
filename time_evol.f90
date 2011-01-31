SUBROUTINE time_evolution
  !! time_evolution - evolves the density matrix forward in time
  use time
  IMPLICIT NONE
  
  REAL (Long) :: dt2  ! half delt

!  integer :: ii

  it=0
  call calcPotDiag
  call output
  call output

  DO it=1,Nt  !Nt  changed for debugging

!     write(*,*)'running step',it    ! moved to output.f90

     ! update current time
     t=it*delt
     dt2=delt*0.5d0

     ! evolve in k-space
     CALL evol_k(dt2)
!     write(*,*)'finished first evol_x'
!     call enforceHermiticityK
!     call output

     ! evolve in x-space
     CALL evol_x(delt)
!     write(*,*)'finished first evol_k'
!     call enforceHermiticityX
!     CALL output

     CALL evol_k(dt2)
!     write(*,*)'finished second evol_x'
!  call howHermitian
!     call enforceHermiticityK
     call output
  ENDDO

END SUBROUTINE time_evolution

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

!subroutine howHermitian
!!! howHermitian - Gives number characterizing Hermiticity of the density matrix.
! use mesh
! use time
! implicit none
!
!! complex*16 :: zdet2d   !complex double determinant function (in bmath.f90)
! complex*16, dimension(:,:), allocatable :: denmatA, denmatZ   !density matrix in human-readable form
! complex*16 :: mag  !subtraction of density matrix and its conjugate
!                      !transpose. Should be zero if Hermitian
! integer :: i1,i2,ii1,ii2  !loop variables
!
! allocate(denmatA(-Nxa2:Nxa2,-Nxr2:Nxr2))
! allocate(denmatZ(-Nxa2:Nxa2,-Nxr2:Nxr2))
!
! do i1=-Nxa2,Nxa2
!  do i2=-Nxr2,Nxr2
!   call getDenPts(i1,i2,ii1,ii2)
!   denmatA(i1,i2)=getDenX(i1,i2)
!  enddo
! enddo
!
! !subtract the conjugate transpose
! do i1=-Nxa2,Nxa2
!  do i2=-Nxr2,Nxr2
!   denmatZ(i1,i2)=denmatA(i1,i2)-dconjg(denmatA(i1,-i2))
!  enddo
! enddo
!
! write(71,*)'# time = ',t,'fm/c'
! do i2=-Nxr2,Nxr2
!  write(71,*) (DBLE(denmatZ(i1,i2)),i1=-Nxa2,Nxa2)
! enddo
!
! write(71,*)
! write(71,*)
!
! !square each element
!! denmatZ=denmatZ**2
!
! mag=0.d0
!
! !sum the squares
! do i1=-Nxa2,Nxa2
!  do i2=-Nxr2,Nxr2
!   mag=mag+denmatZ(i1,i2)**2
!  enddo
! enddo
!
! write(*,*)'square-sums:',mag
! 
!
!!919 FORMAT(1000e16.5)
!
!
!end subroutine howHermitian

!subroutine howHermitianX12
!!! howHermitianX12 - Gives number characterizing Hermiticity of the density matrix.
! ! This code is incorrect for the rotated matrix, as it assumes the (x,x') coordinate system. See howHermitian code for correct subroutine.
! use mesh
! use time
! implicit none
!
!! complex*16 :: zdet2d   !complex double determinant function (in bmath.f90)
! complex*16, dimension(:,:), allocatable :: denmatA, denmatZ   !density matrix in human-readable form
! complex*16 :: mag  !subtraction of density matrix and its conjugate
!                      !transpose. Should be zero if Hermitian
! integer :: i1,i2,ii1,ii2  !loop variables
!
! denmat=cmplx(den_re,den_im,8)
!
! allocate(denmatA(-Nxa2:Nxa2,-Nxr2:Nxr2))
! allocate(denmatZ(-Nxa2:Nxa2,-Nxr2:Nxr2))
!
! do i1=-Nxa2,Nxa2
!  ii1=iNxa2(i1)
!  do i2=-Nxr2,Nxr2
!   ii2=iNxr2(i2)
!   denmatA(i1,i2)=denmat(ii1,ii2)
!  enddo
! enddo
!
! !subtract the conjugate transpose
! denmatZ=denmatA-DCONJG(transpose(denmatA))
!
! !square each element
! denmatZ=denmatZ**2
!
! mag=0.d0
!
! !sum the squares
! do i1=-Nxa2,Nxa2
!  do i2=-Nxr2,Nxr2
!   mag=mag+denmatZ(i1,i2)
!  enddo
! enddo
!
! write(*,*)'square-sums:',mag
! 
! write(71,*)'# time = ',t,'fm/c'
! do i2=-Nxr2,Nxr2
!  write(71,919) (DBLE(denmatZ(i1,i2)),i1=-Nxa2,Nxa2)
! enddo
!
! write(71,*)
! write(71,*)
!
!919 FORMAT(1000e16.5)
!
!! if(Nxa==Nxr)write(*,*)zdet2d(denmatA,Nxa2+Nxa2+1)
!
!end subroutine howHermitianX12

SUBROUTINE evol_k(dtim)
  !! evol_k - evolves density matrix according to the semiclassical kinetic energy
  use mesh
  use phys_cons
  use prec_def
  use time
  IMPLICIT NONE

  real*8, intent(in) :: dtim

  INTEGER ika, ikr !loop variables
  real*8 :: cos2k, sin2k, edt, xre, xim,xre2,xim2   ! cos,sin part of exp, exponent itself, den_re, den_im
  real*8 :: k1,k2

  call setState(MOMENTUM)

!  write(*,*)'dtim=',dtim

!  write(*,*)'starting evol_k loop'
  !loop over all grid points
  DO ikr=-Nkr2,Nkr2-1
     
     DO ika=-Nka,Nka-1
!        call getDenPtsK(ikr,ika,iikr,iika)
        
        call getK12(ika,ikr,k1,k2)

        !time evolution operator = exp(-i(E-E')t/h)
        !                        = exp(-ih/2m(k^2-k'^2))
        edt=-hbc/m0*0.5d0*(k1*k1-k2*k2)*dtim !*8d0 is for testing - BWB 2010-12-24
        cos2k=dcos(edt)
        sin2k=dsin(edt)
!        if(sin2k.gt.1.or.sin2k.lt.-1)write(*,*)cos2k

        xre=DBLE(getDenK(ikr,ika))
        xim=DIMAG(getDenK(ikr,ika))
        
        ! exp(i*edt) = cos2k + i*sin2k
        xre2=xre*cos2k - xim*sin2k
        xim2=xre*sin2k + xim*cos2k

        call setDenK(ikr,ika,cmplx(xre2,xim2,8))

     ENDDO
     
  ENDDO

!  write(*,*)'ending evol_k loop'

END SUBROUTINE evol_k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE evol_x(dtim)
  !! evol_x - evolves density matrix in position space
 use cons_laws
 use mesh
 use osc_pars
 use phys_cons
 use prec_def
 use time
 IMPLICIT NONE

 real*8, intent(in) :: dtim !timestep

 INTEGER :: ixr, ixa !loop variables
 real*8 :: cos2k, sin2k, udt !cos,sin part of exp, exponent itself
 real*8 :: xim,xre,xre2,xim2 !x density matrix, imaginary, real
 real*8 :: x1,x2,ux1,ux2 !position in x,x' basis, potential at x,x'

 integer :: ki !used by LIN_INT

! real*8, dimension(-Nxr2:Nxr2):: tpots !debugging vars

 ! initialize potential (esp. if no potential)
 ux1=0.0d0
 ux2=0.0d0

! tpots=0.d0

 call setState(SPACE)

!  write(*,*)'debug: dtim2=',dtim
! write(101,*)'# time=',t

!  write(*,*)'starting evol_k loop'

 ! calculates and stores the potential at each grid point in x
 call calcPotDiag()

! write(*,*)'finished calcing pot diag'

 !loop over all grid points
 DO ixa=-Nxa2,Nxa2-1

  DO ixr=-Nxr2,Nxr2-1
   call getX12(ixa,ixr,x1,x2)

   ki=1
   call LIN_INT(xa,potDiag,Nxa+1,x1,ux1,ki)
   call LIN_INT(xa,potDiag,Nxa+1,x2,ux2,ki)

! debugging test - straight-calc HO:
!     if(ixr==0) then
!      write(101,*)x1,0.5*m0*(w*x1)**2,ux1
!     endif

     !time evolution operator = exp(-i(U(x)-U(x'))t/h)
   udt=-(ux1-ux2)*dtim/hbc
!   tpots(ixr)=udt  !debugging
   cos2k=dcos(udt)
   sin2k=dsin(udt)
!     if(ixr==0)write(*,*)'ixa,ixr,cos2k,sin2k=',ixa,ixr,cos2k,sin2k

!     if(ixr==0)write(*,*)'ixa,ixr,den_im-pree=',ixa,ixr,den_im(iixa,iixr)     
   xre=DBLE(getDenX(ixa,ixr))
   xim=DIMAG(getDenX(ixa,ixr))

   ! exp(i*edt) = cos2k + i*sin2k
   xre2=xre*cos2k - xim*sin2k
   xim2=xre*sin2k + xim*cos2k
!   if(ixr==0)write(*,*)'ixa,ixr,den_im-post=',ixa,ixr,den_im(iixa,iixr)

   call setDenX(ixa,ixr,cmplx(xre2,xim2,8))

!   if(ixr>0) then
!    write(*,*)'ixr,tpots-diff:',ixr,tpots(ixr)+tpots(-ixr)
!   endif

   ENDDO
 ENDDO

 !copy cells to extra redundant parts
 call copyExtra

! write(*,*)
! write(*,*)

END SUBROUTINE evol_x

subroutine calcPotDiag()
 use mesh
 use prec_def
 use time
 implicit none

 integer :: ixa
 real (Long) :: weight ! weighting for adiabatic switching
 real (Long) :: potI, potF !potentials
 real (Long) :: getWeight !functions

! write(*,*)'debug: starting calcPotDiag'
 if (useAdiabatic) weight=getWeight()

! write(*,*)'time,weight:',t,weight,1.0-weight

 do ixa=-Nxa2,Nxa2
!  write(*,*)'debug: ixa=',ixa
  if (useAdiabatic) then
   call getPotX(potI,potInitial,ixa)
   call getPotX(potF,potFinal,ixa)
   potDiag(ixa)=weight*potI + (1.d0-weight)*potF
!  write(*,*)'sofar:',ixa,potDiag(ixa)
!   write(*,*)potDiag(ixa),getPotX(potInitial,ixa)
  else
!  write(*,*)'debug: ixa,potDiag initial=',ixa,potDiag(ixa)
   call getPotX(potDiag(ixa),potFinal,ixa)
!   write(*,*)'debug: ixa,potDiag:',ixa,potDiag(ixa)
  endif
 enddo
! potDiag(Nxa2)=potDiag(-Nxa2)
end subroutine calcPotDiag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11

real*8 function getWeight()
 use time
 implicit none

! real (Long), intent(out) :: weight

 getWeight=1.d0/(1.d0+dexp((t-tad)/wtad))

end function getWeight

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine getPotX(potX,potType,ix)
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
  case default
   write(*,*)'getPotX: Improper potential type:',potType,'Assuming no potential'
   potX=0.d0
 end select

end subroutine getPotX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine potHO(potX,ix)
!! pot_ho - Harmonic oscillator potential in x,x' basis, centered at x=0
  use osc_pars
  use phys_cons
  use prec_def
  use mesh
  implicit none

  real*8,  intent(out) :: potX
  integer, intent(in)  :: ix

!  write(*,*)'debug: m0,w,ix,delxa',m0,w,ix,delxa

  potX=0.5d0*m0*(w*xa(ix))**2

end subroutine potHO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine potHOmf(potX,ixa1)
  !! calcPotMF - calculates meanfield potential at all diagonal points, stores
  ! in potMF
  use mesh
  use osc_pars
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
                           *(id*delxa)**2
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
 use time
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
