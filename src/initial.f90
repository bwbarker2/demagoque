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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> calculates initial quantities derived from input
SUBROUTINE calcInitial
  use cons_laws
  use input_parameters
  use mesh
  use phys_cons
  implicit none

  call initializeMesh

  write(*,*) '1D to 3D factor:',facd

  ! oscillator data
  if(initState_gaussianNuclear) then
   w=hbar/m0*6e0_Long*(rho0/facd)**2/(Nmax+1e0_Long)!   *0.1
   whm=m0*w/hbar
  endif

  if(potFinal==4) then
!   w=ho_mateo_wz
!   whm=m0*w/hbar
!  endif
!  whm=2.25e-2
!  w=2.5e-5_Long  !hard-coded for Ott2001
  w=0.00016_Long
!  w=0.088
!  w=whm*hbar/m0
  whm=m0*w/hbar
 endif

  write(*,*) 'Parameters of the calculation:'
  write(*,*) 'dxr=',delxr,'dxa=',delxa,'dkr=',delkr, &
             'dka=',delka,'whm=',whm,'w=',w
  
  
  allocate(potx(-Nxa2:Nxa2))
  
end SUBROUTINE calcInitial

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine initialState
  use input_parameters
  use mesh
  use phys_cons
  use prec_def
  implicit none

  real (Long) :: wfnho,x1,x2
  complex (Long)  :: y1,y2,den0
  integer :: ixa, ixr, iin
  integer :: initState_kdelta_index

  if(initState_kdelta) then
   initState_kdelta_index=getNearestIndexX(initState_kdelta_x0)
  else
   initState_kdelta_index=99
  endif

  do ixa=Nxan,Nxax
   do ixr=Nxrn,Nxrx


        !convert to x,x' representation
        x1=xx1(ixa,ixr)
        x2=xx2(ixa,ixr)

        den0=0e0_Long
        y1=0e0_Long
        y2=0e0_Long

        if(initState_gaussianNuclear.OR.initState_gaussian) then
         do iin=0,Nmax
           y1=wfnho(x1,iin,whm)
           y2=wfnho(x2,iin,whm)
           den0=den0+y1*y2
         enddo !iin
        endif

          if(initState_cosine) then
           y1=sqrt(initState_cosine_norm/xLa) &
              *cos(initState_cosine_number*pi*x1/xLa + initState_cosine_shift)
           y2=sqrt(initState_cosine_norm/xLa) &
              *cos(initState_cosine_number*pi*x2/xLa + initState_cosine_shift)
           den0=den0+y1*y2
          endif

          if(initState_kdelta) then
           !if we are on the diagonal and at the delta index
           if(abs(x1-x2)<(0.1_Long*delxr) &
              .AND.ixa==initState_kdelta_index) then
            den0=den0+initState_kdelta_norm*2e0_Long*pi/kLa
           endif
          endif

          if(initState_plane) then
           y1=sqrt(initState_plane_norm*0.5_Long/xLa) &
              *exp(-imagi*initState_plane_number*pi*x1/xLa &
                   +initState_plane_shift)
           y2=sqrt(initState_plane_norm*0.5_Long/xLa) &
              *exp(-imagi*initState_plane_number*pi*x2/xLa &
                   +initState_plane_shift)
           den0=den0+conjg(y1)*y2
          endif

          if(useInitState_SqWell) then
           y1=initState_sqWell%getWavefn(x1)
           y2=initState_sqWell%getWavefn(x2)
           den0=den0+y1*y2
          endif

!           if(useImEvol) then
!            ! multiply by e^(2 E t / hbar) to get tr(rho)=1 after imaginary evolution, see notes BWB 2011-02-22.
!            den0=den0+y1*y2*exp(2d0*w*(iin+0.5d0)*delt*Nimev)
!           else
!            den0=den0+y1*y2
!           endif

!        if(abs(den0).lt.1e-40) den0=0.0d0

        call setDenX(ixa,ixr,den0)
!if(ixr==ixr0) write(*,*)x1,x2,den0
!        if(ixa==2) then
!         if(ixr==1.or.ixr==-1)then
!          write(*,'(I3,I3,O24,O24)')ixa,ixr,den0,dble(denmat(ixa,ixr))
!         endif
!        endif

     enddo !ixr
  enddo !ixa

  if(useMeshXAR2) call mesh_xar2_setZeroesX

  ! set density matrix to X state
  denState=SPACE

end subroutine initialState

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine copyExtra
 !! copyExtra - copies matrix to extra part (>Nxr2 and <-Nxr2)
 use input_parameters
 use mesh
 implicit none

 integer :: ixr,ixa

 !left side
 do ixa=-Nxa2,-1

  !upper left corner
  do ixr=Nxr2,Nxr-1
   denmat(ixa,ixr)=denmat(ixa+Nxa2,ixr-Nxr)
  enddo

  !lower left corner
  do ixr=-Nxr,-Nxr2-1
   denmat(ixa,ixr)=denmat(ixa+Nxa2,ixr+Nxr)
  enddo

 enddo

 !right side
 do ixa=0,Nxa2-1

  !upper right corner
  do ixr=Nxr2,Nxr-1
   denmat(ixa,ixr)=denmat(ixa-Nxa2,ixr-Nxr)
  enddo

  !lower right corner
  do ixr=-Nxr,-Nxr2-1
   denmat(ixa,ixr)=denmat(ixa-Nxa2,ixr+Nxr)
  enddo

 enddo

 if(.not.useMeshShifted) then
  !enforce numerical hermiticity
  denmat(:,Nxr2)=0.5_Long*(denmat(:,-Nxr2)+denmat(:,Nxr2))
  denmat(:,-Nxr2)=denmat(:,Nxr2)
 endif

end subroutine copyExtra

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine boost
  use input_parameters
  use mesh
  use phys_cons
  implicit none

  INTEGER ixa,ixr !loop variables
  real (Long) :: kea !momentum, given ea

  complex (Long) :: epx

  call setState(SPACE)

  ! non-relativistic conversion
  kea=sign(1e0_Long,ea)*sqrt(2e0_Long*m0*abs(ea))/hbar

  !loop over all grid points
  do ixa=Nxan,Nxax
   DO ixr=Nxrn,Nxrx
    epx=exp(imagi*kea*(xx1(ixa,ixr)-xx2(ixa,ixr)))

    denmat(ixa,ixr)=denmat(ixa,ixr)*epx
   enddo
  ENDDO
!call mesh_setReflectedLR(.false.)

!     DO ixa=-Nxa2,Nxa2-1
     
        
!  call mesh_setReflectedLR(.true.)
        !convert into x,x' coordinates
!        call getX12(ixa,ixr,x1,x2)

        !dispacement operator = exp(-iK(X'-X))
!        epx=kea*(x1-x2)
!        epx=kea*xr(ixr)
!        denmat(ixa,ixr)=denmat(ixa,ixr)*exp(imagi*epx)
!  call mesh_setReflectedLR(.false.)
!        udt=m0*w**2*2*xa(ixa)*xr(ixr)*dtim/hbar/2.0_Long
!        udt=0
!        cos2k=cos(epx)
!        sin2k=sin(epx)

!        xre=DBLE(getDenX(ixa,ixr))
!        xim=DIMAG(getDenX(ixa,ixr))
        
        ! exp(i*edt) = cos2k + i*sin2k
!        xre2=xre*cos2k - xim*sin2k
!        xim2=xre*sin2k + xim*cos2k

!        call setDenX(ixa,ixr,cmplx(xre2,xim2,8))

!     ENDDO
     

! call copyExtra


end subroutine boost

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine displaceLeft(nx)
 !! displace - shifts density matrix by 'nx' spatial indices in -xa direction
 use mesh
 implicit none

 integer, intent(in) :: nx

 complex (Long), dimension(0:nx-1,Nxrn:Nxrx) :: wmat  !working matrix

 integer :: i

 !copy cells that are shifted 'off' the left (negative) side of the matrix
 do i=0,nx-1
  wmat(i,:)=denmat(Nxan+i,:)
 enddo !i

 !shift cells within denmat
 do i=Nxan,Nxax-nx
  denmat(i,:)=denmat(i+nx,:)
 enddo !i

 !copy periodically shifted cells to other (positive) side of matrix
 do i=Nxax-(nx-1),Nxax
  denmat(i,:)=wmat(i-Nxax+nx-1,:)
 enddo

end subroutine displaceLeft

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine displaceRight(nx)
 !! displace - shifts density matrix by 'nx' spatial indices in +xa direction
 use mesh
 implicit none

 integer, intent(in) :: nx

 complex (Long), dimension(0:nx-1,-Nxr:Nxr-1) :: wmat  !working matrix

 integer :: i

 !copy cells that are shifted 'off' the right (positive) side of the matrix
 do i=0,nx-1
  wmat(i,:)=denmat(Nxa2-nx+i,:)
 enddo !i

 !shift cells within denmat
 do i=Nxa2-1,-Nxa2+nx,-1
  denmat(i,:)=denmat(i-nx,:)
 enddo !i

 !copy periodically shifted cells to other (negative) side of matrix
 do i=-Nxa2,-Nxa2+nx-1
  denmat(i,:)=wmat(i+Nxa2,:)
 enddo

end subroutine displaceRight

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine flipclone
 !! flipclone - adds the density matrix to its clone, flipped across and conjugated (to reverse momentum) 
 use mesh
 implicit none

 complex (Long), dimension(-Nxa2:Nxa2-1,-Nxr2:Nxr2-1) :: den2

 integer :: ixa

 write(*,*)'flipclone started'

 den2=denmat(:,-Nxr2:Nxr2-1)

 write(*,*)'den2 assigned'

 denmat(-Nxa2,-Nxr2:Nxr2-1)=denmat(-Nxa2,-Nxr2:Nxr2-1) &
                            +conjg(denmat(-Nxa2,-Nxr2:Nxr2-1))
 write(*,*)'first denmat column multiply'
 do ixa=-Nxa2+1,Nxa2-1
  denmat(ixa,-Nxr2:Nxr2-1)=denmat(ixa,-Nxr2:Nxr2-1)+conjg(den2(-ixa,-Nxr2:Nxr2-1))
 enddo
 write(*,*)'do-loop finished'

 call copyExtra
 write(*,*)'copyExtra finished'

end subroutine flipclone

!subroutine displace(xc)
! !! displace - shifts density matrix by 'xc' fm in xa direction
! use mesh
! implicit none
!
! real*8, intent(in) :: xc
!
! integer                       :: ixa, ixr, iixr, ki
! real*8, dimension(Nxa,Nxr)    :: d2_re, d2_im  !working matrix
! real*8, dimension(-Nxa2:Nxa2) :: row_re, row_im
! real*8, dimension(-Nxa2:Nxa2) :: rxc !vector of shifted coordinates
!
! ki=1
!
! if(.not.isDenX) then
!  write(*,*)'displace called on non x-space matrix, wtf?'
!  stop
! endif
!
! do ixa=-Nxa2,Nxa2
!  rxc(ixa)=xa(ixa)-xc
!  if(rxc(ixa)<xa(-Nxa2))then
!   rxc(ixa)=rxc(ixa)+2*xa(Nxa2)
!  elseif(rxc(ixa)>xa(Nxa2)) then
!   rxc(ixa)=rxc(ixa)-2*xa(Nxa2)
!  endif
! enddo
!
!
! do ixr=-Nxr2,Nxr2
!  !construct row vector in -Nxa2,Nxa2 form
!  do ixa=-Nxa2,Nxa2
!   row_re(ixa)=DBLE(getDen(ixa,ixr))
!   row_im(ixa)=DIMAG(getDen(ixa,ixr))
!  enddo
!
!  do ixa=-Nxa2,Nxa2-1
!   call LIN_INT(xa,row_re,Nxa+1,rxc(ixa),d2_re(iNxa2(ixa),iixr),ki)
!   call LIN_INT(xa,row_im,Nxa+1,rxc(ixa),d2_im(iNxa2(ixa),iixr),ki)
!  enddo
! enddo
!
! denmat=cmplx(d2_re,d2_im,8)
!
!end subroutine displace

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine getX12(ixa, ixr, x1, x2)
!! getX12 - converts ixa, ixr grid indices to values of x1,x2 space and
!maintains periodicity.
 use mesh
 use prec_def
 implicit none

 integer,     intent(in)  :: ixa, ixr !human-read index of grid in xa,xr space
 real (Long), intent(out) :: x1,x2    !x,x' position in x,x' space

 real (Long) :: xa1,xr1

 xa1=xa(ixa)
 xr1=xr(ixr)

 call getrX12(xa1,xr1,x1,x2)

! x1=xa(ixa)+0.5d0*xr(ixr)
! x2=xa(ixa)-0.5d0*xr(ixr)
!
! ! if xx is outside the box, move it in periodically
! if(x1.gt.xLa)x1=x1-xLa-xLa
! if(x1.lt.-xLa)x1=x1+xLa+xLa
! if(x2.gt.xLa)x2=x2-xLa-xLa
! if(x2.lt.-xLa)x2=x2+xLa+xLa

end subroutine getX12

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

subroutine getrX12(xa1,xr1,x1,x2)
 !! getrX12 - converst xa,xr coordinates to x1,x2 and maintains periodicity
 use mesh
 implicit none

 real (Long), intent(in)  :: xa1,xr1
 real (Long), intent(out) :: x1,x2

 x1=xa1+0.5_Long*xr1
 x2=xa1-0.5_Long*xr1

! if xx is outside the box, move it in periodically
 if(x1.ge.xLa)x1=x1-xLa-xLa
 if(x1.le.-xLa)x1=x1+xLa+xLa
 if(x2.ge.xLa)x2=x2-xLa-xLa
 if(x2.le.-xLa)x2=x2+xLa+xLa

end subroutine getrX12

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine getK12(ika,ikr,k1,k2)
!! getK12 - converts ixa, ixr grid indices to values of x1,x2 space
 use mesh
 implicit none

 integer,     intent(in)  :: ika, ikr !human-read index of grid in xa,xr space
 real (Long), intent(out) :: k1,k2    !x,x' position in x,x' space

 k1=ka(ika)+0.5_Long*kr(ikr)
 k2=ka(ika)-0.5_Long*kr(ikr)

! k isn't periodic, so don't move periodically - BWB 2010-11-21
! ! if k1,k2 is outside the box, move it in periodically
! if(k1.gt.kLa)k1=k1-kLa-kLa
! if(k1.lt.-kLa)k1=k1+kLa+kLa
! if(k2.gt.kLa)k2=k2-kLa-kLa
! if(k2.lt.-kLa)k2=k2+kLa+kLa

end subroutine getK12

