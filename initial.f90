!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculates initial quantities derived from input
SUBROUTINE calcInitial
  use cons_laws
  use fftw_constants
  use mesh
  use osc_pars
  use phys_cons
  use time
  implicit none


 
  call initializeMesh


  ! oscillator data
  w=hbc/m0*6.d0*(rho0/facd)**2/(Nmax+1.d0)!   *0.1
  whm=m0*w/hbc

    write(*,*) 'Parameters of the calculation:'
  write(*,*) 'dxr=',delxr,'dxa=',delxa,'dkr=',delkr,'dka=',delka,'whm=',whm
  
  
!  allocate(tempdenmat(Nxa,Nxr))
  allocate(potx(-Nxa2:Nxa2))
!  allocate(iNxa(Nxa), iNxr(Nxr), iNxa2(-Nxa2:Nxa2), iNxr2(-Nxr2:Nxr2))
 allocate(arraycos(0:Nxr), arraysin(0:Nxr-1), arraycnum(0:Nxr2))


  
end SUBROUTINE calcInitial

subroutine initialState
  use mesh
  use osc_pars
  use prec_def
  implicit none

  real (Long) :: wfnho,y1,y2,xx1,xx2,den0
  integer :: ixa, ixr, iin

  do ixa=-Nxa2,Nxa2-1
     do ixr=-Nxr2,Nxr2-1
        !convert to x,x' representation
        call getX12(ixa,ixr,xx1,xx2)
!        xx1=(xa(ixa)+xr(ixr)/2.d0)
!        xx2=(xa(ixa)-xr(ixr)/2.d0)

        den0=0.0d0
        do iin=0,Nmax
           y1=wfnho(xx1,iin,whm)
           y2=wfnho(xx2,iin,whm)

           den0=den0+y1*y2
        enddo !in

        if(abs(den0).lt.1e-30) den0=0.0d0
        
        call setDenX(ixa,ixr,cmplx(den0,0.d0,8))

     enddo !ixr
  enddo !ixa

  call copyExtra

  ! set density matrix to X state
  denState=SPACE

end subroutine initialState

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine copyExtra
 !! copyExtra - copies matrix to extra part (>Nxr2 and <-Nxr2)
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

end subroutine copyExtra

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine boost
  use mesh
  use phys_cons
  use time
  implicit none

  INTEGER ixr, ixa !loop variables
  real (Long) :: kea !momentum, given ea
  real (Long) :: cos2k, sin2k, epx !cos,sin part of exp, exponent itself
  real (Long) :: xim,xre,xim2,xre2 !x density matrix, imaginary, real
  real (Long) :: x1,x2   !position in x,x' basis

  call setState(SPACE)

  ! non-relativistic conversion
  kea=dsqrt(2.d0*m0*ea)/hbc

  !loop over all grid points
  DO ixa=-Nxa2,Nxa2
     
     DO ixr=-Nxr2,Nxr2
        
        !convert into x,x' coordinates
        call getX12(ixa,ixr,x1,x2)
!        x1=xa(ixa)+0.5*xr(ixr) 
!        x2=xa(ixa)-0.5*xr(ixr) 

        !dispacement operator = exp(-iK(X'-X))
        epx=kea*(x1-x2)
!        udt=m0*w**2*2*xa(ixa)*xr(ixr)*dtim/hbc/2.0_Long
!        udt=0
        cos2k=dcos(epx)
        sin2k=dsin(epx)

        xre=DBLE(getDenX(ixa,ixr))
        xim=DIMAG(getDenX(ixa,ixr))
        
        ! exp(i*edt) = cos2k + i*sin2k
        xre2=xre*cos2k - xim*sin2k
        xim2=xre*sin2k + xim*cos2k

        call setDenX(ixa,ixr,cmplx(xre2,xim2,8))

     ENDDO
     
  ENDDO
 

end subroutine boost

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

 real*8 :: xa1,xr1

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

 real*8, intent(in)  :: xa1,xr1
 real*8, intent(out) :: x1,x2

 x1=xa1+0.5d0*xr1
 x2=xa1-0.5d0*xr1

! if xx is outside the box, move it in periodically
 if(x1.gt.xLa)x1=x1-xLa-xLa
 if(x1.lt.-xLa)x1=x1+xLa+xLa
 if(x2.gt.xLa)x2=x2-xLa-xLa
 if(x2.lt.-xLa)x2=x2+xLa+xLa

end subroutine getrX12

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine getK12(ika,ikr,k1,k2)
!! getK12 - converts ixa, ixr grid indices to values of x1,x2 space and
!maintains periodicity.
 use mesh
 use prec_def
 implicit none

 integer,     intent(in)  :: ika, ikr !human-read index of grid in xa,xr space
 real (Long), intent(out) :: k1,k2    !x,x' position in x,x' space

 k1=ka(ika)+0.5d0*kr(ikr)
 k2=ka(ika)-0.5d0*kr(ikr)

! k isn't periodic, so don't move periodically - BWB 2010-11-21
! ! if k1,k2 is outside the box, move it in periodically
! if(k1.gt.kLa)k1=k1-kLa-kLa
! if(k1.lt.-kLa)k1=k1+kLa+kLa
! if(k2.gt.kLa)k2=k2-kLa-kLa
! if(k2.lt.-kLa)k2=k2+kLa+kLa

end subroutine getK12

