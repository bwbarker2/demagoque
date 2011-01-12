subroutine fft_initial
 use fftw_constants
 use mesh
 implicit none

! denmat=cmplx(den_re,den_im,8)

 call dfftw_plan_r2r_1d(planwcos,Nxr+1,arraycos,arraycos,FFTW_REDFT00,FFTW_ESTIMATE)
 call dfftw_plan_r2r_1d(planwsin,Nxr-1,arraysin,arraysin,FFTW_RODFT00,FFTW_ESTIMATE)

 call dfftw_plan_dft_r2c_1d(planwx,Nxr+1,arraycos,arraycnum,FFTW_ESTIMATE)

 call dfftw_plan_dft_2d(planf,Nxa,Nxr,denmat,denmat,FFTW_FORWARD,FFTW_MEASURE)
 call dfftw_plan_dft_2d(planb,Nxa,Nxr,denmat,denmat,FFTW_BACKWARD,FFTW_MEASURE)

end subroutine fft_initial

subroutine setState(state)
 !! only implemented for x to wigner and inverse transforms so far
 use mesh
 implicit none

 integer, intent(in) :: state

 if(denState.NE.state) then
  select case (denState)

   case (SPACE)
    select case (state)
     case (WIGNER)
      call transform_x_to_wigner_dumb
     case (MOMENTUM)
      call transform_x_to_wigner_dumb
      call transform_wigner_to_k_dumb
    end select

   case (WIGNER)
    select case (state)
     case (SPACE)
      call transform_wigner_to_x_dumb
     case (MOMENTUM)
      call transform_wigner_to_k_dumb
    end select

   case (MOMENTUM)
    select case (state)
     case (WIGNER)
      call transform_k_to_wigner_dumb
     case (SPACE)
      call transform_k_to_wigner_dumb
      call transform_wigner_to_x_dumb
    end select
  end select
 endif

end subroutine setState

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine transform_x_to_wigner_dumb
 !! brute force method to test the theory
 use mesh
 use fftw_constants
 use phys_cons
 implicit none

 integer :: ixa,ixr,ika,iixa,iixr
 real*8 :: exparg
 complex*16 :: array(-Nxr:Nxr-1)

 do ixa=-Nxa2,Nxa2-1

  ! fill arrays to be transformed
  do ixr=-Nxr,Nxr-1

   if(ixr>=-Nxr2.and.ixr<Nxr2) then
    iixa=ixa
    iixr=ixr
   else
    if(ixa<0) then
     iixa=ixa+Nxa2
    else
     iixa=ixa-Nxa2
    endif
    if(ixr>0) then
     iixr=ixr-Nxr
    else
     iixr=ixr+Nxr
    endif
   endif

   array(ixr)=DBLE(getDenX(iixa,iixr))
!   write(*,*)getDenX(0,-1),getDenX(0,1)
  enddo

  do ika=-Nka2,Nka2-1
   denmat2(ixa,ika)=0d0
   do ixr=-Nxr,Nxr-1
    exparg=delxr*delka*ika*ixr
    denmat2(ixa,ika)=denmat2(ixa,ika)+array(ixr)*exp(-imagi*exparg)
   enddo
   denmat2(ixa,ika)=delxr*denmat2(ixa,ika)
   if(DBLE(denmat2(ixa,ika))>2.d0) write(*,*)ixa,ika,denmat2(ixa,ika)
  enddo

 enddo

 denmat=denmat2

 denState=WIGNER

end subroutine transform_x_to_wigner_dumb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine transform_x_to_wigner
 !! transform_x_to_wigner - transform xr to ka variable. For
 !! derivation, see BWB notes 2010-10-19
 use mesh
 use fftw_constants
 use phys_cons
 implicit none

 integer :: iarray,ixa,iixr,iixa

 !cycle through each ixa column (1D transform)
 do ixa=-Nxa2,Nxa2-1

  arraycos=0.d0
  arraysin=0.d0
  ! fill arrays to be transformed
  do iarray=0,Nxr

   if(iarray<Nxr2) then
    iixa=ixa
    iixr=iarray
   else
    if(ixa<0) then
     iixa=ixa+Nxa2
    else
     iixa=ixa-Nxa2
    endif
    iixr=iarray-Nxr
   endif

   arraycos(iarray)=DBLE(getDenX(iixa,iixr))*cos(0.5*PI*iarray) &
                  -DIMAG(getDenX(iixa,iixr))*sin(0.5*PI*iarray)

   if(iarray.NE.0.AND.iarray.NE.Nxr2) then
    arraysin(iarray-1)=DIMAG(getDenX(iixa,iixr))*cos(0.5*PI*iarray) &
                       +DBLE(getDenX(iixa,iixr))*sin(0.5*PI*iarray)
   endif
  enddo

  ! do sin and cos transform, store in denmat
  call dfftw_execute_dft(planwcos,arraycos,arraycos)
  call dfftw_execute_dft(planwsin,arraysin,arraysin)

  denmat2(ixa,-Nxr2)=delxr*arraycos(0)
  denmat2(ixa,Nxr2)=delxr*arraycos(Nxr)
  do iarray=1,Nxr-1
   denmat2(ixa,iarray-Nxr2)=delxr*(arraycos(iarray)+arraysin(iarray-1))
  enddo

 enddo

 denmat=denmat2

 denState=WIGNER

end subroutine transform_x_to_wigner

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine transform_wigner_to_x_dumb
 ! most straightforward way to compute inverse transform
 use mesh
 use phys_cons
 implicit none

 integer :: ixa,ixr,ika

 denmat2=0d0

 do ixa=-Nxa2,Nxa2-1
  do ixr=-Nxr2,Nxr2-1
   do ika=-Nka2,Nka2-1
    denmat2(ixa,ixr)=denmat2(ixa,ixr)+denmat(ixa,ika)*exp(imagi*delxr*delka*ixr*ika)
   enddo
   denmat2(ixa,ixr)=denmat2(ixa,ixr)*delka/(2d0*pi)
  enddo
 enddo

 denmat=denmat2

 denState=SPACE

end subroutine transform_wigner_to_x_dumb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine transform_wigner_to_x
 !! this is just a standard inverse Fourier transform from ka to xr coordinate.
  ! See notes BWB 2010-11-01 for derivation.
 use fftw_constants
 use mesh
 use phys_cons
 implicit none

 integer :: ixa,iarray,ixr

 do ixa=-Nxa2,Nxa2-1
  arraycos=0.d0
  arraycnum=0.d0

  ! first do the real parts, cos and sin
  do iarray=0,Nxr
   arraycos(iarray)=getDen(ixa,iarray-Nxr2)*cos(pi/2*iarray)
  enddo

  call dfftw_execute_dft_r2c(planwx,arraycos,arraycnum)

  do ixr=0,Nxr2
   arraycnum(ixr)=arraycnum(ixr)*exp(imagi*(-kLa)*(-xLr+ixr*delxr))

   !need to conjugate here because r2c is defaultly signed -1 in FFTW, we
   ! need inverse, so +1. But then we are receiving just the first half of the matrix, which should be for xr<0, and we are storing this in the top (second) half of the matrix, so an additional conjugation should occur, so none overall - BWB, 2010-10-31
   arraycnum(ixr)=arraycnum(ixr)*delka/(2.d0*pi)

   call setDenX(ixa,ixr,arraycnum(ixr))
  enddo
 enddo

 denState=SPACE

end subroutine transform_wigner_to_x

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine transform_wigner_to_k_dumb
 use mesh
 use phys_cons
 implicit none

 integer :: ixa,ika,ikr

 do ika=-Nka2,Nka2-1   !this shouldn't be -1, just for testing BWB 2010-11-17
  do ikr=0,Nkr2
   denmat2(ikr,ika)=0.d0
   do ixa=-Nxa2,Nxa2-1
    denmat2(ikr,ika)=denmat2(ikr,ika)+getDen(ixa,ika)*exp(-imagi*delxa*delkr*ixa*ikr)
   enddo
   denmat2(ikr,ika)=denmat2(ikr,ika)*delxa
  enddo
 enddo

 denmat=denmat2

 denState=MOMENTUM

end subroutine transform_wigner_to_k_dumb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine transform_k_to_wigner_dumb
 use mesh
 use phys_cons
 implicit none

  integer :: ixa,ika,ikr

 do ika=-Nka2,Nka2-1
  do ixa=-Nxa2,Nxa2-1
   denmat2(ixa,ika)=0.d0
   do ikr=-Nkr2,Nka2-1
    denmat2(ixa,ika)=denmat2(ixa,ika)+getDen(ikr,ika)*exp(imagi*delxa*delkr*ixa*ikr)
   enddo
   denmat2(ixa,ika)=denmat2(ixa,ika)*delkr/(2*pi)
  enddo
 enddo

 denmat=denmat2

 denState=WIGNER

end subroutine transform_k_to_wigner_dumb

subroutine FT
 use fftw_constants
 use mesh
 use phys_cons
 use prec_def
 implicit none

 integer :: iixa,iixr,iika,iikr

 IF(denState.NE.1)THEN
  write(6,*)'FT: Error: Tried to FT incorrectly'
  RETURN
 ENDIF

! denmat=cmplx(den_re,den_im,8)

 !multiply by phase factor as in notes from BWB 2010-08-25
 do iixa=0,Nxa-1,2
  do iixr=1,Nxr-1,2
   denmat(iixa,iixr)=denmat(iixa,iixr)*(-1.d0)
  enddo
 enddo

 do iixa=1,Nxa-1,2
  do iixr=0,Nxr-1,2
   denmat(iixa,iixr)=denmat(iixa,iixr)*(-1.d0)
  enddo 
 enddo

 call dfftw_execute_dft(planf,denmat,denmat)

 !multiply by ending phase factor
 do iikr=0,Nkr-1
  do iika=0,Nka-1
   denmat(iikr,iika)=denmat(iikr,iika) &
                     *delxr*exp(imagi*xLr*(ka(-Nka2)+iika*delka)) &
                     *delxa*exp(imagi*xLa*(kr(-Nkr2)+iikr*delkr))
  enddo
 enddo

 denState=2

end subroutine FT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine IFT
 use fftw_constants
 use mesh
 use phys_cons
 use prec_def 
 implicit none

 integer :: iika,iikr, iixa,iixr

 IF(denState.NE.2)THEN
  write(6,*)'IFT: Error: Tried to IFT on wrong matrix'
  RETURN
 ENDIF

! denmat=cmplx(den_re,den_im,8)

 !multiply by phase factor as in notes from BWB 2010-08-25
 do iika=0,Nka-1,2
  do iikr=1,Nkr-1,2
   denmat(iikr,iika)=denmat(iikr,iika)*(-1.d0)
  enddo
 enddo

 do iika=1,Nka-1,2
  do iikr=0,Nkr-1,2
   denmat(iikr,iika)=denmat(iikr,iika)*(-1.d0)
  enddo 
 enddo

 call dfftw_execute_dft(planb,denmat,denmat)

  !multiply by ending phase factor
 do iixa=0,Nxa-1
  do iixr=0,Nxr-1
   denmat(iixa,iixr)=denmat(iixa,iixr)/(4*pi*pi) &
                     *delka*exp(imagi*ka(-Nka2)*(xr(-Nxr2)+iixr*delxr)) &
                     *delkr*exp(imagi*kr(-Nkr2)*(xa(-Nxa2)+iixa*delxa))
  enddo
 enddo


! denmat=denmat/dsqrt(Nxa*Nxr*1.d0)

! den_re=DBLE(denmat)
! den_im=DIMAG(denmat)

 denState=1

end subroutine IFT
