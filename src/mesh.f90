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

MODULE mesh
 use bexception
 use prec_def
 implicit none

 REAL (Long)  :: xLa     ! length of box in xa/2 [fm]
 REAL (Long)  :: xLr     ! length of box in xr/2 [fm]
 real (Long)  :: kLa     ! half length of average momentum box
 INTEGER :: Nxa     ! number of mesh points in xa 
 INTEGER :: Nxr     ! number of mesh points in xr
 INTEGER :: Nxa2    ! Nxa/2
 INTEGER :: Nxr2    ! Nxr/2
 INTEGER :: Nka
 integer :: Nkr
 INTEGER :: Nkr2    ! Nkr/2 = Nxa/2
 INTEGER :: Nka2    ! Nka/2 = Nxr/2

 integer :: Nxam    ! minimum logical index of cell in xa
 integer :: Nxax    ! maximum logical index of cell in xa
 integer :: Nxrm    ! minimum logical index of cell in xr
 integer :: Nxrx    ! maximum logical index of cell in xr
 integer :: Nkam    ! minimum logical index of cell in ka
 integer :: Nkax    ! maximum logical index of cell in ka
 integer :: Nkrm    ! minimum logical index of cell in kr
 integer :: Nkrx    ! maximum logical index of cell in kr

 REAL (Long)    :: delxa   ! interval in x_average
 REAL (Long)    :: delxr   ! interval in x_relative
 REAL (Long)    :: delka   ! in momentum
 REAL (Long)    :: delkr

 real (Long) :: norm_thy  ! theoretical norm (what it 'should' be)

 ! factor to change units to 3D density matrix (calculated in initial.f90)
 REAL (Long) :: facd
 REAL (Long), DIMENSION(:), ALLOCATABLE :: xa,ka,xr,kr   ! coord of grid point  

 ! density matrix, real and imaginary
 REAL (Long) , DIMENSION(:,:) , ALLOCATABLE :: den_re, den_im
 complex (Long), dimension(:,:), allocatable :: denmat !when I need complex, I store here
 complex (Long), dimension(:,:), allocatable :: denmat2 !when I need complex, I store here

 real (Long), dimension(:), allocatable :: denDiagX !< diagonal of spatial density matrix
 real (Long), dimension(:), allocatable :: denDiagK !< diagonal of momentum density matrix

 integer denState  ! gives current coordinate space of density matrix,
                   ! according to the following integer settings:
 integer, parameter :: SPACE = 0
 integer, parameter :: WIGNER = 1
 integer, parameter :: MOMENTUM = 2

 logical :: isDenProcessed !have the diagonals been calculated?

 logical :: isReflectedLR !is the matrix LR reflected?

  ! iNka2 = iNxr2, iNkr2 = iNka2, set here for clarity
 INTEGER , ALLOCATABLE :: iNkr2(:),iNka2(:)
  ! meanfield potential at each (xa,xr=0) point
 real (Long), allocatable :: potDiag(:)

 real (Long) :: maxxim  !maximum imaginary value

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 integer function getNearestIndexX(xx) result(ixx)
  use prec_def
  implicit none

  real (Long), intent(in) :: xx

  real (Long) :: aixx  ! interpolated index for result
  integer, dimension(-Nxa2:Nxa2-1) :: ixen  !array of indices
 
  integer :: ixa,ki
 
  do ixa=-Nxa2,Nxa2-1
   ixen(ixa)=ixa
  enddo

  ki=1
  call lin_int(xa(-Nxa2:Nxa2-1),ixen,Nxa,xx,aixx,ki)

  ixx=nint(aixx)

 end function getNearestIndexX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine initializeMesh
   use input_parameters
   use phys_cons
   use prec_def
   implicit none

   integer :: ixa,ixr,ikr,ika
   real (Long) :: shift   !amount to shift xa,xr,ka,kr definitions
 
   Nxa2=Nxa/2
   Nxr2=Nxr/2
   Nkr2=Nxa2
   Nka2=Nxr2
   Nkr=Nxa
   Nka=Nxr

   ! allocate arrays
   allocate(xa(-Nxa2:Nxa2), kr(-Nkr2:Nkr2), xr(-Nxr:Nxr), ka(-Nka:Nka))
   allocate(denmat(-Nxa2:Nxa2-1,-Nxr:Nxr-1)) !2x size in xr for naive FT - BWB 2011-01-10
   allocate(denmat2(-Nxa2:Nxa2-1,-Nxr:Nxr-1))
   allocate(denDiagX(-Nxa2:Nxa2-1))
   allocate(denDiagK(-Nka:Nka-1))
   allocate(potDiag(-Nxa2:Nxa2))
   allocate(den_re(-Nxa2:Nxa2-1,-Nxr:Nxr-1))

   !facd calc'd here because can't initialize with non-integer exponents
   facd=sqrt(5e0_Long/3e0_Long)*(deg*pi*(rho0**2)/6e0_Long)**(1e0_Long/3e0_Long) 


   ! mesh in x and k
   delxa=(2e0_Long*xLa)/real(Nxa)
   delxr=(2e0_Long*xLr)/real(Nxr)
   delka=pi/(2e0_Long*xLr)
   delkr=pi/xLa

   potDiag=0e0_Long

   Nxam=-Nxa2
   Nxax=Nxa2-1
   Nxrm=-Nxr2
   Nxrx=Nxr2-1
   Nkrm=-Nkr2
   Nkrx=Nkr2-1
   Nkam=-Nka2
   Nkax=Nka2-1

   if(useMeshShifted) then
    shift=0.5_Long
   else
    shift=0e0_Long
   endif

   do ixa=-Nxa2,Nxa2
    xa(ixa)=delxa*(ixa+shift)
   enddo

   do ixr=-Nxr,Nxr
    xr(ixr)=delxr*(ixr+shift)
   enddo

   do ikr=-Nkr2,Nkr2
    kr(ikr)=delkr*(ikr+shift)
   enddo

   do ika=-Nka,Nka
    ka(ika)=delka*(ika+shift)
   enddo

   kLa=-ka(-Nka)+delka*shift
!   write(*,*)'kLa=',kLa

   isReflectedLR=.false.
   isDenProcessed=.false.

  end subroutine initializeMesh

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex (Long) function getDen(i1,i2)
   !! getDen - returns value of density matrix, given the current denState
   use prec_def
   implicit none

   integer, intent(in) :: i1,i2

   select case (denState)
    case (SPACE)
     getDen=getDenX(i1,i2)
    case (WIGNER)
     getDen=getDenW(i1,i2)
    case (MOMENTUM)
     getDen=getDenK(i1,i2)
    case default
     call throwException('getDen: mesh in invalid state, should never happen',BEXCEPTION_FATAL)
     getDen=999 !this will never run, since fatal exception is called above.
   end select

   end function getDen

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> gets diagonal of x-space density matrix
  real (Long) function getDenDiagX(ixa) result(den)
   use input_parameters
   implicit none

   integer, intent(in) :: ixa

   if(.not.useMeshShifted) then
    den=abs(getDenX(ixa,0))
   else
    if(.not.isDenProcessed) then
     call mesh_processDen
    endif
    den=denDiagX(ixa)
   endif

  end function getDenDiagX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> gets diagonal of k-space density matrix. 
  real (Long) function getDenDiagK(ika) result(den)
   use input_parameters
   use prec_def
   implicit none

   integer, intent(in) :: ika

   if(.not.useMeshShifted) then
    den=abs(getDenK(0,ika))
   else
    if(.not.isDenProcessed) then
     call mesh_processDen
    endif
    den=denDiagK(ika)
   endif

  end function getDenDiagK

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex (Long) function getDenX(ixa,ixr)
   !! getDenX - returns value of spatial density matrix at index (ixa,ixr)
   implicit none

   integer, intent(in) :: ixa,ixr

   getDenX=denmat(ixa,ixr)

  end function getDenX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine mesh_processDen
  use phys_cons
  implicit none

  if(denState.ne.WIGNER) then
   call throwException( &
    'getDenDiagX: density matrix has been changed since last processed.' &
    //'Processing now. Suggest re-ordering to avoid extra FFTs.' &
    ,BEXCEPTION_WARNING)
   call setState(WIGNER)
  endif

  ! sum the momentum components to get diagonal in space
  denDiagX=delka*invsqrt2pi*real(sum(denmat(:,0:Nka2-1),2))

  ! sum the space components to get diagonal in momentum
  denDiagK=delxa*invsqrt2pi*real(sum(denmat(0:Nxa2-1,:),1))
  
  isDenProcessed=.true.
 
 end subroutine mesh_processDen

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !> reflects the matrix about the xa=0 axis, exchanging the left and right sides
 !!
 !! \warning does not exchange x and x' or k and k'! (need to reflect up/down as well)
 subroutine mesh_reflectLR()
  use prec_def
  implicit none

  complex (Long) :: wnum

  integer :: ixa,ixr

   do ixa=-Nxa2,Nxa2-1
    do ixr=0,Nxr-1
     wnum=denmat(ixa,ixr)
     denmat(ixa,ixr)=denmat(-ixa,-ixr)
     denmat(-ixa,-ixr)=wnum
    enddo
   enddo

   do ixa=1,Nxa2-1
    wnum=denmat(ixa,-Nxr)
    denmat(ixa,-Nxr)=denmat(-ixa,-Nxr)
    denmat(-ixa,-Nxr)=wnum
   enddo

   do ixr=1,Nxr-1
    wnum=denmat(-Nxa2,ixr)
    denmat(-Nxa2,ixr)=denmat(-Nxa2,-ixr)
    denmat(-Nxa2,-ixr)=wnum
   enddo

   isReflectedLR=.not.isReflectedLR

 end subroutine mesh_reflectLR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine mesh_setReflectedLR(reflect)
  implicit none

  logical, intent(in) :: reflect

  if(reflect.neqv.isReflectedLR) then
   call mesh_reflectLR
  endif

 end subroutine mesh_setReflectedLR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> sets the value of the spatial density matrix at index (ixa,ixr)
  !!
  !! \note Only sets elements that are actually used in getDenX
  subroutine setDenX(ixa,ixr,value)

   integer, intent(in)    :: ixa,ixr
   complex (Long), intent(in) :: value

   denmat(ixa,ixr)=value

   if(isDenProcessed) isDenProcessed=.false.

  end subroutine setDenX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> returns value of Wigner density matrix at index (ixa,ika)
  complex (Long) function getDenW(ixa,ika)
   implicit none

   integer, intent(in) :: ixa,ika

   if(ixa<Nxa2) then
    getDenW=denmat(ixa,ika)
   else
    getDenW=denmat(-Nxa2,ika)
   endif


  end function getDenW

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine setDenW(ixa,ika, this_value)
   implicit none

   integer, intent(in) :: ixa,ika
   complex (Long), intent(in) :: this_value

   denmat(ixa,ika)=this_value

   if(isDenProcessed) isDenProcessed=.false.

  end subroutine setDenW

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> returns value of spectral density matrix at index (ikr,ika)
  complex (Long) function getDenK(ikr,ika)
   implicit none

   integer, intent(in) :: ikr,ika

   getDenK=denmat(ikr,ika)

  end function getDenK

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> sets value of the spectral density matrix at index (ikr,ika)
  subroutine setDenK(ikr,ika,val)
   implicit none

   complex (Long), intent(in) :: val
   integer, intent(in)    :: ikr, ika

   denmat(ikr,ika)=val

   if(isDenProcessed) isDenProcessed=.false.

  end subroutine setDenK

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !> finds eigenvalues, eigenvectors of spatial density matrix.
   !!
   !! \note This currently only works for Nxr=Nxa (square matrices)
  subroutine getDenEigens(evals,evecs)
   use lib_lapack
   implicit none

   complex (Long), dimension(0:Nxa-1), intent(out) :: evals
   complex (Long), dimension(-Nxa2:Nxa2-1,-Nxr2:Nxr2-1), intent(out) :: evecs

   integer :: ixa,ixr

   if(Nxa==Nxr) then
    continue
   else
    write(*,*)'getEigens only works for Nxa=Nxr, exiting subroutine'
    return
   endif

   call setState(SPACE)

   do ixa=-Nxa2,Nxa2-1
    do ixr=-Nxr2,Nxr2-1
     evecs(ixa,ixr)=getDenX(ixa,ixr)
    enddo
   enddo

   call getEigenSq(evecs,Nxa,evals,evecs)

  end subroutine getDenEigens

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Sets Fourier state of system to desired state (position, wigner, or
  !! momentum. If already in desired state, does nothing.
  subroutine setState(state)
   implicit none
   
   integer, intent(in) :: state !< state to set system to

   real :: totalelapsed, elapsed

   if(denState.NE.state) then

   call cpu_time(totalelapsed)

!call mesh_setReflectedLR(.true.)

    select case (denState)
  
     case (SPACE)
      select case (state)
       case (WIGNER)
!call mesh_setReflectedLR(.true.)
!        call transform_x_to_w_norepeat_fft
!call mesh_setReflectedLR(.false.)
        call transform_x_to_wigner_dumb
!        call transform_x_to_w_dumb_kshift
        call cpu_time(elapsed)
        write(*,*)'transform_x_to_w:',elapsed-totalelapsed,'seconds'
       case (MOMENTUM)
!        call transform_x_to_k_norepeat
!call mesh_setReflectedLR(.true.)
!        call transform_x_to_w_norepeat_fft
!call mesh_setReflectedLR(.false.)
        call transform_x_to_wigner_dumb
!call mesh_setReflectedLR(.true.)
!        call transform_wigner_to_k_fft_exp
!call mesh_setReflectedLR(.false.)
        call transform_wigner_to_k_dumb
!        call transform_w_to_k_norepeat
!        write(*,*)'transform_x_to_k:',etime(elapsed)-totalelapsed,'seconds'
      end select
  
     case (WIGNER)
      select case (state)
       case (SPACE)
!        call transform_w_to_x_norepeat_fft
        call transform_wigner_to_x_dumb
        call cpu_time(elapsed)
        write(*,*)'transform_w_to_x:',elapsed-totalelapsed,'seconds'
       case (MOMENTUM)
        call transform_wigner_to_k_dumb
!        call transform_w_to_k_norepeat
!call mesh_setReflectedLR(.true.)
!        call transform_wigner_to_k_fft_exp
!call mesh_setReflectedLR(.false.)
        call cpu_time(elapsed)
        write(*,*)'transform_w_to_k:',elapsed-totalelapsed,'seconds'
      end select
  
     case (MOMENTUM)
      select case (state)
       case (WIGNER)
!        call transform_k_to_w_fft_norepeat
!        call transform_k_to_wigner_fft_exp
        call transform_k_to_wigner_dumb
        call cpu_time(elapsed)
        write(*,*)'transform_k_to_w:',elapsed-totalelapsed,'seconds'
       case (SPACE)
!        call transform_k_to_w_fft_norepeat
!        call transform_k_to_wigner_fft_exp
        call transform_k_to_wigner_dumb
!        write(*,*)'transform_k_to_w__:',etime(elapsed)-totalelapsed,'seconds'
!        call transform_w_to_x_norepeat_fft
        call transform_wigner_to_x_dumb
!        write(*,*)'transform_k_to_x:',etime(elapsed)-totalelapsed,'seconds'
      end select
    end select
!call mesh_setReflectedLR(.false.)
   endif
  
  end subroutine setState
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine transform_x_to_wigner_trig
   !! brute force method to test the theory, using cos/sin transforms. Theory in notes, BWB 2011-02-25p2
   use phys_cons
   implicit none
  
   integer :: ixa,ixr,ika
   real (Long) :: trigarg
   real (Long) :: wigden(-Nxa2:Nxa2-1,-Nka:Nka-1)
   complex (Long) :: array(0:Nxr)
  
   do ixa=-Nxa2,Nxa2-1
  
    array=0e0_Long
  
    ! fill arrays to be transformed
    do ixr=0,Nxr-1
     array(ixr)=getDenX(ixa,ixr)
    enddo

    array(Nxr)=getDenX(ixa,-Nxr)

    do ika=-Nka,Nka-1
     wigden(ixa,ika)=0e0_Long
     do ixr=1,Nxr-1
      trigarg=ka(ika)*xr(ixr)
      wigden(ixa,ika)=wigden(ixa,ika)+REAL(array(ixr))*cos(trigarg) &
                                    +AIMAG(array(ixr))*sin(trigarg)
     enddo
     wigden(ixa,ika)=wigden(ixa,ika)*2e0_Long
     wigden(ixa,ika)=wigden(ixa,ika)+real(array(0))+real(array(Nxr))*(-1)**ika
     wigden(ixa,ika)=wigden(ixa,ika)*delxr*invsqrt2pi
  
     ! if the cell is unreasonably large, write out
!     if(REAL(wigden(ixa,ika))>2.d0) write(*,*)ixa,ika,wigden(ixa,ika)
    enddo
  
   enddo
  
   do ixa=-Nxa2,Nxa2-1
    do ika=-Nka,Nka-1
     call setDenW(ixa,ika,CMPLX(wigden(ixa,ika),0e0_Long,Long))
    enddo
   enddo
  
   denState=WIGNER
  
  end subroutine transform_x_to_wigner_trig
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine transform_x_to_wigner_dumb
   !! brute force method to test the theory
   use phys_cons
   implicit none
  
   integer :: ixa,ixr,ika
   real (Long) :: exparg
   complex (Long) :: array(-Nxr:Nxr-1)

!write(*,*)'x to w dumb go'
  
   do ixa=-Nxa2,Nxa2-1
  
    array=cmplx(0e0_Long,0e0_Long,Long)
  
    ! fill arrays to be transformed
    do ixr=-Nxr,Nxr-1
     array(ixr)=getDenX(ixa,ixr)
    enddo
  
    do ika=-Nka,Nka-1
!     denmat2(ixa,ika)=cmplx(0d0,0d0,8)
     denmat2(ixa,ika)=czero
     do ixr=1,Nxr-1

!      !debug up-down symmetry
!      if(ixr.ne.-Nxr)then
!       if(array(ixr).ne.array(-ixr)) then
!        write(*,*)'ixa,ixr,diff',ixa,ixr,array(ixr)-array(-ixr)
!       endif
!      endif

      exparg=xr(ixr)*ka(ika)
!      exparg=mod(exparg,2*pi)  !this doesn't seem to have any effect
!      denmat2(ixa,ika)=denmat2(ixa,ika)+array(ixr)*exp(-imagi*exparg)
      denmat2(ixa,ika)=denmat2(ixa,ika) &
                       +array(ixr)*(cos(exparg)-imagi*sin(exparg)) &
                       +array(-ixr)*(cos(exparg)+imagi*sin(exparg))
     enddo
     denmat2(ixa,ika)=denmat2(ixa,ika) &
                      +array(-Nxr)*cos(xr(-Nxr)*ka(ika)) &
                      +array(0)

     denmat2(ixa,ika)=delxr*denmat2(ixa,ika)*invsqrt2pi
  
     ! if the cell is unreasonably large, write out
!     if(REAL(denmat2(ixa,ika))>2.d0) write(*,*)ixa,ika,denmat2(ixa,ika)
    enddo
  
   enddo
  
   do ixa=-Nxa2,Nxa2-1
    do ika=-Nka,Nka-1
     call setDenW(ixa,ika,denmat2(ixa,ika))
    enddo
   enddo
  
   denState=WIGNER
  
  end subroutine transform_x_to_wigner_dumb
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine transform_x_to_w_dumb_kshift
   ! transform x to w using explicitly real result, as per notes BWB 2011-08-29
   use phys_cons
   use prec_def
   implicit none

   integer ixa,ixr,ika

   denmat2=0e0_Long

   do ixa=-Nxa2,Nxa2-1
    do ika=-Nka,Nka-1

     denmat2(ixa,ika)=denmat2(ixa,ika) &
                      +REAL(getDen(ixa,-Nxr))*COS(xr(-Nxr)*ka(ika)) &
                      +getDen(ixa,0)

     do ixr=1,Nxr-1

      denmat2(ixa,ika)=denmat2(ixa,ika) &
                       +2e0_Long*(REAL(getDen(ixa,ixr))*cos(xr(ixr)*ka(ika)) &
                                +AIMAG(getDen(ixa,ixr))*sin(xr(ixr)*ka(ika)))

     enddo !ixr

    enddo !ika
   enddo !ixa

   denmat=denmat2*delxr*invsqrt2pi

   denState=WIGNER

  end subroutine transform_x_to_w_dumb_kshift

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine transform_w_to_x_norepeat_fft
   !! implements simple FFT, as per notes BWB 2011-03-28p2
   use lib_fftw
   use phys_cons
   implicit none

   integer :: ixa, ika, ixr, sgnfac
   complex (Long) :: array(-Nxr:Nxr-1)

   !only compute half of it
   do ixa=-Nxa2,-1
    sgnfac=1
    do ika=-Nka,Nka-1
     array(ika)=denmat(ixa,ika)*sgnfac
     sgnfac=-sgnfac
    enddo !ika

    call ift_z2z_1d(array,array,2*Nxr)

    do ixr=-Nxr+1,Nxr-1,2
     array(ixr)=array(ixr)*(-1e0_Long)
    enddo

    do ixr=-Nxr,Nxr-1
     denmat(ixa,ixr)=array(ixr)*delka*invsqrt2pi
    enddo

   enddo !ixa

   !now copy to second half of matrix
   do ixa=0,Nxa2-1
    do ixr=-Nxr,-1
     denmat(ixa,ixr)=denmat(ixa-Nxa2,ixr+Nxr)
    enddo !ixr
    do ixr=0,Nxr-1
     denmat(ixa,ixr)=denmat(ixa-Nxa2,ixr-Nxr)
    enddo !ixr
   enddo !ixa

   denState=SPACE

  end subroutine transform_w_to_x_norepeat_fft

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !> for derivation, see BWB notes 2011-03-25. Essentially this uses 2 sine transforms and 2 cosine transforms, needed because the exponent in the transform has pi/N instead of 2*pi/N.
   ! \warning This is unfinished, because it is much too complicated. See notes BWB 2011-03-28p2
  subroutine transform_w_to_x_norepeat_fft_bad
   use lib_fftw
   use phys_cons
   implicit none

   integer :: ixa,ika,ixr
   real (Long)  :: array(0:Nka)

   denmat2=0e0_Long

   do ixa=0,Nxa2-1

    !create first input array
    array=0e0_Long
    do ika=1,Nka-1
     array(ika)=0.5_Long*(REAL(denmat(ixa,ika))+REAL(denmat(ixa,-ika)))
    enddo !ika
    array(0)=REAL(denmat(ixa,0))
    array(Nka)=REAL(denmat(ixa,-Nka))

    call ft_re_1d(array,array,Nka+1)

    do ixr=0,Nxr-1
     denmat2(ixa,ixr)=denmat2(ixa,ixr)+array(ixr)
    enddo
    denmat2(ixa,-Nxr)=denmat2(ixa,-Nxr)+array(ixr)

    !second input array
    array=0e0_Long
    do ika=1,Nka-1
     array(ika-1)=0.5_Long*(-AIMAG(denmat(ixa,ika))+AIMAG(denmat(ixa,-ika)))
    enddo !ika

    call ft_ro_1d(array,array,Nka-1)

    do ixr=1,Nxr-1
     denmat2(ixa,ixr)=denmat2(ixa,ixr)+array(ixr-1) &
                                      -AIMAG(denmat(ixa,-Nka))*(-1)**ixr &
                                      -AIMAG(denmat(ixa,0))
    enddo
    denmat2(ixa,-Nxr)=denmat2(ixa,-Nxr)+array(Nxr) &
                                       -AIMAG(denmat(ixa,-Nka)) &
                                       -AIMAG(denmat(ixa,0))
 
    !third input array
    array=0e0_Long
    do ika=1,Nka-1
     array(ika)=0.5e0_Long*(AIMAG(denmat(ixa,ika))+AIMAG(denmat(ixa,-ika)))
    enddo !ika
    array(0)=AIMAG(denmat(ixa,0))
    array(Nka)=AIMAG(denmat(ixa,-Nka))

    call ft_re_1d(array,array,Nka+1)

    do ixr=0,Nxr-1
     denmat2(ixa,ixr)=denmat2(ixa,ixr)+imagi*array(ixr)
    enddo
    denmat2(ixa,-Nxr)=denmat2(ixa,-Nxr)-imagi*array(ixr) !minus for conjugation

!Now do the same thing for the other 3 terms.

    array=0e0_Long
    do ika=1,Nka-1
     array(ika)=0.5_Long*(REAL(denmat(ixa,ika))-REAL(denmat(ixa,-ika)))
    enddo !ika
    array(0)=AIMAG(denmat(ixa,0))
    array(Nka)=AIMAG(denmat(ixa,-Nka))

    call ft_ro_1d(array,array,Nka+1)

    do ixr=0,Nxr-1
     denmat2(ixa,ixr)=denmat2(ixa,ixr)+array(ixr) &
                                      -AIMAG(denmat(ixa,-Nka))*(-1)**ixr &
                                      -AIMAG(denmat(ixa,0))
    enddo
    denmat2(ixa,-Nxr)=denmat2(ixa,-Nxr)+array(Nxr) &
                                       -AIMAG(denmat(ixa,-Nka)) &
                                       -AIMAG(denmat(ixa,0))
 

   enddo !ixa

  end subroutine transform_w_to_x_norepeat_fft_bad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> most straightforward way to compute inverse transform
  subroutine transform_wigner_to_x_trig
   use phys_cons
   implicit none

   integer :: ixa,ixr,ika

   denmat2=CMPLX(0e0_Long,0e0_Long,Long)

   do ixa=-Nxa2,Nxa2-1
    do ixr=-Nxr,Nxr-1
     do ika=-Nka,Nka-1
      denmat2(ixa,ixr)=denmat2(ixa,ixr)+real(getDenW(ixa,ika))*exp(imagi*xr(ixr)*ka(ika))
     enddo
     denmat2(ixa,ixr)=denmat2(ixa,ixr)*delka*invsqrt2pi
    enddo
   enddo
  
   denmat=denmat2
  
   denState=SPACE
  
  end subroutine transform_wigner_to_x_trig
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> most straightforward way to compute inverse transform
  subroutine transform_wigner_to_x_dumb
   use phys_cons
   implicit none

   real (Long) :: exparg
   integer :: ixa,ixr,ika

   denmat2=CMPLX(0e0_Long,0e0_Long,Long)

   do ixa=-Nxa2,Nxa2-1
    do ixr=-Nxr,Nxr-1
!     do ika=-Nka,Nka-1
!      denmat2(ixa,ixr)=denmat2(ixa,ixr)+denmat(ixa,ika)*exp(imagi*delxr*delka*ixr*ika)
!     enddo
     do ika=1,Nka-1

      exparg=xr(ixr)*ka(ika)
      denmat2(ixa,ixr)=denmat2(ixa,ixr) &
                       +denmat(ixa,ika)*(cos(exparg)+imagi*sin(exparg)) &
                       +denmat(ixa,-ika)*(cos(exparg)-imagi*sin(exparg))
     enddo
     denmat2(ixa,ixr)=denmat2(ixa,ixr) &
                      +denmat(ixa,-Nka)*cos(xr(ixr)*ka(-Nka)) &
                      +denmat(ixa,0)

     denmat2(ixa,ixr)=denmat2(ixa,ixr)*delka*invsqrt2pi
    enddo
   enddo
  
   denmat=denmat2
  
   denState=SPACE
  
  end subroutine transform_wigner_to_x_dumb
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !> brute force method to test the theory, using cos/sin transforms. Theory in notes, BWB 2011-02-28p1
  subroutine transform_k_to_wigner_trig
   use phys_cons
   implicit none
  
   integer :: ixa,ika,ikr
   real (Long) :: trigarg
   real (Long) :: wigden(-Nxa2:Nxa2-1,-Nka:Nka-1)
   complex (Long) :: array(0:Nkr2)
  
   do ika=-Nka,Nka-1
  
    array=0e0_Long
  
    ! fill arrays to be transformed
    do ikr=0,Nkr2-1
     array(ikr)=getDenK(ikr,ika)
    enddo

    array(Nkr2)=getDenK(-Nkr2,ika)

    do ixa=-Nxa2,Nxa2-1
     wigden(ixa,ika)=0e0_Long
     do ikr=1,Nkr2-1
      trigarg=xa(ixa)*kr(ikr)
      wigden(ixa,ika)=wigden(ixa,ika)+REAL(array(ikr))*cos(trigarg) &
                                    -AIMAG(array(ikr))*sin(trigarg)
     enddo
     wigden(ixa,ika)=wigden(ixa,ika)*2e0_Long
     wigden(ixa,ika)=wigden(ixa,ika)+REAL(array(0))+REAL(array(Nxr))*(-1)**ixa
     wigden(ixa,ika)=wigden(ixa,ika)*delkr*invsqrt2pi
  
     ! if the cell is unreasonably large, write out
!     if(REAL(wigden(ixa,ika))>2.d0) write(*,*)ixa,ika,wigden(ixa,ika)
    enddo
  
   enddo
  
   do ixa=-Nxa2,Nxa2-1
    do ika=-Nka,Nka-1
     call setDenW(ixa,ika,CMPLX(wigden(ixa,ika),0e0_Long,Long))
    enddo
   enddo
  
   denState=WIGNER
  
  end subroutine transform_k_to_wigner_trig
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine transform_wigner_to_k_trig
   use phys_cons
   implicit none

   integer :: ixa,ika,ikr

   do ika=-Nka,Nka-1
    do ikr=-Nkr2,Nkr2-1
     denmat2(ikr,ika)=0e0_Long
     do ixa=-Nxa2,Nxa2-1
      denmat2(ikr,ika)=denmat2(ikr,ika)+real(getDenW(ixa,ika))*exp(-imagi*xa(ixa)*kr(ikr))
     enddo
     denmat2(ikr,ika)=denmat2(ikr,ika)*delxa*invsqrt2pi
    enddo
   enddo
  
   denmat=denmat2
  
   denState=MOMENTUM
  
  end subroutine transform_wigner_to_k_trig
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine transform_wigner_to_k_dumb
   use phys_cons
   use prec_def
   implicit none

   integer :: ixa,ika,ikr
   real(Long) :: trigarg

   do ika=-Nka,Nka-1
    do ikr=-Nkr2,Nkr2-1
     denmat2(ikr,ika)=czero
!     do ixa=-Nxa2,Nxa2-1
     do ixa=1,Nxa2-1
      trigarg=xa(ixa)*kr(ikr)
!      denmat2(ikr,ika)=denmat2(ikr,ika)+getDen(ixa,ika)*exp(-imagi*delxa*delkr*ixa*ikr)
      denmat2(ikr,ika)=denmat2(ikr,ika) &
                       +denmat(ixa,ika)*(cos(trigarg) &
                                         -imagi*sin(trigarg)) &
                       +denmat(-ixa,ika)*(cos(trigarg) &
                                         +imagi*sin(trigarg))
     enddo !ixa
     denmat2(ikr,ika)=denmat2(ikr,ika) &
                      +denmat(-Nxa2,ika)*cos(xa(-Nxa2)*kr(ikr)) &
!                                          -imagi*sin(xa(-Nxa2)*kr(ikr))) &
                      +denmat(0,ika)

     denmat2(ikr,ika)=denmat2(ikr,ika)*delxa*invsqrt2pi

     !Analytically (see notes BWB 2011-01-25p2), for ika+ikr odd, denmat should
     ! be exactly zero. Make it so. This subroutine is not built for speed,
     ! so leave the above computation alone, but add the following:
     if(abs(mod(ika+ikr,2))==1)denmat2(ikr,ika)=0e0_Long

    enddo
   enddo
  
   denmat=denmat2
  
   denState=MOMENTUM
  
  end subroutine transform_wigner_to_k_dumb
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine transform_wigner_to_k_fft_exp
   use lib_fftw
   use phys_cons
   implicit none

   integer :: ixa,ika,ikr
   complex (Long), dimension(-Nxa2:Nxa2-1) :: arrin,arrout

   arrin=cmplx(0e0_Long,0e0_Long,Long)
   arrout=cmplx(0e0_Long,0e0_Long,Long)

   !scan over ika
   do ika=-Nka,Nka-1

    !construct array to transform
    do ixa=-Nxa2,Nxa2-1

     arrin(ixa)=getDenW(ixa,ika)

     !shift indices of transform by multiplying every other value by -1
     ! (see paper notes BWB 2010-08-25)
     arrin(ixa)=arrin(ixa)*(-1)**(ixa+Nxa2)
    enddo

    call ft_z2z_1d(arrin,arrout,Nxa)

    denState=MOMENTUM

    !shift k indices, multiply coefficients, write to density matrix
    do ikr=-Nkr2,Nkr2-1
     arrout(ikr)=arrout(ikr)*(-1)**(ikr+Nkr2)
     arrout(ikr)=arrout(ikr)*delxa*invsqrt2pi
     call setDenK(ikr,ika,arrout(ikr))
    enddo

   enddo

  end subroutine transform_wigner_to_k_fft_exp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine transform_k_to_wigner_dumb
   use phys_cons
   use prec_def
   implicit none
  
   integer :: ixa,ika,ikr
   real(Long) :: trigarg
  
   do ika=-Nka,Nka-1
    do ixa=-Nxa2,Nxa2-1
     denmat2(ixa,ika)=czero
      do ikr=1,Nkr2-1
       trigarg=xa(ixa)*kr(ikr)
!      denmat2(ikr,ika)=denmat2(ikr,ika)+getDen(ixa,ika)*exp(-imagi*delxa*delkr*ixa*ikr)
      denmat2(ixa,ika)=denmat2(ixa,ika) &
                       +denmat(ikr,ika)*(cos(trigarg)+imagi*sin(trigarg)) &
!                       +REAL(denmat(ikr,ika))*cos(trigarg) &
!                       -Dimag(denmat(ikr,ika))*sin(trigarg)
                       +denmat(-ikr,ika)*(cos(trigarg) &
                                         -imagi*sin(trigarg))
     enddo !ixa
     denmat2(ixa,ika)=denmat2(ixa,ika) &
                      +denmat(-Nkr2,ika)*cos(xa(ixa)*kr(-Nkr2)) &
                      +denmat(0,ika)

!     denmat2(ixa,ika)=denmat2(ixa,ika)+getDen(ikr,ika)*exp(imagi*delxa*delkr*ixa*ikr)
!     enddo
!      denmat2(ixa,ika)=denmat2(ixa,ika)+getDen(ikr,ika)*
     denmat2(ixa,ika)=denmat2(ixa,ika)*delkr*invsqrt2pi
    enddo
   enddo
  
   denmat=denmat2
  
   denState=WIGNER
  
  end subroutine transform_k_to_wigner_dumb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine transform_k_to_wigner_fft_exp
   use lib_fftw
   use phys_cons
   implicit none

   integer :: ixa,ika,ikr
   complex (Long), dimension(-Nxa2:Nxa2-1) :: arrin,arrout

   arrin=cmplx(0e0_Long,0e0_Long,Long)
   arrout=cmplx(0e0_Long,0e0_Long,Long)

   !scan over ika
   do ika=-Nka,Nka-1

    !construct array to transform
    do ikr=-Nkr2,Nkr2-1

     arrin(ikr)=getDenK(ikr,ika)

     !shift indices of transform by multiplying every other value by -1
     ! (see paper notes BWB 2010-08-25)
     arrin(ikr)=arrin(ikr)*(-1)**(ikr+Nkr2)
    enddo

    call ift_z2z_1d(arrin,arrout,Nkr)

    denState=WIGNER

    !shift k indices, multiply coefficients, write to density matrix
    do ixa=-Nxa2,Nxa2-1
     arrout(ixa)=arrout(ixa)*(-1)**(ixa+Nxa2)
     arrout(ixa)=arrout(ixa)*delkr*invsqrt2pi
     call setDenW(ixa,ika,arrout(ixa))
    enddo

   enddo

  end subroutine transform_k_to_wigner_fft_exp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine transform_k_to_w_fft_norepeat
   use phys_cons
   use prec_def
   implicit none
  
   integer :: ixa,ika,ikr
   real(Long) :: trigarg
   real (Long), allocatable, dimension(:),save :: arraycos, arraysin
   integer*8, save :: plan_cos, plan_sin  !< plans for FFTW

   include '/usr/include/fftw3.f'

   !do only odd ika, even ika is done below
   do ika=-Nka+1,Nka-1,2
    do ixa=-Nxa2,Nxa2-1
     denmat2(ixa,ika)=0e0_Long
      do ikr=1,Nkr2-1
       trigarg=xa(ixa)*kr(ikr)
!      denmat2(ikr,ika)=denmat2(ikr,ika)+getDen(ixa,ika)*exp(-imagi*delxa*delkr*ixa*ikr)
      denmat2(ixa,ika)=denmat2(ixa,ika) &
                       +denmat(ikr,ika)*(cos(trigarg)+imagi*sin(trigarg)) &
!                       +REAL(denmat(ikr,ika))*cos(trigarg) &
!                       -Dimag(denmat(ikr,ika))*sin(trigarg)
                       +denmat(-ikr,ika)*(cos(trigarg) &
                                         -imagi*sin(trigarg))
     enddo !ixa
     denmat2(ixa,ika)=denmat2(ixa,ika) &
                      +real(denmat(-Nkr2,ika),Long)*cos(xa(ixa)*kr(-Nkr2)) &
                      +denmat(0,ika)

!     denmat2(ixa,ika)=denmat2(ixa,ika)+getDen(ikr,ika)*exp(imagi*delxa*delkr*ixa*ikr)
!     enddo
!      denmat2(ixa,ika)=denmat2(ixa,ika)+getDen(ikr,ika)*
     denmat2(ixa,ika)=denmat2(ixa,ika)*delkr*invsqrt2pi
    enddo
   enddo
  
   denmat=denmat2

   if(.not.allocated(arraycos)) then
!    allocate(arraycos(0:Nkr2/2),arraysin(0:Nkr2/2))
    allocate(arraycos(0:Nkr2/2),arraysin(0:Nkr2/2))

    call dfftw_plan_r2r_1d(plan_cos,size(arraycos),arraycos,arraycos &
                           ,FFTW_REDFT00,FFTW_ESTIMATE)
    call dfftw_plan_r2r_1d(plan_sin,Nkr2/2-1,arraysin(1:Nkr2/2-1),arraysin(1:Nkr2/2-1) &
                           ,FFTW_RODFT00,FFTW_ESTIMATE)

    !these will always be zero, but they are so handy to have 
    !for whole-array operations below
    arraysin(0)=0e0_Long  
    arraysin(Nkr2/2)=0e0_Long

!    arraycos=0

   endif !not allocated

   !now do even ika rows
   do ika=-Nka,Nka-2,2
     arraycos(0)=real(denmat(0,ika),Long)
    do ikr=1,Nkr2/2-1
     arraycos(ikr)=real(denmat(ikr*2,ika),Long)
     arraysin(ikr)=aimag(denmat(ikr*2,ika))
    enddo
    arraycos(Nkr2/2)=aimag(denmat(-Nkr2,ika))

    call dfftw_execute(plan_cos)
    call dfftw_execute(plan_sin)

    arraycos=arraycos*delkr*invsqrt2pi
    arraysin=arraysin*delkr*invsqrt2pi

    ! fill density matrix, starting with 3rd quarter of data (0 to N/4)
    denmat(0:Nxa2/2,ika)=arraycos-arraysin

    ! now do 2nd quarter (-N/4 to 0)
    do ixa=-Nxa2/2,-1
     denmat(ixa,ika)=arraycos(-ixa)+arraysin(-ixa)
    enddo !ixa

    ! copy from 2nd quarter to 4th quarter of matrix
    denmat(Nxa2/2+1:Nxa2-1,ika)=denmat(-Nxa2/2+1:-1,ika)

    ! copy from 3rd quarter to 1st quarter (-N/2 to -N/4)
    denmat(-Nxa2:-Nxa2/2-1,ika)=denmat(0:Nxa2/2-1,ika)

   enddo !ika

   denState=WIGNER
  
  end subroutine transform_k_to_w_fft_norepeat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine transform_x_to_k_norepeat
   !! transform_x_to_k_nozeroes - transforms from coordinate to momentum space, with no redundancy from periodic extensions. See paper notes BWB 2011-03-11.
   use phys_cons
   implicit none

   integer    :: ikr,ikr2,ika,ika2,ixa,ixr
   real (Long)     :: delka2,delkr2, exparg
   complex (Long) :: val
!   complex (Long),dimension(Nxa,Nxr) :: vals

   real :: elapsed

   delka2=2e0_Long*delka
   delkr2=2e0_Long*delkr

   !store in denmat2. Must zero it out first.
   denmat2=cmplx(0e0_Long,0e0_Long,Long)

   call cpu_time(elapsed)
   write(*,*)'time at starting even k ft:',elapsed

   !first transform even ka and kr
   !first calculate rho'(ja,ka')
   do ikr2=-Nkr2/2,Nkr2/2-1
    ikr=ikr2*2

    do ika2=-Nka2,Nka2-1
     ika=ika2*2

     val=cmplx(0e0_Long,0e0_Long,Long)
     do ixa=-Nxa2,Nxa2-1
      do ixr=-Nxr2,Nxr2-1

       exparg=xa(ixa)*delkr2*ikr2 &
             +xr(ixr)*delka2*ika2
       val=val+denmat(ixa,ixr)*exp(-imagi*exparg)

      enddo !ixr
     enddo !ixa

     val=val*delxa*delxr*invpi
     denmat2(ikr,ika)=val

    enddo !ikr2
   enddo !ika2

   call cpu_time(elapsed)
   write(*,*)'time at starting odd k ft:',elapsed

   !now transform odd ka and kr
   do ikr2=-Nkr2/2,Nkr2/2-1
    ikr=ikr2*2+1

    do ika2=-Nka2,Nka2-1
     ika=ika2*2+1

     val=cmplx(0e0_Long,0e0_Long,Long)
     do ixa=-Nxa2,Nxa2-1
      do ixr=-Nxr2,Nxr2-1

       exparg=xa(ixa)*delkr2*0.5_Long &
             +xr(ixr)*delka2*0.5_Long &
             +xa(ixa)*delkr2*ikr2 &
             +xr(ixr)*delka2*ika2

       val=val+getDenX(ixa,ixr)*exp(-imagi*exparg)

      enddo !ixr
     enddo !ixa

     val=val*delxa*delxr*invpi
     denmat2(ikr,ika)=val

    enddo !ikr2
   enddo !ika2

   call cpu_time(elapsed)
   write(*,*)'time at finishing odd k ft:',elapsed

   do ikr=-Nkr2,Nkr2-1
    do ika=-Nka,Nka2-1
     call setDenK(ikr,ika,denmat2(ikr,ika))
    enddo
   enddo

   denState=MOMENTUM

  end subroutine transform_x_to_k_norepeat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

  subroutine transform_x_to_w_norepeat
   !! transform_x_to_w_norepeat - transforms without repeating, as notes BWB 2011-03-11p2
   use phys_cons
   implicit none

   integer :: ixa,ika2,ika,ixr,sgnfac
   real (Long)  :: delka2

   complex (Long), dimension(-Nxr2:Nxr2-1) :: array

   delka2=delka*2e0_Long

!   denmat2=cmplx(1d0,0d0,8)

   !do even ika first
   do ixa=0,Nxa2-1
    !construct array to transform
    array=czero
    do ixr=-Nxr2,Nxr2-1
     array(ixr)=denmat(ixa,ixr)+denmat(ixa-Nxa2,ixr)
    enddo !ika2

    !transform!
    do ika2=-Nka2,Nka2-1
     ika=ika2*2
     denmat2(ixa,ika)=czero
!     val=cmplx(0d0,0d0,8)
     do ixr=-Nxr2,Nxr2-1
      denmat2(ixa,ika)=denmat2(ixa,ika)+array(ixr) &
                                        *exp(-imagi*xr(ixr)*delka2*ika2)
     enddo !ixr
    enddo !ika2
   enddo !ixa

   !do odd ika now
   do ixa=0,Nxa2-1
    !construct array to transform
    array=czero
    do ixr=-Nxr2,Nxr2-1
     array(ixr)=(denmat(ixa,ixr)-denmat(ixa-Nxa2,ixr)) &
                *exp(-imagi*delxr*delka2*0.5_Long*ixr)
    enddo !ika2

    !transform!
    do ika2=-Nka2,Nka2-1
     ika=ika2*2+1
     denmat2(ixa,ika)=czero
     do ixr=-Nxr2,Nxr2-1
      denmat2(ixa,ika)=denmat2(ixa,ika)+array(ixr) &
                                        *exp(-imagi*xr(ixr)*delka2*ika2)
     enddo !ixr
    enddo !ika2
   enddo !ixa

   !copy to ixa<0 half of matrix, with (-1)**ika factor
   do ixa=-Nxa2,-1
    sgnfac=1
    do ika=-Nka,Nka-1
     denmat2(ixa,ika)=sgnfac*denmat2(ixa+Nxa2,ika)
     sgnfac=-sgnfac
    enddo !ika
   enddo !ixa

   denmat=denmat2*delxr*invsqrt2pi

!   do ixa=-Nxa2,Nxa2-1
!    do ika=-Nka,Nka-1
!     if(abs(REAL(denmat(ixa,ika)))<1d-40) then
!      write(*,*)'bad:',ixa,ika,denmat(ixa,ika)
!     endif
!    enddo
!   enddo

   denState=WIGNER

  end subroutine transform_x_to_w_norepeat
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine transform_x_to_w_norepeat_fft
   !! transform_x_to_w_norepeat - transforms without repeating, as notes BWB 2011-03-11p2
   use lib_fftw
   use phys_cons
   implicit none

   integer :: ixa,ika2,ika,ixr,sgnfac
   real (Long) :: delka2

!   integer*8 :: plan

   complex (Long), dimension(-Nxr2:Nxr2-1) :: array

   delka2=delka*2e0_Long

!   denmat2=cmplx(1d0,0d0,8)

!   call dfftw_plan_dft_1d(plan,Nxr,array,array, FFTW_FORWARD, FFTW_MEASURE)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ARRAY,IKA,IKA2,IXR,SGNFAC)

   !do even ika first
!$OMP DO
   do ixa=0,Nxa2-1
    !construct array to transform
    array=czero
    sgnfac=1
    do ixr=-Nxr2,Nxr2-1
     array(ixr)=(denmat(ixa,ixr)+denmat(ixa-Nxa2,ixr))*sgnfac
     sgnfac=-sgnfac
    enddo !ika2

    !transform!
!   call dfftw_execute_dft(plan,array,array)


   call ft_z2z_1d(array,array,Nxr)

    sgnfac=1
    do ika2=-Nka2,Nka2-1
     ika=ika2*2
     array(ika2)=array(ika2)*sgnfac
     denmat2(ixa,ika)=array(ika2)
     sgnfac=-sgnfac
    enddo

   enddo !ixa
!$OMP END DO

!$OMP END PARALLEL


!   call dfftw_destroy_plan(plan)

   !do odd ika now
   do ixa=0,Nxa2-1
    !construct array to transform
    array=czero
    sgnfac=1
    do ixr=-Nxr2,Nxr2-1
     array(ixr)=(denmat(ixa,ixr)-denmat(ixa-Nxa2,ixr)) &
                *exp(-imagi*xr(ixr)*delka2*0.5_Long) &
                *sgnfac
     sgnfac=-sgnfac
    enddo !ika2

    !transform!
    call ft_z2z_1d(array,array,Nxr)

    sgnfac=1
    do ika2=-Nka2,Nka2-1
     ika=ika2*2+1
     array(ika2)=array(ika2)*sgnfac
     denmat2(ixa,ika)=array(ika2)
     sgnfac=-sgnfac
    enddo
   enddo !ixa

   !copy to ixa<0 half of matrix, with (-1)**ika factor
   do ixa=-Nxa2,-1
    sgnfac=1
    do ika=-Nka,Nka-1
     denmat2(ixa,ika)=sgnfac*denmat2(ixa+Nxa2,ika)
     sgnfac=-sgnfac
    enddo !ika
   enddo !ixa

   denmat=denmat2*delxr*invsqrt2pi

!   do ixa=-Nxa2,Nxa2-1
!    do ika=-Nka,Nka-1
!     if(abs(REAL(denmat(ixa,ika)))<1d-40) then
!      write(*,*)'bad:',ixa,ika,denmat(ixa,ika)
!     endif
!    enddo
!   enddo

   denState=WIGNER

  end subroutine transform_x_to_w_norepeat_fft

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine transform_w_to_k_norepeat
   !! transform_x_to_w_norepeat - transforms without repeating, as notes BWB 2011-03-11p2
   use phys_cons
   implicit none

   integer :: ixa,ika2,ika,ikr,ikr2
   real (Long)  :: delkr2

   complex (Long), dimension(-Nxa2:Nxa2-1) :: array

   delkr2=delkr*2e0_Long

   denmat2=czero

   !do even ikr first
   do ika2=-Nka2,Nka2-1
    ika=ika2*2
    !construct array to transform
    array=czero
    do ixa=-Nxa2,Nxa2-1
     array(ixa)=denmat(ixa,ika)
    enddo !ixa

    !transform!
    do ikr2=-Nkr2/2,Nkr2/2-1
     ikr=ikr2*2
     denmat2(ikr,ika)=czero
!     val=cmplx(0d0,0d0,8)
     do ixa=-Nxa2,Nxa2-1
      denmat2(ikr,ika)=denmat2(ikr,ika)+array(ixa) &
                                        *exp(-imagi*xa(ixa)*delkr2*ikr2)
     enddo !ixa
    enddo !ikr2
   enddo !ika2

   !do odd ikr,ika now
   do ika2=-Nka2,Nka2-1
    ika=ika2*2+1
    !construct array to transform
    array=czero
    do ixa=-Nxa2,Nxa2-1
     array(ixa)=denmat(ixa,ika)*exp(-imagi*xa(ixa)*delkr2*0.5_Long)
    enddo !ixa

    !transform!
    do ikr2=-Nkr2/2,Nkr2/2-1
     ikr=ikr2*2+1
     denmat2(ikr,ika)=czero
!     val=cmplx(0d0,0d0,8)
     do ixa=-Nxa2,Nxa2-1
      denmat2(ikr,ika)=denmat2(ikr,ika)+array(ixa) &
                                        *exp(-imagi*xa(ixa)*delkr2*ikr2)
     enddo !ixa
    enddo !ikr2
   enddo !ika2

   denmat=denmat2*delxa*invsqrt2pi

   denState=MOMENTUM

  end subroutine transform_w_to_k_norepeat
 
END MODULE mesh

