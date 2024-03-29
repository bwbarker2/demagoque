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

 integer :: Nxan    ! minimum logical index of cell in xa
 integer :: Nxax    ! maximum logical index of cell in xa
 integer :: Nxrn    ! minimum logical index of cell in xr
 integer :: Nxrx    ! maximum logical index of cell in xr
 integer :: Nkan    ! minimum logical index of cell in ka
 integer :: Nkax    ! maximum logical index of cell in ka
 integer :: Nkrn    ! minimum logical index of cell in kr
 integer :: Nkrx    ! maximum logical index of cell in kr

 REAL (Long)    :: delxa   ! interval in x_average
 REAL (Long)    :: delxr   ! interval in x_relative
 REAL (Long)    :: delka   ! in momentum
 REAL (Long)    :: delkr

 integer        :: ixr0  ! index of xr where xr=0
 integer        :: ikr0  ! index of kr where kr=0
 integer        :: ika0
 integer        :: ixa0

 real (Long) :: norm_thy  ! theoretical norm (what it 'should' be)

 ! factor to change units to 3D density matrix (calculated in initial.f90)
 REAL (Long) :: facd
 REAL (Long), DIMENSION(:), ALLOCATABLE :: xa,ka,xr,kr   ! coord of grid point  

 real (Long), dimension(:,:), allocatable :: xx1,xx2,kk1,kk2 !< (x,x') coordinates as function of (ixa,ixr) indices

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
!
! !> calculates the density at a grid point, gaussian-smeared over the
! !! surrounding cells with a width sigma.
! complex function mesh_getDensitySmooth(ix,sigma)
! 
! integer,    intent(in) :: ix      !< index of density in center of gaussian
! real(Long), intent(in) :: sigma   !< width of gaussian smear
!
! integer :: ii !loop variable
! real(Long), allocatable, dimension(:), save :: erfs
!
! real(Long), allocatable, save, dimension(:) :: sigmas !< array of sigmas that have gaussian distributions calculated
! real(Long), allocatable, save, dimension(:,:) :: gaussians !< array of gaussians, first dim is key, second is value
! integer :: nsigs !< number of sigmas stored
! integer,save :: isig  !< key of current sigma
!
! !if first run, then allocate arrays
! if(.not.allocated(sigmas)) then
!  allocate(sigmas(1),gaussians(1,Nxan:Nxax))
!  allocate(erfs(Nxan-1:Nxax+1))
!  sigmas(1)=sigma
!  call mesh_calcSmoothGauss(gaussians(1,:),sigmas(1))
! endif
!
! !calculate new gaussian list
! do ii=Nxan-1,Nxax+1
!  erfs(ii)=erf(ii/(sqrt(2e0_Long)*sigma))
! enddo
!
! end function mesh_getDensitySmooth
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! !> Calculates the integral of a gaussian, for each cell on the diagonal.
! subroutine mesh_calcSmoothGauss(garray,sig)
!  real(Long), intent(out), dimension(*) :: garray
!  real(Long), intent(in)  :: sig
!
!  real(Long), dimension(:), allocatable :: erfs
!
!  integer :: ii
!
!  allocate(erfs(Nxan-1:Nxax+1))
!
!  do ii=Nxan-1,Nxax+1
!   erfs(ii)=ii/(sqrt(2e0_Long)*sig)
!  enddo
!
!  erfs=erf(erfs)
!
!  do ii=Nxan,Nxax
!   garray(ii)=0.5_Long*(erfs(ii+1)-erfs(ii-1))
!  enddo
!
! end subroutine mesh_calcSmoothGauss
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 integer function getNearestIndexX(xx) result(ixx)
  use prec_def
  implicit none

  real (Long), intent(in) :: xx

  real (Long) :: aixx  ! interpolated index for result
  real (Long), dimension(Nxan:Nxax) :: ixen  !array of indices
 
  integer :: ixa,ki

!  ixx=nint((xx-xan)/delxa)
 
  do ixa=Nxan,Nxax
   ixen(ixa)=ixa
   write(*,*)ixa,ixen(ixa),xa(ixa)
  enddo

  ki=1
  call lin_int(xa(Nxan:Nxax),ixen,Nxa,xx,aixx,ki)

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
   real (Long) :: xan,xrn,krn,kan !minimums in coordinates

 
   Nxa2=Nxa/2
   Nxr2=Nxr/2
   Nkr2=Nxa2
   Nka2=Nxr2
   Nkr=Nxa
   Nka=Nxr


   !facd calc'd here because can't initialize with non-integer exponents
   facd=sqrt(5e0_Long/3e0_Long)*(deg*pi*(rho0**2)/6e0_Long)**(1e0_Long/3e0_Long) 
   !set minimum and maximum indices for the density matrix. These limits will
   !be used every time something loops over the entire logical matrix. These are set in
   !in a way that makes the src/initial.f90 (copyExtra) routine not needed, as
   !the extra sections are looped over as well.

!   Nxan=-Nxa2

   Nxan=0
   Nxrn=0
   Nkrn=0
   Nkan=0

   if(useFrameXXP) then
    Nxax=Nxa-1
    Nxrx=Nxr-1
    Nkrx=Nkr-1
    Nkax=Nka-1
    delxa=2e0_Long*xLa/Nxa
    delxr=2e0_Long*xLr/Nxr
    delkr=pi/xLa
    delka=pi/xLr
    xan=-xLa
    xrn=-xLr
    krn=-delkr*((Nkrx-Nkrn+1)/2)  !integer division is intentional here
    kan=-delka*((Nkax-Nkan+1)/2)   !integer division is intentional here
   else ! if(useMeshXAR2) then
    Nxax=2*Nxa-1
    Nxrx=2*Nxr-1
    Nkrx=2*Nkr-1
    Nkax=2*Nka-1
    delxa=xLa/Nxa
    delxr=2e0*xLr/Nxr
    delkr=pi/xLa
    delka=pi/(2e0_Long*xLr)
    xan=-xLa
    xrn=-2*xLr
!    xrn=-xLr
    krn=-delkr*((Nkrx-Nkrn+1)/2)
    kan=-delka*((Nkax-Nkan+1)/2)
    ixa0=(Nxax-Nxan+1)/2
    ixr0=(Nxrx-Nxrn+1)/2
    ikr0=(Nkrx-Nkrn+1)/2
    ika0=(Nkax-Nkan+1)/2
  endif

   ! allocate arrays
   ! xa has the following bounds so that it can easily be used in evol_x's
   !  linear interpolation
   allocate(xa(Nxan-1:Nxax+1), kr(Nkrn:Nkrx), xr(Nxrn:Nxrx), ka(Nkan:Nkax))
   allocate(denmat(Nxan:Nxax,Nxrn:Nxrx)) 
   allocate(denmat2(Nxan:Nxax,Nxrn:Nxrx))
   allocate(denDiagX(Nxan:Nxax))
   allocate(denDiagK(Nkan:Nkax))
   allocate(potDiag(Nxan-1:Nxax+1))
   allocate(den_re(Nxan:Nxax,Nxrn:Nxrx))
   allocate(xx1(Nxan:Nxax,Nxrn:Nxrx), xx2(Nxan:Nxax,Nxrn:Nxrx))
   allocate(kk1(Nkrn:Nkrx,Nkan:Nkax), kk2(Nkrn:Nkrx,Nkan:Nkax))

   potDiag=0e0_Long

   if(useMeshShifted) then
    shift=0.5_Long
   else
    shift=0.e0_Long
   endif

   do ixa=Nxan-1,Nxax+1
    xa(ixa)=xan+delxa*(ixa-Nxan+shift)
!    write(*,*)'ixa,xa=',ixa,xa(ixa)
   enddo

!   write(*,*)'ixr0=',ixr0
   do ixr=Nxrn,Nxrx
    xr(ixr)=xrn+delxr*(ixr-Nxrn+shift)
!    write(*,*)'ixr,xr=',ixr,xr(ixr)
   enddo

!   write(*,*)'ikr0=',ikr0
   do ikr=Nkrn,Nkrx
    kr(ikr)=krn+delkr*(ikr-Nkrn+shift)
!    write(*,*)'ikr,kr=',ikr,kr(ikr)
   enddo

   do ika=Nkan,Nkax
    ka(ika)=kan+delka*(ika-Nkan+shift)
!    write(*,*)'ika,ka=',ika,ka(ika)
   enddo

   kLa=max(ka(Nkan),ka(Nkax))
!   write(*,*)'kLa=',kLa

   !set conversion from indices to coordinates
   do ixa=Nxan,Nxax
    do ixr=Nxrn,Nxrx
     if(useFrameXXP) then
      xx1(ixa,ixr)=xa(ixa)
      xx2(ixa,ixr)=xr(ixr)
    else
      xx1(ixa,ixr)=xa(ixa)+0.5_Long*xr(ixr)
      xx2(ixa,ixr)=xa(ixa)-0.5_Long*xr(ixr)
      ! if xx is outside the box, move it in periodically
      ! xx1 needs .ge. and xx2 needs .gt., to keep it periodic
      if(xx1(ixa,ixr).ge. xLa)xx1(ixa,ixr)=xx1(ixa,ixr)-xLa-xLa
      if(xx1(ixa,ixr).le.-xLa)xx1(ixa,ixr)=xx1(ixa,ixr)+xLa+xLa
      if(xx2(ixa,ixr).gt. xLa)xx2(ixa,ixr)=xx2(ixa,ixr)-xLa-xLa
      if(xx2(ixa,ixr).lt.-xLa)xx2(ixa,ixr)=xx2(ixa,ixr)+xLa+xLa

    endif
    enddo
   enddo

   do ikr=Nkrn,Nkrx
    do ika=Nkan,Nkax
     if(useFrameXXP) then
      kk1(ikr,ika)=kr(ikr)
      kk2(ikr,ika)=ka(ika)
     else
      kk1(ikr,ika)=ka(ika)+0.5_Long*kr(ikr)
      kk2(ikr,ika)=ka(ika)-0.5_Long*kr(ikr)
     endif
!     write(*,*)'k1,k2,kk1,kk2=',ikr,ika,kk1(ikr,ika),kk2(ikr,ika)
    enddo
   enddo

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

   if(useFrameXXP) then
    den=real(denmat(ixa,ixa))
    return
   endif

   if(.not.useMeshShifted) then
    den=real(getDenX(ixa,ixr0))
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

   if(useFrameXXP) then
    den=real(denmat(ika,ika))
    return
   endif

   if(.not.useMeshShifted) then
    den=abs(getDenK(ikr0,ika))
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Sets edges of MOMENTUM mesh equal to zero. Then renormalizes so that any probability
  !! that is lost this way is spread out over entire density matrix. 
  subroutine mesh_cutoffEdgesK(npoints)
   use cons_laws
   use input_parameters
   use phys_cons

   integer, intent(in) :: npoints !< number of points to cutoff from edges

   integer :: ii &
              ,npka ! number of points to cut in ka direction

   real(Long) :: probLost ! amount of probability that is lost

   ! if using rotated frame, multiply number of points in ka direction
   if(.not.useFrameXXP) then
    npka=npoints*2
   else
    npka=npoints
   endif

   write(*,*)'mesh_cutoffEdgesK=',npoints,npka

   probLost=0e0_Long

   ! calculate amount of probability that is lost
   do ii=Nkan,Nkan+npka-1
    probLost = probLost+getDenDiagK(ii)
   enddo

   do ii=Nkax-npka+1,Nkax
    probLost = probLost+getDenDiagK(ii)
   enddo

   call ener_k

   probLost = probLost*delka

   if(probLost<0e0_Long) probLost=0e0_Long

   write(*,*)'mesh_cutoffEdgesK probLost=',probLost

   denmat(Nkrn:Nkrn+npoints-1,:) = czero
   denmat(Nkrx-npoints+1:Nkrx,:) = czero
   denmat(:,Nkan:Nkan+npka-1) = czero
   denmat(:,Nkax-npka+1:Nkax) = czero

   ! renormalize
   denmat = denmat*knum/(knum-probLost)

  end subroutine mesh_cutoffEdgesK

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine mesh_xar2_setZeroesK
  use phys_cons

  integer :: ika,ikr, ila,ilr

  do ika=Nxan,Nkax
   ila=nint(ka(ika)/delka)  ! how many grid points away from 0?

   do ikr=Nkrn,Nkrx
    ilr=nint(kr(ikr)/delkr) ! how many grid points away from 0?

    !if outside original (k,k') diamond, set to zero
    if(    (Nka-2-ila-ilr<0) &
       .or.(Nka-2-ila+ilr<0) &
       .or.(Nka  +ila+ilr<0) &
       .or.(Nka  +ila-ilr<0) &
      ) then
     call setDenK(ikr,ika,czero)
     cycle
    !set chessboard of zeroes, if ika+ikr is odd
    elseif(abs(mod( ika-ika0 + ikr-ikr0 ,2))==1) then
     call setDenK(ikr,ika,czero)
     cycle
    endif
   enddo
  enddo
 
 end subroutine mesh_xar2_setZeroesK

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine mesh_xar2_setZeroesX
  use phys_cons

  integer :: ixa,ixr,ija,ijr

  do ixa=Nxan,Nxax
   ija=nint(xa(ixa)/delxa)

   do ixr=Nxrn,Nxrx
    ijr=nint(xr(ixr)/delxr)

    !if outside original (x,x') diamond, set to zero
    if(    (Nxa-2-ija-ijr<0) &
       .or.(Nxa-2-ija+ijr<0) &
       .or.(Nxa  +ija+ijr<0) &
       .or.(Nxa  +ija-ijr<0) &
                            ) then
     call setDenX(ixa,ixr,czero)
     cycle
    !set chessboard of zeroes, if ixa+ixr is odd
    elseif(abs(mod( ixa-ixa0 + ixr-ixr0 ,2))==1) then
     call setDenX(ixa,ixr,czero)
     cycle
    endif
   enddo
  enddo
 
 end subroutine mesh_xar2_setZeroesX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine mesh_xar2_transform_k_to_x
  use phys_cons

  integer :: ixa,ikr,ixr,ika
  complex(Long), dimension(:,:), allocatable,save :: expxakr, expkaxr

  call mesh_xar2_setZeroesK

  if(.not.allocated(expxakr)) then
   allocate(expxakr(Nxan:Nxax,Nkrn:Nkrx),expkaxr(Nkan:Nkax,Nxrn:Nxrx))
   do ixa=Nxan,Nxax
    do ikr=Nkrn,Nkrx
     expxakr(ixa,ikr)=exp(imagi*xa(ixa)*kr(ikr))
    enddo
   enddo
   do ika=Nkan,Nkax
    do ixr=Nxrn,Nxrx
     expkaxr(ika,ixr)=exp(imagi*ka(ika)*xr(ixr))
    enddo
   enddo
  endif

  denmat=matmul(denmat,expkaxr)
  denmat=matmul(expxakr,denmat)

  denmat=2*delka*delkr*inv2pi*denmat

  call mesh_xar2_setZeroesX

  denState=SPACE
 end subroutine mesh_xar2_transform_k_to_x

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine mesh_xar2_transform_x_to_k
  use phys_cons

  integer :: ixa,ikr,ixr,ika
  complex(Long), dimension(:,:), allocatable,save :: expxrka, expkrxa

  call mesh_xar2_setZeroesX

  if(.not.allocated(expxrka)) then
   allocate(expxrka(Nxrn:Nxrx,Nkan:Nkax),expkrxa(Nkrn:Nkrx,Nxan:Nxax))
   do ixr=Nxrn,Nxrx
    do ika=Nkan,Nkax
     expxrka(ixr,ika)=exp(-imagi*xr(ixr)*ka(ika))
    enddo
   enddo
   do ikr=Nkrn,Nkrx
    do ixa=Nxan,Nxax
     expkrxa(ikr,ixa)=exp(-imagi*kr(ikr)*xa(ixa))
    enddo
   enddo
  endif

  denmat=matmul(denmat,expxrka)
  denmat=matmul(expkrxa,denmat)

  denmat=2*delxa*delxr*inv2pi*denmat

  call mesh_xar2_setZeroesK

  denState=MOMENTUM
 end subroutine mesh_xar2_transform_x_to_k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 !> Translates density matrix by <nshift> points. In (xa,xr) frame, this is 
 !! less complicated than (x,x') frame, since we are shifting only in xa
 !! direction and not the xr direction.
 subroutine mesh_xar2_displace(nshiftin)
  integer, intent(in) :: nshiftin !< number of whole places to shift

  integer :: nshift ! number of grid spaces to shift

  !> Because of chessboard pattern, must shift by even number of points.
  nshift = nshiftin*2

  do while(nshift<0)
   nshift=nshift+Nxax-Nxan+1
  enddo

  do while(nshift>=Nxa)
   nshift=nshift-Nxax-Nxan+1
  enddo

  if(nshift==0) return

  denmat2(Nxan:Nxan+nshift-1,:)=denmat(Nxax-nshift+1:Nxax,:)

  denmat2(Nxan+nshift:Nxax,:)=denmat(Nxan:Nxax-nshift,:)

  denmat=denmat2
  
 end subroutine mesh_xar2_displace

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !> Translates density matrix by <nshift> points. ix -> ix+nshift
 !! (in positive direction). Only for (x,x') frame.
 !!
 !! It divides the work into 4 blocks: the square block that is simply
 !! translated forward (pos-pos), the square block that is moved periodically
 !! from the end of the matrix to the beginning (neg,neg), and the two 
 !! rectangular blocks that move from one side of the diagonal to the other
 !! (neg,pos and pos,neg). For a 3x3 matrix, with nshift=1, and (a,e,i) being
 !! the diagonal that is shifted to (i,a,e), this is the result:
 !! \f[ 
 !!  \left( \begin{array}{cc|c}
 !!   a & b & c \\
 !!   d & e & f \\
 !!   \hline
 !!   g & h & i
 !!  \end{array} \right)
 !!  \rightarrow
 !!  \left( \begin{array}{c|cc}
 !!   i & g & h \\
 !!  \hline
 !!   c & a & b \\
 !!   f & d & e
 !!  \end{array} \right)
 !! \f]
 subroutine mesh_xxp_displace(nshiftin)
  use input_parameters
  integer, intent(in) :: nshiftin !< number of grid points to shift


  integer :: ix !< loop variable
  integer :: nshift

  if(.not.useFrameXXP) then
   call throwException('mesh_xxp_displacePos: Can only be called in' &
    //' non-rotated frame, and useFrameXXP=.true.', BEXCEPTION_FATAL)
  endif

  nshift=nshiftin

  do while(nshift<0)
   nshift=nshift+Nxa
  enddo

  do while(nshift>=Nxa)
   nshift=nshift-Nxa
  enddo

  if(nshift==0) return

  !negative-negative
  denmat2(Nxan:Nxan+nshift-1,Nxan:Nxan+nshift-1) &
   =denmat(Nxax-nshift+1:Nxax,Nxax-nshift+1:Nxax)

  !positive-positive
  denmat2(Nxan+nshift:Nxax,Nxan+nshift:Nxax) &
   =denmat(Nxan:Nxax-nshift,Nxan:Nxax-nshift)

  !loop over columns (neg-pos) and rows (pos-neg)
  do ix=Nxan,Nxan+nshift-1
   !neg-pos
   denmat2(ix,Nxan+nshift:Nxax)=denmat(ix+Nxa-nshift,Nxan:Nxax-nshift)

   !pos-neg
   denmat2(Nxan+nshift:Nxax,ix)=denmat(Nxan:Nxax-nshift,ix+Nxa-nshift)
  enddo

  denmat=denmat2

 end subroutine mesh_xxp_displace

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 recursive subroutine mesh_processDen
  use phys_cons
  implicit none

  if(denState.ne.WIGNER) then
   call throwException( &
    'getDenDiagX: density matrix has been changed since last processed.' &
    //'Processing now. Suggest re-ordering to avoid extra FFTs.' &
    ,BEXCEPTION_WARNING)
   call setState(WIGNER)
  else
   ! sum the momentum components to get diagonal in space
   denDiagX=delka*invsqrt2pi*real(sum(denmat(:,-Nka:Nka-1),2))

   ! sum the space components to get diagonal in momentum
   denDiagK=delxa*invsqrt2pi*real(sum(denmat(-Nxa2:Nxa2-1,:),1))
  endif  

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
  subroutine setDenX(ixa,ixr,valu)

   integer, intent(in)    :: ixa,ixr
   complex (Long), intent(in) :: valu

   denmat(ixa,ixr)=valu

   if(isDenProcessed) isDenProcessed=.false.

  end subroutine setDenX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> returns value of Wigner density matrix at index (ixa,ika)
  complex (Long) function getDenW(ixa,ika)
   implicit none

   integer, intent(in) :: ixa,ika

   getDenW=denmat(ixa,ika)

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
   complex (Long), dimension(-Nxan:Nxax,-Nxrn+Nxr2:Nxrx-Nxr2), intent(out) :: evecs

   integer :: ixa,ixr

   if(Nxa==Nxr) then
    continue
   else
    call throwException( &
         'getEigens only works for Nxa=Nxr, exiting subroutine' &
         ,BEXCEPTION_FATAL)
    return
   endif

   call setState(SPACE)

   do ixa=-Nxan,Nxax
    do ixr=-Nxrn+Nxr2,Nxrx-Nxr2
     evecs(ixa,ixr)=getDenX(ixa,ixr)
    enddo
   enddo

   call getEigenSq(evecs,Nxa,evals,evecs)

  end subroutine getDenEigens

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Sets Fourier state of system to desired state (position, wigner, or
  !! momentum. If already in desired state, does nothing.
  subroutine setState(state)
   use input_parameters
   implicit none
   
   integer, intent(in) :: state !< state to set system to

   real :: totalelapsed

   if(denState.eq.state) return

   ! if in MOMENTUM and we are cutting the K edges, do so
   if((denState==MOMENTUM).and.(useCutoffK)) then
    call mesh_cutoffEdgesK(cutoffK_ncells)
   endif

   call cpu_time(totalelapsed)

!call mesh_setReflectedLR(.true.)

    if(useFrameXXP) then
     select case (denState)
      case (SPACE)
       select case (state)
        case (MOMENTUM)
         call xxp_transform_x_to_k_dumb
        case default
         call throwException('setState: selected state not available.' &
                             ,BEXCEPTION_FATAL)
       end select !state
      case (MOMENTUM)
       select case (state)
        case (SPACE)
         call xxp_transform_k_to_x_dumb
        case default
         call throwException('setState: selected state not available.' &
                              ,BEXCEPTION_FATAL)
       end select !state
      case default
       call throwException('setState: Current state does not exist. Weird.' &
                           , BEXCEPTION_FATAL)
     end select !denState
    elseif(useMeshShifted) then

     select case (denState)

      case (SPACE)

       select case (state)

        case (WIGNER)

         call transform_x_to_w_shift
         write(*,*)'transform_x_to_w_shift'
       end select !state

      case (WIGNER)

       select case (state)

        case (SPACE)
         call transform_w_to_x_shift
         write(*,*)'transform_w_to_x_shift'
       end select !state
     end select !denState


    elseif(useMeshXAR2) then
     select case (denState)
      case (SPACE)
       select case (state)
        case (MOMENTUM)
         call mesh_xar2_transform_x_to_k
        case default
         call throwException('setState: only SPACE and MOMENTUM states' &
          // ' implemented for useMeshXAR2=.true.',BEXCEPTION_FATAL)
       end select !state
      case (MOMENTUM)
       select case (state)
        case (SPACE)
         call mesh_xar2_transform_k_to_x
        case default
         call throwException('setState: only SPACE and MOMENTUM states' &
         // ' implemented for useMeshXAR2=.true.',BEXCEPTION_FATAL)
       end select !state
     end select !denState
    end if !useMesh selection

!    select case (denState)  
!     case (SPACE)
!      select case (state)
!       case (WIGNER)
!call mesh_setReflectedLR(.true.)
!        call transform_x_to_w_norepeat_fft
!call mesh_setReflectedLR(.false.)
!        call transform_x_to_wigner_dumb
!        call transform_x_to_w_dumb_kshift
!        call cpu_time(elapsed)
!        write(*,*)'transform_x_to_w:',elapsed-totalelapsed,'seconds'
!       case (MOMENTUM)
!        call transform_x_to_k_norepeat
!call mesh_setReflectedLR(.true.)
!        call transform_x_to_w_norepeat_fft
!call mesh_setReflectedLR(.false.)
!        call transform_x_to_wigner_dumb
!call mesh_setReflectedLR(.true.)
!        call transform_wigner_to_k_fft_exp
!call mesh_setReflectedLR(.false.)
!        call transform_wigner_to_k_dumb
!        call transform_w_to_k_norepeat
!        write(*,*)'transform_x_to_k:',etime(elapsed)-totalelapsed,'seconds'
!      end select
!  
!     case (WIGNER)
!      select case (state)
!       case (SPACE)
!        call transform_w_to_x_norepeat_fft
!        call transform_wigner_to_x_dumb
!        call cpu_time(elapsed)
!        write(*,*)'transform_w_to_x:',elapsed-totalelapsed,'seconds'
!       case (MOMENTUM)
!        call transform_wigner_to_k_dumb
!        call transform_w_to_k_norepeat
!call mesh_setReflectedLR(.true.)
!        call transform_wigner_to_k_fft_exp
!call mesh_setReflectedLR(.false.)
!        call cpu_time(elapsed)
!        write(*,*)'transform_w_to_k:',elapsed-totalelapsed,'seconds'
!      end select
!  
!     case (MOMENTUM)
!      select case (state)
!       case (WIGNER)
!        call transform_k_to_w_fft_norepeat
!        call transform_k_to_wigner_fft_exp
!        call transform_k_to_wigner_dumb
!        call cpu_time(elapsed)
!        write(*,*)'transform_k_to_w:',elapsed-totalelapsed,'seconds'
!       case (SPACE)
!        call transform_k_to_w_fft_norepeat
!        call transform_k_to_wigner_fft_exp
!        call transform_k_to_wigner_dumb
!        write(*,*)'transform_k_to_w__:',etime(elapsed)-totalelapsed,'seconds'
!        call transform_w_to_x_norepeat_fft
!        call transform_wigner_to_x_dumb
!        write(*,*)'transform_k_to_x:',etime(elapsed)-totalelapsed,'seconds'
!      end select
!    end select
!call mesh_setReflectedLR(.false.)
!   endif
!   endif !useMeshShifted



   ! if in MOMENTUM and we are cutting the K edges, do so
   if((denState==MOMENTUM).and.(useCutoffK)) then
    call mesh_cutoffEdgesK(cutoffK_ncells)
   endif

  end subroutine setState

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine transform_w_to_x_shift
   use phys_cons
   implicit none

!   real(Long), dimension(0:Nxr-1,0:Nxr-1) :: sins,coss
   real (Long) :: denr,deni

   integer :: ixa, ixr, ika
!   real(Long),dimension(:), allocatable :: iip5  !integer+0.5 array
!   real(Long) :: piN 

   do ixa=-Nxa2,Nxa2-1

    do ixr=-Nxr,Nxr-1

     denr=0e0_Long
     deni=0e0_Long

     do ika=0,Nka-1
      denr = denr &
             + real(denmat(ixa,ika)+denmat(ixa,-ika-1)) &
               * cos(pi*(ixr+0.5_Long)*(ika+0.5_Long)/Nxr)
      deni = deni &
             + real(denmat(ixa,ika)-denmat(ixa,-ika-1)) &
               * sin(pi*(ixr+0.5_Long)*(ika+0.5_Long)/Nxr)
     enddo !ika
     denmat2(ixa,ixr)=cmplx(denr,deni,Long)
    enddo !ixr
   enddo !ixa

   denmat=delka*invsqrt2pi*denmat2

   denState=SPACE

  end subroutine transform_w_to_x_shift 
  
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
!   complex (Long) :: array(Nxrn:Nxrx)
!   real(Long), dimension(:,:), allocatable,save :: coses, sines
   complex(Long), dimension(:,:), allocatable,save :: exps

   if(.not.allocated(exps)) then
!    allocate(coses(1:Nxrx,Nkan:Nkax))
!    allocate(sines(1:Nxrx,Nkan:Nkax))
    allocate(exps(Nxrn:Nxrx,Nkan:Nkax))
    do ika=Nkan,Nkax
     do ixr=Nxrn,Nxrx
!      coses(ixr,ika)=cos(xr(ixr)*ka(ika))
!      sines(ixr,ika)=sin(xr(ixr)*ka(ika))
      exps(ixr,ika)=exp(-imagi*xr(ixr)*ka(ika))
     enddo
    enddo
   endif

!write(*,*)'x to w dumb go'
  
   do ixa=Nxan,Nxax
  
!    array=cmplx(0e0_Long,0e0_Long,Long)
  
    ! fill arrays to be transformed
!    do ixr=Nxrn,Nxrx
!     array(ixr)=getDenX(ixa,ixr)
!    enddo
  
    do ika=Nkan,Nkax
!     denmat2(ixa,ika)=cmplx(0d0,0d0,8)
     denmat2(ixa,ika)=czero
     do ixr=Nxrn,Nxrx

!      !debug up-down symmetry
!      if(ixr.ne.-Nxr)then
!       if(array(ixr).ne.array(-ixr)) then
!        write(*,*)'ixa,ixr,diff',ixa,ixr,array(ixr)-array(-ixr)
!       endif
!      endif

!      exparg=mod(exparg,2*pi)  !this doesn't seem to have any effect
!      denmat2(ixa,ika)=denmat2(ixa,ika)+array(ixr)*exp(-imagi*exparg)
      denmat2(ixa,ika)=denmat2(ixa,ika) &
!                       +2.d0*(REAL(denmat(ixa,ixr))*coses(ixr,ika) &
!                              +AIMAG(denmat(ixa,ixr))*sines(ixr,ika))
          +denmat(ixa,ixr)*exps(ixr,ika)

     enddo
!     denmat2(ixa,ika)=denmat2(ixa,ika) &
!                      +denmat(ixa,Nxr)*coses(Nxrx,ika) &
!                      +denmat(ixa,ixr0)

     denmat2(ixa,ika)=delxr*denmat2(ixa,ika)*invsqrt2pi
  
     ! if the cell is unreasonably large, write out
!     if(REAL(denmat2(ixa,ika))>2.d0) write(*,*)ixa,ika,denmat2(ixa,ika)
    enddo
  
   enddo
  
   do ixa=Nxan,Nxax
    do ika=Nkan,Nkax
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

  subroutine transform_x_to_w_shift
   use phys_cons
   implicit none

!   real(Long), dimension(0:Nxr-1,0:Nxr-1) :: sins,coss

   integer :: ixa, ixr, ika
!   real(Long),dimension(:), allocatable :: iip5  !integer+0.5 array
!   real(Long) :: piN 

   do ixa=-Nxa2,Nxa2-1

    do ika=-Nka,Nka-1

     denmat2(ixa,ika)=czero

     do ixr=0,Nxr-1
      denmat2(ixa,ika) = denmat2(ixa,ika) &
         +  real(denmat(ixa,ixr))*cos(pi*(ixr+0.5_Long)*(ika+0.5_Long)/Nxr) &
         + aimag(denmat(ixa,ixr))*sin(pi*(ixr+0.5_Long)*(ika+0.5_Long)/Nxr)
     enddo !ixr
    enddo !ika
   enddo !ixa

   denmat=2e0_Long*delxr*invsqrt2pi*denmat2

   denState=WIGNER

   call mesh_processDen

!!! The following is for optimizing the algorithm to minimize arithmetic operations and to maximize use of whole array operations

!   piN = pi/Nxr

!   allocate(iip5(-max(Nxa,Nxr)/2:max(Nxa,Nxr)/2-1))

!   do ixr=lbound(iip5),ubound(iip5)
!    iip5=ixr+0.5_Long
!   enddo

!   do ika=-Nka,-1
!    do ixr=0,Nxr-1
!     coss(ixr,ika)=cos(piN*iip5(ixr)*iip5(ika))
!     sins(ixr)=sin(piN*iip5(ixr)*iip5(ika))
!    enddo
!   enddo
!
!   do ika=0,Nka-1
!    coss(ika)=coss(-ika-1)
!    sins(ika)=-sins(-ika-1)
!   enddo

!   sum(real(denmat(ixa,0:Nxr-1)
    
  end subroutine transform_x_to_w_shift

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

!   real (Long) :: exparg
   integer :: ixa,ixr,ika
   complex(Long), dimension(:,:), allocatable,save :: exps

   if(.not.allocated(exps)) then
    allocate(exps(Nxrn:Nxrx,Nkan:Nkax))
    do ika=Nkan,Nkax
     do ixr=Nxrn,Nxrx
      exps(ixr,ika)=exp(imagi*xr(ixr)*ka(ika))
     enddo
    enddo
   endif


   denmat2=CMPLX(0e0_Long,0e0_Long,Long)

   do ixa=Nxan,Nxax
    do ixr=Nxrn,Nxrx
!     do ika=-Nka,Nka-1
!      denmat2(ixa,ixr)=denmat2(ixa,ixr)+denmat(ixa,ika)*exp(imagi*delxr*delka*ixr*ika)
!     enddo
     do ika=Nkan,Nkax

      denmat2(ixa,ixr)=denmat2(ixa,ixr) &
                       +real(denmat(ixa,ika))*exps(ixr,ika)
!                       +denmat(ixa,-ika)*(cos(exparg)-imagi*sin(exparg))
     enddo
!     denmat2(ixa,ixr)=denmat2(ixa,ixr) &
!                      +denmat(ixa,-Nka)*cos(xr(ixr)*ka(-Nka)) &
!                      +denmat(ixa,0)

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
!   real(Long) :: trigarg
   complex(Long), dimension(:,:), allocatable,save :: exps

   if(.not.allocated(exps)) then
    allocate(exps(Nxan:Nxax,Nkrn:Nkrx))
    do ikr=Nkrn,Nkrx
     do ixa=Nxan,Nxax
      exps(ixa,ikr)=exp(-imagi*xa(ixa)*kr(ikr))
     enddo
    enddo
   endif


   do ika=Nkan,Nkax
    do ikr=Nkrn,Nkrx
     denmat2(ikr,ika)=czero
!     do ixa=-Nxa2,Nxa2-1
     do ixa=Nxan,Nxax
!      denmat2(ikr,ika)=denmat2(ikr,ika)+getDen(ixa,ika)*exp(-imagi*delxa*delkr*ixa*ikr)
      denmat2(ikr,ika)=denmat2(ikr,ika) &
                       +denmat(ixa,ika)*exps(ixa,ikr)
!                       +denmat(ixa,ika)*(cos(trigarg) &
!                                         -imagi*sin(trigarg))
!                       +denmat(-ixa,ika)*(cos(trigarg) &
!                                         +imagi*sin(trigarg))
     enddo !ixa
!     denmat2(ikr,ika)=denmat2(ikr,ika) &
!                      +denmat(-Nxa2,ika)*cos(xa(-Nxa2)*kr(ikr)) &
!                                          -imagi*sin(xa(-Nxa2)*kr(ikr))) &
!                      +denmat(0,ika)

     denmat2(ikr,ika)=denmat2(ikr,ika)*delxa*invsqrt2pi

     !Analytically (see notes BWB 2011-01-25p2), for ika+ikr odd, denmat should
     ! be exactly zero. Make it so. This subroutine is not built for speed,
     ! so leave the above computation alone, but add the following:
!     if(abs(mod(ika+ikr,2))==1)denmat2(ikr,ika)=0e0_Long

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
!   real(Long), dimension(:,:), allocatable,save :: coses, sines
   complex(Long), dimension(:,:), allocatable,save :: exps

   if(.not.allocated(exps)) then
!    allocate(coses(Nxan:Nxax,0:Nkrx))
!    allocate(sines(Nxan:Nxax,0:Nkrx))
    allocate(exps(Nxan:Nxax,Nkrn:Nkrx))
    do ikr=Nkrn,Nkrx
     do ixa=Nxan,Nxax
!      coses(ixa,ikr)=cos(xa(ixa)*kr(ikr))
!      sines(ixa,ikr)=sin(xa(ixa)*kr(ikr))
       exps(ixa,ikr)=exp(imagi*xa(ixa)*kr(ikr))
     enddo
    enddo
   endif


   do ika=Nkan,Nkax
    do ixa=Nxan,Nxax
     denmat2(ixa,ika)=czero
      do ikr=Nkrn,Nkrx
!      denmat2(ikr,ika)=denmat2(ikr,ika)+getDen(ixa,ika)*exp(-imagi*delxa*delkr*ixa*ikr)
      denmat2(ixa,ika)=denmat2(ixa,ika) &
!                       +2.d0*(real(denmat(ikr,ika))*coses(ixa,ikr) &
!                              -aimag(denmat(ikr,ika))*sines(ixa,ikr))
        + denmat(ikr,ika)*exps(ixa,ikr)
!                       +REAL(denmat(ikr,ika))*cos(trigarg) &
!                       -Dimag(denmat(ikr,ika))*sin(trigarg)
!                       +denmat(-ikr,ika)*(cos(trigarg) &
!                                         -imagi*sin(trigarg))
     enddo !ixa
!     denmat2(ixa,ika)=denmat2(ixa,ika) &
!                      +denmat(Nkrx,ika)*coses(ixa,Nkrx) &
!                      +denmat(0,ika)

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
   integer(kind=8), save :: plan_cos, plan_sin  !< plans for FFTW

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

!   integer(kind=8) :: plan

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> transforms the density matrix in (x,x') coordinates from SPACE to MOMENTUM.
  !! Uses matrix multiplication:
  !! \f{eqnarray*}{
  !!    \rho(j,j') &=& \rho(k,k') e^{i k j} e^{- i k' j'} \\
  !!               &=& e^{i k j} \rho(k,k') e^{-i k' j'} \\
  !!               &=& A_{j k} B_{k k'} C_{k' j'}
  !! \f}
  !! Note that A is the conjugate transpose of C.
  subroutine xxp_transform_k_to_x_dumb
   use phys_cons

   integer :: ix,ik
   complex(Long), dimension(:,:), allocatable,save :: expxk, expxkp

   if(.not.allocated(expxk)) then
    allocate(expxk(Nxan:Nxax,Nxan:Nxax),expxkp(Nxan:Nxax,Nxan:Nxax))
    do ix=Nxan,Nxax
     do ik=Nxan,Nxax
      expxk(ix,ik)=exp(imagi*xa(ix)*kr(ik))
      expxkp(ik,ix)=exp(-imagi*xr(ix)*ka(ik))
     enddo
    enddo
   endif

   denmat=matmul(denmat,expxkp)
   denmat=matmul(expxk,denmat)

   denmat=delka*delkr*inv2pi*denmat

   denState = SPACE

  end subroutine xxp_transform_k_to_x_dumb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xxp_transform_x_to_k_dumb
   use phys_cons

   integer :: ix,ik
   complex(Long), dimension(:,:), allocatable,save :: expxk, expxkp

   if(.not.allocated(expxk)) then
    allocate(expxk(Nxan:Nxax,Nxan:Nxax),expxkp(Nxan:Nxax,Nxan:Nxax))
    do ix=Nxan,Nxax
     do ik=Nkan,Nkax
      expxk(ix,ik)=exp(-imagi*xa(ix)*kr(ik))
      expxkp(ik,ix)=exp(imagi*xr(ix)*ka(ik))
     enddo
    enddo
   endif

   denmat=matmul(denmat,expxkp)
   denmat=matmul(expxk,denmat)

   denmat=delxa*delxr*inv2pi*denmat

   denState = MOMENTUM

  end subroutine xxp_transform_x_to_k_dumb
 
END MODULE mesh

