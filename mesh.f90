MODULE mesh
 implicit none

 REAL*8  :: xLa     ! length of box in xa/2 [fm]
 REAL*8  :: xLr     ! length of box in xr/2 [fm]
 real*8  :: kLa     ! half length of average momentum box
 INTEGER :: Nxa     ! number of mesh points in xa 
 INTEGER :: Nxr     ! number of mesh points in xr
 INTEGER :: Nxa2    ! Nxa/2
 INTEGER :: Nxr2    ! Nxr/2
 INTEGER :: Nka
 integer :: Nkr
 INTEGER :: Nkr2    ! Nkr/2 = Nxa/2
 INTEGER :: Nka2    ! Nka/2 = Nxr/2
 REAL*8    :: delxa   ! interval in x_average
 REAL*8    :: delxr   ! interval in x_relative
 REAL*8    :: delka   ! in momentum
 REAL*8    :: delkr
 integer :: potInitial    ! potential  with initial state
 integer :: potFinal      ! potential for time evolution
 logical :: useAdiabatic  ! if true, use adiabatic switching
 logical :: useImCutoff   ! use imaginary off-diagonal cutoff?
  ! factor to change units to 3D density matrix (calculated in initial.f90)
 REAL*8 :: facd
  REAL*8, DIMENSION(:), ALLOCATABLE :: xa,ka,xr,kr   ! coord of grid point  
  ! density matrix, real and imaginary
  REAL*8 , DIMENSION(:,:) , ALLOCATABLE :: den_re, den_im
 complex*16, dimension(:,:), allocatable :: denmat !when I need complex, I store here
 complex*16, dimension(:,:), allocatable :: denmat2 !when I need complex, I store here
 integer denState  ! gives current coordinate space of density matrix,
                   ! according to the following integer settings:
 integer, parameter :: SPACE = 0
 integer, parameter :: WIGNER = 1
 integer, parameter :: MOMENTUM = 2
  ! iNka2 = iNxr2, iNkr2 = iNka2, set here for clarity
 INTEGER , ALLOCATABLE :: iNkr2(:),iNka2(:)
  ! meanfield potential at each (xa,xr=0) point
 real*8, allocatable :: potDiag(:)

 real*8 :: maxxim  !maximum imaginary value

 contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine initializeMesh
   use osc_pars
   use phys_cons
   implicit none

   integer :: ixa,ixr,ikr,ika
 
   Nxa2=Nxa/2
   Nxr2=Nxr/2
   Nkr2=Nxa2
   Nka2=Nxr2
   Nkr=Nxa
   Nka=Nxr

   ! allocate arrays
   allocate(xa(-Nxa2:Nxa2), kr(-Nkr2:Nkr2), xr(-Nxr:Nxr), ka(-Nka:Nka))
   allocate(denmat(-Nxa2:Nxa2,-Nxr:Nxr-1)) !2x size in xr for naive FT - BWB 2011-01-10
   allocate(denmat2(-Nxa2:Nxa2,-Nxr:Nxr-1))
   allocate(potDiag(-Nxa2:Nxa2))

   !facd calc'd here because can't initialize with non-integer exponents
   facd=dsqrt(5.d0/3.d0)*(deg*pi*(rho0**2)/6.d0)**(1.d0/3.d0) 

   write(*,*) '1D to 3D factor:',facd

   ! mesh in x and k
   delxa=(2.d0*xLa)/dble(Nxa)
   delxr=(2.d0*xLr)/dble(Nxr)
   delka=pi/(2d0*xLr)
   delkr=pi/xLa

   potDiag=0.0d0
   do ixa=-Nxa2,Nxa2
    xa(ixa)=ixa*delxa
   enddo

   do ixr=-Nxr,Nxr
    xr(ixr)=ixr*delxr
   enddo

   do ikr=-Nkr2,Nkr2
    kr(ikr)=ikr*delkr
   enddo

   do ika=-Nka,Nka
    ka(ika)=ika*delka
   enddo
 
   kLa=ka(Nka2)

  end subroutine initializeMesh

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex*16 function getDen(i1,i2)
   !! getDen - returns value of density matrix, given the current denState
   implicit none

   integer, intent(in) :: i1,i2

   select case (denState)
    case (SPACE)
     getDen=getDenX(i1,i2)
    case (WIGNER)
     getDen=getDenW(i1,i2)
    case (MOMENTUM)
     getDen=getDenK(i1,i2)
   end select

   end function getDen

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex*16 function getDenX(ixa,ixr)
   !! getDenX - returns value of spatial density matrix at index (ixa,ixr)
   implicit none

   integer, intent(in) :: ixa,ixr

   getDenX=denmat(ixa,ixr)

  end function getDenX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine setDenX(ixa,ixr,value)
   !! setDenX - sets the value of the spatial density matrix at index (ixa,ixr)
   !! NOTE: Only sets elements that are actually used in getDenX

   integer, intent(in)    :: ixa,ixr
   complex*16, intent(in) :: value

   denmat(ixa,ixr)=value

  end subroutine setDenX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex*16 function getDenW(ixa,ika)
   !! getDenW - returns value of Wigner density matrix at index (ixa,ika)
   ! NOTE: Only valid for ika<Nka2 (maybe)
   implicit none

   integer, intent(in) :: ixa,ika

   if(ixa<Nxa2) then
    getDenW=denmat(ixa,ika)
   else
    getDenW=denmat(-Nxa2,ika)
   endif

  end function getDenW

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine setDenW(ixa,ika, value)
   implicit none

   integer, intent(in) :: ixa,ika
   complex*16, intent(in) :: value

   denmat(ixa,ika)=value

  end subroutine setDenW

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex*16 function getDenK(ikr,ika)
   !! getDenK - returns value of spectral density matrix at index (ikr,ika)
   implicit none

   integer, intent(in) :: ikr,ika

   getDenK=denmat(ikr,ika)

  end function getDenK

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine setDenK(ikr,ika,value)
   !! setDenK - sets value of the spectral density matrix at index (ikr,ika)
   !! NOTE: Only valid for ikr>=0 and ika<Nka2 (maybe)
   implicit none

   complex*16, intent(in) :: value
   integer, intent(in)    :: ikr, ika

   denmat(ikr,ika)=value

  end subroutine setDenK

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine getDenEigens(evals,evecs)
   !! getEigens - finds eigenvalues, eigenvectors of spatial density matrix. This currently only works for Nxr=Nxa (square matrices)
   use lib_lapack
   implicit none

   complex*16, dimension(0:Nxa-1), intent(out) :: evals
   complex*16, dimension(-Nxa2:Nxa2-1,-Nxr2:Nxr2-1), intent(out) :: evecs

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

  subroutine setState(state)
   !! only implemented for x to wigner and inverse transforms so far
   implicit none
   
   integer, intent(in) :: state

   real :: elapsed(2)
   real :: totalelapsed  

   if(denState.NE.state) then

   totalelapsed=etime(elapsed)

    select case (denState)
  
     case (SPACE)
      select case (state)
       case (WIGNER)
        call transform_x_to_w_norepeat_fft
!        call transform_x_to_wigner_dumb
!        write(*,*)'transform_x_to_w:',etime(elapsed)-totalelapsed,'seconds'
       case (MOMENTUM)
!        call transform_x_to_k_norepeat
        call transform_x_to_w_norepeat_fft
!        call transform_x_to_wigner_dumb
        call transform_wigner_to_k_fft_exp
!        call transform_wigner_to_k_dumb
!        write(*,*)'transform_x_to_k:',etime(elapsed)-totalelapsed,'seconds'
      end select
  
     case (WIGNER)
      select case (state)
       case (SPACE)
        call transform_w_to_x_norepeat_fft
!        call transform_wigner_to_x_dumb
!        write(*,*)'transform_w_to_x:',etime(elapsed)-totalelapsed,'seconds'
       case (MOMENTUM)
!        call transform_wigner_to_k_dumb
        call transform_wigner_to_k_fft_exp
!        write(*,*)'transform_w_to_k:',etime(elapsed)-totalelapsed,'seconds'
      end select
  
     case (MOMENTUM)
      select case (state)
       case (WIGNER)
        call transform_k_to_wigner_fft_exp
!        call transform_k_to_wigner_dumb
!        write(*,*)'transform_k_to_w:',etime(elapsed)-totalelapsed,'seconds'
       case (SPACE)
        call transform_k_to_wigner_fft_exp
!        call transform_k_to_wigner_dumb
!        write(*,*)'transform_k_to_w__:',etime(elapsed)-totalelapsed,'seconds'
        call transform_w_to_x_norepeat_fft
!        call transform_wigner_to_x_dumb
!        write(*,*)'transform_k_to_x:',etime(elapsed)-totalelapsed,'seconds'
      end select
    end select
   endif
  
  end subroutine setState
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine transform_x_to_wigner_trig
   !! brute force method to test the theory, using cos/sin transforms. Theory in notes, BWB 2011-02-25p2
   use fftw_constants
   use phys_cons
   implicit none
  
   integer :: ixa,ixr,ika
   real*8 :: trigarg
   real*8 :: wigden(-Nxa2:Nxa2-1,-Nka:Nka-1)
   complex*16 :: array(0:Nxr)
  
   do ixa=-Nxa2,Nxa2-1
  
    array=0.d0
  
    ! fill arrays to be transformed
    do ixr=0,Nxr-1
     array(ixr)=getDenX(ixa,ixr)
    enddo

    array(Nxr)=getDenX(ixa,-Nxr)

    do ika=-Nka,Nka-1
     wigden(ixa,ika)=0d0
     do ixr=1,Nxr-1
      trigarg=delxr*delka*ika*ixr
      wigden(ixa,ika)=wigden(ixa,ika)+DBLE(array(ixr))*cos(trigarg) &
                                    +DIMAG(array(ixr))*sin(trigarg)
     enddo
     wigden(ixa,ika)=wigden(ixa,ika)*2d0
     wigden(ixa,ika)=wigden(ixa,ika)+array(0)+array(Nxr)*(-1)**ika
     wigden(ixa,ika)=wigden(ixa,ika)*delxr*invsqrt2pi
  
     ! if the cell is unreasonably large, write out
!     if(DBLE(wigden(ixa,ika))>2.d0) write(*,*)ixa,ika,wigden(ixa,ika)
    enddo
  
   enddo
  
   do ixa=-Nxa2,Nxa2-1
    do ika=-Nka,Nka-1
     call setDenW(ixa,ika,CMPLX(wigden(ixa,ika),0d0,8))
    enddo
   enddo
  
   denState=WIGNER
  
  end subroutine transform_x_to_wigner_trig
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
  subroutine transform_x_to_wigner_dumb
   !! brute force method to test the theory
   use fftw_constants
   use phys_cons
   implicit none
  
   integer :: ixa,ixr,ika
   real*8 :: exparg
   complex*16 :: array(-Nxr:Nxr-1)
  
   do ixa=-Nxa2,Nxa2-1
  
    array=0.d0
  
    ! fill arrays to be transformed
    do ixr=-Nxr,Nxr-1
     array(ixr)=getDenX(ixa,ixr)
    enddo
  
    do ika=-Nka,Nka-1
     denmat2(ixa,ika)=0d0
     do ixr=-Nxr,Nxr-1
      exparg=delxr*delka*ika*ixr
      denmat2(ixa,ika)=denmat2(ixa,ika)+array(ixr)*exp(-imagi*exparg)
     enddo
     denmat2(ixa,ika)=delxr*denmat2(ixa,ika)*invsqrt2pi
  
     ! if the cell is unreasonably large, write out
!     if(DBLE(denmat2(ixa,ika))>2.d0) write(*,*)ixa,ika,denmat2(ixa,ika)
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

  subroutine transform_w_to_x_norepeat_fft
   !! implements simple FFT, as per notes BWB 2011-03-28p2
   use lib_fftw
   use phys_cons
   implicit none

   integer :: ixa, ika, ixr, sgnfac
   complex*16 :: array(-Nxr:Nxr-1)

   !only compute half of it
   do ixa=-Nxa2,-1
    sgnfac=1
    do ika=-Nka,Nka-1
     array(ika)=denmat(ixa,ika)*sgnfac
     sgnfac=-sgnfac
    enddo !ika

    call ift_z2z_1d(array,array,2*Nxr)

    do ixr=-Nxr+1,Nxr-1,2
     array(ixr)=array(ixr)*(-1d0)
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

  subroutine transform_w_to_x_norepeat_fft_bad
   !! for derivation, see BWB notes 2011-03-25. Essentially this uses 2 sine transforms and 2 cosine transforms, needed because the exponent in the transform has pi/N instead of 2*pi/N.
   ! NOTE: This is unfinished, because it is much too complicated. See notes BWB 2011-03-28p2
   use lib_fftw
   use phys_cons
   implicit none

   integer :: ixa,ika,ixr
   real*8  :: array(0:Nka)

   denmat2=0.d0

   do ixa=0,Nxa2-1

    !create first input array
    array=0d0
    do ika=1,Nka-1
     array(ika)=0.5d0*(DBLE(denmat(ixa,ika))+DBLE(denmat(ixa,-ika)))
    enddo !ika
    array(0)=DBLE(denmat(ixa,0))
    array(Nka)=DBLE(denmat(ixa,-Nka))

    call ft_re_1d(array,array,Nka+1)

    do ixr=0,Nxr-1
     denmat2(ixa,ixr)=denmat2(ixa,ixr)+array(ixr)
    enddo
    denmat2(ixa,-Nxr)=denmat2(ixa,-Nxr)+array(ixr)

    !second input array
    array=0d0
    do ika=1,Nka-1
     array(ika-1)=0.5d0*(-DIMAG(denmat(ixa,ika))+DIMAG(denmat(ixa,-ika)))
    enddo !ika

    call ft_ro_1d(array,array,Nka-1)

    do ixr=1,Nxr-1
     denmat2(ixa,ixr)=denmat2(ixa,ixr)+array(ixr-1) &
                                      -DIMAG(denmat(ixa,-Nka))*(-1)**ixr &
                                      -DIMAG(denmat(ixa,0))
    enddo
    denmat2(ixa,-Nxr)=denmat2(ixa,-Nxr)+array(Nxr) &
                                       -DIMAG(denmat(ixa,-Nka)) &
                                       -DIMAG(denmat(ixa,0))
 
    !third input array
    array=0d0
    do ika=1,Nka-1
     array(ika)=0.5d0*(DIMAG(denmat(ixa,ika))+DIMAG(denmat(ixa,-ika)))
    enddo !ika
    array(0)=DIMAG(denmat(ixa,0))
    array(Nka)=DIMAG(denmat(ixa,-Nka))

    call ft_re_1d(array,array,Nka+1)

    do ixr=0,Nxr-1
     denmat2(ixa,ixr)=denmat2(ixa,ixr)+imagi*array(ixr)
    enddo
    denmat2(ixa,-Nxr)=denmat2(ixa,-Nxr)-imagi*array(ixr) !minus for conjugation

!Now do the same thing for the other 3 terms.

    array=0d0
    do ika=1,Nka-1
     array(ika)=0.5d0*(DBLE(denmat(ixa,ika))-DBLE(denmat(ixa,-ika)))
    enddo !ika
    array(0)=DIMAG(denmat(ixa,0))
    array(Nka)=DIMAG(denmat(ixa,-Nka))

    call ft_ro_1d(array,array,Nka+1)

    do ixr=0,Nxr-1
     denmat2(ixa,ixr)=denmat2(ixa,ixr)+array(ixr) &
                                      -DIMAG(denmat(ixa,-Nka))*(-1)**ixr &
                                      -DIMAG(denmat(ixa,0))
    enddo
    denmat2(ixa,-Nxr)=denmat2(ixa,-Nxr)+array(Nxr) &
                                       -DIMAG(denmat(ixa,-Nka)) &
                                       -DIMAG(denmat(ixa,0))
 

   enddo !ixa

  end subroutine transform_w_to_x_norepeat_fft_bad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine transform_wigner_to_x_trig
   ! most straightforward way to compute inverse transform
   use phys_cons
   implicit none

   integer :: ixa,ixr,ika

   denmat2=CMPLX(0d0,0d0,8)

   do ixa=-Nxa2,Nxa2-1
    do ixr=-Nxr,Nxr-1
     do ika=-Nka,Nka-1
      denmat2(ixa,ixr)=denmat2(ixa,ixr)+dble(getDenW(ixa,ika))*exp(imagi*delxr*delka*ixr*ika)
     enddo
     denmat2(ixa,ixr)=denmat2(ixa,ixr)*delka*invsqrt2pi
    enddo
   enddo
  
   denmat=denmat2
  
   denState=SPACE
  
  end subroutine transform_wigner_to_x_trig
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine transform_wigner_to_x_dumb
   ! most straightforward way to compute inverse transform
   use phys_cons
   implicit none

   integer :: ixa,ixr,ika

   denmat2=CMPLX(0d0,0d0,8)

   do ixa=-Nxa2,Nxa2-1
    do ixr=-Nxr,Nxr-1
     do ika=-Nka,Nka-1
      denmat2(ixa,ixr)=denmat2(ixa,ixr)+denmat(ixa,ika)*exp(imagi*delxr*delka*ixr*ika)
     enddo
     denmat2(ixa,ixr)=denmat2(ixa,ixr)*delka*invsqrt2pi
    enddo
   enddo
  
   denmat=denmat2
  
   denState=SPACE
  
  end subroutine transform_wigner_to_x_dumb
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine transform_k_to_wigner_trig
   !! brute force method to test the theory, using cos/sin transforms. Theory in notes, BWB 2011-02-28p1
   use fftw_constants
   use phys_cons
   implicit none
  
   integer :: ixa,ika,ikr
   real*8 :: trigarg
   real*8 :: wigden(-Nxa2:Nxa2-1,-Nka:Nka-1)
   complex*16 :: array(0:Nkr2)
  
   do ika=-Nka,Nka-1
  
    array=0.d0
  
    ! fill arrays to be transformed
    do ikr=0,Nkr2-1
     array(ikr)=getDenK(ikr,ika)
    enddo

    array(Nkr2)=getDenK(-Nkr2,ika)

    do ixa=-Nxa2,Nxa2-1
     wigden(ixa,ika)=0d0
     do ikr=1,Nkr2-1
      trigarg=delxa*delkr*ikr*ixa
      wigden(ixa,ika)=wigden(ixa,ika)+DBLE(array(ikr))*dcos(trigarg) &
                                    -DIMAG(array(ikr))*dsin(trigarg)
     enddo
     wigden(ixa,ika)=wigden(ixa,ika)*2d0
     wigden(ixa,ika)=wigden(ixa,ika)+DBLE(array(0))+DBLE(array(Nxr))*(-1)**ixa
     wigden(ixa,ika)=wigden(ixa,ika)*delkr*invsqrt2pi
  
     ! if the cell is unreasonably large, write out
!     if(DBLE(wigden(ixa,ika))>2.d0) write(*,*)ixa,ika,wigden(ixa,ika)
    enddo
  
   enddo
  
   do ixa=-Nxa2,Nxa2-1
    do ika=-Nka,Nka-1
     call setDenW(ixa,ika,CMPLX(wigden(ixa,ika),0d0,8))
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
     denmat2(ikr,ika)=0.d0
     do ixa=-Nxa2,Nxa2-1
      denmat2(ikr,ika)=denmat2(ikr,ika)+dble(getDenW(ixa,ika))*exp(-imagi*delxa*delkr*ixa*ikr)
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
   implicit none

   integer :: ixa,ika,ikr

   do ika=-Nka,Nka-1
    do ikr=-Nkr2,Nkr2-1
     denmat2(ikr,ika)=0.d0
     do ixa=-Nxa2,Nxa2-1
      denmat2(ikr,ika)=denmat2(ikr,ika)+getDen(ixa,ika)*exp(-imagi*delxa*delkr*ixa*ikr)
     enddo
     denmat2(ikr,ika)=denmat2(ikr,ika)*delxa*invsqrt2pi
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
   complex*16, dimension(-Nxa2:Nxa2-1) :: arrin,arrout

   arrin=cmplx(0d0,0d0,8)
   arrout=cmplx(0d0,0d0,8)

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
   implicit none
  
   integer :: ixa,ika,ikr
  
   do ika=-Nka,Nka-1
    do ixa=-Nxa2,Nxa2-1
     denmat2(ixa,ika)=0.d0
     do ikr=-Nkr2,Nkr2-1
      denmat2(ixa,ika)=denmat2(ixa,ika)+getDen(ikr,ika)*exp(imagi*delxa*delkr*ixa*ikr)
     enddo
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
   complex*16, dimension(-Nxa2:Nxa2-1) :: arrin,arrout

   arrin=cmplx(0d0,0d0,8)
   arrout=cmplx(0d0,0d0,8)

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

  subroutine transform_x_to_k_norepeat
   !! transform_x_to_k_nozeroes - transforms from coordinate to momentum space, with no redundancy from periodic extensions. See paper notes BWB 2011-03-11.
   use phys_cons
   implicit none

   integer    :: ikr,ikr2,ika,ika2,ixa,ixr
   real*8     :: delka2,delkr2, exparg
   complex*16 :: val
!   complex*16,dimension(Nxa,Nxr) :: vals

   real :: elapsed(2)

   delka2=2*delka
   delkr2=2*delkr

   !store in denmat2. Must zero it out first.
   denmat2=cmplx(0d0,0d0,8)

   write(*,*)'time at starting even k ft:',etime(elapsed)

   !first transform even ka and kr
   !first calculate rho'(ja,ka')
   do ikr2=-Nkr2/2,Nkr2/2-1
    ikr=ikr2*2

    do ika2=-Nka2,Nka2-1
     ika=ika2*2

     val=cmplx(0d0,0d0,8)
     do ixa=-Nxa2,Nxa2-1
      do ixr=-Nxr2,Nxr2-1

       exparg=delxa*delkr2*ixa*ikr2 &
             +delxr*delka2*ixr*ika2
       val=val+denmat(ixa,ixr)*exp(-imagi*exparg)

      enddo !ixr
     enddo !ixa

     val=val*delxa*delxr*invpi
     denmat2(ikr,ika)=val

    enddo !ikr2
   enddo !ika2


   write(*,*)'time at starting odd k ft:',etime(elapsed)

   !now transform odd ka and kr
   do ikr2=-Nkr2/2,Nkr2/2-1
    ikr=ikr2*2+1

    do ika2=-Nka2,Nka2-1
     ika=ika2*2+1

     val=cmplx(0d0,0d0,8)
     do ixa=-Nxa2,Nxa2-1
      do ixr=-Nxr2,Nxr2-1

       exparg=delxa*delkr2*5d-1*ixa &
             +delxr*delka2*5d-1*ixr &
             +delxa*delkr2*ixa*ikr2 &
             +delxr*delka2*ixr*ika2

       val=val+getDenX(ixa,ixr)*exp(-imagi*exparg)

      enddo !ixr
     enddo !ixa

     val=val*delxa*delxr*invpi
     denmat2(ikr,ika)=val

    enddo !ikr2
   enddo !ika2

   write(*,*)'time at finishing odd k ft:',etime(elapsed)

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
   real*8  :: delka2

   complex*16, dimension(-Nxr2:Nxr2-1) :: array

   delka2=delka*2d0

!   denmat2=cmplx(1d0,0d0,8)

   !do even ika first
   do ixa=0,Nxa2-1
    !construct array to transform
    array=cmplx(0d0,0d0,8)
    do ixr=-Nxr2,Nxr2-1
     array(ixr)=denmat(ixa,ixr)+denmat(ixa-Nxa2,ixr)
    enddo !ika2

    !transform!
    do ika2=-Nka2,Nka2-1
     ika=ika2*2
     denmat2(ixa,ika)=cmplx(0d0,0d0,8)
!     val=cmplx(0d0,0d0,8)
     do ixr=-Nxr2,Nxr2-1
      denmat2(ixa,ika)=denmat2(ixa,ika)+array(ixr) &
                                        *exp(-imagi*delxr*delka2*ixr*ika2)
     enddo !ixr
    enddo !ika2
   enddo !ixa

   !do odd ika now
   do ixa=0,Nxa2-1
    !construct array to transform
    array=cmplx(0d0,0d0,8)
    do ixr=-Nxr2,Nxr2-1
     array(ixr)=(denmat(ixa,ixr)-denmat(ixa-Nxa2,ixr)) &
                *exp(-imagi*delxr*delka2*0.5d0*ixr)
    enddo !ika2

    !transform!
    do ika2=-Nka2,Nka2-1
     ika=ika2*2+1
     denmat2(ixa,ika)=cmplx(0d0,0d0,8)
     do ixr=-Nxr2,Nxr2-1
      denmat2(ixa,ika)=denmat2(ixa,ika)+array(ixr) &
                                        *exp(-imagi*delxr*delka2*ixr*ika2)
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
!     if(abs(dble(denmat(ixa,ika)))<1d-40) then
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
   real*8  :: delka2

!   integer*8 :: plan

   complex*16, dimension(-Nxr2:Nxr2-1) :: array

   delka2=delka*2d0

!   denmat2=cmplx(1d0,0d0,8)

!   call dfftw_plan_dft_1d(plan,Nxr,array,array, FFTW_FORWARD, FFTW_MEASURE)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ARRAY,IKA,IKA2,IXR,SGNFAC)

   !do even ika first
!$OMP DO
   do ixa=0,Nxa2-1
    !construct array to transform
    array=cmplx(0d0,0d0,8)
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
    array=cmplx(0d0,0d0,8)
    sgnfac=1
    do ixr=-Nxr2,Nxr2-1
     array(ixr)=(denmat(ixa,ixr)-denmat(ixa-Nxa2,ixr)) &
                *exp(-imagi*delxr*delka2*0.5d0*ixr) &
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
!     if(abs(dble(denmat(ixa,ika)))<1d-40) then
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
   real*8  :: delkr2

   complex*16, dimension(-Nxa2:Nxa2-1) :: array

   delkr2=delkr*2d0

!   denmat2=cmplx(1d0,0d0,8)

   !do even ikr first
   do ika2=-Nka2,Nka2-1
    ika=ika2*2
    !construct array to transform
    array=cmplx(0d0,0d0,8)
    do ixa=-Nxa2,Nxa2-1
     array(ixa)=denmat(ixa,ika)
    enddo !ixa

    !transform!
    do ikr2=-Nkr2/2,Nkr2/2-1
     ikr=ikr2*2
     denmat2(ikr,ika)=cmplx(0d0,0d0,8)
!     val=cmplx(0d0,0d0,8)
     do ixa=-Nxa2,Nxa2-1
      denmat2(ikr,ika)=denmat2(ikr,ika)+array(ixa) &
                                        *exp(-imagi*delxa*delkr2*ixa*ikr2)
     enddo !ixa
    enddo !ikr2
   enddo !ika2

   !do odd ikr,ika now
   do ika2=-Nka2,Nka2-1
    ika=ika2*2+1
    !construct array to transform
    array=cmplx(0d0,0d0,8)
    do ixa=-Nxa2,Nxa2-1
     array(ixa)=denmat(ixa,ika)*exp(-imagi*delxa*delkr2*0.5d0*ixa)
    enddo !ixa

    !transform!
    do ikr2=-Nkr2/2,Nkr2/2-1
     ikr=ikr2*2+1
     denmat2(ikr,ika)=cmplx(0d0,0d0,8)
!     val=cmplx(0d0,0d0,8)
     do ixa=-Nxa2,Nxa2-1
      denmat2(ikr,ika)=denmat2(ikr,ika)+array(ixa) &
                                        *exp(-imagi*delxa*delkr2*ixa*ikr2)
     enddo !ixa
    enddo !ikr2
   enddo !ika2

   denmat=denmat2*delxa*invsqrt2pi

   denState=MOMENTUM

  end subroutine transform_w_to_k_norepeat
 
END MODULE mesh

