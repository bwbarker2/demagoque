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
   allocate(denmat(-Nxa2:Nxa2,-Nxr:Nxr)) !2x size in xr for naive FT - BWB 2011-01-10
   allocate(denmat2(-Nxa2:Nxa2-1,-Nxr:Nxr-1))
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
  
   if(denState.NE.state) then
    select case (denState)
  
     case (SPACE)
      select case (state)
       case (WIGNER)
        call transform_x_to_wigner_dumb
       case (MOMENTUM)
        call transform_x_to_k_norepeat
!        call transform_x_to_wigner_dumb
!        call transform_wigner_to_k_dumb
      end select
  
     case (WIGNER)
      select case (state)
       case (SPACE)
        call transform_wigner_to_x_dumb
       case (MOMENTUM)
        call transform_wigner_to_k_fft_exp
      end select
  
     case (MOMENTUM)
      select case (state)
       case (WIGNER)
        call transform_k_to_wigner_fft_exp
       case (SPACE)
        call transform_k_to_wigner_fft_exp
        call transform_wigner_to_x_dumb
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

  subroutine transform_x_to_k_norepeat
   !! transform_x_to_k_nozeroes - transforms from coordinate to momentum space, with no redundancy from periodic extensions. See paper notes BWB 2011-03-11.
   use phys_cons
   implicit none

   integer    :: ikr,ikr2,ika,ika2,ixa,ixr
   real*8     :: delka2,delkr2, exparg
   complex*16 :: val

   real :: elapsed(2)

   delka2=2*delka
   delkr2=2*delkr

   !store in denmat2. Must zero it out first.
   denmat2=cmplx(0d0,0d0,8)

   write(*,*)'time at starting even k ft:',etime(elapsed)

   !first transform even ka and kr
   do ikr2=-Nkr2/2,Nkr2/2-1
    ikr=ikr2*2

    do ika2=-Nka2,Nka2-1
     ika=ika2*2

     val=cmplx(0d0,0d0,8)
     do ixa=-Nxa2,Nxa2-1
      do ixr=-Nxr2,Nxr2-1

       exparg=delxa*delkr2*ixa*ikr2 &
             +delxr*delka2*ixr*ika2
       val=val+getDenX(ixa,ixr)*exp(-imagi*exparg)

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
 

END MODULE mesh

