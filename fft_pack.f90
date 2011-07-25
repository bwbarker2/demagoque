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

!! transforms density matrix from physical to spectral space using the
!! FFTPACK5 subroutines
subroutine FT(L,M,xre,xim)
  use mesh
  use prec_def
  implicit none

  integer :: L,M !,iikr0,ika   !dimensions of matrix
  real (Long) :: xre(L,M),xim(L,M)

  IF(isDenK)THEN
     write(6,*)'FT: Error: Tried to FT on spectral density matrix'
     RETURN
  ENDIF

  call FFT2C(L,M,xre,xim,1)

  isDenX=.false.
  isDenK=.true.

! this is a horrible idea, don't do it, BWB 2010-07-21
!  iikr0=iNkr2(0)
!  do ika=-Nka2,Nka2
!   den_im(iikr0,iNka2(ika))=0.0
!  enddo

!  WRITE(*,*)'FT finished'

end subroutine FT


!!transforms matrix from spectral space to physical space
subroutine IFT(L,M,xre,xim)
  use mesh
  use prec_def
  implicit none

  integer :: L,M !,ixa,iixr0   !dimensions of matrix
  real (Long) :: xre(L,M),xim(L,M)

  IF(isDenX)THEN
     write(6,*)'IFT: Error: Tried to IFT on spatial density matrix'
     RETURN
  ENDIF

  call FFT2C(L,M,xre,xim,2)

  isDenK=.false.
  isDenX=.true.

!  write(*,*)'IFT finished'

! this is a horrible idea, don't do it, BWB 2010-07-21
!   iixr0=iNxr2(0)
!  do ixa=-Nxa2,Nxa2
!   den_im(iNxa2(ixa),iixr0)=0.0
!  enddo
 
end subroutine IFT


!! 2D complex FFT, either backward or forward
subroutine FFT2C(L,M,xre,xim,fb)
  use prec_def
  implicit none

  integer :: L,M   !dimensions of matrix
  real (Long) :: xre(L,M),xim(L,M)
  
  ! variables for fft
  integer :: fb    !1 - forward, 2 - backward transform
  integer :: LENSAV,LENWRK,IER           !length of work array, error flag
  real (Long), allocatable :: WSAVE(:), WORK(:)  !work arrays
  complex (Long) :: cmat(L,M)   !complex matrix
!  external zfft2b


  LENSAV = 2*(L+M) + INT(LOG(REAL(L))) + INT(LOG(REAL(M))) + 8
!  lensav = lensav*2 !test
  LENWRK = 2*L*M
!  lenwrk = lenwrk*2 !test
  allocate(WSAVE(1:LENSAV))
  allocate(WORK(1:LENWRK))

  !fill complex matrix
  cmat=CMPLX(xre,xim)

  !initialize work array
  call ZFFT2I(L,M,WSAVE,LENSAV,IER)

  if(fb.EQ.1)then
     call ZFFT2F(L, L, M, cmat, WSAVE, LENSAV, WORK, LENWRK, IER)
     !cmat=cmat/dsqrt(1.0_Long*Nxr*Nxa)
  elseif(fb.EQ.2)then
     call zfft2b(L, L, M, cmat, WSAVE, LENSAV, WORK, LENWRK, IER)
  else
     write(6,*)'fft_pack: bad option passed: fb=',fb
  endif

  !copy back to xre,xim
  xre=REAL(cmat)
  xim=AIMAG(cmat)

!  WHERE(xre.LE.1e-10_long.AND.xre.GT.-1e-10_long)
!     xre=0
!  ENDWHERE
!  WHERE(xim.LE.1e-10_long.AND.xim.GT.-1e-10_long)
!     xim=0
!  ENDWHERE

end subroutine FFT2C



subroutine FFT1(L,M,xre,xim,fb)
!! FFT1 - computes 2D complex FFT using real sine, cosine transforms
  use fft_vars
  use prec_def
  implicit none

  integer :: io,in,o    !loop variables, o - index of 'other' dimension
  integer :: L,M   !dimensions of matrix
  real (Long) :: xre(L,M),xim(L,M)
  
  integer :: N(2) !dimensions of matrix, N(1)=L, N(2)=M
  integer :: fb    !1 - forward, 2 - backward transform
  integer :: IER=0                         !error flag
  real (Long),allocatable :: vecre(:), vecim(:)   !working vectors
  real (Long),allocatable :: sinim(:), sinre(:), &  !sine transforms
                             cosim(:), cosre(:)     !cosine transforms

  N=(/ L, M /)

  !if work arrays not computed, do so now
  if(.not.allocated(wsaves))then
     call fft_initial(N)
!     write(*,*)'allocating work arrays'
  endif

  !do columns, then rows
  do in=1,2
     allocate(vecre(N(in)), vecim(N(in)))
     allocate(sinim(N(in)), sinre(N(in)), cosre(N(in)), cosim(N(in)))
     
     !find index of N() that we are not in ('other' dim)
     o=2/in
     
     ! loop over columns or rows
     DO io=1,N(o)
        IF(io.EQ.1)THEN
           cosre=xre(io,:)
           cosim=xim(io,:)
        ELSE
           cosre=xre(:,io)
           cosim=xim(:,io)
        ENDIF
        sinre=cosre
        sinim=cosim
        CALL DCOST1F(N(in),1,cosre,N(in),wsavec(in,:),lensav(in),work(in,:),lenwrk(in),ier)
        call DCOST1F(N(in),1,cosim,N(in),wsavec(in,:),lensav(in),work(in,:),lenwrk(in),ier)
        call DSINT1F(N(in),1,sinre,N(in),wsaves(in,:),lensav(in),work(in,:),lenwrk(in),ier)
        call DSINT1F(N(in),1,sinim,N(in),wsaves(in,:),lensav(in),work(in,:),lenwrk(in),ier)
        
        !if inverse transform, negate both sine transforms
        if(fb.EQ.2)THEN
           sinre=sinre*(-1.0)
           sinim=sinim*(-1.0)
        ENDIF

        IF(io.EQ.1)THEN
           xre(io,:)=cosre + sinim
           xim(io,:)=cosim - sinre
        ELSE
           xre(:,io)=cosre + sinim
           xim(:,io)=cosim - sinre
        ENDIF
        !     cmat(il,:)=cvec
        !
        !     do im=1,M
        !        if(real(cvec(im)).GT.1.OR.real(cvec(im)).LT.-1)write(*,*)im,real(cvec(im))
        !     enddo
     enddo
     
     deallocate(vecre,vecim,sinim, sinre, cosre, cosim)

     

  enddo

end subroutine FFT1

subroutine fft_initial(N)
  !!fft_initial - initialize work arrays for FFT
  use fft_vars
  implicit none

  integer :: N(2), in, IER, savmx, wrkmx

  do in=1,2
     LENSAV(in) = 2*N(in) + INT(LOG(REAL(N(in)))) + 4
     LENWRK(in) = 2*N(in)+2
  enddo

  savmx=MAX(lensav(1),lensav(2))
  wrkmx=MAX(lenwrk(1),lenwrk(2))
  
  allocate(WSAVEC(2,savmx))
  allocate(WSAVES(2,savmx))

  allocate(WORK(2,wrkmx))
  
  do in=1,2
     call DCOST1I(N(in),WSAVEC(in,:),LENSAV(in),IER)
     call DSINT1I(N(in),WSAVES(in,:),LENSAV(in),IER)
  enddo
     
end subroutine fft_initial
