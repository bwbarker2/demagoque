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

!> outputs various snapshot files at certain timesteps.
SUBROUTINE output
  use cons_laws
  use input_parameters
  USE mesh
  USE time
  IMPLICIT NONE

!  INTEGER :: iixa,iixr,ixr,ii

  ! write output only if we are at beginning, end, or multiple of ntime
  IF(MOD(it,ntime).EQ.0.OR.it.EQ.Nt)THEN
     CONTINUE
  ELSE
     RETURN
  ENDIF

  if(it==Nt)lastOutput=.true.

  write(*,*)'running step',it

!  call mesh_setReflectedLR(.false.)

!  write(*,*)'starting IF(isDenK)'

  ! decide order of output based on whether density matrix is in x or k space
  IF(denState==MOMENTUM)THEN
!     write(*,*)'starting outK'
     CALL outK
!     write(*,*)'ending outK'
!     if(.not.useFrameXXP) call outW
     CALL outX
  ELSE
     CALL outX
!     if(.not.useFrameXXP) call outW
     CALL outK
  ENDIF

!  write(*,*)ekin+epot

  ! write factor that kinetic must be multiplied by to conserve energy
  write(*,*)'timestep,kinetic_correction_factor:',it,(ep0+ek0-epot)/ekin
  !write(*,*)'timestep,kinetic_correction_factor:',it,(epot-ep0)/(ek0-ekin)

 if(firstOutput)firstOutput=.false.
 if(lastOutput)lastOutput=.false.

! call mesh_setReflectedLR(.true.)

END SUBROUTINE output

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE outX
  use input_parameters
  use mesh
  IMPLICIT NONE

  call setState(SPACE)
  CALL outDenMat(61,62)
  CALL outDiagX
  call setState(SPACE)
  call outDenMatXPhys
  call outSpikinessX
  call ener_x

  call outAnalyticX

  ! output analytic oscillator to compare with numeric
  if(potFinal==-1.AND.Nmax==0.AND.EA<1d-5) call outAnalHarmonic
!  call howHermitian

END SUBROUTINE outX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine outW
 use mesh
 implicit none

 call setState(WIGNER)

 call outDenMat(68,69)

end subroutine outW

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE outK
  use mesh
  IMPLICIT NONE

  call setState(MOMENTUM)
  CALL outDenMat(66,67)
  call outDenMatKPhys
  CALL outDiagK
  CALL ener_k
  CALL outEner

END SUBROUTINE outK

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine outAnalyticX
 use formatting
 use input_parameters
 use mesh
 use time
 implicit none

  integer, save :: analyticx_diag

  integer :: ixa

  if(firstOutput)then
  open(newunit=analyticx_diag, file=char(fout_ev_pre)//'analyticx_diag.dat')
  write(analyticx_diag,*)'# xx     probability'
 endif

 do ixa=Nxan,Nxax
  write(analyticx_diag,*)xa(ixa),abs(initSuperWavefunction%getWavefn(xa(ixa),t))**2 
 enddo

 write(analyticx_diag,*)
 write(analyticx_diag,*)

end subroutine outAnalyticX


SUBROUTINE outDenMat(fileim_u, filere_u)
  !! outDenMat - writes real, imaginary density matrices to file unit specified by inputs
  USE mesh
  USE time
  IMPLICIT NONE

  INTEGER, INTENT(in) :: filere_u, fileim_u  !file units of real, imaginary density matrices
  INTEGER :: ixa,ixr

  WRITE(fileim_u,*)'# time=',t,'fm/c'
  WRITE(filere_u,*)'# time=',t,'fm/c'

  DO ixr=Nxrn,Nxrx
     WRITE(fileim_u,919) (AIMAG(getDen(ixa,ixr)),ixa=Nxan,Nxax)
     WRITE(filere_u,919) (REAL(getDen(ixa,ixr)),ixa=Nxan,Nxax)
  ENDDO

  ! two blank lines give proper formatting for gnuplot
  WRITE(fileim_u,*)
  WRITE(fileim_u,*)
  WRITE(filere_u,*)
  WRITE(filere_u,*)


919 FORMAT(1000es23.14e3)

END SUBROUTINE outDenMat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine outDenMatKPhys()
 use input_parameters
 use mesh
 use prec_def
 use time
 implicit none

 integer, parameter :: funit = 74 ! file unit to write to

 integer :: ikr,ika

 call setState(MOMENTUM)

 write(funit,*)'# time=',t,'fm/c'

 do ikr=Nkrn,Nkrx
  do ika=Nkan,Nkax
   if(.not.useFrameXXP) then
    if(abs(mod(ika+ikr,2))==1) cycle
   endif
   write(funit,*)kr(ikr),ka(ika),0.5_Long*REAL(getDenK(ikr,ika)),0.5_Long*AIMAG(getDenK(ikr,ika))
  enddo
  write(funit,*)
 enddo

 write(funit,*)
 write(funit,*)

end subroutine outDenMatKPhys

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine outDenMatXPhys()
 use input_parameters
 use mesh
 use prec_def
 use time
 implicit none

 integer, parameter :: funit = 72 ! file unit to write to

 integer :: ixa,ixr,ixrl,ixru

 call setState(SPACE)

 write(funit,*)'# time=',t,'fm/c'

 !set bounds for ixr loop. Rotated frame is double Nxr
 if(useFrameXXP) then
  ixrl=Nxrn
  ixru=Nxrx
 else !if(useMeshXAR2) then
  ixrl=Nxrn+Nxr2
  ixru=Nxrx-Nxr2
 endif

 do ixa=Nxan,Nxax
  do ixr=ixrl,ixru
   if(useMeshXAR2) then
    if(abs(mod(ixa+ixr,2))==1) cycle
   endif
   write(funit,*)xa(ixa),xr(ixr),REAL(getDenX(ixa,ixr)),AIMAG(getDenX(ixa,ixr))
  enddo
  write(funit,*)
 enddo

 write(funit,*)
 write(funit,*)

end subroutine outDenMatXPhys

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE outDiagK
  !! outDiagK - writes the value of the density matrices along the k=k' line.
  USE mesh
  USE prec_def
  USE time
  IMPLICIT NONE

  INTEGER :: ika 
  REAL(kind=8) :: ddre  !, dddim

  WRITE(42,*)'# time=',t,'fm/c'
  WRITE(42,*)'# k [fm], real, imaginary amplitudes'

  DO ika=Nkan,Nkax
     ddre=getDenDiagK(ika)
!     dddim=AIMAG(getDenK(0,ika))
     WRITE(42,*) ka(ika),ddre
  ENDDO

  WRITE(42,*)
  WRITE(42,*)

!93 format(3e16.8)

END SUBROUTINE outDiagK

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE outDiagX
  USE mesh
  USE time
  IMPLICIT NONE

  INTEGER :: ixa
  REAL (Long) :: ddre  !, dddim

  real (Long) :: sum1, sum2

  WRITE(41,*)'# time=',t,'fm/c'
  WRITE(41,*)'# x [fm], real, imaginary amplitudes'

  DO ixa=Nxan,Nxax
     ddre=getDenDiagX(ixa)
!     dddim=AIMAG(getDen(ixa,0))
!     write(*,*)'den_im,dddim:',den_im(iNxa2(ixa),iixr0),dddim
     WRITE(41,*) xa(ixa),ddre  !,dddim
  ENDDO

  WRITE(41,*)
  WRITE(41,*)

!93 format(3e16.8)

  ! This computes and outputs <|x|>

  sum1=0_Long
  sum2=0_Long

  do ixa=Nxan,Nxax
   sum1=sum1+getDenDiagX(ixa)
   sum2=sum2+abs(xa(ixa))*getDenDiagX(ixa)
  enddo

  if(firstOutput) write(75,*)'# time     <|x|>'
  
  write(75,*)t,sum2/sum1
  write(*,*)'mean_abs_x=',sum2/sum1

END SUBROUTINE outDiagX

SUBROUTINE outEner
  use cons_laws
  use formatting
  use prec_def
  use time
  implicit none

  real (Long) :: relekin, relekerr   ! relative kinetic energy
  real (Long) :: relepot, releperr   ! relative potential energy

  if(firstOutput)then

     ek0=ekin
     ek0err=ekerr
     ep0=epot
     ep0err=eperr

     write(43,*)'# Relative conserved quantities over time'
     write(43,*)'# form: (q-q0)/q0 [MeV]'
     write(43,*)'# time     ekin     ekerr     epot     eperr'

     write(44,*)'# Conserved quantities over time [MeV]'
     write(44,*)'# time     ekin     ekerr     epot     eperr     nnum     knum'

  endif
  
  !write(*,*)ekin,ek0
  relekin=(ekin-ek0)/ek0
  relepot=(epot-ep0)/ep0
  !propagated error - fixed 2010-07-09
  relekerr=(ekerr/ek0)**2+(ek0err/ek0*(1e0_Long+relekin))**2
  relekerr=sqrt(relekerr)
  releperr=(eperr/ep0)**2+(ep0err/ep0*(1e0_Long+relepot))**2
  releperr=sqrt(releperr)

  write(43,*)t,relekin,relekerr,relepot,releperr
  write(44,*)t,ekin,ekerr,epot,eperr,nnum,knum
  
END SUBROUTINE outEner

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine outDenUnf
 !! outDenUnf - writes the density matrix to the unformatted file obdm.dat
 use input_parameters
 use mesh
 implicit none

 integer :: ixa,ixr

 write(54) Nxa,Nxr,Nmax

 do ixa=Nxan,Nxax
  do ixr=Nxrn,Nxrx
   write(54)denmat(ixa,ixr)
  enddo
 enddo

end subroutine outDenUnf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Calculates and outputs an observable I call spikiness.
!!
!! \f[\textrm{spikiness}=\frac{1}{N}\sum_i^N \left|\frac{y(i+1)-y(i)}{x(i+1)-x(i)}\right|\f]
!!
!! This gives the average slope between two adjacent points on the diagonal. Maybe this should be over whole density matrix? Should I care about slope, or just difference in height?
subroutine outSpikinessX
 use formatting
 use mesh
 use time
 implicit none

 integer :: ixa
 integer, save :: spikifile   !file unit number

 real (Long) :: spike

 if(firstOutput)then
  open(newunit=spikifile, file=char(fout_ev_pre)//'spikiness.dat')
  write(spikifile,*)'# time     spikiness'
 endif

 !need to handle last point cyclically
 do ixa=Nxan,Nxax-1
  spike=spike+abs((getDenDiagX(ixa)-getDenDiagX(ixa+1))/delxa)
 enddo

 spike=spike+abs((getDenDiagX(Nxax)-getDenDiagX(Nxan))/delxa)

 spike=spike/Nxa

 write(spikifile,*)t,spike

 if(lastOutput) then
  close(spikifile)
 endif

end subroutine outSpikinessX


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine inDenUnf
 !! inDenUnf - reads the density matrix that is written in previous subroutine
 use input_parameters
 use mesh
 implicit none

 integer :: ixa,ixr,Nxa1,Nxr1,Nmax1

 read(54) Nxa1,Nxr1,Nmax1

 if(Nxa1.NE.Nxa) then
  write(*,*) 'Nxa=',Nxa,'Nxa from obdm=',Nxa1
  write(*,*) 'RUN ADIABATIC SWITCHING?'
  stop
 endif

 if(Nxr1.NE.Nxr) then
  write(*,*) 'Nxr=',Nxr,'Nxr from obdm=',Nxr1
  write(*,*) 'RUN ADIABATIC SWITCHING?'
  stop
 endif

 if(Nmax1.ne.Nmax) then
  write(*,*) 'Nmaxx=',Nmax,'Nmaxx from obdm=',Nmax1
  write(*,*) 'RUN ADIABATIC SWITCHING?'
  stop
 endif

 do ixa=Nxan,Nxax
  do ixr=Nxrn,Nxrx
   read(54)denmat(ixa,ixr)
  enddo
 enddo

end subroutine inDenUnf
