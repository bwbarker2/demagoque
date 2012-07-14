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

  if(firstOutput) call output_ModeFilesOpen

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

 if(lastOutput) call output_ModeFilesClose

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
  call outAnalyticXDiff
  call outAnalyticXDiffness

  call outEigens

  call outDenUnf

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
 use phys_cons
 use time
 implicit none

  integer, save :: analyticx_diag  !< diagonal of analytic wavefunction output

  integer :: ixa

  if(firstOutput)then
  open(newunit=analyticx_diag, file=char(fout_ev_pre)//'analyticx_diag.dat')
  write(analyticx_diag,*)'# xx     probability'
 endif

 write(analyticx_diag,*)'# time = ',char(phys_cons_unit_abbrev_time)

 do ixa=Nxan,Nxax
  write(analyticx_diag,*)xa(ixa),abs(initSuperWavefunction%getWavefn(xa(ixa),t))**2 
 enddo

 write(analyticx_diag,*)
 write(analyticx_diag,*)

end subroutine outAnalyticX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> outputs data of analytic minus density matrix evolution
subroutine outAnalyticXDiff
 use formatting
 use input_parameters
 use mesh
 use phys_cons
 use time
 implicit none

  integer, save :: file_diff_diag  !< diagonal of analytic - denmat evolution output

  integer :: ixa

  if(firstOutput)then
  open(newunit=file_diff_diag, file=char(fout_ev_pre)//'analyticXDiff_diag.dat')
  write(file_diff_diag,*)'# xx     analytic minus denmat evolution density'
 endif

 write(file_diff_diag,*)'# time = ',char(phys_cons_unit_abbrev_time)

 do ixa=Nxan,Nxax
  write(file_diff_diag,*)xa(ixa) &
    ,abs(initSuperWavefunction%getWavefn(xa(ixa),t))**2 &
     -getDenDiagX(ixa)
 enddo

 write(file_diff_diag,*)
 write(file_diff_diag,*)

end subroutine outAnalyticXDiff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine outAnalyticXDiffness
 use formatting
 use input_parameters
 use mesh
 use phys_cons
 use time
 implicit none

  integer, save :: file_diff_diag  !< diagonal of analytic - denmat evolution output

  integer :: ixa,ixr

  real (Long) :: sumxre,sumxim,sum2xre,sum2xim,totaln
  complex (Long) :: diff

  if(firstOutput)then
  open(newunit=file_diff_diag, file=char(fout_ev_pre)//'analyticXDiffness.dat')
  write(file_diff_diag,*) &
   '# time   ave-diff-re   ave-diff-im   rms-diff-re rms-diff-im'
 endif

! write(file_diff_diag,*)'# time = ',char(phys_cons_unit_abbrev_time)

 sumxre = 0._Long
 sum2xre = 0._Long
 sumxim = 0._Long
 sum2xim = 0._Long

 do ixa=Nxan,Nxax
  do ixr=Nxrn,Nxrx
   diff = conjg(initSuperWavefunction%getWavefn(xx1(ixa,ixr),t)) &
            * initSuperWavefunction%getWavefn(xx2(ixa,ixr),t) &
          - getDenX(ixa,ixr)
   sumxre=sumxre+real(diff)
   sumxim=sumxim+aimag(diff)
   sum2xre=sum2xre+(real(diff))**2
   sum2xim=sum2xim+(aimag(diff))**2
  enddo
 enddo

  totaln=1._Long/((Nxax-Nxan+1)*(Nxrx-Nxrn+1))

  sumxre=sumxre*totaln
  sumxim=sumxim*totaln
  sum2xre=sqrt(sum2xre*totaln)
  sum2xim=sqrt(sum2xim*totaln)

 write(file_diff_diag,*)t,sumxre,sumxim,sum2xre,sum2xim

end subroutine outAnalyticXDiffness


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine outEigens
 use bexception
 use formatting, only : fout_ev_pre
 use iso_varying_string, only : char
 use lib_lapack
 use mesh
 use phys_cons, only : phys_cons_unit_abbrev_time
 use time
 implicit none

 integer, save :: filename ! unit of output file
 integer :: ii

 complex(Long), dimension(Nxax-Nxan+1) :: evals
 complex(Long), dimension(Nxax-Nxan+1,Nxax-Nxan+1) :: evecs

 if(Nxax-Nxan/=Nxrx-Nxrn) then
  call throwException('outEigens: Matrix not square, not outputting eigenvalues', BEXCEPTION_WARNING)
  return
 endif

 call getEigenSq(denmat,Nxax-Nxan+1,evals,evecs)

 if(firstOutput) then
  open(newunit=filename, file=char(fout_ev_pre)//'eigens.dat')
  write(filename,*)'# List of eigenvalues at each time'
  write(filename,*)'# real   imaginary'
 endif

 write(filename,*)'# time =',t,char(phys_cons_unit_abbrev_time)

 do ii=1,size(evals)
  write(filename,*)ii,real(evals(ii)),aimag(evals(ii))
 enddo

 write(filename,*)
 write(filename,*)

 if(lastOutput) close(filename)

end subroutine outEigens


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
 use formatting
 use input_parameters
 use mesh
 use time
 implicit none

 integer :: ixa,ixr,fileu
 
 open(newunit=fileu, form='unformatted', &
      file=char(fout_ev_pre//'ufo_'//time_getString(it/ntime)//'.dat'))


! write(fileu) Nxa,Nxr,Nmax

 do ixa=Nxan,Nxax
  do ixr=Nxrn,Nxrx
   write(fileu)denmat(ixa,ixr)
  enddo
 enddo

 close(fileu)

end subroutine outDenUnf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine output_ModeFilesClose
 implicit none

 close(41)
 close(42)
 close(43)
 close(44)
 close(45)
 close(61)
 close(62)
 close(66)
 close(67)
 close(68)
 close(69)
 close(70)
 close(72)
 close(73)
 close(74)
 close(75)
 
end subroutine output_ModeFilesClose

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_ModeFilesOpen
 use formatting
 implicit none

 open(unit=41,file=char(fout_ev_pre)//'denmat_x_t.dat')
 open(unit=42,file=char(fout_ev_pre)//'denmat_k_t.dat')
 open(unit=43,file=char(fout_ev_pre)//'cons_rel.dat')
 open(unit=44,file=char(fout_ev_pre)//'cons_abs.dat')
 open(unit=45,file=char(fout_ev_pre)//'denmatan_x_t.dat') 
 OPEN(unit=61,file=char(fout_ev_pre)//'2dxim.dat')
 OPEN(unit=62,file=char(fout_ev_pre)//'2dxre.dat')
 OPEN(unit=66,file=char(fout_ev_pre)//'2dkim.dat')
 OPEN(unit=67,file=char(fout_ev_pre)//'2dkre.dat')
 open(unit=68,file=char(fout_ev_pre)//'2dwim.dat')
 open(unit=69,file=char(fout_ev_pre)//'2dwre.dat')
 OPEN(unit=70,file=char(fout_ev_pre)//'pk2.dat')
 open(unit=72,file=char(fout_ev_pre)//'2dx.dat')
 open(unit=73,file=char(fout_ev_pre)//'2dw.dat')
 open(unit=74,file=char(fout_ev_pre)//'2dk.dat')
 open(unit=75,file=char(fout_ev_pre)//'mean_abs_x.dat')

end subroutine output_ModeFilesOpen

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

subroutine inDenUnf(filename,indenmat)
 !! inDenUnf - reads the density matrix that is written in previous subroutine
 use input_parameters
 use mesh
 implicit none

 character(len=*), intent(in) :: filename !< name of file to read from
 complex (Long), dimension(Nxan:Nxax,Nxrn:Nxrx), intent(out) &
                :: indenmat !< matrix to save to

 integer :: ixa,ixr,fileu

 open(newunit=fileu, form='unformatted', file=filename)

 do ixa=Nxan,Nxax
  do ixr=Nxrn,Nxrx
   read(fileu)indenmat(ixa,ixr)
  enddo
 enddo

 close(fileu)

end subroutine inDenUnf
