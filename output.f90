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

SUBROUTINE output
  use cons_laws
  USE mesh
  USE out
  USE time
  IMPLICIT NONE

!  INTEGER :: iixa,iixr,ixr,ii

  ! write output only if we are at beginning, end, or multiple of ntime
  IF(MOD(it,ntime).EQ.0.OR.it.EQ.Nt)THEN
     CONTINUE
  ELSE
     RETURN
  ENDIF

  write(*,*)'running step',it


!  write(*,*)'starting IF(isDenK)'

  ! decide order of output based on whether density matrix is in x or k space
  IF(denState==MOMENTUM)THEN
!     write(*,*)'starting outK'
     CALL outK
!     write(*,*)'ending outK'
     call outW
     CALL outX
  ELSE
     CALL outX
     call outW
     CALL outK
  ENDIF

  write(*,*)ekin+epot

  ! write factor that kinetic must be multiplied by to conserve energy
  write(*,*)'timestep,kinetic_correction_factor:',it,(ep0+ek0-epot)/ekin
  !write(*,*)'timestep,kinetic_correction_factor:',it,(epot-ep0)/(ek0-ekin)

 if(firstOutput)firstOutput=.false.

END SUBROUTINE output


SUBROUTINE outX
  use mesh
  use osc_pars
  use time
  IMPLICIT NONE

  call setState(SPACE)

  CALL outDiagX
  CALL outDenMat(61,62)
  call ener_x

  ! output analytic oscillator to compare with numeric
  if(potFinal==0.AND.Nmax==0.AND.EA<1d-5) call outAnalHarmonic
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
  CALL outDiagK
  CALL ener_k
  CALL outEner

END SUBROUTINE outK

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE outDenMat(fileim_u, filere_u)
  !! outDenMat - writes real, imaginary density matrices to file unit specified by inputs
  USE mesh
  USE time
  IMPLICIT NONE

  INTEGER, INTENT(in) :: filere_u, fileim_u  !file units of real, imaginary density matrices
  INTEGER :: ixa,ixr

  WRITE(fileim_u,*)'# time=',t,'fm/c'
  WRITE(filere_u,*)'# time=',t,'fm/c'

  DO ixr=-Nxr,Nxr-1
     WRITE(fileim_u,919) (DIMAG(getDen(ixa,ixr)),ixa=-Nxa2,Nxa2-1)
     WRITE(filere_u,919) (DBLE(getDen(ixa,ixr)),ixa=-Nxa2,Nxa2-1)
  ENDDO

  ! two blank lines give proper formatting for gnuplot
  WRITE(fileim_u,*)
  WRITE(fileim_u,*)
  WRITE(filere_u,*)
  WRITE(filere_u,*)


919 FORMAT(1000es23.14e3)

END SUBROUTINE outDenMat


SUBROUTINE outDiagK
  !! outDiagK - writes the value of the density matrices along the k=k' line.
  USE mesh
  USE prec_def
  USE time
  IMPLICIT NONE

  INTEGER :: ika 
  REAL*8 :: ddre, dddim

  WRITE(42,*)'# time=',t,'fm/c'
  WRITE(42,*)'# k [fm], real, imaginary amplitudes'

  DO ika=-Nka2,Nka2-1
     ddre=DBLE(getDenK(0,ika))
     dddim=DIMAG(getDenK(0,ika))
     WRITE(42,*) ka(ika),ddre,dddim
  ENDDO

  WRITE(42,*)
  WRITE(42,*)

!93 format(3e16.8)

END SUBROUTINE outDiagK


SUBROUTINE outDiagX
  USE mesh
  USE time
  IMPLICIT NONE

  INTEGER :: ixa
  REAL (Long) :: ddre, dddim

  WRITE(41,*)'# time=',t,'fm/c'
  WRITE(41,*)'# x [fm], real, imaginary amplitudes'

  DO ixa=-Nxa2,Nxa2-1
     ddre=DBLE(getDen(ixa,0))
     dddim=DIMAG(getDen(ixa,0))
!     write(*,*)'den_im,dddim:',den_im(iNxa2(ixa),iixr0),dddim
     WRITE(41,*) xa(ixa),ddre,dddim
  ENDDO

  WRITE(41,*)
  WRITE(41,*)

93 format(3e16.8)

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
     write(44,*)'# time     ekin     ekerr     epot     eperr'

  endif
  
  !write(*,*)ekin,ek0
  relekin=(ekin-ek0)/ek0
  relepot=(epot-ep0)/ep0
  !propagated error - fixed 2010-07-09
  relekerr=(ekerr/ek0)**2+(ek0err/ek0*(1+relekin))**2
  relekerr=dsqrt(relekerr)
  releperr=(eperr/ep0)**2+(ep0err/ep0*(1+relepot))**2
  releperr=dsqrt(releperr)

  write(43,*)t,relekin,relekerr,relepot,releperr
  write(44,*)t,ekin,ekerr,epot,eperr
  
END SUBROUTINE outEner

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine outDenUnf
 !! outDenUnf - writes the density matrix to the unformatted file obdm.dat
 use mesh
 use osc_pars
 implicit none

 integer :: ixa,ixr

 write(54) Nxa,Nxr,Nmax

 do ixa=-Nxa2,Nxa2
  do ixr=-Nxr,Nxr-1
   write(54)denmat(ixa,ixr)
  enddo
 enddo

end subroutine outDenUnf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine inDenUnf
 !! inDenUnf - reads the density matrix that is written in previous subroutine
 use mesh
 use osc_pars
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

 do ixa=-Nxa2,Nxa2
  do ixr=-Nxr,Nxr-1
   read(54)denmat(ixa,ixr)
  enddo
 enddo

end subroutine inDenUnf
