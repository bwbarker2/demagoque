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

subroutine ener_k
  !! ener_k - calculates kinetic energy, total momentum from k-space den. mat.
  use cons_laws
  use mesh
  use phys_cons
  IMPLICIT NONE

!  real (Long) :: pk2(-Nka2:Nka2)  !den.mat. times k^2
  real (Long) :: knum, eknum

  integer :: ika
  
  call setState(MOMENTUM)

  ekin=0d0
  knum=0d0
  do ika=-Nka2,Nka2-1
     ekin=ekin+ka(ika)*ka(ika)*DBLE(getDen(0,ika))
     knum=knum+DBLE(getDenK(0,ika))
!     write(*,*)pk2(ixa)
!     write(70,*)'ika,pk2',it,ika,pk2(ika)
  enddo

  ekin=ekin*hbc*hbc/(2d0*m0)/norm_thy*delka
  knum=knum*delka/norm_thy

!  call dint_simp1(Nxr+1, pk2, delka, ekin, ekerr)
!  call dint_simp1(Nxr+1, knn, delka, knum, eknum)

  ! Arnau has a reason for the last part of the next lines - BWB 2010-03-26
  ! I have a better reason to not have it - BWB 2010-09-01
!  ekin=ekin*hbc*hbc/(2.d0*m0)/(Nmax+1.d0)            !/delka*delxa/Nxr*Nxa*Nxr
!  knum=knum/(Nmax+1.d0)                        !/delka*delxa/Nxr*Nxa*Nxr
!  ekin=ekin/delka

!  write(*,*)'m0=',m0

!  ekerr=ekerr*hbc*hbc/(2.d0*m0)/(Nmax+1.d0)    !/delka*delxa*Nxa
  write(*,*)'ekin,ekerr=',ekin,ekerr
  write(*,*)'knum,eknum=',knum,eknum

end subroutine ener_k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ener_x
  !! ener_x - calculates total potential energy
  use cons_laws
  use mesh
  implicit none

  real (Long) :: nerr
!  real (long) :: uoth(-Nxa2:Nxa2)
!  real (Long) :: epotOth,epotOtherr   ! other epot calc

  integer :: ixa
  
  call setState(SPACE)

  call calcPotDiag()

  epot=0d0
  nnum=0d0
  do ixa=-Nxa2,Nxa2-1
     epot=epot+potDiag(ixa)*DBLE(getDenX(ixa,0))
!     uu(ixa)=uu(ixa)*den_re(iNxa2(ixa),iNxr2(0))/(Nmax+1)
     nnum=nnum+DBLE(getDenX(ixa,0))
    !write(*,*)uu(ixa)
!    uoth(ixa)=potMF(ixa)*den_re(iNxa2(ixa),iNxr2(0))/(Nmax+1)/(Nmax+1)/2
!    uoth(ixa)=potMF(ixa)*0.5/(Nmax+1)/(Nmax+1)/2
  enddo
  epot=epot*delxa/norm_thy
  nnum=nnum*delxa/norm_thy

!  call dint_simp1(Nxa+1, uoth, delxa, epotOth, epotOtherr)
!  call dint_simp1(Nxa, uu, delxa, epot, eperr)
!  call dint_simp1(Nxa+1, nn, delxa, nnum, nerr)

  ! prevent double-counting interaction pairs, normalize to
  ! potential per particle
!  epotOth=0.5*epotOth/(Nmax+1)

  write(*,*)'nnum,nerr=',nnum,nerr
  write(*,*)'epot,eperr=',epot,eperr 
!  write(*,*)'epotOth=',epotOth,epotOtherr

end subroutine ener_x