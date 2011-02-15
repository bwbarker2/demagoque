subroutine ener_k
  !! ener_k - calculates kinetic energy, total momentum from k-space den. mat.
  use cons_laws
  use mesh
  use osc_pars
  use phys_cons
  use time
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

  ekin=ekin*hbc*hbc/(2d0*m0)/(Nmax+1)*delka
  knum=knum*delka/(Nmax+1)

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
  use osc_pars
  use cons_laws
  use mesh
  implicit none

  real (Long) :: nerr
!  real (long) :: uoth(-Nxa2:Nxa2)
!  real (Long) :: epotOth,epotOtherr   ! other epot calc

  integer :: ixa
  
  call setState(SPACE)

!  call calcPotDiag()

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
  epot=epot*delxa/(Nmax+1)
  nnum=nnum*delxa/(Nmax+1)

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
