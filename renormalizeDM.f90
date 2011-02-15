subroutine renormalizeDM
 !! renormalizeDM - renormalizes density matrix. This should be run immediately after the output subroutine, so that the current number of particles has been calculated.
 use cons_laws
 use mesh
 implicit none

 integer :: ixa,ixr
 real*8  :: factor

 call setState(SPACE)

 factor=1d0/nnum

 do ixa=-Nxa2,Nxa2-1
  do ixr=-Nxr2,Nxr2-1

   call setDenX(ixa,ixr,factor*getDenX(ixa,ixr))

  enddo
 enddo

end subroutine renormalizeDM

