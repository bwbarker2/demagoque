module cons_laws
  use prec_def
  implicit none

  real (Long) :: ekin   ! total energy calculated in k-space
  real (Long) :: ekerr  ! uncertainty of above
  real (Long) :: ek0    ! initial energy calc'd in k-space
  real (Long) :: ek0err ! error in initial

  real (Long), allocatable :: potx(:) ! potential on diagonal
  real (Long) :: epot   ! total potential energy in x-space
  real (Long) :: eperr  ! uncertainty of above
  real (Long) :: ep0    ! initial pot energy
  real (Long) :: ep0err ! error in initial

  real (Long) :: nnum   ! number of particles/norm_thy calc'd in pos space
  real (Long) :: knum   ! number of particles/norm_thy calc'd in mom space

end module cons_laws


