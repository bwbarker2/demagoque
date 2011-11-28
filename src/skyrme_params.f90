module skyrme_params
 use prec_def
 implicit none

 !mean field interaction parameters
 ! Vautherin, Treiner, Veneroni PB191,6(1987)
! real*8, parameter :: t0  = -2164.1d0
! real*8, parameter :: t3  = 15079.0d0
! real*8, parameter :: sig = 0.25d0

 !Rios variational parameters (give saturation energy 0.16 MeV, saturation density 0.16 fm-3, compressibility 220 MeV)
 real (Long), parameter :: t0  = -2150.1_Long
 real (Long), parameter :: t3  = 14562e0_Long
 real (Long), parameter :: sig = 0.257_Long
 
end module skyrme_params
