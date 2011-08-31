module skyrme_params
 implicit none

 !mean field interaction parameters
 ! Vautherin, Treiner, Veneroni PB191,6(1987)
! real*8, parameter :: t0  = -2164.1d0
! real*8, parameter :: t3  = 15079.0d0
! real*8, parameter :: sig = 0.25d0

 !Rios variational parameters (give saturation energy 0.16 MeV, saturation density 0.16 fm-3, compressibility 220 MeV)
 real*8, parameter :: t0  = -2150.1d0
 real*8, parameter :: t3  = 14562d0
 real*8, parameter :: sig = 0.257d0
 
end module skyrme_params
