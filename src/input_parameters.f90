module input_parameters
 use prec_def
 implicit none

 integer :: potInitial    ! potential  with initial state
 integer :: potFinal      ! potential for time evolution

 real (Long) :: ea  !energy per particle [MeV]

 integer :: ntime ! write data every ntime timesteps

 REAL (Long) :: delt   ! timestep [fm/c]
 integer     :: Nevt   ! number of timesteps for evolution

 !! params_cutoff - parameters for imaginary off-diagonal cutoff. For explanation of formula, etc, see A. Rios et al., Annals of Physics 326 (2011) 1274, specifically page 1298. 
 logical     :: useImCutoff  ! use imaginary off-diagonal cutoff?
 real (Long) :: cutoff_w0    !strength of cutoff
 real (Long) :: cutoff_x0    !size of cutoff
 real (Long) :: cutoff_d0    !steepness of cutoff

 real (Long) :: initialSeparation ! 0 if not used, initial distance in fm between slabs

 logical     :: initState_gaussianNuclear
 REAL (long) :: w
 REAL (long) :: whm
 INTEGER     :: Nmax   ! maximum oscillator shell

 logical     :: initState_cosine !add sine wave to initial state
 integer     :: initState_cosine_number ! wave number n
 real (Long) :: initState_cosine_norm !normalization of sine wave, default 1
 real (Long) :: initState_cosine_shift !phase shift, default 0

 logical     :: initState_plane  !add plane wave to initial state
 integer     :: initState_plane_number
 real (Long) :: initState_plane_norm
 real (Long) :: initState_plane_shift

 logical     :: initState_kdelta  !add kronecker delta function to initial state
 real (Long) :: initState_kdelta_norm
 real (Long) :: initState_kdelta_x0  !center of delta function


! logical :: initState_randomG  !add a random number (Gaussian) to initial state
! real (Long) :: initState_random_norm  !norm
! real (Long) :: initState_random_fwhm  !full width half max

 integer     :: splitOperatorMethod !0 if not used, otherwise order of method

 logical     :: useImEvol   ! use imaginary evolution?
 integer     :: Nimev  ! number of timesteps for imaginary evolution
 
 logical     :: useFlipClone  ! symmetric collision method

 logical     :: useAdiabatic  ! if true, use adiabatic switching
 integer     :: iadib         ! 1 = run adiabatic, 0 = read adiabatic from file
 integer     :: Nad           !adiabatic switching parameters
 real (Long) :: tad,wtad  


end module input_parameters
