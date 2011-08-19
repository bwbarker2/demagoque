module input_parameters
 use prec_def
 implicit none

 logical :: initState_gaussianNuclear

 logical :: initState_cosine !add sine wave to initial state
 integer :: initState_cosine_number ! wave number n
 real (Long) :: initState_cosine_norm !normalization of sine wave, default 1
 real (Long) :: initState_cosine_shift !phase shift, default 0

 logical :: initState_plane  !add plane wave to initial state
 integer :: initState_plane_number
 real (Long) :: initState_plane_norm
 real (Long) :: initState_plane_shift

 logical :: initState_kdelta  !add kronecker delta function to initial state
 real (Long) :: initState_kdelta_norm
 real (Long) :: initState_kdelta_x0  !center of delta function

! logical :: initState_randomG  !add a random number (Gaussian) to initial state
! real (Long) :: initState_random_norm  !norm
! real (Long) :: initState_random_fwhm  !full width half max


end module input_parameters

