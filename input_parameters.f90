module input_parameters
 use prec_def
 implicit none

 logical :: initState_gaussianNuclear

 logical :: initState_cosine !add sine wave to initial state
 integer :: initState_cosine_number ! wave number n
 real (Long) :: initState_cosine_norm !normalization of sine wave, default 1
 real (Long) :: initState_cosine_shift !phase shift, default 0

end module input_parameters

