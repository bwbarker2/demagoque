module input_parameters
 use class_SuperWavefunction
 use class_WfKronigPenney
 use class_WfSquareWell
 use prec_def
 implicit none

 logical     :: initState_gaussian
 real (Long) :: ho_mateo_wz,   &  !< angular freq in z direction
                ho_mateo_wt,   &  !< angular freq in transverse direction
                ho_mateo_scat, &  !< s-wave scattering length
                ho_mateo_Npart    !< total number of particles

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

! logical     :: initState_BEC_q1D  !< use Bose-Einstein Condensate for initial state
! real (Long) :: initState_BEC_q1D_numpart  !< number of atoms in condensate
! real (Long) :: initState_BEC_q1D_scatLength  !< s-wave scattering length

 logical     :: initState_gaussianNuclear
 REAL (long) :: w      !< angular frequency
 REAL (long) :: whm    !< angular wavenumber squared
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

 logical     :: useInitState_KronigPenney !< use Kronig-Penney model eigenstates
 class(WfKronigPenney), pointer :: initState_KronigPenney !< Kronig-Penney eigenstate

 logical     :: useInitState_sqWell      !< use square well potential
 ! \todo try changing this to class instead of type in future versions of gfortran - it works as a class with ifort 12.1.0
 type(WfSquareWell) :: initState_SqWell !< square well state object

 class(SuperWavefunction),pointer :: initSuperWavefunction !< initial state wavefunctions

! logical :: initState_randomG  !add a random number (Gaussian) to initial state
! real (Long) :: initState_random_norm  !norm
! real (Long) :: initState_random_fwhm  !full width half max

 integer     :: splitOperatorMethod !0 if not used, otherwise order of method

 logical     :: unitSystem_bec      !< Use unit system for Bose-Einstein Condensate
 logical     :: unitSystem_nuclear  !< Use unit system for nuclear collision [default]

 logical     :: useCutoffK  !< cutoff in both k and k' edges
 integer     :: cutoffK_ncells  !< number of cells to cut off

 logical     :: useImEvol   ! use imaginary evolution?
 integer     :: Nimev  ! number of timesteps for imaginary evolution
 
 logical     :: useFlipClone  ! symmetric collision method
 logical     :: useMeshShifted !shift mesh by 1/2 of index, (j+1/2)(k+1/2)

 logical     :: useMeshXAR2  !< use rotated mesh with 2N x 2N system

 logical     :: useFrameXXP !< use unrotated (x,x') coordinate frame

 logical     :: useAdiabatic  ! if true, use adiabatic switching
 integer     :: iadib         ! 1 = run adiabatic, 0 = read adiabatic from file
 integer     :: Nad           !adiabatic switching parameters
 real (Long) :: tad,wtad  

end module input_parameters

