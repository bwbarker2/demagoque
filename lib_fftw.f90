module lib_fftw
 implicit none

 include '/usr/include/fftw3.f'

 contains

  subroutine ft_z2z_1d(arrayin, arrayout, num)
   !! ft_z2z - computes Fourier transform for generic complex 1D array.
   implicit none

   integer,    intent(in)  :: num
   complex*16, intent(in)  :: arrayin(0:num-1)
   complex*16, intent(out) :: arrayout(0:num-1)

   integer*8 :: plan

   call dfftw_plan_dft_1d(plan,num,arrayin,arrayout, FFTW_FORWARD, FFTW_ESTIMATE)

   call dfftw_execute(plan)

  end subroutine ft_z2z_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ift_z2z_1d(arrayin, arrayout, num)
   !! ft_z2z - computes inverse Fourier transform for generic complex 1D array.
   implicit none

   integer,    intent(in)  :: num
   complex*16, intent(in)  :: arrayin(0:num-1)
   complex*16, intent(out) :: arrayout(0:num-1)

   integer*8 :: plan

   call dfftw_plan_dft_1d(plan,num,arrayin,arrayout, FFTW_BACKWARD, FFTW_ESTIMATE)

   call dfftw_execute(plan)

  end subroutine ift_z2z_1d

end module lib_fftw
