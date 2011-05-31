module lib_fftw
 implicit none

 include '/usr/include/fftw3.f'

 logical ft_re_1d_init

 contains

  subroutine ft_z2z_1d(arrayin, arrayout, num)
   !! ft_z2z - computes Fourier transform for generic complex 1D array.
   implicit none

   integer,    intent(in)  :: num
   complex*16 :: arrayin(0:num-1)
   complex*16 :: arrayout(0:num-1)

   integer*8 :: plan

!$OMP MASTER

   call dfftw_plan_dft_1d(plan,num,arrayin,arrayout, FFTW_FORWARD, FFTW_ESTIMATE)

!$OMP END MASTER

   call dfftw_execute(plan)

   call dfftw_destroy_plan(plan)

  end subroutine ft_z2z_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ift_z2z_1d(arrayin, arrayout, num)
   !! ft_z2z - computes inverse Fourier transform for generic complex 1D array.
   implicit none

   integer,    intent(in)  :: num
   complex*16 :: arrayin(0:num-1)
   complex*16 :: arrayout(0:num-1)

   integer*8 :: plan

   call dfftw_plan_dft_1d(plan,num,arrayin,arrayout, FFTW_BACKWARD, FFTW_ESTIMATE)

   call dfftw_execute(plan)

   call dfftw_destroy_plan(plan)
  end subroutine ift_z2z_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine ft_re_1d(arrayin, arrayout, num)
   implicit none
   
   integer,    intent(in)  :: num
   real*8 :: arrayin(0:num-1)
   real*8 :: arrayout(0:num-1)

   integer*8 :: plan

   call dfftw_plan_r2r_1d(plan,num,arrayin,arrayout,FFTW_REDFT00,FFTW_ESTIMATE)

   call dfftw_execute(plan)
   
  end subroutine ft_re_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ft_ro_1d(arrayin, arrayout, num)
   implicit none
   
   integer,    intent(in)  :: num
   real*8 :: arrayin(0:num-1)
   real*8 :: arrayout(0:num-1)

   integer*8 :: plan

   call dfftw_plan_r2r_1d(plan,num,arrayin,arrayout,FFTW_RODFT00,FFTW_ESTIMATE)

   call dfftw_execute(plan)
   
  end subroutine ft_ro_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  subroutine ft_z2d_1d(arrayin, arrayout, num)
!   !! ft_z2r_1d - computes complex to real Fourier transform, given a Hermitian symmetric input
!   implicit none
!   
!   integer, intent(in) :: num
!
!   end subroutine ft_z2d_1d

end module lib_fftw

