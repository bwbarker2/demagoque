! Copyright (C) 2011  Brent W. Barker
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program (gpl-3.0.txt).  If not, see
!    <http://www.gnu.org/licenses/>.
!
!    Author: Brent W. Barker
!            barker at nscl dot msu dot edu
!            National Superconducting Cyclotron Laboratory
!            Michigan State University
!            1 Cyclotron, East Lansing, MI 48824-1321

module lib_fftw
 use prec_def
 implicit none

 include '/usr/include/fftw3.f'

 logical ft_re_1d_init

 contains

  subroutine ft_z2z_1d(arrayin, arrayout, num)
   !! ft_z2z - computes Fourier transform for generic complex 1D array.
   implicit none

   integer,    intent(in)  :: num
   complex (Long) :: arrayin(0:num-1)
   complex (Long) :: arrayout(0:num-1)

   integer(kind=8) :: plan

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
   complex (Long) :: arrayin(0:num-1)
   complex (Long) :: arrayout(0:num-1)

   integer(kind=8) :: plan

   call dfftw_plan_dft_1d(plan,num,arrayin,arrayout, FFTW_BACKWARD, FFTW_ESTIMATE)

   call dfftw_execute(plan)

   call dfftw_destroy_plan(plan)
  end subroutine ift_z2z_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine ft_re_1d(arrayin, arrayout, num)
   implicit none
   
   integer,    intent(in)  :: num
   real (Long) :: arrayin(0:num-1)
   real (Long) :: arrayout(0:num-1)

   integer(kind=8) :: plan

   call dfftw_plan_r2r_1d(plan,num,arrayin,arrayout,FFTW_REDFT00,FFTW_ESTIMATE)

   call dfftw_execute(plan)
   
  end subroutine ft_re_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ft_ro_1d(arrayin, arrayout, num)
   implicit none
   
   integer,    intent(in)  :: num
   real (Long) :: arrayin(0:num-1)
   real (Long) :: arrayout(0:num-1)

   integer(kind=8) :: plan

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

