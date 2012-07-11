!> \mainpage demagoque - The Documentation
!!
!!   \section license License
!!    Copyright (C) 2011  Brent W. Barker
!!
!!    This program is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) any later version.
!!
!!    This program is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with this program (gpl-3.0.txt).  If not, see
!!    <http://www.gnu.org/licenses/>.
!!
!! \section author Author
!!            Brent W. Barker<br />
!!            barker at nscl dot msu dot edu<br />
!!            National Superconducting Cyclotron Laboratory<br />
!!            Michigan State University<br />
!!            1 Cyclotron, East Lansing, MI 48824-1321
!!
!! \section coords Coordinate system
!! <p>In the rotated coordinate system, from the original x,x' system, the
!! following coordinate system is used:</p>
!!
!! <code>
!! xa = (x+x')/2 <br />
!! xr = (x-x')   <br />
!! ka = (k+k')/2 <br />
!! kr = (k-k')   <br />
!! </code>
!!
!! <p>In the spatial coordinates, the density matrix is of the form <code>denmat(xa,xr)</code>.
!! In the spectral coordinates, it is of the form <code>denmat(kr,ka)</code>. Note that the
!! order of relative and absolute is switched. This is because the Fourier
!! transform in these coordinates associates the xr coordinate with ka, and xa
!! with kr.
!!
!! \section potentials Potentials
!! <p>These are selected with the potInitial and potFinal variables. The
!! initial state should probably be an eigenstate of <code>potInitial</code>. The <code>potFinal</code> potential is
!! that which is used for the time evolution. If potInitial and potFinal are
!! different, then the system is switched adiabatically from the initial to the
!! final potential. Here is the definition of the different potInitial/Final
!! integers:</p>
!!
!! <table>
!! <tr><td> code </td><td> definition </td></tr>
!! <tr><td> -1 </td><td> no potential at all, free space </td></tr>
!! <tr><td>  0 </td><td>  external harmonic oscillator centered at x=0</td></tr>
!! <tr><td> 1 </td><td>  nonlocal meanfield harmonic oscillator</td></tr>
!! <tr><td> 2 </td><td>  Skyrme-like contact potential (local density dependent) </td></tr>
!! <tr><td> 3 </td><td> same as pot=0, but with exact evolution from Chin, Krotsheck, Phys Rev E72, 036705 (2005) </td></tr>
!! <tr><td> 4 </td><td> external HO with effective 1D nonpolynomial meanfield from Mateo, Delgado, Malomed, Phys Rev A 83, 053610 (2011) </td></tr>
!! </table>
!!
!! \section options Options
!!
!! <p>Options are given at the end of the code. Eventually the potentials will be options as well. Each line consists of a space-separated list of parameters. The first is a string that defines the option to set, followed by a list of parameters for that option. Options can be safely commented out with a starting '!'</p>
!! <table>
!! <tr><td>option </td><td>  definition </td></tr>
!!
!! <tr><td>initialSeparation </td><td> initial separation between center-of-masses of fragments
!!                       in fm. Currently rounds displacement to nearest grid
!!                       point, rather than interpolating.<br />
!!
!!                       Arguments: real*8 initialSeparation</td></tr>
!!
!! <tr><td> splitOperatorMethod </td><td>  time evolve using SOM. Parameter is order of method.
!!                       Available orders are 3 and 5. Formulae from
!!                       A.D.Bandrauk, H. Shen, J. Chem. Phys. 99, 1185 (1993).<br />
!!
!!                       Arguments: integer splitOperatorMethod</td></tr>
!!
!! <tr><td> useFlipClone </td><td> create symmetric system by adding to the system its
!!                       conjugate, reflected about the xa axis. Option
!!                       'initialSeparation' must be set.<br />
!!
!!                       Arguments: None</td></tr>
!!
!! <tr><td> useImCutoff </td><td> imaginary off-diagonal cutoff. <br />
!!
!!                       Arguments: real*8 cutoff_w0, real*8 cutoff_x0, real*8
!!                       cutoff_d0 </td></tr>
!! </table>
PROGRAM dmtdhf
  use input_parameters
  use formatting
  USE mesh
  use time
  IMPLICIT NONE

!  integer :: fftw_init_out
  real :: timeElapsed

  call init_prec_def
  call phys_cons_init

!  call dfftw_init_threads(fftw_init_out)
!  if(fftw_init_out==0)write(*,*)'ERROR: dfftw_init_threads error:',fftw_init_out

!  call dfftw_plan_with_nthreads(2)

  CALL getStdIn

  CALL calcInitial

!  call fft_initial

  CALL initialState

!  call mesh_setReflectedLR(.true.) !debug, reflect it

  ! set output directory
  fout_dir='results/'

  !write(*,*)'finished initialState'
  if (useImEvol) then

   call setModePrefix('imev_')

   open(unit=41,file='results/imev_denmat_x_t.dat')
   open(unit=42,file='results/imev_denmat_k_t.dat')
   open(unit=43,file='results/imev_cons_rel.dat')
   open(unit=44,file='results/imev_cons_abs.dat')
   open(unit=45,file='results/imev_denmatan_x_t.dat') 
   OPEN(unit=61,file='results/imev_2dxim.dat')
   OPEN(unit=62,file='results/imev_2dxre.dat')
   OPEN(unit=66,file='results/imev_2dkim.dat')
   OPEN(unit=67,file='results/imev_2dkre.dat')
   open(unit=68,file='results/imev_2dwim.dat')
   open(unit=69,file='results/imev_2dwre.dat')
   OPEN(unit=70,file='results/imev_pk2.dat')
   open(unit=72,file='results/imev_2dx.dat')
   open(unit=73,file='results/imev_2dw.dat')
   open(unit=74,file='results/imev_2dk.dat')
   open(unit=75,file='results/imev_mean_abs_x.dat')
   open(unit=76,file='results/imev_pot.dat')

   Nt=Nimev

   call time_initialize

   call time_evolution

   call renormalizeDM

   useImEvol=.false.

   close(41)
   close(42)
   close(43)
   close(44)
   close(45)
   close(61)
   close(62)
   close(66)
   close(67)
   close(68)
   close(69)
   close(70)
   close(72)
   close(73)
   close(74)
   close(75)
   close(76)

  endif !useImEvol

  if (potInitial.NE.potFinal) useAdiabatic=.true.
!  useAdiabatic=.false.
  !write(*,*)'useAdiabatic=',useAdiabatic

  if (useAdiabatic.and.iadib==1) then

   call setModePrefix('ad_')

   ! open files to output adiabatic evolution information
   open(unit=41,file='results/ad_denmat_x_t.dat')
   open(unit=42,file='results/ad_denmat_k_t.dat')
   open(unit=43,file='results/ad_cons_rel.dat')
   open(unit=44,file='results/ad_cons_abs.dat')
   open(unit=45,file='results/ad_denmatan_x_t.dat') 
   open(unit=54,file='results/obdm.dat', form='unformatted')
   OPEN(unit=61,file='results/ad_2dxim.dat')
   OPEN(unit=62,file='results/ad_2dxre.dat')
   OPEN(unit=66,file='results/ad_2dkim.dat')
   OPEN(unit=67,file='results/ad_2dkre.dat')
   open(unit=68,file='results/ad_2dwim.dat')
   open(unit=69,file='results/ad_2dwre.dat')
   OPEN(unit=70,file='results/ad_pk2.dat')
   open(unit=72,file='results/ad_2dx.dat')
   open(unit=73,file='results/ad_2dw.dat')
   open(unit=74,file='results/ad_2dk.dat')
   open(unit=75,file='results/ad_mean_abs_x.dat')
   open(unit=76,file='results/ad_pot.dat')

   !set number of timesteps for adiabatic switching
   Nt=Nad

!  call mesh_setReflectedLR(.true.) !debug, reflect it
   call time_evolution
!  call mesh_setReflectedLR(.false.) !debug, reflect it back

   call outDenUnf

   ! close files
   close(41)
   close(42)
   close(43)
   close(44)
   close(45)
   close(54)
   close(61)
   close(62)
   close(66)
   close(67)
   close(68)
   close(69)
   close(70)
   close(72)
   close(73)
   close(74)
   close(75)
   close(76)

   useAdiabatic=.false.

  elseif(useAdiabatic) then

   open(unit=54,file='results/obdm.dat', form='unformatted')

   call inDenUnf

   close(54)

   useAdiabatic=.false.
 
  endif !useAdiabatic

   useAdiabatic=.false.

  call setModePrefix('')
   
  ! open files to output time evolution information
  open(unit=41,file='results/denmat_x_t.dat')
  open(unit=42,file='results/denmat_k_t.dat')
  open(unit=43,file='results/cons_rel.dat')
  open(unit=44,file='results/cons_abs.dat')
  OPEN(unit=61,file='results/2dxim.dat')
  OPEN(unit=62,file='results/2dxre.dat')
  OPEN(unit=66,file='results/2dkim.dat')
  OPEN(unit=67,file='results/2dkre.dat')
  open(unit=68,file='results/2dwim.dat')
  open(unit=69,file='results/2dwre.dat')
  OPEN(unit=70,file='results/pk2.dat')
  OPEN(unit=71,file='results/sumsquarediff.dat')
  open(unit=45,file='results/denmatan_x_t.dat') 
  open(unit=72,file='results/2dx.dat')
  open(unit=73,file='results/2dw.dat')
  open(unit=74,file='results/2dk.dat')
  open(unit=75,file='results/mean_abs_x.dat')
  open(unit=76,file='results/pot.dat')
!  open(unit=101,file='results/pot_diag.dat')

  Nt=Nevt

!  write(*,*)'new Nt=',Nt

!  CALL OUTPUT
!  do it=1,10
!   write(*,*)it
!   call output
!   call outX
!   call outW
!   call outK
!   call outW
!  enddo

!call mesh_setReflectedLR(.true.)
!call mesh_setReflectedLR(.false.)

  if(ea>0e0_Long.or.ea<0e0_Long)call boost
 
!call mesh_setReflectedLR(.false.)

  if(abs(initialSeparation)>epzero) then
   if(useFrameXXP) then
    call mesh_xxp_displace(NINT(initialSeparation*0.25_Long*Nxa/xLa))
   elseif(useMeshXAR2) then
    call mesh_xar2_displace(NINT(initialSeparation*0.25_Long*Nxa/xLa))
   endif
  endif

!  call mesh_setReflectedLR(.false.)

  if(useFlipClone) call throwException('main: flipclone not implemented yet.', BEXCEPTION_FATAL)!call flipclone
!  write(*,*)'flipclone finished'

!  call mesh_setReflectedLR(.false.)
  maxxim=0e0_Long

!  call outX
!  call transform_x_to_w_shift
!  call outW
!  call transform_w_to_x_shift
!  call outX


  CALL time_evolution

  write(*,*)'maxxim:',maxxim

!  do it=1,100
!   write(*,*)it
!   call output
!  enddo

  ! close files
  close(41)
  close(42)
  close(43)
  close(44)
  close(45)
  close(61)
  close(62)
  close(66)
  close(67)
  close(68)
  close(69)
  close(70)
  close(71)
  close(72)
  close(73)
  close(74)
  close(75)
  close(76)
!  close(101)

  call cpu_time(timeElapsed)
  write(*,*)'total time:',timeElapsed,'seconds'

END PROGRAM dmtdhf
