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

   Nt=Nimev

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

  if(ea>0d0.or.ea<0d0)call boost
 
!call mesh_setReflectedLR(.false.)

  if(abs(initialSeparation)>(delxa*0.5d0)) then
   if(initialSeparation>0d0) then
    call displaceLeft(NINT(initialSeparation*0.5d0/delxa))
   else
    call displaceRight(NINT(abs(initialSeparation)*0.5d0/delxa))
   endif
  endif

!  call mesh_setReflectedLR(.false.)

  if(useFlipClone) call flipclone
!  write(*,*)'flipclone finished'

!  call mesh_setReflectedLR(.false.)
  maxxim=0.d0


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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! reads data from standard input
SUBROUTINE getStdIn
 use bexception
 use input_parameters
 use bstring
  USE mesh
  USE phys_cons
  USE prec_def
  IMPLICIT NONE

 character(len=80) :: inline !input line for optional line processing
 integer :: ibeg,iend

  ! ... READ DATA FOR READING AND WRITING
  READ(*,*)

  ! ntime = write data every ntime steps
  READ(*,*) ntime
  WRITE(*,*) 'Write data every',ntime,' iterations'

  ! ... READ PHYSICAL PARAMETERS
  READ(*,*)

  ! delta of time
  READ(*,*) delt
  WRITE(*,*) 'Time increment dt=',delt,' fm/c'

  ! number of timesteps in time evolution
  READ(*,*) Nevt
  WRITE(*,*) 'Number of time increments Nevt=',Nevt


  ! energy per particle in MeV
  READ(*,*) ea
  WRITE(*,*) 'Energy per particle [MeV], ea=',ea

  ! length of box in x_average=(x+x')/2
  READ(*,*) xLa
  WRITE(*,*) 'Length of the box/2 [fm], xLa=',xLa

  ! length of box in x_relative=(x-x')
  READ(*,*) xLr
  WRITE(*,*) 'Width of the box/2 [fm], xLr=',xLr

  ! number of mesh points in x_average
  READ(*,*) Nxa
  WRITE(*,*) 'Number of mesh points in x_average, Nxa=',Nxa

  ! number of mesh points in x_relative
  READ(*,*) Nxr
  WRITE(*,*) 'Number of mesh points in x_relative, Nxr=',Nxr

  ! ... READ POTENTIAL PARAMETERS
  READ(*,*)

  ! type of initial potential
  read(*,*) potInitial
  write(*,*) 'Type # of initial state potential, potInitial=',potInitial

  ! type of final potential
  read(*,*) potFinal
  write(*,*) 'Type # of final state potential, potFinal=',potFinal

  ! adiabatic parameters
  read(*,*) tad,wtad,Nad
  write(*,*) 'adiabatic parameters: tad,wtad,Nad=',tad,wtad,Nad

  ! read from file?
  read(*,*) iadib
  write(*,*) 'read initial state from file, iadib=',iadib

  read(*,*)

  read(*,*) useImEvol
  write(*,*) 'use imaginary evolution? useImEvol=',useImEvol

  read(*,*) Nimev
  write(*,*) 'timesteps for imaginary evolution, Nimev=',Nimev

 !set default options
 initialSeparation=0d0  !don't use initial separation
 initState_gaussian=.false.
 initState_gaussianNuclear=.false.
 initState_cosine=.false.
 initState_kdelta=.false.
 initState_plane=.false.
 Nmax=0
 norm_thy=0d0
 unitSystem_bec=.false.
 unitSystem_nuclear=.true.
 useImCutoff=.false.
 useFlipClone=.false.
 splitOperatorMethod=0  !don't use Split Operator Method

 !read optional lines
 do while(.true.)

  read(*,'(A)') inline

  if(isComment(inline))cycle

  call findFirstWord(inline,' ',ibeg,iend)

  !if it's the sentinel, then exit loop
  if(inline(ibeg:iend)=="END_OF_OPTIONS") then
!   write(*,*)'input sentinel reached, exiting input loop'
   exit
  endif !inline

  call procOptionLine(inline)

 enddo !while true

 !verify option dependencies

 if(useFlipClone.and.abs(initialSeparation)<xLa/Nxa)then
  call throwException( &
   'getStdIn: useFlipClone is set, but not initialSeparation!' &
    //'Setting initialSeparation to half xLa',BEXCEPTION_WARNING)
  initialSeparation=xLa/2d0
 endif

END SUBROUTINE getStdIn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Processes line of options
subroutine procOptionLine(inline)
 use bexception
 use bstring
 use input_parameters
 use mesh
 use phys_cons
 use prec_def
 implicit none

 character(*), intent(in) :: inline !< line to be processed

 integer :: ibeg,iend
 integer, dimension(3) :: iin
 real(Long), dimension(3) :: rin

 call findFirstWord(inline,' ',ibeg,iend)

 select case(inline(ibeg:iend))

  case("initialSeparation")
   read(inline(iend+1:len(inline)),*)initialSeparation
   write(*,*)'Initial Separation, initialSeparation=' &
             ,initialSeparation,'fm'

  case("initState_gaussianNuclear")
    read(inline(iend+1:len(inline)),*)Nmax
    initState_gaussianNuclear=.true.
    norm_thy=norm_thy+Nmax+1
    write(*,*)'Initial state added: initState_gaussianNuclear'
    WRITE(*,*)'  Maximum oscillator shell, Nmax=',Nmax

   case("initState_cosine")
    read(inline(iend+1:len(inline)),*)initState_cosine_number, &
                     initState_cosine_norm, initState_cosine_shift
    initState_cosine=.true.
    norm_thy=norm_thy+initState_cosine_norm
    write(*,*)'Initial state added: initState_cosine'
    write(*,*)'                     initState_cosine_number=',initState_cosine_number
    write(*,*)'                     initState_cosine_norm  =',initState_cosine_norm
    write(*,*)'                     initState_cosine_shift =',initState_cosine_shift

   case("initState_kdelta")
    read(inline(iend+1:len(inline)),*)initState_kdelta_norm, initState_kdelta_x0
    initState_kdelta=.true.
    norm_thy=norm_thy+initState_kdelta_norm
    write(*,*)'Initial state added: initState_kdelta'
    write(*,*)'                     initState_kdelta_norm=',initState_kdelta_norm
    write(*,*)'                     initState_kdelta_x0=',initState_kdelta_x0

   case("initState_plane")
    read(inline(iend+1:len(inline)),*)initState_plane_number, &
                   initState_plane_norm, initState_plane_shift 
    initState_plane=.true.
    norm_thy=norm_thy+initState_plane_norm
    write(*,*)'Initial state added: initState_plane'
    write(*,*)'                     initState_plane_number=',initState_plane_number
    write(*,*)'                     initState_cosine_norm  =',initState_plane_norm
    write(*,*)'                     initState_cosine_shift =',initState_plane_shift

   case("initState_sqWell")
    if(.not.phys_cons_isInitialized) then
     call throwException('procOptionLine: option "unitSystem" must be set before initState_sqWell',BEXCEPTION_FATAL)
    else
     useInitState_sqWell=.true.
     read(inline(iend+1:len(inline)),*)rin(1),iin(1),rin(2),rin(3)
     initState_sqWell=make_SquareWellState(rin(1),iin(1),rin(2),rin(3))
     norm_thy=norm_thy+1
    endif

   case("pot4")
    read(inline(iend+1:len(inline)),*)ho_mateo_wz, ho_mateo_wt, &
                                      ho_mateo_scat, ho_mateo_Npart
    initState_gaussian=.true.
    norm_thy=norm_thy+1
    write(*,*)'Potential 4 options: ho_mateo_wz=',ho_mateo_wz
    write(*,*)'                     ho_mateo_wt=',ho_mateo_wt
    write(*,*)'                     ho_mateo_scat=',ho_mateo_scat
    write(*,*)'                     ho_mateo_Npart=',ho_mateo_Npart
  
   case("splitOperatorMethod")
    read(inline(iend+1:len(inline)),*)splitOperatorMethod
    write(*,*)'Using Split Operator Method, order=',splitOperatorMethod

   case("useImCutoff")
    useImCutoff=.true.
    read(inline(iend+1:len(inline)),*)cutoff_w0,cutoff_x0,cutoff_d0
    write(*,*)'Using imaginary off-diagonal cutoff, w0,x0,d0=' &
              ,cutoff_w0,cutoff_x0,cutoff_d0

   case("unitSystem")
    ! reset unit system declaration so we can specify a new one.
    unitSystem_bec=.false.
    unitSystem_nuclear=.false.

    call findFirstWord(inline(iend+1:len(inline)),' ',ibeg,iend)
    ibeg=ibeg+10  !10 is length of string 'unitSystem' (previous iend)
    iend=iend+10
    select case(inline(ibeg:iend))
     case("bec")
      unitSystem_bec=.true.
      write(*,*)'Defining unit system: bec'
      call phys_cons_initializeBEC

     case("nuclear")
      unitSystem_nuclear=.true.
      write(*,*)'Defining unit system: nuclear'
      call phys_cons_initializeNuclear

     case default
      call throwException('getStdIn: invalid unit system chosen: '//inline(ibeg:iend),BEXCEPTION_FATAL)
     
    end select ! unitSystem

   case("useFlipClone")
    useFlipClone=.true.
    write(*,*)'Making symmetric collision with flipClone'

   case default
    call throwException('getStdIn: Option does not exist: ' // inline(ibeg:iend),BEXCEPTION_FATAL)

  end select ! inline

end subroutine procOptionLine
