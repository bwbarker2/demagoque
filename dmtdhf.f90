PROGRAM dmtdhf
  ! dmtdhf - time evolution in 1D of a 1-body density matrix.
!
! In the rotated coordinate system, from the original x,x' system, the
! following coordinate system is used:
! xa = (x+x')/2
! xr = (x-x')
! ka = (k+k')/2
! kr = (k-k')
!
! in the spatial coordinates, the density matrix is of the form denmat(xa,xr).
! In the spectral coordinates, it is of the form denmat(kr,ka). Note that the
! order of relative and absolute is switched. This is because the Fourier
! transform in these coordinates associates the xr coordinate with ka, and xa
! with kr.
!
! Potentials: These are selected with the potInitial and potFinal variables. The
! initial state must be an eigenstate of potInitial. The potFinal potential is
! that which is used for the time evolution. If potInitial and potFinal are
! different, then the system is switched adiabatically from the initial to the
! final potential. Here is the definition of the different potInitial/Final
! integers:
!
! pot  definition
! -1   no potential at all, free space
! 0    external harmonic oscillator centered at x=0
! 1    nonlocal meanfield harmonic oscillator
! 2    Skyrme-like contact potential (local density dependent)
! 3    same as pot=0, but with exact evolution from Chin, Krotsheck, Phys Rev E72, 036705 (2005)
!
  USE mesh
  USE time
  IMPLICIT NONE

!  integer :: fftw_init_out
  real :: timeElapsed(2)

!  fftw_init_out=dfftw_init_threads()
!  call dfftw_init_threads(fftw_init_out)
!  if(fftw_init_out==0)write(*,*)'ERROR: dfftw_init_threads error:',fftw_init_out

!  call dfftw_plan_with_nthreads(2)

  CALL getStdIn

  CALL calcInitial

!  call fft_initial

  CALL initialState

  !write(*,*)'finished initialState'
  if (useImEvol) then

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

  endif !useImEvol

  if (potInitial.NE.potFinal) useAdiabatic=.true.
!  useAdiabatic=.false.
  !write(*,*)'useAdiabatic=',useAdiabatic

  if (useAdiabatic.and.iadib==1) then
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

   !set number of timesteps for adiabatic switching
   Nt=Nad

   call time_evolution

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

   useAdiabatic=.false.

  elseif(useAdiabatic) then

   open(unit=54,file='results/obdm.dat', form='unformatted')

   call inDenUnf

   close(54)

   useAdiabatic=.false.
 
  endif !useAdiabatic

   useAdiabatic=.false.
 
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

  call boost

  call displace(10)

  call flipclone

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
!  close(101)

  write(*,*)'total time:',etime(timeElapsed)
  write(*,*)'user,walltime:',timeElapsed

END PROGRAM dmtdhf


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! reads data from standard input
SUBROUTINE getStdIn
  USE mesh
  USE osc_pars
  USE out
  USE phys_cons
  USE prec_def
  USE time
  IMPLICIT NONE

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

  ! oscillator maximum shell
  READ(*,*) Nmax
  WRITE(*,*) 'Maximum oscillator shell, Nmax=',Nmax

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

END SUBROUTINE getStdIn
