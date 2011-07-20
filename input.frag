! ........... READ AND WRITE PARAMETERS
10 ! ntime
! ........... PHYSICAL PARAMETERS
0.5 ! delt [fm/c]
200 ! Nevt
1   ! Nmax
25d0 ! EA [MeV]
25d0  ! xLa [fm]
9d0 ! xLr [fm]
100 ! Nxa
36 ! Nxr
! ............. POTENTIAL PARAMETERS
0  ! potInitial
2 ! potFinal
30 70 1200  ! tad,wtad,Nad
1  ! iadib
! ............. IMAGINARY EVOLUTION PARAMETERS
.false.     ! useImEvol - use imaginary evolution?
100   ! nimev - number of timesteps for this evolution
! ............. START OPTIONAL PARAMETERS
initialSeparation 15d0
splitOperatorMethod 3
useImCutoff 1000d0 6d0 2d0  !useImCutoff, cutoff w0,x0,d0
useFlipClone
END_OF_OPTIONS !sentinel for end of optional parameters

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! INPUT FILE FOR THE EVOLUTION PROGRAM
! ntime ->  Controls how often the data is written
! delt -> Time interval                                                    c
! Nt -> Number of time intervals                                           c
! Nmax -> Maximum oscillator shell
! EA  -> initial kinetic energy per particle
! xLa -> Total length of the box in average coordinate xa=(x+x')/2 is 2*xLa
! xLr -> Total length of the box in relative coordinate xr=x-x' is 2*xLr
! Nxa -> number of grid points in average coordinate
! Nxr -> number of grid points in relative coordinate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
