program demagoqueprocess
 use formatting
 use input_parameters
 use mesh
 use phys_cons
 use prec_def
 use time
 implicit none

 integer :: iit

 complex (Long), dimension(:,:), allocatable :: denmatInitial
 
 fout_dir='analysis/'

 call setModePrefix('')
 
 call init_prec_def
 call phys_cons_init

 call getStdIn

 call calcInitial

 denState=SPACE

 allocate(denmatInitial(Nxan:Nxax,Nxrn:Nxrx))

 Nt=Nevt

 call time_initialize

 call inDenUnf('results/ufo_'//time_getString(0)//'.dat',denmatInitial)

 firstOutput=.true.

 do iit=0,Nt/ntime
  call inDenUnf('results/ufo_'//time_getString(iit)//'.dat',denmat)
  denmat=denmat-denmatInitial
  call output
 
 enddo

 ! input zero timestep

 ! loop through all timesteps, subtract, store output in new output file



! do iit=

end program demagoqueprocess

