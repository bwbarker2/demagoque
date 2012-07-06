!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! reads data from standard input
SUBROUTINE getStdIn
 use bexception
 use input_parameters
 use bstring
 use class_PotHO
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
 initialSeparation=0.d0  !don't use initial separation
 initPotentialList => new_PotentialList()
 initState_gaussian=.false.
 initState_gaussianNuclear=.false.
 initState_cosine=.false.
 initState_kdelta=.false.
 initState_plane=.false.
 initSuperWavefunction => new_SuperWavefunction()
 Nmax=0
 norm_thy=0d0
 unitSystem_bec=.false.
 unitSystem_nuclear=.true.
 useCutoffK=.false.
 useImCutoff=.false.
 useFlipClone=.false.
 useFrameXXP=.false.
 useMeshShifted=.false.
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

write(*,*)'procOptionLine: initSuperWavefunction=',initSuperWavefunction%getWavefn(20.2_Long)
 !verify option dependencies

 if(useCutoffK.and.(2*cutoffK_ncells>=Nxr)) then
  call throwException('getStdIn: cutoffK_ncells will cut more than' &
                      // ' entire k-space!',BEXCEPTION_FATAL)
 endif

 if(useFlipClone.and.abs(initialSeparation)<xLa/Nxa)then
  call throwException( &
   'getStdIn: useFlipClone is set, but not initialSeparation!' &
    ,BEXCEPTION_WARNING)
  initialSeparation=xLa/2e0_Long
 endif

 if(useFrameXXP.and.((Nxa/=Nxr).or.abs(xLa-xLr)>epzero)) then
  call throwException( &
   'getStdIn: useFrameXXP=.true, but Nxa/=Nxr or xLa/=xLr. Both of these ' &
   //'need to be equal to each other.', BEXCEPTION_FATAL)
 endif

 if(useFrameXXP.and.useMeshShifted) then
  call throwException( &
   'getStdIn: useFrameXXP and useMeshShifted are both .true. Only one of ' &
   //'them can be true. Do not shift the mesh when in (x,x-prime) ' &
   //'representation.' &
   ,BEXCEPTION_FATAL)
 endif

 if(useMeshShifted.and.(isOdd(Nxa).or.isOdd(Nxr))) then
  call throwException( &
   'getStdIn: useMeshShifted is set, but it requires even Nxa and Nxr!' &
   ,BEXCEPTION_FATAL)
 endif

! if(.not.useMeshShifted.and.(isEven(Nxa).or.isEven(Nxr)).and..not.useFrameXXP) then
!  call throwException( &
!   'getStdIn: useMeshShifted=.false., but this requires odd Nxa and Nxr!' &
!   ,BEXCEPTION_FATAL)
! endif

 if(.not.(useMeshXAR2.neqv.useFrameXXP)) then
  call throwException( &
   'getStdIn: only useMeshXAR2 and useFrameXXP implemented now. Must use 1' &
   //' of these and not both.',BEXCEPTION_FATAL)
 endif

 if(useMeshXAR2.and.(useMeshShifted.or.useFrameXXP)) then
  call throwException( &
   'getStdIn: useMeshXAR2 cannot be used with any other useMesh or useFrame' &
   // ' option!',BEXCEPTION_FATAL)
 endif

END SUBROUTINE getStdIn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Processes line of options
subroutine procOptionLine(inline)
 use bexception
 use bstring
 use input_parameters
 use class_PotHO
 use class_PotSquareWell
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

  case("initPotHO")
   read(inline(iend+1:len(inline)),*)rin(1)
   call initPotentialList%add(new_PotHO(m0,rin(1)))
   write(*,*)'Added initial potential HO, frequency=',rin(1)
write(*,*)

  case("initPotSquareWell")
    read(inline(iend+1:len(inline)),*)rin(1),rin(2)
    call initPotentialList%add(new_PotSquareWell(rin(1),rin(2)))
    write(*,*)'Added initial potential SquareWell. half-width,v0=',rin(1),rin(2)

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

   case("initState_KronigPenney")
    read(inline(iend+1:len(inline)),*)iin(1),rin(1),rin(2)
    useInitState_KronigPenney=.true.
!    call initSuperWavefunction%add(new_WfKronigPenney( &
!                             m0,iin(1),rin(1),rin(2),2._Long*xLa))
    initState_KronigPenney => new_WfKronigPenney( &
                             m0,iin(1),rin(1),rin(2),2._Long*xLa)
!write(ERROR_UNIT,*)'between'
    call initSuperWavefunction%add(initState_KronigPenney)
!write(*,*)'procOptionLine: initState_KronigPenney=',initState_KronigPenney%getWavefn(20.2_Long)
!write(*,*)'procOptionLine: initSuperWavefunction=',initSuperWavefunction%getWavefn(20.2_Long)
    norm_thy=norm_thy+1._Long
    write(*,*)'Initial state added: initState_KronigPenney'

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
     initState_sqWell=make_WfSquareWell(rin(1),iin(1),rin(2),rin(3))
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

   case("useMeshShifted")
    useMeshShifted=.true.
    write(*,*)'Using shifted mesh: xj=dx*(j+1/2), k=dk*(k+1/2)'

   case("useMeshXAR2")
    useMeshXAR2=.true.
    write(*,*)'Using XAR2 mesh (chessboard in SPACE)'

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

   case("useCutoffK")
    useCutoffK=.true.
    read(inline(iend+1:len(inline)),*)cutoffK_ncells
    write(*,*)'Using cutoff on edges of k,k-prime space, cutoffK_ncells=' &
              ,cutoffK_ncells

   case("useFlipClone")
    useFlipClone=.true.
    write(*,*)'Making symmetric collision with flipClone'

   case("useFrameXXP")
    useFrameXXP=.true.
    write(*,*)'Using unrotated frame, (x,x-prime) representation'

   case default
    call throwException('getStdIn: Option does not exist: ' // inline(ibeg:iend),BEXCEPTION_FATAL)

  end select ! inline

end subroutine procOptionLine
