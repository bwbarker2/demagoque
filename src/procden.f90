program procden
 use mesh
 implicit none
 
 integer :: fu_2dxre=1
 integer :: fu_denlrsym=2
 integer :: fu_denudsym=3
 integer :: timesteps=40

 complex*16, allocatable,dimension(:,:) :: denlrsym, denudsym
 real*8 :: lrsymmetry, udsymmetry
 real*8 :: lrsymdiff, udsymdiff

 integer :: it,ixa,ixr

 Nxa=100
 Nxr=100

 call initializeMesh

 allocate(denlrsym(1:Nxa2-1,-Nxr2:Nxr2-1))
 allocate(denudsym(-Nxa2:Nxa2-1,1:Nxr2-1))
 denlrsym=cmplx(0d0,0d0,8)
 denudsym=cmplx(1d0,0d0,8)

 open(unit=fu_2dxre,file='results/2dwim.dat', status='old')
 open(unit=fu_denlrsym,file='results/2dwimlrsym.dat')
 open(unit=fu_denudsym,file='results/2dwimudsym.dat')

 do it=0,timesteps
  write(fu_denlrsym,*)'#timestep=',it
  write(fu_denudsym,*)'#timestep=',it

  read(fu_2dxre,*) !comment line

!!!!! reading brent's file !!!!!!!!!!!!!

  !don't read first quarter of data (redundant info)
  do ixr=1,Nxr2  
   read(fu_2dxre,*)
  enddo
 
  !fill denbrent
  do ixr=-Nxr2,Nxr2-1
   read(fu_2dxre,*)(den_re(ixa,ixr),ixa=-Nxa2,Nxa2-1)
!   write(*,*)den_re(:,ixr)
  enddo

!!!!!!!!!! end reading brent's file !!!!!!


!!!!!!!!!! begin reading arnau's file !!!!!!!!!
!
! !don't read first line, it's a duplicate of the last line
! read(fu_2dxre,*)
! !fill denarnau
! do ixa=Nxa2-1,-Nxa2,-1
!  read(fu_2dxre,*)(den_re(ixa,ixr),ixr=Nxr2,-Nxr2,-1)
! enddo
!
!!!!!!!!! end reading arnau's file !!!!!!!!!!!!!!

  do ixa=1,Nxa2-1
   denlrsym(ixa,:)=den_re(ixa,-Nxr2:Nxr2-1)-den_re(-ixa,-Nxr2:Nxr2-1)
  enddo

  lrsymmetry=SUM(abs(denlrsym))/(Nxa2-1)/Nxr
  lrsymdiff=SUM(denlrsym)/(Nxa2-1)/Nxr

  write(*,*)'lrsymmetry,lrsymdiff:',lrsymmetry,lrsymdiff

  do ixr=-Nxr2,Nxr2-1
   write(fu_denlrsym,*)(DBLE(denlrsym(ixa,ixr)),ixa=1,Nxa2-1)
  enddo

  ! begin up-down symmetry calculation

  do ixr=1,Nxr2-1
   denudsym(:,ixr)=den_re(:,ixr)-den_re(:,-ixr)
  enddo

  udsymmetry=SUM(abs(denudsym))/(Nxr2-1)/Nxa
  udsymdiff=SUM(denudsym)/(Nxr2-1)/Nxa

  write(*,*)'udsymmetry, udsymdiff:',udsymmetry,udsymdiff

  do ixr=1,Nxr2-1
   write(fu_denudsym,*)(DBLE(denudsym(ixa,ixr)),ixa=-Nxa2,Nxa2-1)
  enddo

!!!!!!! begin brent's file reading !!!!! 
  !get to next timestep
  do ixr=1,Nxr2+2  
   read(fu_2dxre,*)
  enddo
!!!!!!!!!! end brent's file reading !!!!!

!!!!!! begin arnau's file reading !!!!
! read(fu_2dxre,*)
!!!!! end arnau's file reading

 enddo !it

 close(fu_2dxre)

end program procden

