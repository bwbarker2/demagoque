program procdenextra
 use mesh
 implicit none
 
 integer :: fu_2dxre=1
 integer :: fu_denlrsym=2
 integer :: fu_denudsym=3
 integer :: fu_diagudsym=50
 integer :: fu_diaglrsym=51
 integer :: timesteps=23

 complex*16, allocatable,dimension(:,:) :: denlrsym, denudsym
 real*8 :: lrsymmetry, udsymmetry
 real*8 :: lrsymdiff, udsymdiff
 character, dimension(3) :: state=(/'x','w','k'/)        ! x, w, or k
 character(len=2), dimension(2) :: reim=(/'re','im'/)  !real or imaginary part

 integer :: it,ixa,ixr,istate,ireim

 Nxa=100
 Nxr=36

 call initializeMesh

 allocate(denlrsym(1:Nxa2-1,-Nxr:Nxr-1))
 allocate(denudsym(-Nxa2:Nxa2-1,1:Nxr-1))
 
 do istate=1,3
 do ireim=1,2

 denlrsym=cmplx(0d0,0d0,8)
 denudsym=cmplx(0d0,0d0,8)

 open(unit=fu_2dxre,file='results/2d'//state(istate)//reim(ireim)//'.dat', status='old')
 open(unit=fu_denlrsym,file='results/2d'//state(istate)//reim(ireim)//'lrsym.dat')
 open(unit=fu_denudsym,file='results/2d'//state(istate)//reim(ireim)//'udsym.dat')
 open(unit=fu_diagudsym,file='results/diag'//state(istate)//reim(ireim)//'udsym.dat')
 open(unit=fu_diaglrsym,file='results/diag'//state(istate)//reim(ireim)//'lrsym.dat')

 write(*,*)'processing 2d'//state(istate)//reim(ireim)

 do it=0,timesteps

  write(*,*)'index:',it

  write(fu_denlrsym,*)'#timestep=',it
  write(fu_denudsym,*)'#timestep=',it

  read(fu_2dxre,*) !comment line

!!!!! reading brent's file !!!!!!!!!!!!!

!  !don't read first quarter of data (redundant info)
!  do ixr=1,Nxr2  
!   read(fu_2dxre,*)
!  enddo
 
  !fill denbrent
  do ixr=-Nxr,Nxr-1
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

 if(reim(ireim).eq.'im')then
  do ixa=0,Nxa2-1
   denlrsym(ixa,:)=den_re(ixa,-Nxr:Nxr-1)+den_re(-ixa,-Nxr:Nxr-1)
  enddo
 else
  do ixa=1,Nxa2-1
   denlrsym(ixa,:)=den_re(ixa,-Nxr:Nxr-1)-den_re(-ixa,-Nxr:Nxr-1)
  enddo
 endif
  

  lrsymmetry=SUM(abs(denlrsym))/(Nxa2-1)/(2*Nxr)
  lrsymdiff=SUM(denlrsym)/(Nxa2-1)/(2*Nxr)

  write(*,*)'lrsymmetry,lrsymdiff:',lrsymmetry,lrsymdiff
  write(fu_diaglrsym,*)it,lrsymmetry,lrsymdiff

  do ixr=-Nxr,Nxr-1
   write(fu_denlrsym,*)(DBLE(denlrsym(ixa,ixr)),ixa=1,Nxa2-1)
  enddo

  ! begin up-down symmetry calculation
  if(reim(ireim).eq.'im')then
   do ixr=1,Nxr-1
    denudsym(:,ixr)=den_re(:,ixr)+den_re(:,-ixr)
   enddo
  else
   do ixr=1,Nxr-1
    denudsym(:,ixr)=den_re(:,ixr)-den_re(:,-ixr)
   enddo
  endif

  udsymmetry=SUM(abs(denudsym))/(Nxr-1)/Nxa
  udsymdiff=SUM(denudsym)/(Nxr-1)/Nxa

  write(*,*)'udsymmetry, udsymdiff:',udsymmetry,udsymdiff
  write(fu_diagudsym,*)it,udsymmetry,udsymdiff

  do ixr=1,Nxr-1
   write(fu_denudsym,*)(DBLE(denudsym(ixa,ixr)),ixa=-Nxa2,Nxa2-1)
  enddo

!!!!!!! begin brent's file reading !!!!! 
  !get to next timestep
  do ixr=1,2  
   read(fu_2dxre,*)
  enddo
!!!!!!!!!! end brent's file reading !!!!!

!!!!!! begin arnau's file reading !!!!
! read(fu_2dxre,*)
!!!!! end arnau's file reading

 enddo !it

 close(fu_2dxre)
 close(fu_denudsym)
 close(fu_denlrsym)
 close(fu_diagudsym)
 close(fu_diaglrsym)

 enddo !istate
 enddo !ireim

end program procdenextra

