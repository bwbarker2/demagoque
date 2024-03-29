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

program compareAB
 !! compares Arnau's and Brent's code in a very non-modular way
 implicit none

 integer, parameter :: Nxa=100
 integer, parameter :: Nxr=100
 integer, parameter :: Nxa2=Nxa/2
 integer, parameter :: Nxr2=Nxr/2

 integer :: fu_brent, fu_arnau, fu_densub !file unit numbers
 real*8, dimension(-Nxa2:Nxa2-1,-Nxr2:Nxr2) :: denbrent, denarnau,densub
 real*8 :: summ,summ2   !to get average difference
 integer :: numm  !to get number to divide by for average
 real*8 :: xava,xavb   ! <|x|> ( expectation of abs(x) )
 real*8 :: xxava, xxavb  ! <x> (expectation of x)

 complex*16, dimension(1:Nxa2-1,-Nxr2:Nxr2-1) :: denlrsym
 real*8 :: lrsymmetry

 integer :: ixa,ixr,it

 fu_brent=1
 fu_arnau=2
 fu_densub=3

 open(unit=fu_brent,file='results.25x25/2dxre.dat',status='old')
 open(unit=fu_arnau,file='results.arnau/mfevol/2dxre.dat',status='old')
 open(unit=fu_densub,file='results/densub.dat')

!loop through timesteps
do it=0,4
 write(fu_densub,*)'# time=',it*10,'fm/c'

 !! read brent's, store in denbrent !!

 read(fu_brent,*) !comment line

 !don't read first quarter of data (redundant info)
 do ixr=1,Nxr2  
  read(fu_brent,*)
 enddo
 
 !fill denbrent
 do ixr=-Nxr2,Nxr2-1
  read(fu_brent,*)(denbrent(ixa,ixr),ixa=-Nxa2,Nxa2-1)
 enddo
! write(*,*)denbrent(:,0)

 !don't read first line, it's a duplicate of the last line
 read(fu_arnau,*)
 !fill denarnau
 do ixa=Nxa2-1,-Nxa2,-1
  read(fu_arnau,*)(denarnau(ixa,ixr),ixr=Nxr2,-Nxr2,-1)
 enddo

! write(*,*)denarnau(:,-1)

! !divide denarnau by the 1D->3D factor he used
! denarnau=denarnau/0.48680747073350783d0
! denarnau=denarnau/0.5069
! denarnau=denarnau*0.25d0

 !calculate left-right symmetry, useful for initially symmetry system. lrsymmetry should be 0 for perfect symmetry. denlrsym is the matrix containing the subtraction of right side minus the left side
! lrsymmetry=0d0
 do ixa=1,Nxa2-1
  denlrsym(ixa,:)=denbrent(ixa,-Nxr2:Nxr2-1)-denbrent(-ixa,-Nxr2:Nxr2-1)
 enddo

 lrsymmetry=SUM(denlrsym**2)

 write(*,*)'lrsymmetry:',lrsymmetry

 !calculate xav for each system
 xava=0d0
 xavb=0d0
 do ixa=-Nxr2,Nxr2-1
  xava=xava+abs(ixa)*denarnau(ixa,0)*0.5d0
  xavb=xavb+abs(ixa)*denbrent(ixa,0)*0.5d0
 enddo

 xxava=0d0
 xxavb=0d0
 do ixa=-Nxr2,Nxr2-1
  xxava=xxava+ixa*denarnau(ixa,0)*0.5d0
  xxavb=xxavb+ixa*denbrent(ixa,0)*0.5d0
 enddo


 write(*,*)'xava,xavb=',xava,xavb,(xava-xavb)/(xava+xavb)*2
 write(*,*)'xxava,xxavb=',xxava,xxavb,(xxava-xxavb)/(xxava+xxavb)*2

 !calculate densub
 do ixa=-Nxa2,Nxa2-1
  do ixr=-Nxr2,Nxr2-1
!   if (denarnau(ixa,ixr).lt.1d-15.AND.denbrent(ixa,ixr).lt.1d-15) then
!    densub(ixa,ixr)=0d0
!   elseif(denarnau(ixa,ixr)+denbrent(ixa,ixr).eq.0d0) then
!    write(*,*)'whoa'
!   else
    densub(ixa,ixr)=(denarnau(ixa,ixr)-denbrent(ixa,ixr)) !&
!                    /(abs(denarnau(ixa,ixr))+abs(denbrent(ixa,ixr)))
!    if(abs(densub(ixa,ixr))>2)write(*,*)densub(ixa,ixr)
!   endif
  enddo
 enddo

 summ=0d0
 summ2=0d0
 numm=0
 do ixa=-Nxa2,Nxa2-1
  do ixr=-Nxr2,Nxr2-1
    summ=summ+densub(ixa,ixr)
!    write(*,*)densub(ixa,ixr)
    summ2=summ2+densub(ixa,ixr)*densub(ixa,ixr)
    numm=numm+1
  enddo
 enddo

! write(*,*)numm,summ
 summ=summ/numm
 summ2=sqrt(summ2)/numm

 write(*,*)summ,summ2

! do ixa=-Nxa2,Nxa2-1
!  do ixr=-Nxr2,Nxr2-1
!   if(densub(ixa,ixr).ne.1d0)then
!    densub(ixa,ixr)=densub(ixa,ixr)/summ
!   endif
!  enddo
! enddo

! densub=denarnau/denbrent

 do ixr=-Nxr2,Nxr2-1
  write(fu_densub,*)(densub(ixa,ixr),ixa=-Nxa2,Nxa2-1)
 enddo

 !get to next timestep
 do ixr=1,Nxr2+2  
  read(fu_brent,*)
 enddo

 read(fu_arnau,*)

 !give 2 blank lines in densub
 write(fu_densub,*)
 write(fu_densub,*)

enddo !cycle through timesteps

 close(fu_brent)
 close(fu_arnau)
 close(fu_densub)

! write(*,*)densub(:,-2)

end program compareAB
