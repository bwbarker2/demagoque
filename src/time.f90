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

!> time parameters
MODULE time
  USE prec_def
  implicit none

  private
  public :: it,t,Nt,firstOutput,lastOutput,time_initialize,time_getString

  INTEGER     :: it     !< current iteration of time
  REAL (Long) :: t      !< current time during time evolution [fm/c]
  INTEGER     :: Nt     !< number of timesteps in current mode (adiabatic or
                        !!  time)

  !> array of strings for each time for use in output filenames
  character(len=20), dimension(:), allocatable :: timestrings
  integer :: timestrlen !< length of timestring

  logical     :: firstOutput  !< first output of a given evolution
  logical     :: lastOutput   !< last output of given evolution

contains

 pure function time_getString(tindex) result(timestr)
  character(len=timestrlen) :: timestr

  integer, intent(in) :: tindex

  timestr=trim(timestrings(tindex))

 end function time_getString

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11

 subroutine time_initialize()
  use bstring
  use input_parameters
  implicit none

  integer :: big,iit,ii,istart,iend
  real(Long) :: ttemp
  character(len=:), allocatable :: timestr

  !> First we initialize the timestrings. This involves parsing all the output
  !! times, padding the beginning of the string with zeroes, and cutting the
  !! extra zeroes off the end.

  allocate(timestrings(0:Nt/ntime))

  big=10
  do while (Nt*delt>big)
   big=big*10
  enddo
  big=big/10

  do iit=0,Nt/ntime
   ttemp = iit*delt*ntime

   write(timestrings(iit),'(F20.10)')ttemp
   timestrings(iit)=adjustl(timestrings(iit))

   if(ttemp<epzero) ttemp=1._Long

   do while(ttemp<big)
    timestrings(iit)='0'//timestrings(iit)(1:len(timestrings(iit))-1)
    ttemp=ttemp*10._Long
   end do

   ii=index(timestrings(iit),'.')
   timestrings(iit)(ii:ii)='_'

!write(*,*)'iit,timestrings=',timestrings(iit)

  enddo

  ! now trim the excess zeroes from the right side of the string

  ! the longest string can be found at t=delt
  allocate(character(len=len(timestrings(1))) :: timestr)

  write(timestr,'(F20.10)')delt
  timestr=adjustl(timestr)

  !iend-istart is the number of decimal places to use
  iend=verify(timestr,' 0',back=.true.)
  istart=index(timestr,'.')

!  write(*,*)'time_initialize:ii,timestr=',trim(timestr),istart,iend

  iend=index(timestrings(1),'_')+iend-istart

  timestrlen=iend

  do ii=0,size(timestrings)-1
   timestrings(ii)=timestrings(ii)(1:iend)
!write(*,*)'time_initialize:ii,timestrings=',ii,timestrings(ii)
  enddo

 end subroutine time_initialize

END MODULE time

