module bstring
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

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine findFirstWord(line,delimiters,istart,iend)
 !! returns indices denoting the extent of the first string found in the string inarr  -- NOTE: not finished yet
 implicit none

 character(len=*), intent(in) :: line !input array
 character(len=*), intent(in) :: delimiters !characters that are not in 'words'
 integer, intent(out) :: istart, iend !index of start, end of word in 'line'

 integer :: i, lenline

 !if there is no word, set these values
 istart=0
 iend=0

 lenline=len(line)

 if(lenline==0)return

 do i=1,lenline
  !if the i'th character is not a delimiter, then start the word here
  if(index(delimiters,line(i:i))==0) then
   istart=i
   exit
  endif
 enddo

 if(istart==0)return  !if no word exists, exit without looking for end of word

 !if no delimiter is found in following do-loop, then the word extends to end of line
 iend=lenline

 if(lenline==1)return  !if word exists and is only 1 character, then iend is set

 do i=istart+1,lenline
  !if the i'th character is a delimiter, then end the word a character before
  if(index(delimiters,line(i:i)).NE.0) then
   iend=i-1
   exit
  endif
 enddo

end subroutine findFirstWord

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

character function getFirstNonBlankChar(chararr)
 implicit none

 character(len=*), intent(in) :: chararr
 
 integer :: i

 do i=1,len(chararr)
  getFirstNonBlankChar=chararr(i:i)
  if(getFirstNonBlankChar.ne.' ')exit
 enddo

end function getFirstNonBlankChar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

logical function isComment(line)
 implicit none

 character(len=*) :: line
 integer :: linelen

 linelen=len(line)

 isComment=(getFirstNonBlankChar(line)=='!')

end function isComment


end module bstring

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!program testbstring
! call test_bstring
!end program testbstring

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_bstring

 write(*,*)'Testing module bstring'

 call test_findFirstWord

end subroutine test_bstring

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_findFirstWord
 use bstring
 implicit none

 character(len=0) :: zerolen
 character(len=1) :: onelen
 character(len=10) :: tenlen

 integer :: ibegin,istop

 write(*,*)'> Testing subroutine findFirstWord'

 call findFirstWord(zerolen,' ',ibegin,istop)
 if(ibegin.ne.0.or.istop.ne.0)then
  write(*,*)'!!! zero-length line fails, istart,iend=',ibegin,istop
 endif

 onelen=' '
 call findFirstWord(onelen,' ',ibegin,istop)
 if(ibegin.ne.0.or.istop.ne.0)then
  write(*,*)'!!! blank one-length line fails, istart,iend=',ibegin,istop
 endif

 onelen='a'
 call findFirstWord(onelen,' ',ibegin,istop)
 if(ibegin.ne.1.or.istop.ne.1)then
  write(*,*)'!!! filled one-length line fails, istart,iend=',ibegin,istop
 endif

 tenlen="abcd      "
 call findFirstWord(tenlen,' ',ibegin,istop)
 if(ibegin.ne.1.or.istop.ne.4)then
  write(*,*)'!!! 10-length line fails:'
  write(*,*)'tenlen=',tenlen
  write(*,*)'istart,iend=',ibegin,istop
 endif

 tenlen=" "
 call findFirstWord(tenlen,' ',ibegin,istop)
 if(ibegin.ne.0.or.istop.ne.0)then
  write(*,*)'!!! 10-length line fails:'
  write(*,*)'tenlen=',tenlen
  write(*,*)'istart,iend=',ibegin,istop
 endif

 tenlen="    abcdef"
 call findFirstWord(tenlen,' ',ibegin,istop)
 if(ibegin.ne.5.or.istop.ne.10)then
  write(*,*)'!!! 10-length line fails:'
  write(*,*)'tenlen=',tenlen
  write(*,*)'istart,iend=',ibegin,istop
 endif

 tenlen="abc defg"
 call findFirstWord(tenlen,' ',ibegin,istop)
 if(ibegin.ne.1.or.istop.ne.3)then
  write(*,*)'!!! 10-length line fails:'
  write(*,*)'tenlen=',tenlen
  write(*,*)'istart,iend=',ibegin,istop
 endif

end subroutine test_findFirstWord

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


