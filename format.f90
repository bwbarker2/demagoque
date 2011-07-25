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

module format
  implicit none

  integer :: dummy
  
  ! ... Format for output files
90 format(i6)
91 format(e16.8)
92 format(2e16.8)
93 format(3e16.8)
94 format(4e16.8)
95 format(5e16.8)
98 format(2e11.3,3x,2e11.3,3x,2e11.3,3x,2e11.3)

! ... FORMAT OF THE WRITING FILES
200 format(/,'# Time=',1d15.5)
202 format(/,/,'# Time=',1d15.5)

201 format('Time',1x,I4,'/',I4,3x,6g16.6)

143 format('# Initial values',/,'#',8x,'EKINX',14x,'EKINK' &
         ,14x,'EPOT',14x,'ETOT',14x,'DENX',14x,'DENK' &
         ,14x,'MOMK',/,'#',7e18.9)
243 format(/,'# Relative differences',/ &
         ,'#',7x,'TIME',11x,'EKINX/EKINX0',6x,'EKINK/EKINK0' &
         ,6x,'EPOT/EPOT0',8x,'ETOT/ETOT0',8x,'DENX/DENX0' &
         ,8x,'DENK/DENK0',8x,'MOMK/MOMK0')
343 format(8e18.9) 

144 format('# Initial values',/,'#',8x,'TIME',14x,'EKINX' &
         ,13x,'EKINK',13x,'EPOT' &
         ,14x,'ETOT',14x,'DENX',14x,'DENK',12x,'TOTMOM',/ &
         ,8e18.9)
344 format(8e18.9) 

145 format('#',5x,'TIME [fm]',8x,'ETOT/A [MeV]' &
         ,8x,'ETOTX/A',11x,'EKIN/A' &
         ,12x,'EKINX/A',11x,'EPOT/A',12x,'TOTMOM')
345 format(8e18.9) 

347 format(3e18.9) 

146 format('#',4x,'TIME [fm]',11x,'XCM [fm]' &
         ,6x,'ELONGATION [fm]',4x,'XCM^2 [fm^2]',8x,'<X> [fm^2]')

919 format(1000e16.8)

end module format
