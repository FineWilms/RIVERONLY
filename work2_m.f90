! Conformal Cubic Atmospheric Model
    
! Copyright 2015 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
! This file is part of the Conformal Cubic Atmospheric Model (CCAM)
!
! CCAM is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! CCAM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with CCAM.  If not, see <http://www.gnu.org/licenses/>.

!------------------------------------------------------------------------------
    
module work2_m

implicit none

private
public zoh,zoq,qsttg,wetfac
public zo,theta
public dgdtg
public work2_init,work2_end

real, dimension(:), allocatable, save :: zoh,zoq,qsttg,wetfac
real, dimension(:), allocatable, save :: zo,theta
real, dimension(:), allocatable, save :: dgdtg

contains

subroutine work2_init(ifull,iextra,kl,nsib)

implicit none

integer, intent(in) :: ifull,iextra,kl,nsib

allocate(zoh(ifull),zoq(ifull),qsttg(ifull),wetfac(ifull))
allocate(zo(ifull),theta(ifull))
zo=0.
zoh=0.
zoq=0.
qsttg=0.
wetfac=1.
if (nsib==3.or.nsib==5) then
  allocate(dgdtg(ifull))
  dgdtg=0.
end if

return
end subroutine work2_init

subroutine work2_end

implicit none

deallocate(zoh,zoq,qsttg,wetfac)
deallocate(zo,theta)
if (allocated(dgdtg)) then
  deallocate(dgdtg)
end if

return
end subroutine work2_end

end module work2_m