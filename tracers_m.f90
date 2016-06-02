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
    
module tracers_m

implicit none

private
public ngas,ntrac,ntracmax,nllp
public ilt,jlt,klt,ngasmax
public tr,traver,acloss_g !,gasmin
public trpm,npm
public tracers_init,tracers_end

! parameters should be controlled from namelist
integer, parameter :: nllp=0
integer, save :: ngas=0
integer, save :: ntrac=0
integer, save :: ntracmax=1
integer, save :: ngasmax=1
integer, save :: ilt=1
integer, save :: jlt=1
integer, save :: klt=1
real, dimension(:,:,:), allocatable, save :: tr,traver
real, dimension(:,:,:), allocatable, save :: trpm
real, dimension(:), allocatable, save :: acloss_g
integer, dimension(:), allocatable, save :: npm

contains

subroutine tracers_init(il,jl,kl,iextra)

implicit none

integer, intent(in) :: il,jl,kl,iextra

ntrac=ngas+nllp
ntracmax=max(ntrac,1) ! ntracmax >= 1
ngasmax=max(ngas,1)   ! ngasmax >= 1

! old trick for common blocks
! now just set ilt to il or 1, etc
!ilt=il**npwr
!jlt=jl**npwr
!klt=kl**npwr

ilt=il
jlt=jl
klt=kl

allocate(tr(ilt*jlt+iextra,klt,ntracmax),traver(ilt*jlt,klt,ntrac))

return
end subroutine tracers_init

subroutine tracers_end

implicit none

deallocate(tr,traver)

return
end subroutine tracers_end

end module tracers_m