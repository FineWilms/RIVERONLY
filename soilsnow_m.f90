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
    
module soilsnow_m

implicit none

private
public tggsn,tgg,wb,wbice,smass,ssdn,ssdnn,snowd 
public osnowd,snage,sno,grpl,gflux,sgflux,snowflx,otgsoil 
public runoff,albvisnir
public fracice,sicedep
public isflag
public soilsnow_init,soilsnow_end

integer, dimension(:), allocatable, save :: isflag
real, dimension(:), allocatable, save :: ssdnn,snowd
real, dimension(:), allocatable, save :: osnowd,snage,sno,grpl,gflux,sgflux,snowflx,otgsoil
real, dimension(:), allocatable, save :: runoff
real, dimension(:), allocatable, save :: fracice,sicedep
real, dimension(:,:), allocatable, save :: tggsn,tgg,wb,wbice,smass,ssdn
real, dimension(:,:), allocatable, save :: albvisnir

contains

subroutine soilsnow_init(ifull,iextra,kl,ms,nsib)

implicit none

integer, intent(in) :: ifull,iextra,kl,ms,nsib

allocate(tggsn(ifull,3),tgg(ifull,ms),wb(ifull,ms),wbice(ifull,ms))
allocate(smass(ifull,3),ssdn(ifull,3),ssdnn(ifull),snowd(ifull))
allocate(snage(ifull),sno(ifull),grpl(ifull),gflux(ifull))
allocate(runoff(ifull),albvisnir(ifull,2))
allocate(fracice(ifull),sicedep(ifull))
allocate(isflag(ifull))
if (nsib==3.or.nsib==5) then
  allocate(sgflux(ifull))
  allocate(osnowd(ifull),otgsoil(ifull),snowflx(ifull))
end if

return
end subroutine soilsnow_init

subroutine soilsnow_end

implicit none

deallocate(tggsn,tgg,wb,wbice,smass,ssdn,ssdnn,snowd)
deallocate(snage,sno,grpl,gflux)
deallocate(runoff,albvisnir)
deallocate(fracice,sicedep)
deallocate(isflag)
if (allocated(sgflux)) then
  deallocate(sgflux)
  deallocate(osnowd,otgsoil,snowflx)
end if

return
end subroutine soilsnow_end

end module soilsnow_m
