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

! CCAM interface for surface flux routines. Includes standard land-surface scheme,
! prescribed SSTs and sea-ice, CABLE interface, urban interface and MLO interface.

! nsib=7              CABLE land-surface scheme with CCAM diagnostics
! nmlo=0              Prescriped SSTs and sea-ice with JLM skin enhancement    
subroutine sflux(nalpha)
      
use arrays_m                       ! Atmosphere dyamics prognostic arrays
use cable_ccam, only : sib4        ! CABLE interface
use cc_mpi                         ! CC MPI routines
use diag_m                         ! Diagnostic routines
use map_m                          ! Grid map arrays
use morepbl_m                      ! Additional boundary layer diagnostics
use nsibd_m                        ! Land-surface arrays
use pbl_m                          ! Boundary layer arrays
use permsurf_m                     ! Fixed surface arrays
use river                          ! River routing
use sigs_m                         ! Atmosphere sigma levels
use soil_m                         ! Soil and surface data
use soilsnow_m                     ! Soil, snow and surface data
use vecsuv_m                       ! Map to cartesian coordinates
use work2_m                        ! Diagnostic arrays
use xyzinfo_m                      ! Grid coordinate arrays
      
implicit none
    
include 'newmpar.h'                ! Grid parameters
include 'const_phys.h'             ! Physical constants
include 'parm.h'                   ! Model configuration
include 'parmgeom.h'               ! Coordinate data
include 'parmsurf.h'               ! Surface parameters
include 'soilv.h'                  ! Soil parameters

integer iq,k,it,ip
integer, intent(in) :: nalpha
real ri_max,zologbgin,ztv,z1onzt,chnsea
real srcp,afrootpan,es,constz,drst
real xx,consea,afroot,fm,con,dtsol,daf
real con1,den,dden,dfm,root,denma,denha
real conh,conw,zminlog,ri_ice,zoice,zologice
real epotice,qtgnet,eg1,eg2,deg,b1
real gbot,deltat,esatf,zobg,zologbg,zologx
real afland,aftlandg,fhbg,rootbg,denhabg
real thnew,thgnew,thnewa,qtgair,aftland
real thgnewa,ri_tmp,fh_tmp,factchice
real, dimension(:), allocatable, save :: taftfh,taftfhg
real, dimension(:), allocatable, save :: plens
real, dimension(ifull) :: vmag,charnck,taftfhg_temp
real, dimension(ifull) :: zonx,zony,zonz,costh
real, dimension(ifull) :: sinth,uzon,vmer,azmin
real, dimension(ifull) :: uav,vav
real, dimension(ifull) :: oldrunoff,newrunoff,rid,fhd
real, dimension(ifull) :: fgf,rgg,fev,af,dirad,dfgdt,factch
real, dimension(ifull) :: degdt,cie,aft,fh,ri,gamm,rho
real, dimension(ifull) :: dumsg,dumrg,dumx,dums,dumw,tv
real, dimension(ifull) :: neta, oldneta
logical, dimension(:), allocatable, save :: outflowmask

integer, parameter :: nblend=0  ! 0 for original non-blended, 1 for blended af
integer, parameter :: ntest=0   ! ntest= 0 for diags off; ntest= 1 for diags on
real, parameter :: bprm=5.,cms=5.,chs=2.6,vkar=.4
real, parameter :: d3=2.5
real, parameter :: cgsoil=1000.,gksoil=.300e-6,rhog=1600.
real, parameter :: d1land=.03
real, parameter :: fmroot=.57735     ! was .4 till 7 Feb 1996

!     stability dependent drag coefficients using Louis (1979,blm) f'
!     n.b. cduv, cdtq are returned as drag coeffs mult by vmod
!          (cduv=cduv*vmod; cdtq=cdtq*vmod)

!     t, u, v, qg are current values
!     tss is surface temperature
!     dw is soil wetness availability (1. for ocean) - not needed here
!     fg is sensible heat flux (was h0)
!     eg is latent heat flux (was wv)
!     dfgdt is dfgdt (was csen in surfupa/b)
!     degdt is degdt (was ceva in surfupa/b)

      
ri_max=(1./fmroot -1.)/bprm    ! i.e. .14641
zologbgin=log(zmin/zobgin)     ! pre-calculated for all except snow points
ztv=exp(vkar/sqrt(chn10))/10.  ! proper inverse of ztsea
z1onzt=300.*rdry*(1.-sig(1))*ztv/grav
chnsea=(vkar/log(z1onzt))**2   ! should give .00085 for csiro9
oldrunoff(:)=runoff(:)
zo=999.        ! dummy value
factch=999.    ! dummy value
!~ taux=0.        ! dummy value
!~ tauy=0.        ! dummy value
gamm=3.471e+05 ! dummy value
root=0.        ! dummy value
denha=0.       ! dummy value
denma=0.       ! dummy value
fm=0.          ! dummy value

!     using av_vmod (1. for no time averaging)
!      *****  check next comment
!       sflux called at beginning of time loop, hence savu, savv

rho(:) = ps(1:ifull)/(rdry*tss(:))

!--------------------------------------------------------------
call START_LOG(sfluxwater_begin)
                                                                            ! sea
  if ( nriver==1 ) then                                                                          ! river
    where ( .not.land(1:ifull) )                                                                 ! river
      watbdy(1:ifull) = 0. ! water enters ocean and is removed from rivers                       ! river
    end where                                                                                    ! river
  end if                                                                                         ! river
                                                                                        
call END_LOG(sfluxwater_end)                                                                     
!--------------------------------------------------------------      
call START_LOG(sfluxland_begin)    

select case(nsib)                                                                                ! land
  case(7)                                                                                        ! cable
    if (nmaxpr==1) then                                                                          ! cable
      if (myid==0) then                                                                          ! cable
        write(6,*) "Before CABLE"                                                                ! cable
      end if                                                                                     ! cable
      call ccmpi_barrier(comm_world)                                                             ! cable
    end if                                                                                       ! cable
    ! call cable                                                                                 ! cable
    call sib4                                                                                    ! cable
    ! update remaining diagnostic arrays                                                         ! cable
    where ( land(1:ifull) )                                                                      ! cable
      sno(1:ifull) = sno(1:ifull) + conds(1:ifull)                                               ! cable
      grpl(1:ifull) = grpl(1:ifull) + condg(1:ifull)                                             ! cable
    end where                                                                                    ! cable
    if (nmaxpr==1) then                                                                          ! cable
      if (myid==0) then                                                                          ! cable
        write(6,*) "After CABLE"                                                                 ! cable
      end if                                                                                     ! cable
      call ccmpi_barrier(comm_world)                                                             ! cable
    end if                                                                                       ! cable
                                                                                                 ! cable
  case DEFAULT                                                                                   ! land
    write(6,*) "ERROR: Unknown land-use option nsib=",nsib                                       ! land
    call ccmpi_abort(-1)                                                                         ! land
                                                                                                 ! land
end select                                                                                       ! land
call END_LOG(sfluxland_end)                                                                      ! land
!----------------------------------------------------------
! Update runoff for river routing
if ( nriver==1 ) then
  newrunoff=runoff-oldrunoff
  watbdy(1:ifull)=watbdy(1:ifull)+newrunoff ! runoff in mm
end if
!***  end of surface updating loop
! ----------------------------------------------------------------------

return
end subroutine sflux