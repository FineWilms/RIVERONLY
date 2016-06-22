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
! nurban>0            Use urban scheme
    
subroutine sflux(nalpha)
      
use arrays_m                       ! Atmosphere dyamics prognostic arrays
use ateb                           ! Urban
use cable_ccam, only : sib4        ! CABLE interface
use cc_mpi                         ! CC MPI routines
use diag_m                         ! Diagnostic routines
use extraout_m                     ! Additional diagnostics
use gdrag_m                        ! Gravity wave drag
!~ use liqwpar_m                      ! Cloud water mixing ratios
use map_m                          ! Grid map arrays
use mlo                            ! Ocean physics and prognostic arrays
use morepbl_m                      ! Additional boundary layer diagnostics
!~ use nharrs_m                       ! Non-hydrostatic atmosphere arrays
use nsibd_m                        ! Land-surface arrays
use pbl_m                          ! Boundary layer arrays
use permsurf_m                     ! Fixed surface arrays
!~ use prec_m                         ! Precipitation
use river                          ! River routing
use savuvt_m                       ! Saved dynamic arrays
use screen_m                       ! Screen level diagnostics
use sigs_m                         ! Atmosphere sigma levels
use soil_m                         ! Soil and surface data
use soilsnow_m                     ! Soil, snow and surface data
use vecsuv_m                       ! Map to cartesian coordinates
use vvel_m                         ! Additional vertical velocity
use work2_m                        ! Diagnostic arrays
use work3_m                        ! Mk3 land-surface diagnostic arrays
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
integer, parameter :: ntss_sh=0 ! 0 for original, 3 for **3, 4 for **4
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

if (.not.allocated(plens).and.nplens/=0) then
  allocate(plens(ifull))
  plens=0.
end if
if (.not.allocated(taftfh).and.(nsib==3.or.nsib==5)) then
  allocate(taftfh(ifull))
  allocate(taftfhg(ifull))
  taftfh(:)=.05        ! just a diag default for sea points
  taftfhg(:)=7.e-4     ! just a diag default for sea points
end if
      
ri_max=(1./fmroot -1.)/bprm    ! i.e. .14641
zologbgin=log(zmin/zobgin)     ! pre-calculated for all except snow points
ztv=exp(vkar/sqrt(chn10))/10.  ! proper inverse of ztsea
z1onzt=300.*rdry*(1.-sig(1))*ztv/grav
chnsea=(vkar/log(z1onzt))**2   ! should give .00085 for csiro9
oldrunoff(:)=runoff(:)
zo=999.        ! dummy value
factch=999.    ! dummy value
taux=0.        ! dummy value
tauy=0.        ! dummy value
gamm=3.471e+05 ! dummy value
root=0.        ! dummy value
denha=0.       ! dummy value
denma=0.       ! dummy value
fm=0.          ! dummy value

!     using av_vmod (1. for no time averaging)
!      *****  check next comment
!       sflux called at beginning of time loop, hence savu, savv

srcp = sig(1)**(rdry/cp)
ga(:) = 0.              !  for ocean points in ga_ave diagnostic
theta(:) = t(1:ifull,1)/srcp
rho(:) = ps(1:ifull)/(rdry*tss(:))
uav(:) = av_vmod*u(1:ifull,1) + (1.-av_vmod)*savu(:,1)   
vav(:) = av_vmod*v(1:ifull,1) + (1.-av_vmod)*savv(:,1)  
vmod(:) = sqrt(uav(:)**2+vav(:)**2)  ! i.e. vmod for tss_sh
vmag(:) = max( vmod(:), vmodmin )    ! vmag used to calculate ri
if ( ntsur/=7 ) vmod(:) = vmag(:)      ! gives usual way

!--------------------------------------------------------------
call START_LOG(sfluxwater_begin)
if (nmlo==0) then                                                                                ! sea
  if ( nriver==1 ) then                                                                          ! river
    where ( .not.land(1:ifull) )                                                                 ! river
      watbdy(1:ifull) = 0. ! water enters ocean and is removed from rivers                       ! river
    end where                                                                                    ! river
  end if                                                                                         ! river
  
elseif (abs(nmlo)>=1.and.abs(nmlo)<=9) then                                                      ! MLO
                                                                                                 ! MLO
  if (nmaxpr==1) then                                                                            ! MLO
    if (myid==0) then                                                                            ! MLO
      write(6,*) "Before MLO mixing"                                                             ! MLO
    end if                                                                                       ! MLO
    call ccmpi_barrier(comm_world)                                                               ! MLO
  end if                                                                                         ! MLO
  if (abs(nmlo)==1) then                                                                         ! MLO
    ! Single column                                                                              ! MLO
    ! set free surface to zero when water is not conserved                                       ! MLO
    neta=0.                                                                                      ! MLO
    call mloimport(4,neta,0,0)                                                                   ! MLO
  end if                                                                                         ! MLO

  ! pan evaporation diagnostic                                                                   ! MLO

  do ip=1,ipland                                                                                 ! MLO
    iq=iperm(ip)                                                                                 ! MLO
    ri(iq)=min(grav*zmin*(1.-tpan(iq)*srcp/t(iq,1))/vmag(iq)**2,ri_max)                          ! MLO
    if(ri(iq)>0.)then                                                                            ! MLO
      fh(iq)=vmod(iq)/(1.+bprm*ri(iq))**2                                                        ! MLO
    else                                                                                         ! MLO
      root=sqrt(-ri(iq)*zmin/panzo)                                                              ! MLO
      denha=1.+chs*2.*bprm*sqrt(panzo*ztv)*chnsea*root                                           ! MLO
      fh(iq)=vmod(iq)-vmod(iq)*2.*bprm *ri(iq)/denha                                             ! MLO
    endif                                                                                        ! MLO
    epan(iq)=rho(iq)*chnsea*hl*fh(iq)*(qsttg(iq)-qg(iq,1))                                       ! MLO
  end do                                                                                         ! MLO
                                                                                                 ! MLO
  ! inflow and outflow model for rivers                                                          ! MLO
  if ( abs(nmlo)>=2 ) then                                                                       ! MLO
    if ( .not.allocated(outflowmask) ) then                                                      ! MLO
      allocate( outflowmask(1:ifull) )                                                           ! MLO
      call riveroutflowmask(outflowmask)                                                         ! MLO
    end if                                                                                       ! MLO
    neta(1:ifull) = 0.                                                                           ! MLO
    call mloexport(4,neta,0,0)                                                                   ! MLO
    oldneta(1:ifull) = neta(1:ifull)                                                             ! MLO
    where ( outflowmask(1:ifull) )                                                               ! MLO
      neta(1:ifull) = min( neta(1:ifull), max( 0.001*watbdy(1:ifull), 0. ) )                     ! MLO
    end where                                                                                    ! MLO
    where ( .not.land(1:ifull) )                                                                 ! MLO
      dumw(1:ifull) = 1000.*(neta(1:ifull)-oldneta(1:ifull))/dt + watbdy(1:ifull)/dt             ! MLO
    elsewhere                                                                                    ! MLO
      dumw(1:ifull) = 0.                                                                         ! MLO
    end where                                                                                    ! MLO
    watbdy(1:ifull) = watbdy(1:ifull) - dt*dumw(1:ifull)                                         ! MLO
  else                                                                                           ! MLO
    dumw(1:ifull) = 0.                                                                           ! MLO
  end if                                                                                         ! MLO

  ! stuff to keep tpan over land working                                                         ! MLO
  rid=min(grav*zmin*(1.-tpan*srcp/t(1:ifull,1))/vmag**2,ri_max)                                  ! MLO
  where (rid>0.)                                                                                 ! MLO
    fhd=vmod/(1.+bprm*rid)**2                                                                    ! MLO
  elsewhere                                                                                      ! MLO
    fhd=vmod-vmod*2.*bprm*rid/(1.+chs*2.*bprm*sqrt(panzo*ztv)*chnsea*sqrt(-rid*zmin/panzo))      ! MLO
  end where                                                                                      ! MLO
                                                                                                 ! MLO
  where ( .not.land(1:ifull) )                                                                   ! MLO
    snowd=snowd*1000.                                                                            ! MLO
    ga=0.                                                                                        ! MLO
    ustar=sqrt(sqrt(taux*taux+tauy*tauy)/rho)                                                    ! MLO
    tpan=tgg(:,1)                                                                                ! MLO
    factch=sqrt(zo/zoh)                                                                          ! MLO
    sno=sno+conds                                                                                ! MLO
    grpl=grpl+condg                                                                              ! MLO
    ! This cduv accounts for a moving surface                                                    ! MLO
    cduv=sqrt(ustar*ustar*cduv) ! cduv=cd*vmod                                                   ! MLO
    cdtq=cdtq*vmod                                                                               ! MLO
  elsewhere                                                                                      ! MLO
    fg=rho*chnsea*cp*fhd*(tpan-theta)                                                            ! MLO
    ga=sgsave-rgsave-5.67e-8*tpan**4-panfg*fg                                                    ! MLO
    tpan=tpan+ga*dt/(4186.*.254*1000.)                                                           ! MLO
  endwhere                                                                                       ! MLO
 
  if (nmaxpr==1) then                                                                            ! MLO
    if (myid==0) then                                                                            ! MLO
      write(6,*) "After MLO mixing"                                                              ! MLO
    end if                                                                                       ! MLO
    call ccmpi_barrier(comm_world)                                                               ! MLO
  end if                                                                                         ! MLO
end if                                                                                           ! PCOM
call END_LOG(sfluxwater_end)                                                                     ! PCOM
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
      factch(1:ifull) = sqrt(zo(1:ifull)/zoh(1:ifull))                                           ! cable 
      !~ qsttg(1:ifull) = qsat(ps(1:ifull),tss(1:ifull))                                            ! cable
      taux(1:ifull) = rho(1:ifull)*cduv(1:ifull)*u(1:ifull,1)                                    ! cable
      tauy(1:ifull) = rho(1:ifull)*cduv(1:ifull)*v(1:ifull,1)                                    ! cable
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
if ( abs(nmlo)>=2 .or. nriver==1 ) then
  newrunoff=runoff-oldrunoff
  watbdy(1:ifull)=watbdy(1:ifull)+newrunoff ! runoff in mm
end if
!***  end of surface updating loop
! ----------------------------------------------------------------------

return
end subroutine sflux
