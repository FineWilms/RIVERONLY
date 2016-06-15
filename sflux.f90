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
      
! nsib=3              Standard land-surface scheme with SIB and Gratz data
! nsib=5              Standard land-surface scheme with MODIS data
! nsib=6              CABLE land-surface scheme with CABLE diagnostics
! nsib=7              CABLE land-surface scheme with CCAM diagnostics
! nmlo=0              Prescriped SSTs and sea-ice with JLM skin enhancement
! nmlo>0 and mlo<=9   KPP ocean mixing
! nmlo>9              Use external PCOM ocean model
! nurban>0            Use urban scheme
    
subroutine sflux(nalpha)
      
use arrays_m                       ! Atmosphere dyamics prognostic arrays
use ateb                           ! Urban
use cable_ccam, only : sib4        ! CABLE interface
use cc_mpi                         ! CC MPI routines
use diag_m                         ! Diagnostic routines
!~ use estab                          ! Liquid saturation function
use extraout_m                     ! Additional diagnostics
use gdrag_m                        ! Gravity wave drag
use liqwpar_m                      ! Cloud water mixing ratios
use map_m                          ! Grid map arrays
use mlo                            ! Ocean physics and prognostic arrays
!~ use mlodynamics                    ! Ocean dynamics routines
use morepbl_m                      ! Additional boundary layer diagnostics
use nharrs_m                       ! Non-hydrostatic atmosphere arrays
use nsibd_m                        ! Land-surface arrays
use pbl_m                          ! Boundary layer arrays
use permsurf_m                     ! Fixed surface arrays
use prec_m                         ! Precipitation
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

if (diag.or.ntest==1) then
  if (mydiag) then
    if (land(idjd)) then
      write(6,*) 'entering sflux ktau,nsib,ivegt,isoilm,land ',ktau,nsib,ivegt(idjd),isoilm(idjd),land(idjd)
      write(6,*) 'idjd,id,jd,slwa,sgsave ',idjd,id,jd,slwa(idjd),sgsave(idjd)
      write(6,*) 'snowd,sicedep,condx ',snowd(idjd),sicedep(idjd),condx(idjd)
      write(6,*) 't1,tss ',t(idjd,1),tss(idjd)
      write(6,*) 'wb ',(wb(idjd,k),k=1,ms)
      write(6,*) 'tgg ',(tgg(idjd,k),k=1,ms)
    endif
  end if
  call maxmin(t,' t',ktau,1.,kl)
endif

!     using av_vmod (1. for no time averaging)
!      *****  check next comment
!       sflux called at beginning of time loop, hence savu, savv

tv(:) = t(1:ifull,1)*(1.+0.61*qg(1:ifull,1)-qlg(1:ifull,1)-qfg(1:ifull,1) &
                     -qrg(1:ifull,1)-qsng(1:ifull,1)-qgrg(1:ifull,1))
azmin(:) = (bet(1)*tv(:)+phi_nh(:,1))/grav
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
call START_LOG(sfluxurban_begin)                                                                 ! urban
if (nmaxpr==1) then                                                                              ! urban
  if (myid==0) then                                                                              ! urban
    write(6,*) "Before urban"                                                                    ! urban
  end if                                                                                         ! urban
  call ccmpi_barrier(comm_world)                                                                 ! urban
end if                                                                                           ! urban
if (nurban/=0) then                                                                              ! urban
  ! calculate zonal and meridonal winds                                                          ! urban
  zonx=real(                       -sin(rlat0*pi/180.)*y(:))                                     ! urban
  zony=real(sin(rlat0*pi/180.)*x(:)+cos(rlat0*pi/180.)*z(:))                                     ! urban
  zonz=real(-cos(rlat0*pi/180.)*y(:)                       )                                     ! urban
  costh= (zonx*ax(1:ifull)+zony*ay(1:ifull)+zonz*az(1:ifull)) &                                  ! urban
        /sqrt( max(zonx**2+zony**2+zonz**2,1.e-7) )                                              ! urban
  sinth=-(zonx*bx(1:ifull)+zony*by(1:ifull)+zonz*bz(1:ifull)) &                                  ! urban
        /sqrt( max(zonx**2+zony**2+zonz**2,1.e-7) )                                              ! urban
  uzon= costh*uav-sinth*vav ! zonal wind                                                         ! urban
  vmer= sinth*uav+costh*vav ! meridonal wind                                                     ! urban
  newrunoff=runoff-oldrunoff ! new runoff since entering sflux                                   ! urban
  ! since ateb will blend non-urban and urban runoff, it is                                      ! urban
  ! easier to remove the new runoff and add it again after the                                   ! urban
  ! urban scheme has been updated                                                                ! urban
  ! call aTEB                                                                                    ! urban
  dumsg=sgsave/(1.-swrsave*albvisnir(:,1)-(1.-swrsave)*albvisnir(:,2))                           ! urban
  dumrg=-rgsave                                                                                  ! urban
  dumx=condx/dt                                                                                  ! urban
  dums=(conds+condg)/dt                                                                          ! urban
  call atebcalc(fg,eg,tss,wetfac,newrunoff,dt,azmin,dumsg,dumrg,dumx,dums,rho,t(1:ifull,1), &    ! urban
                qg(1:ifull,1),ps(1:ifull),uzon,vmer,vmodmin,0)                                   ! urban
  runoff=oldrunoff+newrunoff ! add new runoff after including urban                              ! urban
  ! here we blend zo with the urban part                                                         ! urban
  call atebzo(zo,zoh,zoq,0)                                                                      ! urban
  factch=sqrt(zo/zoh)                                                                            ! urban
  ! calculate ustar                                                                              ! urban
  cduv=cduv/vmag                                                                                 ! urban
  cdtq=cdtq/vmag                                                                                 ! urban
  call atebcd(cduv,cdtq,0)                                                                       ! urban
  cduv=cduv*vmag                                                                                 ! urban
  cdtq=cdtq*vmag                                                                                 ! urban
  ustar=sqrt(vmod*cduv)                                                                          ! urban
  ! calculate screen level diagnostics                                                           ! urban
  !call atebscrnout(tscrn,qgscrn,uscrn,u10,0)                                                    ! urban
  where ( land(1:ifull) )                                                                        ! urban
    !~ qsttg(1:ifull) = qsat(ps(1:ifull),tss(1:ifull))                                              ! urban
    rnet(1:ifull) = sgsave(1:ifull) - rgsave(1:ifull) - stefbo*tss(1:ifull)**4                   ! urban
    taux(1:ifull) = rho(1:ifull)*cduv(1:ifull)*u(1:ifull,1)                                      ! urban
    tauy(1:ifull) = rho(1:ifull)*cduv(1:ifull)*v(1:ifull,1)                                      ! urban
  end where                                                                                      ! urban
end if                                                                                           ! urban
if (nmaxpr==1) then                                                                              ! urban
  if (myid==0) then                                                                              ! urban
    write(6,*) "After urban"                                                                     ! urban
  end if                                                                                         ! urban
  call ccmpi_barrier(comm_world)                                                                 ! urban
end if                                                                                           ! urban
call END_LOG(sfluxurban_end)                                                                     ! urban
! ----------------------------------------------------------------------
      
! scrnout is the standard CCAM screen level diagnostics.
! autoscrn contains the newer diagnostic calculation
!~ if (nmlo==0.and.(nsib==3.or.nsib==5).and.rescrn==0) then
  !~ call scrnout(zo,ustar,factch,wetfac,qsttg,qgscrn,tscrn,uscrn,u10,rhscrn,af,aft,ri,vmod,bprm,cms,chs,chnsea,nalpha)
!~ else
  !~ call autoscrn
!~ end if

! ----------------------------------------------------------------------
evap(:)=evap(:)+dt*eg(:)/hl !time integ value in mm (wrong for snow)

! Update runoff for river routing
if ( abs(nmlo)>=2 .or. nriver==1 ) then
  newrunoff=runoff-oldrunoff
  watbdy(1:ifull)=watbdy(1:ifull)+newrunoff ! runoff in mm
end if

!***  end of surface updating loop
! ----------------------------------------------------------------------

if(diag.or.ntest==1)then
  WRITE(6,*) 'LLL'
  if ( mydiag ) then
    WRITE(6,*) 'MMM'
    write(6,*) 'at end of sflux'
    write(6,*) 'eg,fg ',eg(idjd),fg(idjd)
    write(6,*) 'tscrn,cduv,zolnd ',tscrn(idjd),cduv(idjd),zolnd(idjd)
    write(6,*) 'snowd,sicedep ',snowd(idjd),sicedep(idjd)
    write(6,*) 'u1,v1,qg1 ',u(idjd,1),v(idjd,1),qg(idjd,1)
    write(6,*) 'w,w2,condx ',wb(idjd,1),wb(idjd,ms),condx(idjd)
    write(6,*) 't1,tss,tgg_2,tgg_ms ',t(idjd,1),tss(idjd),tgg(idjd,2),tgg(idjd,ms)
  end if
  call maxmin(tscrn,'tc',ktau,1.,1)
endif
if(ntest==4.and.ktau==10)then
  do iq=1,ifull
    if(.not.land(iq))then
      write(45,'(2g13.4)') sqrt(u(iq,1)**2+v(iq,1)**2),fg(iq)
      write(46,'(2g13.4)') sqrt(u(iq,1)**2+v(iq,1)**2),eg(iq)
    endif
    write(47,'(2g13.4)') sqrt(u(iq,1)**2+v(iq,1)**2),eg(iq)
  enddo
endif

return
end subroutine sflux

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine sib3(nalpha,taftfh,taftfhg,aft,rho)

! This is the standard land-surface scheme
      
use arrays_m                     ! Atmosphere dyamics prognostic arrays
use cc_mpi                       ! CC MPI routines
use extraout_m                   ! Additional diagnostics
use latlong_m                    ! Lat/lon coordinates
use liqwpar_m                    ! Cloud water mixing ratios
use morepbl_m                    ! Additional boundary layer diagnostics
use nsibd_m                      ! Land-surface arrays
use pbl_m                        ! Boundary layer arrays
use permsurf_m                   ! Fixed surface arrays
use screen_m                     ! Screen level diagnostics
use sigs_m                       ! Atmosphere sigma levels
use soil_m                       ! Soil and surface data
use soilsnow_m                   ! Soil, snow and surface data
use vegpar_m                     ! Vegetation arrays
use work2_m                      ! Diagnostic arrays
use work3_m                      ! Mk3 land-surface diagnostic arrays

implicit none
      
include 'newmpar.h'              ! Grid parameters
include 'const_phys.h'           ! Physical constants
include 'dates.h'                ! Date data
include 'parm.h'                 ! Model configuration
include 'soilv.h'                ! Soil parameters

integer nbarewet,nsigmf
common/nsib/nbarewet,nsigmf      ! Land-surface options

integer iq,k,iveg,layer,isoil,ip,icount
integer, intent(in) :: nalpha
real xxx,tgss,esattg,tgss2,fle,frac,conw_fh,qtgnet
real qtgair,eg1,eg2,deg,egg_alph1,sstar,ff,rsi,den
real wbav,f1,f2,f3,f4,esatf,qsatgf,beta,etr,betetrdt
real prz,dirad1,devf,ewwwa,delta_t0
real delta_t,deltat,es,tsoil
real, dimension(ifull), intent(in) :: taftfh,taftfhg,rho
real, dimension(ifull), intent(inout) :: aft
real, dimension(ifull) :: airr,cc,ccs,condxg,condsg,delta_tx,evapfb1,evapfb2,evapfb3,evapfb4
real, dimension(ifull) :: evapfb5,evapfb1a,evapfb2a,evapfb3a,evapfb4a,evapfb5a,otgf,rmcmax
real, dimension(ifull) :: tgfnew,evapfb,dqsttg,tstom,cls,omc
real, dimension(ifull) :: ftsoil,rlai,srlai,res,tsigmf,fgf,egg
real, dimension(ifull) :: evapxf,ewww,fgg,rdg,rgg,residf,fev
real, dimension(ifull) :: extin,dirad,dfgdt,degdt

integer, parameter :: ntest=0 ! ntest= 0 for diags off; ntest= 1 for diags on
!                                      2 for ewww diags      
!                     N.B. may need vsafe for correct diags
integer, parameter :: itnmeth=5
integer, parameter :: nstomata=1  ! 0 for original; 1 for tropical C4
integer, parameter :: newfgf=0    ! 0 for original; 1 with tscrn; 2 with aft too
integer, parameter :: ndiag_arr=0 ! 0 off; 1 for diag arrays on
integer, parameter :: neva=0      ! neva= 0 for diags off; neva= 1 for eva's diags on
     
do iq=1,ifull
  if(land(iq))then
    iveg=ivegt(iq)
    ! evaluate seasonal impact and the snow depth
    ! impact on the fractional vegetation cover
    tstom(iq)=298.
    if(iveg==6+31)tstom(iq)=302.
    if(iveg>=10.and.iveg<=21.and.abs(rlatt(iq)*180./pi)<25.)tstom(iq)=302.
    tsoil=min(tstom(iq), .5*(.3333*tgg(iq,2)+.6667*tgg(iq,3)+.95*tgg(iq,4) + .05*tgg(iq,5)))
    ftsoil(iq)=max(0.,1.-.0016*(tstom(iq)-tsoil)**2)
    ! which is same as:  ftsoil=max(0.,1.-.0016*(tstom-tsoil)**2)
    ! if( tsoil >= tstom ) ftsoil=1.
  endif ! (land)
enddo

if(nsib==3)then
  do iq=1,ifull
    if(land(iq))then
      iveg=ivegt(iq)
      rlai(iq)=max(.1,rlaim44(iveg)-slveg44(iveg)*(1.-ftsoil(iq)))
      srlai(iq)=rlai(iq)+rlais44(iveg)    ! nsib=3  leaf area index
      rsmin(iq) = rsunc44(iveg)/rlai(iq)  ! nsib=3  
      tsigmf(iq)=max(.001, sigmf(iq)-scveg44(iveg)*(1.-ftsoil(iq)))
      vlai(iq)=rlai(iq)
    else
      vlai(iq)=0.
    endif ! (land)
  enddo
else     ! i.e. nsib=5
  where (land)
    rlai=max(.1,vlai)
    srlai=rlai                  ! nsib=5 leaf area index
    tsigmf=max(.001,sigmf)
  end where
endif  !(nsib==3) .. else ..

if(ktau==1)then
  if(mydiag)write(6,*) 'ipland,ipsice,ipsea in sflux: ',ipland,ipsice,ipsea
  do iq=1,ifull     ! gives default over sea too
    tgf(iq)=t(iq,1)  ! was tss(iq)
    cansto(iq)=0.
  enddo 
  do iq=1,ifull
    if(land(iq))then  ! following gives agreement on restarts
      tscrn(iq)=theta(iq)  ! first guess, needed for newfgf=1
      if(nrungcm==3)then
        do layer=2,ms
          wb(iq,layer)=wb(iq,1)   ! w, w2 and wb all same initially 
        enddo
      endif  ! (nrungcm==3)
    endif  ! (land)
  enddo
  if ( mydiag ) then
    if ( land(idjd) ) then ! MJT bugfix
      iveg=ivegt(idjd)
      isoil = isoilm(idjd)
      tsoil=0.5*(0.3333*tgg(idjd,2)+0.6667*tgg(idjd,3)+0.95*tgg(idjd,4) +  0.05*tgg(idjd,5))
      write(6,*) 'nsib,iveg,isoil,nalpha,newfgf,nsigmf,tsigmf ',nsib,iveg,isoil,nalpha,newfgf,nsigmf,tsigmf(idjd)
      write(6,*) 'ftsoil,scveg44,sigmf ',ftsoil(idjd),scveg44(iveg),sigmf(idjd)
      write(6,*) 'swilt,sfc,wb1-6 ',swilt(isoil),sfc(isoil),(wb(idjd,k),k=1,ms)
      write(6,*) 'srlai,rsmin ',srlai(idjd),rsmin(idjd)
    endif
  endif
endif           ! (ktau==1)

if(ntest==1.and.mydiag) then
  iq=idjd
  iveg=ivegt(iq)
  write(6,*) 'in sib3a iq,iveg ',iq,iveg
  write(6,*) 'snowd,zo,zolnd,tstom ',snowd(iq),zo(iq),zolnd(iq),tstom(iq)
  write(6,*) 'in sib3b iq,idjd,iveg ',iq,idjd,iveg
  write(6,*) 'iveg,sigmf(iq),tsigmfa ',iveg,sigmf(iq),tsigmf(iq)
  tsoil=0.5*(0.3333*tgg(idjd,2)+0.6667*tgg(idjd,3)+0.95*tgg(idjd,4) +  0.05*tgg(idjd,5))
  write(6,*) 'rlaim44,tsoil,ftsoil ',rlaim44(iveg),tsoil,ftsoil(iq)
  write(6,*) 'scveg44,snowd,zo,zolnd,tstom ',scveg44(iveg),snowd(iq),zo(iq),zolnd(iq),tstom(iq)
  write(6,*) 'w2,rlai ',wb(iq,ms),rlai(iq)
endif ! ntest
 
do ip=1,ipland  
  iq=iperm(ip)
  tsigmf(iq)=(1.-snowd(iq)/(snowd(iq)+5.*100.*zo(iq)))*tsigmf(iq)
  ! extin(iq)=exp(-0.6*max(1.,rlai(iq)))  ! good approx uses next 2 (jlm)
  xxx=.6*max(1.,rlai(iq))
  extin(iq)=1.-xxx/(1. +.5*xxx +xxx*xxx/12.) 
  if(ntest==1.and.iq==idjd.and.mydiag) then
    write(6,*) 'in sib3c ip,iq,idjd,iveg ',ip,iq,idjd,ivegt(iq)
    write(6,*) 'iveg,sigmf(iq),tsigmf ',ivegt(iq),sigmf(iq),tsigmf(iq)
    write(6,*) 'scveg44,snowd,zo,zolnd,tstom ',scveg44(iveg),snowd(iq),zo(iq),zolnd(iq),tstom(iq)
    write(6,*) 'alb,sgsave ',albvisnir(iq,1),sgsave(iq)
    write(6,*) 'w2,rlai,extin ',wb(iq,ms),rlai(iq),extin(iq)
  endif ! ntest
  ! bare ground calculation
  tgss=isflag(iq)*tggsn(iq,1) + (1-isflag(iq))*tgg(iq,1)
  !~ esattg=establ(tgss)
  qsttg(iq)=.622*esattg/(ps(iq)-esattg)
  tgss2=tgss*tgss
  dqsttg(iq)=qsttg(iq)*ps(iq)*hlars/((ps(iq)-esattg)*tgss2)
  rgg(iq) =  stefbo*tgss2**2   ! i.e. stefbo*tgss**4
  dirad(iq)=4.*rgg(iq)/tgss
  ! sensible heat flux
  dfgdt(iq)=taftfhg(iq)*rho(iq)*cp
  fgg(iq)=dfgdt(iq)*(tgss-theta(iq))
enddo         ! ip=1,ipland
WRITE(6,*) 'NBAREWET is ', nbarewet 
select case(nbarewet)
  case(0) ! original Eva's, same as NCAR
    do ip=1,ipland  ! all land points in this nsib=3 loop
      iq=iperm(ip)
      isoil = isoilm(iq)
      fle=(wb(iq,1)-swilt(isoil))/(sfc(isoil)-swilt(isoil))          
      wetfac(iq)=max( 0.,min(1.,fle) )
    enddo   ! ip=1,ipland

  case(1) ! simplest bare
    do ip=1,ipland  ! all land points in this nsib=3 loop
      iq=iperm(ip)
      isoil = isoilm(iq)
      fle= wb(iq,1)/ssat(isoil)                                   
      wetfac(iq)=max( 0.,min(1.,fle) )
    enddo   ! ip=1,ipland

  case(2) ! jlm suggestion
    do ip=1,ipland  ! all land points in this nsib=3 loop
      iq=iperm(ip)
      isoil = isoilm(iq)
      frac=max(.01,tsigmf(iq))     ! jlm special for fle
      fle= (wb(iq,1)-frac*max( wb(iq,1),swilt(isoil) ))/(ssat(isoil)*(1.-frac))                         
      wetfac(iq)=max( 0.,min(1.,fle) )
    enddo   ! ip=1,ipland

  case(3) ! jlm suggestion
    do ip=1,ipland  ! all land points in this nsib=3 loop
      iq=iperm(ip)
      isoil = isoilm(iq)
      frac=max(.01,tsigmf(iq))     ! jlm special for fle
      fle= (wb(iq,1)-frac*swilt(isoil) )/(sfc(isoil)-frac*swilt(isoil))                    
      wetfac(iq)=max( 0.,min(1.,fle) )
    enddo   ! ip=1,ipland

  case(4)  ! jlm, similar to Noilhan & Planton
    do ip=1,ipland  ! all land points in this nsib=3 loop
      iq=iperm(ip)
      isoil = isoilm(iq)
      fle=min( 1.,wb(iq,1)/sfc(isoil) )         
      wetfac(iq)=fle*fle*(3.-2.*fle)
    enddo   ! ip=1,ipland

  case(5)  ! jlm, similar to Noilhan & Planton with swilt
    do ip=1,ipland  ! all land points in this nsib=3 loop
      iq=iperm(ip)
      isoil = isoilm(iq)
      fle=max( 0.,min( 1.,(wb(iq,1)-swilt(isoil))/(sfc(isoil)-swilt(isoil)) ) )
      wetfac(iq)=fle*fle*(3.-2.*fle)
    enddo   ! ip=1,ipland

  case(6) ! newer jlm
    do ip=1,ipland  ! all land points in this nsib=3 loop
      iq=iperm(ip)
      isoil = isoilm(iq)
      fle= max( 0.,min(1.,wb(iq,1)/ssat(isoil)) )
      wetfac(iq)=fle*fle*(2.2-1.2*fle)  ! .4 for fle=.5
    enddo   ! ip=1,ipland

  case(7) ! newest piecewise jlm (use with nalpha=1, beta)
    do ip=1,ipland  ! all land points in this nsib=3 loop
      iq=iperm(ip)
      isoil = isoilm(iq)
      wetfac(iq)=.2*(wb(iq,1)/swilt(isoil))**2
      if(wb(iq,1)>swilt(isoil))then
        wetfac(iq)=(.2*(sfc(isoil)-wb(iq,1))+.8*(wb(iq,1)-swilt(isoil)))/(sfc(isoil)-swilt(isoil))
      endif
      if(wb(iq,1)>sfc(isoil))then
        wetfac(iq)=(.8*(ssat(isoil)-wb(iq,1))+(wb(iq,1)-sfc(isoil)))/(ssat(isoil)-sfc(isoil))
      endif
    enddo   ! ip=1,ipland

  case(8) ! like NCAR but uses ssat
    do ip=1,ipland  ! all land points in this nsib=3 loop
      iq=iperm(ip)
      isoil = isoilm(iq)
      fle=(wb(iq,1)-swilt(isoil))/(ssat(isoil)-swilt(isoil))        
      wetfac(iq)=max( 0.,min(1.,fle) )
    enddo   ! ip=1,ipland

end select
      
if ( nalpha==1 ) then    ! beta scheme
  do ip = 1,ipland     ! all land points in this nsib=3 loop
    iq = iperm(ip)
    isoil = isoilm(iq)
    conw_fh = rho(iq)*taftfhg(iq)*hl    ! taftfhg(iq)=aftlandg*fhbg
    epot(iq) = conw_fh*(qsttg(iq)-qg(iq,1))
    egg(iq) = wetfac(iq)*epot(iq)
    degdt(iq) = wetfac(iq)*conw_fh*dqsttg(iq)
  end do   ! ip=1,ipland
else
  ! following is alpha scheme
  do ip = 1,ipland  ! all land points in this nsib=3 loop
    iq = iperm(ip)
    isoil = isoilm(iq)
    conw_fh = rho(iq)*taftfhg(iq)*hl    ! taftfhg(iq)=aftlandg*fhbg
    qtgnet = qsttg(iq)*wetfac(iq) - qg(iq,1)
    qtgair = qsttg(iq)*wetfac(iq) - max(qtgnet, .1*qtgnet)
    eg2 = -conw_fh*qtgair
    eg1 =  conw_fh*qsttg(iq)
    ! evaporation from the bare ground
    egg(iq) = eg1*wetfac(iq) +eg2
    epot(iq) = conw_fh*(qsttg(iq)-qg(iq,1))
    deg = wetfac(iq)*conw_fh*dqsttg(iq)
    ! following reduces degdt by factor of 10 for dew
    degdt(iq) = .55*deg + sign(.45*deg, qtgnet)
  end do   ! ip=1,ipland
end if    ! (nalpha==1) .. else ..
if((ntest==1.or.diag).and.mydiag ) then
  if (land(idjd))then ! MJT bugfix
    iq=idjd
    write(6,*) 'epot,egg,tgg1,snowd ',epot(iq),egg(iq),tgg(iq,1),snowd(iq)
    isoil = isoilm(iq)
    conw_fh=rho(iq)*taftfhg(iq)*hl    ! taftfhg(iq)=aftlandg*fhbg
    qtgnet=  qsttg(iq)*wetfac(iq) -qg(iq,1)
    qtgair=  qsttg(iq)*wetfac(iq)-max(qtgnet,.1*qtgnet)
    eg2=   -conw_fh*qtgair
    eg1=    conw_fh*qsttg(iq)
    ! evaporation from the bare ground
    egg_alph1=wetfac(iq)*epot(iq)
    write(6,*) 'then iq,isoil,conw_fh,qsttg,qtgair ',iq,isoil,conw_fh,qsttg(iq),qtgair
    write(6,*) 'eg1,eg2,wetfac ',eg1,eg2,wetfac(iq)
    write(6,*) 'epot,egg,egg_alph1 ',epot(iq),egg(iq),egg_alph1
  endif
endif  ! (ntest==1)

do ip=1,ipland  ! all land points in this nsib=3 loop
  iq=iperm(ip)
  if(snowd(iq)>1.)then
    egg(iq)=epot(iq)
    wetfac(iq)=1.   ! added jlm 18/3/04 to fix qgscrn inconsistency
    cls(iq)=1.+hlf/hl
  else
    egg(iq)=min(egg(iq),wb(iq,1)*zse(1)*1000.*hl/dt)
    cls(iq)=1.
  endif  ! (snowd(iq)>1.)
enddo   ! ip=1,ipland

if(nsigmf==0)then  ! original
  do ip=1,ipland  ! all land points in this nsib=3 loop
    iq=iperm(ip)
    ga(iq)=-slwa(iq)-rgg(iq)-fgg(iq)-cls(iq)*egg(iq)       
    dgdtg(iq)=-dirad(iq)-dfgdt(iq)-cls(iq)*degdt(iq)
  enddo   ! ip=1,ipland
endif     ! (nsigmf==0)

if(nsigmf==1)then  ! jlm preferred
  ! spreads bare-soil flux across whole grid square      
  do ip=1,ipland  ! all land points in this nsib=3 loop
    iq=iperm(ip)
    ga(iq)=(-slwa(iq)-rgg(iq)-fgg(iq)-cls(iq)*egg(iq))*(1.-tsigmf(iq))
    ! dgtdg is used in soilsnow
    dgdtg(iq)=-(dirad(iq)+dfgdt(iq)+cls(iq)*degdt(iq))*(1.-tsigmf(iq))
  enddo   ! ip=1,ipland
endif     ! (nsigmf==1)

if(ntest==1.and.mydiag)then
  iq=idjd
  write(6,*)'dgdtg,dirad,dfgdt,cls,degdt,tsigmf',dgdtg(iq),dirad(iq),dfgdt(iq),cls(iq),degdt(iq),tsigmf(iq) 
endif

! ----------------------------------------------
do ip=1,ipland  ! all land points in this nsib=3 loop
  iq=iperm(ip)
  isoil = isoilm(iq)
  iveg=ivegt(iq)
  ! components of the stomatal resistance
  sstar=90.+sign(60.,.5-zo(iq))  ! equiv to above 2 lines
  ff= 1.1*sgsave(iq)/(rlai(iq)*sstar)
  rsi = rsmin(iq) * rlai(iq)
  f1= (1.+ff)/(ff+rsi/5000.)
  den=sfc(isoil)-swilt(isoil)                          ! sib3
  wbav=(max(0.,froot(1)*(wb(iq,1)-swilt(isoil)))+  &
        max(0.,froot(2)*(wb(iq,2)-swilt(isoil)))+  &
        max(0.,froot(3)*(wb(iq,3)-swilt(isoil)))+  &
        max(0.,froot(4)*(wb(iq,4)-swilt(isoil)))+  &
        max(0.,froot(5)*(wb(iq,5)-swilt(isoil)))   )/den
  f2=max(1. , .5/ max( wbav,1.e-7)) ! N.B. this is equiv to next 2 (jlm)
  f4=max(1.-.0016*(tstom(iq)-t(iq,1))**2 , .05) ! zero for delta_t=25
  airr(iq) = 1./taftfh(iq)
  cc(iq) =min(condx(iq),           4./(1440. *60./dt))  ! jlm speedup for 4 mm/day
  ccs(iq)=min(conds(iq)+condg(iq), 4./(1440. *60./dt))
  ! depth of the reservoir of water on the canopy
  rmcmax(iq) = max(0.5,srlai(iq)) * .1
  omc(iq) = cansto(iq)  ! value from previous timestep as starting guess
  !~ f3=max(1.-.00025*(establ(t(iq,1))-qg(iq,1)*ps(iq)/.622),.05)
  res(iq)=max(30.,rsmin(iq)*f1*f2/(f3*f4))
  if(ntest==1.and.iq==idjd.and.mydiag)then
    write(6,*) 'rlai,srlai,wbav,den ',rlai(iq),srlai(iq),wbav,den
    write(6,*) 'f1,f2,f3,f4 ',f1,f2,f3,f4
    write(6,*) 'ff,f124,rsi,res ',ff,f1*f2/f4,rsi,res(iq)
    write(6,*) 'qg,qfg,qlg ',qg(iq,1),qfg(iq,1),qlg(iq,1)
  endif
  otgf(iq)=tgf(iq)
  tgfnew(iq)=tgf(iq)
  delta_tx(iq)=10.   ! just to supply max change of 5 deg first time
  evapfb1a(iq)=max(0.,wb(iq,1)-swilt(isoil)) *zse(1)*1000.
  evapfb2a(iq)=max(0.,wb(iq,2)-swilt(isoil)) *zse(2)*1000.
  evapfb3a(iq)=max(0.,wb(iq,3)-swilt(isoil)) *zse(3)*1000.
  evapfb4a(iq)=max(0.,wb(iq,4)-swilt(isoil)) *zse(4)*1000.
  evapfb5a(iq)=max(0.,wb(iq,5)-swilt(isoil)) *zse(5)*1000.
enddo   ! ip loop

do icount=1,itnmeth     ! jlm new iteration
  ! transpiration
  do ip=1,ipland  ! all land points in this nsib=3 loop
    iq=iperm(ip)
    !~ esatf = establ(tgfnew(iq))
    qsatgf=.622*esatf/(ps(iq)-esatf)
    ! wet evaporation
    ewwwa = rho(iq) *(qsatgf-qg(iq,1))/airr(iq) ! in W/m**2 /hl
    ! max available dewfall is 
    ! -(qsatgf-qg1)*dsig*1000*ps/grav in mm (mult by hl/dt for W/m**2)
    ewwwa=max(ewwwa,-abs((qsatgf-qg(iq,1))*dsig(1)*ps(iq))/(grav*dt))
    ewww(iq)=min(dt*ewwwa,omc(iq),dt*ewwwa*omc(iq)/rmcmax(iq))
    ! dew(-ve), no_dew ,        no_dew
    ! cansto is reservoir on leaf
    cansto(iq)=(omc(iq)-ewww(iq)) +cc(iq)
    ewww(iq)=ewww(iq)/dt  ! these changes on 19/1/06 jlm

    ! precipitation reaching the ground under the canopy
    ! water interception by the canopy
    condxg(iq)=max(condx(iq)          -cc(iq) +max(0.,cansto(iq)-rmcmax(iq)),0.) ! keep
    condsg(iq)=max(conds(iq)+condg(iq)-ccs(iq)+max(0.,cansto(iq)-rmcmax(iq)),0.)
    cansto(iq) = min( max(0.,cansto(iq)), rmcmax(iq))
    beta =      cansto(iq)/rmcmax(iq)
    Etr=rho(iq)*max(0.,qsatgf-qg(iq,1))/(airr(iq) +res(iq))  ! jlm
    betetrdt =(1.-beta)*Etr*dt*tsigmf(iq)   ! fixed 23/5/01
    evapfb1(iq)=min(betetrdt*froot(1),evapfb1a(iq))
    evapfb2(iq)=min(betetrdt*froot(2),evapfb2a(iq))
    evapfb3(iq)=min(betetrdt*froot(3),evapfb3a(iq))
    evapfb4(iq)=min(betetrdt*froot(4),evapfb4a(iq))
    evapfb5(iq)=min(betetrdt*froot(5),evapfb5a(iq))
    evapfb(iq)=(evapfb1(iq)+evapfb2(iq)+evapfb3(iq)+evapfb4(iq)+evapfb5(iq))/tsigmf(iq)
    evapxf(iq) = (evapfb(iq)/dt + ewww(iq))*hl  ! converting to W/m**2
    prz = rho(iq)*cp*taftfh(iq)
    if(newfgf==0)fgf(iq) = prz*(tgfnew(iq)-theta(iq))  ! original/usual
    if(newfgf==1)fgf(iq) = prz*(tgfnew(iq)-tscrn(iq))
    if(newfgf==2)fgf(iq)=rho(iq)*aft(iq)*cp*(tgfnew(iq)-tscrn(iq))
    ! limit extreme fgf to avoid undue tgf oscillations  June '04
    fgf(iq)=max(-1000.,min(fgf(iq),1000.))
    rdg(iq) =  stefbo*tgfnew(iq)**4
    residf(iq) = -slwa(iq) - rdg(iq) - fgf(iq) - evapxf(iq)
    dirad1 = 4.*rdg(iq)/300.
    ! next 2 expressions can be approximated without effects
    ! dqg=qsatgf*hlars/300.**2
    ! devf= hl*rho(iq)*dqg*( (1.-beta)/(airr(iq) + res(iq))+beta/airr(iq) ) ! re-factored by jlm
    ! according to jlm prints, devf has only small effect
    devf= (hl*hlars/300.**2)*qsatgf*(1.-beta)/res(iq)
    delta_t0=residf(iq)/(dirad1 + devf + prz)
    delta_t=sign(min(abs(delta_t0),.5*abs(delta_tx(iq))),delta_t0)
    tgfnew(iq)=tgfnew(iq)+delta_t
    delta_tx(iq)=tgfnew(iq)-otgf(iq)
  enddo  ! ip loop
  if((ntest==1.or.diag).and.mydiag)then 
    if(land(idjd))then
      iq=idjd
      write(6,*) 'ktau,icount,iq,omc,cc ',ktau,icount,iq,omc(iq),cc(iq)
      write(6,*) 'rmc,rmcmax,ewww ',cansto(iq),rmcmax(iq),ewww(iq)
      !~ esatf = establ(tgfnew(iq))  ! value for next itn
      qsatgf=.622*esatf/(ps(iq)-esatf)
      ewwwa = rho(iq) *(qsatgf-qg(iq,1))/airr(iq)
      write(6,*) 'esatf,qsatgf,ewwwa ',esatf,qsatgf,ewwwa
      prz = rho(iq)*cp*taftfh(iq)
      beta =      cansto(iq)/rmcmax(iq)
      devf= (hl*hlars/300.**2)*qsatgf*(1.-beta)/res(iq)
      dirad1 = 4.*rdg(iq)/300.
      write(6,*) 'beta,airr,res ',beta,airr(iq),res(iq)
      write(6,*) 'dirad1,devf,prz ',dirad1,devf,prz
      write(6,*) 'theta,tscrn,slwa ',theta(iq),tscrn(iq),slwa(iq)
      write(6,*) 'taftfh,condxg ',taftfh(iq),condxg(iq)
      write(6,*) 'rdg,fgf,evapxf,evapfb ',rdg(iq),fgf(iq),evapxf(iq),evapfb(iq)
      write(6,*) 'delta_tx ',delta_tx(iq)
      write(6,*) 'otgf,tgfnew,residf ',otgf(iq),tgfnew(iq),residf(iq)
    endif  ! (land(idjd))
  endif   ! ((ntest==2.or.diag).and.mydiag)
enddo     !  icount=1,5

do ip=1,ipland  ! all land points in this nsib=3 loop
  iq=iperm(ip)
  if(cansto(iq)<1.e-10)cansto(iq)=0.  ! to avoid underflow 24/1/06
  if(tsigmf(iq) <= .0101) then
    condxpr(iq)=condx(iq)
    condspr(iq)=conds(iq)+condg(iq)
    evapfb(iq) = 0.
    evapxf(iq) = egg(iq)
    fgf(iq)  = fgg(iq)
    rdg(iq)=rgg(iq)
    tgf(iq) = tss(iq)
  else
    tgf(iq)=tgfnew(iq)
    wb(iq,1)=wb(iq,1)-evapfb1(iq)/(zse(1)*1000.)
    wb(iq,2)=wb(iq,2)-evapfb2(iq)/(zse(2)*1000.)
    wb(iq,3)=wb(iq,3)-evapfb3(iq)/(zse(3)*1000.)
    wb(iq,4)=wb(iq,4)-evapfb4(iq)/(zse(4)*1000.)
    wb(iq,5)=wb(iq,5)-evapfb5(iq)/(zse(5)*1000.)
    condxpr(iq)=(1.-tsigmf(iq))*condx(iq)+tsigmf(iq)*condxg(iq)
    condspr(iq)=(1.-tsigmf(iq))*(conds(iq)+condg(iq))+tsigmf(iq)*condsg(iq)
    if(ntest==1.and.abs(residf(iq))>10.) then
      write(6,*) 'iq,otgf(iq),tgf,delta_tx,residf ',iq,otgf(iq),tgf(iq),delta_tx(iq),residf(iq)
    end if
  endif          ! tsigmf <= .01   ... else ...
  fev(iq)=evapfb(iq)/dt*hl*tsigmf(iq) ! passed to soilsnow to update wb
  fes(iq)=(1.-tsigmf(iq))*egg(iq)*cls(iq)  ! also passed to soilsnow
  otgsoil(iq)=isflag(iq)*tggsn(iq,1) + (1-isflag(iq))*tgg(iq,1)
enddo  !  ip=1,ipland
if((ntest==1.or.diag).and.mydiag) then
  if(land(idjd))then ! MJT bugfix
    iq=idjd
    isoil = isoilm(iq)
    iveg=ivegt(iq)
    write(6,*) 'in sib3 before soilsnowv'
    write(6,*) 'evapxf,epot,egg,fev,wetfac ',evapxf(iq),epot(iq),egg(iq),fev(iq),wetfac(iq)
    write(6,*) 'fgf,fgg,fes ',fgf(iq),fgg(iq),fes(iq)
    write(6,*) 'isoil,ssat,tsigmf,rlai ',isoil,ssat(isoil),tsigmf(iq),rlai(iq)
    write(6,*) 'tgg1,t1,theta,tscrn ',tgg(iq,1),t(iq,1),theta(iq),tscrn(iq)
    write(6,*) 'qg1,qsttg ',qg(iq,1),qsttg(iq)
    write(6,*) 'dfgdt,taftfhg,rho ',dfgdt(iq),taftfhg(iq),rho(iq)
    write(6,*) 'rmc,rmcmax(iq) ',cansto(iq),rmcmax(iq)
    if (abs(tgf(iq)-otgf(iq))>4.9) then
      write(6,"('ktau,iq,otgf,tgf,dtgf,t1,t2',i4,i6,5f8.2)") ktau,iq,otgf(iq),tgf(iq),tgf(iq)-otgf(iq),t(iq,1),t(iq,2)
    end if
  endif
endif
!-------------------------------------

call soilsnowv

do ip=1,ipland  ! all land points in this nsib=3 loop
  iq=iperm(ip)
  if(isflag(iq)==0) then
    deltat=tgg(iq,1)-otgsoil(iq)
    fgg(iq)=fgg(iq)+deltat*dfgdt(iq)
    egg(iq)=egg(iq)+deltat*degdt(iq)
    egg(iq)=min(egg(iq),wb(iq,1)*zse(1)*1000.*hl/dt)
    rgg(iq)=rgg(iq)+deltat*dirad(iq)
  else
    deltat=tggsn(iq,1)-otgsoil(iq)
    fgg(iq)=fgg(iq)+deltat*dfgdt(iq)
    egg(iq)=egg(iq)+deltat*degdt(iq)
    rgg(iq)=rgg(iq)+deltat*dirad(iq)
  endif
  ! combined fluxes
  if(snowd(iq)>1.)then
    eg(iq)=tsigmf(iq)*evapxf(iq) + egg(iq)
  else
    eg(iq) = tsigmf(iq)*evapxf(iq) + (1. - tsigmf(iq))*egg(iq)
  endif
  if(nsigmf==2)then
    fg(iq)=tsigmf(iq)*fgf(iq)+fgg(iq)
  else
    fg(iq)=tsigmf(iq)*fgf(iq)+(1.-tsigmf(iq))*fgg(iq)
  endif
  rnet(iq)=-slwa(iq)-(1.-tsigmf(iq))*rgg(iq)-tsigmf(iq)*rdg(iq)

  tgss=isflag(iq)*tggsn(iq,1) + (1-isflag(iq))*tgg(iq,1)  ! jlm
  if(tsigmf(iq)<= .01) then
    tss(iq) = tgss
    tgf(iq) = tgss
  else
    tss(iq)=tsigmf(iq)*tgf(iq)+(1.-tsigmf(iq))*tgss
  endif       ! tsigmf<= .01
  !~ es = establ(tss(iq))     !  from 27/12/05
  qsttg(iq)= .622*es/(ps(iq)-es)  ! recal for scrnout, esp. snow    

enddo   ! ip=1,ipland
      
! Calculate fraction of canopy which is wet
fwet=0.
where (land)
  fwet=cansto/rmcmax
end where

if((ntest==1.or.diag).and.mydiag) then
  if (land(idjd))then ! MJT bugfix
    iq=idjd
    write(6,*) 'even further down sib3 after soilsnowv'
    write(6,*) 'tgg ',(tgg(iq,k),k=1,ms)
    write(6,*) 'wb ',(wb(iq,k),k=1,ms)
    write(6,*) 'isflag,snowd ',isflag(iq),snowd(iq)
    write(6,*) 'evapfb,fev,ewww ',evapfb(iq),fev(iq),ewww(iq)
    write(6,*) 'tsigmf,evapxf,egg ',tsigmf(iq),evapxf(iq),egg(iq)
    write(6,*) 'deltat,degdt,wb,zse ',tgg(iq,1)-otgsoil(iq),degdt(iq),wb(iq,1),zse(1)
    write(6,*) 'eg,fg ',eg(iq),fg(iq)
  endif
endif

return
end subroutine sib3
