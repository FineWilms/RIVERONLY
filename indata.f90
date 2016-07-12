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
    
! This subroutine initialises CCAM prognostic variables and
! surface forcings.  Currently, we assume that disk io is
! slower than MPI, so surface forcings are loaded on a
! single processor and then distributed using MPI_distrubute.
! Note that nested files are split over processors (see
! onthefly.f).

module indata

private
public indataf

interface datacheck
  module procedure rdatacheck, idatacheck
end interface datacheck

contains
    
subroutine indataf(hourst,isoth,nsig,io_nest)
     
use arrays_m                                     ! Atmosphere dyamics prognostic arrays
use bigxy4_m                                     ! Grid interpolation
use cable_ccam, only : loadcbmparm,loadtile      ! CABLE interface
use cc_mpi                                       ! CC MPI routines
use diag_m                                       ! Diagnostic routines
use indices_m                                    ! Grid index arrays
use infile                                       ! Input file routines
use latlong_m                                    ! Lat/lon coordinates
use map_m                                        ! Grid map arrays
use morepbl_m                                    ! Additional boundary layer diagnostics
use nsibd_m                                      ! Land-surface arrays
use onthefly_m                                   ! Input interpolation routines
use pbl_m                                        ! Boundary layer arrays
use permsurf_m                                   ! Fixed surface arrays
use river                                        ! River routing
use sigs_m                                       ! Atmosphere sigma levels
use soil_m                                       ! Soil and surface data
use soilsnow_m                                   ! Soil, snow and surface data
use vecsuv_m                                     ! Map to cartesian coordinates
use xyzinfo_m                                    ! Grid coordinate arrays
      
implicit none
      
include 'newmpar.h'                              ! Grid parameters
include 'const_phys.h'                           ! Physical constants
include 'darcdf.h'                               ! Netcdf data
include 'dates.h'                                ! Date data
include 'filnames.h'                             ! Filenames
include 'parm.h'                                 ! Model configuration
include 'parmdyn.h'                              ! Dynamics parmaters
include 'parmgeom.h'                             ! Coordinate data
include 'soilv.h'                                ! Soil parameters
include 'stime.h'                                ! File date data

integer, parameter :: jlmsigmf=1  ! 1 for jlm fixes to dean's data
integer, parameter :: nfixwb=2    ! 0, 1 or 2; wb fixes with nrungcm=1
integer, parameter :: ntest=0
integer, parameter :: klmax=100   ! Maximum vertical levels

!     for the held-suarez test
real, parameter :: delty = 60.    ! pole to equator variation in equal temperature
real, parameter :: deltheta = 10. ! vertical variation
real, parameter :: rkappa = 2./7.


integer, intent(inout) :: io_nest
integer ii, imo, indexi, indexl, indexs, ip, iq, isoil, isoth
integer iveg, iyr, jj, k, kdate_sav, ktime_sav, l
integer nface, nn, nsig, i, j, n
integer ierr, ic, jc, iqg, ig, jg
integer isav, jsav, ier

character(len=160) :: surfin
character(len=80) :: header

real, intent(out) :: hourst
real, dimension(ifull) :: zss, aa, zsmask
real, dimension(ifull) :: dep, depth, rlai
real, dimension(ifull,5) :: duma
real, dimension(ifull,2) :: ocndwn
real, dimension(ifull,kl,9) :: dumb
real, dimension(:,:), allocatable :: glob2d
real, dimension(:), allocatable :: davt_g
real, dimension(3*kl+1) :: dumc
real, dimension(1:9) :: swilt_diag, sfc_diag
real, dimension(1:ms) :: wb_tmpry
real rlonx,rlatx,alf
real c, cent
real coslat, coslong, costh, den, diffb, diffg, dist
real epsmax, fracs, fracwet, ftsoil
real polenx, poleny, polenz, pslavge
real rad, radu, radv, ri, rj, rlat_d, rlon_d
real rlatd, rlongd
real sinlat, sinlong, sinth, snalb,sumdsig, thet, tsoil
real uzon, vmer, wet3, zonx, zony, zonz, zsdiff, tstom
real xbub, ybub, xc, yc, zc, xt, yt, zt, tbubb, emcent
real deli, delj, centi, distnew, distx, rhs, ril2
real newzo,visalb,niralb

! The following look-up tables are for the Mk3 land-surface scheme
real, dimension(44), parameter :: vegpmin = (/                           &
                    .98,.85,.85,.5,.2,.1 ,.85,.5,.2,.5,                  & ! 1-10
                    .2,.1 ,.5,.2,.1 ,.1,.1 ,.85,.5,.2,                   & ! 11-20
                    .1 ,.85,.60,.50,.5 ,.2,.1 ,.5, .0, .0, .4,           & ! 21-31
                    .98,.75,.75,.75,.5,.86,.65,.79,.3, .42,.02,.54,0./)    ! 32-44
real, dimension(44), parameter :: vegpmax = (/                           &
                    .98,.85,.85,.5,.7,.60,.85,.5,.5,.5,                  & ! 1-10
                    .5,.50,.5,.6,.60,.4,.40,.85,.5,.8,                   & ! 11-20
                    .20,.85,.85,.50,.80,.7,.40,.5, .0, .0, .6,           & ! 21-31
                    .98,.75,.75,.75,.5,.86,.65,.79,.3, .42,.02,.54,0./)    ! 32-44
real, dimension(12), parameter :: fracsum =                              &
                     (/-.5,-.5,-.3,-.1,.1, .3, .5, .5, .3, .1,-.1,-.3/)
real, dimension(44), parameter :: fracwets = (/                          &
                    .5, .5, .5, .5, .5, .5, .5, .5, .5, .5,              & !  1-10 summer
                    .5, .5, .5, .5, .5, .5, .5, .5, .5, .5,              & ! 11-20 summer
                    .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5,          & ! 21-31 summer
       .5,.5, .3, .3, .3, .15, .15, .15, .1, .15, .02, .35, .5/)           ! 32-44 summer
real, dimension(44), parameter :: fracwetw = (/                          &
                    .5, .5, .5, .5, .5, .5, .5, .5, .5, .5,              & !  1-10 winter
                    .5, .5, .5, .5, .5, .5, .5, .5, .5, .5,              & ! 11-20 winter
                    .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5,          & ! 21-31 winter
       .5,.5, .6, .6, .6, .25, .3 , .25, .2, .25, .05, .6, .5 /)           ! 32-44 winter
real vegpsig(44)
data vegpsig/ .98,.85,.85,.5,.2,.05,.85,.5,.2,.5,                        & ! 1-10
              .2,.05,.5,.2,.05,.2,.05,.85,.5,.2,                         & ! 11-20
              .05,.85,.85,.55,.65,.2,.05,.5, .0, .0, .5,                 & ! 21-31
              .98,.75,.75,.75,.5,.86,.65,.79,.3, .42,.02,.54,0./           ! 32-44

real, dimension(klmax) :: tbarr,qgin,zbub
real :: pmsl=1.010e5, thlapse=3.e-3, tsea=290., gauss=2.
real :: heightin=2000., hfact=0.1, uin=0., vin=0.

namelist/tin/gauss,heightin,hfact,pmsl,qgin,tbarr,tsea,uin,vin,thlapse

call START_LOG(indata_begin)

!--------------------------------------------------------------
! SET DEFAULT VALUES
tgg(:,:)=280.
tggsn(:,:)=280.
wb(:,:)=.15
snowd(:)=0. 
condx(:)=0.
zolnd(:)=zobgin
eg(:)=0.
fg(:)=0.
cduv(:)=0.
cdtq(:)=0.
hourst=0.
albvissav(:)=-1.
albvisnir(:,:)=0.3
ivegt(:)=1
isoilm(:)=1
zs(:)=0.
zsmask(:)=0.    
land(:)=.false.
kdate=kdate_s
ktime=ktime_s

!--------------------------------------------------------------
!~ ! READ AND PROCESS ATMOSPHERE SIGMA LEVELS
call ccnf_open(vegfile,ncidveg,ierr)
lncveg=1

!--------------------------------------------------------------
! READ OROGRAPHY (io_in and nhstest)
!     read in fresh zs, land-sea mask (land where +ve), variances
if ( io_in<=4 .and. nhstest>=0 ) then
  if ( myid==0 ) then
    allocate( glob2d(ifull_g,3) )
    if ( lnctopo==1 ) then
      write(6,*) 'read zs from topofile'
      call surfread(glob2d(:,1),'zs',netcdfid=ncidtopo)
      glob2d(:,1) = grav*glob2d(:,1)
      write(6,*) 'read land-sea fraction'
      call surfread(glob2d(:,2),'lsm',netcdfid=ncidtopo)
      write(6,*) 'read he'
      call surfread(glob2d(:,3),'tsd',netcdfid=ncidtopo)
      call ccnf_close(ncidtopo)
    else
      write(6,*) 'read zs from topofile'
      read(66,*,iostat=ierr) glob2d(:,1)
      if ( ierr/=0 ) stop 'end-of-file reached on topofile'
      write(6,*) 'read land-sea fraction'
      read(66,*,iostat=ierr) glob2d(:,2)
      if ( ierr/=0 ) stop 'end-of-file reached on topofile'
      write(6,*) 'read he'
      read(66,*,iostat=ierr) glob2d(:,3)
      if ( ierr/=0 ) stop 'end-of-file reached on topofile'
      close(66)
    end if
    call ccmpi_distribute(duma(:,1:3),glob2d(:,1:3))
    deallocate(glob2d)
  else
    call ccmpi_distribute(duma(:,1:3))
  end if
  zs(1:ifull)     = duma(:,1)
  zsmask(1:ifull) = duma(:,2)

  land(1:ifull)=zsmask(1:ifull)>=0.5
 
else                   ! aquaplanet test -1 to -8 or -22
  zs(:)=0.             ! or pgb from June 2003
  zsmask(:)=0.    
  land(:)=.false.
endif  ! (io_in<=4.and.nhstest>=0)  ..else..

if ( mydiag ) then
  write(6,"('zs#_topof ',9f8.1)") diagvals(zs)
  write(6,"('zs#_mask ',9f8.2)") diagvals(zsmask)
end if


!-----------------------------------------------------------------
! The following parameterisations require special input data.
! Therefore, these routines are initialised in indata.f instead of
! initialised in globpe.f


!--------------------------------------------------------------
! READ SURFACE DATA (nsib and nspecial)
! nsib=3 (original land surface scheme with original 1deg+Dean's datasets)
! nsib=5 (original land surface scheme with MODIS datasets)
! nsib=6 (CABLE land surface scheme with internal screen diagnostics)
! nsib=7 (CABLE land surface scheme with CCAM screen diagnostics)
if (nsib>=1) then
  call insoil
  call rdnsib
  if (nsib==6.or.nsib==7) then
    ! albvisnir at this point holds soil albedo for cable initialisation
    call loadcbmparm(vegfile,vegprev,vegnext)
    ! albvisnir at this point holds net albedo
  end if ! (nsib==6.or.nsib==7) ..else..
 
end if   ! nsib>=1


!**************************************************************
!**************************************************************
! No changes to land, isoilm or ivegt arrays after this point
!**************************************************************
!**************************************************************


!--------------------------------------------------------------
! LAND SURFACE ERROR CHECKING
if(nsib>=1)then   !  check here for soil & veg mismatches
  if (mydiag) write(6,*)'idjd,land,isoil,ivegt ',idjd,land(idjd),isoilm(idjd),ivegt(idjd)
  do iq=1,ifull
    if(land(iq))then
      if(ivegt(iq)==0)then
        write(6,*)'stopping because nsib>=1 and veg type not defined for iq = ',iq
        write(6,*)'lat,long ',rlatt(iq)*180./pi,rlongg(iq)*180./pi
        call ccmpi_abort(-1)
      endif  ! (ivegt(iq)==0)
      if(isoilm(iq)==0)then
        write(6,*)'stopping because nsib>=1 and soil type not defined for iq = ',iq
        call ccmpi_abort(-1)
      endif  ! (isoilm(iq)==0)
    endif    ! (land(iq))
  enddo    !  iq loop
endif      ! (nsib>=1)

!-----------------------------------------------------------------
! INTIALISE MIXED LAYER OCEAN (nmlo)
! nmlo<0 same as below, but save all data in history file
! nmlo=0 no mixed layer ocean
! nmlo=1 mixed layer ocean (KPP)
! nmlo=2 same as 1, but with Smag horz diffusion and river routing
! nmlo=3 same as 2, but with horizontal and vertical advection
if (nriver==1 ) then
  if (myid==0) write(6,*) 'Initialising river routing'
  call rvrinit
end if

!--------------------------------------------------------------
! DEFINE FIXED SURFACE ARRAYS
! Note that now only land and water points are allowed, as
! sea-ice can change during the run
indexl=0
do iq=1,ifull
  if(land(iq))then  ! land
    indexl=indexl+1
    iperm(indexl)=iq
  endif ! (land(iq))
enddo   ! iq loop
ipland=indexl
indexi=ipland
ipsea=ifull
indexs=ipsea+1
do iq=1,ifull
  if(.not.land(iq))then
    indexs=indexs-1     ! sea point
    iperm(indexs)=iq    ! sea point
  endif  ! (sicedep(iq)>0.)
enddo   ! iq loop
ipsice=indexs-1
if (mydiag) write(6,*)'ipland,ipsea: ',ipland,ipsea

!-----------------------------------------------------------------
! READ INITIAL CONDITIONS FROM IFILE (io_in)
ncid=-1  ! initialise nc handle with no files open
if ( io_in<4 ) then
  if ( myid==0 ) then
    write(6,*) 'Read initial conditions from ifile'
  end if
  call histopen(ncid,ifile,ier) ! open parallel initial condition files (onthefly will close ncid)
  call ncmsg("ifile",ier)       ! report error messages
  if (myid==0) then
    write(6,*) 'ncid,ifile ',ncid,trim(ifile)
  end if
  kdate_sav=kdate_s
  ktime_sav=ktime_s
  zss=zs(1:ifull)
  if (abs(io_in)==1) then
     call         onthefly(0,kdate,ktime,psl,zss,tss,sicedep,fracice,t,u,v,qg,tgg,wb,wbice,snowd, &
				tggsn,smass,ssdn,ssdnn,snage,isflag)
  endif   ! (abs(io_in)==1)
  call histclose
  if(mydiag)then
    write(6,*)'ds,zss',ds,zss(idjd)
    write(6,*)'kdate_sav,ktime_sav ',kdate_sav,ktime_sav
    write(6,*)'kdate_s,ktime_s >= ',kdate_s,ktime_s
    write(6,*)'kdate,ktime ',kdate,ktime
    write(6,"(' wbice(1-ms)',9f7.3)")(wbice(idjd,k),k=1,ms)
  endif
  if(kdate/=kdate_sav.or.ktime/=ktime_sav)then
    write(6,*) 'stopping in indata, not finding correct kdate/ktime'
    write(6,*) "kdate,    ktime     ",kdate,ktime
    write(6,*) "kdate_sav,ktime_sav ",kdate_sav,ktime_sav
    call ccmpi_abort(-1)
  endif


  tss(1:ifull)=abs(tss(1:ifull)) ! not done in infile because -ve needed for onthefly
  hourst=.01*ktime

  qg(1:ifull,:)=max(qg(1:ifull,:),0.)
  ps(1:ifull)=1.e5*exp(psl(1:ifull))
        
else

  ! read in namelist for uin,vin,tbarr etc. for special runs
  if (myid==0) write(6,*)'Read IC from namelist tinit'
  read (99, tin)
  if (myid==0) write(6, tin)

endif   ! (io_in<4)

!--------------------------------------------------------------
! DEFINE SURFACE DATA PRESETS (nrungcm)

!     nrungcm<0 controls presets for snowd, wb, tgg and other soil variables
!     they can be: preset/read_in_from_previous_run
!                  written_out/not_written_out    after 24 h as follows:
!          nrungcm = -1  preset           | not written to separate file
!                    -2  preset           |     written  
!                    -3  read_in          |     written  (usual for NWP)
!                    -4  read_in          | not written  (usual for netCDF input)
!                    -5  read_in (not wb) |     written  (should be good)
!                    -6 same as -1 bit tapered wb over dry interio of Aust
!                    >5 like -1 but sets most wb percentages

! preset soil data
WRITe(6,*) 'NRUNGCM = ',nrungcm
!~ if ( .not.lrestart ) then
  if((nrungcm.le.-1.and.nrungcm.ge.-2).or.nrungcm==-6.or.nrungcm>5)then
    ! presetting wb when no soil moisture available initially
    iyr=kdate/10000
    imo=(kdate-10000*iyr)/100
    do iq=1,ifull
      if(land(iq))then
        iveg=ivegt(iq)
        if (nsib==6.or.nsib==7) iveg=1
        isoil=isoilm(iq)
        rlonx=rlongg(iq)*180./pi
        rlatx=rlatt(iq)*180./pi
        ! fracsum(imo) is .5 for nh summer value, -.5 for nh winter value
        fracs=sign(1.,rlatt(iq))*fracsum(imo)  ! +ve for local summer
        if(nrungcm>5)then
          fracwet=.01*nrungcm   ! e.g. 50 gives .5
        else
          fracwet=(.5+fracs)*fracwets(iveg)+(.5-fracs)*fracwetw(iveg)
          ! N.B. for all Dean's points, fracwet=fracwets=fracwetw=.5           
        endif
        wb(iq,ms)= (1.-fracwet)*swilt(isoilm(iq))+fracwet*sfc(isoilm(iq)) 
        if(abs(rlatx)<18.)wb(iq,ms)=sfc(isoilm(iq)) ! tropics
        ! following jlm subtropics from Aug 2003 (.1/.9), (.6, .4)
        if(rlatx<20..and.rlatx>8.) then
          wb(iq,ms)=(.35-.5*fracsum(imo))*swilt(isoilm(iq))+(.65+.5*fracsum(imo))*sfc(isoilm(iq)) ! NH
        end if
        if(rlatx>-16..and.rlatx<-8.) then
          wb(iq,ms)=(.35+.5*fracsum(imo))*swilt(isoilm(iq))+(.65-.5*fracsum(imo))*sfc(isoilm(iq)) ! SH
        end if
        if(rlatx>-32..and.rlatx<-22..and.rlonx>117..and.rlonx<146.) then
          if(nrungcm==-6)then
            ! following tapers it over 2 degrees lat/long
            alf=.5*min(abs(rlonx-117.),abs(rlonx-146.),abs(rlatx+22.),abs(rlatx+32.),2.)
            wb(iq,ms)=alf*swilt(isoilm(iq))+(1.-alf)*wb(iq,ms)
          else
            wb(iq,ms)=swilt(isoilm(iq)) ! dry interior of Australia
          endif
        endif
      endif    ! (land(iq))
    enddo     ! iq loop
    do k=1,ms-1
      wb(:,k)=wb(:,ms)
    enddo    !  k loop

    do iq=1,ifull
      ! fix for antarctic snow
      if(land(iq).and.rlatt(iq)*180./pi<-60.)snowd(iq)=max(snowd(iq),400.)
    enddo   ! iq loop

    if ( mydiag ) then
      iveg=ivegt(idjd)
      if (nsib==6.or.nsib==7) iveg=1
      isoil=isoilm(idjd)
      if (isoil>0) then
        write(6,*)'isoil,iveg,month,fracsum,rlatt: ',isoil,iveg,imo,fracsum(imo),rlatt(idjd)
        fracs=sign(1.,rlatt(idjd))*fracsum(imo) ! +ve for local summer
        fracwet=(.5+fracs)*fracwets(iveg)+(.5-fracs)*fracwetw(iveg)
        write(6,*)'fracs,fracwet,initial_wb: ',fracs,fracwet,wb(idjd,ms)
      end if
    end if
  endif       !  ((nrungcm==-1.or.nrungcm==-2.or.nrungcm==-5)

  if(nrungcm==4)then !  wb fix for ncep input 
    ! this is related to eak's temporary fix for soil moisture
    ! - to compensate for ncep drying out, increase minimum value
    do k=1,ms
      do iq=1,ifull     
        isoil=isoilm(iq)
        wb(iq,k)=min( sfc(isoil) , max(.75*swilt(isoil)+.25*sfc(isoil),wb(iq,k)) )
        ! for antarctic snow
        if(land(iq).and.rlatt(iq)*180./pi<-60.)snowd(iq)=max(snowd(iq),400.)
      enddo   ! iq loop
    enddo    !  k loop
  endif      !  (nrungcm==4)

  if(nrungcm==5)then !  tgg, wb fix for mark 3 input
    ! unfortunately mk 3 only writes out 2 levels
    ! wb just saved as excess above wilting; top level & integrated values
    ! tgg saved for levels 2 and ms 
    do iq=1,ifull     
      isoil=isoilm(iq)
      do k=2,3
        wb(iq,k)=wb(iq,ms)
      enddo    !  k loop
      do k=1,ms
        ! wb(iq,k)=min( sfc(isoil) ,wb(iq,k)+swilt(isoil) ) ! till 22/8/02
        wb(iq,k)=wb(iq,k)+swilt(isoil) 
      enddo    !  k loop
      tgg(iq,3)=.75*tgg(iq,2)+.25*tgg(iq,6)
      tgg(iq,4)= .5*tgg(iq,2)+ .5*tgg(iq,6)
      tgg(iq,5)=.25*tgg(iq,2)+.75*tgg(iq,6)
      ! fix for antarctic snow
      if(land(iq).and.rlatt(iq)*180./pi<-60.)snowd(iq)=max(snowd(iq),400.)
    enddo   ! iq loop
    if (mydiag) then
      write(6,*) 'after nrungcm=5 fixup of mk3 soil variables:'
      write(6,*) 'tgg ',(tgg(idjd,k),k=1,ms)
      write(6,*) 'wb ',(wb(idjd,k),k=1,ms)
    end if
  endif      !  (nrungcm==5)

  if(nrungcm==1)then  ! jlm alternative wb fix for nsib runs off early mark 2 gcm
    if (mydiag ) then
      isoil = isoilm(idjd)
      write(6,"('before nrungcm=1 fix-up wb(1-ms)',9f7.3)") (wb(idjd,k),k=1,ms)
      write(6,*)'nfixwb,isoil,swilt,sfc,ssat,alb ',nfixwb,isoil,swilt(isoil),sfc(isoil),ssat(isoil),albvisnir(idjd,1)
    end if
    do ip=1,ipland  ! all land points in this nsib=1+ loop
      iq=iperm(ip)
      isoil = isoilm(iq)

      if(nfixwb==0)then
        ! very dry jlm suggestion. assume vegfrac ~.5, so try to satisfy
        ! wb0/.36=.5*(wb/sfc + (wb-swilt)/(sfc-swilt) )
        wb(iq,1)=( sfc(isoil)*(sfc(isoil)-swilt(isoil))*wb(iq,1)/.36 +.5*sfc(isoil)*swilt(isoil) )/(sfc(isoil)-.5*swilt(isoil))
        do k=2,ms
          wb(iq,k)=( sfc(isoil)*(sfc(isoil)-swilt(isoil))*wb(iq,ms)/.36+.5*sfc(isoil)*swilt(isoil) )/(sfc(isoil)-.5*swilt(isoil))
        enddo   !  k=2,ms
      endif   ! (nfixwb==0)
      if(nfixwb==1.or.nfixwb==2)then
        ! alternative simpler jlm fix-up	
        ! wb0/.36=(wb-swilt)/(sfc-swilt)
        wb(iq,1)=swilt(isoil)+(sfc(isoil)-swilt(isoil))*wb(iq,1)/.36
        do k=2,ms
          wb(iq,k)=swilt(isoil)+(sfc(isoil)-swilt(isoil))*wb(iq,ms)/.36
        enddo   !  k=2,ms
      endif   ! (nfixwb==1.or.nfixwb==2)
      if(ip==1)write(6,*)'kdate ',kdate
      if(nfixwb==2.and.kdate>3210100.and.kdate<3210200)then
        rlon_d=rlongg(iq)*180./pi
        rlat_d=rlatt(iq)*180./pi
        if(ip==1)then
          write(6,*)'kdate in nfixwb=2 ',kdate
          write(6,*)'iq,rlon_d,rlat_d ',rlon_d,rlat_d
        endif
        ! jlm fix-up for tropical oz in january 321
        if(rlon_d>130..and.rlon_d<150..and.rlat_d>-20..and.rlat_d<0.)then
          do k=1,ms
            wb(iq,k)=max(wb(iq,k),.5*(swilt(isoil)+sfc(isoil))) ! tropics
          enddo   !  k=1,ms
        endif
        ! jlm fix-up for dry interior in january 321
        if(rlon_d>117..and.rlon_d<142..and.rlat_d>-32..and.rlat_d<-22.)then
          do k=1,ms
            wb(iq,k)=swilt(isoil)  ! dry interior
          enddo   !  k=1,ms
        endif
      endif   ! (nfixwb==2)
      if(nfixwb==10)then    ! was default for nrungcm=1 till oct 2001
        ! jlm suggestion, assume vegfrac ~.5, so try to satisfy
        ! wb0/.36=.5*(wb/ssat + (wb-swilt)/(ssat-swilt) )
        wb(iq,1)=( ssat(isoil)*(ssat(isoil)-swilt(isoil))*wb(iq,1)/.36+.5*ssat(isoil)*swilt(isoil) )    &
                /(ssat(isoil)-.5*swilt(isoil))
        do k=2,ms                                                       
          wb(iq,k)=( ssat(isoil)*(ssat(isoil)-swilt(isoil))*wb(iq,ms)/.36+.5*ssat(isoil)*swilt(isoil) ) &
                  /(ssat(isoil)-.5*swilt(isoil))
        enddo   !  k=2,ms
      endif   ! (nfixwb.ne.10)

      do k=1,ms
        wb(iq,k)=max( swilt(isoil) , min(wb(iq,k),sfc(isoil)) )
      enddo     !  k=1,ms
    enddo        !  ip=1,ipland
    if (mydiag) then
      write(6,"('after nrungcm=1 fix-up wb(1-ms)',9f7.3)") (wb(idjd,k),k=1,ms)
      write(6,"('wbice(1-ms)',9f7.3)")(wbice(idjd,k),k=1,ms)
      write(6,"('wb3frac#',9f8.2)") (diagvals(wb(:,3)) - swilt(diagvals(isoilm)))     &
                                  / (sfc(diagvals(isoilm)) - swilt(diagvals(isoilm)))
    end if
  endif          !  (nrungcm==1)

  if(nrungcm==2)then  ! for nsib runs off early mark 2 gcm
    if (mydiag) then
      isoil = isoilm(idjd)
      write(6,*)'before nrungcm=2 fix-up wb(1-ms): ',wb(idjd,:)
      write(6,*)'isoil,swilt,ssat,alb ',isoil,swilt(isoil),ssat(isoil),albvisnir(idjd,1)
    end if
    do ip=1,ipland  ! all land points in this nsib=1+ loop
      iq=iperm(ip)
      isoil = isoilm(iq)
      if( albvisnir(iq,1) >= 0.25 ) then
        diffg=max(0. , wb(iq,1)-0.068)*ssat(isoil)/0.395   ! for sib3
        diffb=max(0. , wb(iq,ms)-0.068)*ssat(isoil)/0.395  ! for sib3
      else
        diffg=max(0. , wb(iq,1)-0.175)*ssat(isoil)/0.42    ! for sib3
        diffb=max(0. , wb(iq,ms)-0.175)*ssat(isoil)/0.42   ! for sib3
      endif
      wb(iq,1)=swilt(isoil)+diffg          ! for sib3
      do k=2,ms                            ! for sib3
        wb(iq,k)=swilt(isoil)+diffb         ! for sib3
      enddo     !  k=2,ms
    enddo      !  ip=1,ipland
    if(mydiag) write(6,*)'after nrungcm=2 fix-up wb(1-ms): ',wb(idjd,:)
  endif          !  (nrungcm==2)
!~ end if ! ( .not.lrestart )

!--------------------------------------------------------------
! FINAL FIXES FOR SURFACE DATA

! soil moisture fixes
do iq=1,ifull
  if(.not.land(iq))then
    wb(iq,:)=0.   ! default over ocean (for plotting)
  endif    !  (.not.land(iq))
enddo     ! iq loop

! snow and ice fixes
snalb=.8
do iq=1,ifull

    !~ if(.not.land(iq))then
      !~ ! from June '03 tgg1	holds actual sea temp, tss holds net temp 
      !~ tgg(iq,1)=max(271.3,tss(iq)) 
      !~ tggsn(iq,1)=tss(iq)         ! a default
    !~ endif   ! (.not.land(iq))
    if(sicedep(iq)>0.)then
      ! at beginning of month set sice temperatures
      tggsn(iq,1)=min(271.2,tss(iq),t(iq,1)+.04*6.5) ! for 40 m level 1
      tss(iq)=tggsn(iq,1)*fracice(iq)+tgg(iq,1)*(1.-fracice(iq))
      albvisnir(iq,1)=.8*fracice(iq)+.1*(1.-fracice(iq))
      albvisnir(iq,2)=.5*fracice(iq)+.1*(1.-fracice(iq))
    endif   ! (sicedep(iq)>0.)

  !~ if(isoilm(iq)==9.and.(nsib==3.or.nsib==5))then
    ! also at beg. of month ensure cold deep temps over permanent ice
    !~ do k=2,ms
      !~ tgg(iq,k)=min(tgg(iq,k),273.1) ! was 260
      !~ wb(iq,k)=max(wb(iq,k),sfc(9))  ! restart value may exceed sfc
      !~ if(wbice(iq,k)<=0.)wbice(iq,k)=.8*wb(iq,k)  ! Dec 07       
    !~ enddo
  !~ endif   ! (isoilm(iq)==9)
enddo    ! iq loop
tpan(1:ifull)=t(1:ifull,1) ! default for land_sflux and outcdf

if ( mydiag ) then
  write(6,*)'near end of indata id+-1, jd+-1'
  write(6,"(' tss#    ',9f8.2)") diagvals(tss)
  write(6,"(' tgg(1)# ',9f8.2)") diagvals(tgg(:,1))
  write(6,"(' tgg(2)# ',9f8.2)") diagvals(tgg(:,2))
  write(6,"(' tgg(3)# ',9f8.2)") diagvals(tgg(:,3))
  write(6,"(' tgg(ms)#',9f8.2)") diagvals(tgg(:,ms))
  write(6,"(' land#   ',9l8)")  diagvals(land)
  write(6,"(' sicedep#   ',9f8.2)") diagvals(sicedep)
  write(6,*)'following from rdnsib'
  write(6,"(' zo#     ',9f8.2)") diagvals(zolnd)
  write(6,"(' wb(1)#  ',9f8.3)") diagvals(wb(:,1))
  wb_tmpry(1:ms) = wb(idjd,1:ms)
  write(6,*)' wb(1-ms): ',wb_tmpry(1:ms)
  write(6,"(' wb(ms)# ',9f8.3)") diagvals(wb(:,ms))
  swilt_diag(1:9) = swilt(diagvals(isoilm))
  sfc_diag(1:9) = sfc(diagvals(isoilm))
  write(6,"(' swilt#  ',9f8.3)") swilt_diag(:)
  write(6,"(' wb3frac#',9f8.3)") (diagvals(wb(:,3)) - swilt_diag(:)) / (sfc_diag(:) - swilt_diag(:))
  write(6,"(' snowd#  ',9f8.2)") diagvals(snowd)
  write(6,"(' fracice#',9f8.3)") diagvals(fracice)
end if

! general initial checks for wb and wbice
do k=1,ms
  do iq=1,ifull
    isoil=isoilm(iq)
    wb(iq,k)=min(ssat(isoil),wb(iq,k))
    wbice(iq,k)=min(.99*wb(iq,k),wbice(iq,k)) 
    if(isoil/=9.and.wbice(iq,k)<=0.)wbice(iq,k)=min(.99,max(0.,.99*(273.1-tgg(iq,k))/5.))*wb(iq,k) ! jlm
  enddo  ! iq loop
enddo   ! ms

if ( mydiag ) then
  write(6,*)'nearer end of indata id+-1, jd+-1'
  write(6,"('tgg(2)# ',9f8.2)")  diagvals(tgg(:,2))
  write(6,"('tgg(ms)#',9f8.2)")  diagvals(tgg(:,ms))
  write(6,"('wbice(1-ms)',9f7.3)")(wbice(idjd,k),k=1,ms)
end if

!**************************************************************
!**************************************************************
! No changes to input data allowed after this point
!**************************************************************
!**************************************************************

!-----------------------------------------------------------------
! UPDATE GENERAL MODEL VARIABLES

! orography
call bounds(zs,corner=.true.)

if ( mydiag ) then
  write(6,*)'for idjd get = ',idjd,ie(idjd),iw(idjd),in(idjd),is(idjd)
  write(6,*)'with zs: ',zs(idjd),zs(ie(idjd)),zs(iw(idjd)),zs(in(idjd)),zs(is(idjd))
end if

! saved albedo
albvissav(:)=albvisnir(:,1) ! VIS
albnirsav(:)=albvisnir(:,2) ! NIR
      
! surface pressure
ps(1:ifull)=1.e5*exp(psl(1:ifull))

!--------------------------------------------------------------
! UPDATE BIOSPHERE DATA (nsib)
if (nsib==6.or.nsib==7) then
  ! Load CABLE data
  if (myid==0) write(6,*) 'Importing CABLE data'
  call loadtile
end if

!--------------------------------------------------------------     
! WRITE FORT.22 FILE FOR GRID INFO
if(nproc==1)then
  coslong=cos(rlong0*pi/180.)   
  sinlong=sin(rlong0*pi/180.)
  coslat=cos(rlat0*pi/180.)
  sinlat=sin(rlat0*pi/180.)
  polenx=-coslat
  poleny=0.
  polenz=sinlat
  write(22,920)
920     format(46x,'land            isoilm')
  write(22,921)
921     format('   iq     i    j  rlong    rlat    thet    map'          &
               '   sicedep zs(m) alb   ivegt  tss    t1    tgg2   tgg6'  &
               '   wb1   wb6 ')
  do j=1,jl
    do i=1,il
      iq=i+(j-1)*il
      zonx=real(            -polenz*y(iq))
      zony=real(polenz*x(iq)-polenx*z(iq))
      zonz=real(polenx*y(iq)             )
      thet=atan2(-zonx*bx(iq)-zony*by(iq)-zonz*bz(iq),zonx*ax(iq)+zony*ay(iq)+zonz*az(iq))*180./pi
      if(thet<0.)thet=thet+360.
      write(22,922) iq,i,j,rlongg(iq)*180./pi,rlatt(iq)*180./pi,  &
                    thet,em(iq),land(iq),sicedep(iq),zs(iq)/grav, &
                    albvisnir(iq,1),                              &
                    isoilm(iq),ivegt(iq),                         &
                    tss(iq),t(iq,1),tgg(iq,2),tgg(iq,ms),         &
                    wb(iq,1),wb(iq,ms),                           &
                    0.
922      format(i6,2i5,3f8.3,f8.4,l2,f4.1,f7.1,f5.2,2i3,          &
                4f7.1,2f6.2,f5.2)
    enddo
  enddo
endif  ! (nproc==1)


call END_LOG(indata_end)
return
end subroutine indataf

!--------------------------------------------------------------
! READ BIOSPHERIC FILES
subroutine rdnsib

use arrays_m                 ! Atmosphere dyamics prognostic arrays
use cc_mpi                   ! CC MPI routines
use infile                   ! Input file routines
use map_m                    ! Grid map arrays
use nsibd_m                  ! Land-surface arrays
use pbl_m                    ! Boundary layer arrays
use soil_m                   ! Soil and surface data
use soilsnow_m               ! Soil, snow and surface data

implicit none

include 'newmpar.h'          ! Grid parameters
include 'const_phys.h'       ! Physical constants
include 'darcdf.h'           ! Netcdf data
include 'filnames.h'         ! Filenames
include 'parm.h'             ! Model configuration
include 'soilv.h'            ! Soil parameters
      
integer iq, iernc
integer ivegmin, ivegmax, ivegmax_g
integer :: idatafix = 0
integer, dimension(:,:), allocatable :: iduma
integer, dimension(ifull,2) :: idumb
integer, dimension(2) :: dumc
real sibvegver
real, dimension(:,:), allocatable :: duma
real, dimension(ifull,7) :: dumb
logical mismatch

real, parameter :: sibvegversion = 2015. ! version id for input data
real, parameter :: falbdflt = 0.
real, parameter :: frsdflt  = 990.
real, parameter :: fzodflt  = 1.
integer, parameter :: ivegdflt  = 0
integer, parameter :: isoildflt = 0

! isoilm_in holds the raw soil data.  isoilm is the data used
! by CCAM after the in-land water bodies (-1) and ocean (0)
! have been combined as water (0).  In the future, -1 will be
! used to better initialise lakes and salt emissions for
! aerosols

if ( nsib >= 6 ) then
  if ( myid == 0 ) then
    write(6,*) "Start reading of nsib>=6 (CABLE) surface datafiles"
    allocate( duma(ifull_g,3) )
    if ( lncveg == 1 ) then
      write(6,*) "Reading soil data"
      call surfread(duma(:,3),'soilt', netcdfid=ncidveg)
      write(6,*) "Reading albedo data"
      call surfread(duma(:,1),'albvis',netcdfid=ncidveg)
      call surfread(duma(:,2),'albnir',netcdfid=ncidveg)
    end if
    call ccmpi_distribute(dumb(:,1:3),duma(:,1:3))
    deallocate( duma )
  else
    call ccmpi_distribute(dumb(:,1:3))
  end if
  ! communicate netcdf status to all processors
  albvisnir(:,1) = dumb(:,1)
  albvisnir(:,2) = dumb(:,2)
  isoilm_in = nint(dumb(:,3))
  isoilm = max( isoilm_in, 0 )
  zolnd = zobgin ! updated in cable_ccam2.f90
  ivegt = 1      ! updated in cable_ccam2.f90
end if
      
!--------------------------------------------------------------
! CHECK FOR LAND SEA MISMATCHES      
mismatch = .false.
if ( datacheck(land,albvisnir(:,1),'albv',idatafix,falbdflt) ) mismatch = .true.
if ( datacheck(land,albvisnir(:,2),'albn',idatafix,falbdflt) ) mismatch = .true.
if ( nsib < 6 ) then
  if ( datacheck(land,rsmin,'rsmin',idatafix,frsdflt) ) mismatch = .true.
end if
ivegmin = minval( ivegt, land(1:ifull) )
ivegmax = maxval( ivegt, land(1:ifull) )
if ( ivegmin<1 .or. ivegmax>44 ) then
  write(6,*) 'stopping in indata, as ivegt out of range'
  write(6,*) 'ivegmin,ivegmax ',ivegmin,ivegmax
  call ccmpi_abort(-1)
end if
if ( datacheck(land,isoilm,'isoilm',idatafix,isoildflt) ) mismatch = .true.

! --- rescale and patch up vegie data if necessary
if ( nsib/=6 .and. nsib/=7 ) then
  dumc(1)   = ivegmax
  call ccmpi_allreduce(dumc(1:1),dumc(2:2),"max",comm_world)
  ivegmax_g = dumc(2)
  if ( ivegmax_g < 14 ) then
    if ( mydiag ) write(6,*) '**** in this run veg types increased from 1-13 to 32-44'
    do iq = 1,ifull            ! add offset to sib values so 1-13 becomes 32-44
      if ( ivegt(iq) > 0 ) ivegt(iq) = ivegt(iq) + 31
    end do
  end if
end if
 
zolnd(:) = max( zolnd(:), zobgin )

return
end subroutine rdnsib

!--------------------------------------------------------------
! FUNCTIONS TO CHECK DATA
logical function rdatacheck( mask,fld,lbl,idfix,val )

implicit  none
      
include 'newmpar.h'      ! Grid parameters

integer, intent(in) :: idfix
real, intent(in) :: val
real, dimension(ifull), intent(inout) :: fld
logical, dimension(ifull), intent(in) :: mask
character(len=*), intent(in) :: lbl

rdatacheck=xrdatacheck(mask,fld,lbl,idfix,0.,val)
      
return
end function rdatacheck

    
logical function xrdatacheck(mask,fld,lbl,idfix,to,from)
      
use cc_mpi, only : myid ! CC MPI routines
      
implicit none

include 'newmpar.h'      ! Grid parameters      
      
real, intent(in) :: from,to
real, dimension(ifull), intent(inout) :: fld
integer, intent(in) :: idfix
integer iq
logical, dimension(ifull), intent(in) :: mask
logical err
character(len=*), intent(in) :: lbl

if (myid==0) write(6,*)' datacheck: verifying field ',lbl

err =.false.
do iq=1,ifull
  if( mask(iq) ) then
    if( fld(iq)==from ) then
      err = .true.
      if( idfix==1 ) then
        fld(iq) = to
        write(6,'(a,2i4,2(a,1pe12.4))') '  changing iq=',iq,' from',from,' to',to
      else
        write(6,*) '  mismatch at iq=',iq,', value',from
      end if
    end if
  end if
end do

xrdatacheck=err

return
end function xrdatacheck

    
logical function idatacheck( mask,ifld,lbl,idfix,ival )

implicit  none
      
include 'newmpar.h'      ! Grid parameters
      
integer, intent(in) :: idfix,ival
integer, dimension(ifull), intent(inout) :: ifld
logical, dimension(ifull), intent(in) :: mask
character(len=*), intent(in) :: lbl

idatacheck=xidatacheck(mask,ifld,lbl,idfix,0,ival)
      
return
end function idatacheck

    
logical function xidatacheck(mask,ifld,lbl,idfix,ifrom,ito)

use cc_mpi, only : myid ! CC MPI routines
      
implicit none

include 'newmpar.h'      ! Grid parameters      
      
integer, intent(in) :: idfix,ifrom,ito
integer iq
integer, dimension(ifull), intent(inout) :: ifld
logical, dimension(ifull), intent(in) :: mask
logical err
character(len=*), intent(in) :: lbl
      
if(myid==0) write(6,*)' datacheck: verifying field ',lbl
err =.false.
do iq=1,ifull
  if( mask(iq) ) then
    if( ifld(iq)==ifrom ) then
      err = .true.
      if( idfix==1 ) then
        ifld(iq) = ito
        write(6,'(a,2i4,2(a,i4))') '  changing iq=',iq,' from',ifrom,' to',ito
      else
        write(6,*) '  mismatch at iq=',iq,', value',ifrom
      end if
    end if
  end if
end do

xidatacheck = err

return
end function xidatacheck

!--------------------------------------------------------------
! READ INTEGER TEXT FILES
subroutine readint(filename,itss,ifully)
      
use cc_mpi            ! CC MPI routines
 
implicit none
      
include 'newmpar.h'   ! Grid parameters
include 'parm.h'      ! Model configuration
include 'parmgeom.h'  ! Coordinate data
            
integer ifully,ilx,jlx,ierr
integer, dimension(ifully) :: itss
integer, dimension(ifull_g) :: glob2d
real rlong0x,rlat0x,schmidtx,dsx
character(len=*) filename
character(len=47) header

if ( myid == 0 ) then
  write(6,*) 'reading data via readint from ',filename
  open(87,file=filename,status='old')
  read(87,*,iostat=ierr) ilx,jlx,rlong0x,rlat0x,schmidtx,dsx,header
  if ( ierr == 0 ) then
    write(6,*) ilx,jlx,rlong0x,rlat0x,schmidtx,dsx,header
    if(ilx/=il_g.or.jlx/=jl_g.or.abs(rlong0x-rlong0)>1.e-6.or.abs(rlat0x-rlat0)>1.E-6.or.abs(schmidtx-schmidt)>1.E-6) then
      write(6,*) 'wrong data file supplied'
      call ccmpi_abort(-1)
    end if
    read(87,*) glob2d
    close(87)
  else if ( ierr < 0 ) then ! Error, so really unformatted file
    close(87)
    write(6,*) 'now doing unformatted read'
    open(87,file=filename,status='old',form='unformatted')
    read(87) glob2d
    close(87)
  else ! ierr > 0
    write(6,*) 'End of file occurred in readint'
    call ccmpi_abort(-1)
  end if
  if (ifully==ifull) then
    call ccmpi_distribute(itss, glob2d)
  else if (ifully==ifull_g) then
    itss=glob2d
  else
    write(6,*) "ERROR: Invalid ifully for readint"
    call ccmpi_abort(-1)
  end if
  write(6,*) trim(header), glob2d(id+(jd-1)*il_g)
else
  if (ifully==ifull) then
    call ccmpi_distribute(itss)
  end if
end if
return
end subroutine readint

!--------------------------------------------------------------
! INITALISE SOIL PARAMETERS
subroutine insoil
      
use cc_mpi, only : myid ! CC MPI routines
      
implicit none
      
include 'newmpar.h'     ! Grid parameters
include 'soilv.h'       ! Soil parameters
      
! The following common block is for soilsnow.f
real zshh, ww
common/soilzs/zshh(ms+1),ww(ms)

integer isoil, k

do isoil = 1,mxst
  cnsd(isoil)  = sand(isoil)*0.3+clay(isoil)*0.25+silt(isoil)*0.265
  hsbh(isoil)  = hyds(isoil)*abs(sucs(isoil))*bch(isoil) !difsat*etasat
  ibp2(isoil)  = nint(bch(isoil))+2
  i2bp3(isoil) = 2*nint(bch(isoil))+3
  if ( myid == 0 ) then
    write(6,"('isoil,ssat,sfc,swilt,hsbh ',i2,3f7.3,e11.4)") isoil,ssat(isoil),sfc(isoil),swilt(isoil),hsbh(isoil)
  end if
end do
cnsd(9) = 2.51

zshh(1)    = .5*zse(1)        ! not used (jlm)
zshh(ms+1) = .5*zse(ms)       ! not used (jlm)
ww(1) = 1.
do k = 2,ms
  zshh(k) = .5*(zse(k-1)+zse(k))  ! z(k)-z(k-1) (jlm)
  ww(k)   = zse(k)/(zse(k)+zse(k-1))
end do

return
end subroutine insoil
!=======================================================================
end module indata