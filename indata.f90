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
    
subroutine indataf(hourst,jalbfix,lapsbot,isoth,nsig,io_nest)
     
use arrays_m                                     ! Atmosphere dyamics prognostic arrays
use bigxy4_m                                     ! Grid interpolation
use cable_ccam, only : loadcbmparm,loadtile      ! CABLE interface
use cc_mpi                                       ! CC MPI routines
use diag_m                                       ! Diagnostic routines
use epst_m                                       ! Off-centre terms
use extraout_m                                   ! Additional diagnostics
use indices_m                                    ! Grid index arrays
use infile                                       ! Input file routines
use latlong_m                                    ! Lat/lon coordinates
!~ use liqwpar_m                                    ! Cloud water mixing ratios
use map_m                                        ! Grid map arrays
use morepbl_m                                    ! Additional boundary layer diagnostics
use nharrs_m, only : lrestart                    ! Non-hydrostatic atmosphere arrays
use nsibd_m                                      ! Land-surface arrays
use onthefly_m                                   ! Input interpolation routines
use pbl_m                                        ! Boundary layer arrays
use permsurf_m                                   ! Fixed surface arrays
use river                                        ! River routing
use sigs_m                                       ! Atmosphere sigma levels
use soil_m                                       ! Soil and surface data
use soilsnow_m                                   ! Soil, snow and surface data
use timeseries, only : init_ts                   ! Tracer time series
use vecs_m                                       ! Eigenvectors for atmosphere dynamics
use vecsuv_m                                     ! Map to cartesian coordinates
use vegpar_m                                     ! Vegetation arrays
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
include 'trcom2.h'                               ! Station data

integer, parameter :: jlmsigmf=1  ! 1 for jlm fixes to dean's data
integer, parameter :: nfixwb=2    ! 0, 1 or 2; wb fixes with nrungcm=1
integer, parameter :: ntest=0
integer, parameter :: klmax=100   ! Maximum vertical levels

!     for the held-suarez test
real, parameter :: delty = 60.    ! pole to equator variation in equal temperature
real, parameter :: deltheta = 10. ! vertical variation
real, parameter :: rkappa = 2./7.

integer, intent(in) :: jalbfix
integer, intent(inout) :: io_nest
integer ii, imo, indexi, indexl, indexs, ip, iq, isoil, isoth
integer iveg, iyr, jj, k, kdate_sav, ktime_sav, l
integer nface, nn, nsig, i, j, n
integer ierr, ic, jc, iqg, ig, jg
integer isav, jsav, ier, lapsbot

character(len=160) :: co2in,radonin,surfin
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
real epsmax, fracs, fracwet, ftsoil, gwdfac, hefact
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
swrsave(:)=0.5
hourst=0.
albvissav(:)=-1.
albvisnir(:,:)=0.3
vlai(:)=0.
ivegt(:)=1
isoilm(:)=1
zs(:)=0.
zsmask(:)=0.
!~ he(:)=0.         
land(:)=.false.
kdate=kdate_s
ktime=ktime_s

!--------------------------------------------------------------
! READ AND PROCESS ATMOSPHERE SIGMA LEVELS
if (myid==0) then
  read(28,*)(dumc(k),k=1,kl),(dumc(2*kl+k),k=1,kl),    &
        (bam(k),k=1,kl),((emat(k,l),k=1,kl),l=1,kl),   &
        ((einv(k,l),k=1,kl),l=1,kl),(qvec(k),k=1,kl),  &
        ((tmat(k,l),k=1,kl),l=1,kl)
  write(6,*) 'kl,lapsbot,sig from eigenv file: ',kl,lapsbot,dumc(1:kl)
  ! File has an sigmh(kl+1) which is not required. Causes bounds violation
  ! to read this.
  read(28,*)(dumc(kl+k),k=1,kl)
  close(28) 
  write(6,*) 'tbar: ',dumc(2*kl+1:3*kl)
  write(6,*) 'bam:  ',bam
       
  ! test netcdf for CABLE input
  if (nsib>=6) then
    call ccnf_open(vegfile,ncidveg,ierr)
    if (ierr==0) then
      dumc(3*kl+1)=1.
    else
      dumc(3*kl+1)=0.  
    end if
  else if (nsib==5) then
    call ccnf_open(vegfile,ncidveg,ierr)
    if (ierr==0) then
      dumc(3*kl+1)=1.
    else
      dumc(3*kl+1)=0.  
    end if
  else
    dumc(3*kl+1)=0.
  end if
endif ! (myid==0)

! distribute vertical and vegfile data to all processors
! dumc(1:kl)   = sig,   dumc(kl+1:2*kl) = sigmh, dumc(2*kl+1:3*kl) = tbar
! dumc(3*kl+1) = lncveg
call ccmpi_bcast(dumc(1:3*kl+1),0,comm_world)
sig   =dumc(1:kl)
sigmh =dumc(kl+1:2*kl)
tbar  =dumc(2*kl+1:3*kl)
lncveg=nint(dumc(3*kl+1))

dsig(1:kl-1)=sigmh(2:kl)-sigmh(1:kl-1)
dsig(kl)=-sigmh(kl)
sumdsig=0.
do k=1,kl
  sumdsig=sumdsig-dsig(k)
  tbardsig(k)=0.
enddo
if (myid==0) write(6,*)'dsig,sumdsig ',dsig,sumdsig
if (isoth>=0) then
  dtmax=1./(sig(1)*log(sig(1)/sig(2)))
  tbardsig(1)=dtmax*(tbar(1)-tbar(2))
  do k=2,kl-1
    tbardsig(k)=(tbar(k+1)-tbar(k-1))/(2.*dsig(k))
  enddo
endif
!     rata and ratb are used to interpolate half level values to full levels
!     ratha and rathb are used to interpolate full level values to half levels
ratha(kl) = 0. ! not used
rathb(kl) = 0. ! not used
rata(kl)=(sigmh(kl)-sig(kl))/sigmh(kl)
ratb(kl)=sig(kl)/sigmh(kl)
do k=1,kl-1
  bet(k+1)=rdry*log(sig(k)/sig(k+1))*.5
  rata(k)=(sigmh(k)-sig(k))/(sigmh(k)-sigmh(k+1))
  ratb(k)=(sig(k)-sigmh(k+1))/(sigmh(k)-sigmh(k+1))
  ratha(k)=(sigmh(k+1)-sig(k))/(sig(k+1)-sig(k))
  rathb(k)=(sig(k+1)-sigmh(k+1))/(sig(k+1)-sig(k))
enddo
if (myid==0) then
  write(6,*)'rata ',rata
  write(6,*)'ratb ',ratb
  write(6,*)'ratha ',ratha
  write(6,*)'rathb ',rathb
end if
   
c=grav/stdlapse
bet(1)=c *(sig(1)**(-rdry/c)-1.)
if(lapsbot==1)bet(1)=-rdry*log(sig(1))
betm(1:kl)=bet(1:kl)
if(lapsbot==2)then     ! may need refinement for non-equal spacing
  do k=2,kl
    bet(k)=.5*rdry*(sig(k-1)-sig(k))/sig(k)
    betm(k)=.5*rdry*(sig(k-1)-sig(k))/sig(k-1)
  enddo
  bet(1)=rdry*(1.-sig(1))/sig(1)
elseif(lapsbot==3)then ! possibly suits nh
  betm(:)=0.
  do k=2,kl
    bet(k)=rdry*log(sig(k-1)/sig(k))
  enddo
  bet(1)=-rdry*log(sig(1))
endif

if ( myid==0 ) then
  write(6,*) 'bet  ',bet
  write(6,*) 'betm ',betm
end if
!if (nh/=0) then
! Non-hydrostatic case
!~ if ( nh==2 .and. lapsbot/=3 ) stop 'nh=2 needs lapsbot=3'
!~ if ( abs(epsp)<=1. ) then
  !~ ! exact treatment when epsp is constant
  !~ call eig(sig,sigmh,tbar,lapsbot,isoth,dt,epsp,epsh,nsig,bet,betm,nh)
!~ else
  !~ call eig(sig,sigmh,tbar,lapsbot,isoth,dt,0.,0.,nsig,bet,betm,nh)
!~ end if
!else
!  ! MJT notes - The hydrostatic case could have called
!  ! eig and avoided a ccmpi_bcast.  However, since eig
!  ! does not always exactly reproduce the input file 
!  ! emat, einv and bam, then we keep the bcast for 
!  ! backwards compatibility
!  call ccmpi_bcast(bam,0,comm_world)
!  call ccmpi_bcast(emat,0,comm_world)
!  call ccmpi_bcast(einv,0,comm_world)
!endif  ! (nh/=0)


! zmin here is approx height of the lowest level in the model
zmin = -rdry*280.*log(sig(1))/grav
if ( myid==0 ) write(6,*) 'zmin = ',zmin


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
  !~ he(:)=0.         
  land(:)=.false.
endif  ! (io_in<=4.and.nhstest>=0)  ..else..

if ( mydiag ) then
  write(6,"('zs#_topof ',9f8.1)") diagvals(zs)
  !~ write(6,"('he#_topof ',9f8.1)") diagvals(he)
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
    call loadcbmparm(vegfile,vegprev,vegnext,phenfile,casafile)
    ! albvisnir at this point holds net albedo
  elseif (nsib==3) then
    ! special options for standard land surface scheme
    if(nspecial==35)then      ! test for Andy Cottrill
      do iq=1,ifull
        rlongd=rlongg(iq)*180./pi
        rlatd=rlatt(iq)*180./pi
        if(rlatd>-32..and.rlatd<-23.5)then
          if(rlongd>145..and.rlongd<=150.)ivegt(iq)=4
          if(rlongd>150..and.rlongd<154.)ivegt(iq)=2
        endif
      enddo
    endif  ! (nspecial==35)
    ! zap vegetation over SEQ for Andy
    if(nspecial==41)then
      do iq=1,ifull
        rlongd=rlongg(iq)*180./pi
        rlatd=rlatt(iq)*180./pi
        if(rlatd>-32. .and. rlatd<-23.5)then
          if(rlongd>145. .and. rlongd<=152.)ivegt(iq)=4 
          if(rlongd>152. .and. rlongd< 154.)ivegt(iq)=2 
        endif
      enddo
    endif  ! (nspecial==41)
    do iq=1,ifull
      ! check for littoral veg over Oz      
      rlongd=rlongg(iq)*180./pi
      rlatd=rlatt(iq)*180./pi
      if(rlongd>110.and.rlongd<155.and.rlatd>-45.and.rlatd<-10)then
        if(ivegt(iq)==28)then
          write(6,*)'littoral vegt ',iq,rlongd,rlatd
          if(rlongd>150.and.rlongd<152.and.rlatd>-28.and.rlatd<-26)then
            ivegt(iq)=24   ! fix-up of Graetz data for Andy from July '07
          endif
        endif
      endif
    enddo
    do iq=1,ifull
      if(land(iq))then  
        ! following line from 16/3/06 avoids sand on Himalayas        
        if(zs(iq)>2000.*grav.and.isoilm(iq)<3)isoilm(iq)=3
      endif
    enddo
    ! put in Antarctica ice-shelf fixes 5/3/07
    do iq=1,ifull
      if(zs(iq)<=0.)then
        rlongd=rlongg(iq)*180./pi
        rlatd=rlatt(iq)*180./pi
        if((rlongd>165..and.rlongd<195..and.rlatd<-77.2-(rlongd-165.)/30.).or.      & ! Ross shelf
           (rlongd>300..and.rlongd<330..and.rlatd<-75.2-(rlongd-300.)*2.8/30.).or.  & ! Ronne shelf
           (rlongd>68..and.rlongd<75..and.rlatd<-64.-(rlongd-60.)*5.8/15.))then       ! Amery shelf
             zs(iq)=1.
             land(iq)=.true.
             isoilm(iq)=9
             ivegt(iq)=42
             if(mydiag)write(6,*)'setting sea to ice sheet for iq = ',iq
        endif
      endif  ! (zs(iq)<=0.)
    enddo
  end if ! (nsib==6.or.nsib==7) ..else..
  if (nsib/=6.and.nsib/=7) then
    ! JJK special option for adjusting surface albedo and roughness
    if(nspecial<-10)then
      do iq=1,ifull
        rlongd=rlongg(iq)*180./pi
        rlatd=rlatt(iq)*180./pi
        if((rlatdn<rlatd .and. rlatd<rlatdx).and.(rlongdn<rlongd .and. rlongd<rlongdx))then
          ! assume nspecial = -vvir where vv = % vis alb and ir = % nir alb
          newzo=real(abs(nspecial)/10000)
          visalb=real((abs(nspecial)-nint(newzo)*10000)/100)
          niralb=real(abs(nspecial)-nint(newzo)*10000-nint(visalb)*100)
          if ( newzo  > 1. ) zolnd(iq)=newzo/1000. ! input mm, output m
          if ( visalb > 1. ) albvisnir(iq,1)=visalb/100. ! (to make 0-1)
          if ( niralb > 1. ) albvisnir(iq,2)=niralb/100. ! (to make 0-1)
          if ( .not. land(iq) ) then
            land(iq)=.true.
            ivegt(iq)=42
            isoilm(iq)=7
          end if
        endif
      enddo
    endif  ! (nspecial<-10)	
  end if ! (nsib/=6.and.nsib/=7)
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
if (nmlo/=0.and.abs(nmlo)<=9) then
  if (myid==0) write(6,*) 'Initialising MLO'
  call surfread(dep,'depth',filename=bathfile)
  where (land)
    dep=0.
  end where
end if
if ( abs(nmlo)>=2 .or. nriver==1 ) then
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
				tggsn,smass,ssdn,ssdnn,snage,isflag, ocndwn)
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

  ! adjust input for differences in orography
  if(newtop==2)then
    ! reduce sea tss to mslp      e.g. for qcca in ncep gcm
    do iq=1,ifull
      if(tss(iq)<0.)then
        if (abs(zss(iq))>1000.) then
          write(6,*)'zss,tss_sea in, out',iq,zss(iq),tss(iq),tss(iq)-zss(iq)*stdlapse/grav
        end if
        tss(iq)=tss(iq)-zss(iq)*stdlapse/grav
      endif
    enddo
  endif                  ! (newtop==2)
  tss(1:ifull)=abs(tss(1:ifull)) ! not done in infile because -ve needed for onthefly
  hourst=.01*ktime
  if ( myid==0 ) then
    write(6,*)'rlongg(1),rlongg(ifull) ',rlongg(1),rlongg(ifull)
    write(6,*)'using em: ',(em(ii),ii=1,10)
    write(6,*)'using  f: ',(f(ii),ii=1,10)
    write(6,*)'in indata hourst = ',hourst
    write(6,*)'sigmas: ',sig
    write(6,*)'sigmh: ',sigmh
  end if
  if ( mydiag ) then
    write(6,*)'newtop, zsold, zs,tss_in,land ',newtop,zss(idjd),zs(idjd),tss(idjd),land(idjd)
  end if
  if (newtop>=1.and..not.lrestart) then    
    if (nproc==1) then
      pslavge=sum(psl(1:ifull)*wts(1:ifull))
      write (6,"(' initial pslavge ',f10.6)") pslavge
    endif 
    do iq=1,ifull
      if (land(iq)) then
        tss(iq)=tss(iq)+(zss(iq)-zs(iq))*stdlapse/grav
        do k=1,ms
          tgg(iq,k)=tgg(iq,k)+(zss(iq)-zs(iq))*stdlapse/grav
        enddo
      endif     ! (land(iq))
    enddo        ! iq loop
    if ( mydiag ) then
      write(6,*)'newtop>=1 new_land_tss,zsold,zs: ',tss(idjd),zss(idjd),zs(idjd)
      ! compensate psl, t(,,1), qg as read in from infile
      write(6,"(' zs#  in     ',9f8.1)") diagvals(zs)
      write(6,"(' zss# in     ',9f8.1)") diagvals(zss)
      write(6,"(' 100*psl#  in',9f8.2)") 100.*diagvals(psl)
      write(6,*) 'now call retopo from indata'
    end if ! ( mydiag )
    call retopo(psl,zss,zs,t,qg)
    if(nmaxpr==1.and.mydiag)then
      write(6,"(' 100*psl# out',9f8.2)") 100.*diagvals(psl)
    endif
    if (nproc==1) then
      pslavge=sum(psl(1:ifull)*wts(1:ifull))
      write (6,"(' after retopo pslavge ',f10.6)") pslavge
    endif 
  endif   ! (newtop>=1.and..not.lrestart)

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
if ( .not.lrestart ) then
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
end if ! ( .not.lrestart )




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
  if (nmlo==0.or.abs(nmlo)>9) then
    if(.not.land(iq))then
      ! from June '03 tgg1	holds actual sea temp, tss holds net temp 
      tgg(iq,1)=max(271.3,tss(iq)) 
      tggsn(iq,1)=tss(iq)         ! a default
    endif   ! (.not.land(iq))
    if(sicedep(iq)>0.)then
      ! at beginning of month set sice temperatures
      tggsn(iq,1)=min(271.2,tss(iq),t(iq,1)+.04*6.5) ! for 40 m level 1
      tss(iq)=tggsn(iq,1)*fracice(iq)+tgg(iq,1)*(1.-fracice(iq))
      albvisnir(iq,1)=.8*fracice(iq)+.1*(1.-fracice(iq))
      albvisnir(iq,2)=.5*fracice(iq)+.1*(1.-fracice(iq))
    endif   ! (sicedep(iq)>0.)
  endif    ! (nmlo==0.or.abs(nmlo)>9) 
  if(isoilm(iq)==9.and.(nsib==3.or.nsib==5))then
    ! also at beg. of month ensure cold deep temps over permanent ice
    do k=2,ms
      tgg(iq,k)=min(tgg(iq,k),273.1) ! was 260
      wb(iq,k)=max(wb(iq,k),sfc(9))  ! restart value may exceed sfc
      if(wbice(iq,k)<=0.)wbice(iq,k)=.8*wb(iq,k)  ! Dec 07       
    enddo
  endif   ! (isoilm(iq)==9)
enddo    ! iq loop
tpan(1:ifull)=t(1:ifull,1) ! default for land_sflux and outcdf

! albedo and roughness fixes
select case(nsib)
  case(5)
    ! MODIS input with standard surface scheme
    osnowd = snowd
    zolog=log(zmin/zolnd)   ! for land use in sflux
    sigmf=0.
    where (land)
      sigmf(:)=max(0.01,min(0.98,1.-exp(-0.4*vlai(:))))
    elsewhere
      vlai=0.
    end where
  case(3)
    ! usual input with standard surface scheme
    osnowd = snowd
    sigmf=0.
    do iq=1,ifull
      if(land(iq))then
        isoil = isoilm(iq)
        iveg  = ivegt(iq)
        sigmf(iq)=min(.8,.95*vegpsig(ivegt(iq)))  ! moved from rdnsib
        if(jlmsigmf==1)then  ! fix-up for dean's veg-fraction
          sigmf(iq)=((sfc(isoil)-wb(iq,3))*vegpmin(iveg)+(wb(iq,3)-swilt(isoil))*vegpmax(iveg))/(sfc(isoil)-swilt(isoil)) 
          sigmf(iq)=max(vegpmin(iveg),min(sigmf(iq),.8)) ! for odd wb
        endif   ! (jlmsigmf==1)
        ! following just for rsmin diags for nstn and outcdf	
        tstom=298.
        if(iveg==6+31)tstom=302.
        if(iveg>=10.and.iveg<=21.and.abs(rlatt(iq)*180./pi)<25.)tstom=302.
        tsoil=min(tstom, .5*(.3333*tgg(iq,2)+.6667*tgg(iq,3)+.95*tgg(iq,4) + .05*tgg(iq,5)))
        ftsoil=max(0.,1.-.0016*(tstom-tsoil)**2)
        ! which is same as:  ftsoil=max(0.,1.-.0016*(tstom-tsoil)**2)
        !                    if( tsoil >= tstom ) ftsoil=1.
        rlai(iq)=  max(.1,rlaim44(iveg)-slveg44(iveg)*(1.-ftsoil))
        rsmin(iq) = rsunc44(iveg)/rlai(iq)
      endif   ! (land(iq)) 
    enddo    !  iq loop
    if(jalbfix==1)then ! jlm fix for albedos, esp. for bare sandy soil
      if(mydiag)then
        isoil=isoilm(idjd)
        if (isoil.gt.0) then
          write(6,*)'before jalbfix isoil,sand,alb,rsmin ',isoil,sand(isoil),albvisnir(idjd,1),rsmin(idjd)
        else
          write(6,*)'before jalbfix isoil,sand,alb,rsmin ',isoil,0.,albvisnir(idjd,1),rsmin(idjd)
        end if
      endif
      do ip=1,ipland  
        iq=iperm(ip)
        isoil = isoilm(iq)
        albvisnir(iq,1)=max(albvisnir(iq,1),sigmf(iq)*albvisnir(iq,1)+(1.-sigmf(iq))*(sand(isoil)*.35+(1.-sand(isoil))*.06))
        albvisnir(iq,2)=albvisnir(iq,1)
      enddo                  !  ip=1,ipland
      if(mydiag)then
        write(6,*)'after jalbfix sigmf,alb ',sigmf(idjd),albvisnir(idjd,1)
      endif
    endif  ! (jalbfix==1)
    if(newrough>0)then
      call calczo
      if(mydiag)write(6,*)'after calczo zolnd ',zolnd(idjd)
      if ( mydiag ) then
        write(6,*)'after calczo with newrough = ',newrough
        write(6,"('zo#    ',9f8.2)") diagvals(zolnd)
      end if
    endif ! (newrough>0)
    zolog=log(zmin/zolnd)   ! for land use in sflux 
end select

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
! UPDATE GRAVITY WAVE DRAG DATA (lgwd)
gwdfac=.01*lgwd       ! most runs used .02 up to fri  10-10-1997
hefact=.1*abs(ngwd)   ! hal used hefact=1. (equiv to ngwd=10)
if(myid==0)write(6,*)'hefact,helim,gwdfac: ',hefact,helim,gwdfac

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
               '   wb1   wb6   ico2  radon')
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


!--------------------------------------------------------------
! INITIALISE STATION OUTPUT (nstn)
if(nstn>0)then
  if (myid==0) then
    write(6,*) 'land stations'
    write(*,"(a)") ' lu istn jstn  iq   slon   slat land rlong  rlat' &
                // ' isoil iveg zs(m) alb  wb3  wet3 vlai  zo   he'
  end if
  call ccmpi_barrier(comm_world)
  do nn=1,nstn
    call latltoij(slon(nn),slat(nn),rlong0,rlat0,schmidt,ri,rj,nface,xx4,yy4,il_g)
    ! These are global indices
    ig=nint(ri)
    jg=nint(rj) + nface*il_g
    mystn(nn) = fproc(ig,nint(rj),nface) == myid
    if ( mystn(nn) ) then
      iqg = ig + (jg-1)*il_g
      deli = nint(ri) - ri
      delj = nint(rj) - rj
      ! Local indices on this processor
      call indv_mpi(iqg,ii,jj,n)
      iq = ii + (jj-1)*ipan + (n-1)*ipan*jpan
      if(.not.land(iq))then
        ! simple search for neighbouring land point 
        ! N.B. does not search over panel/processor bndries
        ! N.B. if no land points, just returns closest point
        isav = ii
        jsav = jj
        dist=100.
        if ( isav < ipan ) then
          distnew = (deli+1)**2 + delj**2 
          if(land(iq+1).and.distnew<dist)then
            ii=isav+1
            dist=distnew
          endif
        end if
        if ( isav > 1 ) then
          distnew = (deli-1)**2 + delj**2 
          if(land(iq-1).and.distnew<dist)then
            ii=isav-1
            dist=distnew
          endif
        end if
        if ( jsav < jpan ) then
          distnew = deli**2 + (delj+1)**2 
          if(land(iq+ipan).and.distnew<dist)then
            jj=jsav+1
            dist=distnew
          endif
        end if
        if ( jsav >= 1 ) then
          distnew = deli**2 +(delj-1)**2 
          if(land(iq-ipan).and.distnew<dist)then
            jj=jsav-1
            dist=distnew
          endif
        end if
      endif              ! (.not.land(iq))
      istn(nn) = ii
      jstn(nn) = jj+(n-1)*jpan
      iq = istn(nn) + (jstn(nn)-1)*ipan
      iveg=ivegt(iq)
      isoil = isoilm(iq)
      wet3=(wb(iq,3)-swilt(isoil))/(sfc(isoil)-swilt(isoil))    
    end if               ! mystn
98  format(i3,i4,i5,i6,2f7.2 ,l3,2f7.2, i3,i6,f7.1,f5.2,4f5.2,f7.1,i4)
    ! Put a barrier here to force stations to be printed in the right order
    call ccmpi_barrier(comm_world)
  enddo  ! nn=1,nstn
endif     !  (nstn>0)

    
!--------------------------------------------------------------
! OPEN MESONEST FILE
if ( mbd/=0 .or. nbd/=0 ) then
  if ( myid == 0 ) then
    write(6,*) "Opening mesonest file"
  end if
  io_in = io_nest                  ! Needs to be seen by all processors
  call histopen(ncid,mesonest,ier) ! open parallel mesonest files
  call ncmsg("mesonest",ier)       ! report error messages
  if ( myid == 0 ) then
    write(6,*) "ncid,mesonest ",ncid,trim(mesonest)
  end if
end if    ! (mbd/=0.or.nbd/=0)       


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
!~ use tracers_m                ! Tracer data
use vegpar_m                 ! Vegetation arrays

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
    else
      write(6,*) "Cannot open vegfile as a netcdf file ",vegfile
      write(6,*) "Assuming ASCII file format"
      call surfread(duma(:,3),'soilt', filename=soilfile)
      call surfread(duma(:,1),'albvis',filename=albfile)
      call surfread(duma(:,2),'albnir',filename=albnirfile)
      duma(:,1:2) = 0.01*duma(:,1:2)
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
! READ REAL TEXT FILES
subroutine readreal(filename,tss,ifully)
 
use cc_mpi            ! CC MPI routines
 
implicit none
      
include 'newmpar.h'   ! Grid parameters
include 'parm.h'      ! Model configuration
include 'parmgeom.h'  ! Coordinate data

integer ierr
integer ilx,jlx,ifully
real, dimension(ifully) :: tss
real, dimension(ifull_g) :: glob2d
real rlong0x,rlat0x,schmidtx,dsx
character(len=*) filename
character(len=47) header

if ( myid == 0 ) then
  write(6,*) 'reading data via readreal from ',trim(filename)
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
  else
    write(6,*) "error in readreal",trim(filename),ierr
    call ccmpi_abort(-1)
  end if
  if (ifully==ifull) then
    call ccmpi_distribute(tss, glob2d)
  else if (ifully==ifull_g) then
    tss=glob2d
  else
    write(6,*) "ERROR: Invalid ifully for readreal"
    call ccmpi_abort(-1)
  end if
  write(6,*) trim(header), glob2d(id+(jd-1)*il_g)
else
  if (ifully==ifull) then
    call ccmpi_distribute(tss)
  end if
end if
return
end subroutine readreal


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

      
!--------------------------------------------------------------
! CALCULATE ROUGHNESS LENGTH
subroutine calczo
      
use arrays_m        ! Atmosphere dyamics prognostic arrays
use map_m           ! Grid map arrays
use nsibd_m         ! Land-surface arrays
use soil_m          ! Soil and surface data
use soilsnow_m      ! Soil, snow and surface data
      
implicit none
      
include 'newmpar.h' ! Grid parameters
include 'parm.h'    ! Model configuration
      
integer iq,iveg
real zomax,zomin,tsoil,sdep
      
real xhc(0:44)
!     vegetation height
data xhc    / 0.0,                                                   & ! 0
              30.0,28.0,25.0,17.0,12.0,10.0, 9.0, 7.0, 5.5, 3.0,     & ! 1-10
              2.5, 2.0, 1.0, 0.6, 0.5, 0.5,0.45,0.75, 0.6,0.45,      &
              0.4, 0.6, 0.6,0.24,0.25,0.35, 0.3, 2.5, 0.0, 0.0,      &
              0.0,                                                   & ! 31
              32.,20.,20.,17.,17., 1., 1., 1., 0.5, 0.6, 0., 1.,0./    !sellers 1996 j.climate

zomax=-1.e29
zomin= 1.e29
if(newrough==2)then
  do iq=1,ifull
    if(land(iq))then
      iveg=ivegt(iq)
      zolnd(iq)=max(zobgin , .1*xhc(iveg))
      zomax=max(zomax,zolnd(iq))
      zomin=min(zomin,zolnd(iq))
    endif  ! (land(iq))then
  enddo   ! iq loop
elseif(newrough==3)then
  do iq=1,ifull
    if(land(iq))then
      iveg=ivegt(iq)
      zolnd(iq)=max(zobgin , .13*xhc(iveg))  ! French factor
      zomax=max(zomax,zolnd(iq))
      zomin=min(zomin,zolnd(iq))
    endif  ! (land(iq))then
  enddo   ! iq loop
else
  do iq=1,ifull
    if(land(iq))then
      iveg=ivegt(iq)
      tsoil  = 0.5*(tgg(iq,ms)+tgg(iq,2))
      sdep=0.
      call cruf1 (iveg,tsoil,sdep,zolnd(iq),zobgin)
      zomax=max(zomax,zolnd(iq))
      zomin=min(zomin,zolnd(iq))
    endif  ! (land(iq))then
  enddo   ! iq loop
endif

write(6,*)"calczo zolnd: zomin,zomax=",zomin,zomax
return ! calczo
end subroutine calczo

subroutine cruf1(iv,tsoil,sdep,zolnd,zobgin)
      
implicit none
      
! kf, 1997
! for each vegetation type (= iv), assign veg height, total lai, albedo,
! and computed aerodynamic, radiative and interception properties.
! jmax0 assigned due to table by ray leuning and estimates  21-08-97 
! apply seasonal variations in lai and height. return via /canopy/
! nb: total lai = xrlai, veglai = xvlai, veg cover fraction = xpfc,
!     with xrlai = xvlai*xpfc
! type  0 to 31: 2d version with graetz veg types
! type 32 to 43: 2d version with gcm veg types
! type 44:       stand-alone version
!-----------------------------------------------------------------------
!   name                             symbol  code hc:cm pfc:%  veglai
!   ocean                                o!     0     0     0  0.0
!   tall dense forest                    t4     1  4200   100  4.8
!   tall mid-dense forest                t3     2  3650    85  6.3
!   dense forest                         m4     3  2500    85  5.0  (?)
!   mid-dense forest                     m3     4  1700    50  3.75
!   sparse forest (woodland)             m2     5  1200    20  2.78
!   very sparse forest (woodland)        m1     6  1000     5  2.5
!   low dense forest                     l4     7   900    85  3.9
!   low mid-dense forest                 l3     8   700    50  2.77
!   low sparse forest (woodland)         l2     9   550    20  2.04
!   tall mid-dense shrubland (scrub)     s3    10   300    50  2.6
!   tall sparse shrubland                s2    11   250    20  1.69
!   tall very sparse shrubland           s1    12   200     5  1.9
!   low mid-dense shrubland              z3    13   100    50  1.37
!   low sparse shrubland                 z2    14    60    20  1.5
!   low very sparse shrubland            z1    15    50     5  1.21
!   sparse hummock grassland             h2    16    50    20  1.58
!   very sparse hummock grassland        h1    17    45     5  1.41
!   dense tussock grassland              g4    18    75    85  2.3
!   mid-dense tussock grassland          g3    19    60    50  1.2
!   sparse tussock grassland             g2    20    45    20  1.71
!   very sparse tussock grassland        g1    21    40     5  1.21
!   dense pasture/herbfield (perennial)  f4    22    60    85  2.3
!   dense pasture/herbfield (seasonal)  f4s    23    60    85  2.3
!   mid-dense pasture/herb (perennial)   f3    24    45    50  1.2
!   mid-dense pasture/herb  (seasonal)  f3s    25    45    50  1.2
!   sparse herbfield*                    f2    26    35    20  1.87
!   very sparse herbfield                f1    27    30     5  1.0
!   littoral                             ll    28   250    50  3.0
!   permanent lake                       pl    29     0     0  0
!   ephemeral lake (salt)                sl    30     0     0  0
!   urban                                 u    31     0     0  0
!   stand alone: hc,rlai from param1      -    44     -   100  -

!   above are dean's. below are sib (added 31 to get model iveg)
!  32  1 - broadleaf evergreen trees (tropical forest)
!  33  2 - broadleaf deciduous trees
!  34  3 - broadleaf and needleaf trees
!  35  4 - needleaf evergreen trees
!  36  5 - needleaf deciduous trees 
!  37  6 - broadleaf trees with ground cover (savannah)
!  38  7 - groundcover only (perennial)
!  39  8 - broadleaf shrubs with groundcover
!  40  9 - broadleaf shrubs with bare soil
!  41 10 - dwarf trees and shrubs with groundcover
!  42 11 - bare soil
!  43 12 - agriculture or C3 grassland (newer defn)
 
!                             soil type
!       texture               
!  0   water/ocean
!  1   coarse               sand/loamy_sand
!  2   medium               clay-loam/silty-clay-loam/silt-loam
!  3   fine                 clay
!  4   coarse-medium        sandy-loam/loam
!  5   coarse-fine          sandy-clay
!  6   medium-fine          silty-clay 
!  7   coarse-medium-fine   sandy-clay-loam
!  8   organi!              peat
!  9   land ice
!-----------------------------------------------------------------------

integer iv
real ftsoil,vrlai,hc,tsoil,sdep,zolnd,zobgin
real rlai,usuh,disp,coexp
real xhc(0:44),xpfc(0:44),xvlai(0:44),xslveg(0:44)
! aerodynamic parameters, diffusivities, water density:
real vonk,a33,csw,ctl
parameter(vonk   = 0.40)     ! von karman constant
parameter(a33    = 1.25)     ! inertial sublayer sw/us
parameter(csw    = 0.50)     ! canopy sw decay (weil theory)
parameter(ctl    = 0.40)     ! wagga wheat (rdd 1992, challenges)
! vegetation height
data xhc    / 0.0,                                                   & ! 0
              30.0,28.0,25.0,17.0,12.0,10.0, 9.0, 7.0, 5.5, 3.0,     & ! 1-10
              2.5, 2.0, 1.0, 0.6, 0.5, 0.5,0.45,0.75, 0.6,0.45,      &
              0.4, 0.6, 0.6,0.24,0.25,0.35, 0.3, 2.5, 0.0, 0.0,      &
              0.0,                                                   & ! 31
              32.,20.,20.,17.,17., 1., 1., 1., 0.5, 0.6, 0., 1.,0./    !sellers 1996 j.climate

! vegetation fractional cover
data xpfc   /0.00,                                                   &
             1.00,0.85,0.85,0.50,0.20,0.05,0.85,0.50,0.20,0.50,      &
             0.20,0.05,0.50,0.20,0.05,0.20,0.05,0.85,0.50,0.20,      &
             0.05,0.85,0.85,0.50,0.50,0.20,0.05,0.50,0.00,0.00,      &
             0.00,                                                   &
             .98,.75,.75,.75,.50,.86,.65,.79,.30,.42,.02,.54,  1.0/

! veg lai from graetz table of 283 veg types (iv=0 to 31), and maximum 
! veg lai for gcm veg types (iv=32 to 43)  stand-alone: 44
data xvlai  / 0.0,                                                   &
              4.80,6.30,5.00,3.75,2.78,2.50,3.90,2.77,2.04,2.60,     &
              1.69,1.90,1.37,1.50,1.21,1.58,1.41,2.30,1.20,1.71,     &
              1.21,2.30,2.30,1.20,1.20,1.87,1.00,3.00,0.00,0.00,     &
              0.00,                                                  &
              6.0,5.0,4.0,4.0,4.0,3.0,3.0,3.0,1.0,4.0,0.5,3.0,  0.0/     ! 32-44

! for seasonally varying lai, amplitude of veg lai seasonal change
data xslveg  /0.00,                                                  &
              0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,     & 
              0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,     &
              0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,     &
              0.00,                                                  &
              2.0,2.0,2.0,2.0,2.0,1.5,1.5,1.5,1.0,0.5,0.5,0.5,  0.0/
!-----------------------------------------------------------------------
! assign aerodynamic, radiative, stomatal, interception properties
! assign total lai (xrlai) from veg lai and pfc, and assign seasonal 
!   variation in lai and veg height where necessary. this is controlled
!   by the factor season (0 =< season =< 1).
ftsoil=max(0.,1.-.0016*(298.-tsoil)**2)
if( tsoil >= 298. ) ftsoil=1.
vrlai = max(0.0,(xvlai(iv)-xslveg(iv)*(1.-ftsoil))*xpfc(iv))
hc    = max(0.0,xhc(iv) - sdep)
rlai  = vrlai*hc/max(0.01,xhc(iv))
!   find roughness length zolnd from hc and rlai:
call cruf2(hc,rlai,usuh,zolnd,disp,coexp)
!   set aerodynamic variables for bare soil and vegetated cases:
zolnd=max(zolnd, zobgin)
if (rlai<0.001 .or. hc<.05) then
  zolnd    = zobgin      ! bare soil surface
  hc     = 0.0  
  rlai   = 0.0
endif
return
end subroutine cruf1
!=======================================================================
subroutine cruf2(h,rlai,usuh,z0,d,coexp)
      
implicit none
      
!-----------------------------------------------------------------------
! m.r. raupach, 24-oct-92
! see: raupach, 1992, blm 60 375-395
!      mrr notes "simplified wind model for canopy", 23-oct-92
!      mrr draft paper "simplified expressions...", dec-92
!-----------------------------------------------------------------------
! inputs:
!   h     = roughness height
!   rlai  = leaf area index (assume rl = frontal area index = rlai/2)
! output:
!   usuh  = us/uh (us=friction velocity, uh = mean velocity at z=h)
!   z0    = roughness length
!   d     = zero-plane displacement
!   coexp = coefficient in exponential in-canopy wind profile
!           u(z) = u(h)*exp(coexp*(z/h-1)), found by gradient-matching
!           canopy and roughness-sublayer u(z) at z=h
!-----------------------------------------------------------------------
real h,rlai,usuh,z0,d,coexp
real psih,rl,usuhl,xx,dh,z0h
! preset parameters:
real cr,cs,beta,ccd,ccw,usuhm,vonk
parameter (cr    = 0.3)          ! element drag coefficient
parameter (cs    = 0.003)        ! substrate drag coefficient
parameter (beta  = cr/cs)        ! ratio cr/cs
parameter (ccd   = 15.0)         ! constant in d/h equation
parameter (ccw   = 2.0)          ! ccw=(zw-d)/(h-d)
parameter (usuhm = 0.3)          ! (max of us/uh)
parameter (vonk  = 0.4)          ! von karman constant
psih=alog(ccw)-1.0+1.0/ccw  ! i.e. .19315
rl = rlai*0.5
! find uh/us
usuhl  = sqrt(cs+cr*rl)            ! sqrt(.003 + .15*rlai)
usuh   = min(usuhl,usuhm)
! find d/h and d 
xx     = sqrt(ccd*max(rl,0.0005))  ! sqrt(7.5*rlai)
dh     = 1.0 - (1.0 - exp(-xx))/xx ! .5*xx -.166*xx*xx + .

d      = dh*h                      ! not used
! find z0h and z0:
z0h    = (1.0 - dh) * exp(psih - vonk/usuh)
z0     = z0h*h 
! find coexp: see notes "simplified wind model ..." eq 34a
coexp  = usuh / (vonk*ccw*(1.0 - dh))
return ! ruff
end subroutine cruf2
!=======================================================================
end module indata