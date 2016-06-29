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

! CABLE interface originally developed by the CABLE group
! Subsequently modified by MJT for 5 tile mosaic and SEAESF radiation scheme
  
! - Currently all tiles have the same soil texture, but independent soil temperatures,
!   moisture, etc.
! - LAI can be interpolated between timesteps using a PWCB fit to the LAI integral
!   or LAI can be taken as constant for the month
! - CO2 can be constant or read from the radiation code.  A tracer CO2 is avaliable
!   when tracers are active
! - The code assumes only one month at a time is integrated in RCM mode.  However,
!   since the input files can be modified at runtime (not all months are loaded
!   at once), then we can support off-line evolving/dynamic vegetation, etc.

! The following mappings between IGBP and CSIRO PFT were recommended by RL
    
! ivegt   IGBP type                             CSIRO PFT
! 1       Evergreen Needleleaf Forest           1.  Evergreen Needleleaf
! 2       Evergreen Broadleaf Forest            1.  Evergreen Broadleaf
! 3       Deciduous Needleaf Forest             1.  Deciduous Needleleaf
! 4       Deciduous Broadleaf Forest            1.  Deciduous Broadleaf
! 5       Mixed Forest                          1.  Deciduous Broadlead                              when -25<lat<25
!                                               0.5 Evergreen Needleleaf    0.5 Deciduous Broadlead  when lat<-25 or lat>25
! 6       Closed Shrublands                     0.8 Shrub                   0.2 (Grass)
! 7       Open Shrublands                       0.2 Shrub                   0.8 (Grass)
! 8       Woody Savannas                        0.6 (Grass)                 0.4 Evergreen Needleleaf when lat<-40 or lat>40
!                                               0.6 (Grass)                 0.4 Evergreen Broadleaf  when -40<lat<40
! 9       Savannas                              0.9 (Grass)                 0.1 Evergreen Needleleaf when lat<-40 or lat>40
!                                               0.9 (Grass)                 0.1 Evergreen Broadleaf  when -40<lat<40
! 10      Grasslands                            1.  (Grass)
! 11      Permanent Wetlands                    1.  Wetland
! 12      Croplands                             1.  (Crop)
! 13      Urban and Built-up                    1.  Urban
! 14      Cropland/Natural Vegetation Mosaic    1.  (Crop)
! 15      Snow and Ice                          1.  Ice
! 16      Barren or Sparsely Vegetated          1.  Barren
! 17      Water Bodies                          1.  Lakes

! where:
!   (Grass)   0.9  C3 0.1  C4 0. Tundra   40<lat<50 or -50<lat<-40
!             0.8  C3 0.2  C4 0. Tundra   30<lat<40 or -40<lat<-30
!             0.5  C3 0.5  C4 0. Tundra   25<lat<30 or -30<lat<-25
!             0.05 C3 0.95 C4 0. Tundra  -25<lat<25
!             0.   C3 0.   C4 1. Tundra   lat<-50 or lat>50

!   (Crop)    0.7 C3  0.3 C4   -30<lat<30
!             0.9 C3  0.1 C4   30<lat<40 or -40<lat<-30
!             1.  C3  0.  C4   lat<-40   or lat>40

! *** NOTE MJT SPECIAL
! The PFT evergreen broadleaf's canopy height is reduced for woody savannas for improved roughness length

! CSIRO PFT index
! 1  Evergreen Needleleaf
! 2  Evergreen Broadleaf
! 3  Deciduous Needleaf
! 4  Deciduous Broadleaf
! 5  Shrub
! 6  C3 grass
! 7  C4 grass
! 8  Tundra
! 9  C3 crop
! 10 C4 crop
! 11 Wetland
! 12 Not used
! 13 Not used
! 14 Barren
! 15 Urban
! 16 Lakes
! 17 Ice
  
! isoilm  type
! 0       water/ocean
! 1       coarse               sand/loamy_sand
! 2       medium               clay-loam/silty-clay-loam/silt-loam
! 3       fine                 clay
! 4       coarse-medium        sandy-loam/loam
! 5       coarse-fine          sandy-clay
! 6       medium-fine          silty-clay 
! 7       coarse-medium-fine   sandy-clay-loam
! 8       organi!              peat
! 9       land ice

module cable_ccam


use cable_common_module
use cable_data_module
use cable_def_types_mod, cbm_ms => ms
use cable_roughness_module
use cable_soil_snow_module
use casa_cnp_module
use casadimension
use casaparm, xroot => froot
use casavariable

implicit none

private
public sib4,loadcbmparm,loadtile,cableinflow
public proglai

! The following options will eventually be moved to the globpe.f namelist
integer, save :: proglai             = 0 ! 0 prescribed LAI, 1 prognostic LAI
integer, parameter :: tracerco2      = 0 ! 0 use radiation CO2, 1 use tracer CO2 
real, parameter :: minfrac = 0.01        ! minimum non-zero tile fraction (improves load balancing)

integer, dimension(5,2), save :: pind  
integer, save :: maxnb
real, dimension(:), allocatable, save :: sv,vl1,vl2,vl3
logical, dimension(:,:), allocatable, save :: tmap
type (air_type), save            :: air
type (bgc_pool_type), save       :: bgc
type (met_type), save            :: met
type (balances_type), save       :: bal
type (radiation_type), save      :: rad
type (roughness_type), save      :: rough
type (soil_parameter_type), save :: soil
type (soil_snow_type), save      :: ssnow
type (sum_flux_type), save       :: sum_flux
type (veg_parameter_type), save  :: veg
type (canopy_type), save         :: canopy
type (casa_balance), save        :: casabal
type (casa_biome), save          :: casabiome
type (casa_flux), save           :: casaflux
type (casa_met), save            :: casamet
type (casa_pool), save           :: casapool
type (phen_variable), save       :: phen
type (physical_constants), save  :: c

contains
! ****************************************************************************

! CABLE-CCAM interface
subroutine sib4

use arrays_m
use cc_mpi
use infile
use latlong_m
use soil_m
use soilsnow_m
use work2_m, only : wetfac

implicit none

include 'newmpar.h'
include 'const_phys.h'
include 'dates.h'
include 'parm.h'

real fjd, dhr
real, dimension(ifull) :: tmps
real(r_2), dimension(mp) :: xksoil
integer jyear,jmonth,jday,jhour,jmin
integer k,mins,nb,iq,j
integer idoy,is,ie

! abort calculation if no land points on this processor  
if (mp<=0) return

! calculate zenith angle
dhr = dt/3600.

!--------------------------------------------------------------
! CABLE
ktau_gl          = 900
kend_gl          = 999
ssnow%owetfac    = ssnow%wetfac
canopy%oldcansto = canopy%cansto

ssnow%otss_0     = ssnow%otss
ssnow%otss       = ssnow%tss
ssnow%owetfac    = ssnow%wetfac
call soil_snow(dt,soil,ssnow,canopy,met,bal,veg)
do k=1,ms
  where ( land )
    tgg(:,k)=0.
    wb(:,k)=0.
    wbice(:,k)=0.
  end where
end do
do k=1,3
  where ( land )
    tggsn(:,k)=0.
    smass(:,k)=0.
    ssdn(:,k)=0.
  end where
end do
where ( land )
  wetfac=0.
end where
tmps=0. ! average isflag

do nb=1,maxnb
  is = pind(nb,1)
  ie = pind(nb,2)
  ! soil
  do k=1,ms
    tgg(:,k)  =tgg(:,k)  +unpack(sv(is:ie)*ssnow%tgg(is:ie,k),        tmap(:,nb),0.)
    wb(:,k)   =wb(:,k)   +unpack(sv(is:ie)*real(ssnow%wb(is:ie,k)),   tmap(:,nb),0.)
    wbice(:,k)=wbice(:,k)+unpack(sv(is:ie)*real(ssnow%wbice(is:ie,k)),tmap(:,nb),0.)
  end do
  ! hydrology
  runoff=runoff+unpack(sv(is:ie)*ssnow%runoff(is:ie)*dt,tmap(:,nb),0.) ! convert mm/s to mm
end do
return
end subroutine sib4





! *************************************************************************************
subroutine loadcbmparm(fveg,fvegprev,fvegnext,fphen,casafile)

!~ use carbpools_m
use cc_mpi
use infile
use latlong_m
use nsibd_m
use pbl_m
use sigs_m
use soil_m
use soilsnow_m
!~ use vegpar_m
  
implicit none
  
include 'newmpar.h'
include 'const_phys.h'
include 'parm.h'
include 'soilv.h'

integer, dimension(ifull,5) :: ivs
integer iq,n,k,ipos,iv,ncount
integer, dimension(1) :: pos
integer jyear,jmonth,jday,jhour,jmin,mins
integer, dimension(1) :: lndtst,lndtst_g
real fc3,fc4,ftu,fg3,fg4,clat,nsum
real fjd,xp
real, dimension(ifull,mxvt,0:2) :: newlai
real, dimension(mxvt,ms) :: froot2
real, dimension(mxvt,ncp) :: tcplant
real, dimension(mxvt,ncs) :: tcsoil
real, dimension(mxvt,mplant) :: ratiocnplant
real, dimension(mxvt,msoil) :: ratiocnsoil,ratiocnsoilmax,ratiocnsoilmin
real, dimension(mxvt,2) :: taul,refl  
real, dimension(ifull,mxvt) :: newgrid
real, dimension(ifull,5) :: svs,vlin,vlinprev,vlinnext
real, dimension(ifull,2) :: albsoilsn
real, dimension(12,msoil) :: rationpsoil
real, dimension(ncp) :: ratecp
real, dimension(ncs) :: ratecs
real, dimension(mxvt) :: canst1,leaf_w,leaf_l,ejmax,hc,rp20
real, dimension(mxvt) :: rpcoef,shelrb,vcmax,xfang
real, dimension(mxvt) :: tminvj,tmaxvj,vbeta
real, dimension(mxvt) :: extkn,rootbeta,vegcf,c4frac
real, dimension(mxvt) :: leafage,woodage,frootage,metage
real, dimension(mxvt) :: strage,cwdage,micage,slowage,passage
real, dimension(mxvt) :: xfherbivore,xxkleafcoldmax,xxkleafdrymax
real, dimension(mxvt) :: xratioNPleafmin,xratioNPleafmax,xratioNPwoodmin,xratioNPwoodmax
real, dimension(mxvt) :: xratioNPfrootmin,xratioNPfrootmax,xfNminloss,xfNminleach,xnfixrate
real, dimension(mxvt) :: xnsoilmin,xplab,xpsorb,xpocc
real, dimension(mxvt) :: cleaf,cwood,cfroot,cmet,cstr,ccwd,cmic,cslow,cpass,nleaf
real, dimension(mxvt) :: nwood,nfroot,nmet,nstr,ncwd,nmic,nslow,npass,xpleaf,xpwood
real, dimension(mxvt) :: xpfroot,xpmet,xpstr,xpcwd,xpmic,xpslow,xppass,clabileage
real, dimension(ifull) :: albsoil, savannafrac
real, dimension(12) :: xkmlabp,xpsorbmax,xfPleach
character(len=*), intent(in) :: fveg,fvegprev,fvegnext,fphen,casafile

WRITE(6,*) 'IN loadcbmparam'
if ( myid==0 ) write(6,*) "Initialising CABLE"

if ( cbm_ms/=ms ) then
  write(6,*) "ERROR: CABLE and CCAM soil levels do not match"
  call ccmpi_abort(-1)
end if

! redefine rhos
rhos=(/ 1600., 1600., 1381., 1373., 1476., 1521., 1373., 1537.,  910., 2600., 2600., 2600., 2600. /)

! biophysical parameter tables
hc    =(/   17.,  35.,  15.5,  20.,   0.6, 0.567, 0.567, 0.567, 0.55, 0.55, 0.567,  0.2, 6.017,  0.2,  0.2,  0.2,  0.2 /)
xfang =(/  0.01,  0.1,  0.01, 0.25,  0.01,  -0.3,  -0.3,  -0.3, -0.3, -0.3,  -0.3,  0.1,    0.,   0.,   0.,   0.,   0. /)
leaf_w=(/ 0.001, 0.05, 0.001, 0.08, 0.005,  0.01,  0.01,  0.01, 0.01, 0.01,  0.01, 0.03, 0.015, 0.00,   0.,   0.,   0. /)
leaf_l=(/ 0.055, 0.10, 0.040, 0.15, 0.100,  0.30,  0.30,  0.30, 0.30, 0.30,  0.30, 0.30, 0.242, 0.03, 0.03, 0.03, 0.03 /)
canst1=0.1
shelrb=2.
extkn=0.001 ! new definition for nitrogen (since CABLE v1.9b)
refl(:,1)=(/ 0.062,0.076,0.056,0.092,0.100,0.110,0.100,0.117,0.100,0.090,0.108,0.055,0.091,0.238,0.143,0.143,0.159 /)
refl(:,2)=(/ 0.302,0.350,0.275,0.380,0.400,0.470,0.400,0.343,0.400,0.360,0.343,0.190,0.310,0.457,0.275,0.275,0.305 /)
taul(:,1)=(/ 0.050,0.050,0.045,0.050,0.050,0.070,0.100,0.080,0.100,0.090,0.075,0.023,0.059,0.039,0.023,0.023,0.026 /)
taul(:,2)=(/ 0.100,0.250,0.144,0.250,0.240,0.250,0.150,0.124,0.150,0.225,0.146,0.198,0.163,0.189,0.113,0.113,0.113 /)
vegcf    =(/    9.,  14.,   9.,   8.,   5.,   7.,   7.,   5.,   7.,   1.,   7.,   1.,   1.,   1.,   1.,   1.,   1. /)
vcmax=(/ 40.E-6,55.E-6,40.E-6,60.E-6,40.E-6,60.E-6,10.E-6,40.E-6,80.E-6,80.E-6,60.E-6,17.E-6,1.E-6,17.E-6,17.E-6,17.E-6,17.E-6 /)
ejmax=2.*vcmax
rp20=(/ 3., 0.6, 3., 2.2, 1., 1.5, 2.8, 2.5, 1.5, 1., 1.5, 1., 1., 1., 1., 1., 1. /)
rpcoef=0.0832
rs20  =(/   1.,   1.,  1.,  1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   0.,   1.,   0.,   0.,   0.,   0. /)
tminvj=(/ -15., -15.,  5.,  5., -15., -15., -15., -15., -15., -15., -15., -15., -15., -15., -15., -15., -15. /)
tmaxvj=(/ -10., -10., 10., 15., -10., -10., -10., -10., -10., -10., -10., -10., -10., -10., -10., -10., -10. /)
vbeta=1.
rootbeta=(/ 0.943,0.962,0.966,0.961,0.964,0.943,0.943,0.943,0.961,0.961,0.943,0.975,0.961,0.961,0.961,0.961,0.961 /)
tcplant(:,1)=(/ 200.  , 300.  , 200. , 300.  , 159. , 250., 250., 250., 150., 150., 250., 1., 0.1, 0., 1., 1., 0. /)
tcplant(:,2)=(/ 10217., 16833., 5967., 12000., 5000., 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0., 0.,  0., 0., 0., 0. /)
tcplant(:,3)=(/ 876.  , 1443. , 511. , 1029. , 500. , 500., 500., 500., 607., 607., 500., 1., 0.1, 0., 1., 1., 0. /)
tcsoil(:,1) =(/ 184.  , 303.  , 107. , 216.  , 100. , 275., 275., 275., 149., 149., 275., 1., 0.1, 1., 1., 1., 1. /)
tcsoil(:,2) =(/ 367.  , 606.  , 214. , 432.  , 250. , 314., 314., 314., 300., 300., 314., 1., 0.1, 1., 1., 1., 1. /)
ratecp(1:3)=(/ 1., 0.03, 0.14 /)
ratecs(1:2)=(/ 2., 0.5 /)
c4frac=(/ 0., 0., 0., 0., 0., 0., 1., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0. /)

! read CABLE biome and LAI data
if ( myid==0 ) then
  write(6,*) "Reading tiled surface data for CABLE"
  call vegta(ivs,svs,vlinprev,vlin,vlinnext,fvegprev,fveg,fvegnext)
else
  call vegtb(ivs,svs,vlinprev,vlin,vlinnext,fvegprev,fveg,fvegnext)
end if
do n=1,5
  svs(:,n)=svs(:,n)/sum(svs,2)
end do

icycle=ccycle
cable_user%fwsoil_switch="standard"

if ( myid==0 ) write(6,*) "Define CABLE and CASA CNP arrays"

! default values (i.e., no land)  
ivegt=0
albsoilsn=0.08  
albsoil=0.08
albvisdir=0.08
albvisdif=0.08
albnirdir=0.08
albnirdif=0.08
zolnd=0.
!~ cplant=0.
!~ csoil=0.
pind=ifull+1
mvtype=mxvt
mstype=mxst

if ( myid==0 ) write(6,*) "Mapping IGBP classes to CSIRO PFTs"
mp=0
newgrid=0.
newlai=0.
savannafrac=0.
do iq=1,ifull
  if ( land(iq) ) then
    clat=rlatt(iq)*180./pi
    ! grass
    if (abs(clat)>50.5) then
      fg3=0.
      fg4=0.
    else if (abs(clat)>49.5) then
      xp=abs(clat)-49.5
      fg3=(1.-xp)*0.9
      fg4=(1.-xp)*0.1
    else if (abs(clat)>40.5) then
      fg3=0.9
      fg4=0.1
    else if (abs(clat)>39.5) then
      xp=abs(clat)-39.5
      fg3=(1.-xp)*0.8+xp*0.9
      fg4=(1.-xp)*0.2+xp*0.1
    else if (abs(clat)>30.5) then
      fg3=0.8
      fg4=0.2
    else if (abs(clat)>29.5) then
      xp=abs(clat)-29.5
      fg3=(1.-xp)*0.5+xp*0.8
      fg4=(1.-xp)*0.5+xp*0.2
    else if (abs(clat)>25.5) then
      fg3=0.5
      fg4=0.5
    else if (abs(clat)>24.5) then
      xp=abs(clat)-24.5
      fg3=(1.-xp)*0.05+xp*0.5
      fg4=(1.-xp)*0.95+xp*0.5
    else
      fg3=0.05
      fg4=0.95
    end if
    ftu=1.-fg3-fg4
    ! crops
    if (abs(clat)>40.5) then
      fc3=1.
    else if (abs(clat)>39.5) then
      xp=abs(clat)-39.5
      fc3=(1.-xp)*0.9+xp
    else if (abs(clat)>30.5) then
      fc3=0.9
    else if (abs(clat)>29.5) then
      xp=abs(clat)-29.5
      fc3=(1.-xp)*0.7+xp*0.9
    else
      fc3=0.7
    end if
    fc4=1.-fc3
    do n=1,5
      select case (ivs(iq,n))
        case (1,2,3,4,11)
          newgrid(iq,ivs(iq,n))=newgrid(iq,ivs(iq,n))+svs(iq,n)
          newlai(iq,ivs(iq,n),0)=newlai(iq,ivs(iq,n),0)+svs(iq,n)*vlinprev(iq,n)
          newlai(iq,ivs(iq,n),1)=newlai(iq,ivs(iq,n),1)+svs(iq,n)*vlin(iq,n)
          newlai(iq,ivs(iq,n),2)=newlai(iq,ivs(iq,n),2)+svs(iq,n)*vlinnext(iq,n)
        case (5)
          if (abs(clat)>25.5) then
            newgrid(iq,1)=newgrid(iq,1)+svs(iq,n)*0.5
            newlai(iq,1,0)=newlai(iq,1,0)+svs(iq,n)*0.5*vlinprev(iq,n)
            newlai(iq,1,1)=newlai(iq,1,1)+svs(iq,n)*0.5*vlin(iq,n)
            newlai(iq,1,2)=newlai(iq,1,2)+svs(iq,n)*0.5*vlinnext(iq,n)
            newgrid(iq,4)=newgrid(iq,4)+svs(iq,n)*0.5
            newlai(iq,4,0)=newlai(iq,4,0)+svs(iq,n)*0.5*vlinprev(iq,n)
            newlai(iq,4,1)=newlai(iq,4,1)+svs(iq,n)*0.5*vlin(iq,n)
            newlai(iq,4,2)=newlai(iq,4,2)+svs(iq,n)*0.5*vlinnext(iq,n)
          else if (abs(clat)>24.5) then
            xp=abs(clat)-24.5
            newgrid(iq,1)=newgrid(iq,1)+svs(iq,n)*0.5*xp
            newlai(iq,1,0)=newlai(iq,1,0)+svs(iq,n)*0.5*vlinprev(iq,n)*xp
            newlai(iq,1,1)=newlai(iq,1,1)+svs(iq,n)*0.5*vlin(iq,n)*xp
            newlai(iq,1,2)=newlai(iq,1,2)+svs(iq,n)*0.5*vlinnext(iq,n)*xp
            newgrid(iq,4)=newgrid(iq,4)+svs(iq,n)*(1.-0.5*xp)
            newlai(iq,4,0)=newlai(iq,4,0)+svs(iq,n)*vlinprev(iq,n)*(1.-0.5*xp)
            newlai(iq,4,1)=newlai(iq,4,1)+svs(iq,n)*vlin(iq,n)*(1.-0.5*xp)
            newlai(iq,4,2)=newlai(iq,4,2)+svs(iq,n)*vlinnext(iq,n)*(1.-0.5*xp)
          else
            newgrid(iq,4)=newgrid(iq,4)+svs(iq,n)
            newlai(iq,4,0)=newlai(iq,4,0)+svs(iq,n)*vlinprev(iq,n)
            newlai(iq,4,1)=newlai(iq,4,1)+svs(iq,n)*vlin(iq,n)
            newlai(iq,4,2)=newlai(iq,4,2)+svs(iq,n)*vlinnext(iq,n)
          end if
        case (6)
          newgrid(iq,5)=newgrid(iq,5)+svs(iq,n)*0.8
          newlai(iq,5,0)=newlai(iq,5,0)+svs(iq,n)*0.8*vlinprev(iq,n)
          newlai(iq,5,1)=newlai(iq,5,1)+svs(iq,n)*0.8*vlin(iq,n)
          newlai(iq,5,2)=newlai(iq,5,2)+svs(iq,n)*0.8*vlinnext(iq,n)
          newgrid(iq,6)=newgrid(iq,6)+svs(iq,n)*0.2*fg3
          newlai(iq,6,0)=newlai(iq,6,0)+svs(iq,n)*0.2*fg3*vlinprev(iq,n)
          newlai(iq,6,1)=newlai(iq,6,1)+svs(iq,n)*0.2*fg3*vlin(iq,n)
          newlai(iq,6,2)=newlai(iq,6,2)+svs(iq,n)*0.2*fg3*vlinnext(iq,n)
          newgrid(iq,7)=newgrid(iq,7)+svs(iq,n)*0.2*fg4
          newlai(iq,7,0)=newlai(iq,7,0)+svs(iq,n)*0.2*fg4*vlinprev(iq,n)
          newlai(iq,7,1)=newlai(iq,7,1)+svs(iq,n)*0.2*fg4*vlin(iq,n)
          newlai(iq,7,2)=newlai(iq,7,2)+svs(iq,n)*0.2*fg4*vlinnext(iq,n)
          newgrid(iq,8)=newgrid(iq,8)+svs(iq,n)*0.2*ftu
          newlai(iq,8,0)=newlai(iq,8,0)+svs(iq,n)*0.2*ftu*vlinprev(iq,n)
          newlai(iq,8,1)=newlai(iq,8,1)+svs(iq,n)*0.2*ftu*vlin(iq,n)
          newlai(iq,8,2)=newlai(iq,8,2)+svs(iq,n)*0.2*ftu*vlinnext(iq,n)
        case (7)
          newgrid(iq,5)=newgrid(iq,5)+svs(iq,n)*0.2
          newlai(iq,5,0)=newlai(iq,5,0)+svs(iq,n)*0.2*vlinprev(iq,n)
          newlai(iq,5,1)=newlai(iq,5,1)+svs(iq,n)*0.2*vlin(iq,n)
          newlai(iq,5,2)=newlai(iq,5,2)+svs(iq,n)*0.2*vlinnext(iq,n)
          newgrid(iq,6)=newgrid(iq,6)+svs(iq,n)*0.8*fg3
          newlai(iq,6,0)=newlai(iq,6,0)+svs(iq,n)*0.8*fg3*vlinprev(iq,n)
          newlai(iq,6,1)=newlai(iq,6,1)+svs(iq,n)*0.8*fg3*vlin(iq,n)
          newlai(iq,6,2)=newlai(iq,6,2)+svs(iq,n)*0.8*fg3*vlinnext(iq,n)
          newgrid(iq,7)=newgrid(iq,7)+svs(iq,n)*0.8*fg4
          newlai(iq,7,0)=newlai(iq,7,0)+svs(iq,n)*0.8*fg4*vlinprev(iq,n)
          newlai(iq,7,1)=newlai(iq,7,1)+svs(iq,n)*0.8*fg4*vlin(iq,n)
          newlai(iq,7,2)=newlai(iq,7,2)+svs(iq,n)*0.8*fg4*vlinnext(iq,n)
          newgrid(iq,8)=newgrid(iq,8)+svs(iq,n)*0.8*ftu
          newlai(iq,8,0)=newlai(iq,8,0)+svs(iq,n)*0.8*ftu*vlinprev(iq,n)
          newlai(iq,8,1)=newlai(iq,8,1)+svs(iq,n)*0.8*ftu*vlin(iq,n)
          newlai(iq,8,2)=newlai(iq,8,2)+svs(iq,n)*0.8*ftu*vlinnext(iq,n)
        case (8)
          if (abs(clat)>40.5) then
            newgrid(iq,1)=newgrid(iq,1)+svs(iq,n)*0.4
            newlai(iq,1,0)=newlai(iq,1,0)+svs(iq,n)*0.4*vlinprev(iq,n)
            newlai(iq,1,1)=newlai(iq,1,1)+svs(iq,n)*0.4*vlin(iq,n)
            newlai(iq,1,2)=newlai(iq,1,2)+svs(iq,n)*0.4*vlinnext(iq,n)
          else if (abs(clat)>39.5) then
            xp=abs(clat)-39.5
            newgrid(iq,1)=newgrid(iq,1)+svs(iq,n)*0.4*xp
            newlai(iq,1,0)=newlai(iq,1,0)+svs(iq,n)*vlinprev(iq,n)*0.4*xp
            newlai(iq,1,1)=newlai(iq,1,1)+svs(iq,n)*vlin(iq,n)*0.4*xp
            newlai(iq,1,2)=newlai(iq,1,2)+svs(iq,n)*vlinnext(iq,n)*0.4*xp   
            savannafrac(iq)=savannafrac(iq)+svs(iq,n)*0.4*(1.-xp)
            newgrid(iq,2)=newgrid(iq,2)+svs(iq,n)*0.4*(1.-xp)
            newlai(iq,2,0)=newlai(iq,2,0)+svs(iq,n)*vlinprev(iq,n)*0.4*(1.-xp)
            newlai(iq,2,1)=newlai(iq,2,1)+svs(iq,n)*vlin(iq,n)*0.4*(1.-xp)
            newlai(iq,2,2)=newlai(iq,2,2)+svs(iq,n)*vlinnext(iq,n)*0.4*(1.-xp)
          else
            savannafrac(iq)=savannafrac(iq)+svs(iq,n)*0.4
            newgrid(iq,2)=newgrid(iq,2)+svs(iq,n)*0.4
            newlai(iq,2,0)=newlai(iq,2,0)+svs(iq,n)*0.4*vlinprev(iq,n)
            newlai(iq,2,1)=newlai(iq,2,1)+svs(iq,n)*0.4*vlin(iq,n)
            newlai(iq,2,2)=newlai(iq,2,2)+svs(iq,n)*0.4*vlinnext(iq,n)
          end if
          newgrid(iq,6)=newgrid(iq,6)+svs(iq,n)*0.6*fg3
          newlai(iq,6,0)=newlai(iq,6,0)+svs(iq,n)*0.6*fg3*vlinprev(iq,n)
          newlai(iq,6,1)=newlai(iq,6,1)+svs(iq,n)*0.6*fg3*vlin(iq,n)
          newlai(iq,6,2)=newlai(iq,6,2)+svs(iq,n)*0.6*fg3*vlinnext(iq,n)
          newgrid(iq,7)=newgrid(iq,7)+svs(iq,n)*0.6*fg4
          newlai(iq,7,0)=newlai(iq,7,0)+svs(iq,n)*0.6*fg4*vlinprev(iq,n)
          newlai(iq,7,1)=newlai(iq,7,1)+svs(iq,n)*0.6*fg4*vlin(iq,n)
          newlai(iq,7,2)=newlai(iq,7,2)+svs(iq,n)*0.6*fg4*vlinnext(iq,n)
          newgrid(iq,8)=newgrid(iq,8)+svs(iq,n)*0.6*ftu
          newlai(iq,8,0)=newlai(iq,8,0)+svs(iq,n)*0.6*ftu*vlinprev(iq,n)
          newlai(iq,8,1)=newlai(iq,8,1)+svs(iq,n)*0.6*ftu*vlin(iq,n)
          newlai(iq,8,2)=newlai(iq,8,2)+svs(iq,n)*0.6*ftu*vlinnext(iq,n)
        case (9)
          if (abs(clat)>40.5) then
            newgrid(iq,1)=newgrid(iq,1)+svs(iq,n)*0.1
            newlai(iq,1,0)=newlai(iq,1,0)+svs(iq,n)*0.1*vlinprev(iq,n)
            newlai(iq,1,1)=newlai(iq,1,1)+svs(iq,n)*0.1*vlin(iq,n)
            newlai(iq,1,2)=newlai(iq,1,2)+svs(iq,n)*0.1*vlinnext(iq,n)
          else if (abs(clat)>39.5) then
            xp=abs(clat)-39.5
            newgrid(iq,1)=newgrid(iq,1)+svs(iq,n)*0.1*xp
            newlai(iq,1,0)=newlai(iq,1,0)+svs(iq,n)*vlinprev(iq,n)*0.1*xp
            newlai(iq,1,1)=newlai(iq,1,1)+svs(iq,n)*vlin(iq,n)*0.1*xp
            newlai(iq,1,2)=newlai(iq,1,2)+svs(iq,n)*vlinnext(iq,n)*0.1*xp
            savannafrac(iq)=savannafrac(iq)+svs(iq,n)*0.1*(1.-xp)
            newgrid(iq,2)=newgrid(iq,2)+svs(iq,n)*0.1*(1.-xp)
            newlai(iq,2,0)=newlai(iq,2,0)+svs(iq,n)*vlinprev(iq,n)*0.1*(1.-xp)
            newlai(iq,2,1)=newlai(iq,2,1)+svs(iq,n)*vlin(iq,n)*0.1*(1.-xp)
            newlai(iq,2,2)=newlai(iq,2,2)+svs(iq,n)*vlinnext(iq,n)*0.1*(1.-xp)
          else
            savannafrac(iq)=savannafrac(iq)+svs(iq,n)*0.1
            newgrid(iq,2)=newgrid(iq,2)+svs(iq,n)*0.1
            newlai(iq,2,0)=newlai(iq,2,0)+svs(iq,n)*0.1*vlinprev(iq,n)
            newlai(iq,2,1)=newlai(iq,2,1)+svs(iq,n)*0.1*vlin(iq,n)
            newlai(iq,2,2)=newlai(iq,2,2)+svs(iq,n)*0.1*vlinnext(iq,n)
          end if
          newgrid(iq,6)=newgrid(iq,6)+svs(iq,n)*0.9*fg3
          newlai(iq,6,0)=newlai(iq,6,0)+svs(iq,n)*0.9*fg3*vlinprev(iq,n)
          newlai(iq,6,1)=newlai(iq,6,1)+svs(iq,n)*0.9*fg3*vlin(iq,n)
          newlai(iq,6,2)=newlai(iq,6,2)+svs(iq,n)*0.9*fg3*vlinnext(iq,n)
          newgrid(iq,7)=newgrid(iq,7)+svs(iq,n)*0.9*fg4
          newlai(iq,7,0)=newlai(iq,7,0)+svs(iq,n)*0.9*fg4*vlinprev(iq,n)
          newlai(iq,7,1)=newlai(iq,7,1)+svs(iq,n)*0.9*fg4*vlin(iq,n)
          newlai(iq,7,2)=newlai(iq,7,2)+svs(iq,n)*0.9*fg4*vlinnext(iq,n)
          newgrid(iq,8)=newgrid(iq,8)+svs(iq,n)*0.9*ftu
          newlai(iq,8,0)=newlai(iq,8,0)+svs(iq,n)*0.9*ftu*vlinprev(iq,n)
          newlai(iq,8,1)=newlai(iq,8,1)+svs(iq,n)*0.9*ftu*vlin(iq,n)
          newlai(iq,8,2)=newlai(iq,8,2)+svs(iq,n)*0.9*ftu*vlinnext(iq,n)
        case (10)
          newgrid(iq,6)=newgrid(iq,6)+svs(iq,n)*fg3
          newlai(iq,6,0)=newlai(iq,6,0)+svs(iq,n)*fg3*vlinprev(iq,n)
          newlai(iq,6,1)=newlai(iq,6,1)+svs(iq,n)*fg3*vlin(iq,n)
          newlai(iq,6,2)=newlai(iq,6,2)+svs(iq,n)*fg3*vlinnext(iq,n)
          newgrid(iq,7)=newgrid(iq,7)+svs(iq,n)*fg4
          newlai(iq,7,0)=newlai(iq,7,0)+svs(iq,n)*fg4*vlinprev(iq,n)
          newlai(iq,7,1)=newlai(iq,7,1)+svs(iq,n)*fg4*vlin(iq,n)
          newlai(iq,7,2)=newlai(iq,7,2)+svs(iq,n)*fg4*vlinnext(iq,n)
          newgrid(iq,8)=newgrid(iq,8)+svs(iq,n)*ftu
          newlai(iq,8,0)=newlai(iq,8,0)+svs(iq,n)*ftu*vlinprev(iq,n)
          newlai(iq,8,1)=newlai(iq,8,1)+svs(iq,n)*ftu*vlin(iq,n)
          newlai(iq,8,2)=newlai(iq,8,2)+svs(iq,n)*ftu*vlinnext(iq,n)
        case (12,14)
          newgrid(iq,9)=newgrid(iq,9)+svs(iq,n)*fc3
          newlai(iq,9,0)=newlai(iq,9,0)+svs(iq,n)*fc3*vlinprev(iq,n)
          newlai(iq,9,1)=newlai(iq,9,1)+svs(iq,n)*fc3*vlin(iq,n)
          newlai(iq,9,2)=newlai(iq,9,2)+svs(iq,n)*fc3*vlinnext(iq,n)
          newgrid(iq,10)=newgrid(iq,10)+svs(iq,n)*fc4
          newlai(iq,10,0)=newlai(iq,10,0)+svs(iq,n)*fc4*vlinprev(iq,n)
          newlai(iq,10,1)=newlai(iq,10,1)+svs(iq,n)*fc4*vlin(iq,n)
          newlai(iq,10,2)=newlai(iq,10,2)+svs(iq,n)*fc4*vlinnext(iq,n)
        case (13)
          newgrid(iq,15)=newgrid(iq,15)+svs(iq,n)
          newlai(iq,15,0)=newlai(iq,15,0)+svs(iq,n)*vlinprev(iq,n)
          newlai(iq,15,1)=newlai(iq,15,1)+svs(iq,n)*vlin(iq,n)
          newlai(iq,15,2)=newlai(iq,15,2)+svs(iq,n)*vlinnext(iq,n)
        case (15)
          newgrid(iq,17)=newgrid(iq,17)+svs(iq,n)
          newlai(iq,17,0)=newlai(iq,17,0)+svs(iq,n)*vlinprev(iq,n)
          newlai(iq,17,1)=newlai(iq,17,1)+svs(iq,n)*vlin(iq,n)
          newlai(iq,17,2)=newlai(iq,17,2)+svs(iq,n)*vlinnext(iq,n)
        case (16)
          newgrid(iq,14)=newgrid(iq,14)+svs(iq,n)
          newlai(iq,14,0)=newlai(iq,14,0)+svs(iq,n)*vlinprev(iq,n)
          newlai(iq,14,1)=newlai(iq,14,1)+svs(iq,n)*vlin(iq,n)
          newlai(iq,14,2)=newlai(iq,14,2)+svs(iq,n)*vlinnext(iq,n)
        case (17)
          newgrid(iq,16)=newgrid(iq,16)+svs(iq,n)
          newlai(iq,16,0)=newlai(iq,16,0)+svs(iq,n)*vlinprev(iq,n)
          newlai(iq,16,1)=newlai(iq,16,1)+svs(iq,n)*vlin(iq,n)
          newlai(iq,16,2)=newlai(iq,16,2)+svs(iq,n)*vlinnext(iq,n)
        case DEFAULT
          write(6,*) "ERROR: Land-type/lsmask mismatch at myid,iq,ivs,land=",myid,iq,ivs(iq,n),land(iq)
          call ccmpi_abort(-1)
      end select
    end do
    if (newgrid(iq,2)>0.) then
      savannafrac(iq)=savannafrac(iq)/newgrid(iq,2)
    end if
    where (newgrid(iq,:)>0.)
      newlai(iq,:,0)=newlai(iq,:,0)/newgrid(iq,:)
      newlai(iq,:,1)=newlai(iq,:,1)/newgrid(iq,:)
      newlai(iq,:,2)=newlai(iq,:,2)/newgrid(iq,:)
    end where
    ipos=count(newgrid(iq,:)>0.)
    do while (ipos>5)
      pos=minloc(newgrid(iq,:),newgrid(iq,:)>0.)
      newgrid(iq,pos(1))=0.
      nsum=sum(newgrid(iq,:))
      newgrid(iq,:)=newgrid(iq,:)/nsum
      ipos=count(newgrid(iq,:)>0.)
    end do    
    do while (any(newgrid(iq,:)<minfrac.and.newgrid(iq,:)>0.))
      pos=minloc(newgrid(iq,:),newgrid(iq,:)>0.)
      newgrid(iq,pos(1))=0.
      nsum=sum(newgrid(iq,:))
      newgrid(iq,:)=newgrid(iq,:)/nsum
    end do
    ipos=count(newgrid(iq,:)>0.) 
    mp=mp+ipos
  end if
end do

if (nmaxpr==1) then
  write(6,*) "myid,landtile ",myid,mp

  lndtst(1)=0
  if (mp>0) lndtst(1)=1
  call ccmpi_reduce(lndtst(1:1),lndtst_g(1:1),"sum",0,comm_world)
  if (myid==0) then
    write(6,*) "Processors with land ",lndtst_g(1),nproc
  end if
end if

! if CABLE is present on this processor, then start allocating arrays
! Write messages here in case myid==0 has no land-points (mp==0)
if (myid==0) then
  write(6,*) "Allocating CABLE and CASA CNP arrays"
  if (icycle==0) then
    write(6,*) "Using CABLE carbon cycle"
  else
    write(6,*) "Using CASA CNP"
  end if
end if

if (mp>0) then
  
  allocate(sv(mp))
  allocate(vl1(mp),vl2(mp),vl3(mp))
  allocate(tmap(ifull,5))
  call alloc_cbm_var(air, mp)
  call alloc_cbm_var(bgc, mp)
  call alloc_cbm_var(canopy, mp)
  call alloc_cbm_var(met, mp)
  call alloc_cbm_var(bal, mp)
  call alloc_cbm_var(rad, mp)
  call alloc_cbm_var(rough, mp)
  call alloc_cbm_var(soil, mp)
  call alloc_cbm_var(ssnow, mp)
  call alloc_cbm_var(sum_flux, mp)
  call alloc_cbm_var(veg, mp)

  ! Cable configuration
  cable_user%ssnow_POTEV = ""
  knode_gl = myid
  kwidth_gl = nint(dt) ! MJT notes - what happens when the timestep is less than a second?
  if (kwidth_gl == 0) then
    write(6,*) "ERROR: Timestep too small for CABLE"
    call ccmpi_abort(-1)
  end if
  
  ! soil parameters
  soil%zse        = zse ! soil layer thickness
  soil%zshh(1)    = 0.5 * soil%zse(1)
  soil%zshh(ms+1) = 0.5 * soil%zse(ms)
  soil%zshh(2:ms) = 0.5 * (soil%zse(1:ms-1) + soil%zse(2:ms))
  
  ! froot is now calculated from soil depth and the new parameter rootbeta 
  ! according to Jackson et al. 1996, Oceologica, 108:389-411
  !totdepth = 0.
  !do k=1,ms
  !  totdepth = totdepth + soil%zse(k)*100.
  !  froot2(:,k) = min(1.,1.-rootbeta(:)**totdepth)
  !enddo
  !do k = ms-1, 2, -1
  !  froot2(:,k) = froot2(:,k) - froot2(:,k-1)
  !enddo
  !froot2(:,ms)=1.-sum(froot(:,1:ms-1),2)
  
  ! Eva's method for ACCESS1.3
  froot2(:,1)=0.05
  froot2(:,2)=0.20
  froot2(:,3)=0.20
  froot2(:,4)=0.20
  froot2(:,5)=0.20
  froot2(:,6)=0.15
 
  sv=0.
  vl1=0.
  vl2=0.
  vl3=0.
  tmap=.false.

  ! pack biome data into CABLE vector
  ! prepare LAI arrays for temporal interpolation (PWCB)  
  ! now up to 5 PFT tiles from 5 IGBP classes (need correct order for vectorisation)
  ipos=0
  do n=1,5
    pind(n,1)=ipos+1
    do iq=1,ifull
      if (land(iq)) then
        ncount=0
        do iv=1,mxvt
          if (newgrid(iq,iv)>0.) then
            ncount=ncount+1
            if (ncount==n) exit
          end if
        end do
        if (ncount==n) then
          ipos=ipos+1
          tmap(iq,n)=.true.
          sv(ipos)=newgrid(iq,iv)
          veg%iveg(ipos)=iv
          soil%isoilm(ipos)=isoilm(iq)
          newlai(iq,iv,:)=max(newlai(iq,iv,:),0.01)
          if (fvegprev/=' '.and.fvegnext/=' ') then
            newlai(iq,iv,1)=newlai(iq,iv,1)+newlai(iq,iv,0)
            newlai(iq,iv,2)=newlai(iq,iv,2)+newlai(iq,iv,1)
            vl1(ipos)=0.5*newlai(iq,iv,1)
            vl2(ipos)=4.*newlai(iq,iv,1)-5.*newlai(iq,iv,0)-newlai(iq,iv,2)
            vl3(ipos)=1.5*(newlai(iq,iv,2)+3.*newlai(iq,iv,0)-3.*newlai(iq,iv,1))
          else
            vl1(ipos)=newlai(iq,iv,1)
            vl2(ipos)=0.
            vl3(ipos)=0.
          end if
          if (veg%iveg(ipos)>=14.and.veg%iveg(ipos)<=17) then
            vl1(ipos)=1.E-8
            vl2(ipos)=0.
            vl3(ipos)=0.
          end if
        end if
      end if
    end do
    pind(n,2)=ipos
  end do
  
  if (ipos/=mp) then
    write(6,*) "ERROR: Internal memory allocation error for CABLE set-up"
    call ccmpi_abort(-1)
  end if

  ! Load CABLE arrays
  ivegt=ivs(:,1) ! diagnostic (usually IGBP, not CSIRO pft)
  veg%meth      = 1
  veg%hc        = hc(veg%iveg)
  veg%canst1    = canst1(veg%iveg)
  veg%ejmax     = ejmax(veg%iveg)
  veg%tminvj    = tminvj(veg%iveg)
  veg%tmaxvj    = tmaxvj(veg%iveg)
  veg%vbeta     = vbeta(veg%iveg)
  veg%rp20      = rp20(veg%iveg)
  veg%rpcoef    = rpcoef(veg%iveg)
  veg%shelrb    = shelrb(veg%iveg)
  veg%vcmax     = vcmax(veg%iveg)
  veg%xfang     = xfang(veg%iveg)
  veg%dleaf     = sqrt(max(leaf_w(veg%iveg)*leaf_l(veg%iveg),1.e-20))
  veg%xalbnir   = 1. ! not used
  veg%taul(:,1) = taul(veg%iveg,1)
  veg%taul(:,2) = taul(veg%iveg,2)  
  veg%refl(:,1) = refl(veg%iveg,1)
  veg%refl(:,2) = refl(veg%iveg,2)  
  veg%extkn     = extkn(veg%iveg)
  veg%rs20      = rs20(veg%iveg)
  veg%vegcf     = vegcf(veg%iveg)
  veg%frac4     = c4frac(veg%iveg)
  do k=1,ms
    veg%froot(:,k)=froot2(veg%iveg,k)
  end do

  ! calculate max tile number
  do n=1,5
    if (pind(n,1)<=mp) then
      maxnb=n
    end if
  end do
  
  ! MJT special case for (woody) savannas
  do n=1,maxnb
    where (veg%iveg(pind(n,1):pind(n,2))==2)
      veg%hc(pind(n,1):pind(n,2))=veg%hc(pind(n,1):pind(n,2))+pack((17.-hc(2))*savannafrac,tmap(:,n))
    end where
  end do
  
  ! Calculate LAI and veg fraction diagnostics
  call getzinp(fjd,jyear,jmonth,jday,jhour,jmin,mins)
  !~ call setlai(sigmf,jyear,jmonth,jday,jhour,jmin)
  !~ vlai=0.
  !~ do n=1,maxnb
    !~ vlai=vlai+unpack(sv(pind(n,1):pind(n,2))*veg%vlai(pind(n,1):pind(n,2)),tmap(:,n),0.)
  !~ end do
  
  ! Load CABLE soil data
  soil%bch     = bch(soil%isoilm)
  soil%css     = css(soil%isoilm)
  soil%rhosoil = rhos(soil%isoilm)
  soil%cnsd    = cnsd(soil%isoilm)
  soil%hyds    = hyds(soil%isoilm)
  soil%sucs    = sucs(soil%isoilm)
  soil%hsbh    = hsbh(soil%isoilm)
  soil%sfc     = sfc(soil%isoilm)
  soil%ssat    = ssat(soil%isoilm)
  soil%swilt   = swilt(soil%isoilm)
  soil%ibp2    = ibp2(soil%isoilm)
  soil%i2bp3   = i2bp3(soil%isoilm)
  soil%pwb_min = (soil%swilt/soil%ssat)**soil%ibp2
  soil%clay    = clay(soil%isoilm)
  soil%sand    = sand(soil%isoilm)
  soil%silt    = silt(soil%isoilm)
  bgc%ratecp(:) = ratecp(:)
  bgc%ratecs(:) = ratecs(:)

  !~ ! store bare soil albedo and define snow free albedo
  !~ do n=1,maxnb
    !~ soil%albsoil(pind(n,1):pind(n,2),1)=pack(albvisnir(:,1),tmap(:,n))
    !~ soil%albsoil(pind(n,1):pind(n,2),2)=pack(albvisnir(:,2),tmap(:,n))
  !~ end do
  !~ soil%albsoil(:,3)=0.05
    
  !~ where (land)
    !~ albsoil(:)=0.5*sum(albvisnir,2)
  !~ end where
  !~ where (albsoil<=0.14.and.land)
    !sfact=0.5 for alb <= 0.14
    !~ albsoilsn(:,1)=(1.00/1.50)*albsoil(:)
    !~ albsoilsn(:,2)=(2.00/1.50)*albsoil(:)
  !~ elsewhere ((albsoil(:)<=0.2).and.land)
    !sfact=0.62 for 0.14 < alb <= 0.20
    !~ albsoilsn(:,1)=(1.24/1.62)*albsoil(:)
    !~ albsoilsn(:,2)=(2.00/1.62)*albsoil(:)
  !~ elsewhere (land)
    !sfact=0.68 for 0.2 < alb
    !~ albsoilsn(:,1)=(1.36/1.68)*albsoil(:)
    !~ albsoilsn(:,2)=(2.00/1.68)*albsoil(:)
  !~ end where
  ! MJT suggestion to get an approx inital albedo (before cable is called)
  !~ where (land)
    !~ albvisnir(:,1)=albsoilsn(:,1)*(1.-sigmf)+0.03*sigmf
    !~ albvisnir(:,2)=albsoilsn(:,2)*(1.-sigmf)+0.20*sigmf
  !~ end where
  !~ albvisdir=albvisnir(:,1) ! To be updated by CABLE
  !~ albvisdif=albvisnir(:,1) ! To be updated by CABLE
  !~ albnirdir=albvisnir(:,2) ! To be updated by CABLE
  !~ albnirdif=albvisnir(:,2) ! To be updated by CABLE

  do n=1,maxnb
    ! MJT patch
    soil%albsoil(pind(n,1):pind(n,2),1)   =pack(albsoil,       tmap(:,n))
    soil%albsoil(pind(n,1):pind(n,2),2)   =pack(albsoil,       tmap(:,n))
    ssnow%albsoilsn(pind(n,1):pind(n,2),1)=pack(albsoilsn(:,1),tmap(:,n)) ! overwritten by CABLE
    ssnow%albsoilsn(pind(n,1):pind(n,2),2)=pack(albsoilsn(:,2),tmap(:,n)) ! overwritten by CABLE
    !~ rad%albedo_T(pind(n,1):pind(n,2))     =pack(albsoil,       tmap(:,n))
    !~ rad%trad(pind(n,1):pind(n,2))         =pack(tss,           tmap(:,n))
    !~ rad%latitude(pind(n,1):pind(n,2))     =pack(rlatt,         tmap(:,n))*180./pi
    !~ rad%longitude(pind(n,1):pind(n,2))    =pack(rlongg,        tmap(:,n))*180./pi
  end do
    
  ssnow%albsoilsn(:,3)=0.05    
  ssnow%t_snwlr=0.05
  ssnow%pudsmx=0.

  !~ canopy%oldcansto=0.  
  !~ canopy%ghflux=0.
  !~ canopy%sghflux=0.
  !~ canopy%ga=0.
  !~ canopy%dgdtg=0.
  !~ canopy%fhs_cor=0.
  !~ canopy%fes_cor=0.
  !~ canopy%ga=0.
  !~ canopy%us=0.01
  ssnow%wb_lake=0. ! not used when mlo.f90 is active
  ssnow%fland=1.
  ssnow%ifland=soil%isoilm
    
  ! Initialise sum flux variables
  sum_flux%sumpn=0.
  sum_flux%sumrp=0.
  sum_flux%sumrs=0.
  sum_flux%sumrd=0.
  sum_flux%sumrpw=0.
  sum_flux%sumrpr=0.
  sum_flux%dsumpn=0.
  sum_flux%dsumrp=0.
  sum_flux%dsumrs=0.
  sum_flux%dsumrd=0.
  
  bal%evap_tot=0.
  bal%precip_tot=0.
  bal%ebal_tot=0.
  bal%rnoff_tot=0.
end if
  
if (myid==0) write(6,*) "Finished defining CABLE and CASA CNP arrays"

return
end subroutine loadcbmparm

! *************************************************************************************
! Load CABLE biome and LAI data
! vegta is for myid==0
subroutine vegta(ivs,svs,vlinprev,vlin,vlinnext,fvegprev,fveg,fvegnext)
  
use cc_mpi
use infile

implicit none
  
include 'newmpar.h'
include 'darcdf.h'
include 'parmgeom.h'  ! rlong0,rlat0,schmidt  
  
character(len=*), intent(in) :: fveg,fvegprev,fvegnext
integer, dimension(ifull,5), intent(out) :: ivs
integer, dimension(ifull_g,5) :: ivsg  
integer, dimension(3) :: spos,npos
integer n,iq,ilx,jlx,iad 
integer ncidx,iernc,varid,ndims
real, dimension(ifull,5), intent(out) :: svs,vlinprev,vlin,vlinnext
real, dimension(ifull_g,5) :: svsg,vling
real rlong0x,rlat0x,schmidtx,dsx,ra,rb,cablever
character(len=47) header  
character(len=6) vname
real, parameter :: cableversion = 223. ! version id for input data

write(6,*) "Reading land-use parameters for CABLE"
if (lncveg == 1) then
  ! assume this file grid has been tested when opened
  spos(1:3)=1
  npos(1)=il_g
  npos(2)=6*il_g
  npos(3)=1
  call ccnf_inq_dimlen(ncidveg,'longitude',ilx)
  call ccnf_inq_dimlen(ncidveg,'latitude',jlx)
  call ccnf_get_attg(ncidveg,'lon0',rlong0x)
  call ccnf_get_attg(ncidveg,'lat0',rlat0x)
  call ccnf_get_attg(ncidveg,'schmidt',schmidtx)
  if(ilx/=il_g.or.jlx/=jl_g.or.rlong0x/=rlong0.or.rlat0x/=rlat0.or.schmidtx/=schmidt) then
    write(6,*) 'wrong data file supplied ',trim(fveg)
    call ccmpi_abort(-1)
  end if
  call ccnf_get_attg(ncidveg,'cableversion',cablever,ierr=iernc)
  if (iernc /= 0) then
    write(6,*) "Missing version of CABLE data"
    write(6,*) "Regenerate land-use data with up-to-date version of igbpveg"
    call ccmpi_abort(-1)
  end if
  if (cablever /= cableversion) then
    write(6,*) "Wrong version of CABLE data"
    write(6,*) "Expecting ",cableversion
    write(6,*) "Found     ",cablever
    call ccmpi_abort(-1)
  end if
  do n=1,5
    write(vname,"(A,I1.1)") "lai",n
    call ccnf_inq_varid(ncidveg,vname,varid)
    call ccnf_inq_varndims(ncidveg,varid,ndims)
    call ccnf_get_vara(ncidveg,varid,spos(1:ndims),npos(1:ndims),vling(:,n)) 
    write(vname,"(A,I1.1)") "vegt",n
    call ccnf_inq_varid(ncidveg,vname,varid)
    call ccnf_inq_varndims(ncidveg,varid,ndims)
    call ccnf_get_vara(ncidveg,varid,spos(1:ndims),npos(1:ndims),svsg(:,n)) 
    ivsg(:,n)=nint(svsg(:,n))
    write(vname,"(A,I1.1)") "vfrac",n
    call ccnf_inq_varid(ncidveg,vname,varid)
    call ccnf_inq_varndims(ncidveg,varid,ndims)
    call ccnf_get_vara(ncidveg,varid,spos(1:ndims),npos(1:ndims),svsg(:,n))
  end do
  call ccmpi_distribute(ivs,ivsg)
  call ccmpi_distribute(svs,svsg)
  call ccmpi_distribute(vlin,vling)
  if (fvegprev/=' '.and.fvegnext/=' ') then
    call ccnf_open(fvegprev,ncidx,iernc)
    if (iernc/=0) then
      write(6,*) 'Cannot read netcdf file ',trim(fvegprev)
      call ccmpi_abort(-1)
    end if
    call ccnf_inq_dimlen(ncidx,'longitude',ilx)
    call ccnf_inq_dimlen(ncidx,'latitude',jlx)
    call ccnf_get_attg(ncidx,'lon0',rlong0x)
    call ccnf_get_attg(ncidx,'lat0',rlat0x)
    call ccnf_get_attg(ncidx,'schmidt',schmidtx)
    if(ilx/=il_g.or.jlx/=jl_g.or.rlong0x/=rlong0.or.rlat0x/=rlat0.or.schmidtx/=schmidt) then
      write(6,*) 'wrong data file supplied ',trim(fvegprev)
      call ccmpi_abort(-1)
    end if
    do n=1,5
      write(vname,"(A,I1.1)") "lai",n
      call ccnf_inq_varid(ncidveg,vname,varid)
      call ccnf_inq_varndims(ncidveg,varid,ndims)
      call ccnf_get_vara(ncidveg,varid,spos(1:ndims),npos(1:ndims),vling(:,n)) 
    end do
    call ccnf_close(ncidx)
    call ccmpi_distribute(vlinprev,vling)
    call ccnf_open(fvegnext,ncidx,iernc)
    if (iernc/=0) then
      write(6,*) 'Cannot read netcdf file ',trim(fvegnext)
      call ccmpi_abort(-1)
    end if
    call ccnf_inq_dimlen(ncidx,'longitude',ilx)
    call ccnf_inq_dimlen(ncidx,'latitude',jlx)
    call ccnf_get_attg(ncidx,'lon0',rlong0x)
    call ccnf_get_attg(ncidx,'lat0',rlat0x)
    call ccnf_get_attg(ncidx,'schmidt',schmidtx)
    if(ilx/=il_g.or.jlx/=jl_g.or.rlong0x/=rlong0.or.rlat0x/=rlat0.or.schmidtx/=schmidt) then
      write(6,*) 'wrong data file supplied ',trim(fvegnext)
      call ccmpi_abort(-1)
    end if
    do n=1,5
      write(vname,"(A,I1.1)") "lai",n
      call ccnf_inq_varid(ncidveg,vname,varid)
      call ccnf_inq_varndims(ncidveg,varid,ndims)
      call ccnf_get_vara(ncidveg,varid,spos(1:ndims),npos(1:ndims),vling(:,n)) 
    end do
    call ccnf_close(ncidx)
    call ccmpi_distribute(vlinnext,vling)
  else
    vlinprev=-1.
    vlinnext=-1.    
  end if

else
  open(87,file=fveg,status='old')
  read(87,*) ilx,jlx,rlong0x,rlat0x,schmidtx,dsx,header
  if(ilx/=il_g.or.jlx/=jl_g.or.rlong0x/=rlong0.or.rlat0x/=rlat0.or.schmidtx/=schmidt) then
    write(6,*) 'wrong data file supplied ',trim(fveg)
    call ccmpi_abort(-1)
  end if
  do iq=1,ifull_g
    read(87,*) iad,ra,rb,ivsg(iq,1),svsg(iq,1),vling(iq,1),ivsg(iq,2),svsg(iq,2),vling(iq,2),ivsg(iq,3),svsg(iq,3),vling(iq,3), &
               ivsg(iq,4),svsg(iq,4),vling(iq,4),ivsg(iq,5),svsg(iq,5),vling(iq,5)
  end do
  close(87)
  call ccmpi_distribute(ivs,ivsg)
  call ccmpi_distribute(svs,svsg)
  call ccmpi_distribute(vlin,vling)
  if (fvegprev/=' '.and.fvegnext/=' ') then
    open(87,file=fvegprev,status='old')
    read(87,*) ilx,jlx,rlong0x,rlat0x,schmidtx,dsx,header
    if(ilx/=il_g.or.jlx/=jl_g.or.rlong0x/=rlong0.or.rlat0x/=rlat0.or.schmidtx/=schmidt) then
      write(6,*) 'wrong data file supplied ',trim(fvegprev)
      call ccmpi_abort(-1)
    end if
    do iq=1,ifull_g
      read(87,*) iad,ra,rb,ivsg(iq,1),svsg(iq,1),vling(iq,1),ivsg(iq,2),svsg(iq,2),vling(iq,2),ivsg(iq,3),svsg(iq,3),vling(iq,3), &
                 ivsg(iq,4),svsg(iq,4),vling(iq,4),ivsg(iq,5),svsg(iq,5),vling(iq,5)
    end do
    close(87)
    call ccmpi_distribute(vlinprev,vling)
    open(87,file=fvegnext,status='old')
    read(87,*) ilx,jlx,rlong0x,rlat0x,schmidtx,dsx,header
    if(ilx/=il_g.or.jlx/=jl_g.or.rlong0x/=rlong0.or.rlat0x/=rlat0.or.schmidtx/=schmidt) then
      write(6,*) 'wrong data file supplied ',trim(fvegnext)
      call ccmpi_abort(-1)
    end if
    do iq=1,ifull_g
      read(87,*) iad,ra,rb,ivsg(iq,1),svsg(iq,1),vling(iq,1),ivsg(iq,2),svsg(iq,2),vling(iq,2),ivsg(iq,3),svsg(iq,3),vling(iq,3), &
                 ivsg(iq,4),svsg(iq,4),vling(iq,4),ivsg(iq,5),svsg(iq,5),vling(iq,5)
    end do
    close(87)
    call ccmpi_distribute(vlinnext,vling)
  else
    vlinprev=-1.
    vlinnext=-1.    
  end if
end if
return
end subroutine vegta
  
! vegtb is for myid != 0
subroutine vegtb(ivs,svs,vlinprev,vlin,vlinnext,fvegprev,fveg,fvegnext)
  
use cc_mpi
  
implicit none

include 'newmpar.h'

character(len=*), intent(in) :: fveg,fvegprev,fvegnext
integer, dimension(ifull,5), intent(out) :: ivs
real, dimension(ifull,5), intent(out) :: svs, vlinprev, vlin, vlinnext

call ccmpi_distribute(ivs)
call ccmpi_distribute(svs)
call ccmpi_distribute(vlin)
if (fvegprev/=' '.and.fvegnext/=' ') then
  call ccmpi_distribute(vlinprev)
  call ccmpi_distribute(vlinnext)
else
  vlinprev=-1.
  vlinnext=-1.
end if    
  
return
end subroutine vegtb

! *************************************************************************************  
! This subroutine loads CABLE tile data
subroutine loadtile

!~ use carbpools_m
use cc_mpi
use infile
use soil_m
use soilsnow_m
!~ use vegpar_m
  
implicit none

include 'newmpar.h'
include 'darcdf.h'
include 'parm.h'  
  
integer k, n, ierr, idv
integer, dimension(1) :: dum
real, dimension(ifull) :: dat
real totdepth
logical tst
character(len=11) vname

! check that CABLE data exists in restart file
! and communicate the result to all processors
! as not all processors are assigned an input file
ierr = 1
if ( io_in == 1 ) then
  if ( myid==0 .or. pfall ) then
    call ccnf_inq_varid(ncid,"tgg1_5",idv,tst)
    if ( tst ) then
      ierr = 1
    else
      ierr = 0
    end if
  end if
  if ( .not.pfall ) then
    dum(1) = ierr
    call ccmpi_bcast(dum(1:1),0,comm_world)
    ierr = dum(1)
  end if
end if
  
! Cannot locate tile data, use diagnostic data instead
if ( ierr /= 0 ) then
  if ( myid == 0 ) write(6,*) "Use gridbox averaged data to initialise CABLE"
  if ( mp > 0 ) then
    do n = 1,maxnb
      do k = 1,ms
        ssnow%tgg(pind(n,1):pind(n,2),k)   = pack(tgg(:,k),  tmap(:,n))
        ssnow%wb(pind(n,1):pind(n,2),k)    = pack(wb(:,k),   tmap(:,n))
        ssnow%wbice(pind(n,1):pind(n,2),k) = pack(wbice(:,k),tmap(:,n))
      end do
      do k = 1,3
        ssnow%tggsn(pind(n,1):pind(n,2),k)  = pack(tggsn(:,k),tmap(:,n))
        ssnow%smass(pind(n,1):pind(n,2),k)  = pack(smass(:,k),tmap(:,n))
        ssnow%ssdn(pind(n,1):pind(n,2),k)   = pack(ssdn(:,k), tmap(:,n))
        ssnow%sdepth(pind(n,1):pind(n,2),k) = pack(snowd/3.,  tmap(:,n))
        ssnow%sconds(pind(n,1):pind(n,2),k) = 0.2
      end do      
      ssnow%ssdnn(pind(n,1):pind(n,2))  = pack(ssdnn, tmap(:,n))
      ssnow%isflag(pind(n,1):pind(n,2)) = pack(isflag,tmap(:,n))
      ssnow%snowd(pind(n,1):pind(n,2))  = pack(snowd, tmap(:,n))
      ssnow%snage(pind(n,1):pind(n,2))  = pack(snage, tmap(:,n))
    end do
    ssnow%rtsoil=50.
    canopy%cansto=0.
    canopy%us=0.01
    ssnow%pudsto=0.
    ssnow%wetfac=0.
    ssnow%osnowd=ssnow%snowd
  end if
end if
  
! Some fixes for rounding errors
if ( mp > 0 ) then

  totdepth = 0.
  do k = 1,ms
    totdepth = totdepth + soil%zse(k)*100.
  enddo

  ssnow%wb = max(ssnow%wb,0._r_2)
  ssnow%wbice = max(ssnow%wbice,0._r_2)
  ssnow%smass = max(ssnow%smass,0.)
  ssnow%rtsoil = max(ssnow%rtsoil,0.)
  ssnow%snowd = max(ssnow%snowd,0.)
  ssnow%osnowd = max(ssnow%osnowd,0.)
  ssnow%wetfac = min(max(ssnow%wetfac,0.),1.)
  canopy%cansto = max(canopy%cansto,0.)

  ssnow%wbtot = 0.
  ssnow%wbtot1 = 0.
  ssnow%wbtot2 = 0.
  ssnow%tggav = 0.
  do k = 1,ms
    ssnow%wbtot = ssnow%wbtot+ssnow%wb(:,k)*1000.0*soil%zse(k)
    ssnow%tggav = ssnow%tggav+soil%zse(k)*ssnow%tgg(:,k)/(totdepth/100.)
    ssnow%gammzz(:,k) = max((1.-soil%ssat)*soil%css* soil%rhosoil                     &
        + real(ssnow%wb(:,k)-ssnow%wbice(:,k))*4.218e3* 1000.                         &
        + real(ssnow%wbice(:,k))*2.100e3*1000.*0.9,soil%css*soil%rhosoil)*soil%zse(k) &
        + (1.-ssnow%isflag)*2090.0*ssnow%snowd
  end do

  if ( icycle == 0 ) then
    bgc%cplant = max(bgc%cplant,0.)
    bgc%csoil = max(bgc%csoil,0.)
  else
    casapool%cplant     = max(0._r_2,casapool%cplant)
    casapool%clitter    = max(0._r_2,casapool%clitter)
    casapool%csoil      = max(0._r_2,casapool%csoil)
    casabal%cplantlast(1:mp,1:mplant)   = casapool%cplant(1:mp,1:mplant)
    casabal%clitterlast(1:mp,1:mlitter) = casapool%clitter(1:mp,1:mlitter)
    casabal%csoillast(1:mp,1:msoil)     = casapool%csoil(1:mp,1:msoil)
    casabal%clabilelast = casapool%clabile
    casabal%sumcbal     = 0.
    casabal%FCgppyear   = 0.
    casabal%FCrpyear    = 0.
    casabal%FCnppyear   = 0.
    casabal%FCrsyear    = 0.
    casabal%FCneeyear   = 0.
    casapool%nplant     = max(1.e-6_r_2,casapool%nplant)
    casapool%nlitter    = max(1.e-6_r_2,casapool%nlitter)
    casapool%nsoil      = max(1.e-6_r_2,casapool%nsoil)
    casapool%nsoilmin   = max(1.e-6_r_2,casapool%nsoilmin)
    casabal%nplantlast(1:mp,1:mplant)   = casapool%nplant(1:mp,1:mplant)
    casabal%nlitterlast(1:mp,1:mlitter) = casapool%nlitter(1:mp,1:mlitter)
    casabal%nsoillast(1:mp,1:msoil)     = casapool%nsoil(1:mp,1:msoil)      
    casabal%nsoilminlast= casapool%nsoilmin
    casabal%sumnbal     = 0.
    casabal%FNdepyear   = 0.
    casabal%FNfixyear   = 0.
    casabal%FNsnetyear  = 0.
    casabal%FNupyear    = 0.
    casabal%FNleachyear = 0.
    casabal%FNlossyear  = 0.
    casapool%pplant     = max(1.0e-7_r_2,casapool%pplant)
    casapool%plitter    = max(1.0e-7_r_2,casapool%plitter)
    casapool%psoil      = max(1.0e-7_r_2,casapool%psoil)
    casapool%Psoillab   = max(1.0e-7_r_2,casapool%psoillab)
    casapool%psoilsorb  = max(1.0e-7_r_2,casapool%psoilsorb)
    casapool%psoilocc   = max(1.0e-7_r_2,casapool%psoilocc)
    casabal%pplantlast(1:mp,1:mplant)   = casapool%pplant(1:mp,1:mplant)
    casabal%plitterlast(1:mp,1:mlitter) = casapool%plitter(1:mp,1:mlitter)
    casabal%psoillast(1:mp,1:msoil)     = casapool%psoil(1:mp,1:msoil)       
    casabal%psoillablast= casapool%psoillab
    casabal%psoilsorblast=casapool%psoilsorb
    casabal%psoilocclast= casapool%psoilocc
    casabal%sumpbal     = 0.
    casabal%FPweayear   = 0.
    casabal%FPdustyear  = 0.
    casabal%FPsnetyear  = 0.
    casabal%FPupyear    = 0.
    casabal%FPleachyear = 0.
    casabal%FPlossyear  = 0.
  end if
end if
  
return
end subroutine loadtile

! *************************************************************************************
! Water inflow from river routing
subroutine cableinflow(inflow,rate)

use soil_m

implicit none

include 'newmpar.h'

integer nb, k
real, dimension(ifull), intent(in) :: rate
real, dimension(ifull), intent(inout) :: inflow
real, dimension(ifull) :: delflow
real, dimension(mp) :: xx, ll, delxx, ratepack

if ( mp <= 0 ) return

do nb = 1,maxnb
  xx(pind(nb,1):pind(nb,2)) = pack( inflow(1:ifull), tmap(:,nb) )
  ratepack(pind(nb,1):pind(nb,2)) = pack( rate(1:ifull), tmap(:,nb) )
end do
delxx(1:mp) = 0.
do k = 1,cbm_ms
  ll(1:mp) = max( soil%sfc(1:mp)-real(ssnow%wb(1:mp,k)), 0. )*1000.*soil%zse(k)
  ll(1:mp) = ll(1:mp)*ratepack(1:mp)
  ll(1:mp) = min( xx(1:mp), ll(1:mp) )
  ssnow%wb(1:mp,k) = ssnow%wb(1:mp,k) + ll(1:mp)/(1000.*soil%zse(k))
  delxx(1:mp) = delxx(1:mp) - ll(1:mp)
end do
delflow(1:ifull) = 0.
do nb = 1,maxnb
  delflow(1:ifull) = delflow(1:ifull) + unpack(sv(pind(nb,1):pind(nb,2))*delxx(pind(nb,1):pind(nb,2)),tmap(:,nb),0.)
end do
inflow(1:ifull) = inflow(1:ifull) + delflow(1:ifull)

return
end subroutine cableinflow

end module cable_ccam

