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

! CCAM netCDF output routines

! itype=1     write outfile history file (compressed)
! itype=-1    write restart file (uncompressed)
! localhist=f single processor output 
! localhist=t parallel output for each processor

! Thanks to Paul Ryan for advice on output netcdf routines.
    
module outcdf
    
private
public outfile, mslp

character(len=3), dimension(12), parameter :: month = (/'jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'/)

contains

subroutine outfile(iout,rundate,nwrite,nstagin,jalbfix,nalpha,mins_rad)
      
use arrays_m
use cc_mpi
use pbl_m
use soilsnow_m ! tgg,wb,snowd
!~ use tracers_m
      
implicit none
      
include 'newmpar.h'
include 'dates.h'    ! mtimer
include 'filnames.h' ! list of files, read in once only
include 'parm.h'

integer iout,nwrite,nstagin
integer, intent(in) :: jalbfix,nalpha,mins_rad
character(len=160) :: co2out,radonout,surfout
character(len=20) :: qgout
character(len=8) :: rundate

call START_LOG(outfile_begin)
      
if ( myid==0 ) then
  write(6,*) "ofile written for iout: ",iout
  write(6,*) "kdate,ktime,mtimer:     ",kdate,ktime,mtimer
end if

if ( nrungcm==-2 .or. nrungcm==-3 .or. nrungcm==-5 ) then
  if ( ktau==nwrite/2 .or. ktau==nwrite ) then
!        usually after first 24 hours, save soil variables for next run
    if ( ktau==nwrite ) then  ! 24 hour write
      if ( ktime==1200 ) then
        co2out=co2_12     ! 'co2.1200'
        radonout=radon_12 ! 'radon.1200'
        surfout=surf_12   ! 'current.1200'
        qgout='qg_12'
      else
        co2out=co2_00     !  'co2.0000'
        radonout=radon_00 ! 'radon.0000'
        surfout=surf_00   ! 'current.0000'
        qgout='qg_00'
      endif
    else                    ! 12 hour write
      if(ktime==1200)then
        co2out=co2_00     !  'co2.0000'
        radonout=radon_00 ! 'radon.0000'
        surfout=surf_00   ! 'current.0000'
        qgout='qg_00'
      else
        co2out=co2_12     ! 'co2.1200'
        radonout=radon_12 ! 'radon.1200'
        surfout=surf_12   ! 'current.1200'
        qgout='qg_12'
      endif
    endif               ! (ktau.eq.nwrite)
    if ( myid == 0 ) then
      write(6,*) "writing current soil & snow variables to ",surfout
      open(unit=77,file=surfout,form='formatted',status='unknown')
      write (77,*) kdate,ktime,' ktau = ',ktau
    end if
    call writeglobvar(77, wb, fmt='(14f6.3)')
    call writeglobvar(77, tss, fmt='(12f7.2)')
    call writeglobvar(77, snowd, fmt='(12f7.1)')
    if ( myid == 0 ) close (77)
    if ( nrungcm==-2 .or. nrungcm==-5 ) then
      if ( myid == 0 ) then
        write(6,*) "writing special qgout file: ",qgout
        open(unit=77,file=qgout,form='unformatted',status='unknown')
      end if
      call writeglobvar(77, qg)
      if ( myid == 0 ) close (77)
    endif  ! (nrungcm.eq.-2.or.nrungcm.eq.-5)
  endif    ! (ktau.eq.nwrite/2.or.ktau.eq.nwrite)
endif      ! (nrungcm.eq.-2.or.nrungcm.eq.-3.or.nrungcm.eq.-5)

!---------------------------------------------------------------------------
if ( iout==19 ) then
  select case(io_rest)  
    case(1)  ! for netCDF 
      if ( myid==0 ) write(6,*) "restart write of data to netCDF"
      call cdfout(rundate,-1,nstagin,jalbfix,nalpha,mins_rad)
    case(3)
      write(6,*) "Error, restart binary output not supported"
      call ccmpi_abort(-1)
  end select
else
  select case(io_out)
    case(1)
      call cdfout(rundate,1,nstagin,jalbfix,nalpha,mins_rad)
    case(3)
      write(6,*) "Error, history binary output not supported"
      call ccmpi_abort(-1)
  end select
end if

call END_LOG(outfile_end)
      
return
end subroutine outfile

    
!--------------------------------------------------------------
! CONFIGURE DIMENSIONS FOR OUTPUT NETCDF FILES
subroutine cdfout(rundate,itype,nstagin,jalbfix,nalpha,mins_rad)


use cable_ccam, only : proglai        ! CABLE
use cc_mpi                            ! CC MPI routines
use infile                            ! Input file routines
use liqwpar_m                         ! Cloud water mixing ratios
!~ use mlo, only : mindep              & ! Ocean physics and prognostic arrays
    !~ ,minwater,mxd,zomode,zoseaice   &
    !~ ,factchseaice
use parmhdff_m                        ! Horizontal diffusion parameters
!~ use tkeeps                            ! TKE-EPS boundary layer

implicit none

include 'newmpar.h'                   ! Grid parameters
include 'dates.h'                     ! Date data
include 'filnames.h'                  ! Filenames
include 'kuocom.h'                    ! Convection parameters
include 'parm.h'                      ! Model configuration
include 'parmdyn.h'                   ! Dynamics parameters
include 'parmgeom.h'                  ! Coordinate data
include 'parmhor.h'                   ! Horizontal advection parameters
include 'parmsurf.h'                  ! Surface parameters

integer ixp,iyp,idlev,idnt,idms,idoc
integer leap
common/leap_yr/leap                   ! Leap year (1 to allow leap years)
integer nbarewet,nsigmf
common/nsib/nbarewet,nsigmf

integer, parameter :: nihead=54
integer, parameter :: nrhead=14
integer, dimension(nihead) :: nahead
integer, dimension(4), save :: dima,dims,dimo
integer, intent(in) :: jalbfix,nalpha,mins_rad
integer itype, nstagin
integer xdim,ydim,zdim,tdim,msdim,ocdim
integer icy, icm, icd, ich, icmi, ics, idv
integer namipo3
integer, save :: idnc=0, iarch=0
real, dimension(nrhead) :: ahead
character(len=180) cdffile
character(len=33) grdtim
character(len=20) timorg
character(len=8) rundate

! Determine file names depending on output
if ( myid==0 .or. localhist ) then
  ! File setup follows
  if ( itype==1 ) then
    ! itype=1 outfile
    iarch=iarch+1
    if ( localhist ) then
      write(cdffile,"(a,'.',i6.6)") trim(ofile), myid
    else
      cdffile=ofile
    endif
  else
    ! itype=-1 restfile
    iarch=1
    if ( localhist ) then
      write(cdffile,"(a,'.',i6.6)") trim(restfile), myid
    else
      cdffile=restfile
    endif
    idnc=0
  endif ! ( itype==1)then

  ! Open new file
  if( iarch==1 )then
    if ( myid==0 ) write(6,'(" nccre of itype,cdffile=",i5," ",a80)') itype,cdffile
    call ccnf_create(cdffile,idnc)
    ! Turn off the data filling
    call ccnf_nofill(idnc)
    ! Create dimensions, lon, runtopo.shlat
    if( localhist ) then
      call ccnf_def_dim(idnc,'longitude',il,xdim)
      call ccnf_def_dim(idnc,'latitude',jl,ydim)
    else
      call ccnf_def_dim(idnc,'longitude',il_g,xdim)
      call ccnf_def_dim(idnc,'latitude',jl_g,ydim)
    endif
    call ccnf_def_dim(idnc,'lev',kl,zdim)
    call ccnf_def_dim(idnc,'zsoil',ms,msdim)
    if ( abs(nmlo)>0. .and. abs(nmlo)<=9 ) then
      call ccnf_def_dim(idnc,'olev',ol,ocdim)
    else
      ocdim=0
    end if
    call ccnf_def_dimu(idnc,'time',tdim)
    if ( myid==0 ) then
      write(6,*) "xdim,ydim,zdim,tdim"
      write(6,*)  xdim,ydim,zdim,tdim
    end if

    ! atmosphere dimensions
    dima = (/ xdim, ydim, zdim, tdim /)

    ! soil dimensions
    dims = (/ xdim, ydim, msdim, tdim /)

    ! ocean dimensions
    dimo = (/ xdim, ydim, ocdim, tdim /)

    ! Define coords.
    call ccnf_def_var(idnc,'longitude','float',1,dima(1:1),ixp)
    call ccnf_put_att(idnc,ixp,'point_spacing','even')
    call ccnf_put_att(idnc,ixp,'units','degrees_east')
    call ccnf_def_var(idnc,'latitude','float',1,dima(2:2),iyp)
    call ccnf_put_att(idnc,iyp,'point_spacing','even')
    call ccnf_put_att(idnc,iyp,'units','degrees_north')
    if ( myid==0 ) write(6,*) 'ixp,iyp=',ixp,iyp

    call ccnf_def_var(idnc,'lev','float',1,dima(3:3),idlev)
    call ccnf_put_att(idnc,idlev,'positive','down')
    call ccnf_put_att(idnc,idlev,'point_spacing','uneven')
    call ccnf_put_att(idnc,idlev,'units','sigma_level')
    call ccnf_put_att(idnc,idlev,'long_name','sigma_level')
    if (myid==0) write(6,*) 'idlev=',idlev

    call ccnf_def_var(idnc,'zsoil','float',1,dims(3:3),idms)
    call ccnf_put_att(idnc,idms,'point_spacing','uneven')
    call ccnf_put_att(idnc,idms,'units','m')
    if (myid==0) write(6,*) 'idms=',idms
        
    if (abs(nmlo)>0.and.abs(nmlo)<=9) then
      call ccnf_def_var(idnc,'olev','float',1,dimo(3:3),idoc)
      call ccnf_put_att(idnc,idoc,'point_spacing','uneven')
      call ccnf_put_att(idnc,idoc,'units','sigma_level')
      if (myid==0) write(6,*) 'idoc=',idoc
    end if

    call ccnf_def_var(idnc,'time','float',1,dima(4:4),idnt)
    call ccnf_put_att(idnc,idnt,'point_spacing','even')
    if ( myid==0 ) then
      write(6,*) 'tdim,idnc=',tdim,idnc
      write(6,*) 'idnt=',idnt
      write(6,*) 'kdate,ktime,ktau=',kdate,ktime,ktau
    end if

    icy = kdate/10000
    icm = max(1, min(12, (kdate-icy*10000)/100))
    icd = max(1, min(31, (kdate-icy*10000-icm*100)))
    if ( icy<100 ) icy = icy + 1900
    ich = ktime/100
    icmi = (ktime-ich*100)
    ics = 0
    write(timorg,'(i2.2,"-",a3,"-",i4.4,3(":",i2.2))') icd,month(icm),icy,ich,icmi,ics
    call ccnf_put_att(idnc,idnt,'time_origin',timorg)
    write(grdtim,'("minutes since ",i4.4,"-",i2.2,"-",i2.2," ",2(i2.2,":"),i2.2)') icy,icm,icd,ich,icmi,ics
    call ccnf_put_att(idnc,idnt,'units',grdtim)
    if ( leap==0 ) then
      call ccnf_put_att(idnc,idnt,'calendar','noleap')
    end if
    if ( myid==0 ) then
      write(6,*) 'timorg=',timorg
      write(6,*) 'grdtim=',grdtim
    end if

!   create the attributes of the header record of the file
    nahead(1) = il_g       ! needed by cc2hist
    nahead(2) = jl_g       ! needed by cc2hist
    nahead(3) = kl         ! needed by cc2hist
    nahead(4) = 5
    nahead(5) = 0          ! nsd not used now
    nahead(6) = io_in
    nahead(7) = nbd
    nahead(8) = 0          ! not needed now  
    nahead(9) = mex
    nahead(10) = mup
    nahead(11) = 2 ! nem
    nahead(12) = mtimer
    nahead(13) = 0         ! nmi
    nahead(14) = nint(dt)  ! needed by cc2hist
    nahead(15) = 0         ! not needed now 
    nahead(16) = nhor
    nahead(17) = nkuo
    nahead(18) = khdif
    nahead(19) = kl        ! needed by cc2hist (was kwt)
    nahead(20) = 0  !iaa
    nahead(21) = 0  !jaa
    nahead(22) = -4
    nahead(23) = 0       ! not needed now      
    nahead(24) = 0  !lbd
    nahead(25) = nrun
    nahead(26) = 0
    nahead(27) = khor
    nahead(28) = ksc
    nahead(29) = kountr
    nahead(30) = 1 ! ndiur
    nahead(31) = 0  ! spare
    nahead(32) = nhorps
    nahead(33) = nsoil
    nahead(34) = ms        ! needed by cc2hist
    nahead(35) = ntsur
    nahead(36) = nrad
    nahead(37) = kuocb
    nahead(38) = nvmix
    nahead(40) = 0    
    nahead(41) = nextout

    nahead(44) = nsib
    nahead(45) = nrungcm
    nahead(46) = ncvmix
    nahead(47) = ngwd
    nahead(48) = lgwd
    nahead(49) = mup
    nahead(50) = nritch_t
    nahead(51) = ldr
    nahead(52) = nevapls
    nahead(53) = nevapcc
    nahead(54) = nt_adv
    ahead(1) = ds
    ahead(2) = 0.  !difknbd
    ahead(3) = 0.  ! was rhkuo for kuo scheme
    ahead(4) = 0.  !du
    ahead(5) = rlong0     ! needed by cc2hist
    ahead(6) = rlat0      ! needed by cc2hist
    ahead(7) = schmidt    ! needed by cc2hist
    ahead(8) = 0.  !stl2
    ahead(9) = 0.  !relaxt
    ahead(10) = 0.  !hourbd
    ahead(11) = tss_sh
    ahead(12) = vmodmin
    ahead(13) = av_vmod
    ahead(14) = epsp
    if ( myid==0 ) then
      write(6,'(" nahead=",(20i4))') nahead
      write(6,*) "ahead=",ahead
    end if
    call ccnf_put_attg(idnc,'int_header',nahead)
    call ccnf_put_attg(idnc,'real_header',ahead)
    call ccnf_put_attg(idnc,'date_header',rundate)
    call ccnf_def_var(idnc,'ds','float',idv)
    call ccnf_def_var(idnc,'dt','float',idv)


  else
    if ( myid==0 ) write(6,'(" outcdf itype,idnc,iarch,cdffile=",i5,i8,i5," ",a80)') itype,idnc,iarch,cdffile
  endif ! ( iarch=1 ) ..else..
endif ! (myid==0.or.localhist)
      
! openhist writes some fields so needs to be called by all processes
call openhist(iarch,itype,dima,localhist,idnc,nstagin,ixp,iyp,idlev,idms,idoc)

if ( myid==0 .or. localhist ) then
  if ( ktau==ntau ) then
    if ( myid==0 ) write(6,*) "closing netCDF file idnc=",idnc      
    call ccnf_close(idnc)
  endif
endif    ! (myid==0.or.local)

return
end subroutine cdfout
      
!--------------------------------------------------------------
! CREATE ATTRIBUTES AND WRITE OUTPUT
subroutine openhist(iarch,itype,idim,local,idnc,nstagin,ixp,iyp,idlev,idms,idoc)

use arrays_m                                     ! Atmosphere dyamics prognostic arrays
use ateb, only : atebsave                        ! Urban
use cable_ccam, only : savetile, savetiledef     ! CABLE interface
use cable_def_types_mod, only : ncs, ncp         ! CABLE dimensions
use casadimension, only : mplant, mlitter, msoil ! CASA dimensions
use carbpools_m                                  ! Carbon pools
use cc_mpi                                       ! CC MPI routines
use cfrac_m                                      ! Cloud fraction
use dpsdt_m                                      ! Vertical velocity
use extraout_m                                   ! Additional diagnostics
!~ use gdrag_m                                      ! Gravity wave drag
use histave_m                                    ! Time average arrays
use infile                                       ! Input file routines
use latlong_m                                    ! Lat/lon coordinates
use liqwpar_m                                    ! Cloud water mixing ratios
use map_m                                        ! Grid map arrays
!~ use mlo, only : wlev,mlosave,mlodiag, &          ! Ocean physics and prognostic arrays
                !~ mloexpdep,wrtemp
use morepbl_m                                    ! Additional boundary layer diagnostics
use nharrs_m                                     ! Non-hydrostatic atmosphere arrays
use nsibd_m                                      ! Land-surface arrays
use pbl_m                                        ! Boundary layer arrays
use prec_m                                       ! Precipitation
use raddiag_m                                    ! Radiation diagnostic
use river                                        ! River routing
use savuvt_m                                     ! Saved dynamic arrays
use savuv1_m                                     ! Saved dynamic arrays
!~ use screen_m                                     ! Screen level diagnostics
use sigs_m                                       ! Atmosphere sigma levels
use soil_m                                       ! Soil and surface data
use soilsnow_m                                   ! Soil, snow and surface data
!~ use tkeeps, only : tke,eps,zidry                 ! TKE-EPS boundary layer
use vegpar_m                                     ! Vegetation arrays
use vvel_m                                       ! Additional vertical velocity
use work2_m                                      ! Diagnostic arrays
use xarrs_m, only : pslx                         ! Saved dynamic arrays

implicit none

include 'newmpar.h'                              ! Grid parameters
include 'const_phys.h'                           ! Physical constants
include 'dates.h'                                ! Date data
include 'filnames.h'                             ! Filenames
include 'kuocom.h'                               ! Convection parameters
include 'parm.h'                                 ! Model configuration
include 'parmdyn.h'                              ! Dynamics parameters
include 'soilv.h'                                ! Soil parameters
include 'version.h'                              ! Model version data

integer ixp,iyp,idlev,idms,idoc
integer i, idkdate, idktau, idktime, idmtimer, idnteg, idnter
integer idv, iq, j, k, n, igas, idnc
integer iarch, itype, nstagin, idum
integer, dimension(4), intent(in) :: idim
integer, dimension(3) :: jdim
real, dimension(ms) :: zsoil
real, dimension(il_g) :: xpnt
real, dimension(jl_g) :: ypnt
real, dimension(ifull) :: aa
real, dimension(ifull) :: ocndep,ocnheight
real, dimension(ifull) :: qtot, tv
real, dimension(ifull,kl) :: tmpry,rhoa
!~ real, dimension(ifull,wlev,4) :: mlodwn
real, dimension(ifull,11) :: micdwn
real, dimension(ifull,28) :: atebdwn
character(len=50) expdesc
character(len=50) lname
character(len=21) mnam,nnam
character(len=8) vname
character(len=3) trnum
logical, intent(in) :: local
logical lwrite,lave,lrad,lday
logical l3hr

lwrite=ktau>0
lave=mod(ktau,nperavg)==0.or.ktau==ntau
lave=lave.and.ktau>0
lrad=mod(ktau,kountr)==0.or.ktau==ntau
lrad=lrad.and.ktau>0
lday=mod(ktau,nperday)==0
lday=lday.and.ktau>0
l3hr=(real(nwt)*dt>10800.)

! idim is for 4-D (3 dimensions+time)
! jdim is for 3-D (2 dimensions+time)
jdim(1:2)=idim(1:2)
jdim(3)=idim(4)

if( myid==0 .or. local ) then

! if this is the first archive, set up some global attributes
  if ( iarch==1 ) then

!   Create global attributes
!   Model run number
    if ( myid==0 ) then
      write(6,*) 'idim=',idim
      write(6,*) 'nrun=',nrun
    end if
    call ccnf_put_attg(idnc,'nrun',nrun)

!   Experiment description
    expdesc = 'CCAM model run'
    call ccnf_put_attg(idnc,'expdesc',expdesc)

!   Model version
    call ccnf_put_attg(idnc,'version',version)

    if ( local ) then
      call ccnf_put_attg(idnc,'processor_num',myid)
      call ccnf_put_attg(idnc,'nproc',nproc)
#ifdef uniform_decomp
      call ccnf_put_attg(idnc,'decomp','uniform1')
#else
      call ccnf_put_attg(idnc,'decomp','face')
#endif
    endif           

!       Sigma levels
    if ( myid==0 ) write(6,*) 'sig=',sig
    call ccnf_put_attg(idnc,'sigma',sig)

    lname = 'year-month-day at start of run'
    call ccnf_def_var(idnc,'kdate','int',1,idim(4:4),idkdate)
    call ccnf_put_att(idnc,idkdate,'long_name',lname)

    lname = 'hour-minute at start of run'
    call ccnf_def_var(idnc,'ktime','int',1,idim(4:4),idktime)
    call ccnf_put_att(idnc,idktime,'long_name',lname)

    lname = 'timer (hrs)'
    call ccnf_def_var(idnc,'timer','float',1,idim(4:4),idnter)
    call ccnf_put_att(idnc,idnter,'long_name',lname)

    lname = 'mtimer (mins)'
    call ccnf_def_var(idnc,'mtimer','int',1,idim(4:4),idmtimer)
    call ccnf_put_att(idnc,idmtimer,'long_name',lname)

    lname = 'timeg (UTC)'
    call ccnf_def_var(idnc,'timeg','float',1,idim(4:4),idnteg)
    call ccnf_put_att(idnc,idnteg,'long_name',lname)

    lname = 'number of time steps from start'
    call ccnf_def_var(idnc,'ktau','int',1,idim(4:4),idktau)
    call ccnf_put_att(idnc,idktau,'long_name',lname)

    lname = 'down'
    call ccnf_def_var(idnc,'sigma','float',1,idim(3:3),idv)
    call ccnf_put_att(idnc,idv,'positive',lname)

    lname = 'atm stag direction'
    call ccnf_def_var(idnc,'nstag','int',1,idim(4:4),idv)
    call ccnf_put_att(idnc,idv,'long_name',lname)

    lname = 'atm unstag direction'
    call ccnf_def_var(idnc,'nstagu','int',1,idim(4:4),idv)
    call ccnf_put_att(idnc,idv,'long_name',lname)

    lname = 'atm stag offset'
    call ccnf_def_var(idnc,'nstagoff','int',1,idim(4:4),idv)
    call ccnf_put_att(idnc,idv,'long_name',lname)

    if ( (nmlo<0.and.nmlo>=-9) .or. (nmlo>0.and.nmlo<=9.and.itype==-1) ) then
      lname = 'ocn stag offset'
      call ccnf_def_var(idnc,'nstagoffmlo','int',1,idim(4:4),idv)
      call ccnf_put_att(idnc,idv,'long_name',lname)     
    end if

    if ( myid==0 ) write(6,*) 'define attributes of variables'

!   For time invariant surface fields
    lname = 'Surface geopotential'
    call attrib(idnc,idim(1:2),2,'zht',lname,'m2/s2',-1000.,90.e3,0,-1)
    !~ lname = 'Std Dev of surface height'
    !~ call attrib(idnc,idim(1:2),2,'he',lname,'m',0.,90.e3,0,-1)
    lname = 'Map factor'
    call attrib(idnc,idim(1:2),2,'map',lname,'none',.001,1500.,0,itype)

!   For time varying surface fields
    lname ='Scaled Log Surface pressure'
    call attrib(idnc,jdim(1:3),3,'psf',lname,'none',-1.3,0.2,0,itype)

    lname = 'Runoff'
    call attrib(idnc,jdim(1:3),3,'runoff',lname,'mm/day',0.,1300.,0,-1) ! -1=long
 
    if ( nmlo<=-2 .or. (nmlo>=2.and.itype==-1) .or. nriver==1 ) then
      lname = 'Surface water depth'
      call attrib(idnc,jdim(1:3),3,'swater',lname,'mm',0.,6.5E3,0,-1) ! -1 = long
    end if

    lname = 'Wetness fraction layer 1' ! 5. for frozen sand
    call attrib(idnc,jdim(1:3),3,'wetfrac1',lname,'none',-6.5,6.5,0,itype)
    lname = 'Wetness fraction layer 2'
    call attrib(idnc,jdim(1:3),3,'wetfrac2',lname,'none',-6.5,6.5,0,itype)
    lname = 'Wetness fraction layer 3'
    call attrib(idnc,jdim(1:3),3,'wetfrac3',lname,'none',-6.5,6.5,0,itype)
    lname = 'Wetness fraction layer 4'
    call attrib(idnc,jdim(1:3),3,'wetfrac4',lname,'none',-6.5,6.5,0,itype)
    lname = 'Wetness fraction layer 5'
    call attrib(idnc,jdim(1:3),3,'wetfrac5',lname,'none',-6.5,6.5,0,itype)
    lname = 'Wetness fraction layer 6'
    call attrib(idnc,jdim(1:3),3,'wetfrac6',lname,'none',-6.5,6.5,0,itype)
     
    ! PH - Add wetfac to output for mbase=-19 option
    lname = 'Surface wetness fraction'
    call attrib(idnc,jdim(1:3),3,'wetfac',lname,'none',-6.5,6.5,0,itype)
    
    ! STANDARD 3D VARIABLES -------------------------------------
    if ( myid==0 ) write(6,*) '3d variables'

    ! RESTART ---------------------------------------------------

    if ( myid==0 ) write(6,*) 'finished defining attributes'
!   Leave define mode
    call ccnf_enddef(idnc)
    if ( myid==0 ) write(6,*) 'leave define mode'

    if ( local ) then
      ! Set these to global indices (relative to panel 0 in uniform decomp)
      do i=1,ipan
        xpnt(i) = float(i) + ioff
      end do
      call ccnf_put_vara(idnc,ixp,1,il,xpnt(1:il))
      i=1
      do n=1,npan
        do j=1,jpan
          ypnt(i) = float(j) + joff + (n-noff)*il_g
          i=i+1
        end do
      end do
      call ccnf_put_vara(idnc,iyp,1,jl,ypnt(1:jl))
    else
      do i=1,il_g
        xpnt(i) = float(i)
      end do
      call ccnf_put_vara(idnc,ixp,1,il_g,xpnt(1:il_g))
      do j=1,jl_g
        ypnt(j) = float(j)
      end do
      call ccnf_put_vara(idnc,iyp,1,jl_g,ypnt(1:jl_g))
    endif

    call ccnf_put_vara(idnc,idlev,1,kl,sig)
    call ccnf_put_vara(idnc,'sigma',1,kl,sig)

    zsoil(1)=0.5*zse(1)
    zsoil(2)=zse(1)+zse(2)*0.5
    zsoil(3)=zse(1)+zse(2)+zse(3)*0.5
    zsoil(4)=zse(1)+zse(2)+zse(3)+zse(4)*0.5
    zsoil(5)=zse(1)+zse(2)+zse(3)+zse(4)+zse(5)*0.5
    zsoil(6)=zse(1)+zse(2)+zse(3)+zse(4)+zse(5)+zse(6)*0.5
    call ccnf_put_vara(idnc,idms,1,ms,zsoil)
    call ccnf_put_vara(idnc,'ds',1,ds)
    call ccnf_put_vara(idnc,'dt',1,dt)
  endif ! iarch==1
! -----------------------------------------------------------      

  ! set time to number of minutes since start 
  call ccnf_put_vara(idnc,'time',iarch,real(mtimer))
  call ccnf_put_vara(idnc,'timer',iarch,timer)
  call ccnf_put_vara(idnc,'mtimer',iarch,mtimer)
  call ccnf_put_vara(idnc,'timeg',iarch,timeg)
  call ccnf_put_vara(idnc,'ktau',iarch,ktau)
  call ccnf_put_vara(idnc,'kdate',iarch,kdate)
  call ccnf_put_vara(idnc,'ktime',iarch,ktime)
  call ccnf_put_vara(idnc,'nstag',iarch,nstag)
  call ccnf_put_vara(idnc,'nstagu',iarch,nstagu)
  idum=mod(ktau-nstagoff,max(abs(nstagin),1))
  idum=idum-max(abs(nstagin),1) ! should be -ve
  call ccnf_put_vara(idnc,'nstagoff',iarch,idum)
  if ( (nmlo<0.and.nmlo>=-9) .or. (nmlo>0.and.nmlo<=9.and.itype==-1) ) then
    call ccnf_put_vara(idnc,'nstagoffmlo',iarch,idum)
  end if
  if ( myid==0 ) then
    write(6,*) 'kdate,ktime,ktau=',kdate,ktime,ktau
    write(6,*) 'timer,timeg=',timer,timeg
  end if
       
endif ! myid == 0 .or. local     

!**************************************************************
! WRITE TIME-INVARIANT VARIABLES
!**************************************************************

if ( ktau==0 .or. itype==-1 ) then  ! also for restart file
  call histwrt3(zs,'zht',idnc,iarch,local,.true.)
  !~ call histwrt3(he,'he',idnc,iarch,local,.true.)
  call histwrt3(em,'map',idnc,iarch,local,.true.)
endif ! (ktau==0.or.itype==-1) 

!**************************************************************
! WRITE 3D VARIABLES (2D + Time)
!**************************************************************

! BASIC -------------------------------------------------------
aa(:) = runoff(1:ifull)*real(nperday)/real(min(nwt,max(ktau,1)))
call histwrt3(aa,'runoff',idnc,iarch,local,lwrite)


! MLO ---------------------------------------------------------      
if ( nmlo/=0 .and. abs(nmlo)<=9 ) then
  ocnheight = min(max(ocnheight,-130.),130.)
  where (.not.land(1:ifull))
    snowd   = micdwn(:,7)*1000.
  end where
end if

if ( nmlo<=-2 .or. (nmlo>=2.and.itype==-1) .or. nriver==1 ) then
  call histwrt3(watbdy(1:ifull),'swater',idnc,iarch,local,.true.)
end if

! SOIL --------------------------------------------------------
aa(:)=(wb(:,1)-swilt(isoilm))/(sfc(isoilm)-swilt(isoilm))
call histwrt3(aa,'wetfrac1',idnc,iarch,local,.true.)
aa(:)=(wb(:,2)-swilt(isoilm))/(sfc(isoilm)-swilt(isoilm))
call histwrt3(aa,'wetfrac2',idnc,iarch,local,.true.)
aa(:)=(wb(:,3)-swilt(isoilm))/(sfc(isoilm)-swilt(isoilm))
call histwrt3(aa,'wetfrac3',idnc,iarch,local,.true.)
aa(:)=(wb(:,4)-swilt(isoilm))/(sfc(isoilm)-swilt(isoilm))
call histwrt3(aa,'wetfrac4',idnc,iarch,local,.true.)
aa(:)=(wb(:,5)-swilt(isoilm))/(sfc(isoilm)-swilt(isoilm))
call histwrt3(aa,'wetfrac5',idnc,iarch,local,.true.)
aa(:)=(wb(:,6)-swilt(isoilm))/(sfc(isoilm)-swilt(isoilm))
call histwrt3(aa,'wetfrac6',idnc,iarch,local,.true.)
      
! PH - Add wetfac to output for mbase=-19 option
call histwrt3(wetfac,'wetfac',idnc,iarch,local,.true.)
      
! DIAGNOSTICS -------------------------------------------------
lwrite=(ktau>0)
      

! **************************************************************
! WRITE 4D VARIABLES (3D + Time)
! **************************************************************

if ( myid==0 .or. local ) then
  call ccnf_sync(idnc)
end if

if ( myid==0 ) then
  write(6,*) "finished writing to ofile"    
end if

return
end subroutine openhist

end module outcdf
