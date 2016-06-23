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

!      PE model on conformal-cubic grid
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      input files are :namelist (via file called "input")
!                       "nrun.dat"
!      data input and output file names are specified in namelist 'datafile'
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     sign convention:
!                      u+ve eastwards  (on the panel)
!                      v+ve northwards (on the panel)

program globpe

use arrays_m                               ! Atmosphere dyamics prognostic arrays
use bigxy4_m                               ! Grid interpolation
use carbpools_m, only : carbpools_init   & ! Carbon pools
    ,fpn,frs,frp
use cc_mpi                                 ! CC MPI routines
use diag_m                                 ! Diagnostic routines
use extraout_m                             ! Additional diagnostics
use gdrag_m, only : gdrag_init             ! Gravity wave drag
use histave_m                              ! Time average arrays
use indata                                 ! Data initialisation
use indices_m                              ! Grid index arrays
use infile                                 ! Input file routines
use latlong_m                              ! Lat/lon coordinates
use map_m                                  ! Grid map arrays
use morepbl_m                              ! Additional boundary layer diagnostics
use nharrs_m, only : nharrs_init         & ! Non-hydrostatic atmosphere arrays
   ,lrestart
use nlin_m                                 ! Atmosphere non-linear dynamics
use nsibd_m                                ! Land-surface arrays
use outcdf                                 ! Output file routines
!~ use parmhdff_m                             ! Horizontal diffusion parameters
use pbl_m                                  ! Boundary layer arrays
use permsurf_m, only : permsurf_init       ! Fixed surface arrays
use river                                  ! River routing
use savuvt_m                               ! Saved dynamic arrays
use savuv1_m                               ! Saved dynamic arrays
use sbar_m                                 ! Saved dynamic arrays
use sigs_m                                 ! Atmosphere sigma levels
use soil_m                                 ! Soil and surface data
use soilsnow_m                             ! Soil, snow and surface data
use vecs_m, only : vecs_init               ! Eigenvectors for atmosphere dynamics
use vecsuv_m                               ! Map to cartesian coordinates
use vegpar_m                               ! Vegetation arrays
use work2_m                                ! Diagnostic arrays
use work3_m                                ! Mk3 land-surface diagnostic arrays
use work3f_m                               ! Grid work arrays
use work3sav_m                             ! Water and tracer saved arrays
use workglob_m                             ! Additional grid interpolation
!~ use xarrs_m                                ! Saved dynamic arrays
use xyzinfo_m                              ! Grid coordinate arrays

implicit none

include 'newmpar.h'                        ! Grid parameters
include 'const_phys.h'                     ! Physical constants
include 'darcdf.h'                         ! Netcdf data
include 'dates.h'                          ! Date data
include 'filnames.h'                       ! Filenames
include 'kuocom.h'                         ! Convection parameters
include 'parm.h'                           ! Model configuration
include 'parmdyn.h'                        ! Dynamics parameters
include 'parmgeom.h'                       ! Coordinate data
include 'parmhor.h'                        ! Horizontal advection parameters
include 'parmsurf.h'                       ! Surface parameters
include 'soilv.h'                          ! Soil parameters
include 'stime.h'                          ! File date data
include 'trcom2.h'                         ! Station data
include 'version.h'                        ! Model version data

#ifdef vampir
#include 'vt_user.inc'
#endif
      
integer leap
common/leap_yr/leap                        ! Leap year (1 to allow leap years)
integer nbarewet,nsigmf
common/nsib/nbarewet,nsigmf                ! Land-surface options

integer, dimension(8) :: tvals1, tvals2, nper3hr
#ifdef usempi3
integer, dimension(2) :: shsize
#endif
integer ilx, io_nest, iq, irest, isoil
integer jalbfix, jlx, k, kktau
integer mins_dt, mins_gmt, mspeca, mtimer_in, nalpha
integer nlx, nmaxprsav, npa, npb, n3hr
integer nstagin, nstaguin, nwrite, nwtsav, mins_rad, secs_rad, mtimer_sav
integer nn, i, j, mstn, ierr, nperhr, nversion
integer ierr2, kmax, isoth, nsig, lapsbot, mbd_min
real, dimension(:,:), allocatable, save :: dums, dumliq
real, dimension(:), allocatable, save :: spare1, spare2
real, dimension(:), allocatable, save :: spmean
real, dimension(9) :: temparray, gtemparray
real clhav, cllav, clmav, cltav, dsx, dtds, es
real gke, hourst, hrs_dt, evapavge
real pslavge, pwater, rel_lat, rel_long, spavge, pwatr
real qtot, aa, bb, cc, bb_2, cc_2, rat
real targetlev
real, parameter :: con = 180./pi
character(len=60) comm, comment
character(len=47) header
character(len=10) timeval
character(len=8) rundate
logical odcalc

! version namelist
namelist/defaults/nversion
! main namelist
namelist/cardin/comment,dt,ntau,nwt,npa,npb,nperavg,ia,ib, &
    ja,jb,id,jd,iaero,mex,mbd,nbd,             &
    mbd_maxscale,ndi,ndi2,nlv,nmaxpr,nrad,ntaft,ntsea,ntsur, &
    nvmix,restol,precon,kdate_s,ktime_s,leap,newtop,mup,lgwd,     &
    ngwd,rhsat,nextout,jalbfix,nalpha,nstag,nstagu,ntbar,nwrite,  &
    irest,nrun,nstn,rel_lat,rel_long,nrungcm,nsib,istn,jstn,iunp, &
    slat,slon,zstn,name_stn,mh_bs,nritch_t,nt_adv,mfix,mfix_qg,   &
    namip,amipo3,nh,nhstest,nsemble,nspecial,panfg,panzo,nplens,  &
    rlatdn,rlatdx,rlongdn,rlongdx,newrough,nglacier,newztsea,     &
    epsp,epsu,epsf,epsh,av_vmod,charnock,chn10,snmin,tss_sh,      &
    vmodmin,zobgin,rlong0,rlat0,schmidt,kbotdav,kbotu,nbox,nud_p, &
    nud_q,nud_t,nud_uv,nud_hrs,nudu_hrs,sigramplow,sigramphigh,   &
    nlocal,nbarewet,nsigmf,qgmin,io_in,io_nest,io_out,io_rest,    &
    tblock,tbave,localhist,m_fly,mstn,nqg,nurban,nmr,ktopdav,     &
    nud_sst,nud_sss,mfix_tr,mfix_aero,kbotmlo,ktopmlo,mloalpha,   &
    nud_ouv,nud_sfh,bpyear,rescrn,helmmeth,nmlo,ol, &
    knh,ccycle,kblock,nud_aero,helim,  &
    fc2,sigbot_gwd,alphaj,cgmap_offset,cgmap_scale,nriver
! file namelist
namelist/datafile/ifile,ofile,albfile,co2emfile,eigenv,hfile,     &
    icefile,mesonest,nmifile,o3file,radfile,restfile,rsmfile,     &
    scamfile,scrnfile,snowfile,so4tfile,soilfile,sstfile,surfile, &
    tmaxfile,tminfile,topofile,trcfil,vegfile,zofile,smoistfile,  &
    soil2file,radonemfile,co2_00,radon_00,surf_00,co2_12,         &
    radon_12,surf_12,laifile,albnirfile,urbanfile,bathfile,       &
    vegprev,vegnext,cnsdir,salfile,oxidantfile,casafile,phenfile

data nversion/0/
data comment/' '/,comm/' '/,irest/1/,jalbfix/1/,nalpha/1/
data mins_rad/-1/,nwrite/0/
data lapsbot/0/,io_nest/1/
      
#ifndef stacklimit
! For linux only - removes stacklimit on all processors
call setstacklimit(-1)
#endif

#ifdef i8r8
if ( kind(iq)/=8 .or. kind(es)/=8 ) then
  write(6,*) "ERROR: CCAM configured for double precision"
  stop
end if
#else
if ( kind(iq)/=4 .or. kind(es)/=4 ) then
  write(6,*) "ERROR: CCAM configured for single precision"
  stop
end if
#endif


!--------------------------------------------------------------
! INITALISE MPI ROUTINES
call ccmpi_init


!--------------------------------------------------------------
! INITALISE TIMING LOGS
call log_off()
call log_setup()
call START_LOG(model_begin)


!--------------------------------------------------------------
! READ NAMELISTS AND SET PARAMETER DEFAULTS
ia       = -1   ! diagnostic index
ib       = -1   ! diagnostic index
ntbar    = -1
rel_lat  = 0.
rel_long = 0.
ktau     = 0
ol       = 20   ! default ocean levels
!~ nhor     = -157
!~ nhorps   = -1
!~ khor     = -8
!~ khdif    = 2
!~ nhorjlm  = 1

! All processors read the namelist, so no MPI comms are needed
open(99,file="input",form="formatted",status="old")
read(99, defaults)
if ( myid==0 ) then
  write(6,'(a20," running for nproc =",i7)') version,nproc
  write(6,*) 'Using defaults for nversion = ',nversion
#ifdef usempi3
  write(6,*) 'Using shared memory with number of nodes ',nodecaptian_nproc
#endif
end if

read(99, cardin)
nperday = nint(24.*3600./dt)
nperhr  = nint(3600./dt)
do n3hr = 1,8
  nper3hr(n3hr) = nint(n3hr*3*3600/dt)
end do
if ( nwt==-99 )     nwt = nperday      ! set default nwt to 24 hours
if ( nperavg==-99 ) nperavg = nwt      ! set default nperavg to nwt
if ( nwrite==0 )    nwrite = nperday   ! only used for outfile IEEE
if ( nmlo/=0 .and. abs(nmlo)<=9 ) then
  ol = max( ol, 1 )
else
  ol = 0
end if

if ( abs(nmlo)>=2 ) nriver=1
read(99, datafile)
! try reading boundary layer turbulence namelist

if ( ierr /= 0 ) rewind(99)       ! rewind namelist if turbnml is not found
nagg = 10!max( 10, naero )           ! maximum size of aggregation
nlx        = 0
mtimer_sav = 0


!--------------------------------------------------------------
! READ TOPOGRAPHY FILE TO DEFINE CONFORMAL CUBIC GRID
il_g    = 48 ! default global grid size
rlong0  = 0. ! default longitude
rlat0   = 0. ! default latitude
schmidt = 1. ! default schmidt factor for grid stretching
kl      = 18 ! default number of vertical levels
if ( myid==0 .and. io_in<=4 ) then
  ! open topo file and check its dimensions
  ! here used to supply rlong0,rlat0,schmidt
  ! Remainder of topo file is read in indata.f90
  write(6,*) 'reading topofile header'
  call ccnf_open(topofile,ncidtopo,ierr)
  if ( ierr == 0 ) then
    ! Netcdf format
    lnctopo = 1 ! flag indicating netcdf file
    call ccnf_inq_dimlen(ncidtopo,'longitude',ilx)
    call ccnf_inq_dimlen(ncidtopo,'latitude',jlx)
    call ccnf_get_attg(ncidtopo,'lon0',rlong0)
    call ccnf_get_attg(ncidtopo,'lat0',rlat0)
    call ccnf_get_attg(ncidtopo,'schmidt',schmidt) 
  else
    ! ASCII format      
    lnctopo = 0 ! flag indicating ASCII file
    open(66,file=topofile,recl=2000,status='old',iostat=ierr)
    if ( ierr /= 0 ) then
      write(6,*) "Error opening topofile ",trim(topofile)
      call ccmpi_abort(-1)
    end if
    read(66,*) ilx,jlx,rlong0,rlat0,schmidt,dsx,header
  end if ! (ierr==0) ..else..
  il_g = ilx        
  write(6,*) 'ilx,jlx              ',ilx,jlx
  write(6,*) 'rlong0,rlat0,schmidt ',rlong0,rlat0,schmidt
end if      ! (myid==0.and.io_in<=4)
! store grid dimensions for broadcast below
temparray(1) = rlong0
temparray(2) = rlat0
temparray(3) = schmidt
temparray(4) = real(il_g)


!--------------------------------------------------------------
! READ EIGENV FILE TO DEFINE VERTICAL LEVELS
if ( myid == 0 ) then
  ! Remanded of file is read in indata.f
  open(28,file=eigenv,status='old',form='formatted',iostat=ierr)
  if ( ierr /= 0 ) then
    write(6,*) "Error opening eigenv file ",trim(eigenv)
    call ccmpi_abort(-1)
  end if
  read(28,*)kmax,lapsbot,isoth,nsig
  kl = kmax
  write(6,*)'kl,ol:              ',kl,ol
  write(6,*)'lapsbot,isoth,nsig: ',lapsbot,isoth,nsig
  temparray(5) = real(kl)
  temparray(6) = real(lapsbot)
  temparray(7) = real(isoth)
  temparray(8) = real(nsig)
end if
      
! Broadcast grid data to all processors
! (Since integers are smaller than 1e7, then they can be exactly
!  represented using real*4)
call ccmpi_bcast(temparray(1:8),0,comm_world)
rlong0  = temparray(1)
rlat0   = temparray(2)
schmidt = temparray(3)
il_g    = nint(temparray(4))
kl      = nint(temparray(5))
lapsbot = nint(temparray(6))
isoth   = nint(temparray(7))
nsig    = nint(temparray(8))

      
!--------------------------------------------------------------
! DEFINE newmpar VARIABLES AND DEFAULTS
! CCAM supports face and uniform grid decomposition over processes using preprocessor directives
! Face decomposition reduces MPI message passing, but only works for factors or multiples of six
! processes.  Uniform decomposition is less restrictive on the number of processes, but requires
! more MPI message passing.
#ifdef uniform_decomp
if ( myid == 0 ) then
  write(6,*) "Using uniform grid decomposition"
end if
#else
if ( myid == 0 ) then
  write(6,*) "Using face grid decomposition"
end if
if ( mod(nproc,6)/=0 .and. mod(6,nproc)/=0 ) then
  write(6,*) "ERROR: nproc must be a multiple of 6 or a factor of 6"
  call ccmpi_abort(-1)
end if
#endif
call proctest(npanels,il_g,nproc,nxp,nyp)
jl_g    = il_g + npanels*il_g                 ! size of grid along all panels (usually 6*il_g)
ifull_g = il_g*jl_g                           ! total number of global horizontal grid points
iquad   = 1 + il_g*((8*npanels)/(npanels+4))  ! grid size for interpolation calculations
il      = il_g/nxp                            ! local grid size on process in X direction
jl      = jl_g/nyp                            ! local grid size on process in Y direction
ifull   = il*jl                               ! total number of local horizontal grid points
! The perimeter of the processor region has length 2*(il+jl).
! The first row has 8 possible corner points per panel and the 
! second has 16. In practice these are not all distinct so there could
! be some optimisation.
#ifdef uniform_decomp
npan   = npanels + 1              ! number of panels on this process
iextra = (4*(il+jl)+24)*npan      ! size of halo for MPI message passing
#else      
npan   = max(1,(npanels+1)/nproc) ! number of panels on this process
iextra = 4*(il+jl) + 24*npan      ! size of halo for MPI message passing
#endif
! nrows_rad is a subgrid decomposition for radiation routines
nrows_rad = jl/6
do while( mod(jl, nrows_rad) /= 0 )
  nrows_rad = nrows_rad - 1
end do
if ( myid==0 ) then
  write(6,*) "il_g,jl_g,il,jl   ",il_g,jl_g,il,jl
  write(6,*) "nxp,nyp,nrows_rad ",nxp,nyp,nrows_rad
end if

! some default values for unspecified parameters
if ( ia<0 ) ia = il/2
if ( ib<0 ) ib = ia + 3
if ( ldr==0 ) mbase = 0
dsig4 = max(dsig2+.01, dsig4)
if( mbd/=0 .and. nbd/=0 ) then
  write(6,*) 'setting nbd=0 because mbd/=0'
  nbd = 0
endif
mbd_min = int(20.*112.*90.*schmidt/real(mbd_maxscale))
if ( mbd<mbd_min .and. mbd/=0 ) then
  if ( myid==0 ) then
    write(6,*) "Increasing mbd to satisfy mbd_maxscale ",mbd_maxscale
    write(6,*) "Original mbd and final mbd = ",mbd,mbd_min
  end if
  mbd = mbd_min
end if
nud_hrs = abs(nud_hrs)  ! just for people with old -ves in namelist
if ( nudu_hrs==0 ) nudu_hrs = nud_hrs

! **** do namelist fixes above this ***
      
!--------------------------------------------------------------
! DISPLAY NAMELIST

nstagin  = nstag    ! -ve nstagin gives swapping & its frequency
nstaguin = nstagu   ! only the sign of nstaguin matters (chooses scheme)
if ( nstagin==5 .or. nstagin<0 ) then
  nstag  = 4
  nstagu = 4
  if ( nstagin == 5 ) then  ! for backward compatability
    nstagin  = -1 
    nstaguin = 5  
  endif
endif
if ( kblock<0 ) then
  kblock = max(kl, ol) ! must occur before indata
end if

!~ tke_umin = vmodmin

!--------------------------------------------------------------
! INITIALISE ifull_g ALLOCATABLE ARRAYS
#ifdef usempi3
! Allocate xx4, yy4, em_g, x_g, y_g and z_g as shared
! memory within a node.  The node captian is responsible
! for updating these arrays.
shsize(1:2) = (/ iquad, iquad /)
call ccmpi_allocshdatar8(xx4,shsize(1:2),xx4_win)
call ccmpi_allocshdatar8(yy4,shsize(1:2),yy4_win)
shsize(1) = ifull_g
call ccmpi_allocshdata(em_g,shsize(1:1),em_g_win)
call ccmpi_allocshdatar8(x_g,shsize(1:1),x_g_win)
call ccmpi_allocshdatar8(y_g,shsize(1:1),y_g_win)
call ccmpi_allocshdatar8(z_g,shsize(1:1),z_g_win)
#else
! Allocate xx4, yy4, em_g, x_g, y_g and z_g for
! each process
allocate( xx4(iquad,iquad), yy4(iquad,iquad) )
allocate( em_g(ifull_g) )
allocate( x_g(ifull_g), y_g(ifull_g), z_g(ifull_g) )
#endif
call xyzinfo_init(ifull_g,ifull,iextra,myid)
call indices_init(ifull_g,ifull,iextra,npanels,npan)
call map_init(ifull_g,ifull,iextra,myid)
call latlong_init(ifull_g,ifull,iextra,myid)      
call vecsuv_init(ifull_g,ifull,iextra,myid)


!--------------------------------------------------------------
! SET UP CC GEOMETRY
! Only one process calls setxyz to save memory with large grids
if ( myid==0 ) then
  write(6,*) "Calling setxyz"
  call workglob_init(ifull_g)
  call setxyz(il_g,rlong0,rlat0,schmidt,x_g,y_g,z_g,wts_g,ax_g,ay_g,az_g,bx_g,by_g,bz_g,xx4,yy4)
end if
! Broadcast the following global data
! xx4 and yy4 are used for calculating depature points
! em_g, x_g, y_g and z_g are for the scale-selective filter (1D and 2D versions)
#ifdef usempi3
call ccmpi_shepoch(xx4_win) ! also yy4_win, em_g_win, x_g_win, y_g_win, z_g_win
if ( node_myid==0 ) then
  call ccmpi_bcastr8(xx4,0,comm_nodecaptian)
  call ccmpi_bcastr8(yy4,0,comm_nodecaptian)
  call ccmpi_bcast(em_g,0,comm_nodecaptian)
  call ccmpi_bcastr8(x_g,0,comm_nodecaptian)
  call ccmpi_bcastr8(y_g,0,comm_nodecaptian)
  call ccmpi_bcastr8(z_g,0,comm_nodecaptian)
end if
call ccmpi_shepoch(xx4_win) ! also yy4_win, em_g_win, x_g_win, y_g_win, z_g_win
#else
call ccmpi_bcastr8(xx4,0,comm_world)
call ccmpi_bcastr8(yy4,0,comm_world)
call ccmpi_bcast(em_g,0,comm_world)
call ccmpi_bcastr8(x_g,0,comm_world)
call ccmpi_bcastr8(y_g,0,comm_world)
call ccmpi_bcastr8(z_g,0,comm_world)
#endif
call ccmpi_bcast(ds,0,comm_world)

if ( myid==0 ) then
  write(6,*) "Calling ccmpi_setup"
end if
call ccmpi_setup(kblock)

      
!--------------------------------------------------------------
! DEALLOCATE ifull_g ARRAYS WHERE POSSIBLE
call worklocl_init(ifull)      
if ( myid==0 ) then
  call ccmpi_distribute(rlong4_l,rlong4)
  call ccmpi_distribute(rlat4_l,rlat4)
  call workglob_end
  deallocate( wts_g, emu_g, emv_g )
  deallocate( ax_g, ay_g, az_g )
  deallocate( bx_g, by_g, bz_g )
  deallocate( f_g, fu_g, fv_g )
  deallocate( dmdx_g, dmdy_g )
  deallocate( rlatt_g, rlongg_g )
else
  call ccmpi_distribute(rlong4_l)
  call ccmpi_distribute(rlat4_l)
end if


!--------------------------------------------------------------
! INITIALISE LOCAL ARRAYS
call arrays_init(ifull,iextra,kl)
call carbpools_init(ifull,iextra,kl,nsib,ccycle)
call extraout_init(ifull,iextra,kl,nextout)
call histave_init(ifull,iextra,kl,ms)
call morepbl_init(ifull,iextra,kl)
call nsibd_init(ifull,iextra,kl,nsib)
call pbl_init(ifull,iextra,kl)
call permsurf_init(ifull,iextra,kl)
call savuvt_init(ifull,iextra,kl)
call savuv1_init(ifull,iextra,kl)
call sbar_init(ifull,iextra,kl)
call sigs_init(ifull,iextra,kl)
call soil_init(ifull,iextra,kl,iaero,nsib)
call soilsnow_init(ifull,iextra,kl,ms,nsib)
call vecs_init(ifull,iextra,kl)
call vegpar_init(ifull,iextra,kl)
call work2_init(ifull,iextra,kl,nsib)
call work3_init(ifull,iextra,kl,nsib)
call work3f_init(ifull,iextra,kl)
!~ call xarrs_init(ifull,iextra,kl)

! Remaining arrays are allocated in indata.f90, since their
! definition requires additional input data (e.g, land-surface)

!--------------------------------------------------------------
! DISPLAY DIAGNOSTIC INDEX AND TIMER DATA
if ( mydiag ) then
  write(6,"('id,jd,rlongg,rlatt in degrees: ',2i4,2f8.2)") id,jd,con*rlongg(idjd),con*rlatt(idjd)
end if
call date_and_time(rundate)
call date_and_time(time=timeval)
if ( myid == 0 ) then
  write(6,*)'RUNDATE IS ',rundate
  write(6,*)'Starting time ',timeval
end if


!--------------------------------------------------------------
! READ INITIAL CONDITIONS
if ( myid == 0 ) then
  write(6,*) "Calling indata"
end if
call indataf(hourst,jalbfix,lapsbot,isoth,nsig,io_nest)

!~ if ( mins_rad<0 ) then
  !~ ! automatic estimate for mins_rad
  !~ secs_rad = min(nint((schmidt*112.*90./real(il_g))*8.*60.), 3600)
  !~ secs_rad = min(secs_rad, nint(real(nwt)*dt))
  !~ secs_rad = max(secs_rad, 1)
  !~ kountr   = nint(real(secs_rad)/dt)
  !~ secs_rad = nint(real(kountr)*dt)
  !~ do while ( mod(3600, secs_rad)/=0 .and. kountr>1 )
    !~ kountr = kountr - 1
    !~ secs_rad = nint(real(kountr)*dt)
  !~ end do
!~ else
  !~ ! user specified mins_rad
  !~ kountr   = nint(real(mins_rad)*60./dt)  ! set default radiation to ~mins_rad m
  !~ secs_rad = nint(real(kountr)*dt)        ! redefine to actual value
!~ end if

!--------------------------------------------------------------
! NRUN COUNTER
if ( myid == 0 ) then
  open(11, file='nrun.dat',status='unknown')
  if ( nrun == 0 ) then
    read(11,*,iostat=ierr2) nrun
    nrun = nrun + 1
  endif                  ! nrun==0
  write(6,*)'this is run ',nrun
  rewind 11
  write(11,*) nrun
  write(11,cardin)
  write(11,datafile)
  close(11)
end if


!-------------------------------------------------------------
! SETUP DIAGNOSTIC ARRAYS
runoff(:)      = 0.
wb_ave(:,:)    = 0.

!--------------------------------------------------------------
! OPEN OUTPUT FILES AND SAVE INITAL CONDITIONS
if ( nwt > 0 ) then
  ! write out the first ofile data set
  if ( myid == 0 ) then
    write(6,*)'calling outfile'
  end if
  call outfile(20,rundate,nwrite,nstagin,jalbfix,nalpha,mins_rad)  ! which calls outcdf
end if    ! (nwt>0)


!--------------------------------------------------------------
! INITIALISE DYNAMICS
dtin = dt
n3hr = 1   ! initial value at start of run
if ( myid == 0 ) then
  write(6,*) "number of time steps per day = ",nperday
  write(6,*) "nper3hr,nper6hr .. ",nper3hr(:)
end if
mspeca = 1
if ( mex/=1 .and. .not.lrestart ) then
  mspeca = 2
  dt     = dtin*.5
endif
call gettin(0)              ! preserve initial mass & T fields

nmaxprsav = nmaxpr
nwtsav    = nwt
hrs_dt    = dtin/3600.      ! time step in hours
mins_dt   = nint(dtin/60.)  ! time step in minutes
mtimer_in = mtimer
 
 
!--------------------------------------------------------------
! BEGIN MAIN TIME LOOP
if ( myid == 0 ) then
  call date_and_time(time=timeval,values=tvals1)
  write(6,*) "Start of loop time ", timeval
end if
call log_on()
call START_LOG(maincalc_begin)

do kktau = 1,ntau   ! ****** start of main time loop

  ktau     = kktau
  timer    = timer + hrs_dt                      ! timer now only used to give timeg
  timeg    = mod(timer+hourst,24.)
  mtimer   = mtimer_in + nint(ktau*dtin/60.)     ! 15/6/01 to allow dt < 1 minute
  mins_gmt = mod(mtimer+60*ktime/100,24*60)

  
  ! ***********************************************************************
  ! START OCEAN DYNAMICS
  ! ***********************************************************************
  
  ! nriver=1 allows the rivers to work without the ocean model

    ! RIVER ROUTING ------------------------------------------------------
    call START_LOG(river_begin)
    if ( nmaxpr==1 ) then
      if ( myid==0 ) then
        write(6,*) "Before river"
      end if
      call ccmpi_barrier(comm_world)
    end if
    call rvrrouter
    if ( nmaxpr==1 ) then
      if ( myid==0 ) then
        write(6,*) "After river"
      end if
      call ccmpi_barrier(comm_world)
    end if
    call END_LOG(river_end)

  ! SURFACE FLUXES ---------------------------------------------
  ! (Includes ocean dynamics and mixing, as well as ice dynamics and thermodynamics)
  call START_LOG(sfluxnet_begin)

  if ( ntsur > 1 ) then  ! should be better after convjlm
    call sflux(nalpha)
  endif   ! (ntsur>1)    

  call END_LOG(sfluxnet_end)
  
  ! ***********************************************************************
  ! DIAGNOSTICS AND OUTPUT
  ! ***********************************************************************
  ! DIAGNOSTICS ------------------------------------------------

  if ( mod(ktau,nmaxpr)==0 .or. ktau==ntau ) then
    call maxmin(wb,'wb',ktau,1.,ms)
  endif                  ! (mod(ktau,nmaxpr)==0)

  ! update diag_averages and daily max and min screen temps 
  ! N.B. runoff is accumulated in sflux
  wb_ave(1:ifull,1:ms) = wb_ave(1:ifull,1:ms) + wb


  if ( ktau==ntau .or. mod(ktau,nperavg)==0 ) then
    do k=1,ms
      wb_ave(1:ifull,k) = wb_ave(1:ifull,k)/min(ntau,nperavg)
    end do
  end if    ! (ktau==ntau.or.mod(ktau,nperavg)==0)

  call log_off()
  if ( ktau==ntau .or. mod(ktau,nwt)==0 ) then
    call outfile(20,rundate,nwrite,nstagin,jalbfix,nalpha,mins_rad)  ! which calls outcdf

    if ( ktau==ntau .and. irest==1 ) then
      ! Don't include the time for writing the restart file
      call END_LOG(maincalc_end)
      ! write restart file
      call outfile(19,rundate,nwrite,nstagin,jalbfix,nalpha,mins_rad)
      if ( myid == 0 ) then
        write(6,*)'finished writing restart file in outfile'
      end if
      call START_LOG(maincalc_begin)
    endif  ! (ktau==ntau.and.irest==1)
  endif    ! (ktau==ntau.or.mod(ktau,nwt)==0)
  call log_on()
 
  if ( mod(ktau,nperavg) == 0 ) then   
    wb_ave(:,:)    = 0.
    runoff(:)      = 0.  ! converted to mm/day in outcdf
  endif  ! (mod(ktau,nperavg)==0)

  if ( mod(ktau,nperday) == 0 ) then   ! re-set at the end of each 24 hours
    if ( ntau<10*nperday .and. nstn>0 ) then     ! print stn info
      do nn = 1,nstn
        if ( .not.mystn(nn) ) cycle
        i = istn(nn)
        j = jstn(nn)
        iq = i+(j-1)*il
      end do
    end if  ! (ntau<10*nperday)
  endif   ! (mod(ktau,nperday)==0)
  
#ifdef vampir
  ! Flush vampir trace information to disk to save memory.
  VT_BUFFER_FLUSH()
#endif

end do                  ! *** end of main time loop
call END_LOG(maincalc_end)
call log_off()

! Report timings of run
if ( myid == 0 ) then
  call date_and_time(time=timeval,values=tvals2)
  write(6,*) "End of time loop ", timeval
  write(6,*) "normal termination of run"
  call date_and_time(time=timeval)
  write(6,*) "End time ", timeval
  aa = 3600.*(tvals2(5)-tvals1(5)) + 60.*(tvals2(6)-tvals1(6)) + (tvals2(7)-tvals1(7)) + 0.001*(tvals2(8)-tvals1(8))
  if ( aa <= 0. ) aa = aa + 86400.
  write(6,*) "Model time in main loop",aa
end if
call END_LOG(model_end)

! close mesonest files
if ( mbd/=0 .or. nbd/=0 ) then
  call histclose
end if

#ifdef simple_timer
! report subroutine timings
call simple_timer_finalize
#endif

! finalize MPI comms
call ccmpi_finalize

end


!--------------------------------------------------------------
! INITIAL PARAMETERS
blockdata main_blockdata

implicit none

include 'newmpar.h'          ! Grid parameters
include 'dates.h'            ! Date data
include 'filnames.h'         ! Filenames
include 'kuocom.h'           ! Convection parameters
include 'parm.h'             ! Model configuration
include 'parmdyn.h'          ! Dynamics parmaters
include 'parmgeom.h'         ! Coordinate data
include 'parmhor.h'          ! Horizontal advection parameters
include 'parmsurf.h'         ! Surface parameters
include 'soilv.h'            ! Soil parameters
include 'stime.h'            ! File date data
include 'trcom2.h'           ! Station data

integer leap
common/leap_yr/leap          ! Leap year (1 to allow leap years)
integer nbarewet,nsigmf
common/nsib/nbarewet,nsigmf  ! Land-surface options

! for cardin
data ia/1/,ib/3/,id/2/,ja/1/,jb/10/,jd/5/,nlv/1/
data ndi/1/,ndi2/0/,nmaxpr/99/     
data kdate_s/-1/,ktime_s/-1/,leap/0/
data mbd/0/,mbd_maxscale/3000/,nbd/0/,nbox/1/,kbotdav/4/,kbotu/0/
data nud_p/0/,nud_q/0/,nud_t/0/,nud_uv/1/,nud_hrs/24/,nudu_hrs/0/
data ktopdav/0/,kblock/-1/
data nud_aero/0/
data nud_sst/0/,nud_sss/0/,nud_ouv/0/,nud_sfh/0/
data mloalpha/10/,kbotmlo/-1000/,ktopmlo/1/
data sigramplow/0./,sigramphigh/0./
! Dynamics options A & B      
data mex/30/,mfix/3/,mfix_qg/1/,mup/1/,nh/0/
data nritch_t/300/,epsp/-15./,epsu/0./,epsf/0./,epsh/0.1/
data precon/-2900/,restol/4.e-7/
data schmidt/1./,rlong0/0./,rlat0/90./,nrun/0/
data helmmeth/0/,mfix_tr/0/,mfix_aero/0/
! Horiz advection options
data nt_adv/7/,mh_bs/4/
! Horiz wind staggering options
data nstag/-10/,nstagu/-1/,nstagoff/0/
! Vertical mixing options
data nvmix/3/,nlocal/6/,ncvmix/0/,lgwd/0/,ngwd/-5/
data helim/800./,fc2/1./,sigbot_gwd/0./,alphaj/1.e-6/
data cgmap_offset/0./,cgmap_scale/1./
! Soil, canopy, PBL options
data nbarewet/0/,newrough/0/,nglacier/1/
data nrungcm/-1/,nsib/3/,nsigmf/1/
data ntaft/2/,ntsea/6/,ntsur/6/,av_vmod/.7/,tss_sh/1./
data vmodmin/.2/,zobgin/.02/,charnock/.018/,chn10/.00125/
data newztsea/1/,newtop/1/                
data snmin/.11/  ! 1000. for 1-layer; ~.11 to turn on 3-layer snow
data nurban/0/,ccycle/0/
! I/O options
data m_fly/4/,io_in/1/,io_out/1/,io_rest/1/
data nperavg/-99/,nwt/-99/,tblock/1/,tbave/1/
data nextout/3/,localhist/.false./
data nstn/0/  
data slat/nstnmax*-89./,slon/nstnmax*0./,iunp/nstnmax*6/
data zstn/nstnmax*0./,name_stn/nstnmax*'   '/ 
 
! initialize file names to something
data albfile/' '/,icefile/' '/,maskfile/' '/
data snowfile/' '/,sstfile/' '/,topofile/' '/,zofile/' '/
data rsmfile/' '/,scamfile/' '/,soilfile/' '/,vegfile/' '/
data co2emfile/' '/,so4tfile/' '/
data smoistfile/' '/,soil2file/' '/,restfile/' '/
data radonemfile/' '/,surfile/' '/,surf_00/'s_00a '/
data surf_12/'s_12a '/,co2_00/' '/,co2_12/' '/,radon_00/' '/
data radon_12/' '/,ifile/' '/,ofile/' '/,nmifile/' '/
data eigenv/' '/,radfile/' '/,o3file/' '/,hfile/' '/,mesonest/' '/
data scrnfile/' '/,tmaxfile/' '/,tminfile/' '/,trcfil/' '/
data laifile/' '/,albnirfile/' '/,urbanfile/' '/,bathfile/' '/
data vegprev/' '/,vegnext/' '/,cnsdir/' '/,salfile/' '/
data oxidantfile/' '/,casafile/' '/,phenfile/' '/
! floating point:
data timer/0./,mtimer/0/

! stuff from insoil  for soilv.h
data rlaim44/4.8, 6.3, 5., 3.75, 2.78, 2.5, 3.9, 2.77, 2.04, 2.6,         & ! 1-10
             1.69, 1.9, 1.37, 1.5, 1.21, 1.58, 1.41, 2.3, 1.2, 1.71,      & ! 11-20
             1.21, 2.3, 2.3, 1.2, 1.2, 1.87, 1., 3., .01, .01, 1.2,       & ! 21-31
             6., 5.5, 5., 4.5, 5., 4., 3., 3.5, 1., 4., .5, 4., 0./         ! 32-44
data rlais44/1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,                      & ! 1-10
             1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,                      & ! 11-20
             1., 1., 1., 1., .6, .6, .5, 1., 0., 0., 1.,                  & ! 21-31
             2., 2., 2., 2., 2., 1.5, 1.5, 1.5, 1., .5, .5, .5, 0./         ! 32-44
data rsunc44/370., 330., 260., 200., 150., 130., 200., 150., 110., 160.,  & ! 1-10
             100., 120.,  90.,  90.,  80.,  90.,  90., 150.,  80., 100.,  & ! 11-20
             80.,  80.,  80.,  60.,  60., 120.,  80., 180., 2*995., 80.,  & ! 21-31
             350., 4*300., 3*230., 150., 230., 995., 150., 9900./           ! 32-44
data scveg44/0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,                      & ! 1-10
             0., 0., 0., 0., 0., .1, .1, .1, .1, .1,                      & ! 11-20
             .1, .2, .4, .2, .1, .1, .1, 0., 0., 0., 0.,                  & ! 21-31
             .05, 0., 0., 0., 0., .05, .05, .05, .1, 0., 0., .4, 0./        ! 32-44
data slveg44/0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,                      & ! 1-10
             0., 0., 0., 0., 0., .1, .1, .1, .1, .1,                      & ! 11-20
             .1, .2, .4, .2, .1, .1, .1, 0., 0., 0., 0.,                  & ! 21-31
             1., 5.5, 3., 1., 3., 3., 3.5, 3., .5, 3.5, .1, 3.5, 0./        ! 32-44
data froot/.05, .10, .35, .40, .10/       ! 10/02/99 veg. root distr.

data silt/.08, .33, .17, .2, .06, .25, .15, .70, .33, .2, .33, .33, .17/    ! with mxst=13
data clay/.09, .3, .67, .2, .42, .48, .27, .17, .30, .2, .3, .3, .67/       ! with mxst=13
data sand/.83, .37, .16, .6, .52, .27, .58, .13, .37, .6, .37, .37, .17/    ! with mxst=13
data swilt/0., .072, .216, .286, .135, .219, .283, .175, .395, .216, .1142, .1547, .2864, .2498/
data sfc/1.,  .143, .301, .367, .218, .31 , .37 , .255, .45, .301, .22 , .25 , .367, .294/
data ssat/2., .398, .479, .482, .443, .426, .482, .420, .451, .479, .435, .451, .482, .476/
data bch/4.2, 7.1, 11.4, 5.15, 10.4, 10.4, 7.12, 5.83, 7.1, 4.9, 5.39, 11.4, 8.52/ ! bch for gravity term
data hyds/166.e-6, 4.e-6, 1.e-6, 21.e-6, 2.e-6, 1.e-6, 6.e-6,800.e-6, 1.e-6, 34.e-6, 7.e-6, 1.3e-6, 2.5e-6/
data sucs/-.106, -.591, -.405, -.348, -.153, -.49, -.299,-.356, -.153, -.218, -.478, -.405, -.63/ ! phisat (m)
data rhos/7*2600., 1300.,  910., 4*2600. /     ! soil density
data  css/7* 850., 1920., 2100., 4*850./       ! heat capacity

data zse/.022, .058, .154, .409, 1.085, 2.872/ ! layer thickness
! so depths of centre of layers: .011, .051, .157, .4385, 1.1855, 3.164
! with base at 4.6     

end
      
!--------------------------------------------------------------
! TEST GRID DECOMPOSITION    
subroutine proctest(npanels,il_g,nproc,nxp,nyp)
implicit none
integer, intent(in) :: il_g, nproc, npanels
integer, intent(out) :: nxp, nyp
integer jl_g

#ifdef uniform_decomp
jl_g = il_g + npanels*il_g     ! size of grid along all panels (usually 6*il_g)
nxp = nint(sqrt(real(nproc)))  ! number of processes in X direction
nyp = nproc/nxp                ! number of processes in Y direction
! search for vaild process decomposition.  CCAM enforces the same grid size on each process
do while ( (mod(il_g,max(nxp,1))/=0.or.mod(nproc,max(nxp,1))/=0.or.mod(il_g,nyp)/=0) .and. nxp>0 )
  nxp = nxp - 1
  nyp = nproc/max(nxp,1)
end do
#else
if ( mod(nproc,6)/=0 .and. mod(6,nproc)/=0 ) then
  nxp = -1
else
  jl_g = il_g + npanels*il_g                 ! size of grid along all panels (usually 6*il_g)
  nxp = max( 1, nint(sqrt(real(nproc)/6.)) ) ! number of processes in X direction
  nyp = nproc/nxp                            ! number of processes in Y direction
  ! search for valid process decomposition.  CCAM enforces the same grid size on each process
  do while ( (mod(il_g,max(nxp,1))/=0.or.mod(nproc/6,max(nxp,1))/=0.or.mod(jl_g,max(nyp,1))/=0) .and. nxp>0 )
    nxp = nxp - 1
    nyp = nproc/max(nxp,1)
  end do
end if
#endif

return
end subroutine proctest
    