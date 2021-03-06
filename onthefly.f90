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

! Main NetCDF input routines.  Host grid is automatically
! interpolated to nested model grid.  Three options are
!   nested=0  Initial conditions
!   nested=1  Nudging fields
!   nested=2  Surface data recycling
      
! This version supports the parallel file routines contained
! in infile.f90.  Hence, restart files do not require any
! gathers and scatters.

! In the case where the grid needs to be interpolated, RMA
! is used to distribute host data from processes, which
! reduces the amount of message passing.
    
! When -Dusempi3 is enabled, then host data arrays are
! shared between processes on a node.  The node captian
! is then responsible for obtaining interpolation data
! for all processes on a node.
    
! Thanks to Paul Ryan for advice on input NetCDF routines
    
module onthefly_m
    
implicit none

private
public onthefly
    
integer, parameter :: nord = 3                                ! 1 for bilinear, 3 for bicubic interpolation
integer, save :: ik, jk, kk, ok, nsibx                        ! input grid size
integer dk, fwsize                                            ! size of temporary arrays
integer, dimension(:,:), allocatable, save :: nface4          ! interpolation panel index
integer, dimension(0:5), save :: comm_face                    ! communicator for processes requiring a input panel
real, save :: rlong0x, rlat0x, schmidtx                       ! input grid coordinates
real, dimension(3,3), save :: rotpoles, rotpole               ! vector rotation data
real, dimension(:,:), allocatable, save :: xg4, yg4           ! interpolation coordinate indices
real, dimension(:), allocatable, save :: axs_a, ays_a, azs_a  ! vector rotation data
real, dimension(:), allocatable, save :: bxs_a, bys_a, bzs_a  ! vector rotation data 
real, dimension(:), allocatable, save :: axs_w, ays_w, azs_w  ! vector rotation data
real, dimension(:), allocatable, save :: bxs_w, bys_w, bzs_w  ! vector rotation data
!~ real, dimension(:), allocatable, save :: sigin                ! input vertical coordinates
logical iotest, newfile                                       ! tests for interpolation and new metadata
logical, dimension(0:5), save :: nfacereq = .false.           ! list of panels required for interpolation
logical, save :: bcst_allocated = .false.                     ! Bcast communicator groups have been defined

contains

! *****************************************************************************
! Main interface for input data that reads grid metadata

  


subroutine onthefly(nested,kdate_r,ktime_r,psl,zss,tss,sicedep,fracice,t,u,v,qg,tgg,wb,wbice,snowd, &
                    tggsn,smass,ssdn,ssdnn,snage,isflag)

use cc_mpi           ! CC MPI routines
use infile           ! Input file routines
use soil_m           ! Soil and surface data

implicit none

include 'newmpar.h'  ! Grid parameters
include 'darcdf.h'   ! Netcdf data
include 'parm.h'     ! Model configuration
include 'stime.h'    ! File date data

integer, parameter :: nihead = 54
integer, parameter :: nrhead = 14

integer, intent(in) :: nested
integer, intent(out) :: kdate_r, ktime_r
integer, save :: maxarchi
integer mtimer, ierx, idvkd, idvkt, idvmt
integer, dimension(nihead) :: nahead
integer, dimension(ifull), intent(out) :: isflag
real timer
real, dimension(ifull,ms), intent(out) :: wb, wbice, tgg
real, dimension(ifull,3), intent(out) :: tggsn, smass, ssdn
real, dimension(:,:), intent(out) :: t, u, v, qg
real, dimension(ifull), intent(out) :: psl, zss, tss, fracice, snowd
real, dimension(ifull), intent(out) :: sicedep, ssdnn, snage
real, dimension(nrhead) :: ahead
real, dimension(10) :: rdum
logical ltest, tst

call START_LOG(onthefly_begin)
!--------------------------------------------------------------------
! pfall indicates all processors have a parallel input file and there
! is no need to broadcast metadata (see infile.f90).  Otherwise read
! metadata on myid=0 and broadcast that data to all processors.
if ( myid==0 .or. pfall ) then
  if ( myid==0 ) write(6,*) 'Entering onthefly for nested,ktau = ',nested,ktau
  
  ! Locate new file and read grid metadata --------------------------
  if ( ncid/=ncidold ) then
    if ( myid==0 ) write(6,*) 'Reading new file metadata'
    iarchi=1   ! default time index for input file
    maxarchi=0 ! default number of timesteps in input file
    call ccnf_get_attg(ncid,'int_header',nahead)
    call ccnf_get_attg(ncid,'real_header',ahead)
    ik      =pil_g      ! grid size
    jk      =pjl_g      ! grid size
    kk      =pka_g      ! atmosphere vertical levels
    ok      =pko_g      ! ocean vertical levels
    nsibx   =nahead(44) ! land-surface parameterisation
    rlong0x =ahead(5)   ! longitude
    rlat0x  =ahead(6)   ! latitude
    schmidtx=ahead(7)   ! schmidt factor
    if ( schmidtx<=0. .or. schmidtx>1. ) then
      ! backwards compatibility option
      rlong0x =ahead(6)
      rlat0x  =ahead(7)
      schmidtx=ahead(8)
    endif  ! (schmidtx<=0..or.schmidtx>1.)        
    call ccnf_inq_dimlen(ncid,'time',maxarchi)
    if ( myid==0 ) then
      write(6,*) "Found ik,jk,kk,ok ",ik,jk,kk,ok
      write(6,*) "      maxarchi ",maxarchi
      write(6,*) "      rlong0x,rlat0x,schmidtx ",rlong0x,rlat0x,schmidtx
    end if
  end if
  
  ! search for required date ----------------------------------------
  if ( myid==0 ) write(6,*)'Search for kdate_s,ktime_s >= ',kdate_s,ktime_s
  ltest = .true.       ! flag indicates that the date is not yet found
  iarchi = iarchi - 1  ! move time index back one step to check current position in file
  ierx = 0             ! indicates normal mtimer format or backwards compatibility mode
  call ccnf_inq_varid(ncid,'kdate',idvkd,tst)
  call ccnf_inq_varid(ncid,'ktime',idvkt,tst)
  call ccnf_inq_varid(ncid,'mtimer',idvmt,tst)
  if ( tst ) then
    ! backwards compatability option
    ierx = 1
    call ccnf_inq_varid(ncid,'timer',idvmt,tst)
  end if
  ! start search for required date/time
  do while( ltest .and. iarchi<maxarchi )
    ! could read this as one array, but we only usually need to advance 1 step
    iarchi = iarchi + 1
    call ccnf_get_vara(ncid,idvkd,iarchi,kdate_r)
    call ccnf_get_vara(ncid,idvkt,iarchi,ktime_r)
    if ( ierx==0 ) then
      call ccnf_get_vara(ncid,idvmt,iarchi,mtimer)
      timer = mtimer/60.
    else
      timer = 0.
      call ccnf_get_vara(ncid,idvmt,iarchi,timer)
      mtimer = nint(timer*60.)
    endif
    if ( mtimer>0 ) then
      ! calculate date if mtimer>0
      call datefix(kdate_r,ktime_r,mtimer)
    end if
    ! ltest = .false. when correct date is found
    ltest = (2400*(kdate_r-kdate_s)-1200*nsemble+(ktime_r-ktime_s))<0
  end do
  if ( nsemble/=0 ) then
    kdate_r = kdate_s
    ktime_r = ktime_s
  end if
  if ( ltest ) then
    ! ran out of file before correct date was located
    ktime_r = -1
  end if
  if ( myid==0 ) then
    write(6,*) 'After search ltest,iarchi =',ltest,iarchi
    write(6,*) '             kdate_r,ktime_r =',kdate_r,ktime_r
  end if

endif  ! ( myid==0 .or. pfall )

! if metadata is not read by all processors, then broadcast ---------
if ( .not.pfall ) then
  rdum(1) = rlong0x
  rdum(2) = rlat0x
  rdum(3) = schmidtx
  ! kdate_r is too large to represent as a single real, so
  ! we split kdate_r into year, month and day
  rdum(4) = real(kdate_r/10000)
  rdum(5) = real(kdate_r/100-nint(rdum(4))*100)
  rdum(6) = real(kdate_r-nint(rdum(4))*10000-nint(rdum(5))*100)
  rdum(7) = real(ktime_r)
  if ( ncid/=ncidold ) then
    rdum(8) = 1.
  else
    rdum(8) = 0.
  end if
  rdum(9)  = real(iarchi)
  rdum(10) = real(nsibx)
  call ccmpi_bcast(rdum(1:10),0,comm_world)
  rlong0x  = rdum(1)
  rlat0x   = rdum(2)
  schmidtx = rdum(3)
  kdate_r  = nint(rdum(4))*10000+nint(rdum(5))*100+nint(rdum(6))
  ktime_r  = nint(rdum(7))
  newfile  = (nint(rdum(8))==1)
  iarchi   = nint(rdum(9))
  nsibx    = nint(rdum(10))
  ik       = pil_g      ! grid size
  jk       = pjl_g      ! grid size
  kk       = pka_g      ! atmosphere vertical levels
  ok       = pko_g      ! ocean vertical levels
else
  newfile = (ncid/=ncidold)
end if

! mark current file as read for metadata
ncidold = ncid

! trap error if correct date/time is not located --------------------
if ( ktime_r<0 ) then
  if ( nested==2 ) then
    if ( myid==0 ) then
      write(6,*) "WARN: Cannot locate date/time in input file"
    end if
    return
  else
    write(6,*) "ERROR: Cannot locate date/time in input file"
    call ccmpi_abort(-1)
  end if
end if
!--------------------------------------------------------------------
      
! Here we call ontheflyx with different automatic array
! sizes.  dk is used for global arrays that are defined
! on myid==0.  fwsize is used for MPI RMA in loading 
! files over multiple processors.
      
! Note that if histrd fails to find a variable, it returns
! zero in the output array
      
if ( myid==0 ) then
  dk = ik ! non-zero automatic array size in onthefly_work
else
  dk = 0  ! zero automatic array size in onthefly_work
end if

! memory needed to read input files
fwsize = pil*pjl*pnpan*mynproc 

call onthefly_work(nested,kdate_r,ktime_r,psl,zss,tss,sicedep,fracice,t,u,v,qg,tgg,wb,wbice, &
                   snowd,tggsn,smass,ssdn,ssdnn,snage,isflag)

if ( myid==0 ) write(6,*) "Leaving onthefly"

call END_LOG(onthefly_end)

return
                    end subroutine onthefly


! *****************************************************************************
! Read data from netcdf file
      
! Input usually consists of either a single input file that is
! scattered across processes, or multiple input files that are
! read by many processes and shared by RMA.  In the case of
! restart files, then there is no need for message passing.


subroutine onthefly_work(nested,kdate_r,ktime_r,psl,zss,tss,sicedep,fracice,t,u,v,qg,tgg,wb,wbice, &
                         snowd,tggsn,smass,ssdn,ssdnn,snage,isflag)

use cable_def_types_mod, only : ncs, ncp       ! CABLE dimensions
use cc_mpi                                     ! CC MPI routines
use infile                                     ! Input file routines
use latlong_m                                  ! Lat/lon coordinates
use morepbl_m                                  ! Additional boundary layer diagnostics
use nsibd_m, only : isoilm                     ! Land-surface arrays
use river                                      ! River routing
!~ use sigs_m                                     ! Atmosphere sigma levels
use soil_m                                     ! Soil and surface data
use utilities                                  ! Grid utilities
use vecsuv_m                                   ! Map to cartesian coordinates
use workglob_m                                 ! Additional grid interpolation
use work2_m                                    ! Diagnostic arrays

implicit none

include 'newmpar.h'                            ! Grid parameters
include 'const_phys.h'                         ! Physical constants
include 'darcdf.h'                             ! Netcdf data
include 'parm.h'                               ! Model configuration
include 'parmdyn.h'                            ! Dynamics parmaters
include 'parmgeom.h'                           ! Coordinate data
include 'soilv.h'                              ! Soil parameters
include 'stime.h'                              ! File date data

real, parameter :: iotol = 1.E-5      ! tolarance for iotest grid matching
      
integer, intent(in) :: nested, kdate_r, ktime_r
integer idv, isoil, nud_test
integer levk, levkin, ier, igas, nemi
integer i, j, k, n, mm, iq, numneg
integer, dimension(fwsize) :: isoilm_a
integer, dimension(ifull), intent(out) :: isflag
integer, dimension(7+3*ms) :: ierc
integer, dimension(3), save :: iers
#ifdef usempi3
integer, dimension(3) :: shsize
integer xx4_win, yy4_win
real(kind=8), dimension(:,:), pointer :: xx4, yy4
#else
real(kind=8), dimension(:,:), allocatable, save :: xx4, yy4
#endif
real(kind=8), dimension(dk*dk*6):: z_a, x_a, y_a
real, dimension(ifull,ms), intent(out) :: wb, wbice, tgg
real, dimension(ifull,3), intent(out) :: tggsn, smass, ssdn
real, dimension(:,:), intent(out) :: t, u, v, qg
real, dimension(ifull), intent(out) :: psl, zss, tss, fracice
real, dimension(ifull), intent(out) :: snowd, sicedep, ssdnn, snage
real, dimension(ifull) :: dum6, tss_l, tss_s, pmsl
real, dimension(dk*dk*6) :: wts_a  ! not used here or defined in call setxyz
real, dimension(fwsize) :: ucc
real, dimension(fwsize) :: fracice_a, sicedep_a
real, dimension(fwsize) :: tss_l_a, tss_s_a, tss_a
real, dimension(fwsize) :: t_a_lev, psl_a
real, dimension(:), allocatable, save :: zss_a, ocndep_l
real, dimension(kk+3) :: dumr
character(len=8) vname
character(len=3) trnum
logical tsstest, tst
logical, dimension(:), allocatable, save :: land_a, sea_a


! land-sea mask method (nemi=3 use soilt, nemi=2 use tgg, nemi=1 use zs)
nemi = 3
      
!retopo fields
nud_test = 1

      
! Determine if interpolation is required
iotest = 6*ik*ik==ifull_g .and. abs(rlong0x-rlong0)<iotol .and. abs(rlat0x-rlat0)<iotol .and. &
         abs(schmidtx-schmidt)<iotol .and. nsib==nsibx
if ( iotest ) then
  io_in = 1   ! no interpolation
else
  io_in = -1  ! interpolation
end if
if ( myid==0 ) write(6,*) "Interpolation iotest,io_in =",iotest,io_in

!--------------------------------------------------------------------
! Allocate interpolation, vertical level and mask arrays
! dk is only non-zero on myid==0
if ( .not.allocated(nface4) ) then
  allocate( nface4(ifull,4), xg4(ifull,4), yg4(ifull,4) )
end if
if ( newfile ) then
  !~ if ( allocated(sigin) ) then
    !~ deallocate( sigin, land_a, sea_a )
    !~ deallocate( axs_a, ays_a, azs_a )
    !~ deallocate( bxs_a, bys_a, bzs_a )          
  !~ end if
  !~ allocate( sigin(kk), land_a(fwsize), sea_a(fwsize) )
  allocate( axs_a(dk*dk*6), ays_a(dk*dk*6), azs_a(dk*dk*6) )
  allocate( bxs_a(dk*dk*6), bys_a(dk*dk*6), bzs_a(dk*dk*6) )
end if
      
!--------------------------------------------------------------------
! Determine input grid coordinates and interpolation arrays
if ( newfile .and. .not.iotest ) then
#ifdef usempi3
  shsize(1) = 1 + 4*ik
  shsize(2) = 1 + 4*ik
  call ccmpi_allocshdatar8(xx4,shsize(1:2),xx4_win)
  call ccmpi_allocshdatar8(yy4,shsize(1:2),yy4_win)
#else
  allocate( xx4(1+4*ik,1+4*ik), yy4(1+4*ik,1+4*ik) )
#endif

  if ( m_fly==1 ) then
    rlong4_l(:,1) = rlongg(:)*180./pi
    rlat4_l(:,1)  = rlatt(:)*180./pi
  end if
          
  if ( myid==0 ) then
    write(6,*) "Defining input file grid"
!   following setxyz call is for source data geom    ****   
    do iq = 1,dk*dk*6
      axs_a(iq) = iq
      ays_a(iq) = iq
      azs_a(iq) = iq
    end do 
    call setxyz(ik,rlong0x,rlat0x,-schmidtx,x_a,y_a,z_a,wts_a,axs_a,ays_a,azs_a,bxs_a,bys_a,bzs_a,xx4,yy4)
  end if ! (myid==0)
#ifdef usempi3
  call ccmpi_shepoch(xx4_win) ! also yy4_win
  if ( node_myid==0 ) then
    call ccmpi_bcastr8(xx4,0,comm_nodecaptian)
    call ccmpi_bcastr8(yy4,0,comm_nodecaptian)
  end if
  call ccmpi_shepoch(xx4_win) ! also yy4_win
#else
  call ccmpi_bcastr8(xx4,0,comm_world)
  call ccmpi_bcastr8(yy4,0,comm_world)
#endif
  
  ! calculate the rotated coords for host and model grid
  rotpoles = calc_rotpole(rlong0x,rlat0x)
  rotpole  = calc_rotpole(rlong0,rlat0)
  if ( myid==0 ) then
    write(6,*)'m_fly,nord ',m_fly,nord
    write(6,*)'kdate_r,ktime_r,ktau,ds',kdate_r,ktime_r,ktau,ds
    write(6,*)'rotpoles:'
    do i = 1,3
      write(6,'(3x,2i1,5x,2i1,5x,2i1,5x,3f8.4)') (i,j,j=1,3),(rotpoles(i,j),j=1,3)
    enddo
    if ( nmaxpr==1 ) then
      write(6,*)'in onthefly rotpole:'
      do i = 1,3
        write(6,'(3x,2i1,5x,2i1,5x,2i1,5x,3f8.4)') (i,j,j=1,3),(rotpole(i,j),j=1,3)
      enddo
      write(6,*)'xx4,yy4 ',xx4(id,jd),yy4(id,jd)
      write(6,*)'before latltoij for id,jd: ',id,jd
      write(6,*)'rlong0x,rlat0x,schmidtx ',rlong0x,rlat0x,schmidtx
    end if                ! (nmaxpr==1)
  end if                  ! (myid==0)

  ! setup interpolation arrays
  do mm = 1,m_fly  !  was 4, now may be set to 1 in namelist
    do iq = 1,ifull
      call latltoij(rlong4_l(iq,mm),rlat4_l(iq,mm),       & !input
                    rlong0x,rlat0x,schmidtx,              & !input
                    xg4(iq,mm),yg4(iq,mm),nface4(iq,mm),  & !output (source)
                    xx4,yy4,ik)
    end do
  end do
  
#ifdef usempi3
  call ccmpi_freeshdata(xx4_win)
  call ccmpi_freeshdata(yy4_win)
#else
  deallocate( xx4, yy4 )  
#endif

  ! Identify panels to be processed
  if ( myid==0 ) then
    nfacereq(:) = .true. ! this is the host processor for bcast
  else
    nfacereq(:) = .false.
    do n = 0,npanels
      nfacereq(n) = any( nface4(:,:)==n )
    end do
  end if
  
  ! Define filemap for MPI RMA method
  call file_wininit
  
  ! Define comm_face for MPI IBcast method
  call splitface
       
end if ! newfile .and. .not.iotest
      
! -------------------------------------------------------------------
! read time invariant data when file is first opened
! need global zss_a for (potentially) landsea mask and psl interpolation
! need global isoilm_a for (potentially) landsea mask
if ( newfile ) then

  ! read vertical levels and missing data checks
  !~ if ( myid==0 .or. pfall ) then
    !~ if ( myid==0 ) write(6,*) "Reading time invariant fields"
    !~ call ccnf_inq_varid(ncid,'lev',idv,tst)
    !~ if ( tst ) call ccnf_inq_varid(ncid,'layer',idv,tst)
    !~ if ( tst ) call ccnf_inq_varid(ncid,'sigma',idv,tst)
    !~ if ( tst) then
      !~ if ( myid==0 ) then
        !~ write(6,*) "No sigma data found in input file"
      !~ end if
      !~ if ( kk>1 ) then
        !~ if ( myid==0 ) then
          !~ write(6,*) "ERORR: multiple levels expected but no sigma data found ",kk
        !~ end if
        !~ call ccmpi_abort(-1)
      !~ end if
      !~ sigin(:) = 1.
    !~ else
      !~ call ccnf_get_vara(ncid,idv,1,kk,sigin)
      !~ if ( myid==0 ) then
        !~ write(6,'(" sigin=",(9f7.4))') (sigin(k),k=1,kk)
      !~ end if
    !~ end if
    ! check for missing data
    !~ iers(1:3) = 0
    !~ call ccnf_inq_varid(ncid,'mixr',idv,tst)
    !~ if ( tst ) iers(1) = -1
    !~ call ccnf_inq_varid(ncid,'siced',idv,tst)
    !~ if ( tst ) iers(2) = -1
    !~ call ccnf_inq_varid(ncid,'fracice',idv,tst)
    !~ if ( tst ) iers(3) = -1
  !~ end if
  
  ! bcast data to all processors unless all processes are reading input files
  if ( .not.pfall ) then
    !~ dumr(1:kk)      = sigin(1:kk)
    dumr(kk+1:kk+3) = real(iers(1:3))
    call ccmpi_bcast(dumr(1:kk+3),0,comm_world)
    !~ sigin(1:kk) = dumr(1:kk)
    !~ iers(1:3)   = nint(dumr(kk+1:kk+3))
  end if

  ! determine whether surface temperature needs to be interpolated (tsstest=.false.)
  tsstest = (iers(2)==0) .and. (iers(3)==0) .and. iotest
  if ( myid==0 ) write(6,*) "tsstest, iers ",tsstest,iers(1:3)
  if ( allocated(zss_a) ) deallocate( zss_a )
  if ( tsstest ) then
    ! load local surface temperature
    allocate( zss_a(ifull) )
    call histrd1(iarchi,ier,'zht',ik,zss_a,ifull)
  else if ( fnresid==1 ) then
    ! load global surface temperature using gather
    allocate( zss_a(6*dk*dk) )
    call histrd1(iarchi,ier,'zht',  ik,zss_a,6*ik*ik,nogather=.false.)
    call histrd1(iarchi,ier,'soilt',ik,ucc  ,6*ik*ik,nogather=.false.)
    if ( myid==0 ) then
      isoilm_a(:) = nint(ucc(:))
      if ( all(isoilm_a(:)==0) ) isoilm_a(:) = -1 ! missing value flag
    end if
  else
    ! load global surface temperature using RMA
    allocate( zss_a(fwsize) )
    call histrd1(iarchi,ier,'zht',  ik,zss_a,6*ik*ik,nogather=.true.)
    call histrd1(iarchi,ier,'soilt',ik,ucc  ,6*ik*ik,nogather=.true.)
    if ( fwsize>0 ) then
      isoilm_a(:) = nint(ucc(:))
      if ( all(isoilm_a(:)==0) ) isoilm_a(:) = -1 ! missing value flag
    end if
  end if
  
else
  ! use saved metadata  
  tsstest = (iers(2)==0) .and. (iers(3)==0) .and. iotest
endif ! newfile ..else..

! -------------------------------------------------------------------
! detemine the reference level below sig=0.9 (used to calculate psl)
!~ levk = 0
!~ levkin = 0
!~ if ( nested==0 .or. ( nested==1.and.nud_test/=0 ) ) then
  !~ do while( sig(levk+1)>0.9 ) ! nested grid
    !~ levk = levk + 1
  !~ end do
  !~ do while( sigin(levkin+1)>0.9 ) ! host grid
    !~ levkin = levkin + 1
  !~ end do
  !~ if ( levkin==0 ) then
    !~ write(6,*) "ERROR: Invalid sigma levels in input file"
    !~ write(6,*) "sigin = ",sigin
    !~ call ccmpi_abort(-1)
  !~ end if
!~ end if
! -------------------------------------------------------------------
! read atmospheric fields for nested=0 or nested=1.and.nud/=0

!~ ! winds
!~ ! read for nested=0 or nested=1
!~ u(1:ifull,:) = 0.
!~ v(1:ifull,:) = 0.
!~ ! mixing ratio
!~ qg(1:ifull,:) = qgmin
!**************************************************************
! This is the end of reading the nudging arrays
!**************************************************************

!--------------------------------------------------------------
! The following data is only read for initial conditions
if ( nested/=1 ) then


  !------------------------------------------------------------------
  ! check soil variables
  if ( myid==0 .or. pfall ) then
    ierc(7:7+3*ms) = 0
    if ( ccycle==0 ) then
      ierc(7) = -1
    else
      call ccnf_inq_varid(ncid,'glai',idv,tst)
      if ( tst ) ierc(7) = -1
    end if
    do k = 1,ms
      write(vname,'("tgg",I1.1)') k
      call ccnf_inq_varid(ncid,vname,idv,tst)
      if ( tst ) ierc(7+k) = -1
      write(vname,'("wetfrac",I1.1)') k
      call ccnf_inq_varid(ncid,vname,idv,tst)
      if ( tst ) ierc(7+ms+k) = -1
      write(vname,'("wb",I1.1)') k
      call ccnf_inq_varid(ncid,vname,idv,tst)
      if ( tst ) ierc(7+2*ms+k) = -1
    end do
  end if

  !--------------------------------------------------
  ! Read MLO sea-ice data
    call gethist1('swater',watbdy)

  !------------------------------------------------------------------
  ! Read soil moisture
  wb(:,:) = 20.5
  if ( all(ierc(8+ms:7+2*ms)==0) ) then
    call fillhist4('wetfrac',wb,ms,sea_a)
    wb(:,:) = wb(:,:) + 20. ! flag for fraction of field capacity
  else
    do k = 1,ms
      if ( ierc(7+ms+k)==0 ) then
        write(vname,'("wetfrac",I1.1)') k
      else if ( ierc(7+2*ms+k)==0 ) then
        write(vname,'("wb",I1.1)') k
      else if ( k<2 .and. ierc(7+2*ms+2)==0 ) then
        vname = "wb2"
      else if ( k<2 ) then
        vname = "wfg"
      else if ( ierc(7+2*ms+6)==0 ) then
        vname = "wb6"
      else
        vname = "wfb"
      end if
      if ( iotest ) then
        call histrd1(iarchi,ier,vname,ik,wb(:,k),ifull)
        if ( ierc(7+ms+k)==0 ) then
          wb(:,k) = wb(:,k) + 20. ! flag for fraction of field capacity
        end if
      else if ( fnresid==1 ) then
        call histrd1(iarchi,ier,vname,ik,ucc,6*ik*ik,nogather=.false.)
        if ( ierc(7+ms+k)==0 ) then
          ucc(:) = ucc(:) + 20.   ! flag for fraction of field capacity
        end if
        call fill_cc1_gather(ucc,sea_a)
        call doints1_gather(ucc,wb(:,k))
      else
        call histrd1(iarchi,ier,vname,ik,ucc,6*ik*ik,nogather=.true.)
        if ( ierc(7+ms+k)==0 ) then
          ucc(:) = ucc(:) + 20.   ! flag for fraction of field capacity
        end if
        call fill_cc1_nogather(ucc,sea_a)
        call doints1_nogather(ucc,wb(:,k))
      end if ! iotest
    end do
  end if
  !unpack field capacity into volumetric soil moisture
  if ( any(wb(:,:)>10.) ) then
    if ( mydiag ) write(6,*) "Unpacking wetfrac to wb",wb(idjd,1)
    wb(:,:) = wb(:,:) - 20.
    do iq = 1,ifull
      isoil = isoilm(iq)
      wb(iq,:) = (1.-wb(iq,:))*swilt(isoil) + wb(iq,:)*sfc(isoil)
    end do
    if ( mydiag ) write(6,*) "giving wb",wb(idjd,1)
  end if
  call fillhist1('wetfac',wetfac,sea_a)
  where ( .not.land )
    wetfac(:) = 1.
  end where

  ! -----------------------------------------------------------------
  ! soil ice and snow data
  call gethist4('wbice',wbice,ms) ! SOIL ICE
  
    call gethist4('tggsn',tggsn,3)
    if ( all(tggsn==0.) ) tggsn=280.

  call gethist4('smass',smass,3)
  call gethist4('ssdn',ssdn,3)
  do k=1,3
    if ( all(ssdn(:,k)==0.) ) then
      where ( snowd>100. )
        ssdn(:,k)=240.
      elsewhere
        ssdn(:,k)=140.
      end where
    end if
  end do
  ssdnn=ssdn(:,1)
  call gethist1('snage',snage)
  call gethist1('sflag',dum6)
  isflag=nint(dum6)

        
endif    ! (nested/=1)

!**************************************************************
! This is the end of reading the initial arrays
!**************************************************************         

! -------------------------------------------------------------------
! tgg holds file surface temperature when no MLO

  where ( .not.land )
    tgg(:,1) = tss
  end where


! -------------------------------------------------------------------
! set-up for next read of file
iarchi = iarchi + 1
kdate_s = kdate_r
ktime_s = ktime_r + 1


return
end subroutine onthefly_work


! *****************************************************************************
! INTERPOLATION ROUTINES                         

! Main interface
! Note that sx is a global array for all processors

subroutine doints1_nogather(s,sout)
      
use cc_mpi                 ! CC MPI routines
use infile                 ! Input file routines

implicit none
     
include 'newmpar.h'        ! Grid parameters
include 'parm.h'           ! Model configuration
      
integer mm
real, dimension(:), intent(in) :: s
real, dimension(:), intent(inout) :: sout
real, dimension(ifull,m_fly) :: wrk
real, dimension(pil*pjl*pnpan,size(filemap),fncount) :: abuf
real, dimension(ik+4,ik+4,0:npanels) :: sx

call START_LOG(otf_ints1_begin)

if ( .not.allocated(filemap) ) then
  write(6,*) "ERROR: Mapping for RMA file windows has not been defined"
  call ccmpi_abort(-1)
end if

! This version uses MPI RMA to distribute data
call ccmpi_filewinget(abuf,s)

sx(1:ik+4,1:ik+4,0:npanels) = 0.
call ccmpi_filewinunpack(sx,abuf)
call sxpanelbounds(sx)

if ( nord==1 ) then   ! bilinear
  do mm = 1,m_fly     !  was 4, now may be 1
    call ints_blb(sx,wrk(:,mm),nface4(:,mm),xg4(:,mm),yg4(:,mm))
  end do
else                  ! bicubic
  do mm = 1,m_fly     !  was 4, now may be 1
    call intsb(sx,wrk(:,mm),nface4(:,mm),xg4(:,mm),yg4(:,mm))
  end do
end if   ! (nord==1)  .. else ..
sout(1:ifull) = sum(wrk(:,:), dim=2)/real(m_fly)

call END_LOG(otf_ints1_end)

return
end subroutine doints1_nogather

subroutine doints1_gather(s,sout)
      
use cc_mpi                 ! CC MPI routines
use infile                 ! Input file routines

implicit none
     
include 'newmpar.h'        ! Grid parameters
include 'parm.h'           ! Model configuration
      
integer mm, n, ik2
real, dimension(:), intent(in) :: s
real, dimension(:), intent(inout) :: sout
real, dimension(ifull,m_fly) :: wrk
real, dimension(ik+4,ik+4,0:npanels) :: sx

call START_LOG(otf_ints1_begin)

if ( .not.bcst_allocated ) then
  write(6,*) "ERROR: Bcst commuicators have not been defined"
  call ccmpi_abort(-1)
end if

! This version uses MPI_Bcast to distribute data
sx(1:ik+4,1:ik+4,0:npanels) = 0.
if ( dk>0 ) then
  ik2 = ik*ik
  sx(3:ik+2,3:ik+2,0:npanels) = reshape( s(1:(npanels+1)*ik2), (/ ik, ik, npanels+1 /) )
  call sxpanelbounds(sx)
end if
do n = 0,npanels
  ! send each face of the host dataset to processes that require it
  if ( nfacereq(n) ) then
    call ccmpi_bcast(sx(:,:,n),0,comm_face(n))
  end if
end do  ! n loop

if ( nord==1 ) then   ! bilinear
  do mm = 1,m_fly     !  was 4, now may be 1
    call ints_blb(sx,wrk(:,mm),nface4(:,mm),xg4(:,mm),yg4(:,mm))
  end do
else                  ! bicubic
  do mm = 1,m_fly     !  was 4, now may be 1
    call intsb(sx,wrk(:,mm),nface4(:,mm),xg4(:,mm),yg4(:,mm))
  end do
end if   ! (nord==1)  .. else ..
sout(1:ifull) = sum(wrk(:,:), dim=2)/real(m_fly)

call END_LOG(otf_ints1_end)

return
end subroutine doints1_gather

subroutine doints4_nogather(s,sout)
      
use cc_mpi                 ! CC MPI routines
use infile                 ! Input file routines

implicit none
     
include 'newmpar.h'        ! Grid parameters
include 'parm.h'           ! Model configuration
      
integer mm, k, kx, kb, ke, kn
real, dimension(:,:), intent(in) :: s
real, dimension(:,:), intent(inout) :: sout
real, dimension(ifull,m_fly) :: wrk
real, dimension(pil*pjl*pnpan,size(filemap),fncount,kblock) :: abuf
real, dimension(ik+4,ik+4,0:npanels) :: sx

call START_LOG(otf_ints4_begin)

kx = size(sout, 2)

if ( .not.allocated(filemap) ) then
  write(6,*) "ERROR: Mapping for RMA file windows has not been defined"
  call ccmpi_abort(-1)
end if

do kb = 1,kx,kblock
  ke = min(kb+kblock-1, kx)
  kn = ke - kb + 1

  ! This version uses MPI RMA to distribute data
  call ccmpi_filewinget(abuf(:,:,:,1:kn),s(:,kb:ke))
    
  ! MJT notes - sx can be made into a shared memory array,
  ! although this requires a MPI_Fence when the abuf
  ! arrays are unpacked for each level.
  
  if ( nord==1 ) then   ! bilinear
    do k = 1,kn
      sx(1:ik+4,1:ik+4,0:npanels) = 0.
      call ccmpi_filewinunpack(sx,abuf(:,:,:,k))
      call sxpanelbounds(sx)
      do mm = 1,m_fly     !  was 4, now may be 1
        call ints_blb(sx,wrk(:,mm),nface4(:,mm),xg4(:,mm),yg4(:,mm))
      end do
      sout(1:ifull,k+kb-1) = sum(wrk(:,:), dim=2)/real(m_fly)
    end do
  else                  ! bicubic
    do k = 1,kn
      sx(1:ik+4,1:ik+4,0:npanels) = 0.
      call ccmpi_filewinunpack(sx,abuf(:,:,:,k))
      call sxpanelbounds(sx)
      do mm = 1,m_fly     !  was 4, now may be 1
        call intsb(sx,wrk(:,mm),nface4(:,mm),xg4(:,mm),yg4(:,mm))
      end do
      sout(1:ifull,k+kb-1) = sum(wrk(:,:), dim=2)/real(m_fly)
    end do
  end if   ! (nord==1)  .. else ..

end do
  
call END_LOG(otf_ints4_end)

return
end subroutine doints4_nogather

subroutine doints4_gather(s,sout)
      
use cc_mpi                 ! CC MPI routines
use infile                 ! Input file routines

implicit none
     
include 'newmpar.h'        ! Grid parameters
include 'parm.h'           ! Model configuration
      
integer mm, n, k, kx, ik2
real, dimension(:,:), intent(in) :: s
real, dimension(:,:), intent(inout) :: sout
real, dimension(ifull,m_fly) :: wrk
real, dimension(ik+4,ik+4,0:npanels) :: sx

call START_LOG(otf_ints4_begin)

kx = size(sout,2)

if ( .not.bcst_allocated ) then
  write(6,*) "ERROR: Bcst commuicators have not been defined"
  call ccmpi_abort(-1)
end if

do k = 1,kx

  ! This version uses MPI_Bcast to distribute data
  sx(1:ik+4,1:ik+4,0:npanels) = 0.
  if ( dk>0 ) then
    ik2 = ik*ik
    !     first extend s arrays into sx - this one -1:il+2 & -1:il+2
    sx(3:ik+2,3:ik+2,0:npanels) = reshape( s(1:(npanels+1)*ik2,k), (/ ik, ik, npanels+1 /) )
    call sxpanelbounds(sx)
  end if
  do n = 0,npanels
    if ( nfacereq(n) ) then
      call ccmpi_bcast(sx(:,:,n),0,comm_face(n))
     end if
  end do

  if ( nord==1 ) then   ! bilinear
    do mm = 1,m_fly     !  was 4, now may be 1
      call ints_blb(sx,wrk(:,mm),nface4(:,mm),xg4(:,mm),yg4(:,mm))
    end do
  else                  ! bicubic
    do mm = 1,m_fly     !  was 4, now may be 1
      call intsb(sx,wrk(:,mm),nface4(:,mm),xg4(:,mm),yg4(:,mm))
    end do
  end if   ! (nord==1)  .. else ..
  sout(1:ifull,k) = sum( wrk(:,:), dim=2 )/real(m_fly)
  
end do
  
call END_LOG(otf_ints4_end)

return
end subroutine doints4_gather

subroutine sxpanelbounds(sx_l)

implicit none

include 'newmpar.h'

integer i, n, n_w, n_e, n_n, n_s
real, dimension(-1:ik+2,-1:ik+2,0:npanels), intent(inout) :: sx_l

do n = 0,npanels
  if ( nfacereq(n) ) then
    if ( mod(n,2)==0 ) then
      n_w = mod(n+5, 6)
      n_e = mod(n+2, 6)
      n_n = mod(n+1, 6)
      n_s = mod(n+4, 6)
      do i = 1,ik
        sx_l(0,i,n)    = sx_l(ik,i,n_w)
        sx_l(-1,i,n)   = sx_l(ik-1,i,n_w)
        sx_l(ik+1,i,n) = sx_l(ik+1-i,1,n_e)
        sx_l(ik+2,i,n) = sx_l(ik+1-i,2,n_e)
        sx_l(i,ik+1,n) = sx_l(i,1,n_n)
        sx_l(i,ik+2,n) = sx_l(i,2,n_n)
        sx_l(i,0,n)    = sx_l(ik,ik+1-i,n_s)
        sx_l(i,-1,n)   = sx_l(ik-1,ik+1-i,n_s)
      end do ! i
      sx_l(-1,0,n)      = sx_l(ik,2,n_w)        ! wws
      sx_l(0,-1,n)      = sx_l(ik,ik-1,n_s)     ! wss
      sx_l(0,0,n)       = sx_l(ik,1,n_w)        ! ws
      sx_l(ik+1,0,n)    = sx_l(ik,1,n_e)        ! es  
      sx_l(ik+2,0,n)    = sx_l(ik-1,1,n_e)      ! ees 
      sx_l(-1,ik+1,n)   = sx_l(ik,ik-1,n_w)     ! wwn
      sx_l(0,ik+2,n)    = sx_l(ik-1,ik,n_w)     ! wnn
      sx_l(ik+2,ik+1,n) = sx_l(2,1,n_e)         ! een  
      sx_l(ik+1,ik+2,n) = sx_l(1,2,n_e)         ! enn  
      sx_l(0,ik+1,n)    = sx_l(ik,ik,n_w)       ! wn  
      sx_l(ik+1,ik+1,n) = sx_l(1,1,n_e)         ! en  
      sx_l(ik+1,-1,n)   = sx_l(ik,2,n_e)        ! ess        
    else
      n_w = mod(n+4, 6)
      n_e = mod(n+1, 6)
      n_n = mod(n+2, 6)
      n_s = mod(n+5, 6)
      do i = 1,ik
        sx_l(0,i,n)    = sx_l(ik+1-i,ik,n_w)
        sx_l(-1,i,n)   = sx_l(ik+1-i,ik-1,n_w)
        sx_l(ik+1,i,n) = sx_l(1,i,n_e)
        sx_l(ik+2,i,n) = sx_l(2,i,n_e)
        sx_l(i,ik+1,n) = sx_l(1,ik+1-i,n_n)
        sx_l(i,ik+2,n) = sx_l(2,ik+1-i,n_n)
        sx_l(i,0,n)    = sx_l(i,ik,n_s)
        sx_l(i,-1,n)   = sx_l(i,ik-1,n_s)
      end do ! i
      sx_l(-1,0,n)      = sx_l(ik-1,ik,n_w)    ! wws
      sx_l(0,-1,n)      = sx_l(2,ik,n_s)       ! wss
      sx_l(0,0,n)       = sx_l(ik,ik,n_w)      ! ws
      sx_l(ik+1,0,n)    = sx_l(1,1,n_e)        ! es
      sx_l(ik+2,0,n)    = sx_l(1,2,n_e)        ! ees
      sx_l(-1,ik+1,n)   = sx_l(2,ik,n_w)       ! wwn   
      sx_l(0,ik+2,n)    = sx_l(1,ik-1,n_w)     ! wnn  
      sx_l(ik+2,ik+1,n) = sx_l(1,ik-1,n_e)     ! een  
      sx_l(ik+1,ik+2,n) = sx_l(2,ik,n_e)       ! enn  
      sx_l(0,ik+1,n)    = sx_l(1,ik,n_w)       ! wn  
      sx_l(ik+1,ik+1,n) = sx_l(1,ik,n_e)       ! en  
      sx_l(ik+1,-1,n)   = sx_l(2,1,n_e)        ! ess         
    end if   ! mod(n,2)==0 ..else..
  end if     ! nfacereq(n)
end do       ! n loop

return
end subroutine sxpanelbounds

subroutine intsb(sx_l,sout,nface_l,xg_l,yg_l)
      
!     same as subr ints, but with sout passed back and no B-S      
!     s is input; sout is output array
!     later may wish to save idel etc between array calls
!     this one does linear interp in x on outer y sides
!     doing x-interpolation before y-interpolation
!     This is a global routine 

implicit none
      
include 'newmpar.h'  ! Grid parameters
include 'parm.h'     ! Model configuration

integer, dimension(ifull), intent(in) :: nface_l
integer :: idel, jdel
integer :: n, iq
real, dimension(ifull), intent(inout) :: sout
real xxg, yyg
real cmin, cmax
real, intent(in), dimension(ifull) :: xg_l, yg_l
real, dimension(-1:ik+2,-1:ik+2,0:npanels), intent(in) :: sx_l
real, dimension(2:3) :: dmul
real, dimension(1:4) :: cmul, emul, rmul

do iq = 1,ifull   ! runs through list of target points
  n = nface_l(iq)
  idel = int(xg_l(iq))
  xxg = xg_l(iq) - idel
  jdel = int(yg_l(iq))
  yyg = yg_l(iq) - jdel

  ! bi-cubic
  cmul(1) = (1.-xxg)*(2.-xxg)*(-xxg)/6.
  cmul(2) = (1.-xxg)*(2.-xxg)*(1.+xxg)/2.
  cmul(3) = xxg*(1.+xxg)*(2.-xxg)/2.
  cmul(4) = (1.-xxg)*(-xxg)*(1.+xxg)/6.
  dmul(2) = (1.-xxg)
  dmul(3) = xxg
  emul(1) = (1.-yyg)*(2.-yyg)*(-yyg)/6.
  emul(2) = (1.-yyg)*(2.-yyg)*(1.+yyg)/2.
  emul(3) = yyg*(1.+yyg)*(2.-yyg)/2.
  emul(4) = (1.-yyg)*(-yyg)*(1.+yyg)/6.
  cmin = minval( sx_l(idel:idel+1,jdel:jdel+1,n) )
  cmax = maxval( sx_l(idel:idel+1,jdel:jdel+1,n) )  
  rmul(1) = sum( sx_l(idel:idel+1,  jdel-1,n)*dmul(2:3) )
  rmul(2) = sum( sx_l(idel-1:idel+2,jdel,  n)*cmul(1:4) )
  rmul(3) = sum( sx_l(idel-1:idel+2,jdel+1,n)*cmul(1:4) )
  rmul(4) = sum( sx_l(idel:idel+1,  jdel+2,n)*dmul(2:3) )
  
  sout(iq) = min( max( cmin, sum( rmul(1:4)*emul(1:4) ) ), cmax ) ! Bermejo & Staniforth
end do    ! iq loop

return
end subroutine intsb

subroutine ints_blb(sx_l,sout,nface_l,xg_l,yg_l) 
      
!     this one does bi-linear interpolation only

implicit none
      
include 'newmpar.h'  ! Grid parameters
include 'parm.h'     ! Model configuration

integer :: n, iq, idel, jdel
integer, intent(in), dimension(ifull) :: nface_l
real, dimension(ifull), intent(inout) :: sout
real, intent(in), dimension(ifull) :: xg_l, yg_l
real, dimension(-1:ik+2,-1:ik+2,0:npanels), intent(in) :: sx_l
real :: xxg, yyg

do iq = 1,ifull  ! runs through list of target points
  n = nface_l(iq)
  idel = int(xg_l(iq))
  xxg = xg_l(iq) - idel
  jdel = int(yg_l(iq))
  yyg = yg_l(iq) - jdel
  sout(iq) = yyg*(xxg*sx_l(idel+1,jdel+1,n) + (1.-xxg)*sx_l(idel,jdel+1,n)) + &
          (1.-yyg)*(xxg*sx_l(idel+1,jdel,n) + (1.-xxg)*sx_l(idel,jdel,n))
enddo    ! iq loop

return
end subroutine ints_blb

! *****************************************************************************
! FILL ROUTINES

subroutine fill_cc1_nogather(a_io,land_a)
      
! routine fills in interior of an array which has undefined points
! this version is for multiple input files

use cc_mpi          ! CC MPI routines
use infile          ! Input file routines

implicit none

integer nrem, j, n
integer ncount, cc, ipf
integer, dimension(pil) :: neighc
real, parameter :: value=999.       ! missing value flag
real, dimension(:), intent(inout) :: a_io
real, dimension(0:pil+1,0:pjl+1,pnpan,mynproc) :: c_io
real, dimension(pil,4) :: c
logical, dimension(:), intent(in) :: land_a
logical, dimension(pil,4) :: maskc

! only perform fill on processors reading input files
if ( fwsize==0 ) return
  
where ( land_a(1:fwsize) )
  a_io(1:fwsize) = value
end where
ncount = count( abs(a_io(1:fwsize)-value)<1.E-6 )
call ccmpi_allreduce(ncount,nrem,'sum',comm_ip)
if ( nrem==6*ik*ik ) return
 
do while ( nrem>0 )
  c_io(1:pil,1:pjl,1:pnpan,1:mynproc) = reshape( a_io(1:fwsize), (/ pil, pjl, pnpan, mynproc /) )
  call ccmpi_filebounds(c_io,comm_ip)
  do ipf = 1,mynproc
    do n = 1,pnpan
      do j = 1,pjl
        c(1:pil,1) = c_io(1:pil,j+1,n,ipf)
        c(1:pil,2) = c_io(1:pil,j-1,n,ipf)
        c(1:pil,3) = c_io(2:pil+1,j,n,ipf)
        c(1:pil,4) = c_io(0:pil-1,j,n,ipf)
        maskc(1:pil,1:4) = c(1:pil,1:4)/=value
        neighc(1:pil) = count( maskc(1:pil,1:4), dim=2 )
        cc = (j-1)*pil + (n-1)*pil*pjl + (ipf-1)*pil*pjl*pnpan
        where ( neighc(1:pil)>0 .and. c_io(1:pil,j,n,ipf)==value )
          a_io(1+cc:pil+cc) = sum( c(1:pil,1:4), mask=maskc(1:pil,1:4), dim=2 )/real(neighc(1:pil))
        end where
      end do
    end do
  end do
  ! test for convergence
  ncount = count( abs(a_io(1:fwsize)-value)<1.E-6 )
  call ccmpi_allreduce(ncount,nrem,'sum',comm_ip)
end do
      
return
end subroutine fill_cc1_nogather

subroutine fill_cc1_gather(a_io,land_a)
      
! routine fills in interior of an array which has undefined points
! this version is for a single input file

use cc_mpi          ! CC MPI routines
use infile          ! Input file routines

implicit none

integer nrem, i, iq, j, n
integer iminb, imaxb, jminb, jmaxb
integer is, ie, js, je
integer, dimension(0:5) :: imin, imax, jmin, jmax
integer, dimension(dk) :: neighb
integer, parameter, dimension(0:5) :: npann=(/1,103,3,105,5,101/)
integer, parameter, dimension(0:5) :: npane=(/102,2,104,4,100,0/)
integer, parameter, dimension(0:5) :: npanw=(/5,105,1,101,3,103/)
integer, parameter, dimension(0:5) :: npans=(/104,0,100,2,102,4/)
real, parameter :: value=999.       ! missing value flag
real, dimension(:), intent(inout) :: a_io
real, dimension(6*dk*dk) :: b_io
real, dimension(0:dk+1) :: a
real, dimension(dk) :: b_north, b_south, b_east, b_west
real, dimension(dk,4) :: b
logical, dimension(:), intent(in) :: land_a
logical, dimension(dk,4) :: mask
logical lflag

! only perform fill on myid==0
if ( dk==0 ) return

where ( land_a(1:6*dk*dk) )
  a_io(1:6*dk*dk) = value
end where
if ( all(abs(a_io(1:6*dk*dk)-value)<1.E-6) ) return

imin(0:5) = 1
imax(0:5) = dk
jmin(0:5) = 1
jmax(0:5) = dk
          
nrem = 1    ! Just for first iteration
do while ( nrem>0 )
  nrem = 0
  b_io(1:6*dk*dk) = a_io(1:6*dk*dk)
  ! MJT restricted fill
  do n = 0,5
    
    iminb = dk
    imaxb = 1
    jminb = dk
    jmaxb = 1
    
    ! north
    if (npann(n)<100) then
      do i = 1,dk
        iq=i+npann(n)*dk*dk
        b_north(i) = b_io(iq)
      end do
    else
      do i = 1,dk
        iq=1+(dk-i)*dk+(npann(n)-100)*dk*dk
        b_north(i) = b_io(iq)
      end do
    end if
    ! south
    if (npans(n)<100) then
      do i = 1,dk
        iq=i+(dk-1)*dk+npans(n)*dk*dk
        b_south(i) = b_io(iq)
      end do
    else
      do i = 1,dk
        iq=dk+(dk-i)*dk+(npans(n)-100)*dk*dk
        b_south(i) = b_io(iq)
      end do
    end if
    ! east
    if (npane(n)<100) then
      do j = 1,dk
        iq=1+(j-1)*dk+npane(n)*dk*dk
        b_east(j) = b_io(iq)
      end do
    else
      do j = 1,dk
        iq=dk+1-j+(npane(n)-100)*dk*dk
        b_east(j) = b_io(iq)
      end do
    end if
    ! west
    if (npanw(n)<100) then
      do j = 1,dk
        iq=dk+(j-1)*dk+npanw(n)*dk*dk
        b_west(j) = b_io(iq)
      end do
    else
      do j = 1,dk
        iq=dk+1-j+(dk-1)*dk+(npanw(n)-100)*dk*dk
        b_west(j) = b_io(iq)
      end do
    end if

    is = imin(n)
    ie = imax(n)
    js = jmin(n)
    je = jmax(n)
    
    if ( js==1 ) then
      ! j = 1
      a(0)     = b_west(1)
      a(dk+1)  = b_east(1)
      a(max(is-1,1))  = b_io(max(is-1,1)+n*dk*dk)
      a(min(ie+1,dk)) = b_io(min(ie+1,dk)+n*dk*dk)
      a(is:ie) = b_io(is+n*dk*dk:ie+n*dk*dk)
      b(is:ie,1) = b_io(is+dk+n*dk*dk:ie+dk+n*dk*dk) ! north
      b(is:ie,2) = b_south(is:ie)                    ! south
      b(is:ie,3) = a(is+1:ie+1)                      ! east
      b(is:ie,4) = a(is-1:ie-1)                      ! west
      mask(is:ie,1:4) = b(is:ie,1:4)/=value
      neighb(is:ie) = count( mask(is:ie,1:4), dim=2)
      where ( neighb(is:ie)>0 .and. a(is:ie)==value )
        a_io(is+n*dk*dk:ie+n*dk*dk) = sum( b(is:ie,1:4), mask=mask(is:ie,1:4), dim=2)/real(neighb(is:ie))
      end where
      lflag = .false.
      do i = is,ie
        if ( neighb(i)==0 ) then
          nrem = nrem + 1 ! current number of points without a neighbour
          iminb = min(i, iminb)
          imaxb = max(i, imaxb)
          lflag = .true.
        end if
      end do
      if ( lflag ) then
        jminb = min(1, jminb)
        jmaxb = max(1, jmaxb)
      end if
    end if
    do j = max(js,2),min(je,dk-1)
      a(0)     = b_west(j)
      a(dk+1)  = b_east(j)
      a(max(is-1,1))  = b_io(max(is-1,1)+(j-1)*dk+n*dk*dk)
      a(min(ie+1,dk)) = b_io(min(ie+1,dk)+(j-1)*dk+n*dk*dk)
      a(is:ie) = b_io(is+(j-1)*dk+n*dk*dk:ie+(j-1)*dk+n*dk*dk)
      b(is:ie,1) = b_io(is+j*dk+n*dk*dk:ie+j*dk+n*dk*dk)         ! north
      b(is:ie,2) = b_io(is+(j-2)*dk+n*dk*dk:ie+(j-2)*dk+n*dk*dk) ! south
      b(is:ie,3) = a(is+1:ie+1)                                  ! east
      b(is:ie,4) = a(is-1:ie-1)                                  ! west
      mask(is:ie,1:4) = b(is:ie,1:4)/=value
      neighb(is:ie) = count( mask(is:ie,1:4), dim=2)
      where ( neighb(is:ie)>0 .and. a(is:ie)==value )
        a_io(is+(j-1)*dk+n*dk*dk:ie+(j-1)*dk+n*dk*dk) = sum( b(is:ie,1:4), mask=mask(is:ie,1:4), dim=2)/real(neighb(is:ie))
      end where
      lflag = .false.
      do i = is,ie
        if ( neighb(i)==0 ) then
          nrem = nrem + 1 ! current number of points without a neighbour
          iminb = min(i, iminb)
          imaxb = max(i, imaxb)
          lflag = .true.
        end if
      end do
      if ( lflag ) then
        jminb = min(j, jminb)
        jmaxb = max(j, jmaxb)
      end if
    end do
    if ( je==dk ) then
      ! j = dk
      a(0)     = b_west(dk)
      a(dk+1)  = b_east(dk)
      a(max(is-1,1))  = b_io(max(is-1,1)-dk+(n+1)*dk*dk)
      a(min(ie+1,dk)) = b_io(min(ie+1,dk)-dk+(n+1)*dk*dk)
      a(is:ie) = b_io(is-dk+(n+1)*dk*dk:ie-dk+(n+1)*dk*dk)
      b(is:ie,1) = b_north(is:ie)                                ! north
      b(is:ie,2) = b_io(is-2*dk+(n+1)*dk*dk:ie-2*dk+(n+1)*dk*dk) ! south
      b(is:ie,3) = a(is+1:ie+1)                                  ! east
      b(is:ie,4) = a(is-1:ie-1)                                  ! west
      mask(is:ie,1:4) = b(is:ie,1:4)/=value
      neighb(is:ie) = count( mask(is:ie,1:4), dim=2)
      where ( neighb(is:ie)>0 .and. a(is:ie)==value )
        a_io(is-dk+(n+1)*dk*dk:ie-dk+(n+1)*dk*dk) = sum( b(is:ie,1:4), mask=mask(is:ie,1:4), dim=2)/real(neighb(is:ie))
      end where
      lflag = .false.
      do i = is,ie
        if ( neighb(i)==0 ) then
          nrem = nrem + 1 ! current number of points without a neighbour
          iminb = min(i, iminb)
          imaxb = max(i, imaxb)
          lflag = .true.
        end if
      end do
      if ( lflag ) then
        jminb = min(dk, jminb)
        jmaxb = max(dk, jmaxb)
      end if
    end if
    
    imin(n) = iminb
    imax(n) = imaxb
    jmin(n) = jminb
    jmax(n) = jmaxb
  end do
end do
  
return
end subroutine fill_cc1_gather

subroutine fill_cc4_nogather(a_io,land_a)
      
! routine fills in interior of an array which has undefined points
! this version is distributed over processes with input files

use cc_mpi          ! CC MPI routines
use infile          ! Input file routines

implicit none

integer nrem, j, n, k, kx
integer ncount, cc, ipf
integer, dimension(pil) :: neighc
real, parameter :: value=999.       ! missing value flag
real, dimension(:,:), intent(inout) :: a_io
real, dimension(0:pil+1,0:pjl+1,pnpan,mynproc,size(a_io,2)) :: c_io
real, dimension(pil,4) :: c
logical, dimension(:), intent(in) :: land_a
logical, dimension(pil,4) :: maskc

kx = size(a_io,2)

! only perform fill on processors reading input files
if ( fwsize==0 ) return

do k = 1,kx
  where ( land_a(1:fwsize) )
    a_io(1:fwsize,k) = value
  end where
end do
ncount = count( abs(a_io(1:fwsize,kx)-value)<1.E-6 )
call ccmpi_allreduce(ncount,nrem,'sum',comm_ip)
if ( nrem==6*ik*ik ) return
 
do while ( nrem > 0 )
  c_io(1:pil,1:pjl,1:pnpan,1:mynproc,1:kx) = reshape( a_io(1:fwsize,1:kx), (/ pil, pjl, pnpan, mynproc, kx /) )
  call ccmpi_filebounds(c_io,comm_ip)
  do k = 1,kx
    do ipf = 1,mynproc
      do n = 1,pnpan
        do j = 1,pjl
          c(1:pil,1) = c_io(1:pil,j+1,n,ipf,k)
          c(1:pil,2) = c_io(1:pil,j-1,n,ipf,k)
          c(1:pil,3) = c_io(2:pil+1,j,n,ipf,k)
          c(1:pil,4) = c_io(0:pil-1,j,n,ipf,k)
          maskc(1:pil,1:4) = c(1:pil,1:4)/=value
          neighc(1:pil) = count( maskc(1:pil,1:4), dim=2)
          cc = (j-1)*pil + (n-1)*pil*pjl + (ipf-1)*pil*pjl*pnpan
          where ( neighc(1:pil)>0 .and. c_io(1:pil,j,n,ipf,k)==value )
            a_io(1+cc:pil+cc,k) = sum( c(1:pil,1:4), mask=maskc(1:pil,1:4), dim=2)/real(neighc(1:pil))
          end where
        end do
      end do
    end do
  end do
  ! test for convergence
  ncount = count( abs(a_io(1:fwsize,kx)-value)<1.E-6 )
  call ccmpi_allreduce(ncount,nrem,'sum',comm_ip)
end do
      
return
end subroutine fill_cc4_nogather

subroutine fill_cc4_gather(a_io,land_a)
      
! routine fills in interior of an array which has undefined points
! this version is for a single input file

use cc_mpi          ! CC MPI routines
use infile          ! Input file routines

implicit none

integer nrem, i, iq, j, n, k, kx
integer iminb, imaxb, jminb, jmaxb
integer is, ie, js, je
integer, dimension(0:5) :: imin, imax, jmin, jmax
integer, dimension(dk) :: neighb
integer, parameter, dimension(0:5) :: npann=(/1,103,3,105,5,101/)
integer, parameter, dimension(0:5) :: npane=(/102,2,104,4,100,0/)
integer, parameter, dimension(0:5) :: npanw=(/5,105,1,101,3,103/)
integer, parameter, dimension(0:5) :: npans=(/104,0,100,2,102,4/)
real, parameter :: value=999.       ! missing value flag
real, dimension(:,:), intent(inout) :: a_io
real, dimension(6*dk*dk,size(a_io,2)) :: b_io
real, dimension(0:dk+1) :: a
real, dimension(dk,size(a_io,2)) :: b_north, b_south, b_east, b_west
real, dimension(dk,4) :: b
logical, dimension(:), intent(in) :: land_a
logical, dimension(dk,4) :: mask
logical lflag

kx = size(a_io,2)

! only perform fill on myid==0
if ( dk==0 ) return

do k = 1,kx
  where ( land_a(1:6*dk*dk) )
    a_io(1:6*dk*dk,k)=value
  end where
end do
if ( all(abs(a_io(1:6*dk*dk,kx)-value)<1.E-6) ) return

imin(0:5) = 1
imax(0:5) = dk
jmin(0:5) = 1
jmax(0:5) = dk
          
nrem = 1    ! Just for first iteration
do while ( nrem>0 )
  nrem = 0
  b_io(1:6*dk*dk,1:kx) = a_io(1:6*dk*dk,1:kx)
  ! MJT restricted fill
  do n = 0,5
       
    iminb = dk
    imaxb = 1
    jminb = dk
    jmaxb = 1
    
    ! north
    if (npann(n)<100) then
      do k = 1,kx
        do i = 1,dk
          iq=i+npann(n)*dk*dk
          b_north(i,k) = b_io(iq,k)
        end do
      end do
    else
      do k = 1,kx
        do i = 1,dk
          iq=1+(dk-i)*dk+(npann(n)-100)*dk*dk
          b_north(i,k) = b_io(iq,k)
        end do
      end do
    end if
    ! south
    if (npans(n)<100) then
      do k = 1,kx
        do i = 1,dk
          iq=i+(dk-1)*dk+npans(n)*dk*dk
          b_south(i,k) = b_io(iq,k)
        end do
      end do
    else
      do k = 1,kx
        do i = 1,dk
          iq=dk+(dk-i)*dk+(npans(n)-100)*dk*dk
          b_south(i,k) = b_io(iq,k)
        end do
      end do
    end if
    ! east
    if (npane(n)<100) then
      do k = 1,kx
        do j = 1,dk
          iq=1+(j-1)*dk+npane(n)*dk*dk
          b_east(j,k) = b_io(iq,k)
        end do
      end do
    else
      do k = 1,kx
        do j = 1,dk
          iq=dk+1-j+(npane(n)-100)*dk*dk
          b_east(j,k) = b_io(iq,k)
        end do
      end do
    end if
    ! west
    if (npanw(n)<100) then
      do k = 1,kx
        do j = 1,dk
          iq=dk+(j-1)*dk+npanw(n)*dk*dk
          b_west(j,k) = b_io(iq,k)
        end do
      end do
    else
      do k = 1,kx
        do j = 1,dk
          iq=dk+1-j+(dk-1)*dk+(npanw(n)-100)*dk*dk
          b_west(j,k) = b_io(iq,k)
        end do
      end do
    end if

    is = imin(n)
    ie = imax(n)
    js = jmin(n)
    je = jmax(n)
    
    if ( js==1 ) then
      ! j = 1
      do k = 1,kx
        a(0)     = b_west(1,k)
        a(dk+1)  = b_east(1,k)
        a(max(is-1,1))  = b_io(max(is-1,1)+n*dk*dk,k)
        a(min(ie+1,dk)) = b_io(min(ie+1,dk)+n*dk*dk,k)
        a(is:ie) = b_io(is+n*dk*dk:ie+n*dk*dk,k)
        b(is:ie,1) = b_io(is+dk+n*dk*dk:ie+dk+n*dk*dk,k) ! north
        b(is:ie,2) = b_south(is:ie,k)                    ! south
        b(is:ie,3) = a(is+1:ie+1)                        ! east
        b(is:ie,4) = a(is-1:ie-1)                        ! west
        mask(is:ie,1:4) = b(is:ie,1:4)/=value
        neighb(is:ie) = count( mask(is:ie,1:4), dim=2 )
        where ( neighb(is:ie)>0 .and. a(is:ie)==value )
          a_io(is+n*dk*dk:ie+n*dk*dk,k) = sum( b(is:ie,1:4), mask=mask(is:ie,1:4), dim=2 )/real(neighb(is:ie))
        end where
      end do
      lflag = .false.
      do i = is,ie
        if ( neighb(i)==0 ) then
          nrem = nrem + 1 ! current number of points without a neighbour
          iminb = min(i, iminb)
          imaxb = max(i, imaxb)
          lflag = .true.
        end if
      end do
      if ( lflag ) then
        jminb = min(1, jminb)
        jmaxb = max(1, jmaxb)
      end if
    end if
    do j = max(js,2),min(je,dk-1)
      do k = 1,kx
        a(0)     = b_west(j,k)
        a(dk+1)  = b_east(j,k)
        a(max(is-1,1))  = b_io(max(is-1,1)+(j-1)*dk+n*dk*dk,k)
        a(min(ie+1,dk)) = b_io(min(ie+1,dk)+(j-1)*dk+n*dk*dk,k)
        a(is:ie) = b_io(is+(j-1)*dk+n*dk*dk:ie+(j-1)*dk+n*dk*dk,k)
        b(is:ie,1) = b_io(is+j*dk+n*dk*dk:ie+j*dk+n*dk*dk,k)         ! north
        b(is:ie,2) = b_io(is+(j-2)*dk+n*dk*dk:ie+(j-2)*dk+n*dk*dk,k) ! south
        b(is:ie,3) = a(is+1:ie+1)                                    ! east
        b(is:ie,4) = a(is-1:ie-1)                                    ! west
        mask(is:ie,1:4) = b(is:ie,1:4)/=value
        neighb(is:ie) = count( mask(is:ie,1:4), dim=2)
        where ( neighb(is:ie)>0 .and. a(is:ie)==value )
          a_io(is+(j-1)*dk+n*dk*dk:ie+(j-1)*dk+n*dk*dk,k) = sum( b(is:ie,1:4), mask=mask(is:ie,1:4), dim=2)/real(neighb(is:ie))
        end where
      end do
      lflag = .false.
      do i = is,ie
        if ( neighb(i)==0 ) then
          nrem = nrem + 1 ! current number of points without a neighbour
          iminb = min(i, iminb)
          imaxb = max(i, imaxb)
          lflag = .true.
        end if
      end do
      if ( lflag ) then
        jminb = min(j, jminb)
        jmaxb = max(j, jmaxb)
      end if
    end do
    if ( je==dk ) then
      ! j = dk
      do k = 1,kx
        a(0)     = b_west(dk,k)
        a(dk+1)  = b_east(dk,k)
        a(max(is-1,1))  = b_io(max(is-1,1)-dk+(n+1)*dk*dk,k)
        a(min(ie+1,dk)) = b_io(min(ie+1,dk)-dk+(n+1)*dk*dk,k)
        a(is:ie) = b_io(is-dk+(n+1)*dk*dk:ie-dk+(n+1)*dk*dk,k)
        b(is:ie,1) = b_north(is:ie,k)                                ! north
        b(is:ie,2) = b_io(is-2*dk+(n+1)*dk*dk:ie-2*dk+(n+1)*dk*dk,k) ! south
        b(is:ie,3) = a(is+1:ie+1)                                    ! east
        b(is:ie,4) = a(is-1:ie-1)                                    ! west
        mask(is:ie,1:4) = b(is:ie,1:4)/=value
        neighb(is:ie) = count( mask(is:ie,1:4), dim=2)
        where ( neighb(is:ie)>0 .and. a(is:ie)==value )
          a_io(is-dk+(n+1)*dk*dk:ie-dk+(n+1)*dk*dk,k) = sum( b(is:ie,1:4), mask=mask(is:ie,1:4), dim=2)/real(neighb(is:ie))
        end where
      end do
      lflag = .false.
      do i = is,ie
        if ( neighb(i)==0 ) then
          nrem = nrem + 1 ! current number of points without a neighbour
          iminb = min(i, iminb)
          imaxb = max(i, imaxb)
          lflag = .true.
        end if
      end do
      if ( lflag ) then
        jminb = min(dk, jminb)
        jmaxb = max(dk, jmaxb)
      end if
    end if
    
    imin(n) = iminb
    imax(n) = imaxb
    jmin(n) = jminb
    jmax(n) = jmaxb
  end do
end do
      
return
end subroutine fill_cc4_gather


! *****************************************************************************
! FILE IO ROUTINES

! This version reads and interpolates a surface field
subroutine gethist1(vname,varout)

use cc_mpi             ! CC MPI routines
use infile             ! Input file routines
      
implicit none
      
include 'newmpar.h'    ! Grid parameters
include 'darcdf.h'     ! Netcdf data

integer ier
real, dimension(:), intent(out) :: varout
real, dimension(fwsize) :: ucc
character(len=*), intent(in) :: vname
      
if ( iotest ) then
  ! read without interpolation or redistribution
  call histrd1(iarchi,ier,vname,ik,varout,ifull)
else if ( fnresid==1 ) then
  ! use bcast method for single input file
  ! requires interpolation and redistribution
  call histrd1(iarchi,ier,vname,ik,ucc,6*ik*ik,nogather=.false.)
  call doints1_gather(ucc, varout)
else
  ! use RMA method for multiple input files
  ! requires interpolation and redistribution
  call histrd1(iarchi,ier,vname,ik,ucc,6*ik*ik,nogather=.true.)
  call doints1_nogather(ucc, varout)
end if ! iotest

return
end subroutine gethist1

! This version reads, fills and interpolates a surface field
subroutine fillhist1(vname,varout,mask_a,filllimit)
      
use cc_mpi             ! CC MPI routines
use infile             ! Input file routines
      
implicit none
      
include 'newmpar.h'    ! Grid parameters
include 'darcdf.h'     ! Netcdf data
      
integer ier
real, intent(in), optional :: filllimit
real, dimension(:), intent(out) :: varout
real, dimension(fwsize) :: ucc
logical, dimension(:), intent(in) :: mask_a
character(len=*), intent(in) :: vname
      
if ( iotest ) then
  ! read without interpolation or redistribution
  call histrd1(iarchi,ier,vname,ik,varout,ifull)
else if ( fnresid==1 ) then
  ! use bcast method for single input file
  ! requires interpolation and redistribution
  call histrd1(iarchi,ier,vname,ik,ucc,6*ik*ik,nogather=.false.)
  if ( present(filllimit) ) then
    where ( ucc(:)>=filllimit )
      ucc(:) = 999.
    end where
  end if  
  call fill_cc1_gather(ucc,mask_a)
  call doints1_gather(ucc, varout)
else
  ! use RMA method for multiple input files
  ! requires interpolation and redistribution
  call histrd1(iarchi,ier,vname,ik,ucc,6*ik*ik,nogather=.true.)
  if ( present(filllimit) ) then
    where ( ucc(:)>=filllimit )
      ucc(:) = 999.
    end where
  end if  
  call fill_cc1_nogather(ucc,mask_a)
  call doints1_nogather(ucc, varout)
end if ! iotest
      
return
end subroutine fillhist1

! This version reads 3D fields
subroutine gethist4(vname,varout,kx)

use cc_mpi             ! CC MPI routines
use infile             ! Input file routines
      
implicit none
      
include 'newmpar.h'    ! Grid parameters
include 'darcdf.h'     ! Netcdf data

integer, intent(in) :: kx
integer ier
real, dimension(:,:), intent(out) :: varout
real, dimension(fwsize,kx) :: ucc
character(len=*), intent(in) :: vname

if ( iotest ) then
  ! read without interpolation or redistribution
  call histrd4(iarchi,ier,vname,ik,kx,varout,ifull)
else if ( fnresid==1 ) then
  ! use bcast method for single input file
  ! requires interpolation and redistribution
  call histrd4(iarchi,ier,vname,ik,kx,ucc,6*ik*ik,nogather=.false.)
  call doints4_gather(ucc, varout)
else
  ! use RMA method for multiple input files
  ! requires interpolation and redistribution
  call histrd4(iarchi,ier,vname,ik,kx,ucc,6*ik*ik,nogather=.true.)
  call doints4_nogather(ucc,varout)
end if ! iotest
      
return
end subroutine gethist4   

! This version reads, fills a 3D field for the ocean
subroutine fillhist4(vname,varout,kx,mask_a,filllimit)
  
use cc_mpi             ! CC MPI routines
use infile             ! Input file routines
      
implicit none
      
include 'newmpar.h'    ! Grid parameters
include 'darcdf.h'     ! Netcdf data
      
integer, intent(in) :: kx
integer ier
real, intent(in), optional :: filllimit
real, dimension(:,:), intent(out) :: varout
real, dimension(fwsize,kx) :: ucc
logical, dimension(:), intent(in) :: mask_a
character(len=*), intent(in) :: vname

if ( iotest ) then
  ! read without interpolation or redistribution
  call histrd4(iarchi,ier,vname,ik,kx,varout,ifull)
else if ( fnresid==1 ) then
  ! use bcast method for single input file
  ! requires interpolation and redistribution
  call histrd4(iarchi,ier,vname,ik,kx,ucc,6*ik*ik,nogather=.false.)
  if ( present(filllimit) ) then
    where ( ucc(:,:)>=filllimit )
      ucc(:,:) = 999.
    end where
  end if
  call fill_cc4_gather(ucc,mask_a)
  call doints4_gather(ucc, varout)
else
  ! use RMA method for multiple input files
  ! requires interpolation and redistribution
  call histrd4(iarchi,ier,vname,ik,kx,ucc,6*ik*ik,nogather=.true.)
  if ( present(filllimit) ) then
    where ( ucc(:,:)>=filllimit )
      ucc(:,:) = 999.
    end where
  end if
  call fill_cc4_nogather(ucc,mask_a)
  call doints4_nogather(ucc, varout)
end if ! iotest

return
end subroutine fillhist4


! *****************************************************************************
! FILE DATA MESSAGE PASSING ROUTINES

! Define RMA windows for distributing file data to processors
subroutine file_wininit

use cc_mpi            ! CC MPI routines
use infile            ! Input file routines

implicit none

include 'newmpar.h'   ! Grid parameters

#ifdef usempi3
integer, dimension(:,:,:,:), pointer :: procarray
integer, dimension(4) :: shsize
integer procarray_win
#else
integer, dimension(ik+4,ik+4,npanels+1,2) :: procarray
#endif

if ( allocated(filemap) ) then
  deallocate( filemap )
end if
if ( allocated(axs_w) ) then
  deallocate( axs_w, ays_w, azs_w )
  deallocate( bxs_w, bys_w, bzs_w )
end if

! No RMA window for single input file
if ( fnresid<=1 ) return

if ( myid==0 ) then
  write(6,*) "Create map for file RMA windows"
end if

#ifdef usempi3
shsize(1) = ik + 4
shsize(2) = ik + 4
shsize(3) = npanels + 1
shsize(4) = 2
call ccmpi_allocshdata(procarray,shsize(1:4),procarray_win)
call ccmpi_shepoch(procarray_win)
if ( node_myid==0 ) then
  call file_wininit_defineprocarray(procarray)
end if
call ccmpi_shepoch(procarray_win)
#else
call file_wininit_defineprocarray(procarray)
#endif

call file_wininit_definefilemap(procarray)

! Distribute fields for vector rotation
if ( myid==0 ) then
  write(6,*) "Distribute vector rotation data to processors reading input files"
end if
    
allocate(axs_w(fwsize), ays_w(fwsize), azs_w(fwsize))
allocate(bxs_w(fwsize), bys_w(fwsize), bzs_w(fwsize))
if ( myid==0 ) then
  call file_distribute(axs_w,axs_a)
  call file_distribute(ays_w,ays_a)
  call file_distribute(azs_w,azs_a)
  call file_distribute(bxs_w,bxs_a)
  call file_distribute(bys_w,bys_a)
  call file_distribute(bzs_w,bzs_a)
else if ( fwsize>0 ) then
  call file_distribute(axs_w)
  call file_distribute(ays_w)
  call file_distribute(azs_w)
  call file_distribute(bxs_w)
  call file_distribute(bys_w)
  call file_distribute(bzs_w)
end if

! Define halo indices for ccmpi_filebounds
if ( myid==0 ) then
  write(6,*) "Setup bounds function for processors reading input files"
end if

call ccmpi_filebounds_setup(procarray,comm_ip,ik)

#ifdef usempi3
call ccmpi_freeshdata(procarray_win)
#endif

if ( myid==0 ) then
  write(6,*) "Finished creating control data for file RMA windows"
end if

return
end subroutine file_wininit

subroutine file_wininit_defineprocarray(procarray)

use cc_mpi            ! CC MPI routines
use infile            ! Input file routines

implicit none

include 'newmpar.h'   ! Grid parameters

integer i, n
integer n_n, n_e, n_s, n_w
integer ip, ipf, jpf, no, ca, cb
integer, dimension(-1:ik+2,-1:ik+2,0:npanels,2), intent(inout) :: procarray ! can be a pointer

! define host process of each input file gridpoint
procarray(-1:ik+2,-1:ik+2,0:npanels,1:2) = -1
do ipf = 0,fnproc/fnresid-1
  do jpf = 1,fnresid
    ip = ipf*fnresid + jpf - 1
    do n = 0,pnpan-1
      no = n - pnoff(ip) + 1
      ca = pioff(ip,no)
      cb = pjoff(ip,no)
      procarray(1+ca:pil+ca,1+cb:pjl+cb,no,1) = jpf - 1 ! processor rank
      procarray(1+ca:pil+ca,1+cb:pjl+cb,no,2) = ipf + 1 ! file rank
    end do
  end do
end do

! update boundaries
do n = 0,npanels
  if ( mod(n,2)==0 ) then
    n_w = mod(n+5, 6)
    n_e = mod(n+2, 6)
    n_n = mod(n+1, 6)
    n_s = mod(n+4, 6)
    do i = 1,ik
      procarray(0,i,n,:)    = procarray(ik,i,n_w,:)
      procarray(-1,i,n,:)   = procarray(ik-1,i,n_w,:)
      procarray(ik+1,i,n,:) = procarray(ik+1-i,1,n_e,:)
      procarray(ik+2,i,n,:) = procarray(ik+1-i,2,n_e,:)
      procarray(i,ik+1,n,:) = procarray(i,1,n_n,:)
      procarray(i,ik+2,n,:) = procarray(i,2,n_n,:)
      procarray(i,0,n,:)    = procarray(ik,ik+1-i,n_s,:)
      procarray(i,-1,n,:)   = procarray(ik-1,ik+1-i,n_s,:)
    end do ! i
    procarray(-1,0,n,:)      = procarray(ik,2,n_w,:)        ! wws
    procarray(0,-1,n,:)      = procarray(ik,ik-1,n_s,:)     ! wss
    procarray(0,0,n,:)       = procarray(ik,1,n_w,:)        ! ws
    procarray(ik+1,0,n,:)    = procarray(ik,1,n_e,:)        ! es  
    procarray(ik+2,0,n,:)    = procarray(ik-1,1,n_e,:)      ! ees 
    procarray(-1,ik+1,n,:)   = procarray(ik,ik-1,n_w,:)     ! wwn
    procarray(0,ik+2,n,:)    = procarray(ik-1,ik,n_w,:)     ! wnn
    procarray(ik+2,ik+1,n,:) = procarray(2,1,n_e,:)         ! een  
    procarray(ik+1,ik+2,n,:) = procarray(1,2,n_e,:)         ! enn  
    procarray(0,ik+1,n,:)    = procarray(ik,ik,n_w,:)       ! wn  
    procarray(ik+1,ik+1,n,:) = procarray(1,1,n_e,:)         ! en  
    procarray(ik+1,-1,n,:)   = procarray(ik,2,n_e,:)        ! ess  
  else
    n_w = mod(n+4, 6)
    n_e = mod(n+1, 6)
    n_n = mod(n+2, 6)
    n_s = mod(n+5, 6)
    do i = 1,ik
      procarray(0,i,n,:)    = procarray(ik+1-i,ik,n_w,:)
      procarray(-1,i,n,:)   = procarray(ik+1-i,ik-1,n_w,:)
      procarray(ik+1,i,n,:) = procarray(1,i,n_e,:)
      procarray(ik+2,i,n,:) = procarray(2,i,n_e,:)
      procarray(i,ik+1,n,:) = procarray(1,ik+1-i,n_n,:)
      procarray(i,ik+2,n,:) = procarray(2,ik+1-i,n_n,:)
      procarray(i,0,n,:)    = procarray(i,ik,n_s,:)
      procarray(i,-1,n,:)   = procarray(i,ik-1,n_s,:)
    end do ! i
    procarray(-1,0,n,:)      = procarray(ik-1,ik,n_w,:)    ! wws
    procarray(0,-1,n,:)      = procarray(2,ik,n_s,:)       ! wss
    procarray(0,0,n,:)       = procarray(ik,ik,n_w,:)      ! ws
    procarray(ik+1,0,n,:)    = procarray(1,1,n_e,:)        ! es
    procarray(ik+2,0,n,:)    = procarray(1,2,n_e,:)        ! ees
    procarray(-1,ik+1,n,:)   = procarray(2,ik,n_w,:)       ! wwn   
    procarray(0,ik+2,n,:)    = procarray(1,ik-1,n_w,:)     ! wnn  
    procarray(ik+2,ik+1,n,:) = procarray(1,ik-1,n_e,:)     ! een  
    procarray(ik+1,ik+2,n,:) = procarray(2,ik,n_e,:)       ! enn  
    procarray(0,ik+1,n,:)    = procarray(1,ik,n_w,:)       ! wn  
    procarray(ik+1,ik+1,n,:) = procarray(1,ik,n_e,:)       ! en  
    procarray(ik+1,-1,n,:)   = procarray(2,1,n_e,:)        ! ess          
  end if     ! if mod(n,2)==0 ..else..
end do       ! n

return
end subroutine file_wininit_defineprocarray

subroutine file_wininit_definefilemap(procarray)

use cc_mpi            ! CC MPI routines

implicit none

include 'newmpar.h'   ! Grid parameters
include 'parm.h'      ! Model configuration

integer mm, iq, idel, jdel, n
integer ncount, iproc, rproc
integer, dimension(-1:ik+2,-1:ik+2,0:npanels,2), intent(in) :: procarray
logical, dimension(-1:nproc-1) :: lproc

! calculate which grid points and input files are needed by this processor
lproc(-1:nproc-1) = .false.
do mm = 1,m_fly
  do iq = 1,ifull
    idel = int(xg4(iq,mm))
    jdel = int(yg4(iq,mm))
    n = nface4(iq,mm)
    ! search stencil of bi-cubic interpolation
    lproc(procarray(idel,  jdel+2,n,1)) = .true.
    lproc(procarray(idel+1,jdel+2,n,1)) = .true.
    lproc(procarray(idel-1,jdel+1,n,1)) = .true.
    lproc(procarray(idel  ,jdel+1,n,1)) = .true.
    lproc(procarray(idel+1,jdel+1,n,1)) = .true.
    lproc(procarray(idel+2,jdel+1,n,1)) = .true.
    lproc(procarray(idel-1,jdel,  n,1)) = .true.
    lproc(procarray(idel  ,jdel,  n,1)) = .true.
    lproc(procarray(idel+1,jdel,  n,1)) = .true.
    lproc(procarray(idel+2,jdel,  n,1)) = .true.
    lproc(procarray(idel,  jdel-1,n,1)) = .true.
    lproc(procarray(idel+1,jdel-1,n,1)) = .true.
  end do
end do
if ( lproc(-1) ) then
  write(6,*) "ERROR: Internal error in file_wininit"
  call ccmpi_abort(-1)
end if

! Construct a map of files to be accessed by MPI_Get
ncount = count(lproc)
allocate( filemap(ncount) )
ncount = 0
do iproc = 0,nproc-1
  ! stagger reading of windows - does this make any difference with active RMA?
  rproc = modulo( myid+iproc, nproc )
  if ( lproc(rproc) ) then
    ncount = ncount + 1
    filemap(ncount) = rproc
  end if
end do

return
end subroutine file_wininit_definefilemap

! Define commuication group for broadcasting file panel data
subroutine splitface

use cc_mpi            ! CC MPI routines

implicit none

include 'newmpar.h'   ! Grid parameters

integer n, colour

! Free any existing comm_face
if ( bcst_allocated ) then
  do n = 0,npanels
    call ccmpi_commfree(comm_face(n))
  end do
  bcst_allocated = .false.
end if

! No split face for multiple input files
if ( fnresid>1 ) return

if ( myid == 0 ) then
  write(6,*) "Create communication groups for Bcast method in onthefly"
end if

do n = 0,npanels
  if ( nfacereq(n) ) then
    colour = 1
  else
    colour = -1 ! undefined
  end if
  call ccmpi_commsplit(comm_face(n),comm_world,colour,myid)
end do
bcst_allocated = .true.

if ( myid==0 ) then
  write(6,*) "Finished initalising Bcast method for onthefly"
end if

return
end subroutine splitface

end module onthefly_m
