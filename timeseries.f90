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
    
module timeseries

!     rml 25/08/03 declarations from sflux
implicit none
integer, save :: indextime,ntsfreq,ngrdpts,ngrdpts1,n3d,n2d
real(kind=8), save :: tstime
integer, pointer, dimension(:,:), save :: listijk
logical, pointer, dimension(:), save :: writesurf
integer, pointer, dimension(:), save :: tsid
character(len=10), allocatable, dimension(:), save :: varname3,varname2
!     rml 25/11/03 declarations from sflux
integer, save :: indship,nshippts,inshipid(4),outshipid(3)
integer, save :: nshipout
integer, allocatable, dimension(:), save :: shipdate,shiptime
!     rml 19/09/07 conversion of u,v from ccam grid to zonal/merid
real, dimension(:), save, allocatable :: costh,sinth


contains
! *********************************************************************

subroutine init_ts(ngas,dt)
!~ use tracermodule, only : sitefile,shipfile
implicit none
integer ngas
real dt

!~ if (sitefile.ne.'') call readsitelist(ngas)
!~ if (shipfile.ne.'') call readshiplist(ngas,dt)

return
end subroutine init_ts

! ********************************************************************
subroutine write_ts(ktau,ntau,dt)
!~ use tracermodule, only : sitefile,shipfile
implicit none
include 'dates.h'
integer jyear,jmonth,jday,jhour,jmin
integer ktau,ntau,mstart,mins
real dt
integer ndoy(12)   ! days from beginning of year (1st Jan is 0)
data ndoy/ 0,31,59,90,120,151,181,212,243,273,304,334/

jyear=kdate/10000
jmonth=(kdate-jyear*10000)/100
jday=kdate-jyear*10000-jmonth*100
jhour=ktime/100
jmin=ktime-jhour*100
mstart=24*60*(ndoy(jmonth)+jday-1) + 60*jhour + jmin ! mins from start of year
!     mtimer contains number of minutes since the start of the run.
mins = mtimer + mstart

return

end subroutine write_ts
! *********************************************************************
subroutine readsitelist(ntrac)
!
!     rml 25/08/03 subroutine to read file containing list of sites for 
!     timeseries output and to open netcdf file for output and 
!     write dimensions etc.
!     All processors read this list and select the points in their own
!     region.
!
use cc_mpi, only : myid, indv_mpi, fproc, ipan, ccmpi_abort
use infile
!     rml 19/09/07 add tracname so that can be written to ts output file
use vecsuv_m
use xyzinfo_m
implicit none
integer kount,kountprof,n,kount500
integer ntrac
integer griddim,ijkdim,timedim(1),tracdim,gridid,dims(3)
!     rml 19/09/07 addition of tracer name array in output file
integer lendim, tracnamid
integer surfdim,gridsurfid
integer gridorderid, surforderid
integer, allocatable, save, dimension(:,:) :: templist
integer, allocatable, save, dimension(:) :: gridorder, surforder
integer i,i1,k,ntop
character(len=80) outfile
character(len=8) chtemp
character(len=80) head
integer :: ig, jg, tmpval ! Better name for this
include 'dates.h'
include 'newmpar.h'  ! kl
integer ip, nface, ii, jj, iqg, istn, jstn
!     rml 19/09/07 arrays needed for wind conversion
include 'parmgeom.h'    ! rlong0, rlat0
include 'const_phys.h'  ! pi
logical windconv
real coslong,sinlong,coslat,sinlat,polenx,poleny,polenz
real zonx,zony,zonz,den
integer iq

windconv = .false.

!     read file of site locations for timeseries output

allocate(templist(ngrdpts1,3))
if ( nproc > 1 ) then
  allocate(surforder(ngrdpts1) )
end if
kountprof=0
kount500=0

! Read list of point locations, selecting those that are in
! this processor's region.
! Perhaps need to have an additional netcdf variable for ip, so
! that the original order can be reconstructed from the multiple files.
k = 0


ngrdpts1 = k ! Reset to the actual number I have

!     Read in any additional variables to output besides tracer
n2d=0
n3d=0


!     rml 19/09/07 fill arrays needed for wind conversion
if (windconv) then
  allocate(costh(ifull),sinth(ifull))
!       code taken from cc2hist 
  coslong=cos(rlong0*pi/180.)
  sinlong=sin(rlong0*pi/180.)
  coslat=cos(rlat0*pi/180.)
  sinlat=sin(rlat0*pi/180.)
  polenx=-coslat
  poleny=0.
  polenz=sinlat
  do iq=1,ifull
!         Set up unit zonal vector components
    zonx = real(poleny*z(iq)-polenz*y(iq))
    zony = real(polenz*x(iq)-polenx*z(iq))
    zonz = real(polenx*y(iq)-poleny*x(iq))
!         Allow for poles by taking max
    den = sqrt( max(zonx**2 + zony**2 + zonz**2,1.e-7) )
    costh(iq) =  (zonx*ax(iq)+zony*ay(iq)+zonz*az(iq))/den
    sinth(iq) = -(zonx*bx(iq)+zony*by(iq)+zonz*bz(iq))/den
  enddo
endif

ngrdpts = ngrdpts1 + kountprof*(kl-1) + kount500*8
allocate(listijk(ngrdpts,3))
allocate(writesurf(ngrdpts))
if ( nproc > 1 ) then
  allocate(gridorder(ngrdpts) )
end if
kount = 0
do k=1,ngrdpts1
!     rml 11/11/05 add '98' option for all levels to ~500 hPa
  if (templist(k,3).ge.98) then
    if (templist(k,3).eq.98) then
      ntop = 9
    else
      ntop = kl
    endif
    do n=1,ntop
      kount = kount + 1
      listijk(kount,1:2) = templist(k,1:2)
      listijk(kount,3) = n
      if (n.eq.1) then
        writesurf(kount) = .true.
      else
        writesurf(kount) = .false.
      endif
      if ( nproc > 1 ) gridorder(kount) = surforder(k)
    enddo
  else
    kount = kount + 1
    listijk(kount,:) = templist(k,:)
    writesurf(kount) = .true.
    if ( nproc > 1 ) gridorder(kount) = surforder(k)
  endif
enddo
if (kount.ne.ngrdpts) then
  write(6,*) 'location file: kount.ne.ngrdpts'
  call ccmpi_abort(-1)
end if
!     deallocate(templist)
allocate(tsid(3+n3d+n2d))
!
!     open netcdf file for writing output
write(chtemp,'(i8)') kdate
if ( nproc == 1 ) then
  outfile = 'ts.'//chtemp(1:4)//'.'//chtemp(5:6)//'.nc'
else
!       Include a 6 digid processor number in the file
  write(outfile,'(a,a,a,a,a, i6.6,a)') 'ts.', chtemp(1:4), '.',chtemp(5:6), '.p', myid, '.nc'
end if
!
!     when multiprocessor, possible to have no grid points in an output file
!     in that case, don't write the file
if (ngrdpts>0) then
  call ccnf_create(outfile,tsid(1))
  call ccnf_nofill(tsid(1))
!       define dimensions
  call ccnf_def_dim(tsid(1),'gridpts',ngrdpts,griddim)
  call ccnf_def_dim(tsid(1),'surfpts',ngrdpts1,surfdim)
  call ccnf_def_dim(tsid(1),'ijk',3,ijkdim)
  call ccnf_def_dim(tsid(1),'tracers',ntrac,tracdim)
  call ccnf_def_dimu(tsid(1),'time',timedim(1))
!     rml 19/09/07 add length dimension for tracname character array
  call ccnf_def_dim(tsid(1),'len',13,lendim)
!       define variables
!       rml 19/09/07 add tracer name variable
  dims(1)=lendim; dims(2)=tracdim
  call ccnf_def_var(tsid(1),'tracer_name','char',2,dims,tracnamid)
  dims(1)=griddim; dims(2)=ijkdim
  call ccnf_def_var(tsid(1),'grid','int',2,dims,gridid)
  dims(1)=surfdim; dims(2)=ijkdim
  call ccnf_def_var(tsid(1),'gridsurf','int',2,dims,gridsurfid)

  if ( nproc > 1 ) then
    dims(1)=griddim
    call ccnf_def_var(tsid(1),'gridorder','int',1,dims,gridorderid)
    dims(1)=surfdim
    call ccnf_def_var(tsid(1),'surforder','int',1,dims,surforderid)
  end if

  call ccnf_def_var(tsid(1),'time','double',1,timedim,tsid(2))
  dims(1)=griddim; dims(2)=tracdim; dims(3)=timedim(1)
  call ccnf_def_var(tsid(1),'concts','float',3,dims,tsid(3))
  do n=1,n3d
    dims(1)=griddim; dims(2)=timedim(1)
    call ccnf_def_var(tsid(1),varname3(n),'float',2,dims(1:2),tsid(3+n))
  enddo
  do n=1,n2d
    if (trim(varname2(n))=='flux') then
      dims(1)=surfdim; dims(2)=tracdim; dims(3)=timedim(1)
      call ccnf_def_var(tsid(1),varname2(n),'float',3,dims,tsid(3+n3d+n))
    else
      dims(1)=surfdim; dims(2)=timedim(1)
      call ccnf_def_var(tsid(1),varname2(n),'float',2,dims(1:2),tsid(3+n3d+n))
    endif
  enddo
!
!     leave define mode
  call ccnf_enddef(tsid(1))
!



!     write grid point arrays
  call ccnf_put_vara(tsid(1),gridid,listijk)
! Need explicit section of templist here, because array may have
! been allocated larger
  call ccnf_put_vara(tsid(1),gridsurfid,templist(:ngrdpts1,:))
  if ( nproc > 1 ) then
   call ccnf_put_vara(tsid(1),gridorderid,gridorder(1:ngrdpts))
   call ccnf_put_vara(tsid(1),surforderid,surforder(1:ngrdpts1))
  end if
  call ccnf_sync(tsid(1))
  deallocate(templist)
endif
!
indextime=1
return
end subroutine readsitelist
!
! **********************************************************************
   
subroutine writetimeseries(ktau,ntau,jyear,mins)

!
!     rml: subroutine to write timeseries data to netcdf file
!  rml 10/11/05: added pressure, surface flux and pblh for TC
!
use arrays_m    ! temp, q, ps
use cable_def_types_mod, only : ncs, ncp ! Used in carbpool.h
use cc_mpi, only : ccmpi_abort
use extraout_m  ! cloud arrays
use infile
use morepbl_m   ! rnet,eg,fg
use pbl_m       ! tss
use prec_m      ! precip
use sigs_m      ! sigma levels for pressure
use soil_m      ! albedo
use soilsnow_m  ! soil temp (tgg)
!~ use vegpar_m    ! rlai
use vvel_m      ! vertical velocity
implicit none
real, dimension(:,:), allocatable, save :: cts
real, dimension(:), allocatable, save :: vts
integer start(3),ncount(3),n,iq,kount,m,jyear,mins
integer ktau,ntau,k
logical surfflux
include 'newmpar.h'    ! dimensions for tr array
include 'const_phys.h' ! grav (for height calc)

real temparr2(il*jl,kl),temparr(il*jl)

if (ngrdpts.eq.0) return

#ifdef outsync
  call ccnf_sync(tsid(1))
#endif

if (mod(ktau,ntsfreq).eq.0) then
  tstime = real(jyear,8) + real(mins,8)/real(365.*24.*60.,8)
  call ccnf_put_vara(tsid(1),tsid(2),indextime,tstime)
  do n=1,ngrdpts
    iq = listijk(n,1) + (listijk(n,2)-1)*il
! rml 23/2/10 add in background that removed at start of run
  enddo
  start(1)=1; start(2)=1; start(3)=indextime
  ncount(1)=ngrdpts; ncount(3)=1
  call ccnf_put_vara(tsid(1),tsid(3),start,ncount,cts)
  deallocate(cts)

  do m=1,n3d
    select case(trim(varname3(m)))
    case ('t') ; temparr2=t(1:ifull,:)
!        rml 19/09/07 fix for u, v - convert to zonal/meridional from ccam grid
    case ('u') 
      do k=1,kl
        temparr2(:,k) = costh*u(1:ifull,k) - sinth*v(1:ifull,k)
      enddo
    case ('v') 
      do k=1,kl
        temparr2(:,k) = sinth*u(1:ifull,k) + costh*v(1:ifull,k)
      enddo
! rml 16/7/10 add output of height of levels for TC-methane
    case ('zg')
      temparr2(:,1) = bet(1)*t(1:ifull,1)/grav + zs(1:ifull)/grav     
      do k=2,kl                                                   
        temparr2(:,k) = temparr2(:,k-1) + (bet(k)*t(1:ifull,k)+betm(k)*t(1:ifull,k-1))/grav               
      enddo                                                       
    case ('qg') ; temparr2=qg(1:ifull,:)
    case ('sdotm') ; temparr2=sdot(:,1:kl)
    case ('sdotp') ; temparr2=sdot(:,2:kl+1)
    case ('pressure')
      do k=1,kl
        temparr2(:,k)=ps(1:ifull)*sig(k)
      enddo
    case default
      write(6,*) varname3(m),' not found'
      call ccmpi_abort(-1)
    end select
    allocate(vts(ngrdpts))
    do n=1,ngrdpts
      iq = listijk(n,1) + (listijk(n,2)-1)*il
      vts(n)=temparr2(iq,listijk(n,3))
    enddo
    start(1)=1; start(2)=indextime
    ncount(1)=ngrdpts; ncount(2)=1
    call ccnf_put_vara(tsid(1),tsid(3+m),start(1:2),ncount(1:2),vts)
    deallocate(vts)
  enddo

  do m=1,n2d
    temparr=0.
    surfflux=.false.
    select case(trim(varname2(m)))
    case ('cloudlo') ; temparr=cloudlo
    case ('cloudmi') ; temparr=cloudmi
    case ('cloudhi') ; temparr=cloudhi
    case ('ps')      ; temparr=ps(1:ifull)
    case ('tss')     ; temparr=tss
    case ('rnet')    ; temparr=rnet
    case ('eg')      ; temparr=eg
    case ('fg')      ; temparr=fg
    case ('alb')     ; temparr=swrsave*albvisnir(:,1)+(1.-swrsave)*albvisnir(:,2) ! MJT cable
    case ('sgsave')  ; temparr=sgsave
    case ('rgsave')  ; temparr=rgsave
    case ('precip')  ; temparr=precip
    case ('tgg4')    ; temparr=tgg(:,4)
    case ('tgg5')    ; temparr=tgg(:,5)
    case ('tgg6')    ; temparr=tgg(:,6)
    !~ case ('rlai')    ; temparr=vlai
    case ('pblh')    ; temparr=pblh
    case ('flux')  
      kount=0
      do n=1,ngrdpts
        if (writesurf(n)) then
           kount=kount+1
           iq = listijk(n,1) + (listijk(n,2)-1)*il
        endif
      enddo 
      start(1)=1; start(2)=1; start(3)=indextime
      ncount(1)=ngrdpts1;ncount(3)=1
      call ccnf_put_vara(tsid(1),tsid(3+n3d+m),start,ncount,cts)
      deallocate(cts)
      surfflux=.true.
   case default
      write(6,*) trim(varname2(m)),' not found'
      call ccmpi_abort(-1)
   end select
          
   if (.not.surfflux) then
      allocate(vts(ngrdpts1))
      kount = 0
      do n=1,ngrdpts
        if (writesurf(n)) then
          kount = kount+1
          iq = listijk(n,1) + (listijk(n,2)-1)*il
          vts(kount)=temparr(iq)
        endif
      enddo
      start(1)=1; start(2)=indextime
      ncount(1)=ngrdpts1; ncount(2)=1
      call ccnf_put_vara(tsid(1),tsid(3+n3d+m),start(1:2),ncount(1:2),vts)
      deallocate(vts)
   endif
  enddo

  indextime=indextime+1
endif
if (ktau.eq.ntau) call ccnf_close(tsid(1))

return
end subroutine writetimeseries

end module timeseries
