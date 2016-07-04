!Conformal Cubic Atmospheric Model
    
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


      integer meso,nrungcm,newtop,                                  &
     &        kountr,                                  &
     &        nhstest,nspecial,newrough,newsoilm,nsib,nsoil,       &
     &        ntaft,ntsea,ntsur,ntsur2,lgwd,newztsea,nglacier,       &
     &        kbotu,nbox,   &
     &        ktau,ndi,ndi2,ntau,nperavg,nperday,nmaxpr,nlv,    &
     &        ia,ib,ja,jb,id,jd,idjd,                                    &
     &        io_clim,io_in,io_out,io_rest,io_spec,                      &
     &        nwt,nqg,nrun,nextout,nclim,m_fly,nsemble,tblock,tbave,     &
     &        nud_sss,   &
     &        mloalpha,nud_ouv,nud_sfh,kblock,ccycle,   &
     &        nriver
      real qgmin,                                                        &
     &     aleadfr,snmin,tss_sh,charnock,chn10,zobgin,   &
     &     rlongdn,rlongdx,rlatdn,rlatdx,ds,dt,dtin,timea,panfg,panzo,   &
     &     bpyear,cgmap_offset,cgmap_scale
      logical diag,localhist,amipo3
      common/parm1/meso,nrungcm,newtop,bpyear,      &
     &  qgmin     ! min value, esp. for stratosphere [1.e-6]

      common/parmradn/kountr,amipo3   

      common/parmvmix/cgmap_offset,cgmap_scale

      common/parmtest/nhstest,nsemble,nspecial,rlongdn,rlongdx,    &
     &                rlatdn,rlatdx

      common/parmsfce/newrough,newsoilm,nsib,nsoil,ntsea,ntsur,ntsur2,   &
     &                lgwd,newztsea,aleadfr,snmin,       &
     &                tss_sh,nglacier,charnock,chn10,zobgin,ntaft,       &
     &                panfg,panzo,ccycle,nriver

      common/parmnudg/kbotu,nbox,   &
     &                nud_sss,      &
     &                mloalpha,nud_ouv,nud_sfh,kblock

      common/parmtime/ktau,ntau,nperavg,nperday,ds,dt,dtin,timea,nmaxpr, &
     &                diag,nlv,ia,ib,ja,jb,id,jd,idjd,ndi,ndi2

      common/parmio/io_clim,io_in,io_out,io_rest,io_spec,                &  ! type of I/O
     &            nwt,nqg,nrun,nextout,nclim,m_fly,tblock,tbave,         &
     &            localhist

