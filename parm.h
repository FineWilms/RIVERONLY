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


      integer nrungcm,                                  &
     &        nhstest,newrough,newsoilm,nsib,nsoil,       &
     &        ntsur,nglacier,       &
     &        nbox,   &
     &        ktau,ndi,ntau,nperavg,nperday,nmaxpr,    &
     &        ia,ib,ja,jb,id,jd,idjd,                                    &
     &        io_in,io_out,io_rest,                      &
     &        nwt,nrun,nextout,m_fly,nsemble,     &
     &        kblock,ccycle,   &
     &        nriver
      real snmin,                                                        &
     &     chn10,zobgin,   &
     &     rlongdn,rlongdx,rlatdn,rlatdx,ds,dt,dtin,timea,   &
     &     bpyear
      logical diag,localhist
      common/parm1/nrungcm,bpyear

      common/parmtest/nhstest,nsemble,rlongdn,rlongdx,    &
     &                rlatdn,rlatdx

      common/parmsfce/newrough,newsoilm,nsib,nsoil,ntsur,   &
     &                snmin,       &
     &                nglacier,chn10,zobgin,       &
     &                ccycle,nriver
      common/parmnudg/nbox,   &
     &                kblock

      common/parmtime/ktau,ntau,nperavg,nperday,ds,dt,dtin,timea,nmaxpr, &
     &                diag,ia,ib,ja,jb,id,jd,idjd,ndi

      common/parmio/io_in,io_out,io_rest,                &  ! type of I/O
     &            nwt,nrun,nextout,m_fly,         &
     &            localhist

