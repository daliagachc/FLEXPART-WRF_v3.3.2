!***********************************************************************
!* Copyright 2012,2013                                                *
!* Jerome Brioude, Delia Arnold, Andreas Stohl, Wayne Angevine,       *
!* John Burkhart, Massimo Cassiani, Adam Dingwell, Richard C Easter, Sabine Eckhardt,*
!* Stephanie Evan, Jerome D Fast, Don Morton, Ignacio Pisso,          *
!* Petra Seibert, Gerard Wotawa, Caroline Forster, Harald Sodemann,   *
!*                                                                     *
!* This file is part of FLEXPART WRF                                   *
!*                                                                     *
!* FLEXPART is free software: you can redistribute it and/or modify    *
!* it under the terms of the GNU General Public License as published by*
!* the Free Software Foundation, either version 3 of the License, or   *
!* (at your option) any later version.                                 *
!*                                                                     *
!* FLEXPART is distributed in the hope that it will be useful,         *
!* but WITHOUT ANY WARRANTY; without even the implied warranty of      *
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       *
!* GNU General Public License for more details.                        *
!*                                                                     *
!* You should have received a copy of the GNU General Public License   *
!* along with FLEXPART.  If not, see <http://www.gnu.org/licenses/>.   *
!***********************************************************************

      subroutine readwind_timeav(indj,n,uuh,vvh,wwh)
!**********************************************************************
!                                                                     * 
!             TRAJECTORY MODEL SUBROUTINE READWIND                    *
!                                                                     *
!**********************************************************************
!                                                                     * 
!  April 2012, J. Brioude:                                            *
!             This routine handles the difference in time for         *
!             time-average fields.                                    *
!                                                                     *
!  2015-03-26, A. Dingwell:                                           *
!             Updated calls to read_ncwrfout_gridinfo to match updates*
!             in that subroutine.                                     * 
!                                                                     *
!**********************************************************************
!                                                                     *
! Note:  This is the FLEXPART_WRF version of subroutine readwind.     *
!    The met fields are read from WRF netcdf output files.            *
!    There are many differences from the FLEXPART version.            *
!                                                                     *
! DESCRIPTION:                                                        *
!                                                                     *
! READING OF ECMWF METEOROLOGICAL FIELDS FROM INPUT DATA FILES. THE   *
! INPUT DATA FILES ARE EXPECTED TO BE AVAILABLE IN GRIB CODE          *
!                                                                     *
! INPUT:                                                              *
! indj               indicates number of the wind field to be read in *
! n                  temporal index for meteorological fields (1 to 3)*
!                                                                     *
! IMPORTANT VARIABLES FROM COMMON BLOCK:                              *
!                                                                     *
! wfname             File name of data to be read in                  *
! nx,ny,nuvz,nwz     expected field dimensions                        *
! nlev_ec            number of "T-grid" vertical levels wwf model     *
!                    (the unstaggered "bottom_top" dimension)         *
! uu,vv,ww           wind fields                                      *
! tt,qv              temperature and specific humidity                *
! ps                 surface pressure                                 *
!                                                                     *
!**********************************************************************
!

!      include 'includepar'
!      include 'includecom'
  use par_mod
  use com_mod
  use flexwrf_ncdf_mod
! subr arguments
      integer :: indj, n

      real(kind=4) :: uuh(0:nxmax-1,0:nymax-1,nuvzmax)
      real(kind=4) :: vvh(0:nxmax-1,0:nymax-1,nuvzmax)
      real(kind=4) :: wwh(0:nxmax-1,0:nymax-1,nwzmax)
      real(kind=4) :: uuh2(0:nxmax-1,0:nymax-1,nuvzmax)
      real(kind=4) :: vvh2(0:nxmax-1,0:nymax-1,nuvzmax)
      real(kind=4) :: wwh2(0:nxmax-1,0:nymax-1,nwzmax)

! local variables
      integer,parameter :: ndims_max=4

      integer :: i, idiagaa, ierr, ifn, itime
      integer :: iduma,indj2
      integer :: j, jhhmmss, jyyyymmdd
      integer :: k, kbgn
      integer :: lendim(ndims_max), lendim_exp(ndims_max), &
          lendim_max(ndims_max)
      integer :: levdiff2,deltat,deltat2
      integer :: ndims, ndims_exp
      integer :: n_west_east, n_south_north, n_bottom_top
      integer :: m_grid_id_dum, m_parent_grid_id_dum, &
        m_parent_grid_ratio_dum,  &
        i_parent_start_dum, j_parent_start_dum, &
        map_proj_id_dum,  &
        ext_scalar,pbl_physics,mp_physics_dum,num_land_cat_dum

      real :: dx_met, dy_met
      real :: duma, dumb, dumc, dumd, dume
      real :: dumdz
!      real(kind=4) :: dumarray_aa(nwzmax+1)
!      real(kind=4) :: dumarray_pp(0:nxmax-1,0:nymax-1,nwzmax+1)
      real :: dumarray_aa(nwzmax+1)
      real :: dumarray_pp(0:nxmax-1,0:nymax-1,nwzmax+1)
      real :: ewater_mb, esatwater_mb
      real :: ew      ! this is an external function
      real :: map_stdlon_dum, map_truelat1_dum, map_truelat2_dum
      real :: pint
      real :: toler

!      real(kind=4) :: ewss(0:nxmax-1,0:nymax-1),nsss(0:nxmax-1,0:nymax-1)
      real :: ewss(0:nxmax-1,0:nymax-1),nsss(0:nxmax-1,0:nymax-1)
      real :: plev1,pmean,tv,fu,hlev1,ff10m,fflev1

      real(kind=dp) :: jul,juldate

      character(len=160) :: fnamenc, varname,fnamenc2

      logical :: hflswitch

!
!   get grid info from the wrf netcdf file
!   and check it for consistency against values from gridcheck
!

!        print*,'entering timeav'

      fnamenc = path(2)(1:length(2))//wfname(indj)
      idiagaa = 0

      call read_ncwrfout_gridinfo( ierr, idiagaa, fnamenc, &
        n_west_east, n_south_north, n_bottom_top, & 
        dx_met, dy_met,  &
        m_grid_id_dum, m_parent_grid_id_dum, m_parent_grid_ratio_dum, &
        i_parent_start_dum, j_parent_start_dum, &
        map_proj_id_dum, map_stdlon_dum,  &
        map_truelat1_dum, map_truelat2_dum, &
        ext_scalar,pbl_physics,mp_physics_dum, num_land_cat_dum )
      if (ierr .ne. 0) then
          write(*,9100) 'error getting gridinfor for met file', fnamenc
          stop
      end if

9100  format( / '*** readwind -- ', a )
9110  format( / '*** readwind -- ', a, 1x, i8 / &
        'file = ', a )
9120  format( / '*** readwind -- ', a, 2(1x,i8) / &
        'file = ', a )
9130  format( / '*** readwind -- ', a, 3(1x,i8) / &
        'file = ', a )
9115  format( / '*** readwind -- ', a / a, 1x, i8 / &
        'file = ', a )
9125  format( / '*** readwind -- ', a / a, 2(1x,i8) / &
        'file = ', a )
9135  format( / '*** readwind -- ', a / a, 3(1x,i8) / &
        'file = ', a )

      toler = 2.0e-7

      if (nx .ne. n_west_east) then
          write(*,9100) 'nx not consistent', fnamenc
          stop
      end if
      if (ny .ne. n_south_north) then
          write(*,9100) 'ny not consistent', fnamenc
          stop
      end if
      if (nlev_ec .ne. n_bottom_top) then
          write(*,9100) 'nlev_ec not consistent', fnamenc
          stop
      end if
      if (nwz .ne. n_bottom_top+1) then
          write(*,9100) 'nwz not consistent', fnamenc
          stop
      end if
!     if (nuvz .ne. n_bottom_top+1) then
!         write(*,9100) 'nuvz not consistent', fnamenc
!         stop
!     end if

      if (m_grid_id(0) .ne. m_grid_id_dum) then
          write(*,9100) 'm_grid_id not consistent', fnamenc
          write(*,*) m_grid_id(0), m_grid_id_dum
          stop
      end if
      if (m_parent_grid_id(0) .ne. m_parent_grid_id_dum) then
          write(*,9100) 'm_parent_grid_id not consistent', fnamenc
          stop
      end if
      if (m_parent_grid_ratio(0) .ne. m_parent_grid_ratio_dum) then
          write(*,9100) 'm_parent_grid_ratio not consistent', fnamenc
          stop
      end if
      if (i_parent_start(0) .ne. i_parent_start_dum) then
          write(*,9100) 'i_parent_start not consistent', fnamenc
          stop
      end if
      if (j_parent_start(0) .ne. j_parent_start_dum) then
          write(*,9100) 'j_parent_start not consistent', fnamenc
          stop
      end if

      if (abs(dx - dx_met) .gt. toler*abs(dx)) then
          write(*,9100) 'dx not consistent', fnamenc
          stop
      end if
      if (abs(dy - dy_met) .gt. toler*abs(dy)) then
          write(*,9100) 'dy not consistent', fnamenc
          stop
      end if


! READ THE FIRST FILE
       if (ldirect.eq.1) indj2=indj
       if (ldirect.eq.-1) indj2=indj
       deltat=wfdt(indj2)
      fnamenc = path(2)(1:length(2))//wfname(indj2)
! locate the date/time in the file
      itime = 0
1100  itime = itime + 1
      call read_ncwrfout_1datetime( ierr, fnamenc, &
          itime, jyyyymmdd, jhhmmss )
      if (ierr .eq. -1) then
          write(*,9100) 'error reading time from met file', fnamenc
          stop
      else if (ierr .ne. 0) then
          write(*,9125) 'unable to locate date/time in met file',  &
              'indj, itime =', indj2, itime, fnamenc
          stop
      else 
          jul = juldate( jyyyymmdd, jhhmmss )
          duma = (jul-bdate)*86400.
          iduma = nint(duma)
          if (iduma .ne. wftime(indj2)) goto 1100
      end if
      if (option_verbose.eq.1) then
      write(*,*) 
      write(*,*) 'readwind processing wrfout file ='
      write(*,*) fnamenc
      write(*,*) 'itime, ymd, hms =', itime, jyyyymmdd, jhhmmss
      endif
      kbgn = 1 + add_sfc_level

      if (wind_option.eq.1) varname = 'AVGFLX_RUM'
      do i = 1, ndims_max
          lendim_exp(i) = 0
          lendim_max(i) = 1
      end do
      lendim_exp(1) = nx+1
      lendim_max(1) = nxmax
      lendim_exp(2) = ny
      lendim_max(2) = nymax
      lendim_exp(3) = nuvz-add_sfc_level
      lendim_max(3) = nuvzmax
      ndims_exp = 4
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, uuh(0,0,kbgn), &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of U', fnamenc
      if (wind_option.le.0) print*,'you asked snapshot winds'
      if (wind_option.eq.1) print*,'you asked mean winds'
        print*,'change wind_option'

          stop
      end if


! v wind velocity
!   the wrf output file contains (nuvz-add_sfc_level) levels
!   read the data into k=kbgn,nuvz
!   (interpolate it from "V-grid" to "T-grid" later)
      if (wind_option.eq.1) varname = 'AVGFLX_RVM'
      lendim_exp(1) = nx
      lendim_max(1) = nxmax
      lendim_exp(2) = ny+1
      lendim_max(2) = nymax
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, vvh(0,0,kbgn), &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of V', fnamenc
      if (wind_option.eq.0) print*,'you asked snapshot winds'
      if (wind_option.eq.1) print*,'you asked mean winds'
        print*,'change wind_option'
          stop
      end if


! w wind velocity
!   this is on the "W-grid", and 
!   the wrf output file contains nwz levels, so no shifting needed
      if (wind_option.eq.1) varname = 'AVGFLX_WWM'
!     print*,'varname',varname
      lendim_exp(1) = nx
      lendim_max(1) = nxmax
      lendim_exp(2) = ny
      lendim_max(2) = nymax
      lendim_exp(3) = nwz
      lendim_max(3) = nwzmax
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, wwh, &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of W', fnamenc
      if (wind_option.eq.0) print*,'you asked snapshot winds'
      if (wind_option.eq.1) print*,'you asked mean winds'
        print*,'change wind_option'
          stop
      end if

! READ THE SECOND FILE and average
! deltat must be equal to deltat2, otherwise the time-average wind cannot be
! fixed.
! deltat is assumed to be the same than the WRF output.

       if (ldirect.eq.1) indj2=indj+1
       if (ldirect.eq.-1) indj2=indj-1
       deltat2=wfdt(indj2)
      fnamenc = path(2)(1:length(2))//wfname(indj2)

!        print*,'test timeav',deltat,deltat2,indj2,numbwf
     if (deltat.eq.deltat2 .and. indj2.ge.1 .and.indj2.le.numbwf) then
!        print*,'reading second file in timeav',deltat,deltat2,indj2
! locate the date/time in the file
      itime = 0
1101  itime = itime + 1
      call read_ncwrfout_1datetime( ierr, fnamenc, &
          itime, jyyyymmdd, jhhmmss )
      if (ierr .eq. -1) then
          write(*,9100) 'error reading time from met file', fnamenc
          stop
      else if (ierr .ne. 0) then
          write(*,9125) 'unable to locate date/time in met file',  &
              'indj, itime =', indj2, itime, fnamenc
          stop
      else 
          jul = juldate( jyyyymmdd, jhhmmss )
          duma = (jul-bdate)*86400.
          iduma = nint(duma)
          if (iduma .ne. wftime(indj2)) goto 1101
      end if
      if (option_verbose.eq.1) then
      write(*,*) 
      write(*,*) 'readwind processing wrfout file ='
      write(*,*) fnamenc
      write(*,*) 'itime, ymd, hms =', itime, jyyyymmdd, jhhmmss
       endif
      kbgn = 1 + add_sfc_level

      if (wind_option.eq.1) varname = 'AVGFLX_RUM'
!     print*,'varname',varname
      do i = 1, ndims_max
          lendim_exp(i) = 0
          lendim_max(i) = 1
      end do
      lendim_exp(1) = nx+1
      lendim_max(1) = nxmax
      lendim_exp(2) = ny
      lendim_max(2) = nymax
      lendim_exp(3) = nuvz-add_sfc_level
      lendim_max(3) = nuvzmax
      ndims_exp = 4
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, uuh2(0,0,kbgn), &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of U', fnamenc
      if (wind_option.le.0) print*,'you asked snapshot winds'
      if (wind_option.eq.1) print*,'you asked mean winds'
        print*,'change wind_option'

          stop
      end if


! v wind velocity
!   the wrf output file contains (nuvz-add_sfc_level) levels
!   read the data into k=kbgn,nuvz
!   (interpolate it from "V-grid" to "T-grid" later)
      if (wind_option.eq.1) varname = 'AVGFLX_RVM'
!     print*,'varname',varname
      lendim_exp(1) = nx
      lendim_max(1) = nxmax
      lendim_exp(2) = ny+1
      lendim_max(2) = nymax
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, vvh2(0,0,kbgn), &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of V', fnamenc
      if (wind_option.eq.0) print*,'you asked snapshot winds'
      if (wind_option.eq.1) print*,'you asked mean winds'
        print*,'change wind_option'
          stop
      end if


! w wind velocity
!   this is on the "W-grid", and 
!   the wrf output file contains nwz levels, so no shifting needed
      if (wind_option.eq.1) varname = 'AVGFLX_WWM'
!     print*,'varname',varname
      lendim_exp(1) = nx
      lendim_max(1) = nxmax
      lendim_exp(2) = ny
      lendim_max(2) = nymax
      lendim_exp(3) = nwz
      lendim_max(3) = nwzmax
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, wwh2, &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of W', fnamenc
      if (wind_option.eq.0) print*,'you asked snapshot winds'
      if (wind_option.eq.1) print*,'you asked mean winds'
        print*,'change wind_option'
          stop
      end if

!        print*,'start modif timeav'
         do k=1,nuvzmax
        do j=0,nymax-1
       do i=0,nxmax-1
           uuh(i,j,k)=0.5*(uuh(i,j,k)+uuh2(i,j,k))
           vvh(i,j,k)=0.5*(vvh(i,j,k)+vvh2(i,j,k))
         enddo
        enddo
       enddo
         do k=1,nwzmax
        do j=0,nymax-1
       do i=0,nxmax-1
           wwh(i,j,k)=0.5*(wwh(i,j,k)+wwh2(i,j,k))
         enddo
        enddo
       enddo
!        print*,'end modif timeav'


      endif ! test on deltat

!         print*,'out of entering timeav'

      return    
      end subroutine readwind_timeav

