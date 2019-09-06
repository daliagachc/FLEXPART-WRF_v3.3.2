!***********************************************************************
!* Copyright 2012,2013                                                *
!* Jerome Brioude, Delia Arnold, Andreas Stohl, Wayne Angevine,       *
!* John Burkhart, Massimo Cassiani, Adam Dingwell, Richard C Easter, Sabine Eckhardt,*
!* Stephanie Evan, Jerome D Fast, Don Morton, Ignacio Pisso,          *
!* Petra Seibert, Gerard Wotawa, Caroline Forster, Harald Sodemann,   *
!* Adam Dingwell,                                                      *
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

      subroutine readwind(indj,n,uuh,vvh,wwh,divh)
!**********************************************************************
!                                                                     * 
!             TRAJECTORY MODEL SUBROUTINE READWIND                    *
!                                                                     *
!**********************************************************************
!                                                                     * 
! AUTHOR:      G. WOTAWA                                              *
! DATE:        1997-08-05                                             *
! LAST UPDATE: 2000-10-17, Andreas Stohl                              *
!                                                                     * 
! Bernd C. Krueger, Feb. 2001:  Variables tth and qvh                 *
!                               (on eta coordinates) in common block  *
!                                                                     * 
! Oct-Dec, 2005: R. Easter.  Major changes for WRF.                   *
!    06-nov-2005 rce - change uuh,vvh dimension back to original      * 
!    16-nov-2005 rce - zzh is shifted like pph,tth                    * 
!                                                                     * 
!    11-June-2007,   W.WANG -- read TKE, change ndims_exp=1 for P_TOP
!    19-Oct -2007,             read tcc, RAINC, RAINNC, CLDFRA,W0AVG
!                              Note RAINC, RAINNC are accumulated prec
!    Dec 2011, J Brioude: modifications, notably for the mean wind
!
! D. Arnold May 2012: quick fix - for CLDFRA
! CLDFRA was  for sub-grid variability of precipitation
! within a cell from  Hertel et al., 1995 (grid cells of approx. 150 km)
! WRF will be typically be run at higher resolutions, therefore
! we simply skip this condition and consider the maximum fraction
! no more reading of CLDFRA needed, initialized to 0.
! Same modifications for the nested domains!
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
!     real(kind=4) :: urot(0:nxmax-1,0:nymax-1,nuvzmax)
!     real(kind=4) :: vrot(0:nxmax-1,0:nymax-1,nuvzmax)
      real(kind=4) :: divh(0:nxmax-1,0:nymax-1,nuvzmax)
      real(kind=4) :: mu(0:nxmax-1,0:nymax-1,1),mu2
      real(kind=4) :: mub(0:nxmax-1,0:nymax-1,1)
!     real :: utrue1,vtrue1,utrue2,vtrue2,dumy
!      real(kind=4) :: m_u(0:nxmax-1,0:nymax-1,1)
!      real(kind=4) :: m_v(0:nxmax-1,0:nymax-1,1)
!      real(kind=4) :: m_w(0:nxmax-1,0:nymax-1,1)

!     real :: uuh(0:nxmax-1,0:nymax-1,nuvzmax)
!     real :: vvh(0:nxmax-1,0:nymax-1,nuvzmax)
!     real :: urot(0:nxmax-1,0:nymax-1,nuvzmax)
!     real :: vrot(0:nxmax-1,0:nymax-1,nuvzmax)
!     real :: wwh(0:nxmax-1,0:nymax-1,nwzmax)
!     real :: divh(0:nxmax-1,0:nymax-1,nuvzmax)
!     real :: mu(0:nxmax-1,0:nymax-1,1),mu2
!     real :: mub(0:nxmax-1,0:nymax-1,1)
!     real :: m_u(0:nxmax-1,0:nymax-1,1)
!     real :: m_v(0:nxmax-1,0:nymax-1,1)
!     real :: m_x(0:nxmax-1,0:nymax-1,1)
!     real :: m_y(0:nxmax-1,0:nymax-1,1)
!     real :: m_z(0:nxmax-1,0:nymax-1,1)

! local variables
      integer,parameter :: ndims_max=4

      integer :: i, idiagaa, ierr, ifn, itime, iclass
      integer :: iduma
      integer :: j, jhhmmss, jyyyymmdd
      integer :: k, kbgn
      integer :: lendim(ndims_max), lendim_exp(ndims_max), &
          lendim_max(ndims_max)
      integer :: levdiff2
      integer :: ndims, ndims_exp
      integer :: n_west_east, n_south_north, n_bottom_top
      integer :: m_grid_id_dum, m_parent_grid_id_dum, &
        m_parent_grid_ratio_dum,  &
        i_parent_start_dum, j_parent_start_dum, &
        map_proj_id_dum,  &
        ext_scalar,pbl_physics,mp_physics_dum, num_land_cat

      real :: dx_met, dy_met
      real :: duma, dumb, dumc, dumd, dume
      real :: dumdz
!      real(kind=4) :: dumarray_aa(nwzmax+1)
!      real(kind=4) :: dumarray_pp(0:nxmax-1,0:nymax-1,nwzmax+1)
      real :: dumarray_aa(nwzmax+1)
      real :: dumarray_pp(0:nxmax-1,0:nymax-1,nwzmax+1)
      real(kind=4) :: landuse_wrf(0:nxmax-1,0:nymax-1) !land use cat. from WRF
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
      fnamenc = path(2)(1:length(2))//wfname(indj)
      idiagaa = 0

      call read_ncwrfout_gridinfo( ierr, idiagaa, fnamenc, &
        n_west_east, n_south_north, n_bottom_top, & 
        dx_met, dy_met,  &
        m_grid_id_dum, m_parent_grid_id_dum, m_parent_grid_ratio_dum, &
        i_parent_start_dum, j_parent_start_dum, &
        map_proj_id_dum, map_stdlon_dum,  &
        map_truelat1_dum, map_truelat2_dum, &
        ext_scalar,pbl_physics,mp_physics_dum, num_land_cat)
      if (ierr .ne. 0) then
          write(*,9100) 'error getting gridinfor for met file', fnamenc
          stop
      end if

9100  format( / '*** readwind -- ', a )
9110   format( / '*** readwind -- ', a, 1x, i8 / &
     'file = ', a )
9120 format( / '*** readwind -- ', a, 2(1x,i8) / &
       'file = ', a )
9130 format( / '*** readwind -- ', a, 3(1x,i8) / &
       'file = ', a )
9115 format( / '*** readwind -- ', a / a, 1x, i8 / &
       'file = ', a )
9125 format( / '*** readwind -- ', a / a, 2(1x,i8) / &
       'file = ', a )
9135 format( / '*** readwind -- ', a / a, 3(1x,i8) / &
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
!      if (nuvz .ne. n_bottom_top+1) then
!          write(*,9100) 'nuvz not consistent', fnamenc
!          stop
!      end if

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
              'indj, itime =', indj, itime, fnamenc
          stop
      else 
          jul = juldate( jyyyymmdd, jhhmmss )
          duma = (jul-bdate)*86400.
          iduma = nint(duma)
          if (iduma .ne. wftime(indj)) goto 1100
      end if
      if (option_verbose.eq.1) then
      write(*,*) 
      write(*,*) 'readwind processing wrfout file ='
      write(*,*) fnamenc
      write(*,*) 'itime, ymd, hms =', itime, jyyyymmdd, jhhmmss

      endif
! read eta_w_wrf, eta_u_wrf, p_top_wrf, ylat2d, xlon2d from the 
! netcdf wrfout file and compare to those from the 1st met. file

      varname = 'ZNW'
      do i = 1, ndims_max
          lendim_exp(i) = 0
          lendim_max(i) = 1
      end do
      lendim_exp(1) = nwz
      lendim_max(1) = nwzmax+1
      ndims_exp = 2
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, dumarray_aa, &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
        fnamenc2='wrfout_d03_zn.nc'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc2, &
          varname, dumarray_aa, &
          itime, &
          ndims, ndims_exp, ndims_max, &
          lendim, lendim_exp, lendim_max )  

      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of ZNW', fnamenc
          stop
      end if
      end if
      do k = 1, nwz
          if (abs(eta_w_wrf(k) - dumarray_aa(k))  &
                  .gt. toler*abs(eta_w_wrf(k))) then
              write(*,9100) 'eta_w_wrf not consistent', fnamenc
              stop
          end if
      end do

      varname = 'ZNU'
      lendim_exp(1) = nwz-1
      lendim_max(1) = nwzmax+1
      ndims_exp = 2
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, dumarray_aa, &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
        fnamenc2='wrfout_d03_zn.nc'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc2, &
          varname, dumarray_aa, &
          itime, &
          ndims, ndims_exp, ndims_max, &
          lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then

          write(*,9100) 'error doing ncread of ZNU', fnamenc
          stop
      end if
      end if
      do k = 1, nwz-1
          if (abs(eta_u_wrf(k) - dumarray_aa(k))  &
                  .gt. toler*abs(eta_u_wrf(k))) then
              write(*,9100) 'eta_u_wrf not consistent', fnamenc
              stop
          end if
      end do

!      varname = 'P_TOP'
!      lendim_exp(1) = 1
!      lendim_max(1) = 6
!      ndims_exp = 2
!      if (ext_scalar .lt. 0) ndims_exp = 1
!      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
!      	  varname, duma, &
!      	  itime, &
!      	  ndims, ndims_exp, ndims_max, &
!      	  lendim, lendim_exp, lendim_max )
!      if (ierr .ne. 0) then
!          write(*,9100) 'error doing ncread of P_TOP', fnamenc
!          stop
!      end if
!      if (abs(p_top_wrf - duma) .gt. toler*abs(p_top_wrf)) then
!          write(*,9100) 'p_top_wrf not consistent', fnamenc
!          stop
!      end if

      varname = 'XLAT'
      lendim_exp(1) = nx
      lendim_max(1) = nxmax
      lendim_exp(2) = ny
      lendim_max(2) = nymax
      ndims_exp = 3
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, dumarray_pp, &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,*)
          write(*,9100) 'error doing ncread of XLAT', fnamenc
          stop
      end if
      toler = 1.0e-6
      do j = 0, ny-1
      do i = 0, nx-1
          if (abs(ylat2d(i,j) - dumarray_pp(i,j,1)) .gt.  &
                              toler*abs(ylat2d(i,j))) then
              write(*,9100) 'ylat2d not consistent', fnamenc
              write(*,'(a,2i5,2f16.6)') 'i,j,ylats =', i, j, &
                      ylat2d(i,j), dumarray_pp(i,j,1)
              stop
          end if
      end do
      end do

      varname = 'XLONG'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, dumarray_pp, &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,*)
          write(*,9100) 'error doing ncread of XLONG', fnamenc
          stop
      end if
      do j = 0, ny-1
      do i = 0, nx-1
          if (abs(xlon2d(i,j) - dumarray_pp(i,j,1)) .gt.  &
                              toler*abs(xlon2d(i,j))) then
              write(*,9100) 'xlon2d not consistent', fnamenc
              write(*,'(a,2i5,2f16.6)') 'i,j,xlons =', i, j, &
                      xlon2d(i,j), dumarray_pp(i,j,1)
              stop
          end if
      end do
      end do


!
!
! now read the data fields for current time
! the following are read from ecmwf met files
!       U VELOCITY
!       V VELOCITY
!       W VELOCITY
!       TEMPERATURE
!       SPEC. HUMIDITY  
!       SURF. PRESS.
!       SEA LEVEL PRESS.
!       10 M U VELOCITY
!       10 M V VELOCITY
!       2 M TEMPERATURE
!       2 M DEW POINT  
!       SNOW DEPTH
!       CLOUD COVER
!       LARGE SCALE PREC.
!       CONVECTIVE PREC.
!       SENS. HEAT FLUX
!       SOLAR RADIATION
!       EW SURFACE STRESS
!       NS SURFACE STRESS
!       ECMWF OROGRAPHY
!       STANDARD DEVIATION OF OROGRAPHY
!       ECMWF LAND SEA MASK
!
!
      hflswitch=.false.
      strswitch=.false.
      levdiff2=nlev_ec-nwz+1

      kbgn = 1 + add_sfc_level
! at this point
!   if add_sfc_level=1, then nuvz=nwz   and kbgn=2
!   if add_sfc_level=0, then nuvz=nwz-1 and kbgn=1

! u wind velocity
!   the wrf output file contains (nuvz-add_sfc_level) levels
!   read the data into k=kbgn,nuvz
!   (interpolate it from "U-grid" to "T-grid" later)
      if (wind_option.le.0) varname = 'U'
      if (wind_option.eq.1) varname = 'AVGFLX_RUM'
      if (wind_option.eq.2) varname = 'U'

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
      idiagaa=0
      if (time_option.eq.0) then
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

      endif !test on time_option

! v wind velocity
!   the wrf output file contains (nuvz-add_sfc_level) levels
!   read the data into k=kbgn,nuvz
!   (interpolate it from "V-grid" to "T-grid" later)
      if (wind_option.le.0) varname = 'V'
      if (wind_option.eq.1) varname = 'AVGFLX_RVM'
      if (wind_option.eq.2) varname = 'V'
      lendim_exp(1) = nx
      lendim_max(1) = nxmax
      lendim_exp(2) = ny+1
      lendim_max(2) = nymax
      if (time_option.eq.0) then
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

      endif !test on time_option

! w wind velocity
!   this is on the "W-grid", and 
!   the wrf output file contains nwz levels, so no shifting needed
      if (wind_option.le.0) varname = 'W'
      if (wind_option.eq.1) varname = 'AVGFLX_WWM'
      if (wind_option.eq.2) varname = 'WW'
!     print*,'varname',varname
      lendim_exp(1) = nx
      lendim_max(1) = nxmax
      lendim_exp(2) = ny
      lendim_max(2) = nymax
      lendim_exp(3) = nwz
      lendim_max(3) = nwzmax
      if (time_option.eq.0) then
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

      endif !test on time_option

! pressure - read base state and perturbation pressure,
!     then combine
!   the wrf output file contains (nuvz-add_sfc_level) levels
!   read the data into k=kbgn,nuvz
      varname = 'PB'
      lendim_exp(3) = nuvz-add_sfc_level
      lendim_max(3) = nuvzmax
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, pph(0,0,kbgn,n), &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of PB', fnamenc
          stop
      end if

      varname = 'P'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, dumarray_pp(0,0,kbgn), &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of P', fnamenc
          stop
      end if

      do k = kbgn, nuvz
      do j = 0, nymin1
      do i = 0, nxmin1
          pph(i,j,k,n) = pph(i,j,k,n) + dumarray_pp(i,j,k)
      end do
      end do
      end do


! height - read base state and perturbation geopotential,
!     then combine and divide by gravity
!   these are on the "W-grid", and 
!     the wrf output file contains nwz levels
!   shift them also so they will be consistent with pph
      varname = 'PHB'
      lendim_exp(3) = nwz
      lendim_max(3) = nwzmax+1
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, zzh(0,0,kbgn,n), &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of PB', fnamenc
          stop
      end if

      varname = 'PH'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, dumarray_pp(0,0,kbgn), &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of P', fnamenc
          stop
      end if

      do k = kbgn, nwz+add_sfc_level
      do j = 0, nymin1
      do i = 0, nxmin1
          zzh(i,j,k,n) =  &
                  (zzh(i,j,k,n) + dumarray_pp(i,j,k))/9.81
      end do
      end do
      end do


! temperature - read perturbation potential temperature,
!     add 300. (base value), then add and convert
!   the wrf output file contains (nuvz-add_sfc_level) levels
!   read the data into k=kbgn,nuvz
      varname = 'T'
      lendim_exp(3) = nuvz-add_sfc_level
      lendim_max(3) = nuvzmax
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, tth(0,0,kbgn,n), &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of T', fnamenc
          stop
      end if

      do k = kbgn, nuvz
      do j = 0, nymin1
      do i = 0, nxmin1
! save potential temperature to ptth
        ptth(i,j,k,n)=tth(i,j,k,n)+300.                 
        tth(i,j,k,n) = (tth(i,j,k,n) + 300.) * &
                  (pph(i,j,k,n)/1.0e5)**0.286
  
      end do
      end do
      end do

!-
      if (turb_option .eq. turb_option_tke .or.  &
          turb_option .eq. turb_option_mytke ) then
      print*, 'READ TKE',turb_option
! TKE - read Turbulent Kinetic,
!   the wrf output file contains (nuvz-add_sfc_level) levels
!   read the data into k=kbgn,nuvz
      varname = 'TKE'
      lendim_exp(3) = nuvz-add_sfc_level
      lendim_max(3) = nuvzmax
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
          varname, tkeh(0,0,kbgn,n), &
          itime, &
          ndims, ndims_exp, ndims_max, &
          lendim, lendim_exp, lendim_max )

      if (ierr .ne. 0) then
      print*,'NO TKE available. Try TKE_PBL instead'
      varname = 'TKE_PBL'
      lendim_exp(3) = nuvz-add_sfc_level
      lendim_max(3) = nuvzmax
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
          varname, tkeh(0,0,kbgn,n), &
          itime, &
          ndims, ndims_exp, ndims_max, &
          lendim, lendim_exp, lendim_max )
      endif
      if (ierr .ne. 0) then
      print*,'NO TKE_PBL available. Try QKE instead'
      varname = 'qke'
      lendim_exp(3) = nuvz-add_sfc_level
      lendim_max(3) = nuvzmax
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
          varname, tkeh(0,0,kbgn,n), &
          itime, &
          ndims, ndims_exp, ndims_max, &
          lendim, lendim_exp, lendim_max )
       tkeh=tkeh/2. !conversion of qke
      endif
      if (ierr .ne. 0) then
      varname = 'QKE'
      lendim_exp(3) = nuvz-add_sfc_level
      lendim_max(3) = nuvzmax
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
          varname, tkeh(0,0,kbgn,n), &
          itime, &
          ndims, ndims_exp, ndims_max, &
          lendim, lendim_exp, lendim_max )
       tkeh=tkeh/2. !conversion of qke
      endif

      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of TKE', fnamenc
      write(*,*)'change turb_option NOT to use TKE or change input file'
          print*,'change SFC_OPTION to 0  as well'
          stop
      end if

       endif


!-

! specific humidity - read mixing ratio (kg-water-vapor/kg-dry-air),
!     then convert to (kg-water-vapor/kg-moist-air)
!   the wrf output file contains (nuvz-add_sfc_level) levels
!   read the data into k=kbgn,nuvz
      varname = 'QVAPOR'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, qvh(0,0,kbgn,n), &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of QVAPOR', fnamenc
          stop
      end if

      do k = kbgn, nuvz
      do j = 0, nymin1
      do i = 0, nxmin1
          qvh(i,j,k,n) = max( qvh(i,j,k,n), 0.0 )
          qvh(i,j,k,n) = qvh(i,j,k,n)/(1.0 + qvh(i,j,k,n))
      end do
      end do
      end do


! surface pressure
      varname = 'PSFC'
      lendim_exp(3) = 0
      lendim_max(3) = 1
      ndims_exp = 3
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, ps(0,0,1,n), &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of PSFC', fnamenc
          stop
      end if

! for the mexico city grid 3 simulation, the surface and
!   level 1 pressures are not as consistent as one would like,
!   with the problems occuring near the domain boundaries.
! so do the following
!   -- calculate surface pressure from lowest level pressure, temp, height
!   -- use wrf pressures (pph array) wherever possible
!      (avoid using surface pressure and the akz/bkz, akm/bkm)

!      do j = 0, nymin1
!      do i = 0, nxmin1
!          duma = ps(i,j,1,n)
!          dumdz = 0.5*(zzh(i,j,kbgn+1,n) - zzh(i,j,kbgn,n))
!          tv = tth(i,j,kbgn,n)*(1.+0.61*qvh(i,j,kbgn,n))
!          ps(i,j,1,n) = pph(i,j,kbgn,n)*exp( dumdz*ga/(r_air*tv) )
!      end do
!      end do


! 10 meter u velocity
!   note:  u10 is on the "T-grid" already
      varname = 'U10'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, u10(0,0,1,n), &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of U10', fnamenc
          stop
      end if


! 10 meter v velocity
!   note:  v10 is on the "T-grid" already
      varname = 'V10'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, v10(0,0,1,n), &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of V10', fnamenc
          stop
      end if


! 2 meter temperature
      varname = 'T2'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, tt2(0,0,1,n), &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of T2', fnamenc
          stop
      end if


! 2 meter dew point - read 2 meter water vapor mixing ratio
!   then calculate the dew point
      varname = 'Q2'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, td2(0,0,1,n), &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of Q2'
!, fnamenc
          do j = 0, nymin1
          do i = 0, nxmin1
! 29-nov-2005 - changed qvh(i,j,1,n) to qvh(i,j,kbgn,n) here
              td2(i,j,1,n) = qvh(i,j,kbgn,n)
          end do
          end do
      end if

      if (wind_option.ge.1) then
!      print*,'mean wind from WRF is used'
!      print*,'option ',wind_option
      varname = 'MU'

      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
          varname, mu(0,0,1), &
          itime, &
          ndims, ndims_exp, ndims_max, &
          lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing MU', fnamenc
          stop
      end if

      varname = 'MUB'

      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
          varname, mub(0,0,1), &
          itime, &
          ndims, ndims_exp, ndims_max, &
          lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing MUV', fnamenc
          stop
      end if
      endif

!      varname = 'MAPFAC_UX'
!      varname = 'MAPFAC_UY' !try
!      lendim_exp(1) = nx+1
!      lendim_max(1) = nxmax
!      lendim_exp(2) = ny
!      lendim_max(2) = nymax
!
!      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
!          varname, m_u(0,0,1), &
!          itime, &
!          ndims, ndims_exp, ndims_max, &
!          lendim, lendim_exp, lendim_max )
!      if (ierr .ne. 0) then
!          write(*,9100) 'error doing MAP U', fnamenc
!          stop
!      end if
!
!      varname = 'MAPFAC_VY'
!      varname = 'MAPFAC_VX' !try
!      lendim_exp(1) = nx
!      lendim_max(1) = nxmax
!      lendim_exp(2) = ny+1
!      lendim_max(2) = nymax
!
!      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
!          varname, m_v(0,0,1), &
!          itime, &
!          ndims, ndims_exp, ndims_max, &
!          lendim, lendim_exp, lendim_max )
!      if (ierr .ne. 0) then
!          write(*,9100) 'error doing MAP V', fnamenc
!          stop
!      end if

!      varname = 'MAPFAC_MX'
!      lendim_exp(1) = nx
!      lendim_max(1) = nxmax
!      lendim_exp(2) = ny
!      lendim_max(2) = nymax
!
!      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
!          varname, m_x(0,0,1), &
!          itime, &
!          ndims, ndims_exp, ndims_max, &
!          lendim, lendim_exp, lendim_max )
!      if (ierr .ne. 0) then
!          write(*,9100) 'error doing MAP X', fnamenc
!      varname = 'MAPFAC_U'
!      lendim_exp(1) = nx+1
!      lendim_max(1) = nxmax
!      lendim_exp(2) = ny
!      lendim_max(2) = nymax
!      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
!          varname, m_u(0,0,1), &
!          itime, &
!          ndims, ndims_exp, ndims_max, &
!          lendim, lendim_exp, lendim_max )
!      do j = 0, nymin1
!      do i = 0, nxmin1
!      m_x(i,j,1)=(m_u(i,j,1)+m_u(i+1,j,1))*0.5 
!      enddo
!      enddo
!      if (ierr .ne. 0) then
!          write(*,9100) 'error doing MAP U', fnamenc
!          print*,'NO MAP FACTOR IS GOING TO BE USED.'
!          print*,'LARGE UNCERTAINTIES TO BE EXPECTED'
!      do j = 0, nymin1
!      do i = 0, nxmin1
!      m_x(i,j,1)=1.
!      enddo
!      enddo
!      end if
!      end if

!      varname = 'MAPFAC_M'
!      lendim_exp(1) = nx
!      lendim_max(1) = nxmax
!      lendim_exp(2) = ny
!      lendim_max(2) = nymax
!
!      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
!          varname, m_z(0,0,1), &
!          itime, &
!          ndims, ndims_exp, ndims_max, &
!          lendim, lendim_exp, lendim_max )
!      if (ierr .ne. 0) then
!          write(*,9100) 'error doing MAP W', fnamenc
!          stop
!      end if

!      varname = 'MAPFAC_MY'
!      lendim_exp(1) = nx
!      lendim_max(1) = nxmax
!      lendim_exp(2) = ny
!      lendim_max(2) = nymax
!
!      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
!          varname, m_y(0,0,1), &
!          itime, &
!          ndims, ndims_exp, ndims_max, &
!          lendim, lendim_exp, lendim_max )
!      if (ierr .ne. 0) then
!          write(*,9100) 'error doing MAP Y', fnamenc
!      varname = 'MAPFAC_V'
!      lendim_exp(1) = nx
!      lendim_max(1) = nxmax
!      lendim_exp(2) = ny+1
!      lendim_max(2) = nymax
!      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
!          varname, m_v(0,0,1), &
!          itime, &
!          ndims, ndims_exp, ndims_max, &
!          lendim, lendim_exp, lendim_max )
!      do j = 0, nymin1
!      do i = 0, nxmin1
!      m_y(i,j,1)=(m_v(i,j,1)+m_v(i,j+1,1))*0.5
!      enddo
!      enddo
!      if (ierr .ne. 0) then
!          write(*,9100) 'ERROR doing MAP V', fnamenc
!          print*,'NO MAP FACTOR IS GOING TO BE USED.'
!          print*,'LARGE UNCERTAINTIES TO BE EXPECTED'
!      do j = 0, nymin1
!      do i = 0, nxmin1
!      m_y(i,j,1)=1.
!      enddo
!      enddo
!      end if
!      end if
      lendim_exp(1) = nx
      lendim_max(1) = nxmax
      lendim_exp(2) = ny
      lendim_max(2) = nymax



! calculate water vapor pressure in mb, from sfc pressure
!   and 2 m mixing ratio
      iduma = 0
      do j = 0, nymin1
      do i = 0, nxmin1
! 29-nov-2005 - added this to catch occasional tt2n=0.0 values
          duma = max( 100.0, tth(i,j,kbgn,n)-50.0 )
          if (tt2(i,j,1,n) .le. duma) then
              iduma = iduma + 1
              if (iduma .eq. 1) then
                  write(*,*) 'readwind - bad tt2 at'
                  write(*,*) 'i, j, tt2 =', i, j, tt2(i,j,1,n)
              end if
!             stop
              tt2(i,j,1,n) = tth(i,j,kbgn,n)
              td2(i,j,1,n) = qvh(i,j,kbgn,n)
          end if
          duma = td2(i,j,1,n)/0.622
          ewater_mb = 0.01*( 0.99976*ps(i,j,1,n)*duma/(1.0+duma) )
          esatwater_mb = 0.01*ew(tt2(i,j,1,n))
          ewater_mb = max( 1.0e-10, min( esatwater_mb, ewater_mb ) )
! then use the following, which is from an old 1970's report
!   (reference not available, but the formula works)
!   tdew(in C) = (4318.76/(19.5166 - ln(ewater(in mb)))) - 243.893
          td2(i,j,1,n) = 273.16 + &
                 (4318.76/(19.5166 - log(ewater_mb))) - 243.893
      end do
      end do
      if (iduma .gt. 0) write(*,*) &
          'readwind - bad tt2 count =', iduma


! sea level pressure - calculate it from surface pressure and 
!    ground elevation using standard atmosphere relations
      do j = 0, nymin1
      do i = 0, nxmin1
          msl(i,j,1,n) = ps(i,j,1,n)/ &
                  ((1.0 - 6.5e-3*oro(i,j)/288.0)**5.2553)
      end do
      end do


! large scale precipitation
! convective  precipitation
!   the wrf output files contain these as "accumulated totals"
!   I need to find out if these are accumulated over the output
!       file frequency, or over the total run.
!   For now, set to zero
! total cloud cover
!   Doesn't appear to be any 2-d cloud cover field in the
!       wrf output.
!   For now, set to zero
      do j = 0, nymin1
      do i = 0, nxmin1
          lsprec(i,j,1,n) = 0.0
          convprec(i,j,1,n) = 0.0
          tcc(i,j,1,n) = 0.0
      end do
      end do

!
! Large-scale precipitation (accumulated rain, mm)
!    will convert to mm/h  in interpolat_rain.f
      varname = 'RAINNC'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
          varname, lsprec(0,0,1,n), &
          itime, &
          ndims, ndims_exp, ndims_max, &
          lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
     write(*,9100) 'error doing ncread of RAINNC,set to zero', fnamenc
          do j = 0, nymin1
          do i = 0, nxmin1
              lsprec(i,j,1,n) = 0.0
          end do
          end do
      end if

!
! Convective cumulus precipitation (accumulated rain, mm)

      varname = 'RAINC'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
          varname, convprec(0,0,1,n), &
          itime, &
          ndims, ndims_exp, ndims_max, &
          lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
     write(*,9100) 'error doing ncread of RAINC, set to zero', fnamenc
          do j = 0, nymin1
          do i = 0, nxmin1
              convprec(i,j,1,n) = 0.0
          end do
          end do
      end if

!
! Clound fraction  (cloud cover)
!      varname = 'CLDFRA'
!      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
!          varname, tcc(0,0,1,n), &
!          itime, &
!          ndims, ndims_exp, ndims_max, &
!          lendim, lendim_exp, lendim_max )
!      if (ierr .ne. 0) then
!!     write(*,9100) 'error doing ncread of CLDFRA, set to zero'
!!, fnamenc
!          do j = 0, nymin1
!          do i = 0, nxmin1
!              tcc(i,j,1,n) = 0.0
!          end do
!          end do
!      end if
!!C        write(*,*)'read CLDFRA 0-sucess ',ierr

! land use, added 2015-03-27 //AD
      if (lu_option.eq.1) then  ! Read land-use from WRF?
        varname = 'LU_INDEX'
        call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
          varname, landuse_wrf(0,0), &
          itime, &
          ndims, ndims_exp, ndims_max, &
          lendim, lendim_exp, lendim_max )
        do k = 1, numclass  ! Re-initialize xlandusen
          do j = 0, nymin1
            do i = 0, nxmin1
              xlanduse(i,j,k) = 0.0
            enddo
          enddo
        enddo
        if (num_land_cat.eq.24) then ! Assume USGS landuse data
          do j = 0, nymin1
            do i = 0, nxmin1
              k = nint(landuse_wrf(i,j))  ! Safely convert element to integer (nearest)
              select case (k) ! Translate USGS categories to Wesely-types
                case(1)     ! USGS: Urban and built-up land 
                  xlanduse(i,j,1)  = 1.  ! Wesely: Urban land
                case(2:4)   ! USGS: Any cropland, pasture
                  xlanduse(i,j,2)  = 1.  ! Wesely: Agricultural land
                case(5)     ! USGS: Cropland/grassland mosaic
                  xlanduse(i,j,10) = 1.  ! Wesely: mixed aggricultural and range land
                case(6)     ! USGS: Cropland/woodland mosaic
                  xlanduse(i,j,4)  = 1.  ! Wesely: Deciduous Forest
                case(7)     ! USGS: Grassland
                  xlanduse(i,j,3)  = 1.  ! Wesely: Range land
                case(8)     ! USGS: Shrubland
                  xlanduse(i,j,11) = 1.  ! Wesely: rocky open areas with growing shrubs
                case(9)     ! USGS: Mixed Grassland and Shrubland
                  xlanduse(i,j,3)  = .5  ! Wesely: Range land
                  xlanduse(i,j,11) = .5  ! Wesely: rocky open areas with growing shrubs
                case(10)    ! USGS: Savanna
                  xlanduse(i,j,11) = 1.  ! Wesely: rocky open areas with growing shrubs
                case(11,12) ! USGS: Deciduous broadleaf (or needle leaf) forests
                  xlanduse(i,j,4)  = 1.  ! Wesely: Deciduous forest
                case(13)    ! USGS: Evergreen broadleaf forest
                  xlanduse(i,j,13) = 1.  ! Wesely: Rainforest
                case(14)    ! USGS: Evergreen needleleaf forest
                  xlanduse(i,j,5)  = 1.  ! Wesely: Coniferous Forest
                case(15)    ! USGS: Mixed forest
                  xlanduse(i,j,6)  = 1.  ! Wesely: Mixed forest including wetland
                case(16)    ! USGS: Water bodies
                  xlanduse(i,j,7)  = 1.  ! Wesely: water, both salt and fresh
                case(17,20) ! USGS: Herbaceous wetland or tundra
                  xlanduse(i,j,9)  = 1.  ! Wesely: nonforested wetland
                case(18,21) ! USGS: Wooded wetland or tundra
                  xlanduse(i,j,6)  = 1.  ! Wesely: Mixed Forest including Wetland
                case(19,23) ! USGS: "Barren or sparsely vegetated" or "Bare ground Tundra"
                  xlanduse(i,j,8)  = 1.  ! Wesely: Barren land mostly Desert
                case(22)    ! USGS: Mixed Tundra
                  xlanduse(i,j,6)  = .5  ! Wesely: Mixed Forest including Wetland
                  xlanduse(i,j,8)  = .5  ! Wesely: Barren land mostly Desert
                case(24)    ! USGS: Snow or Ice
                  xlanduse(i,j,12) = 1.  ! Wesely: snow or ice
                case default
                  write(*,9135) 'Unknown land use category',  &
                    'i,j,cat =', i,j,k, fnamenc
              endselect
            enddo
          enddo
        elseif (num_land_cat.eq.20) then ! Assume MODIS data
          do j = 0, nymin1
            do i = 0, nxmin1
              k = nint(landuse_wrf(i,j))  ! Safely convert element to integer (nearest)
              select case (k) ! Translate MODIS categories to Wesely-types
                case(1)   ! IGBP: Evergreen Needleleaf Forest
                  xlanduse(i,j,5)  = 1.  ! Wesely: Coniferous Forest
                case(2)   ! IGBP: Evergreen Broadleaf Forest
                  xlanduse(i,j,13) = 1.  ! Wesely: Rainforest
                case(3,4) ! IGBP: Deciduous Needleleaf or Broadleaf Forest
                  xlanduse(i,j,4)  = 1.  ! Wesely: Deciduous Forest
                case(5)   ! IGBP: Mixed Forest
                  xlanduse(i,j,6)  = 1.  ! Wesely: Mixed Forest including Wetland
                case(6:9) ! IGBP: Shrublands or Savannas
                  xlanduse(i,j,11) = 1.  ! Wesely: rocky open areas with growing shrubs
                case(10)  ! IGBP: Grasslands
                  xlanduse(i,j,3)  = 1.  ! Wesely: Range land
                case(11)  ! IGBP: Permanent Wetlands
                  xlanduse(i,j,9)  = 1.  ! Wesely: Non-forested wetland
                case(12)  ! IGBP: Croplands
                  xlanduse(i,j,2)  = 1.  ! Wesely: Agricultural land
                case(13)  ! IGBP: Urban and Built-up
                  xlanduse(i,j,1)  = 1.  ! Wesely: Urban land
                case(14)  ! IGBP: Cropland/Natural Vegetatation Mosaic
                  xlanduse(i,j,10)  = 1.  ! Wesely: mixed agricultural and range land
                case(15)  ! IGBP: Snow and Ice
                  xlanduse(i,j,12)  = 1.  ! Wesely: Snow and Ice
                case(16)  ! IGBP: Barren or Sparsely Vegetated
                  xlanduse(i,j,8)  = 1.  ! Wesely: Barren land mostly Desert
                case(17)  ! IGBP: Water Bodies
                  xlanduse(i,j,7)  = 1.  ! Wesely: water, both salt and fresh
                case(18)  ! IGBP: Wooded Tundra (WRF only)
                  xlanduse(i,j,6)  = 1.  ! Wesely: Mixed Forest including Wetland
                case(19)  ! IGBP: Mixed Tundra (WRF only)
                  xlanduse(i,j,6)  = .5  ! Wesely: Mixed Forest including Wetland
                  xlanduse(i,j,8)  = .5  ! Wesely: Barren land mostly Desert
                case(20)  ! IGBP: Barren Tundra (WRF only)
                  xlanduse(i,j,8)  = 1.  ! Wesely: Barren land mostly Desert
                case default
                  write(*,9135) 'Unknown land use category',  &
                    'i,j,cat =', i,j,k, fnamenc
              endselect
            enddo
          enddo
        endif
      endif ! lu_option.eq.1

! snow depth
      varname = 'SNOWH'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, sd(0,0,1,n), &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
!         write(*,9100) 'error doing ncread of SNOWH', fnamenc
          do j = 0, nymin1
          do i = 0, nxmin1
              sd(i,j,1,n) = 0.0
          end do
          end do
      end if


! surface sensible heat flux (positive <--> upwards)
      varname = 'HFX'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, sshf(0,0,1,n), &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
           do j = 0, nymin1
           do i = 0, nxmin1
               sshf(i,j,1,n) = -sshf(i,j,1,n)
           end do
           end do

      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of HFX', fnamenc
          do j = 0, nymin1
          do i = 0, nxmin1
              sshf(i,j,1,n) = 0.0
          end do
          end do
          hflswitch=.false.    ! Heat flux is not available
      else
          hflswitch=.true.     ! Heat flux is available
! limit to values to bounds originally used by flexpart?
!        do 1502 j=0,nymin1
!        do 1502 i=0,nxmin1
!           if(sshf(i,j,1,n).gt.200.) sshf(i,j,1,n)=200.
!           if(sshf(i,j,1,n).lt.-400.) sshf(i,j,1,n)=-400.
!1502    continue
      end if

! ustar
      varname = 'UST'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
            varname, ustar(0,0,1,n), &
            itime, &
            ndims, ndims_exp, ndims_max, &
            lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
         write(*,9100) 'error doing ncread of UST', fnamenc
         do j = 0, nymin1
         do i = 0, nxmin1
             ustar(i,j,1,n) = 0.0
         end do
         end do
         strswitch=.false.    ! ustar is not available
      else
         strswitch=.true.     ! ustar is available
         do j=0,nymin1
         do i=0,nxmin1
!           surfstr(i,j,1,n)=ustar(i,j,1,n)/dumarray_pp(i,j,kbgn)
         enddo
         enddo
      end if

! pblh
      if(sfc_option .eq. sfc_option_wrf) then
      varname = 'PBLH'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
          varname, hmix(0,0,1,n), &
          itime, &
          ndims, ndims_exp, ndims_max, &
          lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of PBLH', fnamenc
          stop
      endif

      endif

! surface solar radiation flux (positive <--> downwards)
      varname = 'SWDOWN'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, ssr(0,0,1,n), &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of SWDOWN', fnamenc
          do j = 0, nymin1
          do i = 0, nxmin1
              ssr(i,j,1,n) = 0.0
          end do
          end do
      else
          do j = 0, nymin1
          do i = 0, nxmin1
              ssr(i,j,1,n) = max( ssr(i,j,1,n), 0.0 )
          end do
          end do
      end if


! ew & ns surface stress
!   Doesn't appear to be any 2-d cloud cover field in the
!       wrf output.
!   For now, set to zero
      do j = 0, nymin1
      do i = 0, nxmin1
          ewss(i,j) = 0.0
          nsss(i,j) = 0.0
      end do
      end do
!     strswitch=.false.    ! Surface stress is not available


! orography
! standard deviation of orography
! land sea mask
!    these should be fixed during a simulation
!    so there is no reason to do them again ??


! *** done with reading the wrf output file ***

!  print*,'uu out1',uuh(0,259,1:10)
!  print*,'mu out1',mu(0,259,1),mub(0,259,1)
!  print*,'m_xn out1',m_x(0,259,1),m_y(0,259,1)


! interpolate uuh from the "U-grid" to the "T-grid"
! interpolate vvh from the "V-grid" to the "T-grid"
! new: convert mass weighted wind to wind.
!      print*,'wind_option',wind_option
      if (wind_option.le.0) then
!     if (wind_option.eq.-1) then
!     call calc_uvmet(uuh,vvh,urot,vrot,1)
!     endif
      do k = kbgn, nuvz
      do j = 0, nymin1
      do i = 0, nxmin1
! needs to rotate u and v to work. needs read alpha or something like that.
! if in mercator, no need.
!     divh(i,j,k)=(uuh(i+1,j,k)-uuh(i,j,k))/dx &
!      +(vvh(i,j+1,k)-vvh(i,j,k))/dy   
      if (wind_option.lt.0) then
      divh(i,j,k)=(uuh(i+1,j,k)-uuh(i,j,k))/dx*m_x(i,j,1) &
       +(vvh(i,j+1,k)-vvh(i,j,k))/dy*m_y(i,j,1)   
      endif 
          uuh(i,j,k) = 0.5*(uuh(i,j,k) + uuh(i+1,j,k))
          vvh(i,j,k) = 0.5*(vvh(i,j,k) + vvh(i,j+1,k))
      end do
      end do
      end do
      elseif (wind_option.eq.1) then
      do k = kbgn, nuvz
      do j = 0, nymin1
      do i = 0, nxmin1
!         uuh(i,j,k) = 0.5*(uuh(i,j,k)*m_u(i,j,1) + uuh(i+1,j,k)*m_u(i+1,j,1))
!         vvh(i,j,k) = 0.5*(vvh(i,j,k)*m_v(i,j,1) + vvh(i,j+1,k)*m_v(i,j+1,1))
           uuh(i,j,k) = 0.5*(uuh(i,j,k) + uuh(i+1,j,k))
           vvh(i,j,k) = 0.5*(vvh(i,j,k) + vvh(i,j+1,k))
      mu2=mu(i,j,1)+mub(i,j,1)
!      uuh(i,j,k) = uuh(i,j,k)/mu2!*0.5*(m_u(i,j,1)+m_u(i+1,j,1))
!      vvh(i,j,k) = vvh(i,j,k)/mu2!*0.5*(m_v(i,j,1)+m_v(i,j+1,1))
! m_y=0.5*(m_v(i,j,1)+m_v(i,j+1,1))
      uuh(i,j,k) = uuh(i,j,k)/mu2 !*m_y(i,j,1) !without m=true wind
      vvh(i,j,k) = vvh(i,j,k)/mu2 !*m_x(i,j,1)
      wwh(i,j,k) = wwh(i,j,k)/mu2 !*m_y(i,j,1) 
      end do
      end do
      end do
      elseif (wind_option.eq.2) then
      do k = kbgn, nuvz
      do j = 0, nymin1
      do i = 0, nxmin1
          uuh(i,j,k) = 0.5*(uuh(i,j,k) + uuh(i+1,j,k))
          vvh(i,j,k) = 0.5*(vvh(i,j,k) + vvh(i,j+1,k))
      mu2=mu(i,j,1)+mub(i,j,1)
      wwh(i,j,k) = wwh(i,j,k)/mu2 !*m_y(i,j,1)
      end do
      end do
      end do

      endif


! for ecmwf flexpart, if nwz = nlev_ec+1, then wwh is set
!   to zero at the top level
! for wrf, nlev_ec==n_bottom_top and so nwz = nlev_ec+1.
!   however, it doesn't seem appropriate to zero wwh at 
!   the model top which might be ~100 hPa.
! so deactivate this for now
!      if(levdiff2.eq.0) then
!        iwmax=nlev_ec+1
!        do 60 i=0,nxmin1
!        do 60 j=0,nymin1
!60      wwh(i,j,nlev_ec+1)=0.
!      endif

! For global fields, assign the leftmost data column also to the rightmost
! data column; if required, shift whole grid by nxshift grid points
!
! FLEXPART_WRF - all "global" stuff is turned off
!*************************************************************************

!     if (xglobal) then
!       call shift_field_0(ewss,nxfield,ny)
!       call shift_field_0(nsss,nxfield,ny)
!       call shift_field_0(oro,nxfield,ny)
!       call shift_field_0(excessoro,nxfield,ny)
!       call shift_field_0(lsm,nxfield,ny)
!       call shift_field(ps,nxfield,ny,1,1,2,n)
!       call shift_field(sd,nxfield,ny,1,1,2,n)
!       call shift_field(msl,nxfield,ny,1,1,2,n)
!       call shift_field(tcc,nxfield,ny,1,1,2,n)
!       call shift_field(u10,nxfield,ny,1,1,2,n)
!       call shift_field(v10,nxfield,ny,1,1,2,n)
!       call shift_field(tt2,nxfield,ny,1,1,2,n)
!       call shift_field(td2,nxfield,ny,1,1,2,n)
!       call shift_field(lsprec,nxfield,ny,1,1,2,n)
!       call shift_field(convprec,nxfield,ny,1,1,2,n)
!       call shift_field(sshf,nxfield,ny,1,1,2,n)
!       call shift_field(ssr,nxfield,ny,1,1,2,n)
!       call shift_field(tth,nxfield,ny,nuvzmax,nuvz,2,n)
!       call shift_field(qvh,nxfield,ny,nuvzmax,nuvz,2,n)
!       call shift_field(uuh,nxfield,ny,nuvzmax,nuvz,1,1)
!       call shift_field(vvh,nxfield,ny,nuvzmax,nuvz,1,1)
!       call shift_field(wwh,nxfield,ny,nwzmax,nwz,1,1)
!     endif

! CALCULATE SURFSTR
       if(sfc_option .eq. sfc_option_diagnosed) then
         do  i=0,nxmin1
         do  j=0,nymin1
           surfstr(i,j,1,n)=sqrt(ewss(i,j)**2+nsss(i,j)**2)
         enddo
         enddo
         strswitch=.false.
       endif

      if ((.not.hflswitch).or.(.not.strswitch)) then
        write(*,*) 'WARNING: No (or incomplete) flux data ' //  &
        'contained in WRF output file ', &
        wfname(indj)
 
! CALCULATE USTAR AND SSHF USING THE PROFILE METHOD
!    As ECMWF has increased the model resolution, such that now the first model
!    level is at about 10 m (where 10-m wind is given), use the 2nd ECMWF level
!    (3rd model level in FLEXPART) for the profile method 
!
! FLEXPART_WRF - use k=(2+add_sfc_level) here instead of k=3
!***************************************************************************

        k = 2 + add_sfc_level
        do i=0,nxmin1
          do j=0,nymin1
!           plev1=akz(3)+bkz(3)*ps(i,j,1,n)
            plev1=pph(i,j,k,n)
            pmean=0.5*(ps(i,j,1,n)+plev1)
            tv=tth(i,j,k,n)*(1.+0.61*qvh(i,j,k,n))
            fu=-r_air*tv/ga/pmean
            hlev1=fu*(plev1-ps(i,j,1,n))   ! HEIGTH OF FIRST MODEL LAYER
            ff10m= sqrt(u10(i,j,1,n)**2+v10(i,j,1,n)**2)
            fflev1=sqrt(uuh(i,j,k)**2+vvh(i,j,k)**2)
            call pbl_profile(ps(i,j,1,n),td2(i,j,1,n),hlev1, &
                             tt2(i,j,1,n),tth(i,j,k,n),ff10m,fflev1, &
                             surfstr(i,j,1,n),sshf(i,j,1,n))
            if(sshf(i,j,1,n).gt.200.) sshf(i,j,1,n)=200.
            if(sshf(i,j,1,n).lt.-400.) sshf(i,j,1,n)=-400.
           enddo
           enddo
      endif


! Assign 10 m wind to model level at eta=1.0 to have one additional model
!     level at the ground
! Specific humidity is taken the same as at one level above
! Temperature is taken as 2 m temperature         
!
! Note that the uuh, vvh, tth, & qvh data have already been shifted
!     upwards by one level, when they were read in.
!**************************************************************************

      if (add_sfc_level .eq. 1) then
      do j = 0, nymin1
      do i = 0, nxmin1
          uuh(i,j,1)   = u10(i,j,1,n)
          vvh(i,j,1)   = v10(i,j,1,n)
          tth(i,j,1,n) = tt2(i,j,1,n)
         ptth(i,j,1,n) = ptth(i,j,2,n)
          qvh(i,j,1,n) = qvh(i,j,2,n)
          tkeh(i,j,1,n)=tkeh(i,j,2,n)
! pressure at 2 m AGL
          pph(i,j,1,n) = 0.99976*ps(i,j,1,n)
! height (MSL) at ground level (shift it down)
          zzh(i,j,1,n) = zzh(i,j,2,n)
! height (MSL) at top of the added level
          zzh(i,j,2,n) = zzh(i,j,1,n) + 4.0
      if (hmix(i,j,1,n).lt.hmixmin) hmix(i,j,1,n)=hmixmin

      enddo
      enddo
      end if


       do i=0,nxmax-1
        do j=0,nymax-1
         do k=1,nuvzmax
           u_wrf(i,j,k,n)=uuh(i,j,k)
           v_wrf(i,j,k,n)=vvh(i,j,k)
         enddo
        enddo
       enddo
 
       do i=0,nxmax-1
        do j=0,nymax-1
         do k=1,nwzmax
           w_wrf(i,j,k,n)=wwh(i,j,k)
         enddo
        enddo
       enddo
 




      return    
      end subroutine readwind

