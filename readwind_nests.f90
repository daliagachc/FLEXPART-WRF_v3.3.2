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
      subroutine readwind_nests(indj,n,uuhn,vvhn,wwhn,divhn)
!                                i   i  o    o    o,   o
!*******************************************************************************
!                                                                              *
!     Note:  This is the FLEXPART_WRF version of subroutine readwind_nests.    *
!            The met fields are read from WRF netcdf output files.             *
!            There are many differences from the FLEXPART version.             *
!                                                                              *
!     This routine reads the wind fields for the nested model domains.         *
!     It is similar to subroutine readwind, which reads the mother domain.     *
!                                                                              *
!     Authors: A. Stohl, G. Wotawa                                             *
!                                                                              *
!     8 February 1999                                                          *
!     Last update: 2015-03-26, A. Dingwell                                     *
!                                                                              *
!                                                                              *
!  Feb. 2001, Bernd C. Krueger:                                                *
!             Variables tthn and qvhn (on eta coordinates) in common block     *
!                                                                              *
!  Oct-Dec 2005, R. Easter:                                                    *
!             Major changes for WRF.                                           *
!                                                                              *
!  11 June  2007, input tkehn from WRF                                         *
!  13 JUNE  2007  add ext_scalar, pbl_physics                                  *
!  19 Oct   2007  add RAINC, RAINNC, CLDFRA                                    *
!  Feb 2012, adapt it for mean wind. Jerome Brioude                            *
!                                                                              * 
!  2015-03-26, A. Dingwell:                                                    *
!             Updated calls to read_ncwrfout_gridinfo to match updates in that *
!             subroutine.                                                      * 
!                                                                              * 
!*******************************************************************************

  use par_mod
  use com_mod
  use flexwrf_ncdf_mod

!      include 'includepar'
!      include 'includecom'

! subr arguments
      integer :: indj,n
      real(kind=4) :: uuhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
      real(kind=4) :: vvhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
      real(kind=4) :: wwhn(0:nxmaxn-1,0:nymaxn-1,nwzmax,maxnests)
      real(kind=4) :: divhn(0:nxmaxn-1,0:nymaxn-1,nwzmax,maxnests)
      real(kind=4) ::  mu(0:nxmaxn-1,0:nymaxn-1,1),mu2
      real(kind=4) :: mub(0:nxmaxn-1,0:nymaxn-1,1)
!      real(kind=4) :: m_u(0:nxmaxn-1,0:nymaxn-1,1)
!      real(kind=4) :: m_v(0:nxmaxn-1,0:nymaxn-1,1)
!      real(kind=4) :: m_w(0:nxmaxn-1,0:nymaxn-1,1)
!      real :: uuhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
!      real :: vvhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
!      real :: wwhn(0:nxmaxn-1,0:nymaxn-1,nwzmax,maxnests)
!      real :: mu(0:nxmaxn-1,0:nymaxn-1,1),mu2
!      real :: mub(0:nxmaxn-1,0:nymaxn-1,1)
!      real(kind=4) :: m_un(0:nxmaxn-1,0:nymaxn-1,1,maxnests)
!      real(kind=4) :: m_vn(0:nxmaxn-1,0:nymaxn-1,1,maxnests)
!      real :: m_w(0:nxmaxn-1,0:nymaxn-1,1)


! local variables
      integer,parameter :: ndims_max=4

      integer :: i, idiagaa, ierr, itime, iclass
      integer :: iduma
      integer :: j, jhhmmss, jyyyymmdd
      integer :: k, kbgn
      integer :: l
      integer :: lendim(ndims_max), lendim_exp(ndims_max), &
          lendim_max(ndims_max)
      integer :: m,levdiff2
      integer :: ndims, ndims_exp, &
              ext_scalar,pbl_physics,mp_physics_dum, num_land_cat
      integer :: n_west_east, n_south_north, n_bottom_top
      integer :: m_grid_id_dum, m_parent_grid_id_dum, &
        m_parent_grid_ratio_dum,  &
        i_parent_start_dum, j_parent_start_dum, &
        map_proj_id_dum

      real :: dx_met, dy_met
      real :: duma, dumb, dumc, dumd, dume
      real :: dumdz
      real :: dumarray_aa(nwzmax+1)
      real(kind=4) :: dumarray_pp(0:nxmaxn-1,0:nymaxn-1,nwzmax+1)
      real(kind=4) :: landuse_wrf(0:nxmaxn-1,0:nymaxn-1) !land use cat. from WRF
      real :: ewater_mb, esatwater_mb
      real :: ew      ! this is an external function
      real :: map_stdlon_dum, map_truelat1_dum, map_truelat2_dum
      real :: toler

      real :: ewss(0:nxmaxn-1,0:nymaxn-1),nsss(0:nxmaxn-1,0:nymaxn-1)
      real :: plev1,pmean,tv,fu,hlev1,ff10m,fflev1

      real(kind=dp) :: jul,juldate
      character(len=160) :: fnamenc, varname,fnamenc2

      logical :: hflswitch


!
! main loop -- process each nest
!
      do l=1,numbnests

!
!   get grid info from the wrf netcdf file
!   and check it for consistency against values from gridcheck
!
      m = numpath+2*(l-1)+1
      fnamenc = path(m)(1:length(m)) // wfnamen(l,indj)

      idiagaa = 0

      call read_ncwrfout_gridinfo( ierr, idiagaa, fnamenc, &
        n_west_east, n_south_north, n_bottom_top,  &
        dx_met, dy_met,  &
        m_grid_id_dum, m_parent_grid_id_dum, m_parent_grid_ratio_dum,  &
        i_parent_start_dum, j_parent_start_dum, &
        map_proj_id_dum, map_stdlon_dum,  &
        map_truelat1_dum, map_truelat2_dum, &
        ext_scalar,pbl_physics,mp_physics_dum, num_land_cat)
      if (ierr .ne. 0) then
          write(*,9100) l, 'error getting gridinfo for met file',  &
              fnamenc
          stop
      end if

! subtract 1 here because i & j indexing in flexpart always starts at 0
      i_parent_start_dum = i_parent_start_dum-1
      j_parent_start_dum = j_parent_start_dum-1

9100	format( / '*** readwind_nests, l=', i2, ' -- ',  &
            a / 'file = ', a )
9110	format( / '*** readwind_nests, l=', i2, ' -- ',  &
            a, 1x, i8 / 'file = ', a )
9120	format( / '*** readwind_nests, l=', i2, ' -- ',  &
            a, 2(1x,i8) / 'file = ', a )
9130	format( / '*** readwind_nests, l=', i2, ' -- ',  &
            a, 3(1x,i8) / 'file = ', a )
9115	format( / '*** readwind_nests, l=', i2, ' -- ',  &
            a / a, 1x, i8 / 'file = ', a )
9125	format( / '*** readwind_nests, l=', i2, ' -- ',  &
            a / a, 2(1x,i8) / 'file = ', a )
9135	format( / '*** readwind_nests, l=', i2, ' -- ',  &
            a / a, 3(1x,i8) / 'file = ', a )

      toler = 2.0e-7

      if (nxn(l) .ne. n_west_east) then
          write(*,9100) l, 'nx not consistent', fnamenc
          stop
      end if
      if (nyn(l) .ne. n_south_north) then
          write(*,9100) l, 'ny not consistent', fnamenc
          stop
      end if
      if (nlev_ec .ne. n_bottom_top) then
          write(*,9100) l, 'nlev_ec not consistent', fnamenc
          stop
      end if
      if (nwz .ne. n_bottom_top+1) then
          write(*,9100) l, 'nwz not consistent', fnamenc
          stop
      end if
      if (nuvz .ne. n_bottom_top+add_sfc_level) then
          write(*,9100) l, 'nuvz not consistent', fnamenc
          stop
      end if

      if (m_grid_id(l) .ne. m_grid_id_dum) then
          write(*,9100) l, 'm_grid_id not consistent', fnamenc
          write(*,*) m_grid_id(l), m_grid_id_dum
          stop
      end if
      if (m_parent_grid_id(l) .ne. m_parent_grid_id_dum) then
          write(*,9100) l, 'm_parent_grid_id not consistent', fnamenc
          stop
      end if
      if (m_parent_grid_ratio(l) .ne. m_parent_grid_ratio_dum) then
          write(*,9100) l, 'm_parent_grid_ratio not consistent', fnamenc
          stop
      end if
      if (i_parent_start(l) .ne. i_parent_start_dum) then
          write(*,9100) l, 'i_parent_start not consistent', fnamenc
          stop
      end if
      if (j_parent_start(l) .ne. j_parent_start_dum) then
          write(*,9100) l, 'j_parent_start not consistent', fnamenc
          stop
      end if

      if (abs(dxn(l) - dx_met) .gt. toler*abs(dxn(l))) then
          write(*,9100) l, 'dx not consistent', fnamenc
          stop
      end if
      if (abs(dyn(l) - dy_met) .gt. toler*abs(dyn(l))) then
          write(*,9100) l, 'dy not consistent', fnamenc
          stop
      end if

! locate the date/time in the file
      itime = 0
1100  itime = itime + 1
      call read_ncwrfout_1datetime( ierr, fnamenc, &
          itime, jyyyymmdd, jhhmmss )
      if (ierr .eq. -1) then
          write(*,9100) l, 'error reading time from met file', fnamenc
          stop
      else if (ierr .ne. 0) then
          write(*,9125) l, 'unable to locate date/time in met file',  &
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
      write(*,*) 'readwind_nests processing wrfout file ='
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
      lendim_max(1) = nwzmax
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
          write(*,9100) l, 'error doing ncread of ZNW', fnamenc
          stop
      end if
      end if
      do k = 1, nwz
          if (abs(eta_w_wrf(k) - dumarray_aa(k))  &
                  .gt. toler*abs(eta_w_wrf(k))) then
              write(*,9100) l, 'eta_w_wrf not consistent', fnamenc
              stop
          end if
      end do

      varname = 'ZNU'
      lendim_exp(1) = nwz-1
      lendim_max(1) = nwzmax
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

          write(*,9100) l, 'error doing ncread of ZNU', fnamenc
          stop
      end if
      end if
      do k = 1, nwz-1
          if (abs(eta_u_wrf(k) - dumarray_aa(k))  &
                  .gt. toler*abs(eta_u_wrf(k))) then
              write(*,9100) l, 'eta_u_wrf not consistent', fnamenc
              stop
          end if
      end do

!      varname = 'P_TOP'
!      lendim_exp(1) = 1
!      lendim_max(1) = 1
!      ndims_exp = 2
!      if (ext_scalar .lt. 0) ndims_exp = 1
!      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
!      	  varname, duma, &
!      	  itime, &
!      	  ndims, ndims_exp, ndims_max, &
!      	  lendim, lendim_exp, lendim_max )
!      if (ierr .ne. 0) then
!          write(*,9100) l, 'error doing ncread of P_TOP', fnamenc
!          stop
!      end if
!      if (abs(p_top_wrf - duma) .gt. toler*abs(p_top_wrf)) then
!          write(*,9100) l, 'p_top_wrf not consistent', fnamenc
!          stop
!      end if

      varname = 'XLAT'
      lendim_exp(1) = nxn(l)
      lendim_max(1) = nxmaxn
      lendim_exp(2) = nyn(l)
      lendim_max(2) = nymaxn
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
      do j = 0, nyn(l)-1
      do i = 0, nxn(l)-1
          if (abs(ylat2dn(i,j,l) - dumarray_pp(i,j,1)) .gt.  &
                              toler*abs(ylat2dn(i,j,l))) then
              write(*,9100) l, 'ylat2dn not consistent', fnamenc
              write(*,'(a,2i5,2f16.6)') 'i,j,ylats =', i, j, &
                      ylat2dn(i,j,l), dumarray_pp(i,j,1)
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
      do j = 0, nyn(l)-1
      do i = 0, nxn(l)-1
          if (abs(xlon2dn(i,j,l) - dumarray_pp(i,j,1)) .gt.  &
                              toler*abs(xlon2dn(i,j,l))) then
              write(*,9100) l, 'xlon2dn not consistent', fnamenc
              write(*,'(a,2i5,2f16.6)') 'i,j,xlons =', i, j, &
                      xlon2dn(i,j,l), dumarray_pp(i,j,1)
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
! JB
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
      lendim_exp(1) = nxn(l)+1
      lendim_max(1) = nxmaxn
      lendim_exp(2) = nyn(l)
      lendim_max(2) = nymaxn
      lendim_exp(3) = nuvz-add_sfc_level
      lendim_max(3) = nuvzmax
      ndims_exp = 4
      if (time_option.eq.0) then
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, uuhn(0,0,kbgn,l), &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) l, 'error doing ncread of U', fnamenc
      if (wind_option.le.0) print*,'you asked snapshot winds'
      if (wind_option.eq.1) print*,'you asked mean winds'
        print*,'change wind_option'
          stop
      end if
      endif

! v wind velocity
!   the wrf output file contains (nuvz-add_sfc_level) levels
!   read the data into k=kbgn,nuvz
!   (interpolate it from "V-grid" to "T-grid" later)
      if (wind_option.le.0) varname = 'V'
      if (wind_option.eq.1) varname = 'AVGFLX_RVM'
      if (wind_option.eq.2) varname = 'V'

      lendim_exp(1) = nxn(l)
      lendim_max(1) = nxmaxn
      lendim_exp(2) = nyn(l)+1
      lendim_max(2) = nymaxn
      if (time_option.eq.0) then
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, vvhn(0,0,kbgn,l), &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) l, 'error doing ncread of V', fnamenc
      if (wind_option.le.0) print*,'you asked snapshot winds'
      if (wind_option.eq.1) print*,'you asked mean winds'
        print*,'change wind_option'
          stop
      end if

      endif
! w wind velocity
!   this is on the "W-grid", and 
!   the wrf output file contains nwz levels, so no shifting needed
!     varname = 'W'
      if (wind_option.le.0) varname = 'W'
      if (wind_option.eq.1) varname = 'AVGFLX_WWM'
      if (wind_option.eq.2) varname = 'WW'

      lendim_exp(1) = nxn(l)
      lendim_max(1) = nxmaxn
      lendim_exp(2) = nyn(l)
      lendim_max(2) = nymaxn
      lendim_exp(3) = nwz
      lendim_max(3) = nwzmax
      if (time_option.eq.0) then
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, wwhn(0,0,1,l), &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) l, 'error doing ncread of W', fnamenc
      if (wind_option.eq.0) print*,'you asked snapshot winds'
      if (wind_option.eq.1) print*,'you asked mean winds'
        print*,'change wind_option'

          stop
      end if
      endif

! pressure - read base state and perturbation pressure,
!     then combine
!   the wrf output file contains (nuvz-add_sfc_level) levels
!   read the data into k=kbgn,nuvz
      varname = 'PB'
      lendim_exp(3) = nuvz-1
      lendim_max(3) = nuvzmax
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, pphn(0,0,kbgn,n,l), &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) l, 'error doing ncread of PB', fnamenc
          stop
      end if

      varname = 'P'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, dumarray_pp(0,0,kbgn), &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) l, 'error doing ncread of P', fnamenc
          stop
      end if

      do k = kbgn, nuvz
      do j = 0, nyn(l)-1
      do i = 0, nxn(l)-1
          pphn(i,j,k,n,l) = pphn(i,j,k,n,l) + dumarray_pp(i,j,k)
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
      	  varname, zzhn(0,0,kbgn,n,l), &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) l, 'error doing ncread of PB', fnamenc
          stop
      end if

      varname = 'PH'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, dumarray_pp(0,0,kbgn), &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) l, 'error doing ncread of P', fnamenc
          stop
      end if

      do k = kbgn, nwz+add_sfc_level
      do j = 0, nyn(l)-1
      do i = 0, nxn(l)-1
          zzhn(i,j,k,n,l) =  &
                  (zzhn(i,j,k,n,l) + dumarray_pp(i,j,k))/9.81
      end do
      end do
      end do


! temperature - read perturbation potential temperature,
!     add t00 (base value), then add and convert
!   the wrf output file contains (nuvz-add_sfc_level) levels
!   read the data into k=kbgn,nuvz
      varname = 'T'
      lendim_exp(3) = nuvz-1
      lendim_max(3) = nuvzmax
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, tthn(0,0,kbgn,n,l), &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) l, 'error doing ncread of T', fnamenc
          stop
      end if

      do k = kbgn, nuvz
      do j = 0, nyn(l)-1
      do i = 0, nxn(l)-1
! save pot tempereature to ptthn
         ptthn(i,j,k,n,l) =  tthn(i,j,k,n,l)+300.
          tthn(i,j,k,n,l) = (tthn(i,j,k,n,l) + 300.) * &
                  (pphn(i,j,k,n,l)/1.0e5)**0.286
      end do
      end do
      end do

      if (turb_option .eq. turb_option_tke .or. &
          turb_option .eq. turb_option_mytke) then
!-
! TKE - read TKE
!   the wrf output file contains (nuvz-add_sfc_level) levels
!   read the data into k=kbgn,nuvz
      varname = 'TKE'
      lendim_exp(3) = nuvz-1
      lendim_max(3) = nuvzmax
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
          varname, tkehn(0,0,kbgn,n,l), &
          itime, &
          ndims, ndims_exp, ndims_max, &
          lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
      varname = 'TKE_PBL'
      lendim_exp(3) = nuvz-1
      lendim_max(3) = nuvzmax
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
          varname, tkehn(0,0,kbgn,n,l), &
          itime, &
          ndims, ndims_exp, ndims_max, &
          lendim, lendim_exp, lendim_max )
      endif
      if (ierr .ne. 0) then
!     print*,'NO TKE_PBL available. Try QKE instead'
      varname = 'qke'
      lendim_exp(3) = nuvz-1
      lendim_max(3) = nuvzmax
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
          varname, tkehn(0,0,kbgn,n,l), &
          itime, &
          ndims, ndims_exp, ndims_max, &
          lendim, lendim_exp, lendim_max )
      tkeh=tkeh/2. !conversion of qke
      endif
      if (ierr .ne. 0) then
          write(*,9100) l, 'error doing ncread of TKE', fnamenc
      write(*,*)'Change turb_option NOT to use TKE, or change inputfile'
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
      	  varname, qvhn(0,0,kbgn,n,l), &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) l, 'error doing ncread of QVAPOR', fnamenc
          stop
      end if

      do k = kbgn, nuvz
      do j = 0, nyn(l)-1
      do i = 0, nxn(l)-1
          qvhn(i,j,k,n,l) = max( qvhn(i,j,k,n,l), 0.0 )
          qvhn(i,j,k,n,l) = qvhn(i,j,k,n,l)/(1.0 + qvhn(i,j,k,n,l))
      end do
      end do
      end do


! surface pressure
      varname = 'PSFC'
      lendim_exp(3) = 0
      lendim_max(3) = 1
      ndims_exp = 3
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, psn(0,0,1,n,l), &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) l, 'error doing ncread of PSFC', fnamenc
          stop
      end if

! for the mexico city grid 3 simulation, the surface and
!   level 1 pressures are not as consistent as one would like,
!   with the problems occuring near the domain boundaries.
! so diagnose surface pressure from other variables

!      do j = 0, nyn(l)-1
!      do i = 0, nxn(l)-1
!
!! better fix 
!!   -- calculate surface pressure from lowest level pressure, temp, height
!!   -- use wrf pressures (pph array) wherever possible
!!      (avoid using surface pressure and the akz/bkz, akm/bkm)
!          duma = psn(i,j,1,n,l)
!          dumdz = 0.5*(zzhn(i,j,kbgn+1,n,l) - zzhn(i,j,kbgn,n,l))
!          tv = tthn(i,j,kbgn,n,l)*(1.+0.61*qvhn(i,j,kbgn,n,l))
!          psn(i,j,1,n,l) = pphn(i,j,kbgn,n,l)*exp( dumdz*ga/(r_air*tv) )
!
!      end do
!      end do


! 10 meter u velocity
!   note:  u10 is on the "T-grid" already
      varname = 'U10'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, u10n(0,0,1,n,l), &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) l, 'error doing ncread of U10', fnamenc
          stop
      end if


! 10 meter v velocity
!   note:  v10 is on the "T-grid" already
      varname = 'V10'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, v10n(0,0,1,n,l), &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) l, 'error doing ncread of V10', fnamenc
          stop
      end if


! 2 meter temperature
      varname = 'T2'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, tt2n(0,0,1,n,l), &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) l, 'error doing ncread of T2', fnamenc
          stop
      end if


! 2 meter dew point - read 2 meter water vapor mixing ratio
!   then calculate the dew point
      varname = 'Q2'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, td2n(0,0,1,n,l), &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) l, 'error doing ncread of Q2', fnamenc
          do j = 0, nyn(l)-1
          do i = 0, nxn(l)-1
! 29-nov-2005 - changed qvhn(i,j,1,n,l) to qvhn(i,j,kbgn,n,l) here
              td2n(i,j,1,n,l) = qvhn(i,j,kbgn,n,l)
          end do
          end do
      end if

      if (wind_option.ge.1) then
!      print*,'mean wind from WRF is used'
      varname = 'MU '

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
          write(*,9100) 'error doing MUB', fnamenc
          stop
      end if
       endif

!      varname = 'MAPFAC_MX'
!      lendim_exp(1) = nxn(l)
!      lendim_max(1) = nxmaxn
!      lendim_exp(2) = nyn(l)
!      lendim_max(2) = nymaxn
!
!      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
!          varname, m_xn(0,0,1,l), &
!          itime, &
!          ndims, ndims_exp, ndims_max, &
!          lendim, lendim_exp, lendim_max )
!      if (ierr .ne. 0) then
!          write(*,9100) 'error doing MAP X', fnamenc
!      varname = 'MAPFAC_U'
!      lendim_exp(1) = nxn(l)+1
!      lendim_max(1) = nxmaxn
!      lendim_exp(2) = nyn(l)
!      lendim_max(2) = nymaxn
!      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
!          varname, m_un(0,0,1,l), &
!          itime, &
!          ndims, ndims_exp, ndims_max, &
!          lendim, lendim_exp, lendim_max )
!      do j = 0, nyn(l)-1
!      do i = 0, nxn(l)-1
!      m_xn(i,j,1,l)=(m_un(i,j,1,l)+m_un(i+1,j,1,l))*0.5
!      enddo
!      enddo
!      if (ierr .ne. 0) then
!          write(*,9100) 'error doing MAP U', fnamenc
!          print*,'NO MAP FACTOR IS GOING TO BE USED.'
!          print*,'LARGE UNCERTAINTIES TO BE EXPECTED'
!      do j = 0, nyn(l)-1
!      do i = 0, nxn(l)-1
!      m_xn(i,j,1)=1.
!      enddo
!      enddo
!      end if
!      end if
!
!      varname = 'MAPFAC_MY'
!      lendim_exp(1) = nxn(l)
!      lendim_max(1) = nxmaxn
!      lendim_exp(2) = nyn(l)
!      lendim_max(2) = nymaxn
!
!      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
!          varname, m_yn(0,0,1,l), &
!          itime, &
!          ndims, ndims_exp, ndims_max, &
!          lendim, lendim_exp, lendim_max )
!      if (ierr .ne. 0) then
!          write(*,9100) 'error doing MAP Y', fnamenc
!      varname = 'MAPFAC_V'
!      lendim_exp(1) = nxn(l)
!      lendim_max(1) = nxmaxn
!      lendim_exp(2) = nyn(l)+1
!      lendim_max(2) = nymaxn
!      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
!          varname, m_vn(0,0,1,l), &
!          itime, &
!          ndims, ndims_exp, ndims_max, &
!          lendim, lendim_exp, lendim_max )
!      do j = 0, nyn(l)-1
!      do i = 0, nxn(l)-1
!      m_yn(i,j,1,l)=(m_vn(i,j,1,l)+m_vn(i,j+1,1,l))*0.5
!      enddo
!      enddo
!      if (ierr .ne. 0) then
!          write(*,9100) 'ERROR doing MAP V', fnamenc
!          print*,'NO MAP FACTOR IS GOING TO BE USED.'
!          print*,'LARGE UNCERTAINTIES TO BE EXPECTED'
!      do j = 0, nyn(l)-1
!      do i = 0, nxn(l)-1
!      m_yn(i,j,1,l)=1.
!      enddo
!      enddo
!      end if
!      end if
      lendim_exp(1) = nxn(l)
      lendim_max(1) = nxmaxn
      lendim_exp(2) = nyn(l)
      lendim_max(2) = nymaxn

!      varname = 'MAPFAC_U'
!      lendim_exp(1) = nxn(l)+1
!      lendim_max(1) = nxmaxn
!      lendim_exp(2) = nyn(l)
!      lendim_max(2) = nymaxn
!
!      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
!          varname, m_un(0,0,1,l), &
!          itime, &
!          ndims, ndims_exp, ndims_max, &
!          lendim, lendim_exp, lendim_max )
!      if (ierr .ne. 0) then 
!          write(*,9100) 'error doing MAP U', fnamenc
!          stop 
!      end if
!
!      varname = 'MAPFAC_V'
!      lendim_exp(1) = nxn(l)
!      lendim_max(1) = nxmaxn
!      lendim_exp(2) = nyn(l)+1
!      lendim_max(2) = nymaxn
!
!      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
!          varname, m_vn(0,0,1,l), &
!          itime, &
!          ndims, ndims_exp, ndims_max, &
!          lendim, lendim_exp, lendim_max )
!      if (ierr .ne. 0) then
!          write(*,9100) 'error doing MAP V', fnamenc
!          stop
!      end if
!
!      varname = 'MAPFAC_M'
!      lendim_exp(1) = nxn(l)
!      lendim_max(1) = nxmaxn
!      lendim_exp(2) = nyn(l)
!      lendim_max(2) = nymaxn
!
!      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
!          varname, m_w(0,0,1), &
!          itime, &
!          ndims, ndims_exp, ndims_max, &
!          lendim, lendim_exp, lendim_max )
!      if (ierr .ne. 0) then
!          write(*,9100) 'error doing MAP W', fnamenc
!          stop
!      end if





! calculate water vapor pressure in mb, from sfc pressure
!   and 2 m mixing ratio
      iduma = 0
      do j = 0, nyn(l)-1
      do i = 0, nxn(l)-1
! 29-nov-2005 - added this to catch occasional tt2n=0.0 values
          duma = max( 100.0, tthn(i,j,kbgn,n,l)-50.0 )
          if (tt2n(i,j,1,n,l) .le. duma) then
              iduma = iduma + 1
              if (iduma .eq. 1) then
                  write(*,*) 'readwind_nests - bad tt2n at'
                  write(*,*) 'l, i, j, tt2n =', l, i, j, tt2n(i,j,1,n,l)
              end if
              tt2n(i,j,1,n,l) = tthn(i,j,kbgn,n,l)
              td2n(i,j,1,n,l) = qvhn(i,j,kbgn,n,l)
          end if
          duma = td2n(i,j,1,n,l)/0.622
          ewater_mb = 0.01*( 0.99976*psn(i,j,1,n,l)*duma/(1.0+duma) )
          esatwater_mb = 0.01*ew(tt2n(i,j,1,n,l))
          ewater_mb = max( 1.0e-10, min( esatwater_mb, ewater_mb ) )
! then use the following, which is from an old 1970's report
!   (reference not available, but the formula works)
!   tdew(in C) = (4318.76/(19.5166 - ln(ewater(in mb)))) - 243.893
          td2n(i,j,1,n,l) = 273.16 + &
                 (4318.76/(19.5166 - log(ewater_mb))) - 243.893
      end do
      end do
      if (iduma .gt. 0) write(*,*) &
          'readwind_nests - bad tt2n count =', iduma


! sea level pressure - calculate it from surface pressure and 
!    ground elevation using standard atmosphere relations
      do j = 0, nyn(l)-1
      do i = 0, nxn(l)-1
          msln(i,j,1,n,l) = psn(i,j,1,n,l)/ &
                  ((1.0 - 6.5e-3*oron(i,j,l)/288.0)**5.2553)
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
      do j = 0, nyn(l)-1
      do i = 0, nxn(l)-1
          lsprecn(i,j,1,n,l) = 0.0
          convprecn(i,j,1,n,l) = 0.0
          tccn(i,j,1,n,l) = 0.0
      end do
      end do

!C
! Large scale precipitation, (accumulated value, mm)

      varname = 'RAINNC'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
          varname, lsprecn(0,0,1,n,l), &
          itime, &
          ndims, ndims_exp, ndims_max, &
          lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
      write(*,9100) l, 'error doing ncread of RAINNC, set to zero', fnamenc
          do j = 0, nyn(l)-1
          do i = 0, nxn(l)-1
              lsprecn(i,j,1,n,l) = 0.0
          end do
          end do
      end if

!
! Convective precipitation, (accumulated value, mm)
 
      varname = 'RAINC'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
          varname, convprecn(0,0,1,n,l), &
          itime, &
          ndims, ndims_exp, ndims_max, &
          lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
      write(*,9100) l, 'error doing ncread of RAINC, set to zero', fnamenc
          do j = 0, nyn(l)-1
          do i = 0, nxn(l)-1
              convprecn(i,j,1,n,l) = 0.0
          end do
          end do
      end if

! CLOUD FRACTION (clound cover)
 
      varname = 'CLDFRA'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
          varname, tccn(0,0,1,n,l), &
          itime, &
          ndims, ndims_exp, ndims_max, &
          lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
!     write(*,9100) l, 'error doing ncread of CLDFRA, set to zero', fnamenc
          do j = 0, nyn(l)-1
          do i = 0, nxn(l)-1
              tccn(i,j,1,n,l) = 0.0
          end do
          end do
      end if


! land use, added 2015-03-27 //AD
      if (lu_option.eq.1) then  ! Read land-use from WRF?
        varname = 'LU_INDEX'
        call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
          varname, landuse_wrf(0,0), &
          itime, &
          ndims, ndims_exp, ndims_max, &
          lendim, lendim_exp, lendim_max )
        do k = 1, numclass  ! Re-initialize xlandusen
          do j = 0, nyn(l)-1
            do i = 0, nxn(l)-1
              xlandusen(i,j,k,l) = 0.0
            enddo
          enddo
        enddo
        if (num_land_cat.eq.24) then ! Assume USGS landuse data
          do j = 0, nyn(l)-1
            do i = 0, nxn(l)-1
              k = nint(landuse_wrf(i,j))  ! Safely convert element to integer (nearest)
              select case (k) ! Translate USGS categories to Wesely-types
                case(1)   ! USGS: Urban and built-up land 
                  xlandusen(i,j,1,l)  = 1.  ! Wesely: Urban land
                case(2:4)   ! USGS: Any cropland, pasture
                  xlandusen(i,j,2,l)  = 1.  ! Wesely: Agricultural land
                case(5)     ! USGS: Cropland/grassland mosaic
                  xlandusen(i,j,10,l) = 1.  ! Wesely: mixed aggricultural and range land
                case(6)     ! USGS: Cropland/woodland mosaic
                  xlandusen(i,j,4,l)  = 1.  ! Wesely: Deciduous Forest
                case(7)     ! USGS: Grassland
                  xlandusen(i,j,3,l)  = 1.  ! Wesely: Range land
                case(8)     ! USGS: Shrubland
                  xlandusen(i,j,11,l) = 1.  ! Wesely: rocky open areas with growing shrubs
                case(9)     ! USGS: Mixed Grassland and Shrubland
                  xlandusen(i,j,3,l)  = .5  ! Wesely: Range land
                  xlandusen(i,j,11,l) = .5  ! Wesely: rocky open areas with growing shrubs
                case(10)    ! USGS: Savanna
                  xlandusen(i,j,11,l) = 1.  ! Wesely: rocky open areas with growing shrubs
                case(11,12) ! USGS: Deciduous broadleaf (or needle leaf) forests
                  xlandusen(i,j,4,l)  = 1.  ! Wesely: Deciduous forest
                case(13)    ! USGS: Evergreen broadleaf forest
                  xlandusen(i,j,13,l)  = 1. ! Wesely: Rainforest
                case(14)    ! USGS: Evergreen needleleaf forest
                  xlandusen(i,j,5,l)  = 1.  ! Wesely: Coniferous Forest
                case(15)    ! USGS: Mixed forest
                  xlandusen(i,j,6,l)  = 1.  ! Wesely: Mixed forest including wetland
                case(16)    ! USGS: Water bodies
                  xlandusen(i,j,7,l)  = 1.  ! Wesely: water, both salt and fresh
                case(17,20) ! USGS: Herbaceous wetland or tundra
                  xlandusen(i,j,9,l)  = 1.  ! Wesely: nonforested wetland
                case(18,21) ! USGS: Wooded wetland or tundra
                  xlandusen(i,j,6,l)  = 1.  ! Wesely: Mixed Forest including Wetland
                case(19,23) ! USGS: "Barren or sparsely vegetated" or "Bare ground Tundra"
                  xlandusen(i,j,8,l)  = 1.  ! Wesely: Barren land mostly Desert
                case(22)    ! USGS: Mixed Tundra
                  xlandusen(i,j,6,l)  = .5  ! Wesely: Mixed Forest including Wetland
                  xlandusen(i,j,8,l)  = .5  ! Wesely: Barren land mostly Desert
                case(24)    ! USGS: Snow or Ice
                  xlandusen(i,j,12,l) = 1.  ! Wesely: snow or ice
                case default
                  write(*,9135) l, 'Unknown land use category',  &
                    'i,j,cat =', i,j,k, fnamenc
              endselect
            enddo
          enddo
        elseif (num_land_cat.eq.20) then ! Assume MODIS data
          do j = 0, nyn(l)-1
            do i = 0, nxn(l)-1
              k = nint(landuse_wrf(i,j))  ! Safely convert element to integer (nearest)
              select case (k) ! Translate MODIS categories to Wesely-types
                case(1)   ! IGBP: Evergreen Needleleaf Forest
                  xlandusen(i,j,5,l)  = 1.  ! Wesely: Coniferous Forest
                case(2)   ! IGBP: Evergreen Broadleaf Forest
                  xlandusen(i,j,13,l) = 1.  ! Wesely: Rainforest
                case(3,4) ! IGBP: Deciduous Needleleaf or Broadleaf Forest
                  xlandusen(i,j,4,l)  = 1.  ! Wesely: Deciduous Forest
                case(5)   ! IGBP: Mixed Forest
                  xlandusen(i,j,6,l)  = 1.  ! Wesely: Mixed Forest including Wetland
                case(6:9) ! IGBP: Shrublands or Savannas
                  xlandusen(i,j,11,l) = 1.  ! Wesely: rocky open areas with growing shrubs
                case(10)  ! IGBP: Grasslands
                  xlandusen(i,j,3,l)  = 1.  ! Wesely: Range land
                case(11)  ! IGBP: Permanent Wetlands
                  xlandusen(i,j,9,l)  = 1.  ! Wesely: Non-forested wetland
                case(12)  ! IGBP: Croplands
                  xlandusen(i,j,2,l)  = 1.  ! Wesely: Agricultural land
                case(13)  ! IGBP: Urban and Built-up
                  xlandusen(i,j,1,l)  = 1.  ! Wesely: Urban land
                case(14)  ! IGBP: Cropland/Natural Vegetatation Mosaic
                  xlandusen(i,j,10,l)  = 1.  ! Wesely: mixed agricultural and range land
                case(15)  ! IGBP: Snow and Ice
                  xlandusen(i,j,12,l)  = 1.  ! Wesely: Snow and Ice
                case(16)  ! IGBP: Barren or Sparsely Vegetated
                  xlandusen(i,j,8,l)  = 1.  ! Wesely: Barren land mostly Desert
                case(17)  ! IGBP: Water Bodies
                  xlandusen(i,j,7,l)  = 1.  ! Wesely: water, both salt and fresh
                case(18)  ! IGBP: Wooded Tundra (WRF only)
                  xlandusen(i,j,6,l)  = 1.  ! Wesely: Mixed Forest including Wetland
                case(19)  ! IGBP: Mixed Tundra (WRF only)
                  xlandusen(i,j,6,l)  = .5  ! Wesely: Mixed Forest including Wetland
                  xlandusen(i,j,8,l)  = .5  ! Wesely: Barren land mostly Desert
                case(20)  ! IGBP: Barren Tundra (WRF only)
                  xlandusen(i,j,8,l)  = 1.  ! Wesely: Barren land mostly Desert
                case default
                  write(*,9135) l, 'Unknown land use category',  &
                    'i,j,cat =', i,j,k, fnamenc
              endselect
            enddo
          enddo
        endif
      endif ! lu_option.eq.1

! snow depth
      varname = 'SNOWH'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, sdn(0,0,1,n,l), &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
!         write(*,9100) l, 'error doing ncread of SNOWH', fnamenc
          do j = 0, nyn(l)-1
          do i = 0, nxn(l)-1
              sdn(i,j,1,n,l) = 0.0
          end do
          end do
      end if


! surface sensible heat flux (positive <--> upwards)
      varname = 'HFX'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, sshfn(0,0,1,n,l), &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
          do j = 0, nyn(l)-1
          do i = 0, nxn(l)-1
              sshfn(i,j,1,n,l) = -sshfn(i,j,1,n,l)
          end do
          end do

      if (ierr .ne. 0) then
          write(*,9100) l, 'error doing ncread of HFX', fnamenc
          do j = 0, nyn(l)-1
          do i = 0, nxn(l)-1
              sshfn(i,j,1,n,l) = 0.0
          end do
          end do
          hflswitch=.false.    ! Heat flux is not available
      else
          hflswitch=.true.     ! Heat flux is available
! limit to values to bounds originally used by flexpart?
!         do 1502 j=0,nyn(l)-1
!         do 1502 i=0,nxn(l)-1
!            if(sshfn(i,j,1,n,l).gt.200.) sshfn(i,j,1,n,l)=200.
!            if(sshfn(i,j,1,n,l).lt.-400.) sshfn(i,j,1,n,l)=-400.
!1502     continue
      end if

! ustar
      varname = 'UST'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
          varname, ustarn(0,0,1,n,l), &
          itime, &
          ndims, ndims_exp, ndims_max, &
          lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) l, 'error doing ncread of UST', fnamenc
          do j = 0, nyn(l)
          do i = 0, nxn(l)
              ustarn(i,j,1,n,l) = 0.0
          end do
          end do
          strswitch=.false.    ! ustar is not available
      else
          strswitch=.true.     ! ustar is available
          do j=0,nyn(l)
          do i=0,nxn(l)
            surfstrn(i,j,1,n,l)=ustarn(i,j,1,n,l)/dumarray_pp(i,j,kbgn)
            enddo
            enddo
  
      end if

      if(sfc_option .eq. sfc_option_wrf) then
! pblh
      varname = 'PBLH'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
          varname, hmixn(0,0,1,n,l), &
          itime, &
          ndims, ndims_exp, ndims_max, &
          lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) l, 'error doing ncread of PBLH', fnamenc
          stop
      endif

      endif

! surface solar radiation flux (positive <--> downwards)
      varname = 'SWDOWN'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, ssrn(0,0,1,n,l), &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) l, 'error doing ncread of SWDOWN', fnamenc
          do j = 0, nyn(l)-1
          do i = 0, nxn(l)-1
              ssrn(i,j,1,n,l) = 0.0
          end do
          end do
      else
          do j = 0, nyn(l)-1
          do i = 0, nxn(l)-1
              ssrn(i,j,1,n,l) = max( ssrn(i,j,1,n,l), 0.0 )
          end do
          end do
      end if


! ew & ns surface stress
!   Doesn't appear to be any 2-d cloud cover field in the
!       wrf output.
!   For now, set to zero
      do j = 0, nyn(l)-1
      do i = 0, nxn(l)-1
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
! AD: land sea mask is NOT fixed if the simulation updates sea ice.
!     WRF labels sea ice covered areas as land!


! *** done with reading the wrf output file ***


!  print*,'uu out1',uuhn(0,259,1:10,1)
!  print*,'mu out1',mu(0,259,1),mub(0,259,1)
!  print*,'m_xn out1',m_xn(0,259,1,1),m_yn(0,259,1,1)


! interpolate uuh from the "U-grid" to the "T-grid"
! interpolate vvh from the "V-grid" to the "T-grid"
      if (wind_option.le.0) then
      do k = kbgn, nuvz
      do j = 0, nyn(l)-1
      do i = 0, nxn(l)-1
      if (wind_option.lt.0) then
      divhn(i,j,k,l)=(uuhn(i+1,j,k,l)-uuhn(i,j,k,l))/dxn(l)*m_xn(i,j,1,l) &
       +(vvhn(i,j+1,k,l)-vvhn(i,j,k,l))/dyn(l)*m_yn(i,j,1,l)
      endif
          uuhn(i,j,k,l) = 0.5*(uuhn(i,j,k,l) + uuhn(i+1,j,k,l))
          vvhn(i,j,k,l) = 0.5*(vvhn(i,j,k,l) + vvhn(i,j+1,k,l))
      end do
      end do
      end do
      elseif (wind_option.eq.1) then
      do k = kbgn, nuvz
      do j = 0, nyn(l)-1
      do i = 0, nxn(l)-1
          uuhn(i,j,k,l) = 0.5*(uuhn(i,j,k,l) + uuhn(i+1,j,k,l))
          vvhn(i,j,k,l) = 0.5*(vvhn(i,j,k,l) + vvhn(i,j+1,k,l))
      mu2=mu(i,j,1)+mub(i,j,1)
      uuhn(i,j,k,l) = uuhn(i,j,k,l)/mu2 !*m_yn(i,j,1,l)
      vvhn(i,j,k,l) = vvhn(i,j,k,l)/mu2 !*m_xn(i,j,1,l)
      wwhn(i,j,k,l) = wwhn(i,j,k,l)/mu2 !*m_yn(i,j,1,l)
      end do
      end do
      end do
      elseif (wind_option.eq.2) then
      do k = kbgn, nuvz
      do j = 0, nyn(l)-1
      do i = 0, nxn(l)-1
          uuhn(i,j,k,l) = 0.5*(uuhn(i,j,k,l) + uuhn(i+1,j,k,l))
          vvhn(i,j,k,l) = 0.5*(vvhn(i,j,k,l) + vvhn(i,j+1,k,l))
      mu2=mu(i,j,1)+mub(i,j,1)
      wwhn(i,j,k,l) = wwhn(i,j,k,l)/mu2 !*m_yn(i,j,1,l)
      end do
      end do
      end do

      endif

!  print*,'uu out2',uuhn(0,259,1:10,1)
!  print*,'mu out2',mu(0,259,1),mub(0,259,1)
!  print*,'m_xn out2',m_xn(0,259,1,1),m_yn(0,259,1,1)

! CALCULATE SURFSTR
      if(sfc_option .eq. sfc_option_diagnosed) then
        do j=0,nyn(l)-1
        do i=0,nxn(l)-1
        surfstrn(i,j,1,n,l)=sqrt(ewss(i,j)**2+nsss(i,j)**2) 
       enddo
       enddo
        strswitch=.false.    ! Surface stress is not available
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
        do j=0,nyn(l)-1
          do i=0,nxn(l)-1
!           plev1=akz(3)+bkz(3)*psn(i,j,1,n,l)
            plev1=pphn(i,j,k,n,l)
            pmean=0.5*(psn(i,j,1,n,l)+plev1)
            tv=tthn(i,j,k,n,l)*(1.+0.61*qvhn(i,j,k,n,l))
            fu=-r_air*tv/ga/pmean
            hlev1=fu*(plev1-psn(i,j,1,n,l))   ! HEIGTH OF FIRST MODEL LAYER
            ff10m= sqrt(u10n(i,j,1,n,l)**2+v10n(i,j,1,n,l)**2)
            fflev1=sqrt(uuhn(i,j,k,l)**2+vvhn(i,j,k,l)**2)
            call pbl_profile(psn(i,j,1,n,l),td2n(i,j,1,n,l),hlev1, &
                             tt2n(i,j,1,n,l),tthn(i,j,k,n,l), &
                             ff10m,fflev1, &
                             surfstrn(i,j,1,n,l),sshfn(i,j,1,n,l))
            if(sshfn(i,j,1,n,l).gt.200.) sshfn(i,j,1,n,l)=200.
            if(sshfn(i,j,1,n,l).lt.-400.) sshfn(i,j,1,n,l)=-400.
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
      do  j = 0, nyn(l)-1
      do  i = 0, nxn(l)-1
          uuhn(i,j,1,l)   = u10n(i,j,1,n,l)
          vvhn(i,j,1,l)   = v10n(i,j,1,n,l)
          tthn(i,j,1,n,l) = tt2n(i,j,1,n,l)
         ptthn(i,j,1,n,l) = ptthn(i,j,2,n,l)
          qvhn(i,j,1,n,l) = qvhn(i,j,2,n,l)
         tkehn(i,j,1,n,l) =tkehn(i,j,2,n,l)
! pressure at 2 m AGL
          pphn(i,j,1,n,l) = 0.99976*psn(i,j,1,n,l)
! height (MSL) at ground level (shift it down)
          zzhn(i,j,1,n,l) = zzhn(i,j,2,n,l)
! height (MSL) at top of the added level
          zzhn(i,j,2,n,l) = zzhn(i,j,1,n,l) + 4.0
      if (hmixn(i,j,1,n,l).lt.hmixmin) hmixn(i,j,1,n,l)=hmixmin

      enddo
      enddo
      end if


       do i=0,nxn(L)-1
        do j=0,nyn(L)-1
         do k=1,nuvzmax
           un_wrf(i,j,k,n,l)=uuhn(i,j,k,l)
           vn_wrf(i,j,k,n,l)=vvhn(i,j,k,l)
         enddo
        enddo
       enddo
 
       do i=0,nxn(L)-1
        do j=0,nyn(L)-1
         do k=1,nwzmax
           wn_wrf(i,j,k,n,l)=wwhn(i,j,k,l)
         enddo
        enddo
       enddo




      enddo !loop over the nests




      return    
      end subroutine readwind_nests
