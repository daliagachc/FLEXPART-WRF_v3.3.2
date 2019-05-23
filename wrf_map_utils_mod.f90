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
! file include_map_utils - created 22-nov-2005
!
!   this file contains the "module data" from 
!	file .../wrfsi2.1/src/mod/module_map_utils.F
!-----------------------------------------------------------------------

module wrf_map_utils_mod
      implicit none

      logical :: proj_init, proj_cyclic

      integer :: proj_code, proj_nx, proj_ny
      integer :: proj_latlon, proj_merc, proj_lc, proj_ps, proj_rotlat
      real :: proj_lat1, proj_lon1, proj_dx, &
           proj_dlat, proj_dlon, proj_clat, proj_clon, &
           proj_stdlon, proj_truelat1, proj_truelat2, &
           proj_hemi, proj_cone, proj_polei, proj_polej, &
           proj_rsw, proj_rebydx
      real :: pi, deg_per_rad, rad_per_deg, earth_radius_m

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! define data structures to define various projections
!
! type proj_info
!
!   logical        proj_init     ! flag to indicate if this struct is 
!                                ! ready for use
!   logical        proj_cyclic   ! flag indicating if this grid
!                                ! is cyclic in the longitudinal
!                                ! direction...happens with
!                                ! global lat/lon grids like gfs/avn
!   integer        proj_code     ! integer code for projection type
!   integer        proj_nx
!   integer        proj_ny
!   real           proj_lat1    ! sw latitude (1,1) in degrees (-90->90n)
!   real           proj_lon1    ! sw longitude (1,1) in degrees (-180->180e)
!   real           proj_dx       ! grid spacing in meters at truelats, used
!                                ! only for ps, lc, and merc projections
!   real           proj_dlat     ! lat increment for lat/lon grids
!   real           proj_dlon     ! lon increment for lat/lon grids
!   real           proj_clat     ! center latitude of grid
!   real           proj_clon     ! center longitude of grid
!   real           proj_stdlon   ! longitude parallel to y-axis (-180->180e)
!   real           proj_truelat1 ! first true latitude (all projections)
!   real           proj_truelat2 ! second true lat (lc only)
!   real           proj_hemi     ! 1 for nh, -1 for sh
!   real           proj_cone     ! cone factor for lc projections
!   real           proj_polei    ! computed i-location of pole point
!   real           proj_polej    ! computed j-location of pole point
!   real           proj_rsw      ! computed radius to sw corner
!   real           proj_rebydx   ! earth radius divided by dx
!
! end type proj_info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module wrf_map_utils_mod
