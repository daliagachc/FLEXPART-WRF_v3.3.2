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
!**********************************************************************
! FLEXPART_WRF SOURCE FILE MAP_PROJ_WRF - CONTAINS                    *
!       subroutine xyindex_to_ll_wrf                                  *
!       subroutine xymeter_to_ll_wrf                                  *
!       subroutine ll_to_xyindex_wrf                                  *
!       subroutine ll_to_xymeter_wrf                                  *
!       subroutine test_xyindex_to_ll_wrf                             *
!                                                                     * 
!**********************************************************************
!                                                                     * 
! AUTHOR:      R. Easter & J. Fast, PNNL                              *
! DATE:        Dec 2005                                               *
! LAST UPDATE: same                                                   *
!                                                                     * 
!**********************************************************************
!                                                                     *
! DESCRIPTION:                                                        *
!                                                                     *
! Converts between "grid index" or "grid meter" coordinates           *
! and latitude-longitude (degrees)                                    *
!                                                                     *
!**********************************************************************


!-----------------------------------------------------------------------
subroutine xyindex_to_ll_wrf( lgrid, x_i, y_j, x_lon, y_lat )
!
!   calculates longitude/latitude from xy "grid index" coordinates
!
!   arguments
!	lgrid - input - grid identifier
!	x_i, y_j - input = "grid index" coordinates on grid "lgrid".  
!		x_i ranges from 0 to nx-1.  y_j ranges from 0 to nj-1.
!	x_lon, y_lat - output = longitude and latitude in degrees
!
!   *** note ***
!	if x_i is outside [-grace,  +grace+nx]
!	or y_j is outside [-grace,  +grace+ny]
!   the routine writes an error message and halts
!
!	include 'includepar'
!	include 'includecom'
  use par_mod
  use com_mod

!   arguments
         implicit none
	integer :: lgrid
	real :: x_i, y_j, x_lon, y_lat

!   local variables
	integer :: ia, ib, ja, jb
	real :: dumxi, dumyj
	real :: fx, fy
	real,parameter :: grace=2.01


!   first check that x_i & y_j are "within bounds"
	if (lgrid .le. 0) then
	    if ((x_i .lt. -grace) .or. (x_i .gt. grace+real(nx-1)) .or. &
      	        (y_j .lt. -grace) .or. (y_j .gt. grace+real(ny-1))) then
		write(*,'(/a/a,i4,1p,2e12.3)') &
      		    '*** xyindex_to_ll_wrf error -- bad inputs', &
      		    '    lgrid, x_i, y_j =', lgrid, x_i, y_j
!	stop
	    end if
	else 
	    if ((x_i .lt. -grace                   ) .or.  &
                (x_i .gt. grace+real(nxn(lgrid)-1)) .or. &
      	        (y_j .lt. -grace                   ) .or. &
                (y_j .gt. grace+real(nyn(lgrid)-1))) then
		write(*,'(/a/a,i4,1p,2e12.3)') &
      		    '*** xyindex_to_ll_wrf error -- bad inputs', &
      		    '    lgrid, x_i, y_j =', lgrid, x_i, y_j
!	stop
	    end if
	end if

	if (map_proj_method .gt. 0) then
	    if (map_proj_id .eq. 1) goto 2000
	    if (map_proj_id .eq. 2) goto 2000
	    if (map_proj_id .eq. 3) goto 2000
	    if (map_proj_id .eq. 4) goto 2000
	end if

	if (lgrid .le. 0) then
	    goto 5000
	else
	    goto 6000
	end if

!
!   use map projection routines in map_proj_wrf_subaa.f
!
2000	continue
	if (lgrid .le. 0) then
	    dumxi = 1.0 + x_i
	    dumyj = 1.0 + y_j
	else
	    dumxi = 1.0 + (x_i/xresoln(lgrid)) + xln(lgrid)
	    dumyj = 1.0 + (y_j/yresoln(lgrid)) + yln(lgrid)
	end if
        call ij_to_latlon( dumxi, dumyj, y_lat, x_lon )
	return


!
!   do interpolation using the outer grid xlon2d,ylat2d
!
5000	continue
	if (x_i .le. 0.0) then
	    ia = 0
	else if (x_i .ge. real(nxmin1)) then
	    ia = nxmin1
	else
	    ia = ifix( x_i )
	end if
	fx = x_i - ia
	fx = max( -2.0, min( fx, 3.0 ) )
	ib = ia + 1

	if (y_j .le. 0.0) then
	    ja = 0
	else if (y_j .ge. real(nymin1)) then
	    ja = nymin1
	else
	    ja = ifix( y_j )
	end if
	fy = y_j - ja
	fy = max( -2.0, min( fy, 3.0 ) )
	jb = ja + 1

	x_lon = xlon2d(ia,ja)*(1.0-fx)*(1.0-fy) + &
      	        xlon2d(ia,jb)*(1.0-fx)*fy       + &
      	        xlon2d(ib,ja)*fx*(1.0-fy)       + &
      	        xlon2d(ib,jb)*fx*fy


	y_lat = ylat2d(ia,ja)*(1.0-fx)*(1.0-fy) + &
      	        ylat2d(ia,jb)*(1.0-fx)*fy       + &
               ylat2d(ib,ja)*fx*(1.0-fy)       + &
      	        ylat2d(ib,jb)*fx*fy

	return


!
!   do interpolation using the nested grid xlon2dn,ylat2dn
!
6000	continue
	if (x_i .le. 0.0) then
	    ia = 0
	else if (x_i .ge. real(nxn(lgrid)-1)) then
	    ia = nxn(lgrid)-1
	else
	    ia = ifix( x_i )
	end if
	fx = x_i - ia
	fx = max( -2.0, min( fx, 3.0 ) )
	ib = ia + 1

	if (y_j .le. 0.0) then
	    ja = 0
	else if (y_j .ge. real(nyn(lgrid)-1)) then
	    ja = nyn(lgrid)-1
	else
	    ja = ifix( y_j )
	end if
	fy = y_j - ja
	fy = max( -2.0, min( fy, 3.0 ) )
	jb = ja + 1

	x_lon = xlon2dn(ia,ja,lgrid)*(1.0-fx)*(1.0-fy) + &
      	        xlon2dn(ia,jb,lgrid)*(1.0-fx)*fy       + &
      	        xlon2dn(ib,ja,lgrid)*fx*(1.0-fy)       + &
      	        xlon2dn(ib,jb,lgrid)*fx*fy


	y_lat = ylat2dn(ia,ja,lgrid)*(1.0-fx)*(1.0-fy) + &
      	        ylat2dn(ia,jb,lgrid)*(1.0-fx)*fy       + &
               ylat2dn(ib,ja,lgrid)*fx*(1.0-fy)       + &
      	        ylat2dn(ib,jb,lgrid)*fx*fy

	return

	end subroutine xyindex_to_ll_wrf


!-----------------------------------------------------------------------
        subroutine xyindex_to_ll_wrf_out(lgrid,x_i,y_j,x_lon,y_lat)
!
!   calculates longitude/latitude from xy "grid index" coordinates
!
!   arguments
!	lgrid - input - grid identifier
!	x_i, y_j - input = "grid index" coordinates on grid "lgrid".  
!		x_i ranges from 0 to nx-1.  y_j ranges from 0 to nj-1.
!	x_lon, y_lat - output = longitude and latitude in degrees
!
!   *** note ***
!	if x_i is outside [-grace,  +grace+nx]
!	or y_j is outside [-grace,  +grace+ny]
!   the routine writes an error message and halts
!
!	include 'includepar'
!	include 'includecom'
  use par_mod
  use com_mod

!   arguments
         implicit none
	integer :: lgrid
	real :: x_i, y_j, x_lon, y_lat

!   local variables
	integer :: ia, ib, ja, jb
	real :: dumxi, dumyj
	real :: fx, fy
	real,parameter :: grace=2.01


!   first check that x_i & y_j are "within bounds"
!	if (lgrid .le. 0) then
!	    if ((x_i .lt. -grace) .or. (x_i .gt. grace+real(nx-1)) .or.
!     &	        (y_j .lt. -grace) .or. (y_j .gt. grace+real(ny-1))) then
!		write(*,'(/a/a,i4,1p,2e12.3)')
!     &		    '*** xyindex_to_ll_wrf error -- bad inputs',
!     &		    '    lgrid, x_i, y_j =', lgrid, x_i, y_j
!		stop
!	    end if
!	else 
!	    if ((x_i .lt. -grace                   ) .or. 
!     &          (x_i .gt. grace+real(nxn(lgrid)-1)) .or.
!     &	        (y_j .lt. -grace                   ) .or. 
!     &          (y_j .gt. grace+real(nyn(lgrid)-1))) then
!		write(*,'(/a/a,i4,1p,2e12.3)')
!     &		    '*** xyindex_to_ll_wrf error -- bad inputs',
!     &		    '    lgrid, x_i, y_j =', lgrid, x_i, y_j
!		stop
!	    end if
!	end if

	if (map_proj_method .gt. 0) then
	    if (map_proj_id .eq. 1) goto 2000
	    if (map_proj_id .eq. 2) goto 2000
	    if (map_proj_id .eq. 3) goto 2000
	    if (map_proj_id .eq. 4) goto 2000
	end if

	if (lgrid .le. 0) then
	    goto 5000
	else
	    goto 6000
	end if

!
!   use map projection routines in map_proj_wrf_subaa.f
!
2000	continue
	if (lgrid .le. 0) then
	    dumxi = 1.0 + x_i
	    dumyj = 1.0 + y_j
	else
	    dumxi = 1.0 + (x_i/xresoln(lgrid)) + xln(lgrid)
	    dumyj = 1.0 + (y_j/yresoln(lgrid)) + yln(lgrid)
	end if
        call ij_to_latlon( dumxi, dumyj, y_lat, x_lon )
	return


!
!   do interpolation using the outer grid xlon2d,ylat2d
!
5000	continue
	if (x_i .le. 0.0) then
	    ia = 0
	else if (x_i .ge. real(nxmin1)) then
	    ia = nxmin1
	else
	    ia = ifix( x_i )
	end if
	fx = x_i - ia
	fx = max( -2.0, min( fx, 3.0 ) )
	ib = ia + 1

	if (y_j .le. 0.0) then
	    ja = 0
	else if (y_j .ge. real(nymin1)) then
	    ja = nymin1
	else
	    ja = ifix( y_j )
	end if
	fy = y_j - ja
	fy = max( -2.0, min( fy, 3.0 ) )
	jb = ja + 1

	x_lon = xlon2d(ia,ja)*(1.0-fx)*(1.0-fy) + &
      	        xlon2d(ia,jb)*(1.0-fx)*fy       + &
      	        xlon2d(ib,ja)*fx*(1.0-fy)       + &
      	        xlon2d(ib,jb)*fx*fy


	y_lat = ylat2d(ia,ja)*(1.0-fx)*(1.0-fy) + &
      	        ylat2d(ia,jb)*(1.0-fx)*fy       + &
      	        ylat2d(ib,ja)*fx*(1.0-fy)       + &
      	        ylat2d(ib,jb)*fx*fy

	return


!
!   do interpolation using the nested grid xlon2dn,ylat2dn
!
6000	continue
	if (x_i .le. 0.0) then
	    ia = 0
	else if (x_i .ge. real(nxn(lgrid)-1)) then
	    ia = nxn(lgrid)-1
	else
	    ia = ifix( x_i )
	end if
	fx = x_i - ia
	fx = max( -2.0, min( fx, 3.0 ) )
	ib = ia + 1

	if (y_j .le. 0.0) then
	    ja = 0
	else if (y_j .ge. real(nyn(lgrid)-1)) then
	    ja = nyn(lgrid)-1
	else
	    ja = ifix( y_j )
	end if
	fy = y_j - ja
	fy = max( -2.0, min( fy, 3.0 ) )
	jb = ja + 1

	x_lon = xlon2dn(ia,ja,lgrid)*(1.0-fx)*(1.0-fy) + &
      	        xlon2dn(ia,jb,lgrid)*(1.0-fx)*fy       + &
      	        xlon2dn(ib,ja,lgrid)*fx*(1.0-fy)       + &
      	        xlon2dn(ib,jb,lgrid)*fx*fy


	y_lat = ylat2dn(ia,ja,lgrid)*(1.0-fx)*(1.0-fy) + &
      	        ylat2dn(ia,jb,lgrid)*(1.0-fx)*fy       + &
      	        ylat2dn(ib,ja,lgrid)*fx*(1.0-fy)       + &
      	        ylat2dn(ib,jb,lgrid)*fx*fy

	return

	end subroutine xyindex_to_ll_wrf_out


!-----------------------------------------------------------------------
	subroutine xymeter_to_ll_wrf( xmeter, ymeter, x_lon, y_lat )
!
!   calculates longitude/latitude from xy "grid meter" coordinates
!
!   arguments
!	xmeter, ymeter - input = "grid meter" coordinates on the mother grid.  
!	x_lon, y_lat - output = longitude and latitude in degrees
!
!	include 'includepar'
!	include 'includecom'
  use par_mod
  use com_mod

!   arguments
       implicit none
	real :: xmeter, ymeter, x_lon, y_lat

!   local variables
	real :: x_i, y_j

	 x_i = (xmeter - xmet0)/dx
	 y_j = (ymeter - ymet0)/dy
	call xyindex_to_ll_wrf( 0, x_i, y_j, x_lon, y_lat )

	return
	end subroutine xymeter_to_ll_wrf

	subroutine xymeter_to_ll_wrf_out(xmeter,ymeter,x_lon,y_lat)
!
!   calculates longitude/latitude from xy "grid meter" coordinates
!
!   arguments
!	xmeter, ymeter - input = "grid meter" coordinates on the mother grid.  
!	x_lon, y_lat - output = longitude and latitude in degrees
!
!	include 'includepar'
!	include 'includecom'
  use par_mod
  use com_mod

!   arguments
       implicit none
	real :: xmeter, ymeter, x_lon, y_lat

!   local variables
	real :: x_i, y_j

	 x_i = (xmeter - xmet0)/dx
	 y_j = (ymeter - ymet0)/dy
	call xyindex_to_ll_wrf_out(0,x_i, y_j,x_lon,y_lat )

	return
	end subroutine xymeter_to_ll_wrf_out

!-----------------------------------------------------------------------
	subroutine ll_to_xyindex_wrf( x_lon, y_lat, x_i, y_j )
!
!   calculates xy "grid index" coordinates from longitude/latitude
!
!   arguments
!	x_lon, y_lat - input - longitude and latitude in degrees
!	x_i, y_j - output = "grid index" coordinates on the mother grid.  
!		x_i ranges from 0 to nx-1.  y_j ranges from 0 to nj-1.
!
!	include 'includepar'
!	include 'includecom'
  use par_mod
  use com_mod
       implicit none

!   arguments
	real :: x_i, y_j, x_lon, y_lat

!   local variables
	integer :: i, ii, ia, ib, ip, ipass, iijjlohi
	integer :: j, jj, ja, jb, jp
	real :: dumcos, dumlat, dumlon
	real :: dumr2, dumr2min
	real :: dumxi, dumyj
	real :: flo, fhi
	real :: glo, ghi
	real,parameter :: grace=2.01
	real :: xxcen, xxcenb, yycen, yycenb
	real :: xxyydel, xxyydelmin
	real :: x_lon_sv, y_lat_sv



	x_lon_sv = x_lon
	y_lat_sv = y_lat
	if (map_proj_method .gt. 0) then
	    if (map_proj_id .eq. 1) goto 2000
	    if (map_proj_id .eq. 2) goto 2000
	    if (map_proj_id .eq. 3) goto 2000
	    if (map_proj_id .eq. 4) goto 2000
	end if
	goto 5000

!
!   use map projection routines in map_proj_wrf_subaa.f
!
2000	continue

        call latlon_to_ij( y_lat, x_lon, dumxi, dumyj )
	x_i = dumxi - 1.0
	y_j = dumyj - 1.0
	goto 8000


!
!   do it by search/minimization of distance,
!   using the outer grid xlon2d/ylat2d values from WRF met file
!
5000	continue

!
!   first locate the i,j for which the lon,lat at
!   i+0.5,j+0.5 are closest to x_lon, y_lat
!
	dumr2min = 1.0e30
	dumcos = cos( y_lat )
	do j = 0, ny-1
	do i = 0, nx-1
	    dumlat = 0.25*( ylat2d(i,j  ) + ylat2d(i+1,j  ) +  &
                            ylat2d(i,j+1) + ylat2d(i+1,j+1) )  
	    dumlon = 0.25*( xlon2d(i,j  ) + xlon2d(i+1,j  ) +  &
                            xlon2d(i,j+1) + xlon2d(i+1,j+1) )  
	    dumr2 = (y_lat-dumlat)**2 + ((x_lon-dumlon)*dumcos)**2
	    if (dumr2 .lt. dumr2min) then
		dumr2min = dumr2
		ib = i
		jb = j
	    end if
	end do
	end do
	i = ib
	j = jb
	ip = i+1
	jp = j+1

!
!   next determine the position between i/i+1 & j/j+1
!   that is closest x_lon, y_lat
!
	xxyydelmin = 5.0e-8*max( abs(x_lon), abs(y_lat) )
	xxyydel = 1.00
	xxcen = 0.5
	yycen = 0.5
	iijjlohi = 3
	ipass = 0

!	write(*,9510)
!	write(*,9520) 0, (i+xxcen), (j+yycen), x_lon, y_lat
9510	format( / 'll_to_xyindex' )
9520	format( 'ipass, x_i, y_j, lon, lat', i3, 4f14.7 )

5200	ipass = ipass + 1
	dumr2min = 1.0e30
	do jj = -iijjlohi, iijjlohi
	do ii = -iijjlohi, iijjlohi
	    fhi = xxcen + ii*xxyydel
	    ghi = yycen + jj*xxyydel
	    flo = 1.0 - fhi
	    glo = 1.0 - ghi
	    dumlat = glo*(flo*ylat2d(i,j ) + fhi*ylat2d(ip,j )) +  &
                     ghi*(flo*ylat2d(i,jp) + fhi*ylat2d(ip,jp))  
	    dumlon = glo*(flo*xlon2d(i,j ) + fhi*xlon2d(ip,j )) +  &
                     ghi*(flo*xlon2d(i,jp) + fhi*xlon2d(ip,jp))  
	    dumr2 = (y_lat-dumlat)**2 + ((x_lon-dumlon)*dumcos)**2
	    if (dumr2 .lt. dumr2min) then
		dumr2min = dumr2
		xxcenb = fhi
		yycenb = ghi
	    end if
	end do
	end do
	xxcen = xxcenb
	yycen = yycenb
!	write(*,9520) ipass, (i+xxcen), (j+yycen), dumlon, dumlat
	if (xxyydel .gt. xxyydelmin) then
	    xxyydel = xxyydel*0.5
	    if (ipass .eq. 4) iijjlohi = 2
	    goto 5200
	end if
!	write(*,9520) ipass, (i+xxcen), (j+yycen), dumlon, dumlat

	x_i = i + xxcen
	y_j = j + yycen


!
!   check for x_i, y_j in bounds before returning
!
8000	continue
	if ((x_i .lt. -grace) .or. (x_i .gt. grace+real(nx-1)) .or. &
      	    (y_j .lt. -grace) .or. (y_j .gt. grace+real(ny-1))) then
	    write(*,'(/a/a,1p,2e12.3/a,1p,2e12.3)') &
      		'*** ll_to_xyindex_wrf error -- x_i, y_j out of bounds', &
      		'    x_lon, y_lat =', x_lon_sv, y_lat_sv, &
      		'    x_i,   y_j   =', x_i, y_j
!     pause
!        ylat2d(5000,2000)=0.     
!	    stop
	end if
	return

	end subroutine ll_to_xyindex_wrf


!-----------------------------------------------------------------------
	subroutine ll_to_xymeter_wrf( x_lon, y_lat, xmeter, ymeter )
!
!   calculates xy "grid meter" coordinates from longitude/latitude
!
!   arguments
!	x_lon, y_lat - input - longitude and latitude in degrees
!	xmeter, ymeter - output = "grid meter" coordinates on the mother grid.  
!
!	include 'includepar'
!	include 'includecom'
  use par_mod
  use com_mod
        implicit none

!   arguments
	real :: xmeter, ymeter, x_lon, y_lat

!   local variables
	real :: x_i, y_j

	call ll_to_xyindex_wrf( x_lon, y_lat, x_i, y_j )
	xmeter = xmet0 + dx*x_i
	ymeter = ymet0 + dy*y_j

	return
	end subroutine ll_to_xymeter_wrf


!-----------------------------------------------------------------------
	subroutine test_xyindex_to_ll_wrf( lgrid )
!
!   tests the map projection routines by comparing
!	lat,lon from projection routine against the
!	lat,lon from the WRF met. files
!
!   arguments
!	lgrid - input - grid identifier
!
!	include 'includepar'
!	include 'includecom'
  use par_mod
  use com_mod

!   arguments
         implicit none
	integer :: lgrid
!   local variables
	integer :: idum, ix, jdum, jy
	integer :: map_set_proj_code
	real :: dumdx, dumxi, dumyj
	real :: dumlat, dumlatb, dumlon, dumlonb
	real :: err
	real :: rmserr



!
!   check if map projection is lambert conformal or polar stereographic
!
	if (map_proj_id .eq. 1) then
	    map_set_proj_code = 3    ! lambert conformal
	else if (map_proj_id .eq. 2) then
	    map_set_proj_code = 5    ! polar stereographic
	else if (map_proj_id .eq. 3) then
	    map_set_proj_code = 1    ! mercator
	else if (map_proj_id .eq. 4) then
	    map_set_proj_code = 0    ! lat/lon
	else
	    write(*,'(/ 10(a/) )') &
      	'************************************************************', &
      	'*                                                          *', &
      	'*    WARNING - map projection is not polar sterographic    *', &
        '*              or lambert conformal                        *', &
      	'*                                                          *', &
      	'*              x,y <--> lat,lon conversions will be done   *', &
      	'*              by interpolation & searching, and will      *', &
      	'*              have limited accuracy near poles            *', &
      	'*                                                          *', &
      	'************************************************************'
	    map_proj_method = 0
	    return
	end if
	if (lgrid .gt. 0) goto 2000

!
!   make call to map projection setup routine
!
!   (The 0.999812 factor is due to different earth_radius
!    values in wrfsi code (6370.0 vs 6371.2).)
       dumdx = dx*0.999812  ! WA 2/11/10 is this right?
!      dumdx = dx*0.9987242 ! WA 2/11/10 is this right?
!      dumdx = dx*0.9979203 ! WA 2/11/10 is this right?
       coefdx=0.999812
!	write(*,'(/2a,2i5)')  &
 !     		'test_xyindex_to_ll_wrf calling map_set -- ', &
  !    		'map_proj_id, map_set_proj_code =', &
  !    		map_proj_id, map_set_proj_code
        call map_set( map_set_proj_code, &
      		ylat2d(0,0), xlon2d(0,0), dumdx,    &
      		map_stdlon, map_truelat1, map_truelat2,   &
      		nx, ny )
	map_proj_method = 1

!
!   compute lat,lon from xi,jy at 9 points (center, 4 corners,
!	and 4 boundary midpoints)  
!   compare to lon,lat read from the WRF met. file and report rmserr
!
2000	continue
	rmserr = 0.0
	do jdum = 0, 2
	do idum = 0, 2
	    if (lgrid .le. 0) then
		jy = nint( 0.5*real((ny-1)*jdum) )
		ix = nint( 0.5*real((nx-1)*idum) )
		dumyj = 1 + jy
		dumxi = 1 + ix
		dumlatb = ylat2d(ix,jy)
		dumlonb = xlon2d(ix,jy)
	    else
		jy = nint( 0.5*real((nyn(lgrid)-1)*jdum) )
		ix = nint( 0.5*real((nxn(lgrid)-1)*idum) )
		dumyj = 1.0 + (jy/yresoln(lgrid)) + yln(lgrid)
		dumxi = 1.0 + (ix/xresoln(lgrid)) + xln(lgrid)
		dumlatb = ylat2dn(ix,jy,lgrid)
		dumlonb = xlon2dn(ix,jy,lgrid)
	    end if
            call ij_to_latlon( dumxi, dumyj, dumlat, dumlon )
	    err = (dumlat-dumlatb)**2 + &
      			((dumlon-dumlonb)*cos(dumlatb))**2
!         print*,'err',idum,jdum,dumlat,dumlatb,dumlon,dumlonb
	    rmserr = rmserr + err
!       print*,dumxi, dumyj, dumlat, dumlon,dumlatb,dumlonb
	end do
	end do
	rmserr = 111.2*sqrt( rmserr/9.0 )
!        print*,'rmserr',rmserr
        if (rmserr .le. 0.02*(1.0e-3*dumdx)) then
        write(*,'(/a,i3,1pe10.2)')  &
      	'test_xyindex_to_ll_wrf -- lgrid, rmserr (km) =', lgrid, rmserr
        endif
!
!   if rms error exceeds 0.02*dx, something is wrong
!   do not use the map projection routines in this case
!
        if (rmserr .gt. 0.02*(1.0e-3*dumdx)) then 
!      print*,'pb with rmserr. try an other earth radius'
       dumdx = dx  ! try another earth radius 
       coefdx=1.

	write(*,'(/2a,2i5)')  &
      		'test_xyindex_to_ll_wrf calling map_set -- ', &
      		'map_proj_id, map_set_proj_code =', &
      		map_proj_id, map_set_proj_code
        call map_set( map_set_proj_code, &
      		ylat2d(0,0), xlon2d(0,0), dumdx,    &
      		map_stdlon, map_truelat1, map_truelat2,   &
      		nx, ny )
	map_proj_method = 1

!
!   compute lat,lon from xi,jy at 9 points (center, 4 corners,
!	and 4 boundary midpoints)  
!   compare to lon,lat read from the WRF met. file and report rmserr
!
	rmserr = 0.0
	do jdum = 0, 2
	do idum = 0, 2
	    if (lgrid .le. 0) then
		jy = nint( 0.5*real((ny-1)*jdum) )
		ix = nint( 0.5*real((nx-1)*idum) )
		dumyj = 1 + jy
		dumxi = 1 + ix
		dumlatb = ylat2d(ix,jy)
		dumlonb = xlon2d(ix,jy)
	    else
		jy = nint( 0.5*real((nyn(lgrid)-1)*jdum) )
		ix = nint( 0.5*real((nxn(lgrid)-1)*idum) )
		dumyj = 1.0 + (jy/yresoln(lgrid)) + yln(lgrid)
		dumxi = 1.0 + (ix/xresoln(lgrid)) + xln(lgrid)
		dumlatb = ylat2dn(ix,jy,lgrid)
		dumlonb = xlon2dn(ix,jy,lgrid)
	    end if
            call ij_to_latlon( dumxi, dumyj, dumlat, dumlon )
	    err = (dumlat-dumlatb)**2 + &
      			((dumlon-dumlonb)*cos(dumlatb))**2
	    rmserr = rmserr + err
	end do
	end do
	rmserr = 111.2*sqrt( rmserr/9.0 )

        write(*,'(/a,i3,1pe10.2)')  &
      	'test_xyindex_to_ll_wrf -- lgrid, rmserr (km) =', lgrid, rmserr


        endif

!
!   if rms error exceeds 0.02*dx, something is wrong
!   do not use the map projection routines in this case
!
	dumdx = dx
	if (lgrid .gt. 0) dumdx = dxn(lgrid)
	if (rmserr .gt. 0.02*(1.0e-3*dumdx)) then
    	    map_proj_method = 0
	    write(*,'(/ 9(a/) )') &
      	'************************************************************', &
      	'*                                                          *', &
      	'*  WARNING - the coordinate transfo error exceeds 0.02*dx  *', &
      	'*                                                          *', &
      	'************************************************************'
!     print*,'projection set to',map_proj_method   
      print*,'problem with the projection. wrf+flexpart stops'
       stop
	end if

	return
	end subroutine test_xyindex_to_ll_wrf


