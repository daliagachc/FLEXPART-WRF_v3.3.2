!***********************************************************************
!* Copyright 2012,2013                                                *
!* Jerome Brioude, Delia Arnold, Andreas Stohl, Wayne Angevine,       *
!* John Burkhart, Massimo Cassiani, Adam Dingwell, Richard C Easter, Sabine Eckhardt,*
!* Stephanie Evan, Jerome D Fast, Don Morton, Ignacio Pisso,          *
!* Petra Seibert, Gerard Wotawa, Caroline Forster, Harald Sodemann,   *
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
      subroutine centerofmass(xl,yl,n,xcenter,ycenter)
!                             i  i  i    o       o
!*******************************************************************************
!                                                                              *
!     Note:  This is the FLEXPART_WRF version of subroutine assignland.        *
!            The computational grid is the WRF x-y grid rather than lat-lon.   *
!                                                                              *
!   This routine calculates the center of mass of n points on the Earth.       *
!   Input are the longitudes (xl) and latitudes (yl) of the individual         *
!   points, output is the longitude and latitude of the centre of mass.        *
!                                                                              *
!     Author: A. Stohl                                                         *
!                                                                              *
!     24 January 2002                                                          *
!                                                                              *
!    26 Oct 2005, R. Easter - changes associated with WRF horizontal grid.     *
!                             Since x & y coords are in meters,                *
!                             so just sum/average the xl & yl.                 *
!                                                                              *
!*******************************************************************************

  use par_mod

  implicit none

  integer :: n,l
  real :: xl(n),yl(n),xll,yll,xav,yav,zav,x,y,z,xcenter,ycenter


  xav=0.
  yav=0.
  zav=0.

  do l=1,n

! Convert longitude and latitude from degrees to radians
!*******************************************************

! for FLEXPART_WRF, x & y coords are in meters, 
! so just sum/average the xl & yl
!       xll=xl(l)*pi180
!       yll=yl(l)*pi180

! Calculate 3D coordinates from longitude and latitude
!*****************************************************

! for FLEXPART_WRF, this isn't necessary
!       x = cos(yll)*sin(xll)
!       y = -1.*cos(yll)*cos(xll)
!       z = sin(yll)


! Find the mean location in Cartesian coordinates
!************************************************

!       xav=xav+x
!       yav=yav+y
!       zav=zav+z
        xav=xav+xl(l)
        yav=yav+yl(l)

    enddo

      xav=xav/real(n)
      yav=yav/real(n)
!     zav=zav/float(n)


! Project the point back onto Earth's surface
!********************************************

! for FLEXPART_WRF, this isn't necessary
!     xcenter=atan2(xav,-1.*yav)
!     ycenter=atan2(zav,sqrt(xav*xav+yav*yav))

! Convert back to degrees
!************************

! for FLEXPART_WRF, this isn't necessary
!     xcenter=xcenter/pi180
!     ycenter=ycenter/pi180
      xcenter=xav
      ycenter=yav

end subroutine centerofmass

