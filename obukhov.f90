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
!     real function obukhov(ps,tsurf,tdsurf,tlev,ustar,hf,akm,bkm)
      real function obukhov(ps,tsurf,tdsurf,tlev,ustar,hf,plev)
!********************************************************************
!                                                                   *
!                       Author: G. WOTAWA                           *
!                       Date:   1994-06-27                          *
!                                                                   *
!  Update: A. Stohl, 2000-09-25, avoid division by zero by          *
!  setting ustar to minimum value                                   *
!                                                                   *
!  20 Oct 2005 - R. Easter - pass plev as argument instead of       *
!                            akm,bkm                                *
!                                                                   *
!********************************************************************
!                                                                   *
!     This program calculates Obukhov scale height from surface     *
!     meteorological data and sensible heat flux.                   *
!                                                                   *
!********************************************************************
!                                                                   *
!     INPUT:                                                        *
!                                                                   *
!     ps      surface pressure [Pa]                                 *
!     tsurf   surface temperature [K]                               *
!     tdsurf  surface dew point [K]                                 *
!     tlev    temperature first model level [K]                     *
!     ustar   scale velocity [m/s]                                  *
!     hf      surface sensible heat flux [W/m2]                     *
!     akm     ECMWF vertical discretization parameter               *
!     bkm     ECMWF vertical discretization parameter               *
!                                                                   *
!********************************************************************
 
!      include 'includepar'
!     real akm(nwzmax),bkm(nwzmax)
!      real ps,tsurf,tdsurf,tlev,ustar,hf,e,ew,tv,rhoa,plev
!      real ak1,bk1,theta,thetastar
  use par_mod

  implicit none

  real :: ps,tsurf,tdsurf,tlev,ustar,hf,e,ew,tv,rhoa,plev
  real :: theta,thetastar


      e=ew(tdsurf)                           ! vapor pressure
      tv=tsurf*(1.+0.378*e/ps)               ! virtual temperature
      rhoa=ps/(r_air*tv)                      ! air density
!     ak1=(akm(1)+akm(2))/2.
!     bk1=(bkm(1)+bkm(2))/2.
!     plev=ak1+bk1*ps                        ! Pressure level 1
      theta=tlev*(100000./plev)**(r_air/cpa) ! potential temperature
      if (ustar.le.0.) ustar=1.e-8
      thetastar=hf/(rhoa*cpa*ustar)           ! scale temperature
      if(abs(thetastar).gt.1.e-10) then
         obukhov=theta*ustar**2/(karman*ga*thetastar)
      else
         obukhov=9999                        ! zero heat flux
      endif
      if (obukhov.gt. 9999.) obukhov= 9999.
      if (obukhov.lt.-9999.) obukhov=-9999.

end function obukhov

