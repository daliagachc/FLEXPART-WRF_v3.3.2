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
      subroutine interpol_rain_nests(yy1,yy2,yy3,iy1,iy2,nxmaxn,nymaxn,nzmax, &
      maxnests,ngrid,nxn,nyn,memind,xt,yt,level,itime1,itime2,itime, &
      yint1,yint2,yint3,intiy1,intiy2,icmv)
!                                     i   i   i    i      i      i
!        i       i    i   i    i    i  i    i     i      i      i
!       o     o     o
!****************************************************************************
!                                                                           *
!  Interpolation of meteorological fields on 2-d model layers for nested    *
!  grids. This routine is related to levlin3interpol.f for the mother domain*
!                                                                           *
!  In horizontal direction bilinear interpolation interpolation is used.    *
!  Temporally a linear interpolation is used.                               *
!  Three fields are interpolated at the same time.                          *
!                                                                           *
!  This is a special version of levlininterpol to save CPU time.            *
!                                                                           *
!  1 first time                                                             *
!  2 second time                                                            *
!                                                                           *
!                                                                           *
!     Author: A. Stohl                                                      *
!                                                                           *
!     15 March 2000                                                         *
!                                                                           *
!****************************************************************************
!                                                                           *
! Variables:                                                                *
!                                                                           *
! dt1,dt2              time differences between fields and current position *
! dz1,dz2              z distance between levels and current position       *
! height(nzmax)        heights of the model levels                          *
! indexh               help variable                                        *
! indz                 the level closest to the current trajectory position *
! indzh                help variable                                        *
! itime                current time                                         *
! itime1               time of the first wind field                         *
! itime2               time of the second wind field                        *
! ix,jy                x,y coordinates of lower left subgrid point          *
! level                level at which interpolation shall be done           *
! memind(3)            points to the places of the wind fields              *
! nx,ny                actual field dimensions in x,y and z direction       *
! nxmax,nymax,nzmax    maximum field dimensions in x,y and z direction      *
! xt                   current x coordinate                                 *
! yint                 the final interpolated value                         *
! yt                   current y coordinate                                 *
! yy(0:nxmax,0:nymax,nzmax,3) meteorological field used for interpolation   *
! zt                   current z coordinate                                 *
!                                                                           *
!  Changed 10/22/2007  yy1,yy2 are accumulated rain (mm)
!                      convert them into hourly rain (mm/hr)
!****************************************************************************

    implicit none

    integer :: maxnests,ngrid
    integer :: nxn(maxnests),nyn(maxnests),nxmaxn,nymaxn,nzmax,memind(2)
    integer :: m,ix,jy,ixp,jyp,itime,itime1,itime2,level,indexh,i1,i2
    integer :: ip1,ip2,ip3,ip4
    integer :: intiy1,intiy2,ipsum,icmv
    real :: yy1(0:nxmaxn-1,0:nymaxn-1,nzmax,2,maxnests)
    real :: yy2(0:nxmaxn-1,0:nymaxn-1,nzmax,2,maxnests)
    real :: yy3(0:nxmaxn-1,0:nymaxn-1,nzmax,2,maxnests)
    integer :: iy1(0:nxmaxn-1,0:nymaxn-1,2,maxnests)
    integer :: iy2(0:nxmaxn-1,0:nymaxn-1,2,maxnests)
    real :: ddx,ddy,rddx,rddy,dt1,dt2,dt,y1(2),y2(2),y3(2),yi1(2),yi2(2)
    real :: xt,yt,yint1,yint2,yint3,yint4,p1,p2,p3,p4



! If point at border of grid -> small displacement into grid
!***********************************************************

      if (xt.ge.real(nxn(ngrid)-1)) xt=real(nxn(ngrid)-1)-0.00001
      if (yt.ge.real(nyn(ngrid)-1)) yt=real(nyn(ngrid)-1)-0.00001



!**********************************************************************
! 1.) Bilinear horizontal interpolation
! This has to be done separately for 2 fields (Temporal)
!*******************************************************

! Determine the lower left corner and its distance to the current position
!************************************************************************* 

      ix=int(xt)
      jy=int(yt)
      ixp=ix+1
      jyp=jy+1
      ddx=xt-real(ix)
      ddy=yt-real(jy)
      rddx=1.-ddx
      rddy=1.-ddy
      p1=rddx*rddy
      p2=ddx*rddy
      p3=rddx*ddy
      p4=ddx*ddy
     

! Loop over 2 time steps
!***********************

! y1 and y2 are accumulated rain, need change to hourly rain
 
        i1=memind(1)
        i2=memind(2)
! time interval between two fields, second to hour
        dt=real(itime2-itime1)/3600.0
 
        yint1=p1*(yy1(ix ,jy ,level,i2,ngrid)- & 
                  yy1(ix ,jy ,level,i1,ngrid)) &
            + p2*(yy1(ixp,jy ,level,i2,ngrid)- &
                  yy1(ixp,jy ,level,i1,ngrid)) &
            + p3*(yy1(ix ,jyp,level,i2,ngrid)- &
                  yy1(ix ,jyp,level,i1,ngrid)) &
            + p4*(yy1(ixp,jyp,level,i2,ngrid)- &
                  yy1(ixp,jyp,level,i1,ngrid))
        yint1=yint1/dt
 
        yint2=p1*(yy2(ix ,jy ,level,i2,ngrid)- &
                  yy2(ix ,jy ,level,i1,ngrid)) &
            + p2*(yy2(ixp,jy ,level,i2,ngrid)- &
                  yy2(ixp,jy ,level,i1,ngrid)) &
            + p3*(yy2(ix ,jyp,level,i2,ngrid)- &
                  yy2(ix ,jyp,level,i1,ngrid)) &
            + p4*(yy2(ixp,jyp,level,i2,ngrid)- &
                  yy2(ixp,jyp,level,i1,ngrid))
       yint2=yint2/dt
 
 
! Y3 is cloud fraction in an hour

 
      do m=1,2
        indexh=memind(m)
        
        y3(m)=p1*yy3(ix ,jy ,level,indexh,ngrid) &
            + p2*yy3(ixp,jy ,level,indexh,ngrid) &
            + p3*yy3(ix ,jyp,level,indexh,ngrid) &
            + p4*yy3(ixp,jyp,level,indexh,ngrid)
        enddo


! CDA new clouds

      do m=1,2
        indexh=memind(m)

        ip1=1
        ip2=1
        ip3=1
        ip4=1
        if (iy1(ix ,jy ,indexh,ngrid) .eq. icmv) ip1=0
        if (iy1(ixp,jy ,indexh,ngrid) .eq. icmv) ip2=0
        if (iy1(ix ,jyp,indexh,ngrid) .eq. icmv) ip3=0
        if (iy1(ixp,jyp,indexh,ngrid) .eq. icmv) ip4=0
        ipsum= ip1+ip2+ip3+ip4
        if (ipsum .eq. 0) then
          yi1(m)=icmv
        else
          yi1(m)=(ip1*p1*iy1(ix ,jy ,indexh,ngrid) &
            + ip2*p2*iy1(ixp,jy ,indexh,ngrid) &
            + ip3*p3*iy1(ix ,jyp,indexh,ngrid) &
            + ip4*p4*iy1(ixp,jyp,indexh,ngrid))/ipsum
        endif

        ip1=1
        ip2=1
        ip3=1
        ip4=1
        if (iy2(ix ,jy ,indexh,ngrid) .eq. icmv) ip1=0
        if (iy2(ixp,jy ,indexh,ngrid) .eq. icmv) ip2=0
        if (iy2(ix ,jyp,indexh,ngrid) .eq. icmv) ip3=0
        if (iy2(ixp,jyp,indexh,ngrid) .eq. icmv) ip4=0
        ipsum= ip1+ip2+ip3+ip4
        if (ipsum .eq. 0) then
          yi2(m)=icmv
        else
          yi2(m)=(ip1*p1*iy2(ix ,jy ,indexh,ngrid) &
            + ip2*p2*iy2(ixp,jy ,indexh,ngrid) &
            + ip3*p3*iy2(ix ,jyp,indexh,ngrid) &
            + ip4*p4*iy2(ixp,jyp,indexh,ngrid))/ipsum
        endif
        ip1=1
        ip2=1
        ip3=1
        ip4=1
        if (iy2(ix ,jy ,indexh,ngrid) .eq. icmv) ip1=0
        if (iy2(ixp,jy ,indexh,ngrid) .eq. icmv) ip2=0
        if (iy2(ix ,jyp,indexh,ngrid) .eq. icmv) ip3=0
        if (iy2(ixp,jyp,indexh,ngrid) .eq. icmv) ip4=0
        ipsum= ip1+ip2+ip3+ip4
        if (ipsum .eq. 0) then
          yi2(m)=icmv
        else
          yi2(m)=(ip1*p1*iy2(ix ,jy ,indexh,ngrid) &
            + ip2*p2*iy2(ixp,jy ,indexh,ngrid) &
            + ip3*p3*iy2(ix ,jyp,indexh,ngrid) &
            + ip4*p4*iy2(ixp,jyp,indexh,ngrid))/ipsum
        endif
      enddo
!CPS end clouds


10    continue


!************************************
! 2.) Temporal interpolation (linear)
!************************************

      dt1=real(itime-itime1)
      dt2=real(itime2-itime)
      dt=dt1+dt2

!      yint1=(y1(1)*dt2+y1(2)*dt1)/dt
!      yint2=(y2(1)*dt2+y2(2)*dt1)/dt
      yint3=(y3(1)*dt2+y3(2)*dt1)/dt


!CPS clouds:
      intiy1=(yi1(1)*dt2 + yi1(2)*dt1)/dt
      if (yi1(1) .eq. float(icmv)) intiy1=yi1(2)
      if (yi1(2) .eq. float(icmv)) intiy1=yi1(1)

      intiy2=(yi2(1)*dt2 + yi2(2)*dt1)/dt
      if (yi2(1) .eq. float(icmv)) intiy2=yi2(2)
      if (yi2(2) .eq. float(icmv)) intiy2=yi2(1)

      if (intiy1 .ne. icmv .and. intiy2 .ne. icmv) then
        intiy2 = intiy1 + intiy2 ! convert cloud thickness to cloud top
      else
        intiy1=icmv
        intiy2=icmv
      endif
!CPS end clouds


end subroutine interpol_rain_nests

