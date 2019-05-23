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

subroutine assignland
!*******************************************************************************
!                                                                              *
!     Note:  This is the FLEXPART_WRF version of subroutine assignland.        *
!            The computational grid is the WRF x-y grid rather than lat-lon.   *
!                                                                              *
!     This routine assigns fractions of the 8 landuse classes to each ECMWF    *
!     grid point.                                                              *
!     The landuse inventory of                                                 *
!                                                                              *
!     van de Velde R.J., Faber W.S., van Katwijk V.F., Kuylenstierna J.C.I.,   *
!     Scholten H.J., Thewessen T.J.M., Verspuij M., Zevenbergen M. (1994):     *
!     The Preparation of a European Land Use Database. National Institute of   *
!     Public Health and Environmental Protection, Report nr 712401001,         *
!     Bilthoven, The Netherlands.                                              *
!                                                                              *
!     is used to create a detailed landuse inventory for Europe.               *
!                                                                              *
!     Outside of Europe, the ECMWF land/sea mask is used to distinguish        *
!     between sea (-> ocean) and land (-> grasslands).                         *
!                                                                              *
!     Author: A. Stohl                                                         *
!                                                                              *
!     5 December 1996                                                          *
!     8 February 1999 Additional use of nests, A. Stohl                        *
!                                                                              *
!    14 October  2005 R. Easter -- modified for WRF.                           *
!                     The landuse inventory is not used at all.                *
!                     The land/sea mask is used everywhere.                    *
!                                                                              *
!    Feb 2012 J. Brioude: modified to fortran 90                               *
!                                                                              *
!    2015-03-25,  A. Dingwell: Fixed indentation (FLEXPART coding standard)    *
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
! xlanduse          fractions of numclass landuses for each model grid point   *
! xlandinvent       landuse inventory (1/6 deg resolution)                     *
!                                                                              *
!*******************************************************************************

  use par_mod
  use com_mod

  implicit none
  integer :: ix,jy,k,l,li,nrefine,iix,jjy,n,nxn2,nyn2
  integer,parameter :: lumaxx=1200,lumaxy=600
  real,parameter :: xlon0lu=-180.,ylat0lu=-90.
  real,parameter :: dxlu=0.3
  real :: xlon,ylat,sumperc,p,xi,yj
  real :: xlandusep(lumaxx,lumaxy,numclass)
!  real :: xlanduse(0:nxmax-1,0:nymax-1,numclass)
!  real :: xlandusen(0:nxmaxn-1,0:nymaxn-1,numclass,maxnests)

!      include 'includepar'
!      include 'includecom'
!
!      integer ix,jy,i,j,k,n,l
!      real x,y,xlon,ylat

  do ix=1,lumaxx
    do jy=1,lumaxy
      do k=1,numclass
        xlandusep(ix,jy,k)=0.
      end do
      sumperc=0.
      do li=1,3
        sumperc=sumperc+landinvent(ix,jy,li+3)
      end do
      do li=1,3
        k=landinvent(ix,jy,li)
        if (sumperc.gt.0) then
          p=landinvent(ix,jy,li+3)/sumperc
        else
          p=0
        endif
        ! p has values between 0 and 1
        xlandusep(ix,jy,k)=p
      end do
    end do
  end do

  nrefine=10

  do ix=0,nxmin1
    do jy=0,nymin1
    ! FLEXPART_WRF - use this routine to get lat,lon
      do k=1,numclass
        sumperc=0.
        xlanduse(ix,jy,k)=0.
      enddo
      do iix=1, 1
        do jjy=1, 1
          call xyindex_to_ll_wrf( 0, real(ix), real(jy), xlon, ylat )
          
          if (xlon.ge.(xlon0lu+lumaxx*dxlu))  then
               xlon=xlon-lumaxx*dxlu
          endif
          xi=int((xlon-xlon0lu)/dxlu)+1
          yj=int((ylat-ylat0lu)/dxlu)+1
          if (xi.gt.lumaxx) xi=xi-lumaxx
          if (yj.gt.lumaxy) yj=yj-lumaxy
          if (xi.lt.0) then
            write (*,*) 'problem with landuseinv sampling: ', &
              xlon,xlon0lu,ix,iix,xlon0,dx,nxmax
            stop
          endif
          do k=1,numclass
            xlanduse(ix,jy,k)= &
              xlanduse(ix,jy,k)+xlandusep(int(xi),int(yj),k)
            sumperc=sumperc+xlanduse(ix,jy,k)  ! just for the check if landuseinv. is available
          end do
        end do
      end do
      if (sumperc.gt.0) then                       ! detailed landuse available
        sumperc=0.
        do k=1,numclass
          xlanduse(ix,jy,k)= &
            xlanduse(ix,jy,k)/real(nrefine*nrefine)
          sumperc=sumperc+xlanduse(ix,jy,k)
        end do
        !cc the sum of all categories should be 1 ... 100 percent ... in order to get vdep right!
        if (sumperc.lt.1-1E-5) then
          do k=1,numclass
            xlanduse(ix,jy,k)= &
              xlanduse(ix,jy,k)/sumperc
          end do
        endif
      else
        if (lsm(ix,jy).lt.0.1) then           ! over sea  -> ocean
          xlanduse(ix,jy,3)=1.
        else                                  ! over land -> rangeland
          xlanduse(ix,jy,7)=1.
        endif
      endif
    end do
  end do


!   open(4,file='landusetest',form='formatted')
!   do 56 k=1,13
 !  do 55 ix=0,nxmin1
!  55    write (4,*) (xlanduse(ix,jy,k),jy=0,nymin1)
!  56    continue
!   close(4)
!   write (*,*) 'landuse written'
 ! stop
 !   open(4,file='landseatest'//ck,form='formatted')
 !  do 57 ix=0,nxmin1
!  57       write (4,*) (lsm(ix,jy),jy=0,nymin1)
!    write (*,*) 'landseamask written'

  do l=1,numbnests
    do jy=0,nyn(l)-1
      do ix=0,nxn(l)-1

        do k=1,numclass
          xlandusen(ix,jy,k,l)=0.0
          sumperc=0.
        enddo

        do iix=1, 1
          do jjy=1, 1
!           xlon=(ix+(iix-1)/real(nrefine))*dx+xlon0        ! longitude, should be between -180 and 179
!           ylat=(jy+(jjy-1)/real(nrefine))*dy+ylat0        ! and lat. of each gridpoint
            call xyindex_to_ll_wrf( l, real(ix), real(jy), xlon, ylat )

            xi=int((xlon-xlon0lu)/dxlu)+1
            yj=int((ylat-ylat0lu)/dxlu)+1
            if (xi.gt.lumaxx) xi=xi-lumaxx
            if (yj.gt.lumaxy) yj=yj-lumaxy
            do k=1,numclass
              xlandusen(ix,jy,k,l)=xlandusen(ix,jy,k,l)+ &
                xlandusep(int(xi),int(yj),k)
              sumperc=sumperc+xlandusen(ix,jy,k,l)
            end do
          end do
        end do
        if (sumperc.gt.0) then                     ! detailed landuse available
          sumperc=0.
          do k=1,numclass
            xlandusen(ix,jy,k,l)= &
              xlandusen(ix,jy,k,l)/real(nrefine*nrefine)
            sumperc=sumperc+xlandusen(ix,jy,k,l)
          end do
  !cc the sum of all categories should be 1 ... 100 percent ... in order to get vdep right!
          if (sumperc.lt.1-1E-5) then
            do k=1,numclass
              xlandusen(ix,jy,k,l)=xlandusen(ix,jy,k,l)/sumperc
            end do
          endif
        else                                    ! check land/sea mask
          if (lsmn(ix,jy,l).lt.0.1) then   ! over sea  -> ocean
            xlandusen(ix,jy,3,l)=1.
          else                        ! over land -> grasslands
            xlandusen(ix,jy,7,l)=1.
          endif
        endif
      end do
    end do
  end do

end subroutine assignland
