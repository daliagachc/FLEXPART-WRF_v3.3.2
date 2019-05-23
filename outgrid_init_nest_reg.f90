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

      subroutine outgrid_init_nest_reg
!*******************************************************************************
!                                                                              *
!     Note:  This is the FLEXPART_WRF version of subroutine outgrid_init.      *
!            The computational grid is the WRF x-y grid rather than lat-lon.   *
!                                                                              *
!  This routine calculates, for each grid cell of the output grid, the         *
!  volume, the surface area, and the areas of the northward and eastward       *
!  facing surfaces.                                                            *
!                                                                              *
!     Author: A. Stohl                                                         *
!                                                                              *
!     7 August 2002                                                            *
!                                                                              *
!    26 Oct 2005, R. Easter - changes in gridarea, areaeast, areanorth         *
!                             associated with WRF horizontal grid.             *
!     Dec 2005, R. Easter - changed names of "*lon0*" & "*lat0*" variables     *
!     July 2012, J Brioude - modified for regular output grid.                 *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
!                                                                              *
! area               surface area of all output grid cells                     *
! areaeast           eastward facing wall area of all output grid cells        *
! areanorth          northward facing wall area of all output grid cells       *
! volume             volumes of all output grid cells                          *
!                                                                              *
!*******************************************************************************

  use unc_mod
  use outg_mod
  use par_mod
  use com_mod

!      include 'includepar'
!      include 'includecom'
     implicit none 
      integer :: ix,jy,kz,k,i,nage,l,iix,jjy,ixp,jyp,i1,j1,j,ngrid
      real :: ylat,gridarea,ylatp,ylatm,hzone,cosfact,cosfactm,cosfactp
      real :: ymet,xlon
      real :: xmet,xl,yl,ddx,ddy,rddx,rddy,p1,p2,p3,p4,xtn,ytn,oroh
  integer :: ks,kp,stat
  real,parameter :: eps=nxmax/3.e5



! Compute surface area and volume of each grid cell: area, volume;
! and the areas of the northward and eastward facing walls: areaeast, areanorth
!***********************************************************************

      do jy=0,numygridn-1

        ylat=outlat0n+(real(jy)+0.5)*dyoutln
        ylatp=ylat+0.5*dyoutln
        ylatm=ylat-0.5*dyoutln
        if ((ylatm.lt.0).and.(ylatp.gt.0.)) then
          hzone=dyoutln*r_earth*pi180
        else

!C Calculate area of grid cell with formula M=2*pi*R*h*dx/360,
!C see Netz, Formeln der Mathematik, 5. Auflage (1983), p.90
!*************************************************************
!
          cosfact=cos(ylat*pi180)*r_earth
          cosfactp=cos(ylatp*pi180)*r_earth
          cosfactm=cos(ylatm*pi180)*r_earth
          if (cosfactp.lt.cosfactm) then
            hzone=sqrt(r_earth**2-cosfactp**2)- &
            sqrt(r_earth**2-cosfactm**2)
          else
            hzone=sqrt(r_earth**2-cosfactm**2)- &
            sqrt(r_earth**2-cosfactp**2)
          endif
        endif

!C Surface are of a grid cell at a latitude ylat
!***********************************************
!
         gridarea=2.*pi*r_earth*hzone*dxoutln/360.

! for FLEXPART_WRF, dx & dy are in meters, and no cos(lat) is needed
! ??? maybe should incorporate map factor here,
!     and also for areaeast & areanorth ???
!       gridarea=dxoutn*dyoutn

        do ix=0,numxgridn-1
          arean(ix,jy)=gridarea

! Volume = area x box height
!***************************

          volumen(ix,jy,1)=arean(ix,jy)*outheight(1)

          do kz=2,numzgrid

          volumen(ix,jy,kz)=arean(ix,jy)*(outheight(kz)-outheight(kz-1))
      end do
    end do
  end do




!******************************************************************
! Determine average height of model topography in output grid cells
!******************************************************************

! Loop over all output grid cells
!********************************

      do jjy=0,numygridn-1
        do iix=0,numxgridn-1
          oroh=0.

! Take 100 samples of the topography in every grid cell
!******************************************************

          do j1=1,10
! for FLEXPART_WRF, x & y coords are in meters,
! and the lon & lat variables below are in meters.
            ylat=outlat0n+(float(jjy)+float(j1)/10.-0.05)*dyoutln !in degrees
!           ymet=out_ym0n+(real(jjy)+real(j1)/10.-0.05)*dyoutn
!           yl=(ymet-ymet0)/dy
            do i1=1,10
            xlon=outlon0n+(float(iix)+float(i1)/10.-0.05)*dxoutln !in degrees
!              xmet=out_xm0n+(real(iix)+real(i1)/10.-0.05)*dxoutn
        call ll_to_xymeter_wrf(xlon,ylat,xmet,ymet)
           yl=(ymet-ymet0)/dy
           xl=(xmet-xmet0)/dx
!             xl=(xmet-xmet0)/dx

! Determine the nest we are in
!*****************************

              ngrid=0
              do j=numbnests,1,-1
                if ((xl.gt.xln(j)).and.(xl.lt.xrn(j)).and. &
                (yl.gt.yln(j)).and.(yl.lt.yrn(j))) then
                  ngrid=j
                  goto 43
                endif
          end do
43            continue

! Determine (nested) grid coordinates and auxiliary parameters used for interpolation
!************************************************************************************

              if (ngrid.gt.0) then
                xtn=(xl-xln(ngrid))*xresoln(ngrid)
                ytn=(yl-yln(ngrid))*yresoln(ngrid)
                ix=int(xtn)
                jy=int(ytn)
                ddy=ytn-real(jy)
                ddx=xtn-real(ix)
              else
                ix=int(xl)
                jy=int(yl)
                ddy=yl-real(jy)
                ddx=xl-real(ix)
              endif
              ixp=ix+1
              jyp=jy+1
              rddx=1.-ddx
              rddy=1.-ddy
              p1=rddx*rddy
              p2=ddx*rddy
              p3=rddx*ddy
              p4=ddx*ddy

          if (ngrid.gt.0) then
            oroh=oroh+p1*oron(ix ,jy ,ngrid) &
                 + p2*oron(ixp,jy ,ngrid) &
                 + p3*oron(ix ,jyp,ngrid) &
                 + p4*oron(ixp,jyp,ngrid)
          else
            oroh=oroh+p1*oro(ix ,jy) &
                 + p2*oro(ixp,jy) &
                 + p3*oro(ix ,jyp) &
                 + p4*oro(ixp,jyp)
          endif
        end do
      end do

! Divide by the number of samples taken
!**************************************

          orooutn(iix,jjy)=oroh/100.
    end do
  end do


  ! gridunc,griduncn        uncertainty of outputted concentrations
  allocate(griduncn(0:numxgridn-1,0:numygridn-1,numzgrid,maxspec, &
       maxpointspec_act,nclassunc,maxageclass),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR:could not allocate nested gridunc'

  if (ldirect.gt.0) then
  allocate(wetgriduncn(0:numxgridn-1,0:numygridn-1,maxspec, &
       maxpointspec_act,nclassunc,maxageclass),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR:could not allocate nested gridunc'
  allocate(drygriduncn(0:numxgridn-1,0:numygridn-1,maxspec, &
       maxpointspec_act,nclassunc,maxageclass),stat=stat)
  allocate(drygriduncn2(0:numxgridn-1,0:numygridn-1,maxspec, &
       maxpointspec_act,nclassunc,maxageclass),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR:could not allocate nested gridunc'
  endif

  !write (*,*) 'Dimensions for fields', numxgrid,numygrid, &
  !     maxspec,maxpointspec_act,nclassunc,maxageclass

  ! allocate fields for concoutput with maximum dimension of outgrid
  ! and outgrid_nest
  ! Initial condition field

  !************************
  ! Initialize output grids
  !************************

  do kp=1,maxpointspec_act
  do ks=1,nspec
    do nage=1,nageclass
      do jy=0,numygridn-1
        do ix=0,numxgridn-1
          do l=1,nclassunc
  ! Deposition fields
            if (ldirect.gt.0) then
              wetgriduncn(ix,jy,ks,kp,l,nage)=0.
              drygriduncn(ix,jy,ks,kp,l,nage)=0.
            endif
  ! Concentration fields
            do kz=1,numzgrid
              griduncn(ix,jy,kz,ks,kp,l,nage)=0.
            end do
          end do
        end do
      end do
    end do
  end do
  end do


end subroutine outgrid_init_nest_reg

