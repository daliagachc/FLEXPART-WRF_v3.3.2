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
      subroutine outgrid_init_irreg
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
!                                                                              *
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

  use flux_mod
  use oh_mod
  use unc_mod
  use outg_mod
  use par_mod
  use com_mod

!      include 'includepar'
!      include 'includecom'
     implicit none 
      integer :: ix,jy,kz,k,i,nage,l,iix,jjy,ixp,jyp,i1,j1,j,ngrid
!     real ylat,gridarea,ylatp,ylatm,hzone,cosfact,cosfactm,cosfactp
      real :: ymet,gridarea,xl1,xl2,yl1,yl2,m1,m2,tmpx,tmpy
      real :: xmet,xl,yl,ddx,ddy,rddx,rddy,p1,p2,p3,p4,xtn,ytn,oroh
  integer :: ks,kp,stat,ix2,jy2
  real,parameter :: eps=nxmax/3.e5
      real :: lon2(4),lat2(4)
      real ( kind = 8 ) :: sphere01_polygon_area,haversine,area1



! Compute surface area and volume of each grid cell: area, volume;
! and the areas of the northward and eastward facing walls: areaeast, areanorth
!***********************************************************************
!     do jy=0,10
!     print*,areamet(1,jy),areamet(2,jy),areamet2(2,jy)
!     enddo
      do jy=0,numygrid-1

!        ylat=outlat0+(real(jy)+0.5)*dyout
!        ylatp=ylat+0.5*dyout
!        ylatm=ylat-0.5*dyout
!        if ((ylatm.lt.0).and.(ylatp.gt.0.)) then
!          hzone=dyout*r_earth*pi180
!        else
!
!C Calculate area of grid cell with formula M=2*pi*R*h*dx/360,
!C see Netz, Formeln der Mathematik, 5. Auflage (1983), p.90
!*************************************************************
!
!          cosfact=cos(ylat*pi180)*r_earth
!          cosfactp=cos(ylatp*pi180)*r_earth
!          cosfactm=cos(ylatm*pi180)*r_earth
!          if (cosfactp.lt.cosfactm) then
!            hzone=sqrt(r_earth**2-cosfactp**2)-
!     +      sqrt(r_earth**2-cosfactm**2)
!          else
!            hzone=sqrt(r_earth**2-cosfactm**2)-
!     +      sqrt(r_earth**2-cosfactp**2)
!          endif
!        endif
!
!C Surface are of a grid cell at a latitude ylat
!***********************************************
!
!        gridarea=2.*pi*r_earth*hzone*dxout/360.

! for FLEXPART_WRF, dx & dy are in meters, and no cos(lat) is needed
! ??? maybe should incorporate map factor here,
!     and also for areaeast & areanorth ???

        do ix=0,numxgrid-1
!        gridarea=dxout*dyout ! what is needed is the true area based on map factors
!JB: find a way to get the area between 2 output grid cell using areamet
!        xl1=(real(ix)*dxout+out_xm0)/dx
!        yl1=(real(jy)*dyout+out_ym0)/dy
!        xl2=(real(ix+1)*dxout+out_xm0)/dx
!        yl2=(real(jy+1)*dyout+out_ym0)/dy
!!      xr=out_xm0+real(numxgrid)*dxout
!        m1=0.5*(m_x(int(xl1),int(yl1),1)+m_x(int(xl2),int(yl1),1))
!        m2=0.5*(m_y(int(xl1),int(yl1),1)+m_x(int(xl1),int(yl2),1))
!      area(ix,jy)=dxout*m1*dyout*m2

! A more precise method
          tmpx=out_xm0+(float(ix))*dxout
          tmpy=out_ym0+(float(jy))*dyout
          call xymeter_to_ll_wrf_out(tmpx,tmpy,lon2(1),lat2(1))
          tmpx=out_xm0+(float(ix+1))*dxout
          tmpy=out_ym0+(float(jy))*dyout
          call xymeter_to_ll_wrf_out(tmpx,tmpy,lon2(2),lat2(2))
          tmpx=out_xm0+(float(ix+1))*dxout
          tmpy=out_ym0+(float(jy+1))*dyout
          call xymeter_to_ll_wrf_out(tmpx,tmpy,lon2(3),lat2(3))
          tmpx=out_xm0+(float(ix))*dxout
          tmpy=out_ym0+(float(jy+1))*dyout
          call xymeter_to_ll_wrf_out(tmpx,tmpy,lon2(4),lat2(4))
        area1=sphere01_polygon_area ( 4, real(lat2,kind=8), real(lon2,kind=8) )
        area(ix,jy)=real(area1)*6370000.*6370000./coefdx/coefdx
!
! Volume = area x box height
!***************************

          volume(ix,jy,1)=area(ix,jy)*outheight(1)

! for FLEXPART_WRF, dx & dy are in meters, and no cos(lat) is needed
!          areaeast(ix,jy,1)=dyout*r_earth*pi180*outheight(1)
!          areanorth(ix,jy,1)=cos(ylat*pi180)*dxout*r_earth*pi180*
!     +    outheight(1)
          areaeast(ix,jy,1)=dyout*outheight(1)
          areanorth(ix,jy,1)=dxout*outheight(1)

          do kz=2,numzgrid

!            areaeast(ix,jy,kz)=dyout*r_earth*pi180*
!     +      (outheight(kz)-outheight(kz-1))
!            areanorth(ix,jy,kz)=cos(ylat*pi180)*dxout*r_earth*pi180*
!     +      (outheight(kz)-outheight(kz-1))
            areaeast(ix,jy,kz)=dyout*(outheight(kz)-outheight(kz-1))
            areanorth(ix,jy,kz)=dxout*(outheight(kz)-outheight(kz-1))

          volume(ix,jy,kz)=area(ix,jy)*(outheight(kz)-outheight(kz-1))
      end do
    end do
  end do




!******************************************************************
! Determine average height of model topography in output grid cells
!******************************************************************

! Loop over all output grid cells
!********************************

      do jjy=0,numygrid-1
        do iix=0,numxgrid-1
          oroh=0.

! Take 100 samples of the topography in every grid cell
!******************************************************

          do j1=1,10
! for FLEXPART_WRF, x & y coords are in meters,
! and the lon & lat variables below are in meters.
            ymet=out_ym0+(real(jjy)+real(j1)/10.-0.05)*dyout
            yl=(ymet-ymet0)/dy
            do i1=1,10
              xmet=out_xm0+(real(iix)+real(i1)/10.-0.05)*dxout
              xl=(xmet-xmet0)/dx

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

          oroout(iix,jjy)=oroh/100.
    end do
  end do


  ! if necessary allocate flux fields
  if (iflux.eq.1) then
    allocate(flux(6,0:numxgrid-1,0:numygrid-1,numzgrid, &
         1:nspec,1:maxpointspec_act,1:nageclass),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate flux array '
  endif
  
  !write (*,*) 'allocating: in a sec',OHREA
  if (OHREA.eqv..TRUE.) then
  !   write (*,*) 'allocating: ',maxxOH,maxyOH,maxzOH
    allocate(OH_field(12,0:maxxOH-1,0:maxyOH-1,maxzOH) &
         ,stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate OH array '
    allocate(OH_field_height(7) &
         ,stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate OH array '
  endif
  ! gridunc,griduncn        uncertainty of outputted concentrations
!  print*,'gridunc allocated'
  allocate(gridunc(0:numxgrid-1,0:numygrid-1,numzgrid,maxspec, &
       maxpointspec_act,nclassunc,maxageclass),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'
  if (ldirect.gt.0) then
  allocate(wetgridunc(0:numxgrid-1,0:numygrid-1,maxspec, &
       maxpointspec_act,nclassunc,maxageclass),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'
  allocate(drygridunc(0:numxgrid-1,0:numygrid-1,maxspec, &
       maxpointspec_act,nclassunc,maxageclass),stat=stat)
  allocate(drygridunc2(0:numxgrid-1,0:numygrid-1,maxspec, &
       maxpointspec_act,nclassunc,maxageclass),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'
  endif
  !write (*,*) 'Dimensions for fields', numxgrid,numygrid, &
  !     maxspec,maxpointspec_act,nclassunc,maxageclass

!  print*,'alloc gridunc',numxgrid-1,numygrid-1,numzgrid,maxspec, &
!         maxpointspec_act,nclassunc,maxageclass


  write (*,*) ' Allocating fields for nested and global output (x,y): ', &
       max(numxgrid,numxgridn),max(numygrid,numygridn)

  ! allocate fields for concoutput with maximum dimension of outgrid
  ! and outgrid_nest
  allocate(gridsigma(0:max(numxgrid,numxgridn)-1, &
       0:max(numygrid,numygridn)-1,numzgrid),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'
  allocate(grid(0:max(numxgrid,numxgridn)-1, &
       0:max(numygrid,numygridn)-1,numzgrid),stat=stat)
  allocate(grid2(0:max(numxgrid,numxgridn)-1, &
       0:max(numygrid,numygridn)-1,numzgrid,maxpointspec_act,nageclass),stat=stat)
  allocate(grid3(0:max(numxgrid,numxgridn)-1, &
       0:max(numygrid,numygridn)-1,numzgrid,maxpointspec_act,nageclass),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'
  allocate(densityoutgrid(0:max(numxgrid,numxgridn)-1, &
       0:max(numygrid,numygridn)-1,numzgrid),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'

  allocate(factor3d(0:max(numxgrid,numxgridn)-1, &
       0:max(numygrid,numygridn)-1,numzgrid),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'
  allocate(sparse_dump_r(max(numxgrid,numxgridn)* &
       max(numygrid,numygridn)*numzgrid),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'
  allocate(sparse_dump_i(max(numxgrid,numxgridn)* &
       max(numygrid,numygridn)*numzgrid),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'

  ! deposition fields are only allocated for forward runs
  if (ldirect.gt.0) then
     allocate(wetgridsigma(0:max(numxgrid,numxgridn)-1, &
          0:max(numygrid,numygridn)-1),stat=stat)
     if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'
     allocate(drygridsigma(0:max(numxgrid,numxgridn)-1, &
          0:max(numygrid,numygridn)-1),stat=stat)
     if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'
     allocate(wetgrid(0:max(numxgrid,numxgridn)-1, &
          0:max(numygrid,numygridn)-1),stat=stat)
     allocate(wetgrid2(0:max(numxgrid,numxgridn)-1, &
          0:max(numygrid,numygridn)-1,maxpointspec_act,nageclass),stat=stat)
     if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'
     allocate(drygrid(0:max(numxgrid,numxgridn)-1, &
          0:max(numygrid,numygridn)-1),stat=stat)
     allocate(drygrid2(0:max(numxgrid,numxgridn)-1, &
          0:max(numygrid,numygridn)-1,maxpointspec_act,nageclass),stat=stat)
     if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'
  endif

  ! Initial condition field

  if (linit_cond.gt.0) then
    allocate(init_cond(0:numxgrid-1,0:numygrid-1,numzgrid,maxspec, &
         maxpointspec_act),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate init_cond'
  endif

  !************************
  ! Initialize output grids
  !************************

  do ks=1,nspec
  do kp=1,maxpointspec_act
    do i=1,numreceptor
  ! Receptor points
      creceptor(i,ks)=0.
    end do
    do nage=1,nageclass
      do jy=0,numygrid-1
        do ix=0,numxgrid-1
          do l=1,nclassunc
  ! Deposition fields
            if (ldirect.gt.0) then
            wetgridunc(ix,jy,ks,kp,l,nage)=0.
            drygridunc(ix,jy,ks,kp,l,nage)=0.
!            drygridunc2(ix,jy,ks,kp,l,nage)=0.
            endif
            do kz=1,numzgrid
              if (iflux.eq.1) then
  ! Flux fields
                 do i=1,5
                   flux(i,ix,jy,kz,ks,kp,nage)=0.
                 end do
              endif
  ! Initial condition field
              if ((l.eq.1).and.(nage.eq.1).and.(linit_cond.gt.0)) &
                   init_cond(ix,jy,kz,ks,kp)=0.
  ! Concentration fields
              gridunc(ix,jy,kz,ks,kp,l,nage)=0.
            end do
          end do
        end do
      end do
    end do
  end do
  end do



end subroutine outgrid_init_irreg

function sphere01_polygon_area ( n, lat, lon )

!*****************************************************************************80
!
!! SPHERE01_POLYGON_AREA returns the area of a spherical polygon.
!
!  Discussion:
!
!    On a unit sphere, the area of a spherical polygon with N sides
!    is equal to the spherical excess:
!
!      E = sum ( interior angles ) - ( N - 2 ) * pi.
!
!    On a sphere with radius R, the area is the spherical excess multiplied
!    by R * R.
!
!    The code was revised in accordance with suggestions in Carvalho and
!    Cavalcanti.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 August 2005
!
!  Author:
!
!    Original C version by Robert Miller.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paulo Cezar Pinto Carvalho, Paulo Roma Cavalcanti,
!    Point in Polyhedron Testing Using Spherical Polygons,
!    in Graphics Gems V,
!    edited by Alan Paeth,
!    Academic Press, 1995,
!    ISBN: 0125434553,
!    LC: T385.G6975.
!
!    Robert Miller,
!    Computing the Area of a Spherical Polygon,
!    Graphics Gems, Volume IV, pages 132-138,
!    Edited by Paul Heckbert,
!    Academic Press, 1994, T385.G6974.
!
!    Eric Weisstein,
!    "Spherical Polygon",
!    CRC Concise Encyclopedia of Mathematics,
!    CRC Press, 1999.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of vertices.
!
!    Input, real ( kind = 8 ) LAT[N], LON[N], the latitudes and longitudes 
!    of the vertices of the spherical polygon.
!
!    Output, real ( kind = 8 ) SPHERE01_POLYGON_AREA, the area of the 
!    spherical polygon, measured in spherical radians.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) area
  real ( kind = 8 ) b
  real ( kind = 8 ) beta1
  real ( kind = 8 ) beta2
  real ( kind = 8 ) c
  real ( kind = 8 ) cos_b1
  real ( kind = 8 ) cos_b2
  real ( kind = 8 ) excess
  real ( kind = 8 ) hav_a
  real ( kind = 8 ) haversine
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) lam
  real ( kind = 8 ) lam1
  real ( kind = 8 ) lam2
  real ( kind = 8 ) lat(n)
  real ( kind = 8 ) lon(n)
  real ( kind = 8 ), parameter :: pi_half = 1.5707963267948966192313D+00
  real ( kind = 8 ) s
  real ( kind = 8 ) sphere01_polygon_area
  real ( kind = 8 ) t
  real ( kind = 8 ),parameter :: degrees_to_radians=3.141592653589793D+00 / 180.0D+00

  area = 0.0D+00

  do j=1,n
  lon(j)=lon(j)*degrees_to_radians
  lat(j)=lat(j)*degrees_to_radians
  enddo
  do j = 1, n + 1
! do j = 1, n 

    if ( j == 1 ) then
      lam1 = lon(j)
      beta1 = lat(j)
      lam2 = lon(j+1)
      beta2 = lat(j+1)
      cos_b1 = cos ( beta1 )
      cos_b2 = cos ( beta2 )
    else
!      k = mod ( j + 1, n + 1 )
!    k = mod ( j , n )
      k=j
      if (j.gt.n) k=1
      lam1 = lam2
      beta1 = beta2
      lam2 = lon(k)
      beta2 = lat(k)
! print*,'sphere',n,k,lon(k),lat(k)
      cos_b1 = cos_b2
      cos_b2 = cos ( beta2 )
    end if

    if ( lam1 /= lam2 ) then

      hav_a = haversine ( beta2 - beta1 ) &
        + cos_b1 * cos_b2 * haversine ( lam2 - lam1 )
      a = 2.0D+00 * asin ( sqrt ( hav_a ) )

      b = pi_half - beta2
      c = pi_half - beta1
      s = 0.5D+00 * ( a + b + c )
!
!  Given the three sides of a spherical triangle, we can use a formula
!  to find the spherical excess.
!
      t = tan ( s / 2.0D+00 ) * tan ( ( s - a ) / 2.0D+00 ) &
        * tan ( ( s - b ) / 2.0D+00 ) * tan ( ( s - c ) / 2.0D+00 )

      excess = abs ( 4.0D+00 * atan ( sqrt ( abs ( t ) ) ) )

      if ( lam1 < lam2 ) then
        lam = lam2 - lam1
      else
        lam = lam2 - lam1 + 4.0D+00 * pi_half
      end if

      if ( 2.0D+00 * pi_half < lam ) then
        excess = -excess 
      end if

      area = area + excess

    end if

  end do

  sphere01_polygon_area = abs ( area )

  return
end

function haversine ( a )

!*****************************************************************************80
!
!! HAVERSINE computes the haversine of an angle.
!
!  Discussion:
!
!    haversine(A) = ( 1 - cos ( A ) ) / 2
!
!    The haversine is useful in spherical trigonometry.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the angle.
!
!    Output, real ( kind = 8 ) HAVERSINE, the haversine of the angle.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) haversine

  haversine = ( 1.0D+00 - cos ( a ) ) / 2.0D+00

  return
end
