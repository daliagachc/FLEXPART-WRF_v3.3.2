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

      subroutine outgrid_init_reg()
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
!    July 2012, J. Brioude - modified for regular output grid                  *
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
      real :: ylat,gridarea,ylatp,ylatm,hzone,cosfact,cosfactm,cosfactp
      real :: ymet,xlon    
      real :: xmet,xl,yl,ddx,ddy,rddx,rddy,p1,p2,p3,p4,xtn,ytn,oroh
  integer :: ks,kp,stat
  real,parameter :: eps=nxmax/3.e5



! Compute surface area and volume of each grid cell: area, volume;
! and the areas of the northward and eastward facing walls: areaeast, areanorth
!***********************************************************************

      do jy=0,numygrid-1

        ylat=outlat0+(float(jy)+0.5)*dyoutl
        ylatp=ylat+0.5*dyoutl              
        ylatm=ylat-0.5*dyoutl              
        if ((ylatm.lt.0).and.(ylatp.gt.0.)) then    
          hzone=dyoutl*r_earth*pi180          
        else                 
          cosfact=cos(ylat*pi180)*r_earth          
          cosfactp=cos(ylatp*pi180)*r_earth       
          cosfactm=cos(ylatm*pi180)*r_earth      
          if (cosfactp.lt.cosfactm) then         
            hzone=sqrt(r_earth**2-cosfactp**2)-  &
            sqrt(r_earth**2-cosfactm**2)
          else
            hzone=sqrt(r_earth**2-cosfactm**2)- &
            sqrt(r_earth**2-cosfactp**2)
          endif
        endif
!
!C Surface are of a grid cell at a latitude ylat
!***********************************************
!
         gridarea=2.*pi*r_earth*hzone*dxoutl/360.

! for FLEXPART_WRF, dx & dy are in meters, and no cos(lat) is needed
! ??? maybe should incorporate map factor here,
!     and also for areaeast & areanorth ???
!       gridarea=dxout*dyout

        do ix=0,numxgrid-1
          area(ix,jy)=gridarea
! Volume = area x box height
!***************************

          volume(ix,jy,1)=area(ix,jy)*outheight(1)

! for FLEXPART_WRF, dx & dy are in meters, and no cos(lat) is needed
          areaeast(ix,jy,1)=dyoutl*r_earth*pi180*outheight(1)
          areanorth(ix,jy,1)=cos(ylat*pi180)*dxoutl*r_earth*pi180* &
          outheight(1)
!         areaeast(ix,jy,1)=dyout*outheight(1)
!         areanorth(ix,jy,1)=dxout*outheight(1)

          do kz=2,numzgrid

             areaeast(ix,jy,kz)=dyoutl*r_earth*pi180* &
             (outheight(kz)-outheight(kz-1))
             areanorth(ix,jy,kz)=cos(ylat*pi180)*dxoutl*r_earth*pi180* &
             (outheight(kz)-outheight(kz-1))
!           areaeast(ix,jy,kz)=dyout*(outheight(kz)-outheight(kz-1))
!           areanorth(ix,jy,kz)=dxout*(outheight(kz)-outheight(kz-1))

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
            ylat=outlat0+(float(jjy)+float(j1)/10.-0.05)*dyoutl !in degrees
!            ymet=out_ym0+(real(jjy)+real(j1)/10.-0.05)*dyout
!            yl=(ymet-ymet0)/dy
            do i1=1,10
            xlon=outlon0+(float(iix)+float(i1)/10.-0.05)*dxoutl !in degrees
        call ll_to_xymeter_wrf(xlon,ylat,xmet,ymet)
           yl=(ymet-ymet0)/dy
           xl=(xmet-xmet0)/dx

!              xmet=out_xm0+(real(iix)+real(i1)/10.-0.05)*dxout
!              xl=(xmet-xmet0)/dx

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
  allocate(gridunc(0:numxgrid-1,0:numygrid-1,numzgrid,maxspec, &
       maxpointspec_act,nclassunc,maxageclass),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'

   print*,'alloc gridunc',numxgrid-1,numygrid-1,numzgrid,maxspec, &
          maxpointspec_act,nclassunc,maxageclass
  if (ldirect.gt.0) then
  allocate(wetgridunc(0:numxgrid-1,0:numygrid-1,maxspec, &
       maxpointspec_act,nclassunc,maxageclass),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'
  allocate(drygridunc(0:numxgrid-1,0:numygrid-1,maxspec, &
       maxpointspec_act,nclassunc,maxageclass),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'
  allocate(drygridunc2(0:numxgrid-1,0:numygrid-1,maxspec, &
       maxpointspec_act,nclassunc,maxageclass),stat=stat)
  endif
  !write (*,*) 'Dimensions for fields', numxgrid,numygrid, &
  !     maxspec,maxpointspec_act,nclassunc,maxageclass

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



end subroutine outgrid_init_reg

