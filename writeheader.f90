!***********************************************************************
!* Copyright 2012,2013                                                *
!* Jerome Brioude, Delia Arnold, Andreas Stohl, Wayne Angevine,       *
!* John Burkhart, Massimo Cassiani, Adam Dingwell, Richard C Easter, Sabine Eckhardt,*
!* Stephanie Evan, Jerome D Fast, Don Morton, Ignacio Pisso,          *
!* Petra Seibert, Gerard Wotawa, Caroline Forster, Harald Sodemann,   *
!*                                                                     *
!* This file is part of FLEXPART WRF            
!                                                                     *
! FLEXPART is free software: you can redistribute it and/or modify    *
! it under the terms of the GNU General Public License as published by*
! the Free Software Foundation, either version 3 of the License, or   *
! (at your option) any later version.                                 *
!                                                                     *
! FLEXPART is distributed in the hope that it will be useful,         *
! but WITHOUT ANY WARRANTY; without even the implied warranty of      *
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       *
! GNU General Public License for more details.                        *
!                                                                     *
! You should have received a copy of the GNU General Public License   *
! along with FLEXPART.  If not, see <http://www.gnu.org/licenses/>.   *
!**********************************************************************

subroutine writeheader

  !*****************************************************************************
  !                                                                            *
  !  This routine produces a file header containing basic information on the   *
  !  settings of the FLEXPART run.                                             *
  !  The header file is essential and must be read in by any postprocessing    *
  !  program before reading in the output data.                                *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     7 August 2002                                                          *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! xlon                   longitude                                           *
  ! xl                     model x coordinate                                  *
  ! ylat                   latitude                                            *
  ! yl                     model y coordinate                                  *
  !                                                                            *
  !*****************************************************************************

  use point_mod
  use outg_mod
  use par_mod
  use com_mod

  implicit none

  integer :: jjjjmmdd,ihmmss,i,ix,jy,j
  real :: xp1,yp1,xp2,yp2
  real :: xsw,xne,ysw,yne,tmpx,tmpy,tmplon,tmplat,xl2,yl2


  !************************
  ! Open header output file
  !************************

  open(unitheader,file=path(1)(1:length(1))//'header', &
       form='unformatted',err=998)


  ! Write the header information
  !*****************************

  if (ldirect.eq.1) then
!   write(unitheader) ibdate,ibtime,'FLEXWRF  V2.1'
  if (outgrid_option .eq. 1) then
    write(unitheader) ibdate,ibtime,'FLEXWRF lalo '
  else
  if (map_proj_id.eq.1) write(unitheader) ibdate,ibtime,'FLEXWRF lamb '
  if (map_proj_id.eq.2) write(unitheader) ibdate,ibtime,'FLEXWRF ster '
  if (map_proj_id.eq.3) write(unitheader) ibdate,ibtime,'FLEXWRF merc '
  if (map_proj_id.eq.4) write(unitheader) ibdate,ibtime,'FLEXWRF glob '
  endif
  else
  if (outgrid_option .eq. 1) then
    write(unitheader) iedate,ietime,'FLEXWRF lalo '
  else
  if (map_proj_id.eq.1) write(unitheader) iedate,ietime,'FLEXWRF lamb '
  if (map_proj_id.eq.2) write(unitheader) iedate,ietime,'FLEXWRF ster '
  if (map_proj_id.eq.3) write(unitheader) iedate,ietime,'FLEXWRF merc '
  if (map_proj_id.eq.4) write(unitheader) iedate,ietime,'FLEXWRF glob '
  endif
  endif

  ! Write info on output interval, averaging time, sampling time
  !*************************************************************

  write(unitheader) loutstep,loutaver,loutsample

  ! Write information on output grid setup
  !***************************************

  if (outgrid_option .eq. 1) then
  write(unitheader) outlon0,outlat0,numxgrid,numygrid, &
       dxoutl,dyoutl
  else
  write(unitheader) outlon0,outlat0,numxgrid,numygrid, &
       dxout,dyout
  endif
  write(unitheader) numzgrid,(outheight(i),i=1,numzgrid)

  call caldate(bdate,jjjjmmdd,ihmmss)
  write(unitheader) jjjjmmdd,ihmmss

  ! Write number of species, and name for each species (+extra name for depositions)
  ! Indicate the dimension of the fields (i.e., 1 for deposition fields, numzgrid for
  ! concentration fields
  !*****************************************************************************

  write(unitheader) 3*nspec,maxpointspec_act
  do i=1,nspec
    write(unitheader) 1,'WD_'//species(i)(1:7)
    write(unitheader) 1,'DD_'//species(i)(1:7)
    write(unitheader) numzgrid,species(i)
  end do

  ! Write information on release points: total number, then for each point:
  ! start, end, coordinates, # of particles, name, mass
  !************************************************************************

  write(unitheader) numpoint
  do i=1,numpoint
    write(unitheader) ireleasestart(i),ireleaseend(i),kindz(i)
    xp1=xpoint1(i)*dx+xlon0
    yp1=ypoint1(i)*dy+ylat0
    xp2=xpoint2(i)*dx+xlon0
    yp2=ypoint2(i)*dy+ylat0
    write(unitheader) xp1,yp1,xp2,yp2,zpoint1(i),zpoint2(i)
    write(unitheader) npart(i),1
    if (numpoint.le.2000) then
      write(unitheader) compoint(i)
    else
      write(unitheader) compoint(2001)
    endif
    do j=1,nspec
      write(unitheader) xmass(i,j)
      write(unitheader) xmass(i,j)
      write(unitheader) xmass(i,j)
    end do
  end do

  ! Write information on some model switches
  !*****************************************

  write(unitheader) method,lsubgrid,lconvection, &
       ind_source,ind_receptor

  ! Write age class information
  !****************************

  write(unitheader) nageclass,(lage(i),i=1,nageclass)


  ! Write topography to output file
  !********************************

  do ix=0,numxgrid-1
    write(unitheader) (oroout(ix,jy),jy=0,numygrid-1)
  end do
  close(unitheader)


      open(53,file=path(1)(1:length(1))//'latlon.txt',form='formatted')
          open(54,file=path(1)(1:length(1))//'latlon_corner.txt' &
          ,form='formatted')

      if (outgrid_option.eq.0) then ! irregular
        call ll_to_xymeter_wrf(outgrid_swlon,outgrid_swlat,xsw,ysw)
        call ll_to_xymeter_wrf(outgrid_nelon,outgrid_nelat,xne,yne)
        do jy=1,numygrid
        do ix=1,numxgrid
          tmpx=out_xm0+(float(ix)-0.5)*dxout
          tmpy=out_ym0+(float(jy)-0.5)*dyout
          call xymeter_to_ll_wrf(tmpx,tmpy,tmplon,tmplat)
        write(53,*) tmplon,tmplat
          tmpx=out_xm0+(float(ix)-1.)*dxout
          tmpy=out_ym0+(float(jy)-1.)*dyout
          call xymeter_to_ll_wrf_out(tmpx,tmpy,tmplon,tmplat)
           write(54,*) tmplon,tmplat
        enddo
        enddo
       else ! regular
        call ll_to_xymeter_wrf(outgrid_swlon,outgrid_swlat,xsw,ysw)
        call ll_to_xymeter_wrf(outgrid_nelon,outgrid_nelat,xne,yne)
        do jy=1,numygrid
        do ix=1,numxgrid
          tmpx=xsw+(xne-xsw)*float(ix-1)/float(numxgrid-1)
          tmpy=ysw+(yne-ysw)*float(jy-1)/float(numygrid-1)
          call xymeter_to_ll_wrf(tmpx,tmpy,tmplon,tmplat)
            xl2=outlon0+(float(ix)-0.5)*dxoutl !long  
            yl2=outlat0+(float(jy)-0.5)*dyoutl !lat  
           write(53,*) xl2,yl2
            xl2=outlon0+float(ix-1)*dxoutl !long 
            yl2=outlat0+float(jy-1)*dyoutl !lat    
           write(54,*) xl2,yl2
        enddo
        enddo
      endif


      close(53)
      close(54)

  return

998   write(*,*) ' #### FLEXPART MODEL ERROR!   THE FILE         #### '
  write(*,*) ' #### '//path(1)(1:length(1))//'header'//' #### '
  write(*,*) ' #### CANNOT BE OPENED. IF A FILE WITH THIS    #### '
  write(*,*) ' #### NAME ALREADY EXISTS, DELETE IT AND START #### '
  write(*,*) ' #### THE PROGRAM AGAIN.                       #### '
  stop

end subroutine writeheader
