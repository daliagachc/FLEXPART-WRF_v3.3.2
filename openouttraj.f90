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

      subroutine openouttraj
!*******************************************************************************
!                                                                              *
!   Note:  This is the FLEXPART_WRF version of subroutine openouttraj.         *
!                                                                              *
!   This routine opens the output file for the plume trajectory output         *
!   produced by the cluster analysis.                                          *
!                                                                              *
!     Author: A. Stohl                                                         *
!     27 January 2001                                                          *
!                                                                              *
!     Dec 2005, R. Easter - changed names of "*lon0*" & "*lat0*" variables     *
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
!                                                                              *
!*******************************************************************************

  use point_mod
  use par_mod
  use com_mod
!      include 'includepar'
!      include 'includecom'
      
      integer i
      real xp1,yp1,xp2,yp2
      real xtmp, ytmp


! Open output file for trajectory output
!***************************************

      open(unitouttraj,file=path(1)(1:length(1))//'trajectories.txt', &
      form='formatted',err=998)

      if (ldirect.eq.1) then
      write(unitouttraj,'(i8,1x,i6,1x,a)') ibdate,ibtime,'FLEXWRF  V3.0'
      else
      write(unitouttraj,'(i8,1x,i6,1x,a)') iedate,ietime,'FLEXWRF  V3.0'
      endif
      write(unitouttraj,*) method,lsubgrid,lconvection
      write(unitouttraj,*) numpoint
      do i=1,numpoint
       if (outgrid_option .eq. 0) then
         xp1=xpoint1(i)*dx+xmet0
         yp1=ypoint1(i)*dy+ymet0
         xp2=xpoint2(i)*dx+xmet0
         yp2=ypoint2(i)*dy+ymet0
        endif
       if (outgrid_option .eq. 1) then
          xtmp = xpoint1(i)*dx+xmet0
          ytmp = ypoint1(i)*dy+ymet0
          call xymeter_to_ll_wrf( xtmp, ytmp, xp1, yp1 )
          xtmp = xpoint2(i)*dx+xmet0
          ytmp = ypoint2(i)*dy+ymet0
          call xymeter_to_ll_wrf( xtmp, ytmp, xp2, yp2 )
        endif
!jdf    write(unitouttraj,*) ireleasestart(i),ireleaseend(i),
!jdf +  xp1,yp1,xp2,yp2,zpoint1(i),zpoint2(i),kindz(i),npart(i)
!jdf    write(unitouttraj,'(a)') compoint(i)(1:40)
        write(unitouttraj,*) ireleasestart(i),ireleaseend(i), &
        xp1,yp1,xp2,yp2,zpoint1(i),zpoint2(i),kindz(i),npart(i)
        write(unitouttraj,'(a20)') compoint(i)(1:20)
       enddo
101     format(2i5,4f11.5,2f11.3,2i5)

      return

998   write(*,*) ' #### FLEXPART MODEL ERROR!   THE FILE         #### '
      write(*,*) ' #### trajectories.txt                         #### '
      write(*,*) ' #### CANNOT BE OPENED. IF A FILE WITH THIS    #### '
      write(*,*) ' #### NAME ALREADY EXISTS, DELETE IT AND START #### '
      write(*,*) ' #### THE PROGRAM AGAIN.                       #### '
      stop

end subroutine openouttraj

