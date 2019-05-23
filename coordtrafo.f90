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
      subroutine coordtrafo
!**********************************************************************
!                                                                     * 
! Note:  This is the FLEXPART_WRF version of subroutine coordtrafo.   *
!                                                                     * 
!             FLEXPART MODEL SUBROUTINE COORDTRAFO                    *
!                                                                     *
!**********************************************************************
!                                                                     * 
! AUTHOR:      G. WOTAWA                                              *
! DATE:        1994-02-07                                             *
! LAST UPDATE: 1996-05-18   A. STOHL                                  *
!                                                                     * 
! Dec 2005, R. Easter - changed names of "*lon0*" & "*lat0*" variables*
! July 2012, J Brioude: modification following flexpart 9             * 
!**********************************************************************
!                                                                     *
! DESCRIPTION: This subroutine transforms x and y coordinates of      *
! particle release points to grid coordinates.                        *
!                                                                     *
!**********************************************************************

  use point_mod
  use par_mod
  use com_mod

  implicit none

  integer :: i,j,k

      if (numpoint.eq.0) goto 30

! TRANSFORM X- AND Y- COORDINATES OF STARTING POINTS TO GRID COORDINATES
!***********************************************************************

      do i=1,numpoint
        xpoint1(i)=(xpoint1(i)-xmet0)/dx
        xpoint2(i)=(xpoint2(i)-xmet0)/dx
        ypoint1(i)=(ypoint1(i)-ymet0)/dy
      ypoint2(i)=(ypoint2(i)-ymet0)/dy
   end do

15    continue


! CHECK IF RELEASE POINTS ARE WITHIN DOMAIN
!******************************************

      do i=1,numpoint

      if ((ypoint1(i).lt.1.e-6).or.(ypoint1(i).ge.real(nymin1)-1.e-6)  &
      .or.(ypoint2(i).lt.1.e-6).or.(ypoint2(i).ge.real(nymin1)-1.e-6)  &
      .or.(xpoint1(i).lt.1.e-6).or.(xpoint1(i).ge.real(nxmin1)-1.e-6) &
      .or.(xpoint2(i).lt.1.e-6).or.(xpoint2(i).ge.real(nxmin1)-1.e-6)) then

!      if ((ypoint1(i).lt.1.e-6).or.(ypoint1(i).ge.real(nymin1)-1.e-6) &
!      .or.(ypoint2(i).lt.1.e-6).or.(ypoint2(i).ge.real(nymin1)-1.e-6) &
!      .or.((.not.xglobal).and.((xpoint1(i).lt.1.e-6).or. &
!      (xpoint1(i).ge.real(nxmin1)-1.e-6).or.(xpoint2(i).lt.1.e-6).or. &
!      (xpoint2(i).ge.real(nxmin1)-1.e-6)))) then

          write(*,*) ' NOTICE: RELEASE POINT OUT OF DOMAIN DETECTED.'
          write(*,*) ' IT IS REMOVED NOW ... '
          write(*,*) ' DEBUG INFO:'
          write(*,*) ' Point number:',i
          write(*,*) ' xpoint1,ypoint1:',xpoint1,ypoint1
          write(*,*) ' xpoint2,ypoint2:',xpoint2,ypoint2
!         write(*,*) ' COMMENT: ',compoint(i)

          if (i.lt.numpoint) then
            do j=i+1,numpoint
              xpoint1(j-1)=xpoint1(j)
              ypoint1(j-1)=ypoint1(j)
              xpoint2(j-1)=xpoint2(j)
              ypoint2(j-1)=ypoint2(j)
              zpoint1(j-1)=zpoint1(j)
              zpoint2(j-1)=zpoint2(j)
              npart(j-1)=npart(j)
           kindz(j-1)=kindz(j)
          ireleasestart(j-1)=ireleasestart(j)
           ireleaseend(j-1)=ireleaseend(j)
           if (j.le.2000) compoint(j-1)=compoint(j)
           do k=1,nspec
             xmass(j-1,k)=xmass(j,k)
          end do

         enddo
          endif

          numpoint=numpoint-1
          if (numpoint.gt.0) goto 15

        endif
  end do

30    if (numpoint.eq.0) then
        write(*,*) ' FLEXPART MODEL SUBROUTINE COORDTRAFO: ERROR ! '
        write(*,*) ' NO PARTICLE RELEASES ARE DEFINED!'
        write(*,*) ' CHECK FILE RELEASES...'
        stop
      endif

end subroutine coordtrafo
 
