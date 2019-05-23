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
      subroutine getfields(itime,nstop)
!                            i     o
!*******************************************************************************
!                                                                              *
!  This subroutine manages the 3 data fields to be kept in memory.             *
!  During the first time step of petterssen it has to be fulfilled that the    *
!  first data field must have |wftime|<itime, i.e. the absolute value of wftime*
!  must be smaller than the absolute value of the current time in [s].         *
!  The other 2 fields are the next in time after the first one.                *
!  Pointers (memind) are used, because otherwise one would have to resort the  *
!  wind fields, which costs a lot of computing time. Here only the pointers are*
!  resorted.                                                                   *
!                                                                              *
!     Author: A. Stohl                                                         *
!                                                                              *
!     29 April 1994                                                            *
!                                                                              *
!  Changes, Bernd C. Krueger, Feb. 2001:                                       *
!        Variables tth,qvh,tthn,qvhn (on eta coordinates) in common block.     *
!        Function of nstop extended.                                           *
!                                                                              *
!  Dec 2005, R. Easter -                                                       *
!          When "memtime(2) = itime = wftime(numbwf)", do not read a new file. *
!          This allows the ending date/time of the flexpart run to match       *
!          the date/time of the last met. file.                                *
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
! lwindinterval [s]    time difference between the two wind fields read in     *
! indj                 indicates the number of the wind field to be read in    *
! indmin               remembers the number of wind fields already treated     *
! memind(2)            pointer, on which place the wind fields are stored      *
! memtime(2) [s]       times of the wind fields, which are kept in memory      *
! itime [s]            current time since start date of trajectory calculation *
! nstop                > 0, if trajectory has to be terminated                 *
! nx,ny,nuvz,nwz       field dimensions in x,y and z direction                 *
! uu(0:nxmax,0:nymax,nuvzmax,2)   wind components in x-direction [m/s]         *
! vv(0:nxmax,0:nymax,nuvzmax,2)   wind components in y-direction [m/s]         *
! ww(0:nxmax,0:nymax,nwzmax,2)    wind components in z-direction [deltaeta/s]  *
! tt(0:nxmax,0:nymax,nuvzmax,2)   temperature [K]                              *
! ps(0:nxmax,0:nymax,2)           surface pressure [Pa]                        *
!                                                                              *
! Constants:                                                                   *
! idiffmax             maximum allowable time difference between 2 wind fields *
!                                                                            *
!*******************************************************************************

  use par_mod
  use com_mod

  implicit none

  integer :: indj,itime,nstop,memaux
  integer :: clck_counts_beg, clck_counts_end, clck_rate
  real :: tins

      real(kind=4) :: uuh(0:nxmax-1,0:nymax-1,nuvzmax)
      real(kind=4) :: vvh(0:nxmax-1,0:nymax-1,nuvzmax)
      real(kind=4) :: wwh(0:nxmax-1,0:nymax-1,nwzmax)
   real(kind=4) :: divh(0:nxmax-1,0:nymax-1,nuvzmax)

!  real :: uuh(0:nxmax-1,0:nymax-1,nuvzmax)
!  real :: vvh(0:nxmax-1,0:nymax-1,nuvzmax)
  real :: pvh(0:nxmax-1,0:nymax-1,nuvzmax)
!  real :: wwh(0:nxmax-1,0:nymax-1,nwzmax)
  real(kind=4) :: uuhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
  real(kind=4) :: vvhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
  real :: pvhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
  real(kind=4) :: wwhn(0:nxmaxn-1,0:nymaxn-1,nwzmax,maxnests)
  real :: start, finish

  integer :: indmin = 1


   real(kind=4) :: divhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
   character(len=28) :: name2
   character(len=8)  :: chartime
    integer :: ix,jy,kz


 
! Check, if wind fields are available for the current time step
!**************************************************************

      nstop=0

      if ((ldirect*wftime(1).gt.ldirect*itime).or. &
      (ldirect*wftime(numbwf).lt.ldirect*itime)) then
        write(*,*) 'FLEXPART WARNING: NO WIND FIELDS ARE AVAILABLE.'
        write(*,*) 'A TRAJECTORY HAS TO BE TERMINATED.'
        nstop=4
        return
      endif


      if ((ldirect*memtime(1).le.ldirect*itime).and. &
      (ldirect*memtime(2).gt.ldirect*itime)) then

! The right wind fields are already in memory -> don't do anything
!*****************************************************************

        continue

! FLEXPART_WRF - following change allows the ending date/time 
! of the flexpart run to match that of the last met. file
      else if ( (ldirect*memtime(1).lt.ldirect*itime).and. &
     (memtime(2).eq.itime) .and. (wftime(numbwf).eq.itime) ) then

        continue

      else if ((ldirect*memtime(2).le.ldirect*itime).and. &
      (memtime(2).ne.999999999)) then
 

! Current time is after 2nd wind field
! -> Resort wind field pointers, so that current time is between 1st and 2nd
!***************************************************************************

        memaux=memind(1)
        memind(1)=memind(2)
        memind(2)=memaux
        memtime(1)=memtime(2)


! Read a new wind field and store it on place memind(2)
!******************************************************

        do indj=indmin,numbwf-1
           if (ldirect*wftime(indj+1).gt.ldirect*itime) then
            call system_clock (clck_counts_beg,clck_rate)
           if ((time_option.eq.1).and.(wind_option.eq.1)) &
              call readwind_timeav(indj+1,memind(2),uuh,vvh,wwh)
            call readwind(indj+1,memind(2),uuh,vvh,wwh,divh)
           if ((time_option.eq.1).and.(wind_option.eq.1)) &
              call readwind_nests_timeav(indj+1,memind(2),uuhn,vvhn,wwhn)
            call readwind_nests(indj+1,memind(2),uuhn,vvhn,wwhn,divhn)
!              if (option_verbose.gt.1) then
            call system_clock (clck_counts_end,clck_rate)
            tins=real(clck_counts_end - clck_counts_beg)/real(clck_rate)
           print*,'readwind',tins
            call system_clock (clck_counts_beg,clck_rate)
!             endif
              call calcpar(memind(2),uuh,vvh,pvh)
              call calcpar_nests(memind(2),uuhn,vvhn,pvhn)
!             if (option_verbose.gt.1) then
            call system_clock (clck_counts_end,clck_rate)
            tins=real(clck_counts_end - clck_counts_beg)/real(clck_rate)
           print*,'calcpar',tins
            call system_clock (clck_counts_beg,clck_rate)
!             endif
              call verttransform(memind(2),uuh,vvh,wwh,pvh,divh)
              call verttransform_nests(memind(2),uuhn,vvhn,wwhn,pvhn,divhn)
!             if (option_verbose.gt.1) then
            call system_clock (clck_counts_end,clck_rate)
            tins=real(clck_counts_end - clck_counts_beg)/real(clck_rate)
           print*,'verttran',tins
!             endif
              memtime(2)=wftime(indj+1)
              nstop = 1
              goto 40
           endif
       enddo
 40     indmin=indj

      else

! No wind fields, which can be used, are currently in memory 
! -> read both wind fields
!***********************************************************

         do indj=indmin,numbwf-1
            if ((ldirect*wftime(indj).le.ldirect*itime).and. &
                 (ldirect*wftime(indj+1).gt.ldirect*itime)) then
               memind(1)=1
         if ((time_option.eq.1).and.(wind_option.eq.1)) &
           call readwind_timeav(indj,memind(1),uuh,vvh,wwh)
         call readwind(indj,memind(1),uuh,vvh,wwh,divh)
         if ((time_option.eq.1).and.(wind_option.eq.1)) &
           call readwind_nests_timeav(indj,memind(1),uuhn,vvhn,wwhn)
         call readwind_nests(indj,memind(1),uuhn,vvhn,wwhn,divhn)
         call calcpar(memind(1),uuh,vvh,pvh)
         call calcpar_nests(memind(1),uuhn,vvhn,pvhn)
         call verttransform(memind(1),uuh,vvh,wwh,pvh,divh)
         call verttransform_nests(memind(1),uuhn,vvhn,wwhn,pvhn,divhn)
         memtime(1)=wftime(indj)
         memind(2)=2
        if ((time_option.eq.1).and.(wind_option.eq.1)) &
           call readwind_timeav(indj+1,memind(2),uuh,vvh,wwh)
        call readwind(indj+1,memind(2),uuh,vvh,wwh,divh)
        if ((time_option.eq.1).and.(wind_option.eq.1)) &
           call readwind_nests_timeav(indj+1,memind(2),uuhn,vvhn,wwhn)
        call readwind_nests(indj+1,memind(2),uuhn,vvhn,wwhn,divhn)
        call calcpar(memind(2),uuh,vvh,pvh)
        call calcpar_nests(memind(2),uuhn,vvhn,pvhn)
        call verttransform(memind(2),uuh,vvh,wwh,pvh,divh)
        call verttransform_nests(memind(2),uuhn,vvhn,wwhn,pvhn,divhn)
        memtime(2)=wftime(indj+1)
        nstop = 1
        goto 60
       endif
      end do
 60      indmin=indj

      endif

      lwindinterv=abs(memtime(2)-memtime(1))

      if (lwindinterv.gt.idiffmax) nstop=3

end subroutine getfields
