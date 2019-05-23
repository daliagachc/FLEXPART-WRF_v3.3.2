!***********************************************************************
!* Copyright 1998,1999,2000,2001,2002,2005,2007,2008,2009,2010         *
!* Andreas Stohl, Petra Seibert, A. Frank, Gerhard Wotawa,             *
!* Caroline Forster, Sabine Eckhardt, John Burkhart, Harald Sodemann   *
!*                                                                     *
!* This file is part of FLEXPART.                                      *
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
      subroutine convmix(itime)
!                          i
!**************************************************************
!     handles all the calculations related to convective mixing
!     Petra Seibert, Bernd C. Krueger, Feb 2001
!     nested grids included, Bernd C. Krueger, May 2001
!
!     Changes by Caroline Forster, April 2004 - February 2005:
!       convmix called every lsynctime seconds
!     CHANGES by A. Stohl:
!       various run-time optimizations - February 2005
!**************************************************************

  use flux_mod
  use par_mod
  use com_mod
  use conv_mod

  implicit none

  integer :: igr,igrold, ipart, itime, ix, j, inest
  integer :: ipconv,stat
  integer :: jy, kpart, ktop, ngrid,kz
!  integer :: igrid(maxpart), ipoint(maxpart), igridn(maxpart,maxnests)
  integer,allocatable, dimension (:) :: igrid,ipoint
  integer,allocatable, dimension (:,:) :: igridn

  ! itime [s]                 current time
  ! igrid(maxpart)            horizontal grid position of each particle
  ! igridn(maxpart,maxnests)  dto. for nested grids
  ! ipoint(maxpart)           pointer to access particles according to grid position

  logical :: lconv
  real :: x, y, xtn,ytn, ztold, delt
  real :: dt1,dt2,dtt
  integer :: mind1,mind2
  ! dt1,dt2,dtt,mind1,mind2       variables used for time interpolation
  integer :: itage,nage
  real,parameter :: eps=nxmax/3.e5

  real :: duma, dumz(nuvzmax+1)

      write(*,'(//a,a//)') &
          '*** Stopping in subr. convmix ***', &
          '    This is not implemented for FLEXPART_WRF'
      stop

    allocate(igrid(maxpart) ,stat=stat)
    allocate(ipoint(maxpart) ,stat=stat)
    allocate(igridn(maxpart,maxnests) ,stat=stat)


      return
      end subroutine convmix
