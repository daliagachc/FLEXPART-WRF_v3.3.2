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
!*******************************************************************************
!   Include file for convection
!   This file contains a global common block used by convect
!   and other subroutines
!   Author: P. Ferstl
!
!   Feb 2001
!
!*******************************************************************************

module conv_mod

  use par_mod, only: nconvlevmax, na, nxmax, nymax, nxmaxn, nymaxn, maxnests, nuvzmax

  implicit none

  !integer,parameter :: nconvlevmax = nuvzmax-1, &
  !                     na = nconvlevmax+1
  !these parameters are defined in par_mod now!

  real :: pconv(nconvlevmax),phconv(na),dpr(nconvlevmax)
  real :: pconv_hpa(nconvlevmax),phconv_hpa(na)

  real :: ft(nconvlevmax), fq(nconvlevmax)
  real :: fmass(nconvlevmax,nconvlevmax),sub(nconvlevmax)
  real :: fmassfrac(nconvlevmax,nconvlevmax)
  real :: cbaseflux(0:nxmax-1,0:nymax-1)
  real :: cbasefluxn(0:nxmaxn-1,0:nymaxn-1,maxnests)
  real :: tconv(na),qconv(na),qsconv(na)
  real :: psconv,tt2conv,td2conv
      real :: umf3(0:nxmax-1,0:nymax-1,nuvzmax)
      real :: uer3(0:nxmax-1,0:nymax-1,nuvzmax)
      real :: udr3(0:nxmax-1,0:nymax-1,nuvzmax)
      real :: dmf3(0:nxmax-1,0:nymax-1,nuvzmax)
      real :: der3(0:nxmax-1,0:nymax-1,nuvzmax)
      real :: ddr3(0:nxmax-1,0:nymax-1,nuvzmax)
      real :: cu_top(0:nxmax-1,0:nymax-1)
      real :: cu_bot(0:nxmax-1,0:nymax-1)

      real :: umf3n(0:nxmax-1,0:nymax-1,nuvzmax,maxnests)
      real :: uer3n(0:nxmax-1,0:nymax-1,nuvzmax,maxnests)
      real :: udr3n(0:nxmax-1,0:nymax-1,nuvzmax,maxnests)
      real :: dmf3n(0:nxmax-1,0:nymax-1,nuvzmax,maxnests)
      real :: der3n(0:nxmax-1,0:nymax-1,nuvzmax,maxnests)
      real :: ddr3n(0:nxmax-1,0:nymax-1,nuvzmax,maxnests)
      real :: cu_topn(0:nxmax-1,0:nymax-1,maxnests)
      real :: cu_botn(0:nxmax-1,0:nymax-1,maxnests)

  integer :: nconvlev,nconvtop

end module conv_mod
