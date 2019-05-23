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

    module mpi_mod
  use par_mod, only: dp
  !includes OH concentration field as well as the height information
  !for this field

  implicit none
  integer,allocatable, dimension (:) :: mpi_npoint,mpi_idt,mpi_itra1, &
   mpi_itramem,mpi_nclass
  real,allocatable, dimension (:) :: mpi_uap,mpi_ucp,mpi_uzp,  &
   mpi_us,mpi_vs,mpi_ws,mpi_ztra1
  real,allocatable, dimension (:,:) ::  mpi_xmass1
  real(kind=dp),allocatable, dimension (:) :: mpi_xtra1,mpi_ytra1
  integer(kind=2),allocatable, dimension (:) ::  mpi_cbt
  integer :: chunksize2
  integer, allocatable, dimension (:) :: dummyi2
  real, allocatable, dimension (:) :: dummyr2
  real(kind=dp),allocatable, dimension (:) :: dummyr22
  integer(kind=2),allocatable, dimension (:) :: dummyi22
end module mpi_mod

