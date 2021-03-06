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

  !   routine used to send integer*2 vectors by MPI
  !    Author: J. Brioude                                                      *
  !    March 2012                                                           *
         subroutine sendint2_mpi(tag,numpart2,chunksize,direc)

      use mpi_mod
      use com_mod

          implicit none
      include 'mpif.h'

!      character :: varname*20
       integer :: chunksize,numpart2,jj1
!      integer(kind=2) :: dummyi(numpart2),
!      integer(kind=2) :: dummyi22(numpart2)
       integer :: myid,ierr,ntasks,ii,jdeb,jfin,jj,direc,tag
!      integer :: MPI_COMM_WORLD

       integer :: jj2,from,jj3   
  integer, dimension(MPI_STATUS_SIZE) :: status
       
      call MPI_COMM_RANK ( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE ( MPI_COMM_WORLD, ntasks, ierr )
     if (direc.eq.0) then ! the slaves get

      if (myid.eq.0) then
       do ii=1,ntasks-1
       do jj2=1,chunksize
        jj=(jj2-1)*ntasks+ii+1
!       do jj=ii+1,numpart2+ii,ntasks
!        jj2=(jj-ii-1)/ntasks+1
!        dummyi2(jj2)=dummyi(jj)
        if (tag.eq.13) dummyi22(jj2)=cbt(jj)
        enddo
    call MPI_SEND(dummyi22, chunksize, MPI_INTEGER2, ii,tag, MPI_COMM_WORLD, ierr)
       enddo
!     chunksize2=int((numpart2-1)/ntasks)+1
!     chunksize2=chunksize
      ii=0
      do jj=1,numpart2,ntasks
      ii=ii+1
      jj2=jj
      enddo
      chunksize2=ii+numpart2-jj2

    if (tag.eq.13) then
     do jj=1,numpart2,ntasks
      jj3=(jj-1)/ntasks+1
     mpi_cbt(jj3)=cbt(jj)
     enddo
     mpi_cbt(jj3:chunksize2)=cbt(jj2:numpart2)

    endif

   else ! the slaves receive
    if (tag.eq.13) call MPI_RECV(mpi_cbt, chunksize, MPI_INTEGER2, 0,13, MPI_COMM_WORLD,status, ierr)
   endif

   else ! the master is going to get 
   if (myid.gt.0) then !slaves send
   if (tag.eq.13) call MPI_SEND(mpi_cbt, chunksize, MPI_INTEGER2, 0,13, MPI_COMM_WORLD, ierr)

    else ! the master gets
     do from =1,ntasks-1
    call MPI_RECV(dummyi22, chunksize, MPI_INTEGER2, from,tag, MPI_COMM_WORLD,status,ierr)
          jj1=(from-1)*chunksize+1
          jj2=from*chunksize
       if (tag.eq.13) cbt(jj1:jj2)=dummyi22(1:chunksize)
     enddo
     if (tag.eq.13) cbt(jj2+1:numpart2)=mpi_cbt(1:chunksize2)

    endif

   endif

       end subroutine sendint2_mpi
