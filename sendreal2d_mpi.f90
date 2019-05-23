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

  !   routine used to send 2D real array by MPI
  !    Author: J. Brioude                                                      *
  !    March 2012                                                           *
        subroutine sendreal2d_mpi(tag2,numpart2,nspec2,chunksize,direc)

      use mpi_mod
      use com_mod
          implicit none
      include 'mpif.h'

!      character :: varname*20
       integer :: chunksize,numpart2,ks,jj1
!      real :: dummyr(numpart2,nspec2),
!      real :: dummyr2(chunksize)
       integer :: myid,ierr,ntasks,ii,jdeb,jfin,jj,tag,direc
!      integer :: MPI_COMM_WORLD

       integer :: tag2,nspec2 
       integer :: jj2,from,jj3   
  integer, dimension(MPI_STATUS_SIZE) :: status

      call MPI_COMM_RANK ( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE ( MPI_COMM_WORLD, ntasks, ierr )

      do ks=1,nspec2
      tag=100+ks

     if (direc.eq.0) then ! the slaves get

      if (myid.eq.0) then
       do ii=1,ntasks-1
       do jj2=1,chunksize
        jj=(jj2-1)*ntasks+ii+1
!       do jj=ii+1,numpart2+ii,ntasks
!        jj2=(jj-ii-1)/ntasks+1
!        dummyr2(jj2)=dummyr(jj,ks)
         dummyr2(jj2)=xmass1(jj,ks)
        enddo
       call MPI_SEND(dummyr2, chunksize, MPI_REAL, ii,tag,MPI_COMM_WORLD, ierr)
       enddo

!     chunksize2=int((numpart2-1)/ntasks)+1
!     chunksize2=chunksize
      ii=0
      do jj=1,numpart2,ntasks
      ii=ii+1
      jj2=jj
      enddo
      chunksize2=ii+numpart2-jj2
     do jj=1,numpart2,ntasks
      jj3=(jj-1)/ntasks+1
     mpi_xmass1(jj3,ks)=xmass1(jj,ks)
     enddo
     mpi_xmass1(jj3:chunksize2,ks)=xmass1(jj2:numpart2,ks)

   else ! slave receive
    call MPI_RECV(mpi_xmass1(1,ks), chunksize, MPI_REAL, 0,tag, MPI_COMM_WORLD,status, ierr)
   endif

   else ! the master is going to get 

   if (myid.gt.0) then !slaves send

   call MPI_SEND(mpi_xmass1(1,ks), chunksize, MPI_REAL,0,tag, MPI_COMM_WORLD, ierr)

    else ! the master gets
     do from =1,ntasks-1
    call MPI_RECV(dummyr2, chunksize, MPI_REAL, from,tag, MPI_COMM_WORLD,status,ierr)
          jj1=(from-1)*chunksize+1
          jj2=from*chunksize
        xmass1(jj1:jj2,ks)=dummyr2(1:chunksize)
     enddo
     xmass1(jj2+1:numpart2,ks)=mpi_xmass1(1:chunksize2,ks)

    endif

   endif


    enddo
       end subroutine sendreal2d_mpi
