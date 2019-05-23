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

  !   routine used to send real vectors by MPI
  !    Author: J. Brioude                                                      *
  !    March 2012                                                           *
         subroutine sendreal_mpi(tag,numpart2,chunksize,direc)

      use mpi_mod
      use com_mod
          implicit none
      include 'mpif.h'

!      character :: varname*20
       integer :: chunksize,numpart2,jj1
!      real :: dummyr(numpart),
!      real :: dummyr2(chunksize)
       integer :: myid,ierr,ntasks,ii,jdeb,jfin,jj,direc,tag
!      integer :: MPI_COMM_WORLD
       integer :: jj2,from,jj3
  integer, dimension(MPI_STATUS_SIZE) :: status

      call MPI_COMM_RANK ( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE ( MPI_COMM_WORLD, ntasks, ierr )

     if (direc.eq.0) then ! the slaves get

     if (myid.eq.0) then ! master sends
      do ii=1,ntasks-1
!      jdeb=(ii-1)*chunksize+1
!      jfin=(ii)*chunksize
!      do jj=jdeb,jfin
!      do jj=ii+1,numpart2+ii,ntasks
       do jj2=1,chunksize
!       jj2=jj-jdeb+1
!       jj2=(jj-ii-1)/ntasks+1
        jj=(jj2-1)*ntasks+ii+1
!       dummyr2(jj2)=dummyr(jj)
       if (tag.eq.4)  dummyr2(jj2)=uap(jj)
       if (tag.eq.5)  dummyr2(jj2)=ucp(jj)
       if (tag.eq.6)  dummyr2(jj2)=uzp(jj)
       if (tag.eq.7)  dummyr2(jj2)=us(jj)
       if (tag.eq.8)  dummyr2(jj2)=vs(jj)
       if (tag.eq.9)  dummyr2(jj2)=ws(jj)
       if (tag.eq.10)  dummyr2(jj2)=ztra1(jj)
       enddo
       call MPI_SEND(dummyr2, chunksize, MPI_REAL, ii,tag, MPI_COMM_WORLD, ierr)
      enddo

!     jdeb=(ntasks-1)*chunksize+1
!     jfin=numpart
!     chunksize2=jfin-jdeb+1
!     chunksize2=int((numpart2-1)/ntasks)+1
!      chunksize2=chunksize
      ii=0
      do jj=1,numpart2,ntasks
      ii=ii+1
      jj2=jj
      enddo
      chunksize2=ii+numpart2-jj2

    if (tag.eq.4) then
     do jj=1,numpart2,ntasks
      jj3=(jj-1)/ntasks+1
     mpi_uap(jj3)=uap(jj)
     enddo
     mpi_uap(jj3:chunksize2)=uap(jj2:numpart2)
    elseif (tag.eq.5) then
     do jj=1,numpart2,ntasks
      jj3=(jj-1)/ntasks+1
     mpi_ucp(jj3)=ucp(jj)
     enddo
     mpi_ucp(jj3:chunksize2)=ucp(jj2:numpart2)
    elseif (tag.eq.6) then
     do jj=1,numpart2,ntasks
      jj3=(jj-1)/ntasks+1
     mpi_uzp(jj3)=uzp(jj)
     enddo
     mpi_uzp(jj3:chunksize2)=uzp(jj2:numpart2)
    elseif (tag.eq.7) then
     do jj=1,numpart2,ntasks
      jj3=(jj-1)/ntasks+1
     mpi_us(jj3)=us(jj)
     enddo
     mpi_us(jj3:chunksize2)=us(jj2:numpart2)
    elseif (tag.eq.8) then
     do jj=1,numpart2,ntasks
      jj3=(jj-1)/ntasks+1
     mpi_vs(jj3)=vs(jj)
     enddo
     mpi_vs(jj3:chunksize2)=vs(jj2:numpart2)
    elseif (tag.eq.9) then
     do jj=1,numpart2,ntasks
      jj3=(jj-1)/ntasks+1
     mpi_ws(jj3)=ws(jj)
     enddo
     mpi_ws(jj3:chunksize2)=ws(jj2:numpart2)
    elseif (tag.eq.10) then
     do jj=1,numpart2,ntasks
      jj3=(jj-1)/ntasks+1
     mpi_ztra1(jj3)=ztra1(jj)
     enddo
     mpi_ztra1(jj3:chunksize2)=ztra1(jj2:numpart2)
    endif
    
     else ! the slaves receive
    if (tag.eq.4) call MPI_RECV(mpi_uap, chunksize, MPI_REAL, 0,4,MPI_COMM_WORLD,status, ierr)
    if (tag.eq.5) call MPI_RECV(mpi_ucp, chunksize, MPI_REAL, 0,5,MPI_COMM_WORLD,status, ierr)
    if (tag.eq.6) call MPI_RECV(mpi_uzp, chunksize, MPI_REAL, 0,6,MPI_COMM_WORLD,status, ierr)
    if (tag.eq.7) call MPI_RECV(mpi_us, chunksize, MPI_REAL, 0,7,MPI_COMM_WORLD,status, ierr)
    if (tag.eq.8) call MPI_RECV(mpi_vs, chunksize, MPI_REAL, 0,8,MPI_COMM_WORLD,status, ierr)
    if (tag.eq.9) call MPI_RECV(mpi_ws, chunksize, MPI_REAL, 0,9,MPI_COMM_WORLD,status, ierr)
    if (tag.eq.10) call MPI_RECV(mpi_ztra1, chunksize, MPI_REAL, 0,10,MPI_COMM_WORLD,status, ierr)
   endif

   else ! the master is going to get 

   if (myid.gt.0) then !slaves send

   if (tag.eq.4) call MPI_SEND(mpi_uap, chunksize, MPI_REAL, 0,4, MPI_COMM_WORLD, ierr)
   if (tag.eq.5) call MPI_SEND(mpi_ucp, chunksize, MPI_REAL, 0,5, MPI_COMM_WORLD, ierr)
   if (tag.eq.6) call MPI_SEND(mpi_uzp, chunksize, MPI_REAL, 0,6, MPI_COMM_WORLD, ierr)
   if (tag.eq.7) call MPI_SEND(mpi_us, chunksize, MPI_REAL, 0,7, MPI_COMM_WORLD, ierr)
   if (tag.eq.8) call MPI_SEND(mpi_vs, chunksize, MPI_REAL, 0,8, MPI_COMM_WORLD, ierr)
   if (tag.eq.9) call MPI_SEND(mpi_ws, chunksize, MPI_REAL, 0,9, MPI_COMM_WORLD, ierr)
   if (tag.eq.10) call MPI_SEND(mpi_ztra1, chunksize, MPI_REAL, 0,10, MPI_COMM_WORLD, ierr)

    else ! the master gets

     do from =1,ntasks-1
    call MPI_RECV(dummyr2, chunksize, MPI_REAL, from,tag, MPI_COMM_WORLD,status,ierr)
          jj1=(from-1)*chunksize+1
          jj2=from*chunksize
       if (tag.eq.4) uap(jj1:jj2)=dummyr2(1:chunksize)
       if (tag.eq.5) ucp(jj1:jj2)=dummyr2(1:chunksize)
       if (tag.eq.6) uzp(jj1:jj2)=dummyr2(1:chunksize)
       if (tag.eq.7) us(jj1:jj2)=dummyr2(1:chunksize)
       if (tag.eq.8) vs(jj1:jj2)=dummyr2(1:chunksize)
       if (tag.eq.9) ws(jj1:jj2)=dummyr2(1:chunksize)
       if (tag.eq.10) ztra1(jj1:jj2)=dummyr2(1:chunksize)

     enddo
     if (tag.eq.4) uap(jj2+1:numpart2)=mpi_uap(1:chunksize2)
     if (tag.eq.5) ucp(jj2+1:numpart2)=mpi_ucp(1:chunksize2)
     if (tag.eq.6) uzp(jj2+1:numpart2)=mpi_uzp(1:chunksize2)
     if (tag.eq.7) us(jj2+1:numpart2)=mpi_us(1:chunksize2)
     if (tag.eq.8) vs(jj2+1:numpart2)=mpi_vs(1:chunksize2)
     if (tag.eq.9) ws(jj2+1:numpart2)=mpi_ws(1:chunksize2)
     if (tag.eq.10) ztra1(jj2+1:numpart2)=mpi_ztra1(1:chunksize2)


    endif

   endif
       end subroutine sendreal_mpi
