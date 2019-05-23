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

  !   routine used to send integer vectors by MPI
  !    Author: J. Brioude                                                      *
  !    March 2012                                                           *
         subroutine sendint_mpi(tag,numpart2,chunksize,direc)

      use mpi_mod
      use com_mod
          implicit none
      include 'mpif.h'

!      character :: varname*20
       integer :: chunksize,numpart2,jj1
!      integer :: dummyi(numpart2),
!      integer :: dummyi2(chunksize),
!      integer :: MPI_COMM_WORLD
       integer :: myid,ierr,ntasks,ii,jdeb,jfin,jj,direc,tag
       integer :: jj2,from,jj3
  integer, dimension(MPI_STATUS_SIZE) :: status
       
      call MPI_COMM_RANK ( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE ( MPI_COMM_WORLD, ntasks, ierr )
     if (direc.eq.0) then ! the slaves get

     if (myid.eq.0) then
      do ii=1,ntasks-1
       do jj2=1,chunksize
        jj=(jj2-1)*ntasks+ii+1

!      do jj=ii+1,numpart2+ii+1,ntasks
!       jj2=(jj-ii-1)/ntasks+1
!       dummyi2(jj2)=dummyi(jj)
       if (tag.eq.1)  dummyi2(jj2)=npoint(jj)
       if (tag.eq.2)  dummyi2(jj2)=idt(jj)
       if (tag.eq.3)  dummyi2(jj2)=itra1(jj)
       if (tag.eq.14)  dummyi2(jj2)=itramem(jj)
       if (tag.eq.15)  dummyi2(jj2)=nclass(jj)
       enddo
       call MPI_SEND(dummyi2, chunksize, MPI_INTEGER, ii,tag, MPI_COMM_WORLD, ierr)
      enddo

!     chunksize2=int((numpart2-1)/ntasks)+1
!     chunksize2=chunksize
      ii=0
      do jj=1,numpart2,ntasks
      ii=ii+1 
      jj2=jj
      enddo
      chunksize2=ii+numpart2-jj2
    if (tag.eq.1) then
     do jj=1,numpart2,ntasks
      jj3=(jj-1)/ntasks+1
     mpi_npoint(jj3)=npoint(jj)
     enddo
     mpi_npoint(jj3:chunksize2)=npoint(jj2:numpart2)
    elseif (tag.eq.2) then
     do jj=1,numpart2,ntasks
      jj3=(jj-1)/ntasks+1
     mpi_idt(jj3)=idt(jj)
     enddo
     mpi_idt(jj3:chunksize2)=idt(jj2:numpart2)
    elseif (tag.eq.3) then
     do jj=1,numpart2,ntasks
      jj3=(jj-1)/ntasks+1
     mpi_itra1(jj3)=itra1(jj)
     enddo
     mpi_itra1(jj3:chunksize2)=itra1(jj2:numpart2)
    elseif (tag.eq.14) then
     do jj=1,numpart2,ntasks
      jj3=(jj-1)/ntasks+1
     mpi_itramem(jj3)=itramem(jj)
     enddo
     mpi_itramem(jj3:chunksize2)=itramem(jj2:numpart2)
    elseif (tag.eq.15) then
     do jj=1,numpart2,ntasks
      jj3=(jj-1)/ntasks+1
     mpi_nclass(jj3)=nclass(jj)
     enddo
     mpi_nclass(jj3:chunksize2)=nclass(jj2:numpart2)
    endif 

    else ! the slaves receive
     if (tag.eq.1) call MPI_RECV(mpi_npoint, chunksize, MPI_INTEGER, 0,1,MPI_COMM_WORLD,status, ierr)
     if (tag.eq.2) call MPI_RECV(mpi_idt, chunksize, MPI_INTEGER, 0,2,MPI_COMM_WORLD,status, ierr)
     if (tag.eq.3) call MPI_RECV(mpi_itra1, chunksize, MPI_INTEGER, 0,3,MPI_COMM_WORLD,status, ierr)
     if (tag.eq.14) call MPI_RECV(mpi_itramem, chunksize, MPI_INTEGER, 0,14,MPI_COMM_WORLD,status, ierr)
     if (tag.eq.15) call MPI_RECV(mpi_nclass, chunksize, MPI_INTEGER, 0,15,MPI_COMM_WORLD,status, ierr)
    endif

   else ! the master is going to get 
   if (myid.gt.0) then !slaves send
   if (tag.eq.1) call MPI_SEND(mpi_npoint, chunksize, MPI_INTEGER, 0,1, MPI_COMM_WORLD, ierr)
   if (tag.eq.2) call MPI_SEND(mpi_idt, chunksize, MPI_INTEGER, 0,2, MPI_COMM_WORLD, ierr)
   if (tag.eq.3) call MPI_SEND(mpi_itra1, chunksize, MPI_INTEGER, 0,3, MPI_COMM_WORLD, ierr)
   if (tag.eq.14) call MPI_SEND(mpi_itramem, chunksize, MPI_INTEGER, 0,14, MPI_COMM_WORLD, ierr)
   if (tag.eq.15) call MPI_SEND(mpi_nclass, chunksize, MPI_INTEGER, 0,15, MPI_COMM_WORLD, ierr)

    else ! the master gets

     do from =1,ntasks-1
    call MPI_RECV(dummyi2, chunksize, MPI_INTEGER, from,tag, MPI_COMM_WORLD,status,ierr)
          jj1=(from-1)*chunksize+1
          jj2=from*chunksize
       if (tag.eq.1) npoint(jj1:jj2)=dummyi2(1:chunksize)
       if (tag.eq.2) idt(jj1:jj2)=dummyi2(1:chunksize)
       if (tag.eq.3) itra1(jj1:jj2)=dummyi2(1:chunksize)
       if (tag.eq.14) itramem(jj1:jj2)=dummyi2(1:chunksize)
       if (tag.eq.15) nclass(jj1:jj2)=dummyi2(1:chunksize)
     enddo
     if (tag.eq.1) npoint(jj2+1:numpart2)=mpi_npoint(1:chunksize2)
     if (tag.eq.2) idt(jj2+1:numpart2)=mpi_idt(1:chunksize2)
     if (tag.eq.3) itra1(jj2+1:numpart2)=mpi_itra1(1:chunksize2)
     if (tag.eq.14) itramem(jj2+1:numpart2)=mpi_itramem(1:chunksize2)
     if (tag.eq.15) nclass(jj2+1:numpart2)=mpi_nclass(1:chunksize2)

    endif

   endif

       end subroutine sendint_mpi
