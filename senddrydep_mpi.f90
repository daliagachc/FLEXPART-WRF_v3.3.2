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

  !   routine used to send real array of dry deposition by MPI
  !    Author: J. Brioude                                                      *
  !    January 2013                                                           *
        subroutine senddrydep_mpi(chunksize)

      use mpi_mod
      use com_mod
      use unc_mod, only: drygridunc2 
      use par_mod, only: nclassunc
          implicit none
      include 'mpif.h'
!      character :: varname*20
       integer :: chunksize,numpart2,ks
!      real :: dummyr(numpart2,nspec2),
       real :: dummyr(chunksize)
       integer :: myid,ierr,ntasks,ii,jdeb,jfin,jj,tag,direc
!      integer :: MPI_COMM_WORLD

       integer :: tag2,nspec2,ix,jy,l,nage,kp
       integer :: jj2,from       
  integer, dimension(MPI_STATUS_SIZE) :: status

      call MPI_COMM_RANK ( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE ( MPI_COMM_WORLD, ntasks, ierr )

!    chunksize=numxgrid*numygrid
  do ks=1,nspec
   do kp=1,maxpointspec_act
    do nage=1,nageclass
     do l=1,nclassunc
      tag=ks*1000000+kp*10000+nage*100+l

   if (myid.gt.0) then !slaves send

    call MPI_SEND(drygridunc2(0,0,ks,kp,l,nage), chunksize, MPI_REAL,0,tag, MPI_COMM_WORLD, ierr)

   else ! the master gets
     do from =1,ntasks-1
    call MPI_RECV(dummyr, chunksize, MPI_REAL, from,tag, MPI_COMM_WORLD,status,ierr)
       do jy=0,numygrid-1
        do ix=0,numxgrid-1
         drygridunc2(ix,jy,ks,kp,l,nage)=drygridunc2(ix,jy,ks,kp,l,nage)+ &
           dummyr(ix+1+jy*numxgrid)
        enddo
       enddo
     enddo

   endif

    enddo
    enddo
    enddo
    enddo
       end subroutine senddrydep_mpi
