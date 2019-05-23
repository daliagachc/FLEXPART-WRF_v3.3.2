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
      program flexwrf_mpi
!
!*******************************************************************************
!                                                                              *
!     This is the Lagrangian Particle Dispersion Model FLEXPART_WRF.           *
!                                                                              *
!         FLEXPART uses met. files from the ECMWF model (in grib format),      *
!             and its internal computational grid is latitude-longitude.       *
!                                                                              *
!         FLEXPART_WRF uses met. files from the WRF model (in NetCDF format),  *
!             and its internal computational grid is the WRF x-y grid.         *
!                                                                              *
!     The main program manages the reading of model run specifications, etc.   *
!     All actual computing is done within subroutine timemanager.              *
!                                                                              *
!     Author: A. Stohl                                                         *
!     18 May 1996                                                              *
!                                                                              *
!     Nov 2005, R. Easter - Added the above comments and changed               *
!                           the program name to "flexpart_wrf"                 *
!                                                                              *
!     Feb 2012, J Brioude- modify the name of the pilt_wrf model from PNLL     *
!                          to flexwrf.
!                          start doing versions                                *
!                input information should be put in flexwr.input               *
!     Mar 2012, J Brioude: Hybrid parallelization of v74. everything converted *
!                   in fortran 90, based on version 90.1 of the main stream    *
!                   version of FLEXPART.                                       *
!     Jun 2012, J Brioude: Add tests on arguments to the flexwrf input file.   *
!                                                                              *
!     2015-03-25, A. Dingwell: prepared dry-deposition section for supporting  *
!                   land-use data from WRF                                     *
!     2015-05-27  A. Dingwell: Now reads roughness length when landuse is read *
!                   from WRF.                                                  *
!     2015-06-05  A. Dingwell: Added module netcdf_output_mod, updated calls   *
!                   new/renamed subroutines.                                   *
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
!                                                                              *
! Constants:                                                                   *
!                                                                              *
!*******************************************************************************


  use point_mod
  use par_mod
  use com_mod
  use conv_mod
  use netcdf_output_mod   ! Contains subroutines for netcdf output format

  !-- modules below added by mc for parallel random generaiton using teh RANLUX
  !and Mersenne-Twister generator
  use luxury
  use mt_stream
  implicit none

  include 'mpif.h'

  integer :: i,j,ix,jy,inest,ii
!  integer :: MPI_COMM_WORLD
  integer :: idummy 
  integer :: inext,inextp,ma(55),iff
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer :: myid,ntasks,ierr,islave
! ,len2
! character(len=200) :: inputname
  real, external :: ran3  ! added by mc to use RAN3 using the original JB random number system 
  integer, allocatable :: seed(:) ! here and below further variable used by the MT generator 
  integer(4) :: iseed = 73519232
  integer :: id
  type (mt_state) :: mts (0: MAX_STREAM)
  character :: nummpi_id*2  !for test on pc

!  integer, parameter :: master=0, mstgtag1=11, msgtag2=12


!      if (myid.eq.0) then
!let's comment the line above to let each node reading and making the same
!thing.

!  save inext,inextp,ma,iff
  iff=0
!    print*,'before 1'
  call MPI_INIT( ierr )
  call MPI_COMM_RANK ( MPI_COMM_WORLD, myid, ierr )
  call MPI_COMM_SIZE ( MPI_COMM_WORLD, ntasks, ierr )

!    print*,'after 1'

  print*,'MPI id',myid
! Generate a large number of random numbers
!         write(*,9100) 'degrees longitude,latitude', 
!     +                 'degrees longitude,latitude'
!      else
!         write(*,9100) 'grid-meters', 'grid-meters'
!      endif
!9100  format( '    x,y coordinates in input  files must be in ', a /
!     +        '    x,y coordinates in output files will be in ', a / )


     if (command_argument_count().eq.0) then 
     print*,'the input file used is flexwrf.input in the ' // &
      'local folder of the executable'
     inputname='flexwrf.input' 
     endif
     if (command_argument_count().gt.0) then  
         call get_command_argument(1,inputname,len2,ierr)
     print*,'the input file used is ' // inputname
     endif
! Generate a large number of random numbers
!******************************************
      if (newrandomgen.eq.0) then
       idummy = -320-(myid*4049)
!      idummy = -320
      do j=1,maxomp
      do i=1,maxrand-1,2
      ii=i+(j-1)*maxrand        
        call gasdev1(idummy,rannumb(ii),rannumb(ii+1),inext,inextp,ma,iff)
      enddo
      enddo
      ii=maxrand*maxomp
      call gasdev1(idummy,rannumb(ii),rannumb(ii-1),inext,inextp,ma,iff)
!     print*,'rand',myid
!     print*,rannumb(1:5)
!     call ranlux(uniform_rannumb,maxrandomp)  ! this generate a uniform distribution
      else
      idummy=254+myid*443   !different seed for different mpi processes are produced so indepedent stream for any mpi process suing RANLUX are certain
      call RLUXGO(3,idummy,0,0)  ! this set the luxury level to 3 and initalize the generator for any myid 
      do i=1,maxrand-1,2
      call gasdevlux2R(rannumb(i),rannumb(i+1)) !this will generate a guassian distribution 
      end do
      call gasdevlux2R(rannumb(maxrand),rannumb(maxrand-1))
      ! Generate a stream of uniform deviate random numbers to be used for CBL
      call ranlux(uniform_rannumb,maxrand)  ! this generate a uniform distribution

!----- comment by MC: now initialize the mersenne twister generator for a number
!max_stream of possible streams 
!----- to be called subsequently by any openmp process activated. note RANLUX
!above is suppose to be the best generator
!----- but it is slower than mersenne twister and moreover it would require some
!adaptation for workiong with openmp processes
      !do this on any mpi_process taht will have  a copy of all the MT generator
      !initialization
      ! set parameters

      call set_mt19937

      !  initialize MT state type
      call new (mts(0))

      call init (mts(0),iseed)   !iseed unique and defined above. note that the lenght of the period of the master stream is about 2^19000

      !  initialize additional streams from the master. this is done jumping
      !  between different points in the stream any child stream has period
      !  2^256
      do id=1, MAX_STREAM
      call create_stream (mts(0),mts(id),id)
      end do
      end if

! Read the unified input file - jdf
!***************************

      call readinput

! Get terrain and dry deposition resistances
!*******************************************
  if ( DRYDEP ) then
    if (lu_option.eq.0) then
      if (option_verbose.ge.1) then
        write(*,*) 'Main: Old IGBP data will be used for land use categories'
      endif
      call readlanduse  ! Read the landuse inventory & assign z0
      call assignland   ! Assign fractional cover of landuse to 1st WRF grid
    else
      if(option_verbose.ge.1) then
        write(*,*) 'Main: Land use will be updated at every input interval'
        call read_z0    ! Only read relation between land-use and z0
      endif
      !write(*,*) ' #### FLEXPART MODEL ERROR! LAND USE FROM WRF #### '
      !write(*,*) ' #### IS NOT YET SUPPORTED. YOU SHOULD USE    #### '
      !write(*,*) ' #### LU_OPTION=0 FOR NOW...                  #### '
      !stop
    endif
    call readdepo     ! Read and compute surface resistances to dry deposition of gases
  endif
  
! Convert the release point coordinates from geografical to grid coordinates
!***************************************************************************

      call coordtrafo

! Initialize all particles to non-existent
!*****************************************

!      do j=1,maxpart
!        itra1(j)=-999999999
!      enddo

! For continuation of previous run, read in particle positions
!*************************************************************

      if (ipin.eq.1) then
        call readpartpositions
      else
        numpart=0
    numparticlecount=0

      endif


! Calculate volume, surface area, etc., of all output grid cells
!***************************************************************
!     if (myid.eq.0) then
      if (outgrid_option.eq.0) then 
      call outgrid_init_irreg
      if (nested_output.eq.1) call outgrid_init_nest_irreg() !need to be fixed
      elseif (outgrid_option.eq.1) then
      call outgrid_init_reg
      if (nested_output.eq.1) call outgrid_init_nest_reg() !need to be fixed
      endif
!     endif
!     if (nested_output.eq.1) call outgrid_init_nest() !need to be fixed


  ! Read the OH field
  !******************

  if (OHREA.eqv..TRUE.) &
       call readohfield

! Write basic information on the simulation to a file "header"
! and open files that are to be kept open throughout the simulation
!******************************************************************

      if (myid.eq.0) then
      if (iouttype.eq.2 .and. ipout.gt.0)  call writeheader

      if (iouttype.eq.0 .or. iouttype.eq.1) then ! binary or ascii output
        call writeheader
        if (nested_output.eq.1) call writeheader_nest() !need to be fixed
      else  ! netcdf output
        call nc_create_header_outfile(0,0)
        if (nested_output.eq.1)  call nc_create_header_outfile(0,1)
      endif !iouttype
!  open(unitdates,file=path(2)(1:length(2))//'dates')
      open(unitdates,file=path(1)(1:length(1))//'dates')
      call openreceptors
      if ((iout.eq.4).or.(iout.eq.5)) call openouttraj

      endif
! Releases can only start and end at discrete times (multiples of lsynctime)
!***************************************************************************
      
      do i=1,numpoint
        ireleasestart(i)=nint(real(ireleasestart(i))/ &
        real(lsynctime))*lsynctime
        ireleaseend(i)=nint(real(ireleaseend(i))/ &
        real(lsynctime))*lsynctime
        enddo

! Initialize cloud-base mass fluxes for the convection scheme
!************************************************************

      do jy=0,nymin1
        do ix=0,nxmin1
        cbaseflux(ix,jy)=0.
    end do
  end do

      do  inest=1,numbnests
        do  jy=0,nyn(inest)-1
          do  ix=0,nxn(inest)-1
          cbasefluxn(ix,jy,inest)=0.
    end do
  end do
  end do


! Calculate particle trajectories
!********************************
!     endif !if condition on myid
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)


!     print*,'entering timemanager'
      call timemanager_mpi(mts)


      write(*,'(/a/)') 'CONGRATULATIONS: YOU HAVE SUCCESSFULLY ' //  &
        'COMPLETED A FLEXPART_WRF MODEL RUN!'

   call MPI_FINALIZE ( ierr )

end program flexwrf_mpi 

