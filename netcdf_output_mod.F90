!***********************************************************************
! Copyright 2015                                                       *
!   Wayne Angevine        Delia Arnold          Jerome Brioude         *
!   John Burkhart         Massimo Cassiani      Adam Dingwell          *
!   Richard C Easter      Sabine Eckhardt       Stephanie Evan         *
!   Jerome D Fast         Caroline Forster      Don Morton             *
!   Ignacio Pisso         Petra Seibert         Harald Sodemann        *
!   Andreas Stohl         Gerard Wotawa                                *
!                                                                      * 
!                                                                      *
!  This file is part of FLEXPART WRF                                   *
!                                                                      *
! FLEXPART is free software: you can redistribute it and/or modify     *
! it under the terms of the GNU General Public License as published by *
! the Free Software Foundation, either version 3 of the License, or    *
! (at your option) any later version.                                  *
!                                                                      *
! FLEXPART is distributed in the hope that it will be useful,          *
! but WITHOUT ANY WARRANTY; without even the implied warranty of       *
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        *
! GNU General Public License for more details.                         *
!                                                                      *
! You should have received a copy of the GNU General Public License    *
! along with FLEXPART.  If not, see <http://www.gnu.org/licenses/>.    *
!***********************************************************************

!*******************************************************************************
!                                                                              *
!   This module contains the subroutines used to generate NetCDF output files. *
!                                                                              *
!       Author: A. Dingwell                                                    *
!                                                                              *
!       05 June 2015                                                           *
!                                                                              *
!   Note that the file-extension of this file should be ".F" (upper case)      *
!   This ensure that it is sent to a preprocessor when compiling with any of   *
!   the supported compilers (PGI,Intel or GNU)                                 *
!   The preprocessor is required to determine if netcdf4 compressed output     *
!   be included in the program or not.                                         *
!                                                                              *
!*******************************************************************************

module netcdf_output_mod

  ! uses wrf_map_utils_mod
  use netcdf              ! netcdf module
  use par_mod             ! Global parameters
  use com_mod             ! Global variables
  use point_mod           !
  use outg_mod            ! Output grid definition
  use wrf_map_utils_mod   ! Map projection variables

  
  implicit none
#ifdef NETCDF4_OUTPUT
  integer, parameter :: deflate_level=4   ! compression (1-9, 9 is most intense)
  logical, parameter :: shuffle=.FALSE.   ! For compressing integer arrays
#endif
  
  ! Common names used for defining output variables:
  character, parameter :: vname_t*5   = 'Times'
  character, parameter :: vname_x*5   = 'XLONG'
  character, parameter :: vname_y*4   = 'XLAT'
  
  ! Common names used for naming variable attributes:
  character, parameter :: descr*11    = 'description'
  character, parameter :: units*5     = 'units'
  character, parameter :: coord*11    = 'coordinates'
  character, parameter :: coordxy*10  = 'XLONG XLAT'
  
  
  contains
  subroutine check_ncerror(errcode)

    !*****************************************************************************
    !                                                                            *
    ! This function checks the return value of any call to the netcdf interface. *
    ! The subroutine should be called directly after any call to any nf_*        *
    ! functions, unless some other means of erro handling has been implemented.  *
    !                                                                            *
    !     Author: A. Dingwell                                                    *
    !                                                                            *
    !     27 May 2013                                                            *
    !                                                                            *
    !   2015-06-04: A. Dingwell: Updated subroutine to use NetCDF F90 module.    *
    !*****************************************************************************
    
    integer, intent(in) :: errcode

    if( errcode.ne.NF90_NOERR ) then
      print*, 'Error: ', nf90_strerror(errcode)
      stop
    endif
    return
  end subroutine check_ncerror
  
  subroutine nc_write_global_attributes(ncid,nesting_level)
  !*****************************************************************************
  !                                                                            *
  !  This subroutine writes global attributes to the given output file.        *
  !  It should be called once for each netcdf file (header or data) in the     *
  !  when the respective file is open in define mode                           *
  !                                                                            *
  !      Author: A. Dingwell                                                   *
  !                                                                            *
  !      27 May 2013                                                           *
  !                                                                            *
  !*****************************************************************************
  
    integer, intent(in) :: ncid           ! File handle for output file
    integer, intent(in) :: nesting_level  ! Which grid we should describe
    
    integer ::  ncret ! Return value of calls fot NF*


    integer :: ncgrid_dx,ncgrid_dy        ! dx,dy of current grid in m or latlon
    integer :: ncgrid_nx,ncgrid_ny        ! nx,ny of current grid
    
    ! Suggestion: put these values in an array under the parent module instead 
    ! and simply address the respective element when this function is called:
    if (nesting_level.eq.0) then
      ncgrid_nx = numxgrid
      ncgrid_ny = numygrid
      if (outgrid_option.eq.1) then ! input was in latlon
        ncgrid_dx = dxoutl
        ncgrid_dy = dyoutl
      else  ! input was in metres
        ncgrid_dx = dxout
        ncgrid_dy = dyout
      endif
    elseif (nesting_level.eq.1) then  ! current grid is nested
      ncgrid_nx = numxgridn
      ncgrid_ny = numygridn
      if (outgrid_option.eq.1) then ! input was in latlon
        ncgrid_dx = dxoutln
        ncgrid_dy = dyoutln
      else  ! input was in metres
        ncgrid_dx = dxoutn
        ncgrid_dy = dyoutn
      endif
    endif
  
    if (ldirect.eq.1) then  ! Forward simulation
      if (option_verbose.ge.10) write(*,10) 'forward simulation attributes'
      ncret = NF90_PUT_ATT(ncid,NF90_GLOBAL,'SIMULATION_START_DATE',ibdate)
      call check_ncerror(ncret)
      ncret = NF90_PUT_ATT(ncid,NF90_GLOBAL,'SIMULATION_START_TIME',ibtime)
      call check_ncerror(ncret)

      ncret = NF90_PUT_ATT(ncid,NF90_GLOBAL,'SIMULATION_END_DATE',iedate)
      call check_ncerror(ncret)
      ncret = NF90_PUT_ATT(ncid,NF90_GLOBAL,'SIMULATION_END_TIME',ietime)
      call check_ncerror(ncret)
    else                  ! Backward simulation
      if (option_verbose.ge.10) write(*,10) 'backward simulation attributes'
      ncret = NF90_PUT_ATT(ncid,NF90_GLOBAL,'SIMULATION_START_DATE',iedate)
      call check_ncerror(ncret)
      ncret = NF90_PUT_ATT(ncid,NF90_GLOBAL,'SIMULATION_START_TIME',ietime)
      call check_ncerror(ncret)

      ncret = NF90_PUT_ATT(ncid,NF90_GLOBAL,'SIMULATION_END_DATE',ibdate)
      call check_ncerror(ncret)
      ncret = NF90_PUT_ATT(ncid,NF90_GLOBAL,'SIMULATION_END_TIME',ibtime)
      call check_ncerror(ncret)
    endif

    if (option_verbose.ge.10) write(*,10) 'map projection attributes'
    ncret = NF90_PUT_ATT(ncid,NF90_GLOBAL,'MAP_PROJ',map_proj_id)
    call check_ncerror(ncret)
    if (outgrid_option .eq. 1) then ! native lon-lat grid.
      ncret = &
        NF90_PUT_ATT(ncid,NF90_GLOBAL,'OUTPUT_PROJECTION','Regular Latit/Longit')
      call check_ncerror(ncret)
    else ! Using map projection, add projection info to output:
      ncret = NF90_PUT_ATT(ncid,NF90_GLOBAL,'EARTH_RADIUS_M',earth_radius_m)
      call check_ncerror(ncret)
      ncret = NF90_PUT_ATT(ncid,NF90_GLOBAL,'CEN_LAT',proj_clat)
      call check_ncerror(ncret)
      ncret = NF90_PUT_ATT(ncid,NF90_GLOBAL,'CEN_LON',proj_clon)
      call check_ncerror(ncret)
      
      if (map_proj_id.eq.1) then  ! Lamber conformal projection
        !ncret = &
        !  NF90_PUT_ATT(ncid,NF90_GLOBAL,'OUTPUT_PROJECTION','Lambert conformal')
        !call check_ncerror(ncret)
        ! Attributes on WRF format:
        ncret = NF90_PUT_ATT(ncid,NF90_GLOBAL,'MAP_PROJ_CHAR','Lambert Conformal')
        call check_ncerror(ncret)
        ncret = NF90_PUT_ATT(ncid,NF90_GLOBAL,'STAND_LON',proj_stdlon)
        call check_ncerror(ncret)
        ncret = NF90_PUT_ATT(ncid,NF90_GLOBAL,'TRUELAT1',proj_truelat1)
        call check_ncerror(ncret)
        ncret = NF90_PUT_ATT(ncid,NF90_GLOBAL,'TRUELAT2',proj_truelat2)
        call check_ncerror(ncret)
      elseif (map_proj_id.eq.2) then
        ncret = &
          NF90_PUT_ATT(ncid,NF90_GLOBAL,'OUTPUT_PROJECTION','stereographic')
        call check_ncerror(ncret)
      elseif (map_proj_id.eq.3) then
        ncret = NF90_PUT_ATT(ncid,NF90_GLOBAL,'MAP_PROJ_CHAR','Mercator')
        call check_ncerror(ncret)
        ncret = NF90_PUT_ATT(ncid,NF90_GLOBAL,'TRUELAT1',proj_truelat1)
        call check_ncerror(ncret)
      elseif (map_proj_id.eq.4) then
        ncret = &
          NF90_PUT_ATT(ncid,NF90_GLOBAL,'OUTPUT_PROJECTION','global')
        call check_ncerror(ncret)
      endif
      
      if (nesting_level.eq.0) then  ! main domain
        ncret = NF90_PUT_ATT(ncid,NF90_GLOBAL,'OUTLON0',outlon0)
        call check_ncerror(ncret)
        ncret = NF90_PUT_ATT(ncid,NF90_GLOBAL,'OUTLAT0',outlat0)
        call check_ncerror(ncret)
      elseif (nesting_level.eq.1) then  ! nested domain
        ! Should this be done separately like this?
        ncret = NF90_PUT_ATT(ncid,NF90_GLOBAL,'OUTLON0',outlon0n)
        call check_ncerror(ncret)
        ncret = NF90_PUT_ATT(ncid,NF90_GLOBAL,'OUTLAT0',outlat0n)
        call check_ncerror(ncret)
      endif
    endif

    ! Write info common model settings
    !*********************************
    if (option_verbose.ge.10) write(*,10) 'common model attributes'

    if (option_verbose.ge.10) write(*,10) 'OUTPUT_INTERVAL'
    ncret = NF90_PUT_ATT(ncid,NF90_GLOBAL,'OUTPUT_INTERVAL',loutstep)
    call check_ncerror(ncret)

    if (option_verbose.ge.10) write(*,10) 'AVERAGING_TIME'
    ncret = NF90_PUT_ATT(ncid,NF90_GLOBAL,'AVERAGING_TIME',loutaver)
    call check_ncerror(ncret)

    if (option_verbose.ge.10) write(*,10) 'AVERAGE_SAMPLING'
    ncret = NF90_PUT_ATT(ncid,NF90_GLOBAL,'AVERAGE_SAMPLING',loutsample)
    call check_ncerror(ncret)

    ncret = NF90_PUT_ATT(ncid,NF90_GLOBAL,'NSPEC',nspec)
    call check_ncerror(ncret)
    ncret = NF90_PUT_ATT(ncid,NF90_GLOBAL,'NUMRECEPTOR',numreceptor)
    call check_ncerror(ncret)
    ncret = NF90_PUT_ATT(ncid,NF90_GLOBAL,'NAGECLASS',nageclass)
    call check_ncerror(ncret)

    ncret = NF90_PUT_ATT(ncid,NF90_GLOBAL,'NUMRELEASES',numpoint)
    call check_ncerror(ncret)

    ncret = NF90_PUT_ATT(ncid,NF90_GLOBAL,'DISPERSION_METHOD',method)
    call check_ncerror(ncret)

    ncret = NF90_PUT_ATT(ncid,NF90_GLOBAL,'SUBGRID_TOPOGRAPHY',lsubgrid)
    call check_ncerror(ncret)

    ncret = NF90_PUT_ATT(ncid,NF90_GLOBAL,'CONVECTION_PARAM',lconvection)
    call check_ncerror(ncret)

    ncret = NF90_PUT_ATT(ncid,NF90_GLOBAL,'SUBGRID_TOPOGRAPHY',lsubgrid)
    call check_ncerror(ncret)
    
    ncret = NF90_PUT_ATT(ncid,NF90_GLOBAL,'LU_OPTION',lu_option)
    call check_ncerror(ncret)

    ! Write information on output grid setup
    !***************************************
    if (option_verbose.ge.10) write(*,10) 'WEST-EAST_GRID_DIMENSION'
    ncret = NF90_PUT_ATT(ncid,NF90_GLOBAL,'WEST-EAST_GRID_DIMENSION', &
      ncgrid_nx)
    call check_ncerror(ncret)

    if (option_verbose.ge.10) write(*,10) 'SOUTH-NORTH_GRID_DIMENSION'
    ncret = NF90_PUT_ATT(ncid,NF90_GLOBAL,'SOUTH-NORTH_GRID_DIMENSION', &
      ncgrid_ny)
    call check_ncerror(ncret)
    
    if (option_verbose.ge.10) write(*,10) 'BOTTOM-TOP_GRID_DIMENSION'
    ncret = NF90_PUT_ATT(ncid,NF90_GLOBAL,'BOTTOM-TOP_GRID_DIMENSION', &
      numzgrid)

    if (option_verbose.ge.10) write(*,10) 'DX and DY'
    ncret = NF90_PUT_ATT(ncid,NF90_GLOBAL,'DX',ncgrid_dx)
    call check_ncerror(ncret)

    ncret = NF90_PUT_ATT(ncid,NF90_GLOBAL,'DY',ncgrid_dy)
    call check_ncerror(ncret)
    
    10 format('nc_write_global_attributes: Setting up ',A)
  end subroutine nc_write_global_attributes
  
  subroutine nc_create_main_outfile(itime,nesting_level)
  !*****************************************************************************
  !                                                                            *
  !  This routine perdefines a netcdf ouput file with information on flexpart  *
  !  settings, releases and topography.                                        *
  !                                                                            *
  !      Author: A. Dingwell                                                   *
  !                                                                            *
  !      27 May 2013                                                           *
  !                                                                            *
  ! Modifications:                                                             *
  ! June 5 2013: J. Brioude: Create and write attributes to netcdf output only *
  !  2015-04-28: A. Dingwell: Cleaned up comments and indentation              *
  !  2015-05-04: A. Dingwell: Added receptor points to netcdf output           *
  !  2015-06-04: A. Dingwell: Updated subroutine to use NetCDF F90 module.     *
  !  2015-06-05: A. Dingwell: Moved global attributes to separate subroutine   *
  !*****************************************************************************
  
    integer, intent(in) :: itime          ! seconds since simulation start
    integer, intent(in) :: nesting_level  ! 0=main grid (mother) 1=nest (child)
                            ! this is written to be easy to expand if additional 
                            ! nests are desired in the future

    real(kind=dp) :: jul          ! Julian date
    integer   :: jjjjmmdd,ihmmss  ! date & time as integer
    character :: adate*8,atime*6  ! date and time strings, used for filename

    ! Grid related variables
    real    :: xp1,yp1,xp2,yp2  ! temporary coordinates
    real    :: xsw,xne,ysw,yne,tmpx,tmpy,tmplon,tmplat,xl2,yl2
    integer :: ncgrid_nx,ncgrid_ny        ! nx,ny of current grid
    real    :: ncgrid_swlon,ncgrid_swlat  ! SW corner of current grid in latlon
    real    :: ncgrid_nelon,ncgrid_nelat  ! NE corner of current grid in latlon
    real    :: ncgrid_xm0,ncgrid_ym0      ! lower-left grid coord in metres
    real    :: ncgrid_lon0,ncgrid_lat0    ! lower-left grid coord in latlon

    ! Iterators
    integer i,j,ix,jy

    ! NETCDF file related variables
    integer nclvlid,nclonid,nclatid,ncrecid,ncspcid,ncageid !outgrid dimension ids
    integer ncrelid,ncrseid                                 ! release points dimension ids
    integer ncrepid                             ! receptor points dimension ids
    integer ncrnvid,ncrmvid,ncspvid             ! release points: number,mass,species ids
    integer ncrtvid,ncrxvid,ncryvid,ncrzvid     ! release points: t,x,y,z min/max limits
    integer nctovid,ncarvid                     ! Topography and grid area variable-ids
    integer ncstr1id,ncstr2id,ncstr3id          ! description string length dimid
    integer nclvlvid,nclonvid,nclatvid,ncspcvid,ncagevid  ! outgrid dimension variables
    integer ncdimsid3(6),ncdimsid2(5) ! arrays of dimension ids for outgrid 3D & 2D
    integer ncdimsid32(7),ncdimsid22(6) ! arrays of dimension ids for outgrid 3D & 2D
    !integer ncdimsid_times(2) ! Used to defin time dimension variable

    ! NETCDF filename & attribute related variables
    !character descr*11,units*5,ncname*29,coord*11,coordxy*10
    character ncname*29
    !integer coordxylen
    character unit2d*10   ! unit for deposition fields
    !integer   unit2dlen   ! length of character string

    ! NETCDF misc variables
    integer ncid    ! local container for netcdf file-id (either ncout or ncoutn)
    integer ncret   ! Return-value of calls to nf_* utils
    
    ! integer :: chunks(2) ! shuffle
    ! Attribute notation:

    !coordxylen = 10

    ! Determine current calendar date, needed for the file name
    !**********************************************************
    jul=bdate+real(itime,kind=dp)/86400._dp
    call caldate(jul,jjjjmmdd,ihmmss)
    write(adate,'(i8.8)') jjjjmmdd
    write(atime,'(i6.6)') ihmmss

    !************************
    ! Open header output file
    !************************
    write(ncname,'(A8,I2.2,A1,I8.8,A1,I6.6,A3)') &
      'flxout_d',nesting_level+1,'_',jjjjmmdd,'_',ihmmss,'.nc' ! filename

  ! print*,'step0',itime,jjjjmmdd,ihmmss

  !  call nf_set_log_level(3)
    if (option_verbose.ge.1) write(*,*) &
      'write_ncinfo: creating file: ',path(1)(1:length(1))//ncname
  ! call nf_set_chunk_cache(32000000)
  ! ncret = nf_create(path(1)(1:length(1))//ncname, nf_clobber,ncid)
    ncret = NF90_CREATE(path=path(1)(1:length(1))//ncname, cmode=NF90_NETCDF4,ncid=ncid)
    call check_ncerror(ncret)
    

    ! Determine which nest/outfile we just created so we can set up the grid
    !***********************************************************************
    ! Suggestion: put these values in an array under the parent module instead 
    ! and simply address the respective element when this function is called:
    if (nesting_level.eq.0) then  ! current grid is main grid
      ncout   = ncid  ! copy current file handle to ncout
      ncgrid_nx = numxgrid
      ncgrid_ny = numygrid
      ncgrid_nelon = outgrid_nelon
      ncgrid_nelat = outgrid_nelat
      ncgrid_swlon = outgrid_swlon
      ncgrid_swlat = outgrid_swlat
  !   allocate(ncgrid_oro(ncgrid_nx,ncgrid_ny),stat=stat)
  !   allocate(ncgrid_area(ncgrid_nx,ncgrid_ny),stat=stat)
  !   ncgrid_oro   = oroout(0:ncgrid_nx-1,0:ncgrid_ny-1)
  !   ncgrid_area  = area(0:ncgrid_nx-1,0:ncgrid_ny-1)
  !   print*,'step2'
      if (outgrid_option.eq.1) then ! input was in latlon
        ncgrid_lon0 = outlon0
        ncgrid_lat0 = outlat0
      else  ! input was in metres
        ncgrid_xm0  = out_xm0
        ncgrid_ym0  = out_ym0
      endif
    elseif (nesting_level.eq.1) then  ! current grid is nested
      ncoutn  = ncid  ! copy current file handle to ncoutn
      ncgrid_nx = numxgridn
      ncgrid_ny = numygridn
      ncgrid_nelon = outgridn_nelon
      ncgrid_nelat = outgridn_nelat
      ncgrid_swlon = outgridn_swlon
      ncgrid_swlat = outgridn_swlat
  !   allocate(ncgrid_oro(ncgrid_nx,ncgrid_ny),stat=stat)
  !   allocate(ncgrid_area(ncgrid_nx,ncgrid_ny),stat=stat)
  !   ncgrid_oro   = orooutn(0:ncgrid_nx-1,0:ncgrid_ny-1)
  !   ncgrid_area  = arean(0:ncgrid_nx-1,0:ncgrid_ny-1)
      if (outgrid_option.eq.1) then ! input was in latlon
        ncgrid_lon0 = outlon0n
        ncgrid_lat0 = outlat0n
      else  ! input was in metres
        ncgrid_xm0  = out_xm0n
        ncgrid_ym0  = out_ym0n
      endif
    endif
  !   print*,'step3'

    if (option_verbose.ge.10) &
      write(*,*) 'write_ncinfo: ncout,ncoutn=',ncout,ncoutn

    ! Write global attributes
    !*****************************
    call nc_write_global_attributes(ncid,nesting_level)

    ! Set up netcdf dimensions
    !*************************
    if (option_verbose.ge.10) write(*,10) 'main grid dimensions'
    
    ncret = NF90_DEF_DIM(ncid,'Time',NF90_UNLIMITED,ncrecid)
    call check_ncerror(ncret)

    ncret = NF90_DEF_DIM(ncid,'DateStrLen',15,ncstr3id) !TODO: WRF format
    call check_ncerror(ncret)

    ncret = NF90_DEF_DIM(ncid,'west_east',ncgrid_nx,nclonid)
    call check_ncerror(ncret)

    ncret = NF90_DEF_DIM(ncid,'south_north',ncgrid_ny,nclatid)
    call check_ncerror(ncret)

    ncret = NF90_DEF_DIM(ncid,'bottom_top',numzgrid,nclvlid)
    call check_ncerror(ncret)

    ncret = NF90_DEF_DIM(ncid,'species',nspec,ncspcid)
    call check_ncerror(ncret)

    ncret = NF90_DEF_DIM(ncid,'SpeciesStrLen',10,ncstr1id)
    call check_ncerror(ncret)

    ncret = NF90_DEF_DIM(ncid,'ageclass',nageclass,ncageid)
    call check_ncerror(ncret)

    if (option_verbose.ge.10) write(*,10) 'release point dimensions'
    ncret = NF90_DEF_DIM(ncid,'releases',numpoint,ncrelid)
    call check_ncerror(ncret)
    
    ncret = NF90_DEF_DIM(ncid,'ReleaseStrLen',45,ncstr2id)
    call check_ncerror(ncret)
    
    ncret = NF90_DEF_DIM(ncid,'ReleaseStartEnd',2,ncrseid)
    call check_ncerror(ncret)
    
    if(option_verbose.ge.10) write(*,10) 'receptor point dimensions'
    ncret = NF90_DEF_DIM(ncid,'receptors',numreceptor,ncrepid)
    call check_ncerror(ncret)

    ! Select which dimensions to use for main output grids
    if ((ldirect.eq.1).and.(maxpointspec_act.gt.1)) then
      ncdimsid32(1) = nclonid ! X
      ncdimsid32(2) = nclatid ! Y
      ncdimsid32(3) = nclvlid ! Z
      ncdimsid32(4) = ncrelid ! points
      ncdimsid32(5) = ncspcid ! species
      ncdimsid32(6) = ncageid ! ageclass
      ncdimsid32(7) = ncrecid ! t

      ncdimsid22(1) = nclonid ! X
      ncdimsid22(2) = nclatid ! Y
      ncdimsid22(3) = ncrelid ! points
      ncdimsid22(4) = ncspcid ! species
      ncdimsid22(5) = ncageid ! ageclass
      ncdimsid22(6) = ncrecid ! t
    else
      ncdimsid3(1) = nclonid ! X
      ncdimsid3(2) = nclatid ! Y
      ncdimsid3(3) = nclvlid ! Z
      if (ldirect.eq.1) ncdimsid3(4) = ncspcid ! species
      if (ldirect.eq.-1) ncdimsid3(4) = ncrelid ! points
      ncdimsid3(5) = ncageid ! ageclass
      ncdimsid3(6) = ncrecid ! t

      ncdimsid2(1) = nclonid ! X
      ncdimsid2(2) = nclatid ! Y
      if (ldirect.eq.1) ncdimsid2(3) = ncspcid ! species
      if (ldirect.eq.-1) ncdimsid2(3) = ncrelid ! points
      ncdimsid2(4) = ncageid ! ageclass
      ncdimsid2(5) = ncrecid ! t
    endif

    ! TIMES
    if (option_verbose.ge.10) write(*,10) 'TIMES dimension variable'
#ifdef NETCDF4_OUTPUT
      !ncret = NF90_DEF_VAR(ncid,'Times',NF90_CHAR,dimids=(/ncstr3id,ncrecid/),varid=ncrecvid, &
      ncret = NF90_DEF_VAR(ncid,vname_t,NF90_CHAR,(/ncstr3id,ncrecid/),ncrecvid, &
          shuffle = shuffle, deflate_level=deflate_level)
#else
      ncret = NF90_DEF_VAR(ncid,vname_t,NF90_CHAR,(/ncstr3id,ncrecid/),ncrecvid)
#endif
    call check_ncerror(ncret)
    ncret = NF90_PUT_ATT(ncid,ncrecvid,descr, &
        'TIME OF OUTPUT (END OF AVERAGING INTERVAL)')

    ! RECEPTOR POINTS
    ! only write receptors to domain 1 (maybe change this)
    if (nesting_level.eq.0 .and. numreceptor.gt.0) then  
      if (ldirect.eq.1) then  ! Forward run
        if (option_verbose.ge.10) write(*,10) 'RECEPTORS variable'
        if ((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)) then ! CONCENTRATION OUTPUT?
#ifdef NETCDF4_OUTPUT
          ncret = NF90_DEF_VAR(ncid,'RECEPTOR_CONCENTRATION', &
              NF90_FLOAT,(/ncrepid,ncspcid,ncrecid/),ncrcovid, &
              shuffle = shuffle, deflate_level=deflate_level)
#else
          ncret = NF90_DEF_VAR(ncid,'RECEPTOR_CONCENTRATION', &
            NF90_FLOAT,(/ncrepid,ncspcid,ncrecid/),ncrcovid)
#endif
          call check_ncerror(ncret)
          ncret = NF90_PUT_ATT(ncid,ncrcovid,descr, &
              'CONCENTRATION AT RECEPTOR POINTS')
          call check_ncerror(ncret)
          ncret = NF90_PUT_ATT(ncid,ncrcovid,units,'ng m-3')
          call check_ncerror(ncret)
        endif
        
        if ((iout.eq.2).or.(iout.eq.3)) then  ! MIXING RATIO OUTPUT? 
#ifdef NETCDF4_OUTPUT
          ncret = NF90_DEF_VAR(ncid,'RECEPTOR_MIXINGRATIO', &
              NF90_FLOAT,(/ncrepid,ncspcid,ncrecid/),ncrmivid,&
              shuffle = shuffle, deflate_level=deflate_level)
#else
          ncret = NF90_DEF_VAR(ncid,'RECEPTOR_MIXINGRATIO', &
              NF90_FLOAT,(/ncrepid,ncspcid,ncrecid/),ncrmivid)
          call check_ncerror(ncret)
#endif
          ncret = NF90_PUT_ATT(ncid,ncrmivid,descr, &
              'MASS MIXING RATIO AT RECEPTOR POINTS')
          call check_ncerror(ncret)
          ncret = NF90_PUT_ATT(ncid,ncrmivid,units,'ppt by mass')
          call check_ncerror(ncret)
        endif

      endif ! ldirect.eq.1
    end if ! nesting_level.eq.0

    ! MAIN OUTPUT VARIABLES
    if ((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)) then ! OUTPUT CONCENTRATION
      if (option_verbose.ge.10) write(*,10) 'CONC variable'
  !   print*,ncdimsid3
      if ((ldirect.eq.1).and.(maxpointspec_act.gt.1)) then
#ifdef NETCDF4_OUTPUT
        ncret = NF90_DEF_VAR(ncid,'CONC',NF90_FLOAT,ncdimsid32,nccovid,&
            shuffle = shuffle, deflate_level=deflate_level)
#else
        ncret = NF90_DEF_VAR(ncid,'CONC',NF90_FLOAT,ncdimsid32,nccovid)
#endif
      else
#ifdef NETCDF4_OUTPUT
        ncret = NF90_DEF_VAR(ncid,'CONC',NF90_FLOAT,ncdimsid3,nccovid,&
            shuffle = shuffle, deflate_level=deflate_level)
#else
        ncret = NF90_DEF_VAR(ncid,'CONC',NF90_FLOAT,ncdimsid3,nccovid)
#endif
      endif
      call check_ncerror(ncret)
      ncret = NF90_PUT_ATT(ncid,nccovid,descr, &
        'CONCENTRATION OF AIRBORNE SPECIES')
      call check_ncerror(ncret)
      ncret = NF90_PUT_ATT(ncid,nccovid,coord,coordxy)  ! Coordinate string
      call check_ncerror(ncret)
    endif

    if ((iout.eq.2).or.(iout.eq.3)) then  ! MIXING RATIO
      if (option_verbose.ge.10) write(*,10) 'MIXINGRATIO variable'
    if ((ldirect.eq.1).and.(maxpointspec_act.gt.1)) then
#ifdef NETCDF4_OUTPUT
      ncret = NF90_DEF_VAR(ncid,'MIXINGRATIO',NF90_FLOAT,ncdimsid32,ncravid,&
            shuffle = shuffle, deflate_level=deflate_level)
#else
      ncret = NF90_DEF_VAR(ncid,'MIXINGRATIO',NF90_FLOAT,ncdimsid32,ncravid)
#endif
    else
#ifdef NETCDF4_OUTPUT
      ncret = NF90_DEF_VAR(ncid,'MIXINGRATIO',NF90_FLOAT,ncdimsid3,ncravid,&
            shuffle = shuffle, deflate_level=deflate_level)
#else
      ncret = NF90_DEF_VAR(ncid,'MIXINGRATIO',NF90_FLOAT,ncdimsid3,ncravid)
#endif
    endif
    ncret = NF90_PUT_ATT(ncid,ncravid,descr, &
        'MASS MIXING RATIO OF AIRBORNE SPECIES')
    call check_ncerror(ncret)
    ncret = NF90_PUT_ATT(ncid,ncravid,coord,coordxy)  ! Coordinate string
    call check_ncerror(ncret)
    endif

    if (ldirect.eq.1) then  ! Forward run
      unit2d = 'pg m-2'
      !unit2dlen = 6
      if (option_verbose.ge.10) write(*,10) 'DRYDEP variable'
  !   write(*,*) ncdimsid2
      if ((ldirect.eq.1).and.(maxpointspec_act.gt.1)) then
#ifdef NETCDF4_OUTPUT
        ncret = NF90_DEF_VAR(ncid,'DRYDEP',NF90_FLOAT,ncdimsid2,ncddvid,&
            shuffle = shuffle, deflate_level=deflate_level)
#else
        ncret = NF90_DEF_VAR(ncid,'DRYDEP',NF90_FLOAT,ncdimsid22,ncddvid)
#endif
      else
#ifdef NETCDF4_OUTPUT
        ncret = NF90_DEF_VAR(ncid,'DRYDEP',NF90_FLOAT,ncdimsid2,ncddvid,&
            shuffle = shuffle, deflate_level=deflate_level)
#else
        ncret = NF90_DEF_VAR(ncid,'DRYDEP',NF90_FLOAT,ncdimsid2,ncddvid)
#endif
      endif
      call check_ncerror(ncret)
      ncret = NF90_PUT_ATT(ncid,ncddvid,descr, &
          'ACCUMULATED TOTAL DRY DEPOSITION')
      call check_ncerror(ncret)
      ncret = NF90_PUT_ATT(ncid,ncddvid,units,unit2d)   ! unit string
      call check_ncerror(ncret)
      ncret = NF90_PUT_ATT(ncid,ncddvid,coord,coordxy)  ! coordinate string
      call check_ncerror(ncret)

      if (option_verbose.ge.10) write(*,10) 'WETDEP variable'
      if ((ldirect.eq.1).and.(maxpointspec_act.gt.1)) then
#ifdef NETCDF4_OUTPUT
        ncret = NF90_DEF_VAR(ncid,'WETDEP',NF90_FLOAT,ncdimsid22,ncwdvid,&
            shuffle = shuffle, deflate_level=deflate_level)
#else
        ncret = NF90_DEF_VAR(ncid,'WETDEP',NF90_FLOAT,ncdimsid22,ncwdvid)
#endif
      else
#ifdef NETCDF4_OUTPUT
        ncret = NF90_DEF_VAR(ncid,'WETDEP',NF90_FLOAT,ncdimsid2,ncwdvid,&
            shuffle = shuffle, deflate_level=deflate_level)
#else
        ncret = NF90_DEF_VAR(ncid,'WETDEP',NF90_FLOAT,ncdimsid2,ncwdvid)
#endif
      endif
      call check_ncerror(ncret)
      ncret = NF90_PUT_ATT(ncid,ncwdvid,descr, &
          'ACCUMULATED TOTAL WET DEPOSITION')
      call check_ncerror(ncret)
      ncret = NF90_PUT_ATT(ncid,ncwdvid,units,unit2d)   ! unit string
      call check_ncerror(ncret)
      ncret = NF90_PUT_ATT(ncid,ncwdvid,coord,coordxy)  ! coordinate string
      call check_ncerror(ncret)
      
      ! Add unit attr to mixing ratio and concentration fields
      if (ind_samp.eq.0) then   ! AD: Not sure if this should be done like this in forward runs...
        ncret = NF90_PUT_ATT(ncid,nccovid,units,'ng m-3') !CONC
        call check_ncerror(ncret)
        ncret = NF90_PUT_ATT(ncid,ncravid,units,'ppt by mass')  !MIX
        call check_ncerror(ncret)
      else
        ncret = NF90_PUT_ATT(ncid,nccovid,units,'ppt by mass') !CONC
        call check_ncerror(ncret)
        !    ncret = nf_put_att_text(ncid,ncravid,units,3,'???')  !MIX
        !    call check_ncerror(ncret)
      endif
    else                    ! Backward run
      if ((ind_rel.eq.0).and.(ind_samp.eq.0)) then
        ncret = NF90_PUT_ATT(ncid,nccovid,units,'s') !CONC
      else if ((ind_rel.eq.0).and.(ind_samp.eq.-1)) then
        ncret = NF90_PUT_ATT(ncid,nccovid,units,'s m3 kg-1') !CONC
      else if ((ind_rel.eq.1).and.(ind_samp.eq.0)) then
        ncret = NF90_PUT_ATT(ncid,nccovid,units,'s kg m-3') !CONC
      else if ((ind_rel.eq.1).and.(ind_samp.eq.-1)) then
        ncret = NF90_PUT_ATT(ncid,nccovid,units,'s') !CONC
      end if
      call check_ncerror(ncret)
    endif ! Backward/Forward run

    ! EXIT DEFINE MODE, ENTER DATA MODE
    ncret = NF90_ENDDEF(ncid)
    call check_ncerror(ncret)
    
    return

  10 format('nc_create_main_outfile: Setting up ',A)
  end subroutine nc_create_main_outfile

  subroutine nc_create_header_outfile(itime,nesting_level)
  !*****************************************************************************
  !                                                                            *
  !  This routine create the netcdf header file for FLEXPART-WRF               *
  !                                                                            *
  !      Author: A. Dingwell                                                   *
  !                                                                            *
  !      27 May 2013                                                           *
  !                                                                            *
  !  Modifications                                                             *
  !  June 5 2013: J. Brioude: generate a header*nc                             * 
  !  2015-05-04: A. Dingwell: Added receptor points to netcdf output           *
  !  2015-05-07: H.M.J. Barbosa: Receptor coordinates should be in latlon      *
  !                 AD: Use method fix regardless of outgrid_option            *
  !  2015-06-04: A. Dingwell: Updated subroutine to use NetCDF F90 module.     *
  !  2015-06-05: A. Dingwell: Moved global attributes to separate subroutine   *
  !*****************************************************************************
  integer, intent(in) :: itime          ! seconds since simulation start
  integer, intent(in) :: nesting_level  ! 0=main grid (mother) 1=nest (child)
                            ! this is written to be easy to expand if additional 
                            ! are desired in the future
  integer   :: stat
  real(kind=dp) :: jul          ! Julian date
  integer   :: jjjjmmdd,ihmmss  ! date & time as integer
  character :: adate*8,atime*6  ! date and time strings, used for filename

  ! Grid related variables
  real    :: xp1,yp1,xp2,yp2  ! temporary coordinates
  real    :: xsw,xne,ysw,yne,tmpx,tmpy,tmplon,tmplat,xl2,yl2
  integer :: ncgrid_nx,ncgrid_ny        ! nx,ny of current grid
  integer :: ncgrid_dx,ncgrid_dy        ! dx,dy of current grid in m or latlon
  real    :: ncgrid_swlon,ncgrid_swlat  ! SW corner of current grid in latlon
  real    :: ncgrid_nelon,ncgrid_nelat  ! NE corner of current grid in latlon
  real    :: ncgrid_xm0,ncgrid_ym0      ! lower-left grid coord in metres
  real    :: ncgrid_lon0,ncgrid_lat0    ! lower-left grid coord in latlon

  ! Grid related 2D-variables (reassigning these here is a bit inefficient but
  ! it lets us keep a consistent structure of the code, besides it's only once
  ! per output
  real,allocatable,dimension (:,:)  :: ncgrid_oro,ncgrid_area ! of current grid

  ! Iterators
  integer i,j,ix,jy

  ! NETCDF file related variables
  integer nclvlid,nclonid,nclatid,ncrecid,ncspcid,ncageid !outgrid dimension ids
  integer ncrelid,ncrseid                                 ! release points dimension ids
  integer ncrnvid,ncrmvid,ncspvid             ! release points: number,mass,species ids
  integer ncrtvid,ncrxvid,ncryvid,ncrzvid     ! release points: t,x,y,z min/max limits
  integer ncrexvid,ncreyvid,ncrenvid,ncrepid,ncresid   ! Receptor variable ids (move to header?)
  integer nctovid,ncarvid                     ! Topography and grid area variable-ids
  integer ncstr1id,ncstr2id,ncstr3id          ! description string length dimid
  integer nclvlvid,nclonvid,nclatvid,ncspcvid,ncagevid  ! outgrid dimension variables
  integer nclonvid2,nclatvid2
  integer ncdimsid3(6),ncdimsid2(5) ! arrays of dimension ids for outgrid 3D & 2D

  ! NETCDF filename & attribute related variables
  !character descr*11,units*5,coord*11,coordxy*10
  character ncname*29

  ! NETCDF misc variables
  integer ncid    ! local container for netcdf file-id (either ncout or ncoutn)
  integer ncret   ! Return-value of calls to nf_* utils
  
! integer :: chunks(2) ! shuffle

  !hmjb: avoid creation of temporary vectors while calling netcdf
  integer :: ones2(2) = (/1, 1/)
  integer :: xyvec(2)

  !Receptor coordinates
  real :: x_rep_tmp(maxreceptor),y_rep_tmp(maxreceptor)
  real :: x_rep_tmpb,y_rep_tmpb

  ! Attribute notation:
  !descr = 'description'
  !units = 'units'
  !coord = 'coordinates'
  !coordxy = 'XLONG XLAT'

  ! Determine current calendar date, needed for the file name
  !**********************************************************
  jul=bdate+real(itime,kind=dp)/86400._dp
  call caldate(jul,jjjjmmdd,ihmmss)
  write(adate,'(i8.8)') jjjjmmdd
  write(atime,'(i6.6)') ihmmss

  !************************
  ! Open header output file
  !************************
! write(ncname,'(A8,I2.2,A1,I8.8,A1,I6.6,A3)') &
!   'flxout_d',nesting_level+1,'_',jjjjmmdd,'_',ihmmss,'.nc' ! filename
  write(ncname,'(A8,I2.2,A3)') &
    'header_d',nesting_level+1,'.nc' ! filename


!  call nf_set_log_level(3)
  if (option_verbose.ge.1) write(*,*) &
    'nc_create_header_outfile: creating file: ',path(1)(1:length(1))//ncname
! call nf_set_chunk_cache(32000000)
! ncret = nf_create(path(1)(1:length(1))//ncname, nf_clobber,ncid)
  ncret = NF90_create(path(1)(1:length(1))//ncname, cmode=NF90_NETCDF4,ncid=ncid)
  call check_ncerror(ncret)

  ! Determine which nest/outfile we just created so we can set up the grid
  !***********************************************************************
    ! Suggestion: put these values in an array under the parent module instead 
    ! and simply address the respective element when this function is called:
  if (nesting_level.eq.0) then  ! current grid is main grid
    ncout   = ncid  ! copy current file handle to ncout
    ncgrid_nx = numxgrid
    ncgrid_ny = numygrid
    ncgrid_nelon = outgrid_nelon
    ncgrid_nelat = outgrid_nelat
    ncgrid_swlon = outgrid_swlon
    ncgrid_swlat = outgrid_swlat
    allocate(ncgrid_oro(ncgrid_nx,ncgrid_ny),stat=stat)
    allocate(ncgrid_area(ncgrid_nx,ncgrid_ny),stat=stat)
    ncgrid_oro   = oroout(0:ncgrid_nx-1,0:ncgrid_ny-1)
    ncgrid_area  = area(0:ncgrid_nx-1,0:ncgrid_ny-1)
    if (outgrid_option.eq.1) then ! input was in latlon
      ncgrid_dx = dxoutl
      ncgrid_dy = dyoutl
      ncgrid_lon0 = outlon0
      ncgrid_lat0 = outlat0
    else  ! input was in metres
      ncgrid_dx = dxout
      ncgrid_dy = dyout
      ncgrid_xm0  = out_xm0
      ncgrid_ym0  = out_ym0
    endif
  elseif (nesting_level.eq.1) then  ! current grid is nested
    ncoutn  = ncid  ! copy current file handle to ncoutn
    ncgrid_nx = numxgridn
    ncgrid_ny = numygridn
    ncgrid_nelon = outgridn_nelon
    ncgrid_nelat = outgridn_nelat
    ncgrid_swlon = outgridn_swlon
    ncgrid_swlat = outgridn_swlat
    allocate(ncgrid_oro(ncgrid_nx,ncgrid_ny),stat=stat)
    allocate(ncgrid_area(ncgrid_nx,ncgrid_ny),stat=stat)
    ncgrid_oro   = orooutn(0:ncgrid_nx-1,0:ncgrid_ny-1)
    ncgrid_area  = arean(0:ncgrid_nx-1,0:ncgrid_ny-1)
    if (outgrid_option.eq.1) then ! input was in latlon
      ncgrid_dx = dxoutln
      ncgrid_dy = dyoutln
      ncgrid_lon0 = outlon0n
      ncgrid_lat0 = outlat0n
    else  ! input was in metres
      ncgrid_dx = dxoutn
      ncgrid_dy = dyoutn
      ncgrid_xm0  = out_xm0n
      ncgrid_ym0  = out_ym0n
    endif
  endif

  if (option_verbose.ge.10) &
    write(*,*) 'write_ncheader: ncout,ncoutn=',ncout,ncoutn

  ! Write global attributes
  !*****************************
  call nc_write_global_attributes(ncid,nesting_level)

  ! Set up netcdf dimensions
  !*************************
  if (option_verbose.ge.10) write(*,10) 'main grid dimensions'

  ncret = NF90_DEF_DIM(ncid,'Time',NF90_UNLIMITED,ncrecid)
  call check_ncerror(ncret)

  ncret = NF90_DEF_DIM(ncid,'DateStrLen',15,ncstr3id) !TODO: WRF format
  call check_ncerror(ncret)

  ncret = NF90_DEF_DIM(ncid,'west_east',ncgrid_nx,nclonid)
  call check_ncerror(ncret)

  ncret = NF90_DEF_DIM(ncid,'south_north',ncgrid_ny,nclatid)
  call check_ncerror(ncret)

  ncret = NF90_DEF_DIM(ncid,'bottom_top',numzgrid,nclvlid)
  call check_ncerror(ncret)

  ncret = NF90_DEF_DIM(ncid,'species',nspec,ncspcid)
  call check_ncerror(ncret)

  ncret = NF90_DEF_DIM(ncid,'SpeciesStrLen',10,ncstr1id)
  call check_ncerror(ncret)

  ncret = NF90_DEF_DIM(ncid,'ageclass',nageclass,ncageid)
  call check_ncerror(ncret)

  if (option_verbose.ge.10) write(*,10) 'release point dimensions'
  ncret = NF90_DEF_DIM(ncid,'releases',numpoint,ncrelid)
  call check_ncerror(ncret)
  
  ncret = NF90_DEF_DIM(ncid,'ReleaseStrLen',45,ncstr2id)
  call check_ncerror(ncret)

  ncret = NF90_DEF_DIM(ncid,'ReleaseStartEnd',2,ncrseid)
  call check_ncerror(ncret)

  if(option_verbose.ge.10) write(*,10) 'receptor point dimensions'
  ncret = NF90_DEF_DIM(ncid,'receptors',numreceptor,ncrepid)
  call check_ncerror(ncret)
  
  ncret = NF90_DEF_DIM(ncid,'ReceptorStrLen',16,ncresid)
  call check_ncerror(ncret)

  ! Select which dimensions to use for main output grids
  ncdimsid3(1) = nclonid ! X
  ncdimsid3(2) = nclatid ! Y
  ncdimsid3(3) = nclvlid ! Z
  if (ldirect.eq.1) ncdimsid3(4) = ncspcid ! species
  if (ldirect.eq.-1) ncdimsid3(4) = ncrelid ! points
  ncdimsid3(5) = ncageid ! ageclass
  ncdimsid3(6) = ncrecid ! t

  ncdimsid2(1) = nclonid ! X
  ncdimsid2(2) = nclatid ! Y
  if (ldirect.eq.1) ncdimsid2(3) = ncspcid ! species
  if (ldirect.eq.-1) ncdimsid2(3) = ncrelid ! points
  ncdimsid2(4) = ncageid ! ageclass
  ncdimsid2(5) = ncrecid ! t

  ! Set up dimension variables
  !***************************

  ! XLONG
  if (option_verbose.ge.10) write(*,10) 'XLONG dimension variable'
#ifdef NETCDF4_OUTPUT
  ncret = NF90_DEF_VAR(ncid,vname_x,NF90_FLOAT,ncdimsid2(1:2),nclonvid, &
      shuffle = shuffle, deflate_level=deflate_level)
#else
  ncret = NF90_DEF_VAR(ncid,vname_x,NF90_FLOAT,ncdimsid2(1:2),nclonvid)
#endif
  call check_ncerror(ncret)
  ncret = NF90_PUT_ATT(ncid,nclonvid,descr, &
      'Longitude of center grid, west is negative')
  call check_ncerror(ncret)
  ncret = NF90_PUT_ATT(ncid,nclonvid,units,'degree_east')
  call check_ncerror(ncret)
  
  if (option_verbose.ge.10) write(*,10) 'XLONG CORNER variable'
#ifdef NETCDF4_OUTPUT
  ncret = NF90_DEF_VAR(ncid,'XLONG_CORNER',NF90_FLOAT,ncdimsid2(1:2),nclonvid2, &
      shuffle = shuffle, deflate_level=deflate_level)
#else
  ncret = NF90_DEF_VAR(ncid,'XLONG_CORNER',NF90_FLOAT,ncdimsid2(1:2),nclonvid2)
#endif
  call check_ncerror(ncret)
  ncret = NF90_PUT_ATT(ncid,nclonvid2,descr, &
      'Longitude of lower left corner of grids, west is negative')
  call check_ncerror(ncret)
  ncret = NF90_PUT_ATT(ncid,nclonvid2,units,'degree_east')
  call check_ncerror(ncret)

  ! XLAT
  if (option_verbose.ge.10) write(*,10) 'XLAT dimension variable'
  write(*,*) vname_y
  write(*,*) vname_x
#ifdef NETCDF4_OUTPUT
    ncret = NF90_DEF_VAR(ncid,vname_y,NF90_FLOAT,ncdimsid2(1:2),nclatvid, &
        shuffle = shuffle, deflate_level=deflate_level)
#else
    ncret = NF90_DEF_VAR(ncid,vname_y,NF90_FLOAT,ncdimsid2(1:2),nclatvid)
#endif
  call check_ncerror(ncret)
  ncret = NF90_PUT_ATT(ncid,nclatvid,descr,'Latitude of center grid, south is negative')
  call check_ncerror(ncret)
  ncret = NF90_PUT_ATT(ncid,nclatvid,units,'degree_north')
  call check_ncerror(ncret)
  
  if (option_verbose.ge.10) write(*,10) 'XLAT_CORNER variable'
#ifdef NETCDF4_OUTPUT
  ncret = NF90_DEF_VAR(ncid,'XLAT_CORNER',NF90_FLOAT,ncdimsid2(1:2),nclatvid2, &
      shuffle = shuffle, deflate_level=deflate_level)
#else
  ncret = NF90_DEF_VAR(ncid,'XLAT_CORNER',NF90_FLOAT,ncdimsid2(1:2),nclatvid2)
#endif
  call check_ncerror(ncret)
  ncret = NF90_PUT_ATT(ncid,nclatvid2,descr, &
      'Latitude of lower left corner of grids, south is negative')
  call check_ncerror(ncret)
  ncret = NF90_PUT_ATT(ncid,nclatvid2,units,'degree_north')
  call check_ncerror(ncret)

  ! ZTOP
  if (option_verbose.ge.10) write(*,10) 'ZTOP dimension variable'
#ifdef NETCDF4_OUTPUT
  ncret = NF90_DEF_VAR(ncid,'ZTOP',NF90_FLOAT,ncdimsid3(3),nclvlvid, &
      shuffle = shuffle, deflate_level=deflate_level)
#else
    ncret = NF90_DEF_VAR(ncid,'ZTOP',NF90_FLOAT,ncdimsid3(3),nclvlvid)
#endif
  call check_ncerror(ncret)
  ncret = NF90_PUT_ATT(ncid,nclvlvid,descr, &
      'UPPER BOUNDARY OF MODEL LAYER')
  call check_ncerror(ncret)
  ncret = NF90_PUT_ATT(ncid,nclvlvid,units,'m')
  call check_ncerror(ncret)

  ! SPECIES
  if (option_verbose.ge.10) write(*,10) 'SPECIES dimension variable'
#ifdef NETCDF4_OUTPUT
    ncret = NF90_DEF_VAR(ncid,'SPECIES',NF90_CHAR,(/ncstr1id,ncspcid/),ncspcvid, &
        shuffle = shuffle, deflate_level=deflate_level)
#else
    ncret = NF90_DEF_VAR(ncid,'SPECIES',NF90_CHAR,(/ncstr1id,ncspcid/),ncspcvid)
#endif
  call check_ncerror(ncret)
  ncret = NF90_PUT_ATT(ncid,ncspcvid,descr,'NAME OF SPECIES')
  call check_ncerror(ncret)

  ! AGECLASSES
  if (option_verbose.ge.10) write(*,10) 'AGECLASSES dimension variable'
#ifdef NETCDF4_OUTPUT
    ncret = NF90_DEF_VAR(ncid,'AGECLASS',NF90_INT,ncageid,ncagevid, &
        shuffle = shuffle, deflate_level=deflate_level)
#else
    ncret = NF90_DEF_VAR(ncid,'AGECLASS',NF90_INT,ncageid,ncagevid)
#endif
  call check_ncerror(ncret)
  ncret = NF90_PUT_ATT(ncid,ncagevid,descr,'MAX AGE OF SPECIES IN CLASS')
  call check_ncerror(ncret)
  ncret = NF90_PUT_ATT(ncid,ncagevid,units,'s')
  call check_ncerror(ncret)

  ! TIMES
  if (option_verbose.ge.10) write(*,10) 'TIMES dimension variable'
#ifdef NETCDF4_OUTPUT
    ncret = NF90_DEF_VAR(ncid,vname_t,NF90_CHAR,(/ncstr3id,ncrecid/),ncrecvid, &
        shuffle = shuffle, deflate_level=deflate_level)
#else
    ncret = NF90_DEF_VAR(ncid,vname_t,NF90_CHAR,(/ncstr3id,ncrecid/),ncrecvid)
#endif
  call check_ncerror(ncret)
  ncret = NF90_PUT_ATT(ncid,ncrecvid,descr, &
    'TIME OF OUTPUT (END OF AVERAGING INTERVAL)')

  ! Release related variables
  if (option_verbose.ge.10) write(*,10) 'ReleaseName variable'
#ifdef NETCDF4_OUTPUT
    ncret = NF90_DEF_VAR(ncid,'ReleaseName', &
        NF90_CHAR,(/ncstr2id,ncrelid/),ncrnvid, &
        shuffle = shuffle, deflate_level=deflate_level)
#else
    ncret = NF90_DEF_VAR(ncid,'ReleaseName', &
        NF90_CHAR,(/ncstr2id,ncrelid/),ncrnvid)
#endif
  call check_ncerror(ncret)
  ncret = NF90_PUT_ATT(ncid,ncrnvid,descr,'RELEASE IDENTIFIER/COMMENT')
  call check_ncerror(ncret)
  ncret = NF90_PUT_ATT(ncid,ncrnvid,units,'-')
  call check_ncerror(ncret)

  if (option_verbose.ge.10) write(*,10) 'ReleaseTstart_end variable'
#ifdef NETCDF4_OUTPUT
    ncret = NF90_DEF_VAR(ncid,'ReleaseTstart_end', &
      NF90_INT,(/ncrseid,ncrelid/),ncrtvid, &
      shuffle = shuffle, deflate_level=deflate_level)
#else
    ncret = NF90_DEF_VAR(ncid,'ReleaseTstart_end', &
      NF90_INT,(/ncrseid,ncrelid/),ncrtvid)
#endif
  call check_ncerror(ncret)
  ncret = NF90_PUT_ATT(ncid,ncrtvid,descr, &
    'BEGINNING/ENDING TIME OF RELEASE (SECONDS SINCE RUN START)')
  call check_ncerror(ncret)
  ncret = NF90_PUT_ATT(ncid,ncrtvid,units,'s')
  call check_ncerror(ncret)

  if (option_verbose.ge.10) write(*,10) 'ReleaseXstart_end variable'
#ifdef NETCDF4_OUTPUT
    ncret = NF90_DEF_VAR(ncid,'ReleaseXstart_end',  &
        NF90_FLOAT,(/ncrseid,ncrelid/),ncrxvid, &
        shuffle = shuffle, deflate_level=deflate_level)
#else
    ncret = NF90_DEF_VAR(ncid,'ReleaseXstart_end',  &
        NF90_FLOAT,(/ncrseid,ncrelid/),ncrxvid)
#endif
  call check_ncerror(ncret)
  ncret = NF90_PUT_ATT(ncid,ncrxvid,descr, &
    'WEST/EAST BOUNDARIES OF SOURCE')
  call check_ncerror(ncret)
  ncret = NF90_PUT_ATT(ncid,ncrxvid,units,'degree_north')
  call check_ncerror(ncret)

  if (option_verbose.ge.10) write(*,10) 'ReleaseYstart_end variable'
#ifdef NETCDF4_OUTPUT
    ncret = NF90_DEF_VAR(ncid,'ReleaseYstart_end',  &
      NF90_FLOAT,(/ncrseid,ncrelid/),ncryvid, &
      shuffle = shuffle, deflate_level=deflate_level)
#else
    ncret = NF90_DEF_VAR(ncid,'ReleaseYstart_end',  &
      NF90_FLOAT,(/ncrseid,ncrelid/),ncryvid)
#endif
  call check_ncerror(ncret)
  ncret = NF90_PUT_ATT(ncid,ncryvid,descr, &
    'SOUTH/NORTH BOUNDARIES OF SOURCE')
  call check_ncerror(ncret)
  ncret = NF90_PUT_ATT(ncid,ncryvid,units,'degree_north')
  call check_ncerror(ncret)

  if (option_verbose.ge.10) write(*,10) 'ReleaseZstart_end variable'
#ifdef NETCDF4_OUTPUT
  ncret = NF90_DEF_VAR(ncid,'ReleaseZstart_end',  &
    NF90_FLOAT,(/ncrseid,ncrelid/),ncrzvid, &
    shuffle = shuffle, deflate_level=deflate_level)
#else
    ncret = NF90_DEF_VAR(ncid,'ReleaseZstart_end',  &
    NF90_FLOAT,(/ncrseid,ncrelid/),ncrzvid)
#endif
  call check_ncerror(ncret)
  ncret = NF90_PUT_ATT(ncid,ncrzvid,descr, &
    'BOTTOM/TOP BOUNDARIES OF SOURCE')
  call check_ncerror(ncret)
  ncret = NF90_PUT_ATT(ncid,ncrzvid,units,'m')
  call check_ncerror(ncret)

  if (option_verbose.ge.10) write(*,10) 'ReleaseNP variable'
#ifdef NETCDF4_OUTPUT
    ncret = NF90_DEF_VAR(ncid,'ReleaseNP',NF90_INT,ncrelid,ncspvid, &
        shuffle = shuffle, deflate_level=deflate_level)
#else
    ncret = NF90_DEF_VAR(ncid,'ReleaseNP',NF90_INT,ncrelid,ncspvid)
#endif
  call check_ncerror(ncret)
  ncret = NF90_PUT_ATT(ncid,ncspvid,descr, &
    'TOTAL NUMBER OF PARTICLES RELEASED')
  ncret = NF90_PUT_ATT(ncid,ncspvid,units,'-')
  call check_ncerror(ncret)

  if (option_verbose.ge.10) write(*,10) 'ReleaseXMass variable'
#ifdef NETCDF4_OUTPUT
    ncret = NF90_DEF_VAR(ncid,'ReleaseXMass', &
        NF90_FLOAT,(/ncspcid,ncrelid/),ncrmvid, &
        shuffle = shuffle, deflate_level=deflate_level)
#else
    ncret = NF90_DEF_VAR(ncid,'ReleaseXMass', &
        NF90_FLOAT,(/ncspcid,ncrelid/),ncrmvid)
#endif
  call check_ncerror(ncret)
  ncret = NF90_PUT_ATT(ncid,ncrmvid,descr,'TOTAL MASS RELEASED')
  call check_ncerror(ncret)
  ncret = NF90_PUT_ATT(ncid,ncrmvid,units,'kg')
  call check_ncerror(ncret)

  if (option_verbose.ge.10) write(*,10) 'ReceptorLon variable'
#ifdef NETCDF4_OUTPUT
    ncret = NF90_DEF_VAR(ncid,'ReceptorLon',NF90_FLOAT,ncrepid,ncrexvid, &
        shuffle = shuffle, deflate_level=deflate_level)
#else
    ncret = NF90_DEF_VAR(ncid,'ReceptorLon',NF90_FLOAT,ncrepid,ncrexvid)
#endif
  call check_ncerror(ncret)
  ncret = NF90_PUT_ATT(ncid,ncrexvid,descr,'Longitude of receptor points')
  call check_ncerror(ncret)
  ncret = NF90_PUT_ATT(ncid,ncrexvid,units,'degree_east')
  call check_ncerror(ncret)
  
  if (option_verbose.ge.10) write(*,10) 'ReceptorLat variable'
#ifdef NETCDF4_OUTPUT
    ncret = NF90_DEF_VAR(ncid,'ReceptorLat',NF90_FLOAT,ncrepid,ncreyvid, &
        shuffle = shuffle, deflate_level=deflate_level)
#else
    ncret = NF90_DEF_VAR(ncid,'ReceptorLat',NF90_FLOAT,ncrepid,ncreyvid)
#endif
  call check_ncerror(ncret)
  ncret = NF90_PUT_ATT(ncid,ncreyvid,descr,'Latitude of receptor points')
  call check_ncerror(ncret)
  ncret = NF90_PUT_ATT(ncid,ncreyvid,units,'degree_north')
  call check_ncerror(ncret)

  if (option_verbose.ge.10) write(*,10) 'ReceptorName variable'
#ifdef NETCDF4_OUTPUT
    ncret = NF90_DEF_VAR(ncid,'ReceptorName', &
        NF90_CHAR,(/ncresid,ncrepid/),ncrenvid, &
        shuffle = shuffle, deflate_level=deflate_level)
#else
    ncret = NF90_DEF_VAR(ncid,'ReceptorName', &
        NF90_CHAR,(/ncresid,ncrepid/),ncrenvid)
#endif
  call check_ncerror(ncret)
  ncret = NF90_PUT_ATT(ncid,ncrenvid,descr,'RECEPTOR IDENTIFIER/COMMENT')
  call check_ncerror(ncret)
  ncret = NF90_PUT_ATT(ncid,ncrenvid,units,'-')
  call check_ncerror(ncret)

  ! Since we need to exit define mode before we can insert
  ! variable data, we will include the last file attributes and
  ! define the last variables here.

  ! DIRECTION INDEPENDENT OUTPUT VARIABLES
  if (option_verbose.ge.10) write(*,10) 'TOPOGRAPHY variable'
#ifdef NETCDF4_OUTPUT
    ncret = NF90_DEF_VAR(ncid,'TOPOGRAPHY', &
        NF90_FLOAT,ncdimsid2(1:2),nctovid, &
        shuffle = shuffle, deflate_level=deflate_level)
#else
    ncret = NF90_DEF_VAR(ncid,'TOPOGRAPHY', &
        NF90_FLOAT,ncdimsid2(1:2),nctovid)
#endif
  call check_ncerror(ncret)
  ncret = NF90_PUT_ATT(ncid,nctovid,descr,  &
      'TERRAIN ELEVATION ABOVE SEA LEVEL')
  call check_ncerror(ncret)
  ncret = NF90_PUT_ATT(ncid,nctovid,units,'m')
  call check_ncerror(ncret)
  ncret = NF90_PUT_ATT(ncid,nctovid,coord,coordxy)
  call check_ncerror(ncret)

  if (option_verbose.ge.10) write(*,10) 'GRIDAREA variable'
#ifdef NETCDF4_OUTPUT
  ncret = NF90_DEF_VAR(ncid,'GRIDAREA', &
      NF90_FLOAT,ncdimsid2(1:2),ncarvid, &
      shuffle = shuffle, deflate_level=deflate_level)
#else
  ncret = NF90_DEF_VAR(ncid,'GRIDAREA', &
      NF90_FLOAT,ncdimsid2(1:2),ncarvid)
#endif
  call check_ncerror(ncret)
  ncret = NF90_PUT_ATT(ncid,ncarvid,descr, &
    'SURFACE AREA OF EACH GRID CELL')
  call check_ncerror(ncret)
  ncret = NF90_PUT_ATT(ncid,ncarvid,units,'m2')
  call check_ncerror(ncret)
  ncret = NF90_PUT_ATT(ncid,ncarvid,coord,coordxy)
  call check_ncerror(ncret)

  ! EXIT DEFINE MODE, ENTER DATA MODE
  ncret = NF90_ENDDEF(ncid)
  call check_ncerror(ncret)

  ! DIMENSION VARIABLES
  if (option_verbose.ge.10) write(*,10) 'ZTOP data'
  !ncret = nf_put_var_real(ncid,nclvlvid,outheight)
  ncret = NF90_PUT_VAR(ncid,nclvlvid,outheight)
  call check_ncerror(ncret)

  ! X,Y - Lon,Lat
  if (option_verbose.ge.10) write(*,10) 'XLAT,XLONG data'

  if (outgrid_option.eq.0) then ! irregular
    do jy=1,ncgrid_ny
    xyvec(2)=jy
      do ix=1,ncgrid_nx
        xyvec(1)=ix
        tmpx=ncgrid_xm0+(float(ix)-0.5)*ncgrid_dx
        tmpy=ncgrid_ym0+(float(jy)-0.5)*ncgrid_dy
        call xymeter_to_ll_wrf(tmpx,tmpy,tmplon,tmplat)
        !ncret = nf_put_vara_real(ncid,nclonvid,xyvec,ones2,tmplon)
        ncret = NF90_PUT_VAR(ncid,nclonvid,tmplon,start=xyvec)
        call check_ncerror(ncret)
        !ncret = nf_put_vara_real(ncid,nclatvid,xyvec,ones2,tmplat)
        ncret = NF90_PUT_VAR(ncid,nclatvid,tmplat,start=xyvec)
        call check_ncerror(ncret)
        tmpx=ncgrid_xm0+(float(ix)-1.)*ncgrid_dx
        tmpy=ncgrid_ym0+(float(jy)-1.)*ncgrid_dy
        call xymeter_to_ll_wrf(tmpx,tmpy,tmplon,tmplat)
        !ncret = nf_put_vara_real(ncid,nclonvid2,xyvec,ones2,tmplon)
        ncret = NF90_PUT_VAR(ncid,nclonvid2,tmplon,start=xyvec)
        call check_ncerror(ncret)
        !ncret = nf_put_vara_real(ncid,nclatvid2,xyvec,ones2,tmplat)
        ncret = NF90_PUT_VAR(ncid,nclatvid2,tmplat,start=xyvec)
        call check_ncerror(ncret)
      enddo
    enddo
  else
    do jy=1,ncgrid_ny
      xyvec(2)=jy
      do ix=1,ncgrid_nx
        xyvec(1)=ix
        call ll_to_xymeter_wrf(ncgrid_swlon,ncgrid_swlat,xsw,ysw)
        call ll_to_xymeter_wrf(ncgrid_nelon,ncgrid_nelat,xne,yne)
        tmpx=xsw+(xne-xsw)*float(ix-1)/float(ncgrid_nx-1)
        tmpy=ysw+(yne-ysw)*float(jy-1)/float(ncgrid_ny-1)
        call xymeter_to_ll_wrf(tmpx,tmpy,tmplon,tmplat)
        xl2=ncgrid_lon0+(float(ix)-0.5)*dxoutl !long
        yl2=ncgrid_lat0+(float(jy)-0.5)*dyoutl !lat
        !ncret = nf_put_vara_real(ncid,nclonvid,xyvec,ones2,xl2)
        ncret = NF90_PUT_VAR(ncid,nclonvid,xl2,start=xyvec)
        call check_ncerror(ncret)
        !ncret = nf_put_vara_real(ncid,nclatvid,xyvec,ones2,yl2)
        ncret = NF90_PUT_VAR(ncid,nclatvid,yl2,start=xyvec)
        call check_ncerror(ncret)
        xl2=ncgrid_lon0+(float(ix)-1.)*dxoutl !long
        yl2=ncgrid_lat0+(float(jy)-1.)*dyoutl !lat
        !ncret = nf_put_vara_real(ncid,nclonvid2,xyvec,ones2,xl2)
        ncret = NF90_PUT_VAR(ncid,nclonvid2,xl2,start=xyvec)
        call check_ncerror(ncret)
        !ncret = nf_put_vara_real(ncid,nclatvid2,xyvec,ones2,yl2)
        ncret = NF90_PUT_VAR(ncid,nclatvid2,yl2,start=xyvec)
        call check_ncerror(ncret)
      enddo
    enddo
  endif ! outgrid_option


  ! Write information on release points: total number, then for each point:
  ! start, end, coordinates, # of particles, name, mass
  !************************************************************************
  do i=1,numpoint
    xp1=xpoint1(i)*dx+xlon0 ! This is probably wrong, but it seems to be 
    yp1=ypoint1(i)*dy+ylat0 ! the same in writeheader*.f90, so I'll leave
    xp2=xpoint2(i)*dx+xlon0 ! it for now... //AD
    yp2=ypoint2(i)*dy+ylat0 !

    if (option_verbose.ge.10) write(*,10) 'ReleaseTstart_end data'
    !ncret = nf_put_vara_int(ncid,ncrtvid,   & ! ReleaseTstart_end
    !  (/1,i/),(/2,1/),(/ireleasestart(i),ireleaseend(i)/))
    ncret = NF90_PUT_VAR(ncid,ncrtvid, &        ! ReleaseTstart_end
        (/ireleasestart(i),ireleaseend(i)/), &
        start=(/1,i/), count=(/2,1/) )
    call check_ncerror(ncret)

    if (option_verbose.ge.10) write(*,10) 'ReleaseXstart_end data'
    !ncret = nf_put_vara_real(ncid,ncrxvid,  & ! ReleaseXstart_end
    !  (/1,i/),(/2,1/),(/xp1,xp2/))
    ncret = NF90_PUT_VAR(ncid,ncrxvid, &        ! ReleaseXstart_end
        (/xp1,xp2/), &
        start=(/1,i/), count=(/2,1/) )
    call check_ncerror(ncret)

    if (option_verbose.ge.10) write(*,10) 'ReleaseYstart_end data'
    !ncret = nf_put_vara_real(ncid,ncryvid,  & !ReleaseYstart_end
    !  (/1,i/),(/2,1/),(/yp1,yp2/))
    ncret = NF90_PUT_VAR(ncid,ncryvid, &        ! ReleaseYstart_end
        (/yp1,yp2/), &
        start=(/1,i/), count=(/2,1/) )
    call check_ncerror(ncret)

    if (option_verbose.ge.10) write(*,10) 'ReleaseZstart_end data'
    !ncret = nf_put_vara_real(ncid,ncrzvid,  & !ReleaseZstart_end
    !  (/1,i/),(/2,1/),(/zpoint1(i),zpoint2(i)/))
    ncret = NF90_PUT_VAR(ncid,ncrzvid, &        ! ReleaseZstart_end
        (/zpoint1,zpoint2/), &
        start=(/1,i/), count=(/2,1/) )
    call check_ncerror(ncret)

    if (option_verbose.ge.10) write(*,10) 'ReleaseXMass data'
    !ncret = nf_put_vara_real(ncid,ncrmvid,  & !ReleaseXMass
    !  (/1,i/),(/nspec,1/),xmass(i,1:nspec))
    ncret = NF90_PUT_VAR(ncid,ncrmvid, &        ! ReleaseXMass
        xmass(i,1:nspec), &
        start=(/1,i/), count=(/nspec,1/) )
    call check_ncerror(ncret)

    if (option_verbose.ge.10) write(*,10) 'ReleaseNP data'
    !ncret = nf_put_vara_int(ncid,ncspvid,   & !ReleaseNP
    !  i,1,npart(i))
    ncret = NF90_PUT_VAR(ncid,ncspvid, npart(:)) !ReleaseNP
    ! For some reason I cannot get this to work with indices, they aren't need anyway..
    call check_ncerror(ncret)

    !Release Name/Comment
    j=1 ! Find the length of each release point comment/name
    do while( j.lt.45.and.compoint(i)(j+1:j+1).ne." ")
      j=j+1
    enddo
    if (option_verbose.ge.10) write(*,10) 'ReleaseName data'
    !ncret = nf_put_vara_text(ncid,ncrnvid,(/1,i/),(/j,1/),compoint(i)(1:j))
    ncret = NF90_PUT_VAR(ncid,ncrnvid,compoint(i)(1:j),start=(/1,i/),count=(/j,1/))
    call check_ncerror(ncret)
  enddo

  ! Get output coordinates in units of lat,lon
  do j = 1, numreceptor
    x_rep_tmp(j) = xreceptor(j)*dx + xmet0
    y_rep_tmp(j) = yreceptor(j)*dy + ymet0
    !if (outgrid_option .eq. 1) then  ! AD: it makes more sense to always use lat/lon in netcdf
    x_rep_tmpb = x_rep_tmp(j)
    y_rep_tmpb = y_rep_tmp(j)
    call xymeter_to_ll_wrf( x_rep_tmpb, y_rep_tmpb, x_rep_tmp(j), y_rep_tmp(j) )
    !endif
  enddo
  
  if (option_verbose.ge.10) write(*,10) 'ReceptorLon data'
  !ncret = nf_put_vara_real(ncid,ncrexvid, & !ReceptorLon
  !  1,numreceptor,x_rep_tmp(1:numreceptor))
  ncret = NF90_PUT_VAR(ncid,ncrexvid, &     ! ReceptorLon
      x_rep_tmp(1:numreceptor) )
  call check_ncerror(ncret)
  
  if (option_verbose.ge.10) write(*,10) 'ReceptorLat data'
  !ncret = nf_put_vara_real(ncid,ncreyvid, & ! ReceptorLat
  !  1,numreceptor,y_rep_tmp(1:numreceptor))
  ncret = NF90_PUT_VAR(ncid,ncreyvid, &     ! ReceptorLat
      y_rep_tmp(1:numreceptor) )
  call check_ncerror(ncret)
  
  if (option_verbose.ge.10) write(*,10) 'ReceptorName data'
  do i=1,numreceptor
    j=1
    do while( j.lt.16 .and. receptorname(i)(j+1:j+1).ne." ")
      j=j+1
    enddo
    !ncret = nf_put_vara_text(ncid,ncrenvid,(/1,i/),(/j,1/),receptorname(i)(1:j))
    ncret = NF90_PUT_VAR(ncid,ncrenvid,&
        receptorname(i)(1:j), &
        start=(/1,i/), count=(/j,1/) )
    call check_ncerror(ncret)
  enddo

  ! Write age class information
  !****************************
  if (option_verbose.ge.10) write(*,10) 'AGECLASSES data'
  !ncret = nf_put_var_int(ncid,ncagevid,lage(1:nageclass))
  ncret = NF90_PUT_VAR(ncid,ncagevid,lage(1:nageclass))
  call check_ncerror(ncret)

  ! Write topography to output file
  !********************************
  if (option_verbose.ge.10) write(*,10) 'TOPOGRAPHY data'
! do ix=0,ncgrid_nx-1
  do ix=1,ncgrid_nx
  !  ncret = nf_put_vara_real(ncid,nctovid,  &
  !    (/ix,1/),(/1,ncgrid_ny/),ncgrid_oro(ix,1:ncgrid_ny))
    ncret = NF90_PUT_VAR(ncid,nctovid, &
        ncgrid_oro(ix,1:ncgrid_ny), &
        start=(/ix,1/), count=(/1,ncgrid_ny/) )
    call check_ncerror(ncret)
  enddo

  ! Write grid cell surface area
  !*****************************
  if (option_verbose.ge.10) write(*,10) 'GRIDAREA data'
  do ix=1,ncgrid_nx
    !ncret = nf_put_vara_real(ncid,ncarvid,  &
    !  (/ix,1/),(/1,ncgrid_ny/),ncgrid_area(ix,1:ncgrid_ny))
    ncret = NF90_PUT_VAR(ncid,ncarvid, &
        ncgrid_area(ix,1:ncgrid_ny), &
        start=(/ix,1/), count=(/1,ncgrid_ny/) )
    call check_ncerror(ncret)
  enddo

  ! SAVE CREATED NETCDF TO FILE
  !****************************
  if (option_verbose.ge.1) write(*,*) 'nc_create_header_outfile: writing to disk'
  !ncret = nf_sync(ncid)
  ncret = NF90_SYNC(ncid)
  call check_ncerror(ncret)

  return

10 format('nc_create_header_outfile: Setting up ',A)

    !ncret=nf_close(ncid)
    ncret=NF90_CLOSE(ncid)
    deallocate(ncgrid_oro,ncgrid_area)
end subroutine nc_create_header_outfile

  subroutine nc_write_output(itime,outnum,ks,kp,nage,tot_mu_scalar,rho_recept,grid_vol,nesting_level)
  !*****************************************************************************
  !                                                                            *
  !  This routine writes concentration, mixing ratio and deposition fields     *
  !  to a netcdf file defined by nc_create_main_outfile.                       *
  !                                                                            *
  !  flex_ncheader is called from within write_ncconc when it's time for a new *
  !  output file.                                                              *
  !                                                                            *
  !  write_ncconc should be called by concoutput_irreg and concoutput_reg      *
  !  it is separate from the binary and ascii output routines to avoid mixing  *
  !  of sparse and full grid approaches.  Netcdf will output the full grid.    *
  !                                                                            *
  !      Author: A. Dingwell                                                   *
  !                                                                            *
  !      29 May 2013                                                           *
  !                                                                            *
  ! Modifications:                                                             *
  ! June 5 2013: J. Brioude: compression using deflate level, optimization of  *
  !  the writing procedure. bug fixes for backtrajectory mode                  *
  !  Feb 2014:   A. Dingwell: bug fix                                          *
  !  2015-05-04: A. Dingwell: Added receptor points to netcdf output           *
  !  2015-06-04: A. Dingwell: Updated subroutine to use NetCDF F90 module.     *
  !  2015-06-05: A. Dingwell: moved to module                                  *
  !  2015-07-02: A. Dingwell: Added grid_vol to list of input variables. Could *
  !                 have been avoided with a pointer, but it's tricky (and not *
  !                 very safe) to target an allocatable variable...            *
  !*****************************************************************************
    real,intent(in)    :: outnum         ! Number of samples for each concentration calculation
    integer,intent(in) :: itime          ! Current simulation time [s]
    integer,intent(in) :: ks,kp,nage     ! species, maxpointspec_act and nageclass indices resp.
    real,intent(in)    :: tot_mu_scalar  ! total mass per source and species (backward)
                              ! or unity (forward).  Should probably be sent as
                              ! tot_mu(ks,kp) from concoutput*.f90
    real,intent(in)     :: rho_recept(maxreceptor)     ! Air density at each receptor point (array)
    real,intent(in)     :: grid_vol(:,:,:)  ! Volume of each grid cell (for current nest).
    integer,intent(in) :: nesting_level  ! 0 for main (parent) grid, 1 for nest (child)

    real(kind=dp) :: jul          ! Julian date
    integer   :: jjjjmmdd,ihmmss  ! date & time as integer
    character :: adate*8,atime*6  ! date and time strings, used for filename

    integer :: ncid           ! Pointer to netcdf file, depends on nesting level
    integer :: grid_nx,grid_ny! outgrid dimensions, depend on the nesting level
    integer :: ncret          ! Netcdf:  return code
    integer :: ix,jy,kz       ! iterators
    character :: datestr*15   ! For the Times variable
    integer :: deflate_level=5 ! compression level
      
    if (option_verbose.ge.1) then
      write(*,*)'write_ncconc: writing netcdf output for: domain,kp,nage =',&
        nesting_level+1,kp,nage
    endif

    ! Determine which nest/outfile we are writing to
    !***********************************************
    if (nesting_level.eq.0) then
      ncid    = ncout
      grid_nx = numxgrid
      grid_ny = numygrid
    elseif (nesting_level.eq.1) then
      ncid    = ncoutn
      grid_nx = numxgridn
      grid_ny = numygridn
    else
      write(*,*) '***write_ncconc error: nesting level  must be 0 or 1'
      ! Note for future development: If additional output nests are to be
      ! supported for netcdf output, modification must be made here as well as in
      ! the respective nesting_level if-block in write_ncheader
    endif
    ! Update/Initialize record index
    !*******************************
     if ((ks.eq.1).and.(kp.eq.1).and.(nage.eq.1)) then
  !   print*,'ncirec',ncirec,ncnumrec
    if (nesting_level.eq.0) then  ! Only update for first domain
      if (itime.eq.loutstep) then  ! first output
        ncirec = 1  ! initialize record index
      elseif (ncirec.eq.ncnumrec) then  ! new file
  !      print*,'file is closing'
        ncirec = 1  ! reset record index
        ncret=NF90_CLOSE(ncid)      ! close the old file
        call check_ncerror(ncret)
  !      print*,'file is closed'
      else
        ncirec=ncirec+1 ! move on to next record
      endif
    endif
  !   print*,'ncirec',ncirec,ncnumrec
    endif

    ! Check if it's time to create a new netcdf file
    !***********************************************
    if (ncirec.eq.1) then         ! First output in current file?
  !   write(*,*) 'itime=',itime
      if ((ks.eq.1).and.(kp.eq.1).and.(nage.eq.1)) then
        if (option_verbose.ge.1) &
          write(*,*)'write_ncconc: calling write_ncinfo'
        call nc_create_main_outfile(itime,nesting_level)  ! Create new file
      endif
      ! Reassign file handle to the newly created file:
      if (nesting_level.eq.0) ncid=ncout
      if (nesting_level.eq.1) ncid=ncoutn
    endif

    if (option_verbose.ge.10) &
      write(*,*) 'ncid,nccovid=',ncid,nccovid


    ! Create output for the current record index
    !*******************************************
    jul=bdate+real(itime,kind=dp)/86400._dp
    call caldate(jul,jjjjmmdd,ihmmss)
    write(adate,'(i8.8)') jjjjmmdd
    write(atime,'(i6.6)') ihmmss
    
    if ((ks.eq.1).and.(kp.eq.1).and.(nage.eq.1)) then
      if (option_verbose.ge.10) write(*,*)'write_ncconc: record index',ncirec
      write(datestr,'(I8.8,A1,I6.6)') jjjjmmdd,'_',ihmmss
      ncret = NF90_PUT_VAR(ncid,ncrecvid,datestr,start=(/1,ncirec/),count=(/15,1/))
      call check_ncerror(ncret)
    endif

    if ((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)) then ! concentration
      if (option_verbose.ge.1)  &
        write(*,*)'write_ncconc: concentration output',kp,nage,ncirec,nccovid,ncid
      do kz=1,numzgrid
        do jy=0,grid_ny-1
          do ix=0,grid_nx-1
            grid2(ix,jy,kz,kp,nage)= grid(ix,jy,kz)*factor3d(ix,jy,kz)/tot_mu_scalar
          enddo ! ix=1,grid_nx-1
        enddo ! jy=1,grid_ny-1
      enddo ! kz=1,numzgrid

      if (kp.eq.maxpointspec_act .and. nage.eq.nageclass) then
        if (ldirect.eq.-1) then 
          ncret = NF90_PUT_VAR(ncid,nccovid, &
              grid2(0:grid_nx-1,0:grid_ny-1,1:numzgrid,1:kp,1:nage), &
              start = (/1,1,1,1,1,ncirec/), &
              count = (/grid_nx,grid_ny,numzgrid,kp,nage,1/) )
          call check_ncerror(ncret)
        else
          if (kp.gt.1) then
           ncret = NF90_PUT_VAR(ncid,nccovid, &
              grid2(0:grid_nx-1,0:grid_ny-1,1:numzgrid,1:kp,1:nage), &
              start = (/1,1,1,1,ks,1,ncirec/), &
              count = (/grid_nx,grid_ny,numzgrid,kp,1,nage,1/) ) 
          else
           ncret = NF90_PUT_VAR(ncid,nccovid, &
              grid2(0:grid_nx-1,0:grid_ny-1,1:numzgrid,1,1:nage), &
              start = (/1,1,1,ks,1,ncirec/), &
              count = (/grid_nx,grid_ny,numzgrid,1,nage,1/) )
          endif
          call check_ncerror(ncret)
        endif
      endif
    endif ! concentraion

    if ((iout.eq.2).or.(iout.eq.3)) then  ! mixing ratio
      if (option_verbose.ge.1) &
        write(*,*)'write_ncconc: mixing ratio output',kp,nage,ncirec,ncravid,ncid
      do kz=1,numzgrid
        do jy=0,grid_ny-1
          do ix=0,grid_nx-1
            grid3(ix,jy,kz,kp,nage)= 1.e12*grid(ix,jy,kz)/grid_vol(ix,jy,kz)/outnum*  &
              weightair/weightmolar(ks)/densityoutgrid(ix,jy,kz)
          enddo ! ix=1,grid_nx-1
        enddo ! jy=1,grid_ny-1
      enddo ! kz=1,numzgrid
      if (option_verbose.ge.10) &
        write(*,*)'write_ncconc: ldirect,ks,kp = ',ldirect,ks,kp
      if (kp.eq.maxpointspec_act .and. nage.eq.nageclass) then
        if (ldirect.eq.-1) then
          ncret = NF90_PUT_VAR(ncid,ncravid, &
              grid3(0:grid_nx-1,0:grid_ny-1,1:numzgrid,1:kp,1:nage), &
              start = (/1,1,1,1,1,ncirec/), &
              count = (/grid_nx,grid_ny,numzgrid,kp,nage,1/) )
          call check_ncerror(ncret)
        else
          if (kp.gt.1) then
           ncret = NF90_PUT_VAR(ncid,ncravid, &
              grid3(0:grid_nx-1,0:grid_ny-1,1:numzgrid,1:kp,1:nage), &
              start = (/1,1,1,1,ks,1,ncirec/), &
              count = (/grid_nx,grid_ny,numzgrid,kp,1,nage,1/) )
          else
           ncret = NF90_PUT_VAR(ncid,ncravid, &
              grid3(0:grid_nx-1,0:grid_ny-1,1:numzgrid,1,1:nage), &
              start = (/1,1,1,ks,1,ncirec/), &
              count = (/grid_nx,grid_ny,numzgrid,1,nage,1/) )
          endif
          call check_ncerror(ncret)
        endif
      endif
    endif ! mixing ratio

    if ((ldirect.eq.1).and.(WETDEP)) then ! WETDEP
      if (option_verbose.ge.1)write(*,*)'write_ncconc: wet deposition output'
      do jy=0,grid_ny-1
      do ix=0,grid_nx-1
      if (nesting_level.eq.0)  wetgrid2(ix,jy,kp,nage)=1.e12*wetgrid(ix,jy)/area(ix,jy)
      if (nesting_level.eq.1)  wetgrid2(ix,jy,kp,nage)=1.e12*wetgrid(ix,jy)/arean(ix,jy)
      enddo ! ix=1,grid_nx-1
      enddo ! jy=1,grid_ny-1
      if (kp.eq.maxpointspec_act .and. nage.eq.nageclass) then
    if (ldirect.eq.-1) then
        ncret = NF90_PUT_VAR(ncid,ncwdvid, &
            wetgrid2(0:grid_nx-1,0:grid_ny-1,1:kp,1:nage), &
            start = (/1,1,1,1,ncirec/), &
            count = (/grid_nx,grid_ny,kp,nage,1/) )
        call check_ncerror(ncret)
    else
      if (kp.gt.1) then
        ncret = NF90_PUT_VAR(ncid,ncwdvid, &
            wetgrid2(0:grid_nx-1,0:grid_ny-1,1:kp,1:nage), &
            start = (/1,1,1,ks,1,ncirec/), &
            count = (/grid_nx,grid_ny,kp,1,nage,1/) )
      else
        ncret = NF90_PUT_VAR(ncid,ncwdvid, &
            wetgrid2(0:grid_nx-1,0:grid_ny-1,1,1:nage), &
            start = (/1,1,ks,1,ncirec/), &
            count = (/grid_nx,grid_ny,1,nage,1/) )
      endif
        call check_ncerror(ncret)
      endif
    endif
  !    do jy=0,grid_ny-1
  !    do ix=0,grid_nx-1
  !      ncret = nf_put_vara_real(ncid,ncwdvid, &
  !        (/ix+1,jy+1,kp,nage,ncirec/),(/1,1,1,1,1/), &
  !        1.e12*wetgrid(ix,jy)/area(ix,jy))
  !      call check_ncerror(ncret)
  !    enddo ! ix=1,grid_nx-1
  !    enddo ! jy=1,numygrid-1
    endif ! WETDEP

    if ((ldirect.eq.1).and.(DRYDEP)) then ! DRYDEP
      if (option_verbose.ge.1)write(*,*)'write_ncconc: dry deposition output'
      do jy=0,grid_ny-1
      do ix=0,grid_nx-1
      if (nesting_level.eq.0)  drygrid2(ix,jy,kp,nage)=1.e12*drygrid(ix,jy)/area(ix,jy)
      if (nesting_level.eq.1)  drygrid2(ix,jy,kp,nage)=1.e12*drygrid(ix,jy)/arean(ix,jy)
      enddo ! ix=1,grid_nx-1
      enddo ! jy=1,grid_ny-1
      if (kp.eq.maxpointspec_act .and. nage.eq.nageclass) then
    if (ldirect.eq.-1) then
        ncret = NF90_PUT_VAR(ncid,ncddvid, &
            drygrid2(0:grid_nx-1,0:grid_ny-1,1:kp,1:nage), &
            start = (/1,1,1,1,ncirec/), &
            count = (/grid_nx,grid_ny,kp,nage,1/) )
        call check_ncerror(ncret)
    else
      if (kp.gt.1) then
        ncret = NF90_PUT_VAR(ncid,ncddvid, &
          drygrid2(0:grid_nx-1,0:grid_ny-1,1:kp,1:nage), &
          start = (/1,1,1,ks,1,ncirec/), &
          count = (/grid_nx,grid_ny,kp,1,nage,1/) )
      else
        ncret = NF90_PUT_VAR(ncid,ncddvid, &
          drygrid2(0:grid_nx-1,0:grid_ny-1,1,1:nage), &
          start = (/1,1,ks,1,ncirec/), &
          count = (/grid_nx,grid_ny,1,nage,1/) )
      endif
        call check_ncerror(ncret)
      endif
    endif

  !    do jy=0,grid_ny-1
  !    do ix=0,grid_nx-1
  !      ncret = nf_put_vara_real(ncid,ncddvid, &
  !        (/ix+1,jy+1,kp,nage,ncirec/),(/1,1,1,1,1/), &
  !        1.e12*drygrid(ix,jy)/area(ix,jy))
  !      call check_ncerror(ncret)
  !    enddo ! ix=1,grid_nx-1
  !    enddo ! jy=1,numygrid-1
    endif ! DRYDEP

    if (numreceptor.gt.0) then ! RECEPTOR OUTPUT
      if (option_verbose.ge.10) &
        write(*,*) "write_ncconc: Receptor output, ks =",ks
      if (iout.eq.2 .or. iout.eq.3) then  ! mixing ratio output
        ncret = NF90_PUT_VAR(ncid,ncrmivid, &
          1.e12*creceptor(:,ks)/outnum*weightair/weightmolar(ks)/rho_recept, &
          start = (/1,ks,ncirec/), &
          count = (/numreceptor,1,1/) )
      endif ! mixing ratio output
      
      if ((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)) then ! concentration output
        ncret = NF90_PUT_VAR(ncid,ncrcovid, &
          1.e12*creceptor(:,ks)/outnum, &
          start = (/1,ks,ncirec/), &
          count = (/numreceptor,1,1/) )
             
      endif ! concentration output
    endif ! RECEPTOR OUTPUT


    ncret=NF90_SYNC(ncid)
    call check_ncerror(ncret)
    
    10 format('nc_write_output: Setting up ',A)

  end subroutine nc_write_output

end module netcdf_output_mod
