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

      subroutine readwind_nests_timeav(indj,n,uuhn,vvhn,wwhn)
!**********************************************************************
!                                                                     * 
!             TRAJECTORY MODEL SUBROUTINE READWIND                    *
!                                                                     *
!**********************************************************************
!                                                                     * 
!  April 2012, J. Brioude:                                            *
!             This routine handles the difference in time for         *
!             time-average fields.                                    *
!                                                                     *
!  Feb 2014, A. Griffiths:                                            *
!             bug fix in the dimension of the final loops             *
!                                                                     *
!  2015-03-26, A. Dingwell:                                           *
!             Updated calls to read_ncwrfout_gridinfo to match updates*
!             in that subroutine.                                     * 
!                                                                     *
!**********************************************************************
!                                                                     *
! Note:  This is the FLEXPART_WRF version of subroutine readwind.     *
!    The met fields are read from WRF netcdf output files.            *
!    There are many differences from the FLEXPART version.            *
!                                                                     *
! DESCRIPTION:                                                        *
!                                                                     *
! READING OF ECMWF METEOROLOGICAL FIELDS FROM INPUT DATA FILES. THE   *
! INPUT DATA FILES ARE EXPECTED TO BE AVAILABLE IN GRIB CODE          *
!                                                                     *
! INPUT:                                                              *
! indj               indicates number of the wind field to be read in *
! n                  temporal index for meteorological fields (1 to 3)*
!                                                                     *
! IMPORTANT VARIABLES FROM COMMON BLOCK:                              *
!                                                                     *
! wfname             File name of data to be read in                  *
! nx,ny,nuvz,nwz     expected field dimensions                        *
! nlev_ec            number of "T-grid" vertical levels wwf model     *
!                    (the unstaggered "bottom_top" dimension)         *
! uu,vv,ww           wind fields                                      *
! tt,qv              temperature and specific humidity                *
! ps                 surface pressure                                 *
!                                                                     *
!**********************************************************************
!

!      include 'includepar'
!      include 'includecom'
  use par_mod
  use com_mod
  use flexwrf_ncdf_mod

! subr arguments
      integer :: indj, n

      real(kind=4) :: uuhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
      real(kind=4) :: vvhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
      real(kind=4) :: wwhn(0:nxmaxn-1,0:nymaxn-1,nwzmax,maxnests)
      real(kind=4) :: uuhn2(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
      real(kind=4) :: vvhn2(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
      real(kind=4) :: wwhn2(0:nxmaxn-1,0:nymaxn-1,nwzmax,maxnests)

! local variables
      integer,parameter :: ndims_max=4

      integer :: i, idiagaa, ierr, ifn, itime
      integer :: iduma,indj2,l
      integer :: j, jhhmmss, jyyyymmdd
      integer :: k, kbgn
      integer :: lendim(ndims_max), lendim_exp(ndims_max), &
          lendim_max(ndims_max)
      integer :: levdiff2,deltat,deltat2
      integer :: ndims, ndims_exp
      integer :: n_west_east, n_south_north, n_bottom_top
      integer :: m_grid_id_dum, m_parent_grid_id_dum, &
        m_parent_grid_ratio_dum,  &
        i_parent_start_dum, j_parent_start_dum, &
        map_proj_id_dum,  &
        ext_scalar,pbl_physics,mp_physics_dum,num_land_cat_dum

      real :: dx_met, dy_met
      real :: duma, dumb, dumc, dumd, dume
      real :: dumdz
!      real(kind=4) :: dumarray_aa(nwzmax+1)
!      real(kind=4) :: dumarray_pp(0:nxmax-1,0:nymax-1,nwzmax+1)
      real :: dumarray_aa(nwzmax+1)
      real :: dumarray_pp(0:nxmax-1,0:nymax-1,nwzmax+1)
      real :: ewater_mb, esatwater_mb
      real :: ew      ! this is an external function
      real :: map_stdlon_dum, map_truelat1_dum, map_truelat2_dum
      real :: pint
      real :: toler

!      real(kind=4) :: ewss(0:nxmax-1,0:nymax-1),nsss(0:nxmax-1,0:nymax-1)
      real :: ewss(0:nxmax-1,0:nymax-1),nsss(0:nxmax-1,0:nymax-1)
      real :: plev1,pmean,tv,fu,hlev1,ff10m,fflev1

      real(kind=dp) :: jul,juldate

      character(len=160) :: fnamenc, varname,fnamenc2

      logical :: hflswitch

!
!   get grid info from the wrf netcdf file
!   and check it for consistency against values from gridcheck
!

!
! main loop -- process each nest
!
      do l=1,numbnests

      m = numpath+2*(l-1)+1
      fnamenc = path(m)(1:length(m)) // wfnamen(l,indj)

!      fnamenc = path(2)(1:length(2))//wfname(indj)
      idiagaa = 0

      call read_ncwrfout_gridinfo( ierr, idiagaa, fnamenc, &
        n_west_east, n_south_north, n_bottom_top, & 
        dx_met, dy_met,  &
        m_grid_id_dum, m_parent_grid_id_dum, m_parent_grid_ratio_dum, &
        i_parent_start_dum, j_parent_start_dum, &
        map_proj_id_dum, map_stdlon_dum,  &
        map_truelat1_dum, map_truelat2_dum, &
        ext_scalar,pbl_physics,mp_physics_dum, num_land_cat_dum )
      if (ierr .ne. 0) then
          write(*,9100) 'error getting gridinfor for met file', fnamenc
          stop
      end if

      i_parent_start_dum = i_parent_start_dum-1
      j_parent_start_dum = j_parent_start_dum-1


9100  format( / '*** readwind -- ', a )
9110  format( / '*** readwind -- ', a, 1x, i8 / &
        'file = ', a )
9120  format( / '*** readwind -- ', a, 2(1x,i8) / &
        'file = ', a )
9130  format( / '*** readwind -- ', a, 3(1x,i8) / &
        'file = ', a )
9115  format( / '*** readwind -- ', a / a, 1x, i8 / &
        'file = ', a )
9125  format( / '*** readwind -- ', a / a, 2(1x,i8) / &
        'file = ', a )
9135  format( / '*** readwind -- ', a / a, 3(1x,i8) / &
        'file = ', a )


! READ THE FIRST FILE
       if (ldirect.eq.1) indj2=indj
       if (ldirect.eq.-1) indj2=indj
       deltat=wfdt(indj2)
      fnamenc = path(m)(1:length(m)) // wfnamen(l,indj2)
! locate the date/time in the file
      itime = 0
1100  itime = itime + 1
      call read_ncwrfout_1datetime( ierr, fnamenc, &
          itime, jyyyymmdd, jhhmmss )
      if (ierr .eq. -1) then
          write(*,9100) 'error reading time from met file', fnamenc
          stop
      else if (ierr .ne. 0) then
          write(*,9125) 'unable to locate date/time in met file',  &
              'indj, itime =', indj2, itime, fnamenc
          stop
      else 
          jul = juldate( jyyyymmdd, jhhmmss )
          duma = (jul-bdate)*86400.
          iduma = nint(duma)
          if (iduma .ne. wftime(indj2)) goto 1100
      end if
      write(*,*) 
      write(*,*) 'readwind_nests processing wrfout file ='
      write(*,*) fnamenc
      write(*,*) 'itime, ymd, hms =', itime, jyyyymmdd, jhhmmss

      kbgn = 1 + add_sfc_level

      if (wind_option.eq.1) varname = 'AVGFLX_RUM'
      do i = 1, ndims_max
          lendim_exp(i) = 0
          lendim_max(i) = 1
      end do
      lendim_exp(1) = nxn(l)+1
      lendim_max(1) = nxmaxn
      lendim_exp(2) = nyn(l)
      lendim_max(2) = nymaxn
      lendim_exp(3) = nuvz-add_sfc_level
      lendim_max(3) = nuvzmax
      ndims_exp = 4
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, uuhn(0,0,kbgn,l), &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of U', fnamenc
      if (wind_option.le.0) print*,'you asked snapshot winds'
      if (wind_option.eq.1) print*,'you asked mean winds'
        print*,'change wind_option'

          stop
      end if


! v wind velocity
!   the wrf output file contains (nuvz-add_sfc_level) levels
!   read the data into k=kbgn,nuvz
!   (interpolate it from "V-grid" to "T-grid" later)
      if (wind_option.eq.1) varname = 'AVGFLX_RVM'
      lendim_exp(1) = nxn(l)
      lendim_max(1) = nxmaxn
      lendim_exp(2) = nyn(l)+1
      lendim_max(2) = nymaxn
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, vvhn(0,0,kbgn,l), &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of V', fnamenc
      if (wind_option.eq.0) print*,'you asked snapshot winds'
      if (wind_option.eq.1) print*,'you asked mean winds'
        print*,'change wind_option'
          stop
      end if


! w wind velocity
!   this is on the "W-grid", and 
!   the wrf output file contains nwz levels, so no shifting needed
      if (wind_option.eq.1) varname = 'AVGFLX_WWM'
!      print*,'varname',varname
      lendim_exp(1) = nxn(l)
      lendim_max(1) = nxmaxn
      lendim_exp(2) = nyn(l)
      lendim_max(2) = nymaxn
      lendim_exp(3) = nwz
      lendim_max(3) = nwzmax
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, wwhn(0,0,1,l), &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of W', fnamenc
      if (wind_option.eq.0) print*,'you asked snapshot winds'
      if (wind_option.eq.1) print*,'you asked mean winds'
        print*,'change wind_option'
          stop
      end if

! READ THE SECOND FILE and average
! deltat must be equal to deltat2, otherwise the time-average wind cannot be
! fixed.
! deltat is assumed to be the same than the WRF output.

       if (ldirect.eq.1) indj2=indj+1
       if (ldirect.eq.-1) indj2=indj-1
       deltat2=wfdt(indj2)
      fnamenc = path(m)(1:length(m)) // wfnamen(l,indj2)

     if (deltat.eq.deltat2 .and. indj2.ge.1 .and.indj2.le.numbwf) then
! locate the date/time in the file
      itime = 0
1101  itime = itime + 1
      call read_ncwrfout_1datetime( ierr, fnamenc, &
          itime, jyyyymmdd, jhhmmss )
      if (ierr .eq. -1) then
          write(*,9100) 'error reading time from met file', fnamenc
          stop
      else if (ierr .ne. 0) then
          write(*,9125) 'unable to locate date/time in met file',  &
              'indj, itime =', indj2, itime, fnamenc
          stop
      else 
          jul = juldate( jyyyymmdd, jhhmmss )
          duma = (jul-bdate)*86400.
          iduma = nint(duma)
          if (iduma .ne. wftime(indj2)) goto 1101
      end if
      if (option_verbose.eq.1) then
      write(*,*) 
      write(*,*) 'readwind processing wrfout file ='
      write(*,*) fnamenc
      write(*,*) 'itime, ymd, hms =', itime, jyyyymmdd, jhhmmss
      endif
      kbgn = 1 + add_sfc_level

      if (wind_option.eq.1) varname = 'AVGFLX_RUM'
      do i = 1, ndims_max
          lendim_exp(i) = 0
          lendim_max(i) = 1
      end do
      lendim_exp(1) = nxn(l)+1
      lendim_max(1) = nxmaxn
      lendim_exp(2) = nyn(l)
      lendim_max(2) = nymaxn
      lendim_exp(3) = nuvz-add_sfc_level
      lendim_max(3) = nuvzmax
      ndims_exp = 4
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, uuhn2(0,0,kbgn,l), &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of U', fnamenc
      if (wind_option.le.0) print*,'you asked snapshot winds'
      if (wind_option.eq.1) print*,'you asked mean winds'
        print*,'change wind_option'

          stop
      end if


! v wind velocity
!   the wrf output file contains (nuvz-add_sfc_level) levels
!   read the data into k=kbgn,nuvz
!   (interpolate it from "V-grid" to "T-grid" later)
      if (wind_option.eq.1) varname = 'AVGFLX_RVM'
      lendim_exp(1) = nxn(l)
      lendim_max(1) = nxmaxn
      lendim_exp(2) = nyn(l)+1
      lendim_max(2) = nymaxn
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, vvhn2(0,0,kbgn,l), &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of V', fnamenc
      if (wind_option.eq.0) print*,'you asked snapshot winds'
      if (wind_option.eq.1) print*,'you asked mean winds'
        print*,'change wind_option'
          stop
      end if


! w wind velocity
!   this is on the "W-grid", and 
!   the wrf output file contains nwz levels, so no shifting needed
      if (wind_option.eq.1) varname = 'AVGFLX_WWM'
!     print*,'varname',varname
      lendim_exp(1) = nxn(l)
      lendim_max(1) = nxmaxn
      lendim_exp(2) = nyn(l)
      lendim_max(2) = nymaxn
      lendim_exp(3) = nwz
      lendim_max(3) = nwzmax
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, wwhn2(0,0,1,l), &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of W', fnamenc
      if (wind_option.eq.0) print*,'you asked snapshot winds'
      if (wind_option.eq.1) print*,'you asked mean winds'
        print*,'change wind_option'
          stop
      end if


         do k=1,nuvzmax
        do j=0,nymaxn-1
       do i=0,nxmaxn-1
           uuhn(i,j,k,l)=0.5*(uuhn(i,j,k,l)+uuhn2(i,j,k,l))
           vvhn(i,j,k,l)=0.5*(vvhn(i,j,k,l)+vvhn2(i,j,k,l))
         enddo
        enddo
       enddo
         do k=1,nwzmax
        do j=0,nymaxn-1
       do i=0,nxmaxn-1
           wwhn(i,j,k,l)=0.5*(wwhn(i,j,k,l)+wwhn2(i,j,k,l))
         enddo
        enddo
       enddo


      endif ! test on deltat

      enddo ! loop over the nests

      return    
      end subroutine readwind_nests_timeav

