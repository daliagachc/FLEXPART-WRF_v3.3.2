      subroutine gridcheck
!**********************************************************************
!                                                                     * 
!             TRAJECTORY MODEL SUBROUTINE GRIDCHECK                   *
!             FLEXPART VERSION -> DO NOT USE IN FLEXTRA, PRAEPRO      *
!                                                                     *
!**********************************************************************
!                                                                     * 
! AUTHOR:      R. Easter & J. Fast, PNNL                              *
! DATE:        2005-autumn-??                                         *
!                                                                     * 
!  Dec 2005, R. Easter:                                               *
!             Changed names of "*lon0*" & "*lat0*" variables*         *
!                                                                     *
!  2015-03-26, A. Dingwell:                                           *
!             Updated calls to read_ncwrfout_gridinfo to match updates*
!             in that subroutine.                                     * 
!                                                                     * 
!**********************************************************************
!                                                                     *
! DESCRIPTION:                                                        *
!                                                                     *
! Note:  This is the FLEXPART_WRF version of subroutine gridcheck.    *
!    The computational grid is the WRF x-y grid rather than lat-lon.  *
!    There are many differences from the FLEXPART version.            *
!                                                                     *
! This subroutine determines the grid specifications                  *
! of the WRF model run from the first met. input file.                *
!     (longitudes & latitudes, number of grid points, grid distance,  *
!     vertical discretization)                                        *
! The consistancy of met. input files is checked in the routine       *
! "readwind" each call.  (Grid info must not change.)                 *
!                                                                     *
!                                                                     *
! XMET0,YMET0 -- x,y coordinates (in m) of the lower left             *
!                "T-grid" point                                       *
!           NOTE: These replace XLON0,YLAT0 of FLEXPART (ECMWF),      *
!           which uses longitude,longitude (degrees) as coordinates.  *
!                                                                     *
! NX        number of grid points x-direction                         *
! NY        number of grid points y-direction                         *
! DX        grid distance x-direction (in m)                          *
! DY        grid distance y-direction (in m)                          *
! NUVZ      number of grid points for horizontal wind                 *
!           components in z direction                                 *
! NWZ       number of grid points for vertical wind                   *
!           component in z direction                                  *
!                                                                     *
! For WRF, the "W" grid has NWZ == "bottom_top_stag" levels           *
! For WRF, the "U", "V", and "T" grids have                           *
!     NWZ-1 == "bottom_top"  levels.                                  *
! In the ecmwf FLEXPART, the "U", "V", and "T" grids are              *
!     "augmented" with an additional "surface" layer,                 *
!     and thus have NUVZ==NWZ levels                                  *
! Because of the high vertical resolution often used in WRF,          *
!     it may be desirable to eliminate this "surface layer".          *
!                                                                     *
!                                                                     *
! UVHEIGHT(1)-         heights of gridpoints where u and v are        *
! UVHEIGHT(NUVZ)       given                                          *
! WHEIGHT(1)-          heights of gridpoints where w is given         *
! WHEIGHT(NWZ)                                                        *
!                                                                     *
!**********************************************************************
!
!      include 'includepar'
!      include 'includecom'
!      include 'includeconv'
!  use grib_api
  use par_mod
  use com_mod
  use conv_mod
  use cmapf_mod, only: stlmbr,stcm2p
  use flexwrf_ncdf_mod

  integer,parameter :: ndims_max=4
  integer :: i, idiagaa, ierr, ifn, itime, ix
  integer :: jy
  integer :: kz
  integer :: lendim(ndims_max), lendim_exp(ndims_max), &
      lendim_max(ndims_max)
  integer :: ndims, ndims_exp, &
    ext_scalar, pbl_physics, mp_physics_dum, num_land_cat_dum
  integer :: n_west_east, n_south_north, n_bottom_top

  real :: dx_met, dy_met
  real :: duma, dumb, dumc
  real :: dump1, dump2, dumdz
  real :: pint,m

  character(len=160) :: fnamenc, varname,fnamenc2

  real(kind=4) :: m_u(0:nxmax-1,0:nymax-1,1)
  real(kind=4) :: m_v(0:nxmax-1,0:nymax-1,1)

  real :: tmp_vardata(nzmax)


! set default value of z-dimension to 1 (indicating 2d field)
!  -- perhaps it would improve readability to do this each time
!  -- the other elements of lendim_max are set below
lendim_max(3) = 1

!
!   get grid info from the wrf netcdf file
!
!     write(*,'(//a)') 'gridcheck output'

      if(ideltas.gt.0) then
        ifn=1
      else
        ifn=numbwf
      endif
      fnamenc = path(2)(1:length(2))//wfname(ifn)
      idiagaa = 0

      call read_ncwrfout_gridinfo( ierr, idiagaa, fnamenc, &
        n_west_east, n_south_north, n_bottom_top,  &
        dx_met, dy_met, &
        m_grid_id(0), m_parent_grid_id(0), m_parent_grid_ratio(0), &
        i_parent_start(0), j_parent_start(0), &
        map_proj_id, map_stdlon, map_truelat1, map_truelat2, &
        ext_scalar,pbl_physics,mp_physics_dum, num_land_cat_dum)

       mp_physics=mp_physics_dum

      if (ierr.ne.0) goto 999

!
! set grid dimension and size variables
!
      nx = n_west_east
      nxfield = nx
      ny = n_south_north

      nxmin1=nx-1
      nymin1=ny-1
!    print*,'nxo',nxmin1,nymin1
      nuvz = n_bottom_top
      nwz = n_bottom_top + 1
      nlev_ec = n_bottom_top

! for FLEXPART_WRF, x & y coords are in meters, and
! we define lower-left corner of outermost (mother) grid == (0.0,0.0)
      xmet0 = 0.0
      ymet0 = 0.0

      dx = dx_met
      dy = dy_met

      dxconst = 1.0/dx
      dyconst = 1.0/dy

      l_parent_nest_id(0) = -1

      xglobal=.false.
      if (nxshift.ne.0) stop  &
        'nxshift (par_mod) must be zero for non-global domain'

      sglobal=.false.
      switchsouthg=999999.

      nglobal=.false.
      switchnorthg=999999.

      if (nxshift.lt.0) stop 'nxshift (par_mod) must not be negative'
      if (nxshift.ge.nxfield) stop 'nxshift (par_mod) too large'
     
!
! check that grid dimensions are not too big
!
! FLEXPART_WRF 07-nov-2005 - require (nx+1 .le. nxmax) and (ny+1 .le. nymax)
! because u,v in met. files are on staggered grid
!     if (nx.gt.nxmax) then                         
      if (nx+1 .gt. nxmax) then                         
       write(*,*) 'FLEXPART error: Too many grid points in x direction.'
        write(*,*) 'Reduce resolution of wind fields.'
        write(*,*) 'Or change parameter settings in file par_mod.'
        write(*,*) 'nx+1, nxmax =', nx,nxmax
        stop
      endif

!     if (ny.gt.nymax) then                         
      if (ny+1 .gt. nymax) then                         
       write(*,*) 'FLEXPART error: Too many grid points in y direction.'
        write(*,*) 'Reduce resolution of wind fields.'
        write(*,*) 'Or change parameter settings in file par_mod.'
        write(*,*) 'ny+1, nymax =', ny,nymax
        stop
      endif

      if (nuvz+1.gt.nuvzmax) then                         
        write(*,*) 'FLEXPART error: Too many u,v grid points in z '// &
        'direction.'
        write(*,*) 'Reduce resolution of wind fields.'
        write(*,*) 'Or change parameter settings in file par_mod.'
        write(*,*) nuvz+1,nuvzmax
        stop
      endif

      if (nwz.gt.nwzmax) then                         
        write(*,*) 'FLEXPART error: Too many w grid points in z '// &
        'direction.'
        write(*,*) 'Reduce resolution of wind fields.'
        write(*,*) 'Or change parameter settings in file par_mod.'
        write(*,*) nwz,nwzmax
        stop
      endif

!
!   read latitude and longitude
!   read oro, lsm, and excessoro

      varname = 'MAPFAC_MX'
      lendim_exp(1) = nx
      lendim_max(1) = nxmax
      lendim_exp(2) = ny
      lendim_max(2) = nymax
      ndims_exp = 3
      itime=1
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
          varname, m_x(0,0,1), &
          itime, &
          ndims, ndims_exp, ndims_max, &
          lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
      varname = 'MAPFAC_M'
      lendim_exp(1) = nx
      lendim_max(1) = nxmax
      lendim_exp(2) = ny
      lendim_max(2) = nymax

      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
          varname, m_x(0,0,1), &
          itime, &
          ndims, ndims_exp, ndims_max, &
          lendim, lendim_exp, lendim_max )
      endif 
      if (ierr .ne. 0) then
          print*,'error doing MAP X'
      varname = 'MAPFAC_UX'
      lendim_exp(1) = nx+1
      lendim_max(1) = nxmax
      lendim_exp(2) = ny
      lendim_max(2) = nymax
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
          varname, m_u(0,0,1), &
          itime, &
          ndims, ndims_exp, ndims_max, &
          lendim, lendim_exp, lendim_max )
      do j = 0, nymin1
      do i = 0, nxmin1
      m_x(i,j,1)=(m_u(i,j,1)+m_u(i+1,j,1))*0.5
      enddo
      enddo
      if (ierr .ne. 0) then
          print*,'error doing MAP U'
          print*,'NO MAP FACTOR IS GOING TO BE USED.'
          print*,'LARGE UNCERTAINTIES TO BE EXPECTED'
      do j = 0, nymin1
      do i = 0, nxmin1
      m_x(i,j,1)=1.
      enddo
      enddo
      end if
      end if
!     do j = 0, nymin1
!     do i = 0, nxmin1
!     m_x(i,j,1)=1.
!     enddo
!     enddo

      varname = 'MAPFAC_MY'
      lendim_exp(1) = nx
      lendim_max(1) = nxmax
      lendim_exp(2) = ny
      lendim_max(2) = nymax

      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
          varname, m_y(0,0,1), &
          itime, &
          ndims, ndims_exp, ndims_max, &
          lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
      varname = 'MAPFAC_M'
      lendim_exp(1) = nx
      lendim_max(1) = nxmax
      lendim_exp(2) = ny
      lendim_max(2) = nymax

      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
          varname, m_y(0,0,1), &
          itime, &
          ndims, ndims_exp, ndims_max, &
          lendim, lendim_exp, lendim_max )
      endif
      if (ierr .ne. 0) then
          print*,'error doing MAP Y'
      varname = 'MAPFAC_VY'
      lendim_exp(1) = nx
      lendim_max(1) = nxmax
      lendim_exp(2) = ny+1
      lendim_max(2) = nymax
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
          varname, m_v(0,0,1), &
          itime, &
          ndims, ndims_exp, ndims_max, &
          lendim, lendim_exp, lendim_max )
      do j = 0, nymin1
      do i = 0, nxmin1
      m_y(i,j,1)=(m_v(i,j,1)+m_v(i,j+1,1))*0.5
      enddo
      enddo
      if (ierr .ne. 0) then
          print*,'ERROR doing MAP V'
          print*,'NO MAP FACTOR IS GOING TO BE USED.'
          print*,'LARGE UNCERTAINTIES TO BE EXPECTED'
      do j = 0, nymin1
      do i = 0, nxmin1
      m_y(i,j,1)=1.
      enddo
      enddo
      end if
      end if
      lendim_exp(1) = nx
      lendim_max(1) = nxmax
      lendim_exp(2) = ny
      lendim_max(2) = nymax
!     do j = 0, nymin1
!     do i = 0, nxmin1
!     m_y(i,j,1)=1.
!     enddo
!     enddo

!

      varname = 'XLAT'
      do i = 1, ndims_max
          lendim_exp(i) = 0
          lendim_max(i) = 1
      end do
      itime = 1
      lendim_exp(1) = nx
      lendim_max(1) = nxmax
      lendim_exp(2) = ny
      lendim_max(2) = nymax
      ndims_exp = 3
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, ylat2d, &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,*)
          write(*,*) '*** checkgrid -- error doing ncread of XLAT'
          stop
      end if
!       print*,'values',ylat2d(309,101),ylat2d(300,101),ylat2d(250,101)
!       print*,'values',ylat2d(309,101),ylat2d(300,101),ylat2d(250,101)
      varname = 'XLONG'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, & 
      	  varname, xlon2d, &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,*)
          write(*,*) '*** checkgrid -- error doing ncread of XLONG'
          stop
      end if

      varname = 'HGT'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, oro, &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,*)
          write(*,*) '*** checkgrid -- error doing ncread of HGT'
          stop
      end if

! lsm = landsea mask == land fraction (or non-ocean fraction)
! for now, set lsm=1.0 which means land
      do jy=0,ny-1
      do ix=0,nxfield-1
!         lsm(ix,jy)=zsec4(nxfield*(ny-jy-1)+ix+1)
          lsm(ix,jy)=1.0
      end do
      end do

! for now, set excessoro=0.0
      do jy=0,ny-1
      do ix=0,nxfield-1
          excessoro(ix,jy)=0.0
      end do
      end do
      do jy=1,ny-2
      do ix=1,nxfield-2
      m=4.*oro(ix,jy)+oro(ix-1,jy)+oro(ix+1,jy)+oro(ix,jy-1)+oro(ix,jy+1)
      m=m/8.
      excessoro(ix,jy)=4.*(oro(ix,jy)-m)**2.+(oro(ix-1,jy)-m)**2. &
      +(oro(ix+1,jy)-m)**2.+(oro(ix,jy-1)-m)**2.+(oro(ix,jy+1)-m)**2.
      excessoro(ix,jy)=(excessoro(ix,jy)/8.)**0.5
      end do
      end do


! check that the map projection routines are working
      call test_xyindex_to_ll_wrf( 0 )


! If desired, shift all grids by nxshift grid cells
!**************************************************

!     if (xglobal) then
!       call shift_field_0(oro,nxfield,ny)
!       call shift_field_0(lsm,nxfield,ny)
!       call shift_field_0(excessoro,nxfield,ny)
!     endif


! CALCULATE VERTICAL DISCRETIZATION OF WRF MODEL
! PARAMETER akm,bkm DESCRIBE THE HYBRID "ETA" COORDINATE SYSTEM

! read eta_w_wrf, eta_u_wrf, and p_top_wrf from the netcdf wrfout file
      itime = 1

      varname = 'ZNW'
      do i = 1, ndims_max
          lendim_exp(i) = 0
          lendim_max(i) = 1
      end do
      lendim_exp(1) = nwz
      lendim_max(1) = nwzmax
      ndims_exp = 2
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, eta_w_wrf, &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
        fnamenc2='wrfout_d03_zn.nc'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc2, & 
      	  varname, eta_w_wrf, &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then

          write(*,*)
          write(*,*) '*** checkgrid -- error doing ncread of ZNW'
!         stop
      end if
      end if

      varname = 'ZNU'
      do i = 1, ndims_max
          lendim_exp(i) = 0
          lendim_max(i) = 1
      end do
      lendim_exp(1) = nwz-1
      lendim_max(1) = nwzmax
      ndims_exp = 2
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, eta_u_wrf, &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
        fnamenc2='wrfout_d03_zn.nc'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc2, &
      	  varname, eta_u_wrf, &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
        
      if (ierr .ne. 0) then
          write(*,*)
          write(*,*) '*** checkgrid -- error doing ncread of ZNU'
          stop
      end if
      end if

      varname = 'P_TOP'
      do i = 1, ndims_max
          lendim_exp(i) = 0
          lendim_max(i) = 1
      end do
      lendim_exp(1) = 1
      lendim_max(1) = 6
      ndims_exp = 2
      if (ext_scalar .lt. 0) ndims_exp = 1
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, tmp_vardata, &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max ) 
      p_top_wrf = tmp_vardata(1)
      if (ierr .ne. 0) then
          write(*,*)
          write(*,*) '*** checkgrid -- error doing ncread of P_TOP'
          stop
      end if

! diagnostics for testing
      if (idiagaa .gt. 0) then
          write(*,*)
          write(*,*) 'k, eta_w_wrf, eta_u_wrf ='
          write(*,'(i3,2f11.6)')  &
              (kz, eta_w_wrf(kz), eta_u_wrf(kz), kz=1,nwz-1)
          kz = nwz
          write(*,'(i3,2f11.6)') kz, eta_w_wrf(kz)
          write(*,*)
          write(*,*) 'p_top_wrf =', p_top_wrf
          write(*,*)

          duma = 0.0
          dumb = 1.0e30
          dumc = -1.0e30
          do jy = 0, ny-1
          do ix = 0, nx-1
              duma = duma + oro(ix,jy)
              dumb = min( dumb, oro(ix,jy) )
              dumc = max( dumc, oro(ix,jy) )
          end do
          end do
          duma = duma/(nx*ny)
          write(*,*) 'oro avg, min, max =', duma, dumb, dumc
          write(*,*)
      end if


!
! the wrf eta vertical grid at layer boundaries (w grid) and 
! layer centers (u grid) is defined by
!	eta_w_wrf(kz) = (pdh_w(kz) - p_top_wrf)/(pdh_surface - p_top_wrf)
!	eta_u_wrf(kz) = (pdh_u(kz) - p_top_wrf)/(pdh_surface - p_top_wrf)
! where "pdh_" refers to the dry hydrostatic component of the pressure
!
! so
!	pdh_w(kz) = ((1.0 - eta_w_wrf(kz))*p_top_wrf) + eta_w_wrf(kz)*pdh_surface
!	pdh_u(kz) = ((1.0 - eta_u_wrf(kz))*p_top_wrf) + eta_u_wrf(kz)*pdh_surface
!
! the ecmwf eta vertical grid is defined by
!	p_w(kz) = akm(kz) + bkm(kz)*p_surface
!	p_u(kz) = akz(kz) + bkz(kz)*p_surface
!
! the following definitions of akm, bkm, akz, bkz for wrf would be roughly 
! consistent those for ecmwf EXCEPT that for wrf, they involve the 
! dry hydrostatic component of the pressure
!     do kz = 1, nwz
!         akm(kz) = (1.0 - eta_w_wrf(kz))*p_top_wrf
!         bkm(kz) = eta_w_wrf(kz)
!     end do
!     do kz = 1, nuvz
!         akz(kz) = (1.0 - eta_u_wrf(kz))*p_top_wrf
!         bkz(kz) = eta_u_wrf(kz)
!     end do
!
! *** in FLEXPART_WRF we decided to used pressure from the met. files
!     and drop the akz/bkz/akm/bkm entirely ***
!


! in FLEXPART_ECMWF, the U, V, & T-grid levels are always shifted up by 1, 
!    and an extra near-surface level is defined at kz=1 
!    which is loaded with the 10 m winds and 2 m temperature & humidity
! for FLEXPART_WRF, this is optional, and is done when add_sfc_level=1
! (Note -- it will take a lot of effort to get rid of this augmented
!    level because many of the surface & boundary layer routines
!    are expecting it.  so for now, always augment.)

      dump1 = (101325.0-p_top_wrf)*eta_w_wrf(1) + p_top_wrf
      dump2 = (101325.0-p_top_wrf)*eta_w_wrf(2) + p_top_wrf
      dumdz = log(dump1/dump2)*8.4e3
      write(*,*)
!     write(*,*) 'add_sfc_level =', add_sfc_level
!     write(*,*) 'WRF layer 1 approx. thickness =', dumdz

      if (add_sfc_level .eq. 1) then
          nuvz = nuvz + 1
      else
          write(*,'(/a/a/)') '*** gridcheck fatal error ***', &
              '    add_sfc_level=0 is not yet implemented'
!         stop
      end if


!*******************************************************************************
! following comments are from FLEXPART_ECMWF.  This options for doubled vertical
! resolution has not been tried in FLEXPART_WRF, but it probably could be done
! with little effort.
! ------------------------------------------------------------------------------
! NOTE: In FLEXPART versions up to 4.0, the number of model levels was doubled
! upon the transformation to z levels. In order to save computer memory, this is
! not done anymore in the standard version. However, this option can still be
! switched on by replacing the following lines with those below, that are
! currently commented out. For this, similar changes are necessary in
! verttransform.f and verttranform_nests.f
!*******************************************************************************

      nz=nuvz
      if (nz.gt.nzmax) stop 'nzmax too small'

!     do 100 i=1,nuvz
!       aknew(i)=akz(i)
!100    bknew(i)=bkz(i)

! Switch on following lines to use doubled vertical resolution
!*************************************************************
!     nz=nuvz+nwz-1
!     if (nz.gt.nzmax) stop 'nzmax too small'
!     do 100 i=1,nwz
!       aknew(2*(i-1)+1)=akm(i)
!00     bknew(2*(i-1)+1)=bkm(i)
!     do 110 i=2,nuvz
!       aknew(2*(i-1))=akz(i)
!10     bknew(2*(i-1))=bkz(i)
! End doubled vertical resolution



! Determine the uppermost level for which the convection scheme shall be applied
! by assuming that there is no convection above 50 hPa (for standard SLP)
!
! FLEXPART_WRF - use approx. pressures to set nconvlev, and limit it to nuvz-2
!*******************************************************************************

      do i=1,nuvz-2
!     do i=1,nuvz-1
!       pint=akz(i)+bkz(i)*101325.
        pint = (101325.0-p_top_wrf)*eta_u_wrf(i) + p_top_wrf
        if (pint.lt.5000.0) goto 96
      enddo
96    nconvlev=i
      nconvlev = min( nconvlev, nuvz-2 )
      if (nconvlev.gt.nconvlevmax-1) then
        nconvlev=nconvlevmax-1
        pint = (101325.0-p_top_wrf)*eta_u_wrf(nconvlev) + p_top_wrf
        write(*,*) 'Attention, convection only calculated up to ', & 
            pint*0.01, ' hPa'
      endif


! Output of grid info
!********************

      write(*,*)
      write(*,*)
      write(*,'(a/a,2i7/a,2i7//a,3i7/a,2i7/a,4i7)')  &
        '# of vertical levels in WRF data',  &
        '    n_bottom_top & "true" nuvz:', n_bottom_top,  &
                                           (nuvz-add_sfc_level), &
        '    nwz &     "augmented" nuvz:', nwz, nuvz, &
        '    nwzmax, nuvzmax, nzmax    :', nwzmax, nuvzmax, nzmax, &
        '    nconvlevmax, nconvlev     :', nconvlevmax, nconvlev, &
        '    nx, ny, nxmax, nymax      :', nx, ny, nxmax, nymax  
      write(*,*)
      write(*,'(a)') 'Mother domain:'
      write(*,'(a,f10.1,a1,f10.1,a,f10.1)') '  east-west   range: ', &
        xmet0,' to ',xmet0+(nx-1)*dx,'   Grid distance: ',dx
      write(*,'(a,f10.1,a1,f10.1,a,f10.1)') '  south-north range: ', &
        ymet0,' to ',ymet0+(ny-1)*dy,'   Grid distance: ',dy
      write(*,*)


!
! all done
!
      return


! file open error
999   write(*,*)  
      write(*,*) ' ###########################################'// &
                 '###### '
      write(*,*) '       TRAJECTORY MODEL SUBROUTINE GRIDCHECK:'
      write(*,*) ' CAN NOT OPEN INPUT DATA FILE = '
      write(*,*) wfname(ifn)
      write(*,*) ' ###########################################'// & 
                 '###### '
      write(*,*)
      stop

end subroutine gridcheck



