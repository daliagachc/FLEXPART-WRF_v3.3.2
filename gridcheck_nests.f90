      subroutine gridcheck_nests
!*******************************************************************************
!                                                                              *
!     This routine checks the grid specification for the nested model domains. *
!     It is similar to subroutine gridcheck, which checks the mother domain.   *
!                                                                              *
!     Note:  This is the FLEXPART_WRF version of subroutine gridcheck.         *
!            The computational grid is the WRF x-y grid rather than lat-lon.   *
!            There are many differences from the FLEXPART version.             *
!                                                                              *
!     Authors: A. Stohl, G. Wotawa                                             *
!     8 February 1999                                                          *
!                                                                              *
!  Nov 2005, R. Easter:                                                        *
!             MAJOR revisions for FLEXPART_WRF                                 *
!  Dec 2005, R. Easter:                                                        *
!             changed names of "*lon0*" & "*lat0*" variables                   *
!                                                                              *
!  2015-03-26, A. Dingwell:                                                    *
!             Updated calls to read_ncwrfout_gridinfo to match updates in that *
!             subroutine.                                                      * 
!                                                                              *
!*******************************************************************************

!  use grib_api
  use par_mod
  use com_mod
  use flexwrf_ncdf_mod



      integer,parameter :: ndims_max=4
      integer :: i, ierr, ifn, itime, ix
      integer :: idiagaa, idiagaa_1, idiagaa_2, idiagbb
      integer :: iduma, idumb
      integer :: jy
      integer :: k
      integer :: l, lp
      integer :: lendim(ndims_max), lendim_exp(ndims_max), & 
          lendim_max(ndims_max)
      integer :: m
      integer :: map_proj_id_dum
      integer :: ndims, ndims_exp, &
              ext_scalar,pbl_physics,mp_physics_dum, num_land_cat_dum
      integer :: n_west_east, n_south_north, n_bottom_top
      integer :: nuvzn, nwzn

      real :: dx_met, dy_met
      real :: duma, dumb, dumc, dumx, dumy
      real :: dump1, dump2, dumdz
      real :: eta_w_wrf_nest(nwzmax), eta_u_wrf_nest(nwzmax) 
      real :: map_stdlon_dum, map_truelat1_dum, map_truelat2_dum
      real :: pint, p_top_wrf_nest
      real :: xaux1, xaux2, yaux1, yaux2

      character(len=160) :: fnamenc, varname

      real(kind=4) :: m_un(0:nxmaxn-1,0:nymaxn-1,1,maxnests)
      real(kind=4) :: m_vn(0:nxmaxn-1,0:nymaxn-1,1,maxnests)

      real :: tmp_vardata(nzmax)


      ! set default value of z-dimension to 1 (indicating 2d field)
      ! this prevents an error 
      ! At line 935 of file read_ncwrfout.f90
      ! Fortran runtime error: Index '1' of dimension 3 of array 'vardata_out' above upper bound of 0
      lendim_max(3) = 1

! Loop about all nesting levels
!******************************
!     idiagaa_1 = 1
      idiagaa_1 = 0
      idiagaa_2 = 0
!     idiagbb = 10
      idiagbb = 0

      do l=1,numbnests

      write(*,'(//a,i3)') 'gridcheck_nests output for grid #', l

!
!   get grid info from the wrf netcdf file
!
      if(ideltas.gt.0) then
        ifn=1
      else
        ifn=numbwf
      endif
      m = numpath+2*(l-1)+1
      fnamenc = path(m)(1:length(m)) // wfnamen(l,ifn)

      idiagaa = idiagaa_1

      call read_ncwrfout_gridinfo( ierr, idiagaa, fnamenc, &
        n_west_east, n_south_north, n_bottom_top,  &
        dx_met, dy_met,  &
        m_grid_id(l), m_parent_grid_id(l), m_parent_grid_ratio(l),  &
        i_parent_start(l), j_parent_start(l), &
        map_proj_id_dum, map_stdlon_dum,  &
        map_truelat1_dum, map_truelat2_dum, &
        ext_scalar,pbl_physics,mp_physics_dum, num_land_cat_dum)
      if (ierr .ne. 0) goto 999


      mp_physicsn(l)=mp_physics_dum

! subtract 1 because i & j indexing in flexpart always starts at 0
      i_parent_start(l) = i_parent_start(l)-1
      j_parent_start(l) = j_parent_start(l)-1

!
! set grid dimension and size variables
!
      nxn(l) = n_west_east
      nyn(l) = n_south_north

      nuvzn = n_bottom_top
      nwzn = n_bottom_top + 1

! for FLEXPART_WRF, x & y coords are in meters
      dxn(l) = dx_met
      dyn(l) = dy_met


!
! check that grid dimensions are not too big
!
! flexpart_wrf 07-nov-2005 - require (nxn+1 .le. nxmaxn) and (nyn+1 .le. nymaxn)
! because u,v in met. files are on staggered grid
      if (nxn(l)+1 .gt. nxmaxn) then                         
        write(*,*) 'FLEXPART gridcheck_nests error: ' // &
                   'Too many grid points in x direction.'
        write(*,*) 'Change parameter settings in file par_mod.'
        write(*,*) 'l, nxn(l)+1, nxmaxn =', l, nxn(l), nxmaxn
        stop
      endif

      if (nyn(l)+1 .gt. nymaxn) then                         
        write(*,*) 'FLEXPART gridcheck_nests error: ' // &
                   'Too many grid points in y direction.'
        write(*,*) 'Change parameter settings in file par_mod.'
        write(*,*) 'l, nyn(l)+1, nymaxn =', l, nyn(l), nymaxn
        stop
      endif

      nuvzn = nuvzn+add_sfc_level
      if (nuvzn .ne. nuvz) then                         
        write(*,*) 'FLEXPART gridcheck_nests error: ' // &
                   'nuvzn and nuvz differ'
        write(*,*) 'l, nuvzn, nuvz =', l, nuvzn, nuvz
        stop
      endif

      if (nwzn .ne. nwz) then                         
        write(*,*) 'FLEXPART gridcheck_nests error: ' // &
                   'nwzn and nwz differ'
        write(*,*) 'l, nwzn, nwz =', l, nwzn, nwz
        stop
      endif

! check that map projection info matches parent
      duma = 3.0e-7*max( abs(map_stdlon),   1.0e-30 )
      dumb = 3.0e-7*max( abs(map_truelat1), 1.0e-30 )
      dumc = 3.0e-7*max( abs(map_truelat2), 1.0e-30 )
      iduma = 0
      if (map_proj_id .ne. map_proj_id_dum)             iduma = 1
      if (abs(map_stdlon-map_stdlon_dum)     .gt. duma) iduma = 2
      if (abs(map_truelat1-map_truelat1_dum) .gt. dumb) iduma = 3
      if (abs(map_truelat2-map_truelat2_dum) .gt. dumc) iduma = 4
      if (iduma .ne. 0) then
        write(*,*) 'FLEXPART gridcheck_nests error: ' // &
                   'map projection parameters differ'
        write(*,*) 'l, map param #=', l, iduma
        stop
      end if

      varname = 'MAPFAC_MX'
      lendim_exp(1) = nxn(l)
      lendim_max(1) = nxmaxn
      lendim_exp(2) = nyn(l)
      lendim_max(2) = nymaxn
      ndims_exp = 3
      itime=1
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
          varname, m_xn(0,0,1,l), &
          itime, &
          ndims, ndims_exp, ndims_max, &
          lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
      varname = 'MAPFAC_M'
      lendim_exp(1) = nxn(l)
      lendim_max(1) = nxmaxn
      lendim_exp(2) = nyn(l)
      lendim_max(2) = nymaxn
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
          varname, m_xn(0,0,1,l), &
          itime, &
          ndims, ndims_exp, ndims_max, &
          lendim, lendim_exp, lendim_max )
      endif
      if (ierr .ne. 0) then
          print*,'error doing MAP X'
      varname = 'MAPFAC_UX'
      lendim_exp(1) = nxn(l)+1
      lendim_max(1) = nxmaxn
      lendim_exp(2) = nyn(l)
      lendim_max(2) = nymaxn
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
          varname, m_un(0,0,1,l), &
          itime, &
          ndims, ndims_exp, ndims_max, &
          lendim, lendim_exp, lendim_max )
      do j = 0, nyn(l)-1
      do i = 0, nxn(l)-1
      m_xn(i,j,1,l)=(m_un(i,j,1,l)+m_un(i+1,j,1,l))*0.5
      enddo
      enddo
      if (ierr .ne. 0) then
          print*,'error doing MAP U'
          print*,'NO MAP FACTOR IS GOING TO BE USED.'
          print*,'LARGE UNCERTAINTIES TO BE EXPECTED'
      do j = 0, nyn(l)-1
      do i = 0, nxn(l)-1
      m_xn(i,j,1,l)=1.
      enddo
      enddo
      end if
      end if

      varname = 'MAPFAC_MY'
      lendim_exp(1) = nxn(l)
      lendim_max(1) = nxmaxn
      lendim_exp(2) = nyn(l)
      lendim_max(2) = nymaxn

      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
          varname, m_yn(0,0,1,l), &
          itime, &
          ndims, ndims_exp, ndims_max, &
          lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
      varname = 'MAPFAC_M'
      lendim_exp(1) = nxn(l)
      lendim_max(1) = nxmaxn
      lendim_exp(2) = nyn(l)
      lendim_max(2) = nymaxn
     call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
          varname, m_yn(0,0,1,l), &
          itime, &
          ndims, ndims_exp, ndims_max, &
          lendim, lendim_exp, lendim_max )
      endif
      if (ierr .ne. 0) then
          print*,'error doing MAP Y'
      varname = 'MAPFAC_VY'
      lendim_exp(1) = nxn(l)
      lendim_max(1) = nxmaxn
      lendim_exp(2) = nyn(l)+1
      lendim_max(2) = nymaxn
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
          varname, m_vn(0,0,1,l), &
          itime, &
          ndims, ndims_exp, ndims_max, &
          lendim, lendim_exp, lendim_max )
      do j = 0, nyn(l)-1
      do i = 0, nxn(l)-1
      m_yn(i,j,1,l)=(m_vn(i,j,1,l)+m_vn(i,j+1,1,l))*0.5
      enddo
      enddo
      if (ierr .ne. 0) then
          print*,'ERROR doing MAP V'
          print*,'NO MAP FACTOR IS GOING TO BE USED.'
          print*,'LARGE UNCERTAINTIES TO BE EXPECTED'
      do j = 0, nyn(l)-1
      do i = 0, nxn(l)-1
      m_yn(i,j,1,l)=1.
      enddo
      enddo
      end if
      end if
      lendim_exp(1) = nxn(l)
      lendim_max(1) = nxmaxn
      lendim_exp(2) = nyn(l)
      lendim_max(2) = nymaxn

!
!   read latitude and longitude
!   read oro, lsm, and excessoro
!

      idiagaa = idiagaa_2

      varname = 'XLAT'
      do i = 1, ndims_max
          lendim_exp(i) = 0
          lendim_max(i) = 1
      end do
      itime = 1
      lendim_exp(1) = nxn(l)
      lendim_max(1) = nxmaxn
      lendim_exp(2) = nyn(l)
      lendim_max(2) = nymaxn
      ndims_exp = 3
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, ylat2dn(0,0,l), &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max ) 
      if (ierr .ne. 0) then
          write(*,*)
          write(*,*) '*** checkgrid -- error doing ncread of XLAT'
          stop
      end if

      varname = 'XLONG'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, xlon2dn(0,0,l), &
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
      	  varname, oron(0,0,l), &
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
      do jy=0,nyn(l)-1
      do ix=0,nxn(l)-1
          lsmn(ix,jy,l)=1.0
      end do
      end do

! for now, set excessoro=0.0
      do jy=0,nyn(l)-1
      do ix=0,nxn(l)-1
          excessoron(ix,jy,l)=0.0
      end do
      end do
      do jy=1,nyn(l)-2
      do ix=1,nxn(l)-2
      m=oron(ix,jy,l)+oron(ix-1,jy,l)+oron(ix+1,jy,l)+ &
        oron(ix,jy-1,l)+oron(ix,jy+1,l)
      m=m/5.
      excessoron(ix,jy,l)=(oron(ix,jy,l)-m)**2.+(oron(ix-1,jy,l)-m)**2. &
      +(oron(ix+1,jy,l)-m)**2.+(oron(ix,jy-1,l)-m)**2.+(oron(ix,jy+1,l)-m)**2.
      excessoron(ix,jy,l)=(excessoron(ix,jy,l)/5.)**0.5
      end do
      end do


!
! identify the parent grid (which is probably "l-1", so this code 
!    may be more complicated that necessary)
! set xmet0n, ymet0n, which  are the x,y coords of lower-left corner 
!   of a nested grid (in meters on the mother grid)
!
! note on dumc:
!   the lower-left corner of the nested cell (0,0) coincides with the 
!       lower-left corner of the parent cell (i_parent_start,j_parent_start)
!   this being the case, the center of nested cell (0,0) is shifted
!       by (-dumc*parent_gridsize) relative to the center of the parent cell
!   (for m_parent_grid_ratio = 2, 3, 4, 5; dumc = 1/4, 1/3, 3/8, 2/5)
!
      l_parent_nest_id(l) = -1
      if (m_parent_grid_id(l) .gt. 0) then
      do lp = 0, l-1
         if ( (l_parent_nest_id(l) .eq. -1) .and. &
              (m_parent_grid_id(l) .eq. m_grid_id(lp)) ) then

            l_parent_nest_id(l) = lp
            m = m_parent_grid_ratio(l)
            dumc = real(m-1)/real(m*2)
            if (lp .eq. 0) then
               xmet0n(l) = xmet0 + dx*(i_parent_start(l)-dumc)
               ymet0n(l) = ymet0 + dy*(j_parent_start(l)-dumc)
            else
               xmet0n(l) = xmet0n(lp) + dxn(lp)*(i_parent_start(l)-dumc)
               ymet0n(l) = ymet0n(lp) + dyn(lp)*(j_parent_start(l)-dumc)
            end if
         end if
      end do
      end if
      if (idiagbb .gt. 0) write(*,'(/a,3i8)')  &
            'l, m_grid_id(l), m_parent_grid_id(l)', &
             l, m_grid_id(l), m_parent_grid_id(l)

      if (l_parent_nest_id(l) .eq. -1) then
         write(*,'(/a,i3/)')  &
            'gridcheck_nests fatal error -- ' // &
            'parent grid not found for l =', l
         stop
      end if

!
! diagnostics for testing the nesting calculations
! (set idiagbb=0 to turn it off)
!
      lp = l_parent_nest_id(l)
      if (idiagbb .gt. 0) then
         write(*,'(a,2i8)') 'l_parent, m_grid_id(l_parent)       ', &
                             lp,       m_grid_id(lp)
         write(*,'(a,2i8)') 'm_parent_grid_ratio(l)              ', &
                             m_parent_grid_ratio(l)
         write(*,'(a,i8,f11.2)') &
                            'i_parent_start(l), xi_...           ', &
                             i_parent_start(l), i_parent_start(l)-dumc
         write(*,'(a,i8,f11.2)') &
                            'j_parent_start(l), yj_...           ', &
                             j_parent_start(l), j_parent_start(l)-dumc
      end if
!23456789012345678901234567890123456789012345678901234567890123456789012

      if (idiagbb .gt. 0) then
         write(*,*)
         do jy = j_parent_start(l)-1, j_parent_start(l)
         do ix = i_parent_start(l)-1, i_parent_start(l)
            if (lp .eq. 0) then
               write(*,'(a,2i7,2f11.4)') 'parent i,j,lon,lat', &
                  ix, jy, xlon2d(ix,jy), ylat2d(ix,jy)
            else
               write(*,'(a,2i7,2f11.4)') 'parent i,j,lon,lat', &
                  ix, jy, xlon2dn(ix,jy,lp), ylat2dn(ix,jy,lp) 
            end if
         end do
         end do

         dumc = real( (m_parent_grid_ratio(l)-1) )/ &
                real( (m_parent_grid_ratio(l)*2) )
         dumx = i_parent_start(l) - dumc
         dumy = j_parent_start(l) - dumc
         call xyindex_to_ll_wrf( lp, dumx, dumy, xaux1, yaux1 )
         write(*,'(a,2f7.2,2f11.4)') 'par. xi,yj,lon,lat', &
            dumx, dumy, xaux1, yaux1

         write(*,*)
         write(*,'(a,2i7,2f11.4)') 'nest   i,j,lon,lat', &
            0, 0, xlon2dn(0,0,l), ylat2dn(0,0,l)
         write(*,*)

         dumx = (xmet0n(l) - xmet0)/dx
         dumy = (ymet0n(l) - ymet0)/dy
         call xyindex_to_ll_wrf( 0, dumx, dumy, xaux1, yaux1 )
         write(*,'(a,2f7.2,2f11.4)') 'mot. xi,yj,lon,lat', &
            dumx, dumy, xaux1, yaux1

         iduma = max( 0, ifix(dumx) )
         idumb = max( 0, ifix(dumy) )
         do jy = idumb, idumb+1
         do ix = iduma, iduma+1
               write(*,'(a,2i7,2f11.4)') 'mother i,j,lon,lat', &
                  ix, jy, xlon2d(ix,jy), ylat2d(ix,jy)
         end do
         end do
      end if



! Output of grid info
!********************

      write(*,'(/a,i2)')  &
          'gridcheck_nests -- nested domain #: ',l
      write(*,'(a,f12.1,a1,f12.1,a,f10.1)')  &
          '  X coordinate range: ', &
          xmet0n(l),' to ',xmet0n(l)+(nxn(l)-1)*dxn(l), &
          '   Grid distance: ',dxn(l)
      write(*,'(a,f12.1,a1,f12.1,a,f10.1)')  &
          '  Y coordinate range: ', &
          ymet0n(l),' to ',ymet0n(l)+(nyn(l)-1)*dyn(l), &
          '   Grid distance: ',dyn(l)
      write(*,*)


! Determine, how much the resolutions in the nests are enhanced as
! compared to the mother grid
!*****************************************************************

        xresoln(l)=dx/dxn(l)
        yresoln(l)=dy/dyn(l)


! Determine the mother grid coordinates of the corner points of the
! nested grids
! Convert first to geographical coordinates, then to grid coordinates
!********************************************************************

        xaux1=xmet0n(l)
        xaux2=xmet0n(l)+real(nxn(l)-1)*dxn(l)
        yaux1=ymet0n(l)
        yaux2=ymet0n(l)+real(nyn(l)-1)*dyn(l)

        xln(l)=(xaux1-xmet0)/dx
        xrn(l)=(xaux2-xmet0)/dx
        yln(l)=(yaux1-ymet0)/dy
        yrn(l)=(yaux2-ymet0)/dy


        if ((xln(l).lt.0.).or.(yln(l).lt.0.).or. &
        (xrn(l).gt.real(nxmin1)).or.(yrn(l).gt.real(nymin1))) then
          write(*,*) 'gridcheck_nests error'
          write(*,*) 'Nested domain does not fit into mother domain'
          write(*,*) 'Execution is terminated.'
          stop
        endif


! check that the map projection routines are working
      call test_xyindex_to_ll_wrf( l )


! Check, whether the heights of the model levels of the nested
! wind fields are consistent with those of the mother domain.
! If not, terminate model run.
!*************************************************************

! first read eta_w_wrf_nest, eta_u_wrf_nest, and p_top_wrf_nest 
! from the netcdf wrfout file
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
      	  varname, eta_w_wrf_nest, &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,*)
          write(*,*) '*** checkgrid -- error doing ncread of ZNW'
          stop
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
      	  varname, eta_u_wrf_nest, &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,*)
          write(*,*) '*** checkgrid -- error doing ncread of ZNU'
          stop
      end if

      varname = 'P_TOP'
      do i = 1, ndims_max
          lendim_exp(i) = 0
          lendim_max(i) = 1
      end do
      lendim_exp(1) = 1
      lendim_max(1) = 6
      ndims_exp = 2
      if (ext_scalar .lt. 0) ndims_exp=1
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc, &
      	  varname, tmp_vardata, &
      	  itime, &
      	  ndims, ndims_exp, ndims_max, &
      	  lendim, lendim_exp, lendim_max )
      p_top_wrf_nest = tmp_vardata(1)
      if (ierr .ne. 0) then
          write(*,*)
          write(*,*) '*** checkgrid -- error doing ncread of P_TOP'
          stop
      end if


      do k = 1, nwz
          duma = 3.0e-7*max( abs(eta_w_wrf(k)), 1.0e-30 )
          if (abs(eta_w_wrf(k)-eta_w_wrf_nest(k)) .gt. duma) then
              write(*,*)  &
              'FLEXPART gridcheck_nests error for nesting level',l
              write(*,*)  &
              'eta_w_wrf are not consistent with the mother domain'
              write(*,*) 'k, eta_w_wrf(k), eta_w_wrf_nest(k) =', &
                          k, eta_w_wrf(k), eta_w_wrf_nest(k)
              stop
          endif
      end do

      do k = 1, nwz-1
          duma = 3.0e-7*max( abs(eta_u_wrf(k)), 1.0e-30 )
          if (abs(eta_u_wrf(k)-eta_u_wrf_nest(k)) .gt. duma) then
              write(*,*)  &
              'FLEXPART gridcheck_nests error for nesting level',l
              write(*,*)  &
              'eta_u_wrf are not consistent with the mother domain'
              write(*,*) 'k, eta_u_wrf(k), eta_u_wrf_nest(k) =', &
                          k, eta_u_wrf(k), eta_u_wrf_nest(k)
              stop
          endif
      end do

      duma = 3.0e-7*max( abs(p_top_wrf), 1.0e-30 )
      if (abs(p_top_wrf-p_top_wrf_nest) .gt. duma) then
          write(*,*)  &
          'FLEXPART gridcheck_nests error for nesting level',l
          write(*,*)  &
          'p_top_wrf are not consistent with the mother domain'
          write(*,*) 'p_top_wrf, p_top_wrf_nest', &
                      p_top_wrf, p_top_wrf_nest
          stop
      endif


!   done with nest l
      enddo

      return


999   write(*,*)  
      write(*,*) ' ###########################################'// &
                 '###### '
      write(*,*) ' FLEXPART_WRF subroutine gridcheck_nests:  ' // &
                 'nesting level ', l
      write(*,*) ' can not open input data file '
      write(*,*) '     '//fnamenc
      write(*,*) ' or, an error occured in subr. read_ncwrfout_gridinfo'
      write(*,*) '     with ierr =', ierr
      write(*,*) ' ###########################################'// &
                 '###### '
      stop

end subroutine gridcheck_nests


