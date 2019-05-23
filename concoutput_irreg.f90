!***********************************************************************
!* Copyright 2012,2013                                                 *
!* Jerome Brioude, Jerome Fast, 
!* Andreas Stohl, Petra Seibert, A. Frank, Gerhard Wotawa *
!* Caroline Forster, Sabine Eckhardt, John Burkhart, Harald Sodemann   *
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
      subroutine concoutput_irreg(itime,outnum,gridtotalunc,wetgridtotalunc, &
      drygridtotalunc)
!                             i     i          o             o
!            o
!*******************************************************************************
!                                                                              *
!     Note:  This is the FLEXPART_WRF version of subroutine concoutput         *
!                                                                              *
!     Output of the concentration grid and the receptor concentrations.        *
!                                                                              *
!     Author: A. Stohl                                                         *
!                                                                              *
!     24 May 1995                                                              *
!                                                                              *
!     13 April 1999, Major update: if output size is smaller, dump output      *
!                    in sparse matrix format; additional output of uncertainty *
!                                                                              *
!     05 April 2000, Major update: output of age classes; output for backward  *
!                    runs is time spent in grid cell times total mass of       *
!                    species.                                                  *
!                                                                              *
!     17 February 2002, Appropriate dimensions for backward and forward runs   *
!                       are now specified in file includepar                   *
!                                                                              *
!     Dec 2005, J. Fast - Output files can be either binary or ascii.          *
!                         Sparse output option is turned off.                  *
!     Dec 2005, R. Easter - changed names of "*lon0*" & "*lat0*" variables     *
!     2012, J. Brioude- modify output format to flexpart 8*, latlon output     *
!     2014-05-07  A. Dingwell: Suppressed output of receptor points to separate*
!                   files when netcdf output is set (bug reported by Barbosa)  *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
! outnum          number of samples                                            *
! ncells          number of cells with non-zero concentrations                 *
! sparse          .true. if in sparse matrix format, else .false.              *
! nspeciesdim     either nspec (forward runs), or numpoint (backward runs)     *
! tot_mu          1 for forward, initial mass mixing ration for backw. runs    *
! maxpointspec    maxspec for forward runs, maxpoint for backward runs         *
!                                                                              *
!*******************************************************************************

!      include 'includepar'
!      include 'includecom'
!
!      double precision jul
!      integer itime,i,ix,jy,kz,k,l,iix,jjy,kzz,nage,jjjjmmdd,ihmmss
!      integer ncells(maxpointspec,maxageclass)
!      integer ncellsd(maxpointspec,maxageclass)
!      integer ncellsw(maxpointspec,maxageclass),nspeciesdim
!      real outnum,weightair,densityoutrecept(maxreceptor),xl,yl
!      real densityoutgrid(0:maxxgrid-1,0:maxygrid-1,maxzgrid),
!     +grid(0:maxxgrid-1,0:maxygrid-1,maxzgrid,maxpointspec,maxageclass)
!      real wetgrid(0:maxxgrid-1,0:maxygrid-1,maxpointspec,maxageclass)
!      real drygrid(0:maxxgrid-1,0:maxygrid-1,maxpointspec,maxageclass)
!      real gridsigma(0:maxxgrid-1,0:maxygrid-1,maxzgrid,maxpointspec,
!     +maxageclass),
!     +drygridsigma(0:maxxgrid-1,0:maxygrid-1,maxpointspec,maxageclass),
!     +wetgridsigma(0:maxxgrid-1,0:maxygrid-1,maxpointspec,maxageclass)
!      real auxgrid(nclassunc),gridtotal,gridsigmatotal,gridtotalunc
!      real wetgridtotal,wetgridsigmatotal,wetgridtotalunc
!      real drygridtotal,drygridsigmatotal,drygridtotalunc
!      real factor(0:maxxgrid-1,0:maxygrid-1,maxzgrid)
!      real halfheight,dz,dz1,dz2,tot_mu(maxpointspec)
!      real xnelat,xnelon
!      real xsw,xne,ysw,yne,tmpx,tmpy,tmplon,tmplat
!      parameter(weightair=28.97)
!      logical sparse(maxpointspec,maxageclass)
!      logical sparsed(maxpointspec,maxageclass)
!      logical sparsew(maxpointspec,maxageclass)
!      character adate*8,atime*6
  use unc_mod
  use point_mod
  use outg_mod
  use par_mod
  use com_mod
  use netcdf_output_mod

  implicit none

  real(kind=dp) :: jul
  integer :: itime,i,ix,jy,kz,ks,kp,l,iix,jjy,kzz,nage,jjjjmmdd,ihmmss
  integer :: sp_count_i,sp_count_r
  real :: sp_fact
  real :: outnum,densityoutrecept(maxreceptor),xl,yl
  real :: auxgrid(nclassunc),gridtotal,gridsigmatotal,gridtotalunc
  real :: wetgridtotal,wetgridsigmatotal,wetgridtotalunc
  real :: drygridtotal,drygridsigmatotal,drygridtotalunc
  real :: halfheight,dz,dz1,dz2,tot_mu(maxspec,maxpointspec_act)
  real :: xsw,xne,ysw,yne,tmpx,tmpy,tmplon,tmplat
  real :: start, finish

  real,parameter :: smallnum = tiny(0.0) ! smallest number that can be handled
  ! real,parameter :: weightair=28.97 !AD: moved this to par_mod.f90
  logical :: sp_zer
  character :: adate*8,atime*6
  character(len=3) :: anspec


! Determine current calendar date, needed for the file name
!**********************************************************

  jul=bdate+real(itime,kind=dp)/86400._dp

      call caldate(jul,jjjjmmdd,ihmmss)
      write(adate,'(i8.8)') jjjjmmdd
      write(atime,'(i6.6)') ihmmss
      write(unitdates,'(a)') adate//atime


! For forward simulations, output fields have dimension MAXSPEC,
! for backward simulations, output fields have dimension MAXPOINT.
! Thus, make loops either about nspec, or about numpoint
!*****************************************************************

  if (ldirect.eq.1) then
    do ks=1,nspec
      do kp=1,maxpointspec_act
        tot_mu(ks,kp)=1
      end do
    end do
  else
    do ks=1,nspec
      do kp=1,maxpointspec_act
        tot_mu(ks,kp)=xmass(kp,ks)
      end do
    end do
  endif

!*******************************************************************
! Compute air density: sufficiently accurate to take it
! from coarse grid at some time
! Determine center altitude of output layer, and interpolate density
! data to that altitude
!*******************************************************************

  do kz=1,numzgrid
    if (kz.eq.1) then
      halfheight=outheight(1)/2.
    else
      halfheight=(outheight(kz)+outheight(kz-1))/2.
    endif
    do kzz=2,nz
      if ((height(kzz-1).lt.halfheight).and. &
           (height(kzz).gt.halfheight)) goto 46
    end do
46   kzz=max(min(kzz,nz),2)
    dz1=halfheight-height(kzz-1)
    dz2=height(kzz)-halfheight
    dz=dz1+dz2
    do jy=0,numygrid-1
      do ix=0,numxgrid-1
!        xl=outlon0+real(ix)*dxout
!        yl=outlat0+real(jy)*dyout
        xl=out_xm0+float(ix)*dxout
        yl=out_ym0+float(jy)*dyout
!        xl=(xl-xlon0)/dx
!        yl=(yl-ylat0)/dx
        xl=(xl-xmet0)/dx
        yl=(yl-ymet0)/dy
        iix=max(min(nint(xl),nxmin1),0)
        jjy=max(min(nint(yl),nymin1),0)
        densityoutgrid(ix,jy,kz)=(rho(iix,jjy,kzz,2)*dz1+ &
             rho(iix,jjy,kzz-1,2)*dz2)/dz
      end do
    end do
  end do

  do i=1,numreceptor
    xl=xreceptor(i)
    yl=yreceptor(i)
    iix=max(min(nint(xl),nxmin1),0)
    jjy=max(min(nint(yl),nymin1),0)
    densityoutrecept(i)=rho(iix,jjy,1,2)
  end do

  ! Output is different for forward and backward simulations
  do kz=1,numzgrid
    do jy=0,numygrid-1
      do ix=0,numxgrid-1
        if (ldirect.eq.1) then
          factor3d(ix,jy,kz)=1.e12/volume(ix,jy,kz)/outnum
        else
          factor3d(ix,jy,kz)=real(abs(loutaver))/outnum
        endif
      end do
    end do
  end do

  !*********************************************************************
  ! Determine the standard deviation of the mean concentration or mixing
  ! ratio (uncertainty of the output) and the dry and wet deposition
  !*********************************************************************

  gridtotal=0.
  gridsigmatotal=0.
  gridtotalunc=0.
  wetgridtotal=0.
  wetgridsigmatotal=0.
  wetgridtotalunc=0.
  drygridtotal=0.
  drygridsigmatotal=0.
  drygridtotalunc=0.

!*******************************************************************
! Generate output: may be in concentration (ng/m3) or in mixing
! ratio (ppt) or both
! Output either in full grid dump or sparse matrix format
! For backward simulations, the unit is seconds, stored in grid_conc
!*******************************************************************

! Concentration output
!*********************

!      open(53,file=path(1)(1:length(1))//'latlon.txt',form='formatted')
!          open(54,file=path(1)(1:length(1))//'latlon_corner.txt' &
!          ,form='formatted')
!
!!        xnelat=outgrid_nelat
!!        xnelon=outgrid_nelon
!         print*,'before ll_to',outgrid_swlon,outgrid_swlat,outgrid_nelon,outgrid_nelat
!        call ll_to_xymeter_wrf(outgrid_swlon,outgrid_swlat,xsw,ysw)
!        call ll_to_xymeter_wrf(outgrid_nelon,outgrid_nelat,xne,yne)
!         print*,'after ll_to'
!        do jy=1,numygrid
!        do ix=1,numxgrid
!!         tmpx=out_xm0+(ix-1)*dxout
!!         tmpy=out_ym0+(jy-1)*dyout
!          tmpx=out_xm0+(float(ix)-0.5)*dxout
!          tmpy=out_ym0+(float(jy)-0.5)*dyout
!!          print*,'jb','tmpx','tmpy',dxout,dyout,ix,jy
!          call xymeter_to_ll_wrf(tmpx,tmpy,tmplon,tmplat)
!!jb          if(iouttype.eq.0) write(unitoutgrid) tmplon,tmplat
!!          if(iouttype.eq.1) write(unitoutgrid,*) tmplon,tmplat
!        write(53,*) tmplon,tmplat
!!         tmpx=out_xm0+(ix-1-0.5)*dxout
!!         tmpy=out_ym0+(jy-1-0.5)*dyout
!          tmpx=out_xm0+(float(ix)-1.)*dxout
!          tmpy=out_ym0+(float(jy)-1.)*dyout
!!         tmpx=xsw+(xne-xsw)*float(ix-1)/float(numxgrid-1)
!!         tmpy=ysw+(yne-ysw)*float(jy-1)/float(numygrid-1)
!!          print*,'jb2','tmpx','tmpy',dxout,dyout,ix,jy
!!         call xymeter_to_ll_wrf(tmpx,tmpy,tmplon,tmplat)
!          call xymeter_to_ll_wrf_out(tmpx,tmpy,tmplon,tmplat)
!           write(54,*) tmplon,tmplat
!        enddo
!        enddo
!      close(53)
!      close(54)

!     print*,'in grid conc',nspec,iout,adate,atime
  do ks=1,nspec
    if (iouttype.ne.2) then ! Not netcdf output, open the standard files
      ! AD: I don't think this is right, there is no distinction between ascii
      ! or binary files as there is in concoutput_reg...
      write(anspec,'(i3.3)') ks
      if ((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)) then
        if (ldirect.eq.1) then
          open(unitoutgrid,file=path(1)(1:length(1))//'grid_conc_'//adate// &
               atime//'_'//anspec,form='unformatted')
        else
          open(unitoutgrid,file=path(1)(1:length(1))//'grid_time_'//adate// &
               atime//'_'//anspec,form='unformatted')
        endif
        if (iouttype.eq.0) write(unitoutgrid) itime
        if (iouttype.eq.1) write(unitoutgrid,*) itime
      endif
    endif ! iouttype.ne.2

    if ((iout.eq.2).or.(iout.eq.3)) then      ! mixing ratio
      if(iouttype.ne.2) then ! Not netcdf output, open standard file
        ! AD: still the same issue as my previous comment...
        open(unitoutgridppt,file=path(1)(1:length(1))//'grid_pptv_'//adate// &
              atime//'_'//anspec,form='unformatted')
        if (iouttype.eq.0) write(unitoutgridppt) itime
        if (iouttype.eq.1) write(unitoutgridppt,*) itime
      endif
    endif


!     print*,'in grid conc step 2',maxpointspec_act,nageclass,numygrid,numxgrid,numzgrid
  do kp=1,maxpointspec_act
  do nage=1,nageclass

    do jy=0,numygrid-1
      do ix=0,numxgrid-1

  ! WET DEPOSITION
        if ((WETDEP).and.(ldirect.gt.0)) then
            do l=1,nclassunc
              auxgrid(l)=wetgridunc(ix,jy,ks,kp,l,nage)
            end do
            call mean(auxgrid,wetgrid(ix,jy), &
                 wetgridsigma(ix,jy),nclassunc)
  ! Multiply by number of classes to get total concentration
            wetgrid(ix,jy)=wetgrid(ix,jy) &
                 *nclassunc
            wetgridtotal=wetgridtotal+wetgrid(ix,jy)
  ! Calculate standard deviation of the mean
            wetgridsigma(ix,jy)= &
                 wetgridsigma(ix,jy)* &
                 sqrt(real(nclassunc))
            wetgridsigmatotal=wetgridsigmatotal+ &
                 wetgridsigma(ix,jy)
        endif

  ! DRY DEPOSITION
        if ((DRYDEP).and.(ldirect.gt.0)) then
            do l=1,nclassunc
              auxgrid(l)=drygridunc(ix,jy,ks,kp,l,nage)
            end do
            call mean(auxgrid,drygrid(ix,jy), &
                 drygridsigma(ix,jy),nclassunc)
  ! Multiply by number of classes to get total concentration
            drygrid(ix,jy)=drygrid(ix,jy)* &
                 nclassunc
            drygridtotal=drygridtotal+drygrid(ix,jy)
  ! Calculate standard deviation of the mean
            drygridsigma(ix,jy)= &
                 drygridsigma(ix,jy)* &
                 sqrt(real(nclassunc))
125         drygridsigmatotal=drygridsigmatotal+ &
                 drygridsigma(ix,jy)
        endif
  ! CONCENTRATION OR MIXING RATIO
        do kz=1,numzgrid
            do l=1,nclassunc
              auxgrid(l)=gridunc(ix,jy,kz,ks,kp,l,nage)
            end do
            call mean(auxgrid,grid(ix,jy,kz), &
                 gridsigma(ix,jy,kz),nclassunc)
  ! Multiply by number of classes to get total concentration
            grid(ix,jy,kz)= &
                 grid(ix,jy,kz)*nclassunc
!        if (grid(ix,jy,kz).gt.0. ) print*,grid(ix,jy,kz)
            gridtotal=gridtotal+grid(ix,jy,kz)
  ! Calculate standard deviation of the mean
            gridsigma(ix,jy,kz)= &
                 gridsigma(ix,jy,kz)* &
                 sqrt(real(nclassunc))
            gridsigmatotal=gridsigmatotal+ &
                 gridsigma(ix,jy,kz)
        end do
      end do
    end do

  !*******************************************************************
  ! Generate output: may be in concentration (ng/m3) or in mixing
  ! ratio (ppt) or both
  ! Output the position and the values alternated multiplied by
  ! 1 or -1, first line is number of values, number of positions
  ! For backward simulations, the unit is seconds, stored in grid_time
  !*******************************************************************

  if (iouttype.eq.2) then   ! netcdf output
    if (option_verbose.ge.1) then
      write(*,*) 'concoutput_irreg: Calling write_ncconc for main outgrid'
    endif
    call nc_write_output(itime,outnum,ks,kp,nage,tot_mu(ks,kp),densityoutrecept,volume,0) ! 0= nest level
  else  ! binary or ascii output

    ! Concentration output
    !*********************
    if ((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)) then

  ! Wet deposition
         sp_count_i=0
         sp_count_r=0
         sp_fact=-1.
         sp_zer=.true.
         if ((ldirect.eq.1).and.(WETDEP)) then
         do jy=0,numygrid-1
            do ix=0,numxgrid-1
  !oncentraion greater zero
              if (wetgrid(ix,jy).gt.smallnum) then
                 if (sp_zer.eqv..true.) then ! first non zero value
                    sp_count_i=sp_count_i+1
                    sparse_dump_i(sp_count_i)=ix+jy*numxgrid
                    sp_zer=.false.
                    sp_fact=sp_fact*(-1.)
                 endif
                 sp_count_r=sp_count_r+1
                 sparse_dump_r(sp_count_r)= &
                      sp_fact*1.e12*wetgrid(ix,jy)/area(ix,jy)
  !                sparse_dump_u(sp_count_r)=
  !+                1.e12*wetgridsigma(ix,jy,ks,kp,nage)/area(ix,jy)
              else ! concentration is zero
                  sp_zer=.true.
              endif
            end do
         end do
         else
            sp_count_i=0
            sp_count_r=0
         endif
     if (iouttype.eq.0) then
         write(unitoutgrid) sp_count_i
         write(unitoutgrid) (sparse_dump_i(i),i=1,sp_count_i)
         write(unitoutgrid) sp_count_r
         write(unitoutgrid) (sparse_dump_r(i),i=1,sp_count_r)
     endif
     if (iouttype.eq.1) then
         write(unitoutgrid,*) sp_count_i
         write(unitoutgrid,*) (sparse_dump_i(i),i=1,sp_count_i)
         write(unitoutgrid,*) sp_count_r
         write(unitoutgrid,*) (sparse_dump_r(i),i=1,sp_count_r)
     endif
  !       write(unitoutgrid) sp_count_u
  !       write(unitoutgrid) (sparse_dump_u(i),i=1,sp_count_r)

  ! Dry deposition
         sp_count_i=0
         sp_count_r=0
         sp_fact=-1.
         sp_zer=.true.
         if ((ldirect.eq.1).and.(DRYDEP)) then
          do jy=0,numygrid-1
            do ix=0,numxgrid-1
              if (drygrid(ix,jy).gt.smallnum) then
                 if (sp_zer.eqv..true.) then ! first non zero value
                    sp_count_i=sp_count_i+1
                    sparse_dump_i(sp_count_i)=ix+jy*numxgrid
                    sp_zer=.false.
                    sp_fact=sp_fact*(-1.)
                 endif
                 sp_count_r=sp_count_r+1
                 sparse_dump_r(sp_count_r)= &
                      sp_fact* &
                      1.e12*drygrid(ix,jy)/area(ix,jy)
  !                sparse_dump_u(sp_count_r)=
  !+                1.e12*drygridsigma(ix,jy,ks,kp,nage)/area(ix,jy)
              else ! concentration is zero
                  sp_zer=.true.
              endif
            end do
          end do
         else
            sp_count_i=0
            sp_count_r=0
         endif
     if (iouttype.eq.0) then
         write(unitoutgrid) sp_count_i
         write(unitoutgrid) (sparse_dump_i(i),i=1,sp_count_i)
         write(unitoutgrid) sp_count_r
         write(unitoutgrid) (sparse_dump_r(i),i=1,sp_count_r)
     endif
     if (iouttype.eq.1) then
         write(unitoutgrid,*) sp_count_i
         write(unitoutgrid,*) (sparse_dump_i(i),i=1,sp_count_i)
         write(unitoutgrid,*) sp_count_r
         write(unitoutgrid,*) (sparse_dump_r(i),i=1,sp_count_r)
     endif

  !       write(*,*) sp_count_u
  !       write(unitoutgrid) (sparse_dump_u(i),i=1,sp_count_r)

  ! Concentrations
         sp_count_i=0
         sp_count_r=0
         sp_fact=-1.
         sp_zer=.true.
          do kz=1,numzgrid
            do jy=0,numygrid-1
              do ix=0,numxgrid-1
                if (grid(ix,jy,kz).gt.smallnum) then
                  if (sp_zer.eqv..true.) then ! first non zero value
                    sp_count_i=sp_count_i+1
                    sparse_dump_i(sp_count_i)= &
                         ix+jy*numxgrid+kz*numxgrid*numygrid
                    sp_zer=.false.
                    sp_fact=sp_fact*(-1.)
                   endif
                   sp_count_r=sp_count_r+1
                   sparse_dump_r(sp_count_r)= &
                        sp_fact* &
                        grid(ix,jy,kz)* &
                        factor3d(ix,jy,kz)/tot_mu(ks,kp)
  !                 if ((factor(ix,jy,kz)/tot_mu(ks,kp)).eq.0)
  !    +              write (*,*) factor(ix,jy,kz),tot_mu(ks,kp),ks,kp
  !                sparse_dump_u(sp_count_r)=
  !+               ,gridsigma(ix,jy,kz,ks,kp,nage)*
  !+               factor(ix,jy,kz)/tot_mu(ks,kp)
              else ! concentration is zero
                  sp_zer=.true.
              endif
              end do
            end do
          end do
     if (iouttype.eq.0) then
         write(unitoutgrid) sp_count_i
         write(unitoutgrid) (sparse_dump_i(i),i=1,sp_count_i)
         write(unitoutgrid) sp_count_r
         write(unitoutgrid) (sparse_dump_r(i),i=1,sp_count_r)
     endif
     if (iouttype.eq.1) then
         write(unitoutgrid,*) sp_count_i
         write(unitoutgrid,*) (sparse_dump_i(i),i=1,sp_count_i)
         write(unitoutgrid,*) sp_count_r
         write(unitoutgrid,*) (sparse_dump_r(i),i=1,sp_count_r)
     endif
  !       write(unitoutgrid) sp_count_u
  !       write(unitoutgrid) (sparse_dump_u(i),i=1,sp_count_r)



    endif !  concentration output

    ! Mixing ratio output
    !********************

    if ((iout.eq.2).or.(iout.eq.3)) then      ! mixing ratio

    ! Wet deposition
         sp_count_i=0
         sp_count_r=0
         sp_fact=-1.
         sp_zer=.true.
         if ((ldirect.eq.1).and.(WETDEP)) then
          do jy=0,numygrid-1
            do ix=0,numxgrid-1
                if (wetgrid(ix,jy).gt.smallnum) then
                  if (sp_zer.eqv..true.) then ! first non zero value
                    sp_count_i=sp_count_i+1
                    sparse_dump_i(sp_count_i)= &
                         ix+jy*numxgrid
                    sp_zer=.false.
                    sp_fact=sp_fact*(-1.)
                 endif
                 sp_count_r=sp_count_r+1
                 sparse_dump_r(sp_count_r)= &
                      sp_fact* &
                      1.e12*wetgrid(ix,jy)/area(ix,jy)
  !                sparse_dump_u(sp_count_r)=
  !    +            ,1.e12*wetgridsigma(ix,jy,ks,kp,nage)/area(ix,jy)
              else ! concentration is zero
                  sp_zer=.true.
              endif
            end do
          end do
         else
           sp_count_i=0
           sp_count_r=0
         endif
     if (iouttype.eq.0) then
         write(unitoutgridppt) sp_count_i
         write(unitoutgridppt) (sparse_dump_i(i),i=1,sp_count_i)
         write(unitoutgridppt) sp_count_r
         write(unitoutgridppt) (sparse_dump_r(i),i=1,sp_count_r)
     endif
     if (iouttype.eq.1) then
         write(unitoutgridppt,*) sp_count_i
         write(unitoutgridppt,*) (sparse_dump_i(i),i=1,sp_count_i)
         write(unitoutgridppt,*) sp_count_r
         write(unitoutgridppt,*) (sparse_dump_r(i),i=1,sp_count_r)
     endif
  !       write(unitoutgridppt) sp_count_u
  !       write(unitoutgridppt) (sparse_dump_u(i),i=1,sp_count_r)

  ! Dry deposition
         sp_count_i=0
         sp_count_r=0
         sp_fact=-1.
         sp_zer=.true.
         if ((ldirect.eq.1).and.(DRYDEP)) then
          do jy=0,numygrid-1
            do ix=0,numxgrid-1
                if (drygrid(ix,jy).gt.smallnum) then
                  if (sp_zer.eqv..true.) then ! first non zero value
                    sp_count_i=sp_count_i+1
                    sparse_dump_i(sp_count_i)= &
                         ix+jy*numxgrid
                    sp_zer=.false.
                    sp_fact=sp_fact*(-1)
                 endif
                 sp_count_r=sp_count_r+1
                 sparse_dump_r(sp_count_r)= &
                      sp_fact* &
                      1.e12*drygrid(ix,jy)/area(ix,jy)
  !                sparse_dump_u(sp_count_r)=
  !    +            ,1.e12*drygridsigma(ix,jy,ks,kp,nage)/area(ix,jy)
              else ! concentration is zero
                  sp_zer=.true.
              endif
            end do
          end do
         else
           sp_count_i=0
           sp_count_r=0
         endif
     if (iouttype.eq.0) then
         write(unitoutgridppt) sp_count_i
         write(unitoutgridppt) (sparse_dump_i(i),i=1,sp_count_i)
         write(unitoutgridppt) sp_count_r
         write(unitoutgridppt) (sparse_dump_r(i),i=1,sp_count_r)
     endif
     if (iouttype.eq.1) then
         write(unitoutgridppt,*) sp_count_i
         write(unitoutgridppt,*) (sparse_dump_i(i),i=1,sp_count_i)
         write(unitoutgridppt,*) sp_count_r
         write(unitoutgridppt,*) (sparse_dump_r(i),i=1,sp_count_r)
     endif
  !       write(unitoutgridppt) sp_count_u
  !       write(unitoutgridppt) (sparse_dump_u(i),i=1,sp_count_r)


  ! Mixing ratios
         sp_count_i=0
         sp_count_r=0
         sp_fact=-1.
         sp_zer=.true.
          do kz=1,numzgrid
            do jy=0,numygrid-1
              do ix=0,numxgrid-1
                if (grid(ix,jy,kz).gt.smallnum) then
                  if (sp_zer.eqv..true.) then ! first non zero value
                    sp_count_i=sp_count_i+1
                    sparse_dump_i(sp_count_i)= &
                         ix+jy*numxgrid+kz*numxgrid*numygrid
                    sp_zer=.false.
                    sp_fact=sp_fact*(-1.)
                 endif
                 sp_count_r=sp_count_r+1
                 sparse_dump_r(sp_count_r)= &
                      sp_fact* &
                      1.e12*grid(ix,jy,kz) &
                      /volume(ix,jy,kz)/outnum* &
                      weightair/weightmolar(ks)/densityoutgrid(ix,jy,kz)
  !                sparse_dump_u(sp_count_r)=
  !+              ,1.e12*gridsigma(ix,jy,kz,ks,kp,nage)/volume(ix,jy,kz)/
  !+              outnum*weightair/weightmolar(ks)/
  !+              densityoutgrid(ix,jy,kz)
              else ! concentration is zero
                  sp_zer=.true.
              endif
              end do
            end do
          end do
     if (iouttype.eq.0) then
         write(unitoutgridppt) sp_count_i
         write(unitoutgridppt) (sparse_dump_i(i),i=1,sp_count_i)
         write(unitoutgridppt) sp_count_r
         write(unitoutgridppt) (sparse_dump_r(i),i=1,sp_count_r)
     endif
     if (iouttype.eq.1) then
         write(unitoutgridppt,*) sp_count_i
         write(unitoutgridppt,*) (sparse_dump_i(i),i=1,sp_count_i)
         write(unitoutgridppt,*) sp_count_r
         write(unitoutgridppt,*) (sparse_dump_r(i),i=1,sp_count_r)
     endif
  !       write(unitoutgridppt) sp_count_u
  !       write(unitoutgridppt) (sparse_dump_u(i),i=1,sp_count_r)

    endif ! output for ppt

  endif ! iouttype.eq.2

  end do
  end do

    if((iouttype.eq.0).or.(iouttype.eq.1)) then ! binary or ascii output
      close(unitoutgridppt)
      close(unitoutgrid)
    endif

  end do

  if (gridtotal.gt.0.) gridtotalunc=gridsigmatotal/gridtotal
  if (wetgridtotal.gt.0.) wetgridtotalunc=wetgridsigmatotal/ &
       wetgridtotal
  if (drygridtotal.gt.0.) drygridtotalunc=drygridsigmatotal/ &
       drygridtotal

  ! Dump values at receptors (binary or ascii only):
  ! When iouttype==2, this is handled by the default netcdf output routines, 
  ! no extra call is needed here.
  if ( (iouttype.eq.0 .or. iouttype.eq.1) .and. numreceptor.gt.0 ) then
    ! Dump of receptor mixing ratios
    if (iout.eq.2 .or. iout.eq.3) then ! mix. rat.
      write(unitoutreceptppt) itime
      do ks=1,nspec
        write(unitoutreceptppt) (1.e12*creceptor(i,ks)/outnum* &
              weightair/weightmolar(ks)/densityoutrecept(i),i=1,numreceptor)
      end do
    endif
    ! Dump of receptor concentrations
    if (iout.eq.1 .or. iout.eq.3 .or. iout.eq.5) then
      write(unitoutrecept) itime
      do ks=1,nspec
        write(unitoutrecept) (1.e12*creceptor(i,ks)/outnum, i=1,numreceptor)
      end do
    endif
  endif ! Dump values at receptors

  do ks=1,nspec
  do kp=1,maxpointspec_act
    do i=1,numreceptor
      creceptor(i,ks)=0.
    end do
    do jy=0,numygrid-1
      do ix=0,numxgrid-1
        do l=1,nclassunc
          do nage=1,nageclass
            do kz=1,numzgrid
              gridunc(ix,jy,kz,ks,kp,l,nage)=0.
            end do
          end do
        end do
      end do
    end do
  end do
  end do


end subroutine concoutput_irreg

