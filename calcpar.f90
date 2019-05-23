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
      subroutine calcpar(n,uuh,vvh,pvh)
!                        i  i   i   o
!*******************************************************************************
!                                                                              *
!     Computation of several boundary layer parameters needed for the          *
!     dispersion calculation and calculation of dry deposition velocities.     *
!     All parameters are calculated over the entire grid.                      *
!                                                                              *
!     Note:  This is the FLEXPART_WRF version of subroutine calcpar.           *
!                                                                              *
!     Author: A. Stohl                                                         *
!                                                                              *
!     21 May 1995                                                              *
!                                                                              *
! ------------------------------------------------------------------           *
!     Petra Seibert, Feb 2000:                                                 *
!     convection scheme:                                                       *
!     new variables in call to richardson                                      *
!                                                                              *
!     Changes, Bernd C. Krueger, Feb. 2001:
!        Variables tth and qvh (on eta coordinates) in common block
!                                                                              *
!     17 Oct 2005 - R. Easter - added ierr in call to richardson               *
!     18 Oct 2005 - J. Fast - limit ustar to < 5.0 m/s                         *
!     -- Oct 2005 - R. Easter - use xy_to_ll_wrf to get latitude               *
!             use pph for calculating zlev                                     *
!             pass level-2 pph directly to obukhov                             *
!     11 Nov 2005 - R. Easter - changed name of "xy to latlon" routine         *
!     15 Nov 2005 - R. Easter - pass pplev to richardson instead of akz,bkz    *
!    July 2012: J. Brioude: coded in fortran 90 and parallelized               *
!   2015-05-27  A. Dingwell - Bugfix: Charnock's relation should apply to water*
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
! n                  temporal index for meteorological fields (1 to 3)         *
!                                                                              *
! Constants:                                                                   *
!                                                                              *
!                                                                              *
! Functions:                                                                   *
! scalev             computation of ustar                                      *
! obukhov            computatio of Obukhov length                              *
!                                                                              *
!*******************************************************************************

  use par_mod
  use com_mod

  implicit none

  integer :: n,ix,jy,i,kz,lz,kzmin,ierr
  real :: ttlev(nuvzmax),qvlev(nuvzmax),obukhov,scalev,ol,hmixplus
  real :: ulev(nuvzmax),vlev(nuvzmax),ew,rh,vd(maxspec),subsceff,ylat
  real :: altmin,tvold,pold,zold,pint,tv,zlev(nuvzmax)
   real(kind=4) :: uuh(0:nxmax-1,0:nymax-1,nuvzmax)
   real(kind=4) :: vvh(0:nxmax-1,0:nymax-1,nuvzmax)

! real :: uuh(0:nxmax-1,0:nymax-1,nuvzmax)
! real :: vvh(0:nxmax-1,0:nymax-1,nuvzmax)
  real :: pvh(0:nxmax-1,0:nymax-1,nuvzmax)
  real,parameter :: const=r_air/ga

  real :: xlon,dumx,dumy,dumxb,dumyb,pplev(nuvzmax),hmixdummy

! Loop over entire grid
!**********************
!      ientry = ientry + 1

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(i,ix,jy,kz,lz,kzmin,tvold,pold,zold,zlev,tv,pint, &
!$OMP rh,ierr,subsceff,ulev,vlev,pplev,ttlev,qvlev,ol,altmin,xlon,ylat )
!$OMP DO
      do jy=0,nymin1

        do ix=0,nxmin1

! Set minimum height for tropopause
!**********************************

! FLEXPART_WRF - use this routine to get lat,lon
!       ylat=ylat0+real(jy)*dy
        call xyindex_to_ll_wrf( 0, real(ix), real(jy), xlon, ylat )

!       if ( ((ix.eq.0) .or. (ix.eq.nxmin1) .or. 
!    &                       (ix.eq.nxmin1/2)) .and.
!    &       ((jy.eq.0) .or. (jy.eq.nymin1) .or. 
!    &                       (jy.eq.nymin1/2)) ) then
!           if (ientry .eq. 1) then
!               write(*,'(a,2i4,2f12.5)') 
!    &              'calcpar i,j, xlon,ylat', ix, jy, xlon, ylat
!               write(*,'(a, 8x,2f12.5)') 
!    &              '             dlon,dlat', 
!    &              (xlon-xlon2d(ix,jy)), (ylat-ylat2d(ix,jy))
!               call ll_to_xyindex_wrf(
!    &              xlon2d(ix,jy), ylat2d(ix,jy), dumx, dumy )
!               write(*,'(a, 8x,2f12.5)') 
!    &              '             dxkm,dykm', 
!    &              ((dumx-ix)*dx*1.0e-3), ((dumy-jy)*dy*1.0e-3) 
!
!               if ((ix .eq. 0) .and. (jy .eq. 0)) then
!                  dumxb = 2.33
!                  dumyb = 3.44
!                  call xyindex_to_ll_wrf( 0, dumxb, dumyb, dumx, dumy )
!                  call ll_to_xyindex_wrf( dumx, dumy, dumx, dumy )
!                  write(*,'(a,2f5.2,2f12.5)') 
!    &                'xi,yj,     dxkm,dykm', dumxb, dumyb,
!    &                ((dumx-dumxb)*dx*1.0e-3), ((dumy-dumyb)*dy*1.0e-3)
!                  dumxb = 4.55
!                  dumyb = 6.77
!                  call xyindex_to_ll_wrf( 0, dumxb, dumyb, dumx, dumy )
!                  call ll_to_xyindex_wrf( dumx, dumy, dumx, dumy )
!                  write(*,'(a,2f5.2,2f12.5)') 
!    &                'xi,yj,     dxkm,dykm', dumxb, dumyb,
!    &                ((dumx-dumxb)*dx*1.0e-3), ((dumy-dumyb)*dy*1.0e-3)
!               end if
!
!           end if
!       end if

        if ((ylat.ge.-20.).and.(ylat.le.20.)) then
          altmin = 5000.
        else
          if ((ylat.gt.20.).and.(ylat.lt.40.)) then
            altmin=2500.+(40.-ylat)*125.
          else if ((ylat.gt.-40.).and.(ylat.lt.-20.)) then
            altmin=2500.+(40.+ylat)*125.
          else
            altmin=2500.
          endif
        endif

! 1) Calculation of friction velocity
!************************************
          if ( (.not.strswitch)) then
          ustar(ix,jy,1,n)=scalev(ps(ix,jy,1,n),tt2(ix,jy,1,n), &
          td2(ix,jy,1,n),surfstr(ix,jy,1,n))
          endif
          if (ustar(ix,jy,1,n).le.1.e-8) ustar(ix,jy,1,n)=1.e-8
! FLEXPART_WRF - limit ustar
          if (ustar(ix,jy,1,n).ge.5.0)   ustar(ix,jy,1,n)=5.0

! 2) Calculation of inverse Obukhov length scale
!***********************************************

! FLEXPART_WRF - pass k=2 pressure directly
!         ol=obukhov(ps(ix,jy,1,n),tt2(ix,jy,1,n),td2(ix,jy,1,n),
!    +    tth(ix,jy,2,n),ustar(ix,jy,1,n),sshf(ix,jy,1,n),akm,bkm)
          ol=obukhov(ps(ix,jy,1,n),tt2(ix,jy,1,n),td2(ix,jy,1,n), &
          tth(ix,jy,2,n),ustar(ix,jy,1,n),sshf(ix,jy,1,n), &
          pph(ix,jy,2,n) )
          if (ol.ne.0.) then
            oli(ix,jy,1,n)=1./ol
          else
            oli(ix,jy,1,n)=99999.
          endif


! 3) Calculation of convective velocity scale and mixing height
!**************************************************************
            
          do i=1,nuvz
            ulev(i) =uuh(ix,jy,i)
            vlev(i) =vvh(ix,jy,i)
            pplev(i)=pph(ix,jy,i,n)
            ttlev(i)=tth(ix,jy,i,n)
            qvlev(i)=qvh(ix,jy,i,n)
            zlev(i)=0.5*(zzh(ix,jy,i+1,n)+zzh(ix,jy,i,n))-zzh(ix,jy,1,n)
           enddo
! FLEXPART_WRF - use  & check ierr argument
! FLEXPART_WRF - pass pplev instead of akz,bkz
!         call richardson(ps(ix,jy,1,n),ustar(ix,jy,1,n),ttlev,qvlev,
!    +    ulev,vlev,nuvz,akz,bkz,sshf(ix,jy,1,n),tt2(ix,jy,1,n),
!    +    td2(ix,jy,1,n),hmix(ix,jy,1,n),wstar(ix,jy,1,n),hmixplus)
          call richardson(ps(ix,jy,1,n),ustar(ix,jy,1,n),ttlev,qvlev, &
          ulev,vlev,nuvz,  pplev,sshf(ix,jy,1,n),tt2(ix,jy,1,n), &
          td2(ix,jy,1,n),hmix(ix,jy,1,n),wstar(ix,jy,1,n),hmixplus, &
!         td2(ix,jy,1,n),hmixdummy,wstar(ix,jy,1,n),hmixplus, &
          ierr,sfc_option )
!JB
! no reflec
!         hmix(ix,jy,1,n)=5000.

          if (ierr .gt. 0) then
              write(*,9500) 'warning', ix, jy
          else if (ierr .lt. 0) then
              write(*,9500) 'failure', ix, jy
              stop
          end if
9500      format( 'calcpar - richardson ', a, ' - ix,jy=', 2i5 )


          if(lsubgrid.eq.1) then
            subsceff=min(excessoro(ix,jy),hmixplus)
!           subsceff=hmixplus
          else
            subsceff=0.
          endif
!
! CALCULATE HMIX EXCESS ACCORDING TO SUBGRIDSCALE VARIABILITY AND STABILITY
!
          hmix(ix,jy,1,n)=hmix(ix,jy,1,n)+subsceff
!         print*,'hmix',hmix(ix,jy,1,n),subsceff
          hmix(ix,jy,1,n)=max(hmixmin,hmix(ix,jy,1,n)) ! set minimum PBL height
          hmix(ix,jy,1,n)=min(hmixmax,hmix(ix,jy,1,n)) ! set maximum PBL height

! 4) Calculation of dry deposition velocities
!********************************************

          if (DRYDEP) then
            ! AD: Bugfix: Land use indexes have changed, water is now class 7
            ! (indexing starts at 1)
            !z0(4)=0.016*ustar(ix,jy,1,n)*ustar(ix,jy,1,n)/ga
            !z0(9)=0.016*ustar(ix,jy,1,n)*ustar(ix,jy,1,n)/ga
            z0(7)=0.016*ustar(ix,jy,1,n)*ustar(ix,jy,1,n)/ga

! Calculate relative humidity at surface
!***************************************
            rh=ew(td2(ix,jy,1,n))/ew(tt2(ix,jy,1,n))

            call getvdep(n,ix,jy,ustar(ix,jy,1,n), &
            tt2(ix,jy,1,n),ps(ix,jy,1,n),1./oli(ix,jy,1,n), &
            ssr(ix,jy,1,n),rh,lsprec(ix,jy,1,n)+convprec(ix,jy,1,n),vd)

            do i=1,nspec
              vdep(ix,jy,i,n)=vd(i)
            enddo
          endif

!******************************************************
! Calculate height of thermal tropopause (Hoinka, 1997)
!******************************************************

! 1) Calculate altitudes of ECMWF model levels
!*********************************************

          tvold=tt2(ix,jy,1,n)*(1.+0.378*ew(td2(ix,jy,1,n))/ &
                                         ps(ix,jy,1,n))
          pold=ps(ix,jy,1,n)
          zold=0.
! FLEXPART_WRF - set zlev(1)
          zlev(1)=zold
          do kz=2,nuvz
! FLEXPART_WRF - use pph for pressure
!           pint=akz(kz)+bkz(kz)*ps(ix,jy,1,n)  ! pressure on model layers
            pint=pph(ix,jy,kz,n)  ! pressure on model layers
            tv=tth(ix,jy,kz,n)*(1.+0.608*qvh(ix,jy,kz,n))

            if (abs(tv-tvold).gt.0.2) then
             zlev(kz)=zold+const*log(pold/pint)*(tv-tvold)/log(tv/tvold)
            else
              zlev(kz)=zold+const*log(pold/pint)*tv
            endif
            tvold=tv
            pold=pint
            zold=zlev(kz)
            enddo

! 2) Define a minimum level kzmin, from which upward the tropopause is
!    searched for. This is to avoid inversions in the lower troposphere
!    to be identified as the tropopause
!************************************************************************

      do kz=1,nuvz
        if (zlev(kz).ge.altmin) then
          kzmin=kz
          goto 45
        endif
      end do
45    continue

! 3) Search for first stable layer above minimum height that fulfills the
!    thermal tropopause criterion
!************************************************************************

      do kz=kzmin,nuvz
        do lz=kz+1,nuvz
          if ((zlev(lz)-zlev(kz)).gt.2000.) then
            if (((tth(ix,jy,kz,n)-tth(ix,jy,lz,n))/ &
                 (zlev(lz)-zlev(kz))).lt.0.002) then
              tropopause(ix,jy,1,n)=zlev(kz)
              goto 51
            endif
            goto 50
          endif
        end do
50      continue
      end do
51    continue


    end do
  end do
!$OMP END DO
!$OMP END PARALLEL

! Calculation of potential vorticity on 3-d grid, if plume trajectory mode is used
!*********************************************************************************

      if ((iout.eq.4).or.(iout.eq.5).or.(mdomainfill.eq.2)) then
        call calcpv(n,uuh,vvh,pvh)
      endif


end subroutine calcpar

