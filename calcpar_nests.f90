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
      subroutine calcpar_nests(n,uuhn,vvhn,pvhn)
!                              i  i    i    o
!*******************************************************************************
!                                                                              *
!     Computation of several boundary layer parameters needed for the          *
!     dispersion calculation and calculation of dry deposition velocities.     *
!     All parameters are calculated over the entire grid.                      *
!     This routine is similar to calcpar, but is used for the nested grids.    *
!                                                                              *
!     Note:  This is the FLEXPART_WRF version of subroutine calcpar.           *
!                                                                              *
!     Author: A. Stohl                                                         *
!                                                                              *
!     8 February 1999                                                          *
!                                                                              *
! ------------------------------------------------------------------           *
!     Petra Seibert, Feb 2000:                                                 *
!     convection scheme:                                                       *
!     new variables in call to richardson                                      *
!                                                                              *
!     Changes, Bernd C. Krueger, Feb. 2001:                                    *
!        Variables tth and qvh (on eta coordinates) in common block            *
!                                                                              *
!     14 Nov 2005 - R. Easter -                                                *
!          use xyindex_to_ll_wrf to get latitude                               *
!          limit ustar to < 5.0 m/s                                            *
!          added ierr in call to richardson                                    *
!          use pph for calculating zlev                                        *
!          pass level-2 pph directly to obukhov                                *
!     15 Nov 2005 - R. Easter - pass pplev to richardson instead of akz,bkz    *
!    Jul 2012: J. Brioude: coded in fortran90 and parallelized                 *
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

  integer :: n,ix,jy,i,l,kz,lz,kzmin,ierr
  real :: ttlev(nuvzmax),qvlev(nuvzmax),obukhov,scalev,ol,hmixplus
  real :: pplev(nuvzmax),xlon
  real :: ulev(nuvzmax),vlev(nuvzmax),ew,rh,vd(maxspec),subsceff,ylat
  real :: altmin,tvold,pold,zold,pint,tv,zlev(nuvzmax)
  real(kind=4) :: uuhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
  real(kind=4) :: vvhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
  real :: pvhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
  real,parameter :: const=r_air/ga


! Loop over all nests
!********************
!     ientry = ientry + 1

      do l=1,numbnests

! Loop over entire grid
!**********************
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(i,ix,jy,kz,lz,kzmin,tvold,pold,zold,zlev,tv,pint, &
!$OMP rh,ierr,subsceff,ulev,vlev,pplev,ttlev,qvlev,ol,altmin,xlon,ylat )
!$OMP DO
      do jy=0,nyn(l)-1

        do ix=0,nxn(l)-1

! Set minimum height for tropopause
!**********************************

! FLEXPART_WRF - use this routine to get lat,lon
!       ylat=ylat0n(l)+real(jy)*dyn(l)
        call xyindex_to_ll_wrf( l, real(ix), real(jy), xlon, ylat )

!       if ( ((ix.eq.0) .or. (ix.eq.(nxn(l)-1)) .or. 
!    &                       (ix.eq.(nxn(l)-1)/2)) .and.
!    &       ((jy.eq.0) .or. (jy.eq.(nyn(l)-1)) .or. 
!    &                       (jy.eq.(nyn(l)-1)/2)) ) then
!           if (ientry .eq. 1) then
!               write(*,'(a,3i4,2f12.5)') 
!    &              'calcpar_nests l,i,j, xlon,ylat', 
!    &              l, ix, jy, xlon, ylat
!               write(*,'(a,12x,2f12.5)') 
!    &              '                     dlon,dlat', 
!    &              (xlon-xlon2dn(ix,jy,l)), (ylat-ylat2dn(ix,jy,l))
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

          ustarn(ix,jy,1,n,l)=scalev(psn(ix,jy,1,n,l),tt2n(ix,jy,1,n,l), &
          td2n(ix,jy,1,n,l),surfstrn(ix,jy,1,n,l))
          endif
! FLEXPART_WRF - limit ustar
          if (ustarn(ix,jy,1,n,l).le.1.e-8) ustarn(ix,jy,1,n,l)=1.e-8
          if (ustarn(ix,jy,1,n,l).ge.5.0)   ustarn(ix,jy,1,n,l)=5.0

! 2) Calculation of inverse Obukhov length scale
!***********************************************

! FLEXPART_WRF - pass k=2 pressure directly
!         ol=obukhov(psn(ix,jy,1,n,l),tt2n(ix,jy,1,n,l),
!    +    td2n(ix,jy,1,n,l),tthn(ix,jy,2,n,l),ustarn(ix,jy,1,n,l),
!    +    sshfn(ix,jy,1,n,l),akm,bkm)
          ol=obukhov(psn(ix,jy,1,n,l),tt2n(ix,jy,1,n,l), &
          td2n(ix,jy,1,n,l),tthn(ix,jy,2,n,l),ustarn(ix,jy,1,n,l), &
          sshfn(ix,jy,1,n,l),pphn(ix,jy,2,n,l) )
          if (ol.ne.0.) then
            olin(ix,jy,1,n,l)=1./ol
          else
            olin(ix,jy,1,n,l)=99999.
          endif


! 3) Calculation of convective velocity scale and mixing height
!**************************************************************

          do i=1,nuvz
            ulev(i) =uuhn(ix,jy,i,l)
            vlev(i) =vvhn(ix,jy,i,l)
            pplev(i)=pphn(ix,jy,i,n,l)
            ttlev(i)=tthn(ix,jy,i,n,l)
            qvlev(i)=qvhn(ix,jy,i,n,l)
      end do

! FLEXPART_WRF - use  & check ierr argument
! FLEXPART_WRF - pass pplev instead of akz,bkz
!         call richardson(psn(ix,jy,1,n,l),ustarn(ix,jy,1,n,l),ttlev,
!    +    qvlev,ulev,vlev,nuvz,akz,bkz,sshfn(ix,jy,1,n,l),
!    +    tt2n(ix,jy,1,n,l),td2n(ix,jy,1,n,l),hmixn(ix,jy,1,n,l),
!    +    wstarn(ix,jy,1,n,l),hmixplus)
          call richardson(psn(ix,jy,1,n,l),ustarn(ix,jy,1,n,l),ttlev, &
          qvlev,ulev,vlev,nuvz,  pplev,sshfn(ix,jy,1,n,l), &
          tt2n(ix,jy,1,n,l),td2n(ix,jy,1,n,l),hmixn(ix,jy,1,n,l), &
          wstarn(ix,jy,1,n,l),hmixplus,ierr,sfc_option)

          if (ierr .gt. 0) then
              write(*,9500) 'warning', l, ix, jy
          else if (ierr .lt. 0) then
              write(*,9500) 'failure', l, ix, jy
              stop
          end if
9500      format( 'calcpar_nests - richardson ', a, ' - l,ix,jy=', 3i5 )

          if(lsubgrid.eq.1) then
            subsceff=min(excessoron(ix,jy,l),hmixplus)
          else
            subsceff=0
          endif
!
! CALCULATE HMIX EXCESS ACCORDING TO SUBGRIDSCALE VARIABILITY AND STABILITY
!
          hmixn(ix,jy,1,n,l)=hmixn(ix,jy,1,n,l)+subsceff
          hmixn(ix,jy,1,n,l)=max(hmixmin,hmixn(ix,jy,1,n,l)) ! minim PBL height
          hmixn(ix,jy,1,n,l)=min(hmixmax,hmixn(ix,jy,1,n,l)) ! maxim PBL height


! 4) Calculation of dry deposition velocities
!********************************************

          if (DRYDEP) then
            z0(4)=0.016*ustarn(ix,jy,1,n,l)*ustarn(ix,jy,1,n,l)/ga
            z0(9)=0.016*ustarn(ix,jy,1,n,l)*ustarn(ix,jy,1,n,l)/ga

! Calculate relative humidity at surface
!***************************************
            rh=ew(td2n(ix,jy,1,n,l))/ew(tt2n(ix,jy,1,n,l))

            call getvdep_nests(n,ix,jy,ustarn(ix,jy,1,n,l), &
            tt2n(ix,jy,1,n,l),psn(ix,jy,1,n,l),1./olin(ix,jy,1,n,l), &
            ssrn(ix,jy,1,n,l),rh,lsprecn(ix,jy,1,n,l)+ &
            convprecn(ix,jy,1,n,l),vd,l)

            do i=1,nspec
              vdepn(ix,jy,i,n,l)=vd(i)
            enddo
          endif

!******************************************************
! Calculate height of thermal tropopause (Hoinka, 1997)
!******************************************************

! 1) Calculate altitudes of ECMWF model levels
!*********************************************

          tvold=tt2n(ix,jy,1,n,l)*(1.+0.378*ew(td2n(ix,jy,1,n,l))/ &
                                         psn(ix,jy,1,n,l))
          pold=psn(ix,jy,1,n,l)
          zold=0.
! FLEXPART_WRF - set zlev(1)
          zlev(1)=zold
      do kz=2,nuvz
        pint=akz(kz)+bkz(kz)*psn(ix,jy,1,n,l)  ! pressure on model layers
        tv=tthn(ix,jy,kz,n,l)*(1.+0.608*qvhn(ix,jy,kz,n,l))

        if (abs(tv-tvold).gt.0.2) then
         zlev(kz)=zold+const*log(pold/pint)*(tv-tvold)/log(tv/tvold)
        else
          zlev(kz)=zold+const*log(pold/pint)*tv
        endif
        tvold=tv
        pold=pint
        zold=zlev(kz)
      end do

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
            if (((tthn(ix,jy,kz,n,l)-tthn(ix,jy,lz,n,l))/ &
                 (zlev(lz)-zlev(kz))).lt.0.002) then
              tropopausen(ix,jy,1,n,l)=zlev(kz)
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

        if ((iout.eq.4).or.(iout.eq.5)) then
          call calcpv_nests(l,n,uuhn,vvhn,pvhn)
        endif

        enddo


end subroutine calcpar_nests

