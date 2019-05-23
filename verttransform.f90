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

      subroutine verttransform(n,uuh,vvh,wwh,pvh,divh)
!                              i  i   i   i   i
!*******************************************************************************
!                                                                              *
!     This subroutine transforms temperature, dew point temperature and        *
!     wind components from eta to meter coordinates.                           *
!     The vertical wind component is transformed from Pa/s to m/s using        *
!     the conversion factor pinmconv.                                          *
!     In addition, this routine calculates vertical density gradients          *
!     needed for the parameterization of the turbulent velocities.             *
!                                                                              *
!     Note:  This is the FLEXPART_WRF version of subroutine assignland.        *
!            The computational grid is the WRF x-y grid rather than lat-lon.   *
!                                                                              *
!     Author: A. Stohl, G. Wotawa                                              *
!                                                                              *
!     12 August 1996                                                           *
!     Update: 16 January 1998                                                  *
!                                                                              *
!     Major update: 17 February 1999                                           *
!     by G. Wotawa                                                             *
!                                                                              *
!     - Vertical levels for u, v and w are put together                        *
!     - Slope correction for vertical velocity: Modification of calculation    *
!       procedure                                                              *
!                                                                              *
!     Changes, Bernd C. Krueger, Feb. 2001:                                    *
!        Variables tth and qvh (on eta coordinates) from common block          *
!                                                                              *
!     Oct-Nov 2005 - R. Easter - conversion to wrf                             *
!     17 Nov 2005 - R. Easter - terrain correction applied to ww.  There are   *
!            now 3 options, controlled by "method_w_terrain_correction"        *
!                                                                              *
!     11 June 2007,  conversion of tkeh to tke
!     25 June 2007   conversion of ptth to ptt
!     Jan 2012, J Brioude:  modified to handle different wind options and openmp 
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
! nx,ny,nz                        field dimensions in x,y and z direction      *
! uu(0:nxmax,0:nymax,nzmax,2)     wind components in x-direction [m/s]         *
! vv(0:nxmax,0:nymax,nzmax,2)     wind components in y-direction [m/s]         *
! ww(0:nxmax,0:nymax,nzmax,2)     wind components in z-direction [deltaeta/s]  *
! tt(0:nxmax,0:nymax,nzmax,2)     temperature [K]                              *
! pv(0:nxmax,0:nymax,nzmax,2)     potential voriticity (pvu)                   *
! ps(0:nxmax,0:nymax,2)           surface pressure [Pa]                        *
!                                                                              *
!*******************************************************************************

  use par_mod
  use com_mod
!      include 'includepar'
!      include 'includecom'

  implicit none

      integer :: ix,jy,kz,iz,n,kmin,kl,klp,ix1,jy1,ixp,jyp,ixm,jym
! CDA 
      integer :: icloudtop,OMP_GET_THREAD_NUM

      integer :: method_z_compute,aa
      real :: uvzlev(nuvzmax),rhoh(nuvzmax),pinmconv(nzmax)
      real :: ew,pint,tv,tvold,pold,dz1,dz2,dz,ui,vi
      real :: xlon,ylat,xlonr,dzdx,dzdy
      real :: dzdx1,dzdx2,dzdy1,dzdy2
      real :: uuaux,vvaux,uupolaux,vvpolaux,ddpol,ffpol,wdummy
      real(kind=4) :: uuh(0:nxmax-1,0:nymax-1,nuvzmax)
      real(kind=4) :: vvh(0:nxmax-1,0:nymax-1,nuvzmax)
      real(kind=4) :: wwh(0:nxmax-1,0:nymax-1,nwzmax)
   real(kind=4) :: divh(0:nxmax-1,0:nymax-1,nuvzmax)

!     real :: uuh(0:nxmax-1,0:nymax-1,nuvzmax)
!     real :: divh(0:nxmax-1,0:nymax-1,nuvzmax)
      real :: div(0:nxmax-1,0:nymax-1,nuvzmax)
!     real :: vvh(0:nxmax-1,0:nymax-1,nuvzmax)
      real :: pvh(0:nxmax-1,0:nymax-1,nuvzmax)
!     real :: wwh(0:nxmax-1,0:nymax-1,nwzmax)
      real :: wzlev(nwzmax),uvwzlev(0:nxmax-1,0:nymax-1,nzmax)
      real :: wwh_svaa(nwzmax), wtc_stat(4,nzmax),u,v
      real,parameter :: const=r_air/ga
! CDA cloud commented
!      integer :: rain_cloud_above,kz_inv

      integer :: kz_inv
      real :: f_qvsat,pressure
! CDA some new declarations and mods
!      real :: rh,lsp,convp
      real :: rh,lsp,convp,prec,rhmin
      real,parameter :: precmin = 0.002



      logical :: init = .true.
! CDA
      logical :: lconvectprec = .true.



! set method_w_terrain_correction  & method_z_compute
      method_w_terrain_correction = 20
      method_z_compute = 10
       aa=0
      do iz = 1, nz
      do ix = 1, 4
          wtc_stat(ix,iz) = 0.0
      end do
      end do


!*************************************************************************
! If verttransform is called the first time, initialize heights of the   *
! z levels in meter. The heights are the heights of model levels, where  *
! u,v,T and qv are given, and of the interfaces, where w is given. So,   *
! the vertical resolution in the z system is doubled. As reference point,*
! the lower left corner of the grid is used.                             *
! Unlike in the eta system, no difference between heights for u,v and    *
! heights for w exists.                                                  *
!*************************************************************************

      if (init) then

! Search for a point with high surface pressure (i.e. not above significant topography)
! Then, use this point to construct a reference z profile, to be used at all times
!
! FLEXPART_WRF - use grid point with highest surface pressure
!**************************************************************************************

        pint = -1.0
        ixm = -999888777
        jym = -999888777
        do jy=0,nymin1
          do ix=0,nxmin1
!           if (ps(ix,jy,1,n).gt.100000.) then
            if (ps(ix,jy,1,n).gt.pint) then
              pint = ps(ix,jy,1,n)
              ixm=ix
              jym=jy
!             goto 3
            endif
         enddo
         enddo
3       continue
!       write(*,'(/a,2i4,1pe11.2)') 
!    &          'verttransform -- ixm,jym,ps() =', ixm, jym, pint


        tvold=tt2(ixm,jym,1,n)*(1.+0.378*ew(td2(ixm,jym,1,n))/ &
        ps(ixm,jym,1,n))
        pold=ps(ixm,jym,1,n)
        height(1)=0.

        do kz=2,nuvz
! use pressure from wrf met file
!         pint=akz(kz)+bkz(kz)*ps(ixm,jym,1,n)
          pint=pph(ixm,jym,kz,n)
          tv=tth(ixm,jym,kz,n)*(1.+0.608*qvh(ixm,jym,kz,n))


! NOTE: In FLEXPART versions up to 4.0, the number of model levels was doubled
! upon the transformation to z levels. In order to save computer memory, this is
! not done anymore in the standard version. However, this option can still be
! switched on by replacing the following lines with those below, that are
! currently commented out.
! Note that two more changes are necessary in this subroutine below.
! One change is also necessary in gridcheck.f, and another one in verttransform_nests.
!*************************************************************************************

          if (abs(tv-tvold).gt.0.2) then
            height(kz)= &
            height(kz-1)+const*log(pold/pint)* &
            (tv-tvold)/log(tv/tvold)
          else
            height(kz)=height(kz-1)+ &
            const*log(pold/pint)*tv
          endif

! 
! *** NOTE -- the doubled vertical resolution has not been tested in FLEXPART_WRF
! 
! Switch on following lines to use doubled vertical resolution
!*************************************************************
!         if (abs(tv-tvold).gt.0.2) then
!           height((kz-1)*2)=
!    +      height(max((kz-2)*2,1))+const*log(pold/pint)*
!    +      (tv-tvold)/log(tv/tvold)
!         else
!           height((kz-1)*2)=height(max((kz-2)*2,1))+
!    +      const*log(pold/pint)*tv
!         endif
! End doubled vertical resolution
 
! FLEXPART_WRF - get height from zzh
          if (method_z_compute .eq. 10) then
             if ((add_sfc_level .eq. 1) .and. (kz .eq. 2)) then
                height(kz) = 0.5*(zzh(ixm,jym,   3,n)+zzh(ixm,jym, 1,n)) &
                           - zzh(ixm,jym,1,n)
             else
                height(kz) = 0.5*(zzh(ixm,jym,kz+1,n)+zzh(ixm,jym,kz,n)) &
                           - zzh(ixm,jym,1,n)
             end if
          end if

          tvold=tv
          pold=pint
         enddo
           do kz=1,nz-1
         heightmid(kz)=0.5*(height(kz)+height(kz+1))
          enddo
! 
! *** NOTE -- the doubled vertical resolution has not been tested in FLEXPART_WRF
! 
! Switch on following lines to use doubled vertical resolution
!*************************************************************
!       do 7 kz=3,nz-1,2
!         height(kz)=0.5*(height(kz-1)+height(kz+1))
!       height(nz)=height(nz-1)+height(nz-1)-height(nz-2)
! End doubled vertical resolution


! Determine highest levels that can be within PBL
!************************************************

    do kz=1,nz
      if (height(kz).gt.hmixmax) then
        nmixz=kz
        goto 9
      endif
    end do
9   continue

! Do not repeat initialization of the Cartesian z grid
!*****************************************************

        init=.false.

      endif


! Loop over the whole grid
!*************************

!$OMP PARALLEL &
!$OMP PRIVATE(ix,jy,ixm,jym,tvold,pold,pint,tv,rhoh,uvzlev,wzlev, &
!$OMP uvwzlev,pinmconv,kz,iz,kmin,dz1,dz2,dz,ix1,jy1,ixp,jyp, &
!$OMP dzdy,dzdx,aa,u,v,wwh_svaa )
!$OMP DO SCHEDULE(STATIC,5)  
!!!!$OMP COLLAPSE(2) &
      do jy=0,nymin1
        do ix=0,nxmin1
          tvold=tt2(ix,jy,1,n)*(1.+0.378*ew(td2(ix,jy,1,n))/ &
                                         ps(ix,jy,1,n))
          pold=ps(ix,jy,1,n)
          uvzlev(1)=0.
          wzlev(1)=0.
          rhoh(1)=pold/(r_air*tvold)


! Compute heights of eta levels
!******************************

          do kz=2,nuvz
! use pressure from wrf met file
!           pint=akz(kz)+bkz(kz)*ps(ix,jy,1,n)
            pint=pph(ix,jy,kz,n)
            tv=tth(ix,jy,kz,n)*(1.+0.608*qvh(ix,jy,kz,n))
            rhoh(kz)=pint/(r_air*tv)

            if (abs(tv-tvold).gt.0.2) then
              uvzlev(kz)=uvzlev(kz-1)+const*log(pold/pint)* &
              (tv-tvold)/log(tv/tvold)
            else
              uvzlev(kz)=uvzlev(kz-1)+const*log(pold/pint)*tv
            endif
           
            tvold=tv
          pold=pint
      end do

          do kz=2,nwz-1
          wzlev(kz)=(uvzlev(kz+1)+uvzlev(kz))/2.
      end do

          wzlev(nwz)=wzlev(nwz-1)+ &
          uvzlev(nuvz)-uvzlev(nuvz-1)

! FLEXPART_WRF - get uvzlev & wzlev from zzh
          if (method_z_compute .eq. 10) then
            do kz = 2, nuvz
              if ((add_sfc_level .eq. 1) .and. (kz .eq. 2)) then
                uvzlev(kz) = 0.5*(zzh(ix,jy,   3,n) + zzh(ix,jy, 1,n)) &
                           - zzh(ix,jy,1,n)
              else
                uvzlev(kz) = 0.5*(zzh(ix,jy,kz+1,n) + zzh(ix,jy,kz,n)) &
                           - zzh(ix,jy,1,n)
              end if
            end do
            do kz = 2, nwz
              wzlev(kz) = zzh(ix,jy,kz+add_sfc_level,n)  &
                        - zzh(ix,jy,1,n)
            end do
          end if

          uvwzlev(ix,jy,1)=0.0
          do kz=2,nuvz
          uvwzlev(ix,jy,kz)=uvzlev(kz)
      end do

!     if ((ix .eq. ixm) .and. (jy .eq. jym)) then
!        write(*,'(/a)') 
!    &     'kz, height, uvzlev, wzlev, zzh-zzh(1) at ixm,jym  (in km)'
!        write(*,'(i3,4f8.3)') (kz, height(kz)*1.0e-3, 
!    &     uvzlev(kz)*1.0e-3, wzlev(kz)*1.0e-3,
!    &     (zzh(ix,jy,kz,n)-zzh(ix,jy,1,n))*1.0e-3, kz=nz,1,-1)
!        ixm = -9
!     end if

! Switch on following lines to use doubled vertical resolution
! Switch off the three lines above.
!*************************************************************
!22          uvwzlev(ix,jy,(kz-1)*2)=uvzlev(kz)
!          do 23 kz=2,nwz
!23          uvwzlev(ix,jy,(kz-1)*2+1)=wzlev(kz)
! End doubled vertical resolution

! pinmconv=(h2-h1)/(p2-p1)
!
! in flexpart_ecmwf, pinmconv is used to convert etadot to w
! in FLEXPART_WRF, vertical velocity is already m/s, so pinmconv=1.0

          if (wind_option.le.0) then
          pinmconv(1)=1.0
          do kz=2,nz-1
             pinmconv(kz)=1.0
          enddo
          pinmconv(nz)=1.0
          elseif (wind_option.ge.1) then

!         pinmconv(1)=(uvzlev(1+add_sfc_level)-0.) & 
!         /(eta_u_wrf(1)-1.)
          pinmconv(1)=(wzlev(2)-0.) &
          /(eta_w_wrf(2)-1.)
          do kz=2,nz-1

!          pinmconv(kz)=(uvzlev(kz+add_sfc_level)-uvzlev(kz-1+add_sfc_level)) & 
!          /(eta_u_wrf(kz)-eta_u_wrf(kz-1))
!          /(pph(ix,jy,kz+add_sfc_level,n)-pph(ix,jy,kz-1+add_sfc_level,n)) &
!          *(pph(ix,jy,1,n)-pph(ix,jy,nz,n))
!          *(ps(ix,jy,1,n)-p_top_wrf)

          pinmconv(kz)=(wzlev(kz+1)-wzlev(kz-1)) &
          /(eta_w_wrf(kz+1)-eta_w_wrf(kz-1))
           enddo   

          pinmconv(nwz)=pinmconv(nwz-1) ! 
          endif
! Levels, where u,v,t and q are given
!************************************

          uu(ix,jy,1,n)=uuh(ix,jy,1)
          vv(ix,jy,1,n)=vvh(ix,jy,1)
          div(ix,jy,1)=divh(ix,jy,1)
          tt(ix,jy,1,n)=tth(ix,jy,1,n)
          qv(ix,jy,1,n)=qvh(ix,jy,1,n)
          pv(ix,jy,1,n)=pvh(ix,jy,1)
          rho(ix,jy,1,n)=rhoh(1)
          uu(ix,jy,nz,n)=uuh(ix,jy,nuvz)
          vv(ix,jy,nz,n)=vvh(ix,jy,nuvz)
          tt(ix,jy,nz,n)=tth(ix,jy,nuvz,n)
          qv(ix,jy,nz,n)=qvh(ix,jy,nuvz,n)
          pv(ix,jy,nz,n)=pvh(ix,jy,nuvz)
          rho(ix,jy,nz,n)=rhoh(nuvz)
          tke(ix,jy,1,n)=tkeh(ix,jy,1,n)
          tke(ix,jy,nz,n)=tkeh(ix,jy,nuvz,n)
          ptt(ix,jy,1,n)=ptth(ix,jy,1,n)
          ptt(ix,jy,nz,n)=ptth(ix,jy,nuvz,n)


           kmin=2
          do iz=2,nz-1
            do kz=kmin,nuvz
              if(heightmid(iz).gt.uvzlev(nuvz)) then
               div(ix,jy,iz)=div(ix,jy,nz)
                goto 230
              endif
              if ((heightmid(iz).gt.uvzlev(kz-1)).and. &
                  (heightmid(iz).le.uvzlev(kz))) then
               dz1=heightmid(iz)-uvzlev(kz-1)
               dz2=uvzlev(kz)-heightmid(iz)
               dz=dz1+dz2
          div(ix,jy,iz)=(divh(ix,jy,kz-1)*dz2+divh(ix,jy,kz)*dz1)/dz
               kmin=kz
           goto 230
          endif
        end do
230      continue
      end do

           kmin=2
          do iz=2,nz-1
            do kz=kmin,nuvz
              if(height(iz).gt.uvzlev(nuvz)) then
                uu(ix,jy,iz,n)=uu(ix,jy,nz,n)
                vv(ix,jy,iz,n)=vv(ix,jy,nz,n)
                tt(ix,jy,iz,n)=tt(ix,jy,nz,n)
                qv(ix,jy,iz,n)=qv(ix,jy,nz,n)
                pv(ix,jy,iz,n)=pv(ix,jy,nz,n)
                rho(ix,jy,iz,n)=rho(ix,jy,nz,n)
                tke(ix,jy,iz,n)=tke(ix,jy,nz,n)
                ptt(ix,jy,iz,n)=ptt(ix,jy,nz,n)

                goto 30
              endif
              if ((height(iz).gt.uvzlev(kz-1)).and. &
                  (height(iz).le.uvzlev(kz))) then
               dz1=height(iz)-uvzlev(kz-1)
               dz2=uvzlev(kz)-height(iz)
               dz=dz1+dz2
               uu(ix,jy,iz,n)=(uuh(ix,jy,kz-1)*dz2+uuh(ix,jy,kz)*dz1)/dz
               vv(ix,jy,iz,n)=(vvh(ix,jy,kz-1)*dz2+vvh(ix,jy,kz)*dz1)/dz
               tt(ix,jy,iz,n)=(tth(ix,jy,kz-1,n)*dz2 &
                    +tth(ix,jy,kz,n)*dz1)/dz
               qv(ix,jy,iz,n)=(qvh(ix,jy,kz-1,n)*dz2 &
                    +qvh(ix,jy,kz,n)*dz1)/dz
               pv(ix,jy,iz,n)=(pvh(ix,jy,kz-1)*dz2+pvh(ix,jy,kz)*dz1)/dz
               rho(ix,jy,iz,n)=(rhoh(kz-1)*dz2+rhoh(kz)*dz1)/dz
               tke(ix,jy,iz,n)=(tkeh(ix,jy,kz-1,n)*dz2 &
                    +tkeh(ix,jy,kz,n)*dz1)/dz
               ptt(ix,jy,iz,n)=(ptth(ix,jy,kz-1,n)*dz2 &
                    +ptth(ix,jy,kz,n)*dz1)/dz


               kmin=kz
           goto 30
          endif
        end do
30      continue
      end do


! Levels, where w is given
!*************************

!          ww(ix,jy,1,n)=wwh(ix,jy,1)*pinmconv(1)
!          ww(ix,jy,nz,n)=wwh(ix,jy,nwz)*pinmconv(nz)
!          kmin=2
!          do iz=2,nz
!            do  kz=kmin,nwz
!              if ((height(iz).gt.wzlev(kz-1)).and. &
!                  (height(iz).le.wzlev(kz))) then
!               dz1=height(iz)-wzlev(kz-1)
!               dz2=wzlev(kz)-height(iz)
!               dz=dz1+dz2
!!              ww(ix,jy,iz,n)=(wwh(ix,jy,kz-1)*dz2*pinmconv(kz-1)
!!    +  +wwh(ix,jy,kz)*dz1*pinmconv(kz))/dz
!               kmin=kz
!               goto 40
!              endif
!        end do
!40      continue
!      end do

          if (method_w_terrain_correction .eq. 20) then
! apply w correction assuming that the WRF w is "absolute w";
! apply it here to wwh; set wwh=0 at iz=1
!            do iz = 1, nz
!               wtc_stat(1,iz) = wtc_stat(1,iz) + ww(ix,jy,iz,n)
!               wtc_stat(2,iz) = wtc_stat(2,iz) + abs(ww(ix,jy,iz,n))
!            end do

!             if ((ix.eq.0) .and. (jy.eq.0)) write(*,*) 
!     &            'verttransform doing method_w_terrain_correction =', 
!     &            method_w_terrain_correction
             ix1 = max( ix-1, 0 )
             jy1 = max( jy-1, 0 )
             ixp = min( ix+1, nx-1 )
             jyp = min( jy+1, ny-1 )
          if (wind_option.eq.0) then
            dzdx=(oro(ixp,jy) - oro(ix1,jy))/(dx*(ixp-ix1)*m_x(ix,jy,1))
            dzdy=(oro(ix,jyp) - oro(ix,jy1))/(dy*(jyp-jy1)*m_y(ix,jy,1))
             do kz = 1, nwz-1
                wwh_svaa(kz) = wwh(ix,jy,kz)
                wwh(ix,jy,kz) = wwh(ix,jy,kz)*pinmconv(kz) &
!               wwh(ix,jy,kz) = 
              - (uuh(ix,jy,kz)*dzdx + vvh(ix,jy,kz)*dzdy)   !this is correct. term of variation of geopot not necessary

                if (kz .eq. 1) wwh(ix,jy,kz) = 0.0
              aa=aa+1
             end do
          elseif (wind_option.ge.1) then
             do kz = 2, nwz-1
                wwh_svaa(kz) = wwh(ix,jy,kz)
!             dzdx=(zzh(ixp,jy,kz,n) - zzh(ix1,jy,kz,n))
!     +  /(dx*(ixp-ix1))
!             dzdy=(zzh(ix,jyp,kz,n) - zzh(ix,jy1,kz,n))
!     +  /(dy*(jyp-jy1))
!             dzdx=(zzh(ixp,jy,kz,n) - zzh(ix1,jy,kz,n)-zzh(ixp,jy,1,n) &
!        +zzh(ix1,jy,1,n)) &
!        /(dx*(ixp-ix1)*m_x(ix,jy,1))
!             dzdy=(zzh(ix,jyp,kz,n) - zzh(ix,jy1,kz,n)-zzh(ix,jyp,1,n) &
!        +zzh(ix,jy1,1,n)) &
!        /(dy*(jyp-jy1)*m_y(ix,jy,1))

        dzdx=(zzh(ixp,jy,kz+add_sfc_level,n) - zzh(ix1,jy,kz+add_sfc_level,n) & 
             -zzh(ixp,jy,1,n)+zzh(ix1,jy,1,n))/(dx*(ixp-ix1)*m_x(ix,jy,1))
        dzdy=(zzh(ix,jyp,kz+add_sfc_level,n) - zzh(ix,jy1,kz+add_sfc_level,n)  &
             -zzh(ix,jyp,1,n)+zzh(ix,jy1,1,n))/(dy*(jyp-jy1)*m_y(ix,jy,1))
           u=0.5*(uuh(ix,jy,kz+add_sfc_level)+uuh(ix,jy,kz-1+add_sfc_level))
           v=0.5*(vvh(ix,jy,kz+add_sfc_level)+vvh(ix,jy,kz-1+add_sfc_level))
                wwh(ix,jy,kz) = wwh(ix,jy,kz)*pinmconv(kz)  &
            + (u*dzdx + v*dzdy) ! variation of geopot on sigma is necessary

!                wwh(ix,jy,kz) = wwh(ix,jy,kz)*pinmconv(kz) &
!           + (uuh(ix,jy,kz)*dzdx + vvh(ix,jy,kz)*dzdy) ! variation of geopot on sigma is necessary
!            if (kz .eq. 1) wwh(ix,jy,kz) = 0.0
             if (kz .eq. 1) wwh(ix,jy,kz) = wwh(ix,jy,kz)*pinmconv(kz)
!             aa=aa+1
             end do
         endif
         if (wind_option.eq.-1) then
!        ww(ix,jy,1,n)=wwh(ix,jy,1)
         ww(ix,jy,1,n)=0.
         do iz=2,nz
         ww(ix,jy,iz,n)=ww(ix,jy,iz-1,n)-(height(iz)-height(iz-1))* &
        div(ix,jy,iz-1)
         enddo
         else

             ww(ix,jy,1,n)=wwh(ix,jy,1)
             ww(ix,jy,nz,n)=wwh(ix,jy,nwz)
             kmin=2
             do iz=2,nz
               do kz=kmin,nwz
                 if ((height(iz).gt.wzlev(kz-1)).and. &
                     (height(iz).le.wzlev(kz))) then
                  dz1=height(iz)-wzlev(kz-1)
                  dz2=wzlev(kz)-height(iz)
                  dz=dz1+dz2
                  ww(ix,jy,iz,n)=(wwh(ix,jy,kz-1)*dz2 &
        +wwh(ix,jy,kz)*dz1) &
                    /dz
                  kmin=kz
           goto 4000
          endif
        end do
4000      continue
      end do
           endif

!             do kz = 1, nwz
!                wwh(ix,jy,kz) = wwh_svaa(kz)
!             end do

!            do iz = 1, nz
!               wtc_stat(3,iz) = wtc_stat(3,iz) + ww(ix,jy,iz,n)
!               wtc_stat(4,iz) = wtc_stat(4,iz) + abs(ww(ix,jy,iz,n))
!            end do
          end if

! Compute density gradients at intermediate levels
!*************************************************

          drhodz(ix,jy,1,n)=(rho(ix,jy,2,n)-rho(ix,jy,1,n))/ &
            (height(2)-height(1))
          do kz=2,nz-1
          drhodz(ix,jy,kz,n)=(rho(ix,jy,kz+1,n)-rho(ix,jy,kz-1,n))/ &
            (height(kz+1)-height(kz-1))
      end do

          drhodz(ix,jy,nz,n)=drhodz(ix,jy,nz-1,n)

    end do
  end do
!$OMP END DO
!$OMP END PARALLEL

!****************************************************************
! Compute slope of eta levels in windward direction and resulting
! vertical wind correction
!
! The ECMWF model uses a hybrid-pressure vertical coordinate, "eta"
!    The "eta" coordinate transitions from terrain-following near
!    the surface to constant pressure in the stratosphere.
!    The vertical velocities in the ECMWF grib files are "eta_dot"
! FLEXPART uses a "height above ground" vertical coordinate
!    which we will call "hag".
!    The vertical velocity is uses (in ww array) is "hag_dot".
! Converting from eta_dot to hag_dot involves
!    >> multiplying by pinmconv = [d(hag)/d(eta)]
!    >> adding a term that accounts for the fact that
!       "eta" varies on constant "hag" surfaces.
!       This term is [u*d(hag)/dx + v*d(hag)/dy], with the
!       partial derivatives taken with "eta" being constant
!
! The WRF model uses a similar (to ECMWF) vertical coordinate.
!    HOWEVER, the vertical velocities in the WRF output files
!    are the "true/absolute w" in m/s.  (Is this true?)
! Converting from "absolute w" to hag_dot involves
!    adding a term that accounts for the fact that
!    "absolute z" varies on constant "hag" surfaces.
!    This term is [- u*d(oro)/dx - v*d(oro)/dy]
!
! The FLEXPART code did not apply the terrain corrections
!    at jy=0 & ny-1; ix=0 & nx-1; iz=1 & nz.
! FLEXPART_WRF applies the correction at all grid points
!****************************************************************


! If north pole is in the domain, calculate wind velocities in polar
! stereographic coordinates
!*******************************************************************
 
      if (nglobal) then
        write(*,*)
        write(*,*) '*** stopping in verttransform ***'
        write(*,*) '    the nglobal code section should not be active'
        write(*,*)
        stop
!        do 74 jy=int(switchnorthg)-2,nymin1
!          ylat=ylat0+real(jy)*dy
!          do 74 ix=0,nxmin1
!            xlon=xlon0+real(ix)*dx
!            do 74 iz=1,nz
!74            call cc2gll(northpolemap,ylat,xlon,uu(ix,jy,iz,n),
!     +        vv(ix,jy,iz,n),uupol(ix,jy,iz,n),
!     +        vvpol(ix,jy,iz,n))
!
!
!        do 76 iz=1,nz
!
!* CALCULATE FFPOL, DDPOL FOR CENTRAL GRID POINT
!          xlon=xlon0+real(nx/2-1)*dx
!          xlonr=xlon*pi/180.
!          ffpol=sqrt(uu(nx/2-1,nymin1,iz,n)**2+
!     &               vv(nx/2-1,nymin1,iz,n)**2)
!          if(vv(nx/2-1,nymin1,iz,n).lt.0.) then
!            ddpol=atan(uu(nx/2-1,nymin1,iz,n)/
!     &                 vv(nx/2-1,nymin1,iz,n))-xlonr
!          else
!            ddpol=pi+atan(uu(nx/2-1,nymin1,iz,n)/
!     &                    vv(nx/2-1,nymin1,iz,n))-xlonr
!          endif
!          if(ddpol.lt.0.) ddpol=2.0*pi+ddpol
!          if(ddpol.gt.2.0*pi) ddpol=ddpol-2.0*pi
! 
!* CALCULATE U,V FOR 180 DEG, TRANSFORM TO POLAR STEREOGRAPHIC GRID
!          xlon=180.0
!          xlonr=xlon*pi/180.
!          ylat=90.0
!          uuaux=-ffpol*sin(xlonr+ddpol)
!          vvaux=-ffpol*cos(xlonr+ddpol)
!          call cc2gll(northpolemap,ylat,xlon,uuaux,vvaux,uupolaux,
!     +      vvpolaux)
!      
!          jy=nymin1
!          do 76 ix=0,nxmin1
!            uupol(ix,jy,iz,n)=uupolaux
!            vvpol(ix,jy,iz,n)=vvpolaux
!76      continue
! 
! 
!* Fix: Set W at pole to the zonally averaged W of the next equator-
!* ward parallel of latitude
! 
!      do 85 iz=1,nz
!          wdummy=0.
!          jy=ny-2
!          do 80 ix=0,nxmin1
!80          wdummy=wdummy+ww(ix,jy,iz,n)
!          wdummy=wdummy/real(nx)
!          jy=nymin1
!          do 85 ix=0,nxmin1
!85          ww(ix,jy,iz,n)=wdummy
 
      endif 

 
! If south pole is in the domain, calculate wind velocities in polar
! stereographic coordinates
!*******************************************************************
 
      if (sglobal) then
        write(*,*)
        write(*,*) '*** stopping in verttransform ***'
        write(*,*) '    the sglobal code section should not be active'
        write(*,*)
        stop
!        do 77 jy=0,int(switchsouthg)+3
!          ylat=ylat0+real(jy)*dy
!          do 77 ix=0,nxmin1
!            xlon=xlon0+real(ix)*dx
!            do 77 iz=1,nz
!77            call cc2gll(southpolemap,ylat,xlon,uu(ix,jy,iz,n),
!     +        vv(ix,jy,iz,n),uupol(ix,jy,iz,n),
!     +        vvpol(ix,jy,iz,n))
!      
!        do 79 iz=1,nz
! 
!* CALCULATE FFPOL, DDPOL FOR CENTRAL GRID POINT
!          xlon=xlon0+real(nx/2-1)*dx
!          xlonr=xlon*pi/180.
!          ffpol=sqrt(uu(nx/2-1,0,iz,n)**2+
!     &               vv(nx/2-1,0,iz,n)**2)
!          if(vv(nx/2-1,0,iz,n).lt.0.) then
!            ddpol=atan(uu(nx/2-1,0,iz,n)/
!     &                 vv(nx/2-1,0,iz,n))+xlonr
!          else
!            ddpol=pi+atan(uu(nx/2-1,0,iz,n)/
!     &                    vv(nx/2-1,0,iz,n))+xlonr
!          endif
!          if(ddpol.lt.0.) ddpol=2.0*pi+ddpol
!          if(ddpol.gt.2.0*pi) ddpol=ddpol-2.0*pi
! 
!* CALCULATE U,V FOR 180 DEG, TRANSFORM TO POLAR STEREOGRAPHIC GRID
!          xlon=180.0
!          xlonr=xlon*pi/180.
!          ylat=-90.0
!          uuaux=+ffpol*sin(xlonr-ddpol)
!          vvaux=-ffpol*cos(xlonr-ddpol)
!          call cc2gll(northpolemap,ylat,xlon,uuaux,vvaux,uupolaux,
!     +      vvpolaux)
!      
!          jy=0
!          do 79 ix=0,nxmin1
!            uupol(ix,jy,iz,n)=uupolaux
!79          vvpol(ix,jy,iz,n)=vvpolaux
! 
! 
!* Fix: Set W at pole to the zonally averaged W of the next equator-
!* ward parallel of latitude
! 
!        do 95 iz=1,nz
!          wdummy=0.
!          jy=1
!          do 90 ix=0,nxmin1
!90          wdummy=wdummy+ww(ix,jy,iz,n)
!          wdummy=wdummy/real(nx)
!          jy=0
!          do 95 ix=0,nxmin1
!95          ww(ix,jy,iz,n)=wdummy
      endif

  !write (*,*) 'initializing clouds, n:',n,nymin1,nxmin1,nz^M
  !   create a cloud and rainout/washout field, clouds occur where rh>80%^M
  !   total cloudheight is stored at level 0^M
    do 100 jy=0,nymin1
      do 100 ix=0,nxmin1
!        rain_cloud_above=0
        lsp=lsprec(ix,jy,1,n)  
        convp=convprec(ix,jy,1,n)
!        cloudsh(ix,jy,n)=0

          prec=lsp+convp
          if (lsp.gt.convp) then !  prectype='lsp'
            lconvectprec = .false.
          else ! prectype='cp '
            lconvectprec = .true.
          endif
          rhmin = 0.90 ! standard condition for presence of clouds

!CPS       note that original by Sabine Eckhart was 80%
!CPS       however, for T<-20 C we consider saturation over ice
!CPS       so I think 90% should be enough


          icloudbot(ix,jy,n)=icmv
          icloudtop=icmv ! this is just a local variable
98        do kz=1,nz
            pressure=rho(ix,jy,kz,n)*r_air*tt(ix,jy,kz,n)
            rh=qv(ix,jy,kz,n)/f_qvsat(pressure,tt(ix,jy,kz,n))
!cps            if (prec.gt.0.01) print*,'relhum',prec,kz,rh,height(kz)
            if (rh .gt. rhmin) then
              if (icloudbot(ix,jy,n) .eq. icmv) then
                icloudbot(ix,jy,n)=nint(height(kz))
              endif
              icloudtop=nint(height(kz)) ! use int to save memory
            endif
          enddo


!CPS try to get a cloud thicker than 50 m
!CPS if there is at least .01 mm/h  - changed to 0.002 and put into
!CPS parameter precpmin
          if ((icloudbot(ix,jy,n) .eq. icmv .or. &
               icloudtop-icloudbot(ix,jy,n) .lt. 50) .and. &
               prec .gt. precmin) then
            rhmin = rhmin - 0.05
            if (rhmin .ge. 0.30) goto 98 ! give up for <= 25% rel.hum.
          endif
!CPS implement a rough fix for badly represented convection
!CPS is based on looking at a limited set of comparison data
          if (lconvectprec .and. icloudtop .lt. 6000 .and. &
              prec .gt. precmin) then
            if (convp .lt. 0.1) then
              icloudbot(ix,jy,n) = 500
              icloudtop =         8000
            else
              icloudbot(ix,jy,n) = 0
              icloudtop =      10000
            endif
          endif
          if (icloudtop .ne. icmv) then
            icloudthck(ix,jy,n) = icloudtop-icloudbot(ix,jy,n)
          else
            icloudthck(ix,jy,n) = icmv
          endif
!CPS  get rid of too thin clouds
          if (icloudthck(ix,jy,n) .lt. 50) then
            icloudbot(ix,jy,n)=icmv
            icloudthck(ix,jy,n)=icmv
          endif

100   continue






!       do kz_inv=1,nz-1
!           kz=nz-kz_inv+1
!           pressure=rho(ix,jy,kz,n)*r_air*tt(ix,jy,kz,n)
!           rh=qv(ix,jy,kz,n)/f_qvsat(pressure,tt(ix,jy,kz,n))
!           clouds(ix,jy,kz,n)=0
!           if (rh.gt.0.8) then ! in cloud
!              if ((lsp.gt.0.01).or.(convp.gt.0.01)) then ! cloud and precipitation
!                 rain_cloud_above=1
!                 cloudsh(ix,jy,n)=cloudsh(ix,jy,n)+ &
!                      height(kz)-height(kz-1)
!                 if (lsp.ge.convp) then
!                    clouds(ix,jy,kz,n)=3 ! lsp dominated rainout
!                 else
!                    clouds(ix,jy,kz,n)=2 ! convp dominated rainout
!                 endif 
!              else ! no precipitation
!                    clouds(ix,jy,kz,n)=1 ! cloud
!              endif
!           else ! no cloud
!              if (rain_cloud_above.eq.1) then ! scavenging
!                 if (lsp.ge.convp) then
!                    clouds(ix,jy,kz,n)=5 ! lsp dominated washout
!                 else
!                    clouds(ix,jy,kz,n)=4 ! convp dominated washout
!                 endif
!              endif
!           endif
!        end do
!      end do
!    end do


end subroutine verttransform
