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
      subroutine verttransform_nests(n,uuhn,vvhn,wwhn,pvhn,divhn)
!                                    i   i    i    i   i
!*******************************************************************************
!                                                                              *
! Note:  This is the FLEXPART_WRF version of subroutine verttransform_nests.   *
!     The computational grid is the WRF x-y grid rather than lat-lon.          *
!                                                                              *
!     This subroutine transforms temperature, dew point temperature and        *
!     wind components from eta to meter coordinates.                           *
!     The vertical wind component is transformed from Pa/s to m/s using        *
!     the conversion factor pinmconv.                                          *
!     In addition, this routine calculates vertical density gradients          *
!     needed for the parameterization of the turbulent velocities.             *
!     It is similar to verttransform, but makes the transformations for        *
!     the nested grids.                                                        *
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
!     Changes, Bernd C. Krueger, Feb. 2001:       (marked "C-cv")              *
!        Variables tthn and qvhn (on eta coordinates) from common block        *
!                                                                              *
!     16 Nov 2005, R. Easter - changes for FLEXPART_WRF                        *
!     17 Nov 2005 - R. Easter - terrain correction applied to ww.  There are   *
!            now 3 options, controlled by "method_w_terrain_correction"        *
!                                                                              *
!     11 June 2007  --  convert TKEhn to tken
!     25 June 2007  --  convert ptthn to pttn
!     Jan 2012, J Brioude:  modified to handle different wind options and openmp 
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
! nxn,nyn,nuvz,nwz                field dimensions in x,y and z direction      *
! uun                             wind components in x-direction [m/s]         *
! vvn                             wind components in y-direction [m/s]         *
! wwn                             wind components in z-direction [deltaeta/s]  *
! ttn                             temperature [K]                              *
! pvn                             potential vorticity (pvu)                    *
! psn                             surface pressure [Pa]                        *
!                                                                              *
!*******************************************************************************

  use par_mod
  use com_mod

!      include 'includepar'
!      include 'includecom'
    implicit none
     integer :: ix,jy,kz,iz,n,l,kmin,kl,klp,ix1,jy1,ixp,jyp,ixm,jym
     integer :: icloudtop
     real :: rh,lsp,convp,prec,rhmin
     integer :: method_z_compute,aa,dimx,dimy
     real :: uvzlev(nuvzmax),wzlev(nwzmax),rhoh(nuvzmax),pinmconv(nzmax)
     real :: uvwzlev(0:nxmaxn-1,0:nymaxn-1,nzmax)
     real :: ew,pint,tv,tvold,pold,const,dz1,dz2,dz,ui,vi
     real :: dzdx,dzdy
     real :: dzdx1,dzdx2,dzdy1,dzdy2
     real :: pvhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
     real :: divn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
     real(kind=4) :: uuhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
     real(kind=4) :: vvhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
     real(kind=4) :: wwhn(0:nxmaxn-1,0:nymaxn-1,nwzmax,maxnests)
     real(kind=4) :: divhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
     real :: wwhn_svaa(nwzmax),u,v
     parameter(const=r_air/ga)
!      integer :: rain_cloud_above,kz_inv

      real :: f_qvsat,pressure
!      real :: rh,lsp,convp
      real,parameter :: precmin = 0.002

! CDA
      logical :: lconvectprec = .true.


! set method_z_compute
      method_z_compute = 10


! Loop over all nests
!********************

      do l=1,numbnests
        dimy=nyn(l)-1
        dimx=nxn(l)-1
!      print*,'start omp '
! Loop over the whole grid
!*************************
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ix,jy,ixm,jym,tvold,pold,pint,tv,rhoh,uvzlev,wzlev, &
!$OMP uvwzlev,pinmconv,kz,iz,kmin,dz1,dz2,dz,ix1,jy1,ixp,jyp, &
!$OMP dzdy,dzdx,aa,u,v )
!$OMP DO
      do jy=0,dimy
        do ix=0,dimx

           tvold=tt2n(ix,jy,1,n,l)*(1.+0.378*ew(td2n(ix,jy,1,n,l))/ &
                                         psn(ix,jy,1,n,l))
          pold=psn(ix,jy,1,n,l)
          uvzlev(1)=0.
          wzlev(1)=0.
          rhoh(1)=pold/(r_air*tvold)


! Compute heights of eta levels
!******************************

          do kz=2,nuvz
! FLEXPART_WRF - pphn hold pressure
!           pint=akz(kz)+bkz(kz)*psn(ix,jy,1,n,l)
            pint=pphn(ix,jy,kz,n,l)
            tv=tthn(ix,jy,kz,n,l)*(1.+0.608*qvhn(ix,jy,kz,n,l))
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


!      print*,'etap 1',ix,jy

          do kz=2,nwz-1
            wzlev(kz)=(uvzlev(kz+1)+uvzlev(kz))/2.
      end do
          wzlev(nwz)=wzlev(nwz-1)+ &
          uvzlev(nuvz)-uvzlev(nuvz-1)

! FLEXPART_WRF - get uvzlev & wzlev from zzh
          if (method_z_compute .eq. 10) then
            do kz = 2, nuvz
              if ((add_sfc_level .eq. 1) .and. (kz .eq. 2)) then
                uvzlev(kz) = 0.5*(zzhn(ix,jy,3,n,l) +  &
                                  zzhn(ix,jy,1,n,l)) &
                           - zzhn(ix,jy,1,n,l)
              else
                uvzlev(kz) = 0.5*(zzhn(ix,jy,kz+1,n,l) + &
                                  zzhn(ix,jy,kz  ,n,l)) &
                           - zzhn(ix,jy,1,n,l)
              end if
            end do
            do kz = 2, nwz
              wzlev(kz) = zzhn(ix,jy,kz+add_sfc_level,n,l)  &
                        - zzhn(ix,jy,1,n,l)
            end do
          end if

!      print*,'etap 2',ix,jy
! NOTE: In FLEXPART versions up to 4.0, the number of model levels was doubled
! upon the transformation to z levels. In order to save computer memory, this is
! not done anymore in the standard version. However, this option can still be
! switched on by replacing the following lines with those below, that are
! currently commented out.  
! Note that one change is also necessary in gridcheck.f,
! and three changes in verttransform.f
!
! *** NOTE -- the doubled vertical resolution has not been tested in FLEXPART_WRF
!*******************************************************************************
          uvwzlev(ix,jy,1)=0.0
          do kz=2,nuvz
          uvwzlev(ix,jy,kz)=uvzlev(kz)
         enddo
! Switch on following lines to use doubled vertical resolution
! Switch off the three lines above.
!
! *** NOTE -- the doubled vertical resolution has not been tested in FLEXPART_WRF
!*************************************************************
!22          uvwzlev(ix,jy,(kz-1)*2)=uvzlev(kz)
!          do 23 kz=2,nwz
!23          uvwzlev(ix,jy,(kz-1)*2+1)=wzlev(kz)
! End doubled vertical resolution

! pinmconv=(h2-h1)/(p2-p1)
!
! in flexpart_ecmwf, pinmconv is used to convert etadot to w
! in FLEXPART_WRF, vertical velocity is already m/s, so pinmconv=1.0
!
!         pinmconv(1)=(uvwzlev(ix,jy,2)-uvwzlev(ix,jy,1))/
!    +    ((aknew(2)+bknew(2)*psn(ix,jy,1,n,l))-
!    +    (aknew(1)+bknew(1)*psn(ix,jy,1,n,l)))
          if (wind_option.eq.0) then
          pinmconv(1)=1.0
          do kz=2,nz-1
             pinmconv(kz)=1.0
          enddo
          pinmconv(nz)=1.0
          elseif (wind_option.eq.1) then
!          pinmconv(1)=(uvzlev(2)-uvzlev(1)) &
!          /(eta_u_wrf(1)-1.)
          pinmconv(1)=(wzlev(2)-0.) &
          /(eta_w_wrf(2)-1.)
          do kz=2,nz-1
!          pinmconv(kz)=(uvzlev(kz)-uvzlev(kz-1)) &
!          /(eta_u_wrf(kz)-eta_u_wrf(kz-1))
          pinmconv(kz)=(wzlev(kz+1)-wzlev(kz-1)) &
          /(eta_w_wrf(kz+1)-eta_w_wrf(kz-1))
           enddo
          pinmconv(nwz)=pinmconv(nwz-1)
          endif


!      print*,'etap 3',ix,jy

! Levels, where u,v,t and q are given
!************************************

          uun(ix,jy,1,n,l)=uuhn(ix,jy,1,l)
          vvn(ix,jy,1,n,l)=vvhn(ix,jy,1,l)
          divn(ix,jy,1,l)=divhn(ix,jy,1,l)
          ttn(ix,jy,1,n,l)=tthn(ix,jy,1,n,l)
          qvn(ix,jy,1,n,l)=qvhn(ix,jy,1,n,l)
          pvn(ix,jy,1,n,l)=pvhn(ix,jy,1,l)
          rhon(ix,jy,1,n,l)=rhoh(1)
          uun(ix,jy,nz,n,l)=uuhn(ix,jy,nuvz,l)
          vvn(ix,jy,nz,n,l)=vvhn(ix,jy,nuvz,l)
          ttn(ix,jy,nz,n,l)=tthn(ix,jy,nuvz,n,l)
          qvn(ix,jy,nz,n,l)=qvhn(ix,jy,nuvz,n,l)
          pvn(ix,jy,nz,n,l)=pvhn(ix,jy,nuvz,l)
          rhon(ix,jy,nz,n,l)=rhoh(nuvz)
          tken(ix,jy,1,n,l)=tkehn(ix,jy,1,n,l)
          tken(ix,jy,nz,n,l)=tkehn(ix,jy,nuvz,n,l)
          pttn(ix,jy,1,n,l)=ptthn(ix,jy,1,n,l)
          pttn(ix,jy,nz,n,l)=ptthn(ix,jy,nuvz,n,l)


!      print*,'etap 3.5',ix,jy
           kmin=2
          do iz=2,nz-1
            do kz=kmin,nuvz
              if(heightmid(iz).gt.uvzlev(nuvz)) then
               divn(ix,jy,iz,l)=divn(ix,jy,nz,l)
                goto 230
              endif
!!      print*,'etap 3.7',kz,iz,heightmid(iz),uvzlev(kz-1),uvzlev(kz)
              if ((heightmid(iz).gt.uvzlev(kz-1)).and. &
                  (heightmid(iz).le.uvzlev(kz))) then
               dz1=heightmid(iz)-uvzlev(kz-1)
               dz2=uvzlev(kz)-heightmid(iz)
               dz=dz1+dz2
          divn(ix,jy,iz,l)=(divhn(ix,jy,kz-1,l)*dz2+divhn(ix,jy,kz,l)*dz1)/dz
               kmin=kz
           goto 230
          endif
        end do
230      continue
      end do

!      print*,'etap 4',ix,jy

          kmin=2
          do iz=2,nz-1
            do kz=kmin,nuvz
              if(height(iz).gt.uvzlev(nuvz)) then
                uun(ix,jy,iz,n,l)=uun(ix,jy,nz,n,l)
                vvn(ix,jy,iz,n,l)=vvn(ix,jy,nz,n,l)
                ttn(ix,jy,iz,n,l)=ttn(ix,jy,nz,n,l)
                qvn(ix,jy,iz,n,l)=qvn(ix,jy,nz,n,l)
                pvn(ix,jy,iz,n,l)=pvn(ix,jy,nz,n,l)
                rhon(ix,jy,iz,n,l)=rhon(ix,jy,nz,n,l)
                tken(ix,jy,iz,n,l)=tken(ix,jy,nz,n,l)
                pttn(ix,jy,iz,n,l)=pttn(ix,jy,nz,n,l)
                goto 30
              endif
              if ((height(iz).gt.uvzlev(kz-1)).and. &
                  (height(iz).le.uvzlev(kz))) then
               dz1=height(iz)-uvzlev(kz-1)
               dz2=uvzlev(kz)-height(iz)
               dz=dz1+dz2
               uun(ix,jy,iz,n,l)=(uuhn(ix,jy,kz-1,l)*dz2+ &
               uuhn(ix,jy,kz,l)*dz1)/dz
               vvn(ix,jy,iz,n,l)=(vvhn(ix,jy,kz-1,l)*dz2+ &
               vvhn(ix,jy,kz,l)*dz1)/dz
               ttn(ix,jy,iz,n,l)=(tthn(ix,jy,kz-1,n,l)*dz2+ &
               tthn(ix,jy,kz,n,l)*dz1)/dz
               qvn(ix,jy,iz,n,l)=(qvhn(ix,jy,kz-1,n,l)*dz2+ &
               qvhn(ix,jy,kz,n,l)*dz1)/dz
               pvn(ix,jy,iz,n,l)=(pvhn(ix,jy,kz-1,l)*dz2+ &
               pvhn(ix,jy,kz,l)*dz1)/dz
               rhon(ix,jy,iz,n,l)=(rhoh(kz-1)*dz2+rhoh(kz)*dz1)/dz
               tken(ix,jy,iz,n,l)=(tkehn(ix,jy,kz-1,n,l)*dz2+ &
               tkehn(ix,jy,kz,n,l)*dz1)/dz
               pttn(ix,jy,iz,n,l)=(ptthn(ix,jy,kz-1,n,l)*dz2+ &
               ptthn(ix,jy,kz,n,l)*dz1)/dz


               kmin=kz
           goto 30
          endif
        end do
30      continue
      end do

!      print*,'continue to ww in nests'
! Levels, where w is given
!*************************

          if (method_w_terrain_correction .eq. 20) then
! apply w correction assuming that the WRF w is "absolute w";
! apply it here to wwh; set wwh=0 at iz=1
             ix1 = max( ix-1, 0 )
             jy1 = max( jy-1, 0 )
             ixp = min( ix+1, nxn(l)-1 )
             jyp = min( jy+1, nyn(l)-1 )
          if (wind_option.eq.0) then
         dzdx=(oron(ixp,jy,l)-oron(ix1,jy,l))/(dxn(l)*(ixp-ix1)*m_xn(ix,jy,1,l))
         dzdy=(oron(ix,jyp,l)-oron(ix,jy1,l))/(dyn(l)*(jyp-jy1)*m_yn(ix,jy,1,l))

             do kz = 1, nwz-1
!               wwhn_svaa(kz) = wwhn(ix,jy,kz,l)
                wwhn(ix,jy,kz,l) = wwhn(ix,jy,kz,l) &
                     - (uuhn(ix,jy,kz,l)*dzdx + vvhn(ix,jy,kz,l)*dzdy)  
!               if (kz .eq. 1) wwhn(ix,jy,kz,l) = 0.0
             end do
          elseif (wind_option.ge.1) then
             do kz = 2, nwz-1
!               wwhn_svaa(kz) = wwhn(ix,jy,kz,l)
             dzdx=(zzhn(ixp,jy,kz+add_sfc_level,n,l) - zzhn(ix1,jy,kz+add_sfc_level,n,l) &
        -zzhn(ixp,jy,1,n,l) &
        +zzhn(ix1,jy,1,n,l)) &
        /(dxn(l)*(ixp-ix1)*m_xn(ix,jy,1,l))
             dzdy=(zzhn(ix,jyp,kz+add_sfc_level,n,l) - zzhn(ix,jy1,kz+add_sfc_level,n,l) &
        -zzhn(ix,jyp,1,n,l) &
        +zzhn(ix,jy1,1,n,l)) &
        /(dyn(l)*(jyp-jy1)*m_yn(ix,jy,1,l))

        dzdx=(zzhn(ixp,jy,kz+add_sfc_level,n,l) - zzhn(ix1,jy,kz+add_sfc_level,n,l) &
             -zzhn(ixp,jy,1,n,l)+zzhn(ix1,jy,1,n,l))/(dxn(l)*(ixp-ix1)*m_xn(ix,jy,1,l))
        dzdy=(zzhn(ix,jyp,kz+add_sfc_level,n,l) - zzhn(ix,jy1,kz+add_sfc_level,n,l)  &
             -zzhn(ix,jyp,1,n,l)+zzhn(ix,jy1,1,n,l))/(dyn(l)*(jyp-jy1)*m_yn(ix,jy,1,l))
           u=0.5*(uuhn(ix,jy,kz+add_sfc_level,l)+uuhn(ix,jy,kz-1+add_sfc_level,l))
           v=0.5*(vvhn(ix,jy,kz+add_sfc_level,l)+vvhn(ix,jy,kz-1+add_sfc_level,l))

                wwhn(ix,jy,kz,l) = wwhn(ix,jy,kz,l)*pinmconv(kz) &
!            + (uuhn(ix,jy,kz,l)*dzdx + vvhn(ix,jy,kz,l)*dzdy) ! variation of geopot on sigma is necessary
            + (u*dzdx + v*dzdy) ! variation of geopot on sigma is necessary
             if (kz .eq. 1) wwhn(ix,jy,kz,l) = wwhn(ix,jy,kz,l)*pinmconv(kz)

!             if (kz .eq. 1) wwhn(ix,jy,kz,l) = 0.0
              end do
         endif
         if (wind_option.eq.-1) then
!!        ww(ix,jy,1,n)=wwh(ix,jy,1)
         wwn(ix,jy,1,n,l)=0.
         do iz=2,nz
         wwn(ix,jy,iz,n,l)=wwn(ix,jy,iz-1,n,l)-(height(iz)-height(iz-1))* &
        divn(ix,jy,iz-1,l)
         enddo
         else
!       print*,'converting ww in nest'
          wwn(ix,jy,1,n,l)=wwhn(ix,jy,1,l)
          wwn(ix,jy,nz,n,l)=wwhn(ix,jy,nwz,l)
          kmin=2
          do iz=2,nz
            do kz=kmin,nwz
              if ((height(iz).gt.wzlev(kz-1)).and. &
                  (height(iz).le.wzlev(kz))) then
               dz1=height(iz)-wzlev(kz-1)
               dz2=wzlev(kz)-height(iz)
               dz=dz1+dz2
               wwn(ix,jy,iz,n,l)=(wwhn(ix,jy,kz-1,l)*dz2+ &
         wwhn(ix,jy,kz,l)*dz1)/dz
               kmin=kz
           goto 40
          endif
        end do
40      continue
      end do
         endif 

!          if (method_w_terrain_correction .eq. 20) then
!             do kz = 1, nwz
!                wwhn(ix,jy,kz,l) = wwhn_svaa(kz)
!             end do
!          end if
          end if

! Compute density gradients at intermediate levels
!*************************************************

          drhodzn(ix,jy,1,n,l)=(rhon(ix,jy,2,n,l)-rhon(ix,jy,1,n,l))/ &
            (height(2)-height(1))
          do kz=2,nz-1
          drhodzn(ix,jy,kz,n,l)=(rhon(ix,jy,kz+1,n,l)- &
            rhon(ix,jy,kz-1,n,l))/(height(kz+1)-height(kz-1))
      end do
          drhodzn(ix,jy,nz,n,l)=drhodzn(ix,jy,nz-1,n,l)

    end do
  end do
!$OMP END DO
!$OMP END PARALLEL

!         print*,'end of ww, now clouds, nests'
!****************************************************************
! Compute slope of eta levels in windward direction and resulting
! vertical wind correction
!
! See notes in verttransform.f about the w correction done here.
!****************************************************************
  !write (*,*) 'initializing clouds, n:',n,nymin1,nxmin1,nz^M
  !   create a cloud and rainout/washout field, clouds occur where rh>80%^M
  !   total cloudheight is stored at level 0^M

     do 100 jy=0,nyn(l)-1
       do 100 ix=0,nxn(l)-1
!        rain_cloud_above=0
         lsp=lsprecn(ix,jy,1,n,l)
         convp=convprecn(ix,jy,1,n,l)

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

          icloudbotn(ix,jy,n,l)=icmv
          icloudtop=icmv ! this is just a local variable
98        do kz=1,nz
            pressure=rhon(ix,jy,kz,n,l)*r_air*ttn(ix,jy,kz,n,l)
            rh=qvn(ix,jy,kz,n,l)/f_qvsat(pressure,ttn(ix,jy,kz,n,l))
!cps            if (prec.gt.0.01) print*,'relhum',prec,kz,rh,height(kz)
            if (rh .gt. rhmin) then
              if (icloudbotn(ix,jy,n,l) .eq. icmv) then
                icloudbotn(ix,jy,n,l)=nint(height(kz))
              endif
              icloudtop=nint(height(kz)) ! use int to save memory
            endif
          enddo

!CPS try to get a cloud thicker than 50 m
!CPS if there is at least .01 mm/h  - changed to 0.002 and put into
!CPS parameter precpmin
          if ((icloudbotn(ix,jy,n,l) .eq. icmv .or. &
               icloudtop-icloudbotn(ix,jy,n,l) .lt. 50) .and. &
               prec .gt. precmin) then
            rhmin = rhmin - 0.05
            if (rhmin .ge. 0.30) goto 98 ! give up for <= 25% rel.hum.
          endif
!CPS implement a rough fix for badly represented convection
!CPS is based on looking at a limited set of comparison data
          if (lconvectprec .and. icloudtop .lt. 6000 .and. &
              prec .gt. precmin) then
            if (convp .lt. 0.1) then
              icloudbotn(ix,jy,n,l) = 500
              icloudtop =         8000
            else
              icloudbotn(ix,jy,n,l) = 0
              icloudtop =      10000
            endif
          endif
          if (icloudtop .ne. icmv) then
            icloudthckn(ix,jy,n,l) = icloudtop-icloudbotn(ix,jy,n,l)
          else
            icloudthckn(ix,jy,n,l) = icmv
          endif
!CPS  get rid of too thin clouds
          if (icloudthckn(ix,jy,n,l) .lt. 50) then
            icloudbotn(ix,jy,n,l)=icmv
            icloudthckn(ix,jy,n,l)=icmv
          endif

100   continue
      enddo ! nests    

      return
      end
