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

      subroutine tke_partition_hanna(z, &
   ust,wst,ol,h,zeta,sigu,sigv,tlu,tlv,tlw, &
   sigw,dsigwdz,dsigw2dz,uprof,vprof,tkeprof,pttprof,indz,indzp)
!                      i
!*******************************************************************************
!                                                                              *
!     Computation of \sigma_u,v,w,dsigwdz, and dsigw2dz based on TKE from WRF  *
!                                                                              *
!     Author: A. Stohl                                                         *
!                                                                              *
!     4 December 1997                                                          *
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
! dsigwdz [1/s]     vertical gradient of sigw                                  *
! ol [m]            Obukhov length                                             *
! sigu, sigv, sigw  standard deviations of turbulent velocity fluctuations     *
! tlu [s]           Lagrangian time scale for the along wind component.        *
! tlv [s]           Lagrangian time scale for the cross wind component.        *
! tlw [s]           Lagrangian time scale for the vertical wind component.     *
! ust, ustar [m/s]  friction velocity                                          *
! wst, wstar [m/s]  convective velocity scale                                  *
!                                                                              *
!*******************************************************************************
! 12 JUNE 2007  compute sigu,sigv,sigw,dsigwdz, and dsigw2dz from TKE          *
! fu2               fraction for u2
! fv2               fraction for v2
! fw2               fraction for w2
!
! 25 JUNE 2007 merged from hanna_tke and tke_partition.f

!      include 'includepar'
!      include 'includecom'
!      include 'includehanna'
!      include 'includeinterpol'
  use par_mod
  use com_mod
!  use hanna_mod
!  use interpol_mod
   implicit none 
      real :: corr,z,zzz,hanna_sigu,hanna_sigv,hanna_sigw,dz,dz1,dz2
      real :: fu2(2),fv2(2),fw2(2)
      real :: siguprof(2),sigvprof(2),sigwprof(2)
      integer :: k,indz,indzp
     real :: uprof(nzmax),vprof(nzmax),tkeprof(nzmax),pttprof(nzmax)
    real :: ust,wst,ol,h,zeta,sigu,sigv,tlu,tlv,tlw
    real :: sigw,dsigwdz,dsigw2dz

       do k=1,2 
         fu2(k)=0.333
         fv2(k)=0.333
         fw2(k)=0.333        
         siguprof(k)=0.0
         sigvprof(k)=0.0
         sigwprof(k)=0.0
       enddo

! determine fraction
       do k=1,2                                   ! k fraction      
!**********************
! 1. Neutral conditions
!**********************
        if (k .eq. 1) zzz=height(indz)
        if (k .eq. 2) zzz=height(indzp)
      if (h/abs(ol).lt.1.) then
        ust=max(1.e-4,ust)
        corr=height(indz+k-1)/ust
        hanna_sigu=1.e-2+2.0*ust*exp(-3.e-4*corr)
        hanna_sigw=1.e-2+1.3*ust*exp(-2.e-4*corr)
        hanna_sigv=hanna_sigw

!***********************
! 2. Unstable conditions
!***********************

      else if (ol.lt.0.) then


! Determine sigmas
!*****************

        hanna_sigu=1.e-2+ust*(12-0.5*h/ol)**0.33333
        hanna_sigv=hanna_sigu
        hanna_sigw=sqrt(1.2*wst**2*(1.-.9*zeta)*height(indz+k-1)/h &
       **0.66666+(1.8-1.4*height(k)/h)*ust**2)+1.e-2

!*********************
! 3. Stable conditions
!*********************
      else
        hanna_sigu=1.e-2+2.*ust*(1.-height(indz+k-1)/h)
        hanna_sigv=1.e-2+1.3*ust*(1.-height(indz+k-1)/h)
        hanna_sigw=hanna_sigv
      endif
        fu2(k)=hanna_sigu**2/(hanna_sigu**2+hanna_sigv**2+hanna_sigw**2)
        fv2(k)=hanna_sigv**2/(hanna_sigu**2+hanna_sigv**2+hanna_sigw**2)
        fw2(k)=hanna_sigw**2/(hanna_sigu**2+hanna_sigv**2+hanna_sigw**2)
      enddo                                          !k fraction

!- compute sigu,v,w
      do k=1,2             ! siguprof
        siguprof(k)=max(sqrt(2.0*tkeprof(indz+k-1)*fu2(k)),1.e-2)
        sigvprof(k)=max(sqrt(2.0*tkeprof(indz+k-1)*fv2(k)),1.e-2)
        sigwprof(k)=max(sqrt(2.0*tkeprof(indz+k-1)*fw2(k)),1.e-2)
!C         write(*,*)'z=',height(indz+k-1),'tke=', tkeprof(indz+k-1)
      enddo                 ! siguprof
!         write(*,*)'tkeprof=',(tkeprof(k),k=1,nz)
!- interpolate sigu,sigv, sigw
          dz=1./(height(indzp)-height(indz))
          dz1=(z - height(indz))*dz
          dz2=(height(indzp)-z)*dz 
           
          sigu=dz1*siguprof(2)+dz2*siguprof(1)
          sigv=dz1*sigvprof(2)+dz2*sigvprof(1)
          sigw=dz1*sigwprof(2)+dz2*sigwprof(1)
          
          dsigwdz=max(1.e-10,(sigwprof(2)-sigwprof(1))*dz )          
          dsigw2dz=max(1.e-10,(sigwprof(2)**2-sigwprof(1)**2)*dz )

!-- compute length scales based on hanna(1982)
!*************************
! 1. Neutral conditions
!**********************
 
      if (h/abs(ol).lt.1.) then
        tlu=0.5*z/sigw/(1.+1.5e-3*corr)
        tlv=tlu
        tlw=tlu
 
!***********************
! 2. Unstable conditions
!***********************
 
      else if (ol.lt.0.) then
 
! Determine average Lagrangian time scale
!****************************************
 
        tlu=0.15*h/sigu
        tlv=tlu
        if (z.lt.abs(ol)) then
          tlw=0.1*z/(sigw*(0.55-0.38*abs(z/ol)))
        else if (zeta.lt.0.1) then
          tlw=0.59*z/sigw
        else
          tlw=0.15*h/sigw*(1.-exp(-5*zeta))
        endif
 
 
!*********************
! 3. Stable conditions
!*********************
 
      else
        tlu=0.15*h/sigu*(sqrt(zeta))
        tlv=0.467*tlu
        tlw=0.1*h/sigw*zeta**0.8
      endif
 
      tlu=max(10.,tlu)
      tlv=max(10.,tlv)
      tlw=max(30.,tlw)
 




      end subroutine tke_partition_hanna
