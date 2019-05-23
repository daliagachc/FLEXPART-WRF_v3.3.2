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
      subroutine interpol_all_nests(itime,xt,yt,zt, &
    uprof,vprof,wprof, usigprof,vsigprof,wsigprof, &
    rhoprof,rhogradprof, tkeprof,pttprof, &
    u,v,w,usig,vsig,wsig,pvi, &
    p1,p2,p3,p4,ddx,ddy,rddx,rddy,dtt,dt1,dt2, &
    ix,jy,ixp,jyp,ngrid,indz,indzp, depoindicator, &
    indzindicator, &
   ust,wst,ol,h,zeta,sigu,sigv,tlu,tlv,tlw, &
   sigw,dsigwdz,dsigw2dz,mu,mv)

!                                     i   i  i  i
!*******************************************************************************
!                                                                              *
!  This subroutine interpolates everything that is needed for calculating the  *
!  dispersion.                                                                 *
!  Version for interpolating nested grids.                                     *
!                                                                              *
!    Author: A. Stohl                                                          *
!                                                                              *
!    9 February 1999                                                           *
!    16 December 1997                                                          *
!                                                                              *
!  Revision March 2005 by AST : all output variables in common block           *
!                               calculation of standard deviation done in this *
!                               routine rather than subroutine call in order   *
!                               to save computation time                       *
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
! itime [s]          current temporal position                                 *
! memtime(3) [s]     times of the wind fields in memory                        *
! xt,yt,zt           coordinates position for which wind data shall be calculat*
!                                                                              *
! Constants:                                                                   *
!                                                                              *
!*******************************************************************************
!     12 JUNE 2007  compute tkeprof,  add a variable y4(2)
!     25 June 2007  compute pttprof, y5
!                   compute tkeprof for all levels
!*******************************************************************************
  use par_mod
  use com_mod
!  use interpol_mod
!  use hanna_mod

  implicit none

  integer :: itime
  real :: xt,yt,zt

  ! Auxiliary variables needed for interpolation
  real :: ust1(2),wst1(2),oli1(2),oliaux
  real :: y1(2),y2(2),y3(2),rho1(2),rhograd1(2),y4(2),y5(2)
  real :: usl,vsl,wsl,usq,vsq,wsq,xaux
  integer :: i,m,n,indexh
  real,parameter :: eps=1.0e-30
  real :: uprof(nzmax),vprof(nzmax),wprof(nzmax)
  real :: usigprof(nzmax),vsigprof(nzmax),wsigprof(nzmax)
  real :: rhoprof(nzmax),rhogradprof(nzmax)
  real :: tkeprof(nzmax),pttprof(nzmax)
  real :: u,v,w,usig,vsig,wsig,pvi,mu,mv

  real :: p1,p2,p3,p4,ddx,ddy,rddx,rddy,dtt,dt1,dt2
  integer :: ix,jy,ixp,jyp,ngrid,indz,indzp
  logical :: depoindicator(maxspec)
  logical :: indzindicator(nzmax)
  real :: ust,wst,ol,h,zeta,sigu,sigv,tlu,tlv,tlw
  real :: sigw,dsigwdz,dsigw2dz


!********************************************
! Multilinear interpolation in time and space
!********************************************

! Determine the lower left corner and its distance to the current position
!*************************************************************************

      ddx=xt-real(ix)                     
      ddy=yt-real(jy)
      rddx=1.-ddx
      rddy=1.-ddy
      p1=rddx*rddy
      p2=ddx*rddy
      p3=rddx*ddy
      p4=ddx*ddy

! Calculate variables for time interpolation
!*******************************************

      dt1=real(itime-memtime(1))
      dt2=real(memtime(2)-itime)
      dtt=1./(dt1+dt2)


!*****************************************
! 1. Interpolate u*, w* and Obukhov length
!*****************************************

! a) Bilinear horizontal interpolation
  do m=1,2
    indexh=memind(m)

    ust1(m)=p1*ustarn(ix ,jy ,1,indexh,ngrid) &
         + p2*ustarn(ixp,jy ,1,indexh,ngrid) &
         + p3*ustarn(ix ,jyp,1,indexh,ngrid) &
         + p4*ustarn(ixp,jyp,1,indexh,ngrid)
    wst1(m)=p1*wstarn(ix ,jy ,1,indexh,ngrid) &
         + p2*wstarn(ixp,jy ,1,indexh,ngrid) &
         + p3*wstarn(ix ,jyp,1,indexh,ngrid) &
         + p4*wstarn(ixp,jyp,1,indexh,ngrid)
    oli1(m)=p1*olin(ix ,jy ,1,indexh,ngrid) &
         + p2*olin(ixp,jy ,1,indexh,ngrid) &
         + p3*olin(ix ,jyp,1,indexh,ngrid) &
         + p4*olin(ixp,jyp,1,indexh,ngrid)
  end do
       mu =p1*m_xn(ix ,jy ,1,ngrid) &
         + p2*m_xn(ixp,jy ,1,ngrid) &
         + p3*m_xn(ix ,jyp,1,ngrid) &
         + p4*m_xn(ixp,jyp,1,ngrid)
       mv =p1*m_yn(ix ,jy ,1,ngrid) &
         + p2*m_yn(ixp,jy ,1,ngrid) &
         + p3*m_yn(ix ,jyp,1,ngrid) &
         + p4*m_yn(ixp,jyp,1,ngrid)


! b) Temporal interpolation

      ust=(ust1(1)*dt2+ust1(2)*dt1)*dtt
      wst=(wst1(1)*dt2+wst1(2)*dt1)*dtt
      oliaux=(oli1(1)*dt2+oli1(2)*dt1)*dtt

      if (oliaux.ne.0.) then
        ol=1./oliaux
      else
        ol=99999.
      endif


!*****************************************************
! 2. Interpolate vertical profiles of u,v,w,rho,drhodz
!*****************************************************


! Determine the level below the current position
!***********************************************

  do i=2,nz
    if (height(i).gt.zt) then
      indz=i-1
      indzp=i
      goto 6
    endif
  end do
6   continue

!**************************************
! 1.) Bilinear horizontal interpolation
! 2.) Temporal interpolation (linear)
!**************************************

! Loop over 2 time steps and indz levels
!***************************************

  do n=indz,indz+1
    usl=0.
    vsl=0.
    wsl=0.
    usq=0.
    vsq=0.
    wsq=0.
    do m=1,2
      indexh=memind(m)
      y1(m)=p1*uun(ix ,jy ,n,indexh,ngrid) &
           +p2*uun(ixp,jy ,n,indexh,ngrid) &
           +p3*uun(ix ,jyp,n,indexh,ngrid) &
           +p4*uun(ixp,jyp,n,indexh,ngrid)
      y2(m)=p1*vvn(ix ,jy ,n,indexh,ngrid) &
           +p2*vvn(ixp,jy ,n,indexh,ngrid) &
           +p3*vvn(ix ,jyp,n,indexh,ngrid) &
           +p4*vvn(ixp,jyp,n,indexh,ngrid)
      y3(m)=p1*wwn(ix ,jy ,n,indexh,ngrid) &
           +p2*wwn(ixp,jy ,n,indexh,ngrid) &
           +p3*wwn(ix ,jyp,n,indexh,ngrid) &
           +p4*wwn(ixp,jyp,n,indexh,ngrid)
      rhograd1(m)=p1*drhodzn(ix ,jy ,n,indexh,ngrid) &
           +p2*drhodzn(ixp,jy ,n,indexh,ngrid) &
           +p3*drhodzn(ix ,jyp,n,indexh,ngrid) &
           +p4*drhodzn(ixp,jyp,n,indexh,ngrid)
      rho1(m)=p1*rhon(ix ,jy ,n,indexh,ngrid) &
           +p2*rhon(ixp,jy ,n,indexh,ngrid) &
           +p3*rhon(ix ,jyp,n,indexh,ngrid) &
           +p4*rhon(ixp,jyp,n,indexh,ngrid)


     usl=usl+uun(ix ,jy ,n,indexh,ngrid)+uun(ixp,jy ,n,indexh,ngrid) &
          +uun(ix ,jyp,n,indexh,ngrid)+uun(ixp,jyp,n,indexh,ngrid)
     vsl=vsl+vvn(ix ,jy ,n,indexh,ngrid)+vvn(ixp,jy ,n,indexh,ngrid) &
          +vvn(ix ,jyp,n,indexh,ngrid)+vvn(ixp,jyp,n,indexh,ngrid)
     wsl=wsl+wwn(ix ,jy ,n,indexh,ngrid)+wwn(ixp,jy ,n,indexh,ngrid) &
          +wwn(ix ,jyp,n,indexh,ngrid)+wwn(ixp,jyp,n,indexh,ngrid)

    usq=usq+uun(ix ,jy ,n,indexh,ngrid)*uun(ix ,jy ,n,indexh,ngrid)+ &
         uun(ixp,jy ,n,indexh,ngrid)*uun(ixp,jy ,n,indexh,ngrid)+ &
         uun(ix ,jyp,n,indexh,ngrid)*uun(ix ,jyp,n,indexh,ngrid)+ &
         uun(ixp,jyp,n,indexh,ngrid)*uun(ixp,jyp,n,indexh,ngrid)
    vsq=vsq+vvn(ix ,jy ,n,indexh,ngrid)*vvn(ix ,jy ,n,indexh,ngrid)+ &
         vvn(ixp,jy ,n,indexh,ngrid)*vvn(ixp,jy ,n,indexh,ngrid)+ &
         vvn(ix ,jyp,n,indexh,ngrid)*vvn(ix ,jyp,n,indexh,ngrid)+ &
         vvn(ixp,jyp,n,indexh,ngrid)*vvn(ixp,jyp,n,indexh,ngrid)
    wsq=wsq+wwn(ix ,jy ,n,indexh,ngrid)*wwn(ix ,jy ,n,indexh,ngrid)+ &
         wwn(ixp,jy ,n,indexh,ngrid)*wwn(ixp,jy ,n,indexh,ngrid)+ &
         wwn(ix ,jyp,n,indexh,ngrid)*wwn(ix ,jyp,n,indexh,ngrid)+ &
         wwn(ixp,jyp,n,indexh,ngrid)*wwn(ixp,jyp,n,indexh,ngrid)
    end do
        uprof(n)=(y1(1)*dt2+y1(2)*dt1)*dtt
        vprof(n)=(y2(1)*dt2+y2(2)*dt1)*dtt
        wprof(n)=(y3(1)*dt2+y3(2)*dt1)*dtt
        rhoprof(n)=(rho1(1)*dt2+rho1(2)*dt1)*dtt
        rhogradprof(n)=(rhograd1(1)*dt2+rhograd1(2)*dt1)*dtt
        indzindicator(n)=.false.

! Compute standard deviations
!****************************

        xaux=usq-usl*usl/8.
        if (xaux.lt.eps) then
          usigprof(n)=0.
        else
          usigprof(n)=sqrt(xaux/7.)
        endif

        xaux=vsq-vsl*vsl/8.
        if (xaux.lt.eps) then
          vsigprof(n)=0.
        else
          vsigprof(n)=sqrt(xaux/7.)
        endif


        xaux=wsq-wsl*wsl/8.
        if (xaux.lt.eps) then
          wsigprof(n)=0.
        else
          wsigprof(n)=sqrt(xaux/7.)
        endif

  end do

! compute tke profile for all levels

        do n=1,nz
         do m=1,2
             indexh=memind(m)
                    y4(m)=p1*tken(ix ,jy ,n,indexh,ngrid) &
                 +p2*tken(ixp,jy ,n,indexh,ngrid) &
                 +p3*tken(ix ,jyp,n,indexh,ngrid) &
                 +p4*tken(ixp,jyp,n,indexh,ngrid)
          y5(m)=p1*pttn(ix ,jy ,n,indexh,ngrid) &
                 +p2*pttn(ixp,jy ,n,indexh,ngrid) &
                 +p3*pttn(ix ,jyp,n,indexh,ngrid) &
                 +p4*pttn(ixp,jyp,n,indexh,ngrid)

        enddo
       tkeprof(n)=(y4(1)*dt2+y4(2)*dt1)*dtt
         pttprof(n)=(y5(1)*dt2+y5(2)*dt1)*dtt
        enddo


     end subroutine interpol_all_nests

