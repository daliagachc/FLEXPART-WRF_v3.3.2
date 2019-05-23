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
      subroutine initialize(itime,ldt,up,vp,wp, &
      usigold,vsigold,wsigold,xt,yt,zt,icbt, &
      ngrid,depoindicator,indzindicator,cpt2,ompid,myid,n_threads,mts )
!    uprof,vprof,wprof, usigprof,vsigprof,wsigprof, &
!    rhoprof,rhogradprof, tkeprof,pttprof, &
!    u,v,w,usig,vsig,wsig,pvi, &
!!   p1,p2,p3,p4,ddx,ddy,rddx,rddy,dtt,dt1,dt2, &
!    ix,jy,ixp,jyp,ngrid,indz,indzp, depoindicator, &
!    indzindicator, &
!    ust,wst,ol,h,zeta,sigu,sigv,tlu,tlv,tlw, &
!    sigw,dsigwdz,dsigw2dz,cpt,ompid)

!                             i    i   o  o  o 
!        o       o       o    i  i  i   o
!*******************************************************************************
!                                                                              *
!  Calculation of trajectories utilizing a zero-acceleration scheme.           *
!  The time step is determined by the Courant-Friedrichs-Lewy (CFL) criterion. *
!  This means that the time step must be so small that the displacement within *
!  this time step is smaller than 1 grid distance. Additionally, a temporal    *
!  CFL criterion is introduced: the time step must be smaller than the time    *
!  interval of the wind fields used for interpolation.                         *
!  For random walk simulations, these are the only time step criteria.         *
!  For the other options, the time step is also limited by the Lagrangian time *
!  scale.                                                                      *
!                                                                              *
!     Author: A. Stohl                                                         *
!                                                                              *
!     16 December 1997                                                         *
!                                                                              *
!  Literature:                                                                 *
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
! h [m]              Mixing height                                             *
! lwindinterv [s]    time interval between two wind fields                     *
! itime [s]          current temporal position                                 *
! ldt [s]            Suggested time step for next integration                  *
! ladvance [s]       Total integration time period                             *
! rannumb(maxrand)   normally distributed random variables                     *
! up,vp,wp           random velocities due to turbulence                       *
! usig,vsig,wsig     uncertainties of wind velocities due to interpolation     *
! usigold,vsigold,wsigold  like usig, etc., but for the last time step         *
! xt,yt,zt           Next time step's spatial position of trajectory           *
!                                                                              *
!                                                                              *
! Constants:                                                                   *
! cfl                factor, by which the time step has to be smaller than the *
!                    spatial CFL-criterion                                     *
! cflt               factor, by which the time step has to be smaller than the *
!                    temporal CFL-criterion                                    *
!                                                                              *
!*******************************************************************************
! 12 JUNE 2007 W. Wang
!              use WRF TKE option to compute turbulence 
!  Mar 2012: J. Brioude modification to handle openmp.                         *
! Jan 2013 M. Cassiani modification to use CBL scheme                            
!*******************************************************************************
  use par_mod
  use com_mod
  use mt_stream

!  use interpol_mod
!  use hanna_mod
!  use ran_mod
  implicit none

  integer :: itime
  integer :: ldt,nrand,ompid
!OMP_GET_THREAD_NUM
  integer(kind=2) :: icbt
  real :: zt,dz,dz1,dz2,up,vp,wp,usigold,vsigold,wsigold,ran3
  real(kind=dp) :: xt,yt
!  save idummy

  integer :: idummy = -7

  real :: uprof(nzmax),vprof(nzmax),wprof(nzmax)
  real :: usigprof(nzmax),vsigprof(nzmax),wsigprof(nzmax)
  real :: rhoprof(nzmax),rhogradprof(nzmax)
  real :: tkeprof(nzmax),pttprof(nzmax)
  real :: u,v,w,usig,vsig,wsig,pvi,mu,mv

  real :: p1,p2,p3,p4,ddx,ddy,rddx,rddy,dtt,dt1,dt2
  integer :: ix,jy,ixp,jyp,ngrid,indz,indzp,cpt2,maxrand2
  logical :: depoindicator(maxspec)
  logical :: indzindicator(nzmax)

  real :: ust,wst,ol,h,zeta,sigu,sigv,tlu,tlv,tlw
  real :: sigw,dsigwdz,dsigw2dz

  real :: dcas,dcas1,dcas2  !modified by  by mc, random number needed in initialize_cbl_vel
  integer ::  myid,n_threads !added by mc for parallel random number generation
  integer(4) :: rannum
  real(4) :: real_rannum
  type (mt_state) :: mts (0: MAX_STREAM)


  idummy=7
  icbt=1           ! initialize particle to "no reflection"

       if (newrandomgen.eq.0) then
       cpt2=cpt2+1
!      cpt=cpt+1000367
      cpt2=mod(cpt2,maxrand)+1;

!      nrand=int(ran3(idummy,inext,inextp,ma,iff)*real(maxrand-1))+1
!      nrand=int(ran3(idummy)*real(maxrand-1))+1
      nrand=cpt2+ompid*maxrand
       maxrand2=maxrandomp
       else
!mc
!      rannum=genrand_int32(mts(ompid+1+(myid*n_threads)))  !integer random
!      number at 32 bit resolution
       rannum=genrand_int32(mts(ompid+1))  !integer random number at 32 bit resolution
       real_rannum = sngl(0.5_DP + 0.2328306e-9_DP * rannum) !conversion to single precision 32bit real between 0-1
       nrand=int(real_rannum*real(maxrand-1))+1
       maxrand2=maxrand
       endif

!******************************
! 2. Interpolate necessary data
!******************************

! Compute maximum mixing height around particle position
!*******************************************************

      ix=int(xt)
      jy=int(yt)
      ixp=ix+1
      jyp=jy+1

  h=max(hmix(ix ,jy ,1,memind(1)), &
       hmix(ixp,jy ,1,memind(1)), &
       hmix(ix ,jyp,1,memind(1)), &
       hmix(ixp,jyp,1,memind(1)), &
       hmix(ix ,jy ,1,memind(2)), &
       hmix(ixp,jy ,1,memind(2)), &
       hmix(ix ,jyp,1,memind(2)), &
       hmix(ixp,jyp,1,memind(2)))

! JB 
      zeta=zt/h

!*************************************************************
! If particle is in the PBL, interpolate once and then make a
! time loop until end of interval is reached
!*************************************************************

      if (zeta.le.1.) then

        call interpol_all(itime,real(xt),real(yt),zt, &
    uprof,vprof,wprof, usigprof,vsigprof,wsigprof, &
    rhoprof,rhogradprof, tkeprof,pttprof, &
    u,v,w,usig,vsig,wsig,pvi, &
    p1,p2,p3,p4,ddx,ddy,rddx,rddy,dtt,dt1,dt2, &
    ix,jy,ixp,jyp,ngrid,indz,indzp, depoindicator, &
    indzindicator, &
    ust,wst,ol,h,zeta,sigu,sigv,tlu,tlv,tlw, &
    sigw,dsigwdz,dsigw2dz,mu,mv)


! Vertical interpolation of u,v,w,rho and drhodz
!***********************************************

! Vertical distance to the level below and above current position
! both in terms of (u,v) and (w) fields
!****************************************************************

        dz1=zt-height(indz)
        dz2=height(indzp)-zt
        dz=1./(dz1+dz2)

        u=(dz1*uprof(indzp)+dz2*uprof(indz))*dz
        v=(dz1*vprof(indzp)+dz2*vprof(indz))*dz
        w=(dz1*wprof(indzp)+dz2*wprof(indz))*dz


! Compute the turbulent disturbances

! Determine the sigmas and the timescales
!****************************************
! FLEXPART WRF
!          write(*,*)'initial.f','turb_option=',turb_option
!          write(*,*)'turb_option_mytke=',turb_option_mytke
       if (turb_option .eq. turb_option_mytke) then
!           write(*,*)'initial.f'
            call tke_partition_my(zt, &
   ust,wst,ol,h,zeta,sigu,sigv,tlu,tlv,tlw, &
   sigw,dsigwdz,dsigw2dz,uprof,vprof,tkeprof,pttprof,indz,indzp)
       elseif (turb_option .eq. turb_option_tke) then
              call tke_partition_hanna(zt, &
   ust,wst,ol,h,zeta,sigu,sigv,tlu,tlv,tlw, &
   sigw,dsigwdz,dsigw2dz,uprof,vprof,tkeprof,pttprof,indz,indzp)
       else

         if (turbswitch) then
           call hanna(zt, &
   ust,wst,ol,h,zeta,sigu,sigv,tlu,tlv,tlw, &
   sigw,dsigwdz,dsigw2dz)

         else
           call hanna1(zt, &
   ust,wst,ol,h,zeta,sigu,sigv,tlu,tlv,tlw, &
   sigw,dsigwdz,dsigw2dz)
         endif
       endif       

! Determine the new diffusivity velocities
!*****************************************

        if (nrand+2.gt.maxrand2) nrand=1
        up=rannumb(nrand)*sigu
        vp=rannumb(nrand+1)*sigv
        wp=rannumb(nrand+2)
        nrand=nrand+2 ! added by mc: it was missing previously a bug I think here and in original flexpart

        if (.not.turbswitch) then     ! modified by mc
        wp=wp*sigw
       else if (cblflag.eq.1) then   ! modified by mc
        if (-h/ol.gt.5) then          !unstable conditions from -h/ol >5 
        !if (ol.lt.0.) then
        !if (ol.gt.0.) then !by mc : gt.0 is only for test the correct is lt.0^M
            dcas=uniform_rannumb(nrand) !uniform^M
            dcas1=rannumb(nrand)        !gaussian^M
            nrand=nrand+3
            call initialize_cbl_vel(idummy,zt,ust,wst,h,sigw,wp,dcas,dcas1,ol)
        else
            wp=wp*sigw
        end if
       end if
! Determine time step for next integration
!*****************************************

        if (turbswitch) then
          ldt=int(min(tlw,h/max(2.*abs(wp*sigw),1.e-5), &
          0.5/abs(dsigwdz),600.)*ctl)
        else
          ldt=int(min(tlw,h/max(2.*abs(wp),1.e-5),600.)*ctl)
        endif
        ldt=max(ldt,mintime)


        usig=(usigprof(indzp)+usigprof(indz))/2.
        vsig=(vsigprof(indzp)+vsigprof(indz))/2.
        wsig=(wsigprof(indzp)+wsigprof(indz))/2.

      else



!**********************************************************
! For all particles that are outside the PBL, make a single
! time step. Only horizontal turbulent disturbances are
! calculated. Vertical disturbances are reset.
!**********************************************************


! Interpolate the wind
!*********************

        call interpol_wind(itime,real(xt),real(yt),zt, &
    uprof,vprof,wprof, usigprof,vsigprof,wsigprof, &
    rhoprof,rhogradprof, tkeprof,pttprof, &
    u,v,w,usig,vsig,wsig,pvi, &
    p1,p2,p3,p4,ddx,ddy,rddx,rddy,dtt,dt1,dt2, &
    ix,jy,ixp,jyp,ngrid,indz,indzp, depoindicator, &
    indzindicator,mu,mv)


! Compute everything for above the PBL

! Assume constant turbulent perturbations
!****************************************

        ldt=abs(lsynctime)

        if (nrand+1.gt.maxrand2) nrand=1
        up=rannumb(nrand)*0.3
        vp=rannumb(nrand+1)*0.3
        nrand=nrand+2
        wp=0.
        sigw=0.

      endif

!****************************************************************
! Add mesoscale random disturbances
! This is done only once for the whole lsynctime interval to save
! computation time
!****************************************************************


! It is assumed that the average interpolation error is 1/2 sigma
! of the surrounding points, autocorrelation time constant is
! 1/2 of time interval between wind fields
!****************************************************************

      if (nrand+2.gt.maxrand2) nrand=1
      usigold=rannumb(nrand)*usig
      vsigold=rannumb(nrand+1)*vsig
      wsigold=rannumb(nrand+2)*wsig

end subroutine initialize

