!***********************************************************************
!* Copyright 2012,2013                                                *
!* Jerome Brioude, Delia Arnold, Andreas Stohl, Wayne Angevine,       *
!* John Burkhart, Massimo Cassiani, Adam Dingwell, Richard C Easter, Sabine Eckhardt,*
!* Stephanie Evan, Jerome D Fast, Don Morton, Ignacio Pisso,          *
!* Petra Seibert, Gerard Wotawa, Caroline Forster, Harald Sodemann,   *
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

      subroutine advance(itime,nrelpoint,ldt,up,vp,wp, &
      usigold,vsigold,wsigold,nstop,xt,yt,zt,prob,icbt, &
      ngrid,depoindicator,indzindicator,cpt2,ompid,myid,n_threads,mts)  !comment by mc: ...,ompid,myid,n_threads) added  n_threads for MT parallel random number generator

!                          i    i  i/oi/oi/o
!       i/o     i/o     i/o     o  i/oi/oi/o i/o  i/o
!*******************************************************************************
!                                                                              *
!  Note:  This is the FLEXPART_WRF version of subroutine gridcheck.            *
!    The computational grid is the WRF x-y grid rather than lat-lon.           *
!                                                                              *
!  Calculation of turbulent particle trajectories utilizing a                  *
!  zero-acceleration scheme, which is corrected by a numerically more accurate *
!  Petterssen scheme whenever possible.                                        *
!                                                                              *
!  Particle positions are read in, incremented, and returned to the calling    *
!  program.                                                                    *
!                                                                              *
!  In different regions of the atmosphere (PBL vs. free troposphere),          *
!  different parameters are needed for advection, parameterizing turbulent     *
!  velocities, etc. For efficiency, different interpolation routines have      *
!  been written for these different cases, with the disadvantage that there    *
!  exist several routines doing almost the same. They all share the            *
!  included file 'includeinterpol'. The following                              *
!  interpolation routines are used:                                            *
!                                                                              *
!  interpol_all(_nests)     interpolates everything (called inside the PBL)    *
!  interpol_misslev(_nests) if a particle moves vertically in the PBL,         *
!                           additional parameters are interpolated if it       *
!                           crosses a model level                              *
!  interpol_wind(_nests)    interpolates the wind and determines the           *
!                           standard deviation of the wind (called outside PBL)*
!                           also interpolates potential vorticity              *
!  interpol_wind_short(_nests) only interpolates the wind (needed for the      *
!                           Petterssen scheme)                                 *
!  interpol_vdep(_nests)    interpolates deposition velocities                 *
!                                                                              *
!                                                                              *
!     Author: A. Stohl                                                         *
!                                                                              *
!     16 December 1997                                                         *
!                                                                              *
!  Changes:                                                                    *
!                                                                              *
!  8 April 2000: Deep convection parameterization                              *
!                                                                              *
!  May 2002: Petterssen scheme introduced                                      *
!                                                                              *
!  26 Oct 2005, R. Easter - changes for horizontal grid in m instead of lat,lon*
!  10 Nov 2005, R. Easter - zero turbulent wind components is                  *
!                           turbulence is turned off                           *
!  Mar 2012, J. Brioude: modification to handle openmp.                        *
!  turbulence option 3 is not going   *
!  to work. but it shouldn't be used anyway ^M
!  Jan 2013  M. Cassiani (look for mc or MC in the code):
!  *^M
!               introduction of CBL skewed turbulence model
!               *^M
!               & parallel random number generation
!               *                       *
!********************************************************************************
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
! cbt                1 if particle not transferred to forbidden state, else -1 *
! dawsave            accumulated displacement in along-wind direction          *
! dcwsave            accumulated displacement in cross-wind direction          *
! dxsave             accumulated displacement in longitude                     *
! dysave             accumulated displacement in latitude                      *
! h [m]              Mixing height                                             *
! lwindinterv [s]    time interval between two wind fields                     *
! itime [s]          time at which this subroutine is entered                  *
! itimec [s]         actual time, which is incremented in this subroutine      *
! href [m]           height for which dry deposition velocity is calculated    *
! ladvance [s]       Total integration time period                             *
! ldirect            1 forward, -1 backward                                    *
! ldt [s]            Time step for the next integration                        *
! lsynctime [s]      Synchronisation interval of FLEXPART                      *
! ngrid              index which grid is to be used                            *
! nrand              index for a variable to be picked from rannumb            *
! nstop              if > 1 particle has left domain and must be stopped       *
! prob               probability of absorption due to dry deposition           *
! rannumb(maxrand)   normally distributed random variables                     *
! rhoa               air density                                               *
! rhograd            vertical gradient of the air density                      *
! up,vp,wp           random velocities due to turbulence (along wind, cross    *
!                    wind, vertical wind                                       *
! usig,vsig,wsig     mesoscale wind fluctuations                               *
! usigold,vsigold,wsigold  like usig, etc., but for the last time step         *
! vdepo              Deposition velocities for all species                     *
! xt,yt,zt           Particle position                                         *
!                                                                              *
!*******************************************************************************

  use point_mod
  use par_mod
  use com_mod
  use mt_stream  !added by mc for random number generation^M
!  use test_well_mod !added by mc for testting well mixed
!  use interpol_mod
!  use hanna_mod
  use cmapf_mod
!  use ieee_arithmetic
! include 'sprng_f.h'
!  use ran_mod
!      include 'includepar'
!      include 'includecom'
!      include 'includeinterpol'
!      include 'includehanna'

  implicit none
  real(kind=dp) :: xt,yt
  real :: zt,xts,yts,weight
  integer :: itime,itimec,nstop,ldt,i,j,k,nrand,loop,memindnext
  integer :: ngr,nix,njy,ks,nsp,nrelpoint,ii,ompid,myid,nombre
  real :: dz,dz1,dz2,xlon,ylat,xpol,ypol,gridsize
  real :: ru,rv,rw,dt,ux,vy,cosfact,xtn,ytn,tropop
  real :: prob(maxspec),up,vp,wp,dxsave,dysave,dawsave
  real :: dcwsave,mu,mv
  real :: usigold,vsigold,wsigold,r,rs
  real :: uold,vold,wold,vdepo(maxspec)
  !real uprof(nzmax),vprof(nzmax),wprof(nzmax)
  !real usigprof(nzmax),vsigprof(nzmax),wsigprof(nzmax)
  !real rhoprof(nzmax),rhogradprof(nzmax)
  real :: rhoa,rhograd,ran3,delz,dtf,rhoaux,dtftlw,uxscale,wpscale
  real :: ptot_lhh,Q_lhh,phi_lhh,ath,bth   !modified by mc added for CBL scheme
  real :: old_wp_buf,del_test              !modified by mc added for CBL scheme re-initlization fo particle after NaN
  integer(kind=2) :: icbt
  real,parameter :: eps=nxmax/3.e5,eps2=1.e-9
  integer :: flagrein                      !re-initialization flag for particles: modified by mc

  real :: uprof(nzmax),vprof(nzmax),wprof(nzmax)
  real :: usigprof(nzmax),vsigprof(nzmax),wsigprof(nzmax)
  real :: rhoprof(nzmax),rhogradprof(nzmax)
  real :: tkeprof(nzmax),pttprof(nzmax)
  real :: u,v,w,usig,vsig,wsig,pvi
  real(kind=dp) :: xtold
  real :: p1,p2,p3,p4,ddx,ddy,rddx,rddy,dtt,dt1,dt2
  integer :: ix,jy,ixp,jyp,ngrid,indz,indzp,cpt2,maxrand2
  logical :: depoindicator(maxspec)
  logical :: indzindicator(nzmax)

  real :: ust,wst,ol,h,zeta,sigu,sigv,tlu,tlv,tlw
  real :: sigw,dsigwdz,dsigw2dz

  real :: wp2, zt2, ust2, wst2, h2, rhoa2, rhograd2, sigw2, &
    dsigwdz2, tlw2, ptot_lhh2, Q_lhh2, phi_lhh2, ath2, bth2, ol2

  logical :: isnan2

!!! CHANGE: TEST OF THE WELL-MIXED CRITERION
!        integer iclass
!        parameter(iclass=10)
!        double precision zacc,tacc,t(iclass),th(0:iclass),hsave
!        logical dump
!        save zacc,tacc,t,th,hsave,dump
!c        itimeod=0
!!! CHANGE

  integer :: idummy = 7
  real    :: settling = 0.
  !added by mc for random number generation ---------------------  
  integer ::  n_threads !added by mc for parallel random number generation
  integer(4) :: rannum
  real(4) :: real_rannum
  type (mt_state) :: mts (0: MAX_STREAM)
  integer,SAVE :: nan_count(max_STREAM)=0
!-------------------------------------------------------


!!! CHANGE: TEST OF THE WELL-MIXED CRITERION
!      if (idummy.eq.-7) then
!      open(550,file='WELLMIXEDTEST')
!      do 17 i=0,iclass
!17      th(i)=real(i)/real(iclass)
!      endif
!!! CHANGE

!    if (nombre.eq.103) print*,'usig -1',usig,xts,yts,zt
     if (xt.ne.xt .or. abs(xt).gt.1000.) print*,'problem 0', xt,yt,zt,itime,myid,ompid
     xtold=xt
!   print *,'aa',xt,yt,zt,u,v,w
      nstop=0
      do i=1,nmixz
        indzindicator(i)=.true.
       enddo

      if (DRYDEP) then    ! reset probability for deposition
        do  ks=1,nspec
          depoindicator(ks)=.true.
          prob(ks)=0.
       enddo
      endif

      dxsave=0.           ! reset position displacements
      dysave=0.           ! due to mean wind
      dawsave=0.          ! and turbulent wind
      dcwsave=0.
      usig=0.
      vsig=0.
      wsig=0.
       ust=0.
      wst=0.
      ol=0.
      h=0.
      zeta=0.
      sigu=0.
      sigv=0.
      tlu=0.
      tlv=0.
      tlw=0.
      sigw=0.
      dsigwdz=0.
      dsigw2dz=0. 
!     wp=0.
        itimec=itime
      idummy=7
       if (newrandomgen.eq.0) then
!      cpt2=cpt2+ompid+1
      cpt2=cpt2+1
!     cpt2=cpt2+1000367
      cpt2=mod(cpt2,maxrand)+1;

!      nrand=int(ran3(idummy,inext,inextp,ma,iff)*real(maxrand-1))+1
!     nrand=cpt
      nrand=cpt2+ompid*maxrand
!     print*,cpt2,maxrand,maxrandomp,maxomp
!       print*, rannumb(nrand),nrelpoint
!     print*,rannumb(nrand),myid,ompid 
!      if (nrelpoint.ge.993 .and. nrelpoint.le.998) then
!      write(22,*),itime,cpt2,rannumb(nrand),nrelpoint
!!,myid,OMP_GET_THREAD_NUM()
 !!      write(22,*),itime,cpt2,rannumb(nrand),nrelpoint
!         endif


!        if (nrand+2.gt.maxrand) nrand=1
!    print*,rannumb(nrand)
       maxrand2=maxrandomp

        else
!-------------------------------------------------------------------------------------------------
!----- added by MC: parallel random nuymber generation using MT generator ------------------------
       !print *,'varie3',ompid,myid,n_threads
!      rannum=genrand_int32(mts(ompid+1+(myid*n_threads)))  !integer random number at 32 bit resolution
       rannum=genrand_int32(mts(ompid+1))  !integer random number at 32 bit resolution
       real_rannum = sngl(0.5_DP + 0.2328306e-9_DP * rannum) !conversion to single precision 32bit real
       nrand=int(real_rannum*real(maxrand-1))+1
        !print *,'varie4',rannum,real_rannum,nrand
!--------------------------------------------------------------------------------------------------
        maxrand2=maxrand

        end if
       
! Determine whether lat/long grid or polarstereographic projection
! is to be used
! Furthermore, determine which nesting level to be used
!*****************************************************************

      if (nglobal.and.(yt.gt.switchnorthg)) then
        write(*,*)
        write(*,*) '*** stopping in advance ***'
        write(*,*) '    the n-pole code section should not be active'
        write(*,*)
        ngrid=-1
      else if (sglobal.and.(yt.lt.switchsouthg)) then
        write(*,*)
        write(*,*) '*** stopping in advance ***'
        write(*,*) '    the s-pole code section should not be active'
        write(*,*)
        ngrid=-2
      else
        ngrid=0
        do j=numbnests,1,-1
          if ((xt.gt.xln(j)+eps).and.(xt.lt.xrn(j)-eps).and. &
          (yt.gt.yln(j)+eps).and.(yt.lt.yrn(j)-eps)) then
            ngrid=j
            goto 23
          endif
        enddo 
23      continue
      endif


!***************************
! Interpolate necessary data
!***************************

      if (abs(itime-memtime(1)).lt.abs(itime-memtime(2))) then
        memindnext=1
      else
        memindnext=2
      endif

! Determine nested grid coordinates
!**********************************

      if (ngrid.gt.0) then
        xtn=(xt-xln(ngrid))*xresoln(ngrid)
        ytn=(yt-yln(ngrid))*yresoln(ngrid)
        ix=int(xtn)
        jy=int(ytn)
    nix=nint(xtn)
    njy=nint(ytn)

      else
        ix=int(xt)
        jy=int(yt)
    nix=nint(xt)
    njy=nint(yt)
      endif
      ixp=ix+1
      jyp=jy+1

     if (ix.lt.0) print*,'problem', xt,xtold,yt,zt,itime,myid,ompid,nrelpoint

! Compute maximum mixing height around particle position
!*******************************************************

      h=0.
      if (ngrid.le.0) then
        do k=1,2
          do j=jy,jyp
            do i=ix,ixp
              if (hmix(i,j,1,k).gt.h) h=hmix(i,j,1,k)
        end do
      end do
    end do
    tropop=tropopause(nix,njy,1,1)
  else
    do k=1,2
      do j=jy,jyp
        do i=ix,ixp
          if (hmixn(i,j,1,k,ngrid).gt.h) h=hmixn(i,j,1,k,ngrid)
        end do
      end do
    end do
    tropop=tropopausen(nix,njy,1,1,ngrid)
  endif

  zeta=zt/h


!*************************************************************
! If particle is in the PBL, interpolate once and then make a
! time loop until end of interval is reached
!*************************************************************
!      print*,'zeta',zeta,h,zt,xt
      if (zeta.le.1.) then

! BEGIN TIME LOOP
!================

        loop=0
100       loop=loop+1
          if (method.eq.1) then
            ldt=min(ldt,abs(lsynctime-itimec+itime))
            itimec=itimec+ldt*ldirect
          else
            ldt=abs(lsynctime)
            itimec=itime+lsynctime
          endif
          dt=real(ldt)

          zeta=zt/h


!   print *,'xx0',OMP_GET_THREAD_NUM(),loop,xt,yt,zt,xts,yts,u,v,w
          if (loop.eq.1) then
!    if (nombre.eq.103) print*,'usig 0',usig,xt,yt,zt
            if (ngrid.le.0) then
              xts=real(xt)
              yts=real(yt)
!    if (nombre.eq.103) print*,'usig 0',usig,xts,yts,zt
              call interpol_all(itime,xts,yts,zt, &
    uprof,vprof,wprof, usigprof,vsigprof,wsigprof, &
    rhoprof,rhogradprof, tkeprof,pttprof, &
    u,v,w,usig,vsig,wsig,pvi, &
    p1,p2,p3,p4,ddx,ddy,rddx,rddy,dtt,dt1,dt2, &
    ix,jy,ixp,jyp,ngrid,indz,indzp, depoindicator, &
    indzindicator, &
    ust,wst,ol,h,zeta,sigu,sigv,tlu,tlv,tlw, &
    sigw,dsigwdz,dsigw2dz,mu,mv)

            else
              call interpol_all_nests(itime,xtn,ytn,zt, &
    uprof,vprof,wprof, usigprof,vsigprof,wsigprof, &
    rhoprof,rhogradprof, tkeprof,pttprof, &
    u,v,w,usig,vsig,wsig,pvi, &
    p1,p2,p3,p4,ddx,ddy,rddx,rddy,dtt,dt1,dt2, &
    ix,jy,ixp,jyp,ngrid,indz,indzp, depoindicator, &
    indzindicator, &
    ust,wst,ol,h,zeta,sigu,sigv,tlu,tlv,tlw, &
    sigw,dsigwdz,dsigw2dz,mu,mv)
            endif
!    if (nombre.eq.103) print*,'usig 1',usig,xts,yts,zt

          else


!   print *,'xx',OMP_GET_THREAD_NUM(),xt,yt,zt,xts,yts,u,v,w
! Determine the level below the current position for u,v,rho
!***********************************************************

            do i=2,nz
              if (height(i).gt.zt) then
                indz=i-1
                indzp=i
                goto 6
              endif
              enddo
6           continue

! If one of the levels necessary is not yet available,
! calculate it
!*****************************************************

            do i=indz,indzp
              if (indzindicator(i)) then
                if (ngrid.le.0) then
!    if (nombre.eq.103) print*,'in usig 2',usig
                  call interpol_misslev(i,xt,yt,zt, &
    uprof,vprof,wprof, usigprof,vsigprof,wsigprof, &
    rhoprof,rhogradprof, tkeprof,pttprof, &
    u,v,w,usig,vsig,wsig,pvi, &
    p1,p2,p3,p4,ddx,ddy,rddx,rddy,dtt,dt1,dt2, &
    ix,jy,ixp,jyp,ngrid,indz,indzp, depoindicator, &
    indzindicator, &
    ust,wst,ol,h,zeta,sigu,sigv,tlu,tlv,tlw, &
    sigw,dsigwdz,dsigw2dz)
!JB mw not needed here
                else
                  call interpol_misslev_nests(i,xt,yt,zt, &
    uprof,vprof,wprof, usigprof,vsigprof,wsigprof, &
    rhoprof,rhogradprof, tkeprof,pttprof, &
    u,v,w,usig,vsig,wsig,pvi, &
    p1,p2,p3,p4,ddx,ddy,rddx,rddy,dtt,dt1,dt2, &
    ix,jy,ixp,jyp,ngrid,indz,indzp, depoindicator, &
    indzindicator, &
    ust,wst,ol,h,zeta,sigu,sigv,tlu,tlv,tlw, &
    sigw,dsigwdz,dsigw2dz)
                endif
              endif
            enddo
          endif

!    if (nombre.eq.103) print*,'usig 2',usig

! Vertical interpolation of u,v,w,rho and drhodz
!***********************************************

! Vertical distance to the level below and above current position
! both in terms of (u,v) and (w) fields
!****************************************************************

          dz=1./(height(indzp)-height(indz))
          dz1=(zt-height(indz))*dz
          dz2=(height(indzp)-zt)*dz

          u=dz1*uprof(indzp)+dz2*uprof(indz)
          v=dz1*vprof(indzp)+dz2*vprof(indz)
          w=dz1*wprof(indzp)+dz2*wprof(indz)
          rhoa=dz1*rhoprof(indzp)+dz2*rhoprof(indz)
          rhograd=dz1*rhogradprof(indzp)+dz2*rhogradprof(indz)


! Compute the turbulent disturbances
! Determine the sigmas and the timescales
!****************************************
          if (turb_option .eq. turb_option_mytke) then 
! FLEXPART-WRF

!               write(*,*)'itime=',itime,'xt=',xt,'yt=',yt
              call tke_partition_my(zt, &
   ust,wst,ol,h,zeta,sigu,sigv,tlu,tlv,tlw, &
   sigw,dsigwdz,dsigw2dz,uprof,vprof,tkeprof,pttprof,indz,indzp)
          else if (turb_option .eq. turb_option_tke) then
              call tke_partition_hanna(zt, &
   ust,wst,ol,h,zeta,sigu,sigv,tlu,tlv,tlw, &
   sigw,dsigwdz,dsigw2dz,uprof,vprof,tkeprof,pttprof,indz,indzp)
          else
             if (turbswitch) then
               call hanna(zt,  &
    ust,wst,ol,h,zeta,sigu,sigv,tlu,tlv,tlw, &
    sigw,dsigwdz,dsigw2dz)
                       
             else
               call hanna1(zt, &
    ust,wst,ol,h,zeta,sigu,sigv,tlu,tlv,tlw, &
    sigw,dsigwdz,dsigw2dz)

            endif
          endif   
!      print*,   ust,wst,ol,h,zeta,sigu,sigv,tlu,tlv,tlw, &
!  sigw,dsigwdz,dsigw2dz,indz,indzp
!JB
!          if (h/abs(ol).lt.1.) then
!c        print*,itime,'h and ol',h,ol,'neutral'
!           reflect_switch=0
!       else if (ol.lt.0.) then
!c        print*,itime,'h and ol',h,ol,'unstable'
!c          reflect_switch=1
!           reflect_switch=0
!        else
!c        print*,itime,'h and ol',h,ol,'stable'
!           reflect_switch=0
!         endif 
!*****************************************
! Determine the new diffusivity velocities
!*****************************************

! Horizontal components
!**********************

!         if (nrand+1.gt.maxrandomp) nrand=1
          if (nrand+1.gt.maxrand2) nrand=1
          if (dt/tlu.lt..5) then
            up=(1.-dt/tlu)*up+rannumb(nrand)*sigu*sqrt(2.*dt/tlu)
          else
            ru=exp(-dt/tlu)
            up=ru*up+rannumb(nrand)*sigu*sqrt(1.-ru**2)
          endif
          if (dt/tlv.lt..5) then
            vp=(1.-dt/tlv)*vp+rannumb(nrand+1)*sigv*sqrt(2.*dt/tlv)
          else
            rv=exp(-dt/tlv)
            vp=rv*vp+rannumb(nrand+1)*sigv*sqrt(1.-rv**2)
          endif
          nrand=nrand+2


!         if (nrand+ifine.gt.maxrandomp) nrand=1
          if (nrand+ifine.gt.maxrand2) nrand=1
          rhoaux=rhograd/rhoa
          dtf=dt*fine

          dtftlw=dtf/tlw

! Loop over ifine short time steps for vertical component
!********************************************************

          do i=1,ifine

! Determine the drift velocity and density correction velocity
! Determine the drift velocity and density correction velocity
!*************************************************************
!--------------- lines below are teh original FLEXPART  and are commented out to insert teh cbl options comment by mc
!            if (turbswitch) then
!              if (dtftlw.lt..5) then
!                wp=((1.-dtftlw)*wp+rannumb(nrand+i)*sqrt(2.*dtftlw) &
!                +dtf*(dsigwdz+rhoaux*sigw))*real(icbt)
!              else
!                rw=exp(-dtftlw)
!                wp=(rw*wp+rannumb(nrand+i)*sqrt(1.-rw**2) &
!                +tlw*(1.-rw)*(dsigwdz+rhoaux*sigw))*real(icbt)
!              endif
!              delz=wp*sigw*dtf
!            else
!              rw=exp(-dtftlw)
!              wp=(rw*wp+rannumb(nrand+i)*sqrt(1.-rw**2)*sigw &
!              +tlw*(1.-rw)*(dsigw2dz+rhoaux*sigw**2))*real(icbt)
!              delz=wp*dtf
!            endif
!************ CBL scheme integrated in FLEXPART: added by mc **********!
         if (turbswitch) then
          if (dtftlw.lt..5) then
            if (cblflag.eq.1) then
!     wp2=wp
!     zt2=zt
!    ust2=ust
!    wst2=wst
!    h2=h
!    rhoa2=rhoa
!    rhograd2=rhograd
!    sigw2=sigw
!    dsigwdz2=dsigwdz
!    tlw2=tlw
!    ptot_lhh2=ptot_lhh
!    Q_lhh2=Q_lhh
!    phi_lhh2=phi_lhh
!    ath2=ath
!    bth2=bth
!    ol2=ol
               if (-h/ol.gt.5) then  !modified by mc 
              !if (ol.lt.0.) then   !modified by mc  
              !if (ol.gt.0.) then   !modified by mc : for test
                  !print  *,zt,wp,ath,bth,tlw,dtf,'prima'
                  flagrein=0
                  nrand=nrand+1
                  old_wp_buf=wp
                  del_test=(1.-old_wp_buf)/old_wp_buf

                  !rhoa=1.  !for testing vertical well mixed state, by mc
                  !rhograd=0.  !for testing vertical well mixed state, by mc
                  call cbl(wp,zt,ust,wst,h,rhoa,rhograd,sigw,dsigwdz,tlw,ptot_lhh,Q_lhh,phi_lhh,ath,bth,ol,flagrein)
! see inside the routine for inverse time
                  wp=(wp+ath*dtf+bth*rannumb(nrand)*sqrt(dtf))*real(icbt)
                  delz=wp*dtf
!                  if ((ieee_is_nan(zt).or.ieee_is_nan(wp)).and.(flagrein.eq.0)) print*,'pb4',wp2,zt2,ust2,wst2,h2,rhoa2,rhograd2,sigw2,dsigwdz2,tlw2,ptot_lhh2,Q_lhh2,phi_lhh2,ath2,bth2,ol2,flagrein,i
                  if (abs(wp).gt.50.) flagrein=1
                  if (flagrein.eq.1) then  !added for re-initlization of particle vertical velcoity based on condition inside routine cbl.f90
                      call re_initialize_particle(zt,ust,wst,h,sigw,old_wp_buf,nrand,ol)
!                  if (ieee_is_nan(old_wp_buf)) print*,"PROBLEM WP",wp,old_wp_buf,nrand,ol,zt,ust,wst,h,sigw
                      wp=old_wp_buf
                      delz=wp*dtf
                      !nan_count(myid)=nan_count(myid)+1
                       nan_count(ompid+1)=nan_count(ompid+1)+1
                  else       
                  del_test=(1.-wp)/wp !catch infinity value
!                 if (ieee_is_nan(wp) .or. ieee_is_nan(del_test)) then
                   if (isnan2(wp).or.isnan2(del_test)) then !note that, given the test on particle velocity inside the routine cbl.f90, this condition should never be true!! 
!                 if (isnan(wp).or.isnan(del_test)) then !note that, given the test on particle velocity inside the routine cbl.f90, this condition should never be true!! 
                      nrand=nrand+1
                      call re_initialize_particle(zt,ust,wst,h,sigw,old_wp_buf,nrand,ol)
                      wp=old_wp_buf
                      delz=wp*dtf
                      nan_count(ompid+1)=nan_count(ompid+1)+1
                      !nan_count(myid)=nan_count(myid)+1
                      print *,'NaN counter equal to:',nan_count(ompid+1),'omp',ompid,'mpi',myid
!                     ,'increase
!                     !ifine if this number became a non-negligible fraction of
!                     !the particle number'
                  end if
                  end if
              else 
                  !rhoa=1.  !for testing vertical well mixed state, by mc
                  !rhograd=0. !for testing vertical well mixed state, by mc
                  nrand=nrand+1
                  old_wp_buf=wp                  
                  ath=-wp/tlw+sigw*dsigwdz+wp*wp/sigw*dsigwdz+sigw*sigw/rhoa*rhograd  !1-note for inverse time should be -wp/tlw*ldirect+... calculated for wp=-wp
                                                                                      !2-but since ldirect =-1 for inverse time and this must be calculated for (-wp) and
                                                                                      !3-the gaussian pdf is symmetric (i.e. pdf(w)=pdf(-w) ldirect can be discarded
                  bth=sigw*rannumb(nrand)*sqrt(2.*dtftlw)
                  wp=(wp+ath*dtf+bth)*real(icbt)  
                  delz=wp*dtf
                  del_test=(1.-wp)/wp !catch infinity value
!                    if (ieee_is_nan(wp).or.ieee_is_nan(del_test)) then 
!                print*,'PB',wp2,zt2,ust2,wst2,h2,rhoa2,rhograd2,sigw2,dsigwdz2,tlw2,ptot_lhh2,Q_lhh2,phi_lhh2,ath2,bth2,ol2,flagrein,i
!                print*,'PB2',ath,old_wp_buf,bth,wp,sigw
                  if (isnan2(wp).or.isnan2(del_test).or.abs(wp).gt.50.) then 
!                  if (wp.ne.wp .or. del_test.ne.del_test) then
!                  if (ieee_is_nan(wp) .or. ieee_is_nan(del_test).or.abs(wp).gt.50.) then            
                      nrand=nrand+1                      
                      wp=sigw*rannumb(nrand)
                      delz=wp*dtf
                      nan_count(ompid+1)=nan_count(ompid+1)+1
                      print *,'NaN counter equal to:',nan_count(ompid+1),'omp',ompid,'mpi',myid &
            ,'increase ifine if this number became a non-negligible fraction of the particle number'
                  end if  
              end if
            else
                 wp=((1.-dtftlw)*wp+rannumb(nrand+i)*sqrt(2.*dtftlw) &
                 +dtf*(dsigwdz+rhoaux*sigw))*real(icbt) 
                 delz=wp*sigw*dtf
            end if
          else
            rw=exp(-dtftlw)
            wp=(rw*wp+rannumb(nrand+i)*sqrt(1.-rw**2) &
                 +tlw*(1.-rw)*(dsigwdz+rhoaux*sigw))*real(icbt)
            delz=wp*sigw*dtf
          endif
          
        else
          rw=exp(-dtftlw)
          wp=(rw*wp+rannumb(nrand+i)*sqrt(1.-rw**2)*sigw &
               +tlw*(1.-rw)*(dsigw2dz+rhoaux*sigw**2))*real(icbt)
          delz=wp*dtf
        endif

!***************** end turbulent options : comemnt by mc *********************************!

!        if (ieee_is_nan(wp)) then
!         print*,"PROBLEM WP OUT",wp,old_wp_buf,nrand,ol,zt,ust,wst,h,sigw
!        endif
! FLEXPART_WRF - zero up,vp,wp if turbulence is turned off
            if (turb_option .eq. turb_option_none) then
             up=0.0
             vp=0.0
              wp=0.0
              delz=0.
            end if
!          print*,'delz',delz,zt 
!****************************************************
! Compute turbulent vertical displacement of particle
!****************************************************
!     if (nrelpoint.eq.970) then
!      write(15,*),rannumb(nrand),nrand,nrelpoint,OMP_GET_THREAD_NUM()
!     endif

            if (abs(delz).gt.h) delz=mod(delz,h)

! Determine if particle transfers to a "forbidden state" below the ground
! or above the mixing height
!************************************************************************

            if (delz.lt.-zt) then         ! reflection at ground
              icbt=-1
              zt=-zt-delz
               else if (delz.gt.(h-zt)) then ! reflection at h
                 icbt=-1
                 zt=-zt-delz+2.*h
!            else if (delz.gt.(h-zt) .and. reflect_switch==1) then ! reflection at h
!             else if (delz.gt.(h-zt)) then ! reflection at h
!               icbt=-1
!               zt=-zt-delz+2.*h
            else                         ! no reflection
              icbt=1
              zt=zt+delz
            endif

            if (i.ne.ifine) then
! FLEXPART_WRF, TKE option
              if (turb_option .gt. 1) then 
                  do ii=2,nz
                      if (height(ii).gt.zt) then
                       indz=ii-1
                       indzp=ii
                       goto 69
                      endif
                   enddo
69                 continue
 
! If one of the levels necessary is not yet available,
! calculate it
!*****************************************************
 
            do ii=indz,indzp                !i
              if (indzindicator(ii)) then
                if (ngrid.le.0) then
                  call interpol_misslev(ii,xt,yt,zt, &
    uprof,vprof,wprof, usigprof,vsigprof,wsigprof, &
    rhoprof,rhogradprof, tkeprof,pttprof, &
    u,v,w,usig,vsig,wsig,pvi, &
    p1,p2,p3,p4,ddx,ddy,rddx,rddy,dtt,dt1,dt2, &
    ix,jy,ixp,jyp,ngrid,indz,indzp, depoindicator, &
    indzindicator, &
    ust,wst,ol,h,zeta,sigu,sigv,tlu,tlv,tlw, &
    sigw,dsigwdz,dsigw2dz)
                else
                  call interpol_misslev_nests(ii,xt,yt,zt, &
    uprof,vprof,wprof, usigprof,vsigprof,wsigprof, &
    rhoprof,rhogradprof, tkeprof,pttprof, &
    u,v,w,usig,vsig,wsig,pvi, &
    p1,p2,p3,p4,ddx,ddy,rddx,rddy,dtt,dt1,dt2, &
    ix,jy,ixp,jyp,ngrid,indz,indzp, depoindicator, &
    indzindicator, &
    ust,wst,ol,h,zeta,sigu,sigv,tlu,tlv,tlw, &
    sigw,dsigwdz,dsigw2dz)
                endif
              endif
            enddo                          !i
!              write(*,*)'after reflection'
            if(turb_option .eq. turb_option_mytke) &
               call tke_partition_my(zt, &
   ust,wst,ol,h,zeta,sigu,sigv,tlu,tlv,tlw, &
   sigw,dsigwdz,dsigw2dz,uprof,vprof,tkeprof,pttprof,indz,indzp)
            if(turb_option .eq. turb_option_tke)   &
               call tke_partition_hanna(zt, &
   ust,wst,ol,h,zeta,sigu,sigv,tlu,tlv,tlw, &
   sigw,dsigwdz,dsigw2dz,uprof,vprof,tkeprof,pttprof,indz,indzp)
            else
                zeta=zt/h
                 call hanna_short(zt, &
    ust,wst,ol,h,zeta,sigu,sigv,tlu,tlv,tlw, &
    sigw,dsigwdz,dsigw2dz)

            endif
           endif

          enddo
          
          if (cblflag.ne.1) nrand=nrand+i  !------>>>>>>>>>>>>>>>> modified by mc for accounting of different increment of nrand in cbl flag
!    if (nombre.eq.103) print*,'usig 3',usig

! Determine time step for next integration
!*****************************************

          if (turbswitch) then
            ldt=int(min(tlw,h/max(2.*abs(wp*sigw),1.e-5), &
            0.5/abs(dsigwdz))*ctl)
          else
            ldt=int(min(tlw,h/max(2.*abs(wp),1.e-5))*ctl)
          endif
          ldt=max(ldt,mintime)


! If particle represents only a single species, add gravitational settling
! velocity. The settling velocity is zero for gases, or if particle
! represents more than one species
!*************************************************************************

      if (mdomainfill.eq.0) then
        do nsp=1,nspec
!         print*,nrelpoint,nsp
          if (xmass(nrelpoint,nsp).gt.eps2) goto 887
        end do
887     nsp=min(nsp,nspec)
       if (density(nsp).gt.0.) then
!        print*,'settle'
!       print*,'settle 1'
               call get_settling(itime,real(xt),real(yt),zt,nsp,settling)  !bugfix
        endif
        w=w+settling
      endif


! Horizontal displacements during time step dt are small real values compared
! to the position; adding the two, would result in large numerical errors.
! Thus, displacements are accumulated during lsynctime and are added to the
! position at the end
!****************************************************************************

          dxsave=dxsave+u*dt
!     if (nombre.eq.103) print*,'xt-2',dxsave,u,dt
          dysave=dysave+v*dt
          dawsave=dawsave+up*dt
          dcwsave=dcwsave+vp*dt
          zt=zt+w*dt*real(ldirect)  ! comment out and put zt=zt for testing equation based on the well_mixed conditin comemnt by mc 


          if (zt.gt.h) then
            if (itimec.eq.itime+lsynctime) goto 99
            goto 700    ! complete the current interval above PBL
          endif
          if (zt.lt.0.) zt=-1.*zt    ! if particle below ground -> refletion


!!!! CHANGE: TEST OF THE WELL-MIXED CRITERION
!!!! These lines may be switched on to test the well-mixed criterion
!      if (zt.le.h) then
!        zacc=zacc+zt/h*dt
!        hsave=hsave+h*dt
!        tacc=tacc+dt
!        do 67 i=1,iclass
!          if ((zt/h.gt.th(i-1)).and.(zt/h.le.th(i)))
!     +    t(i)=t(i)+dt
!67        continue
!      endif
!c       print*,'itime',itime
!c     if ((mod(abs(itime),3600).eq.0)) then
!c     if ((mod(abs(itime),3600).eq.0).and.dump) then
!      if (itime<itimeold) then
!       print*,'dump well',itime,itimeold
!       dump=.false.
!        itimeold=itimeold-3600
!       write(550,'(i8,12f10.3)') itime,hsave/tacc,zacc/tacc,
!c      write(550,'(i8,22f10.3)') itime,hsave/tacc,zacc/tacc,
!     + (t(i)/tacc*real(iclass),i=1,iclass)
!        flush(550)
!        zacc=0.
!        tacc=0.
!        do 68 i=1,iclass
!68        t(i)=0.
!        hsave=0.
!      endif
!      if (mod(abs(itime),3600).ne.0) dump=.true.
!c       print*,'itime',itime,3600,mod(abs(itime),3600),dump
!!!! CHANGE
!!!******************  NEW test for THE WELL MIXED CRITERION by mc ***************
!$OMP CRITICAL
!if (zt.lt.h) then
!          i_well=int(zt/h*num_lvl*1.)+1                                          !per fare il test qui devo considerare OMP and MPI...
!          well_mixed_vector(i_well,ompid+1)=well_mixed_vector(i_well,ompid+1)+dt
!          well_mixed_norm(ompid+1)=well_mixed_norm(ompid+1)+dt
!          avg_air_dens(i_well,ompid+1)=avg_air_dens(i_well,ompid+1)+rhoa*dt
!          
!      end if
!      h_well(ompid+1)=h
!$OMP END CRITICAL    
!!*********************************************************************************

! Determine probability of deposition
!************************************

      if ((DRYDEP).and.(zt.lt.2.*href)) then
        do ks=1,nspec
          if (DRYDEPSPEC(ks)) then
            if (depoindicator(ks)) then
              if (ngrid.le.0) then
                   call interpol_vdep(ks,vdepo(ks),ix,jy,ixp,jyp, &
                                      p1,p2,p3,p4,dt1,dt2,dtt,depoindicator)
              else
                call interpol_vdep_nests(ks,vdepo(ks),ix,jy,ixp,jyp, &
                                      p1,p2,p3,p4,dt1,dt2,dtt,depoindicator,ngrid)
              endif
            endif
  ! correction by Petra Seibert, 10 April 2001
  !   this formulation means that prob(n) = 1 - f(0)*...*f(n)
  !   where f(n) is the exponential term
               prob(ks)=1.+(prob(ks)-1.)* &
                    exp(-vdepo(ks)*abs(dt)/(2.*href))
          endif
        end do
      endif

      if (zt.lt.0.) zt=min(h-eps2,-1.*zt)    ! if particle below ground -> reflection

      if (itimec.eq.(itime+lsynctime)) then
!   if (nombre.eq.103) print*,'usig',usig,usigprof(indzp)+usigprof(indz),indz
        usig=0.5*(usigprof(indzp)+usigprof(indz))
        vsig=0.5*(vsigprof(indzp)+vsigprof(indz))
        wsig=0.5*(wsigprof(indzp)+wsigprof(indz))
        goto 99  ! finished
      endif
      goto 100

! END TIME LOOP
!==============


      endif



!**********************************************************
! For all particles that are outside the PBL, make a single
! time step. Only horizontal turbulent disturbances are
! calculated. Vertical disturbances are reset.
!**********************************************************


! Interpolate the wind
!*********************

!JB needs to define mu and mv
700   continue
      if (ngrid.le.0) then
        xts=real(xt)
        yts=real(yt)
        call interpol_wind(itime,xts,yts,zt, &
    uprof,vprof,wprof, usigprof,vsigprof,wsigprof, &
    rhoprof,rhogradprof, tkeprof,pttprof, &
    u,v,w,usig,vsig,wsig,pvi, &
    p1,p2,p3,p4,ddx,ddy,rddx,rddy,dtt,dt1,dt2, &
    ix,jy,ixp,jyp,ngrid,indz,indzp, depoindicator, &
    indzindicator,mu,mv)
 !JB mw not needed here
      else
        call interpol_wind_nests(itime,xtn,ytn,zt, &
    uprof,vprof,wprof, usigprof,vsigprof,wsigprof, &
    rhoprof,rhogradprof, tkeprof,pttprof, &
    u,v,w,usig,vsig,wsig,pvi, &
    p1,p2,p3,p4,ddx,ddy,rddx,rddy,dtt,dt1,dt2, &
    ix,jy,ixp,jyp,ngrid,indz,indzp, depoindicator, &
    indzindicator,mu,mv)
      endif

!    if (nombre.eq.103) print*,'usig 4',usig
! Compute everything for above the PBL

! Assume constant, uncorrelated, turbulent perturbations
! In the stratosphere, use a small vertical diffusivity d_strat,
! in the troposphere, use a larger horizontal diffusivity d_trop.
! Turbulent velocity scales are determined based on sqrt(d_trop/dt)
!******************************************************************

      ldt=abs(lsynctime-itimec+itime)
      dt=real(ldt)
      
  if (zt.lt.tropop) then  ! in the troposphere
    uxscale=sqrt(2.*d_trop/dt)
!   if (nrand+1.gt.maxrandomp) nrand=1
    if (nrand+1.gt.maxrand2) nrand=1
    ux=rannumb(nrand)*uxscale
    vy=rannumb(nrand+1)*uxscale
    nrand=nrand+2
    wp=0.
  else if (zt.lt.tropop+1000.) then     ! just above the tropopause: make transition
    weight=(zt-tropop)/1000.
    uxscale=sqrt(2.*d_trop/dt*(1.-weight))
!   if (nrand+2.gt.maxrandomp) nrand=1
    if (nrand+2.gt.maxrand2) nrand=1
    ux=rannumb(nrand)*uxscale
    vy=rannumb(nrand+1)*uxscale
    wpscale=sqrt(2.*d_strat/dt*weight)
    wp=rannumb(nrand+2)*wpscale+d_strat/1000.
    nrand=nrand+3
  else                 ! in the stratosphere
!   if (nrand.gt.maxrandomp) nrand=1
    if (nrand.gt.maxrand2) nrand=1
    ux=0.
    vy=0.
    wpscale=sqrt(2.*d_strat/dt)
    wp=rannumb(nrand)*wpscale
    nrand=nrand+1
  endif

! FLEXPART_WRF - zero ux,vy,wp if turbulence is turned off
      if (turb_option .eq. turb_option_none) then
        ux=0.0
        vy=0.0
        wp=0.0
      end if


! If particle represents only a single species, add gravitational settling
! velocity. The settling velocity is zero for gases
!*************************************************************************

    if (mdomainfill.eq.0) then
      do nsp=1,nspec
        if (xmass(nrelpoint,nsp).gt.eps2) goto 888
      end do
888   nsp=min(nsp,nspec)
    if (density(nsp).gt.0.) then
!       print*,'settle 2, bef',real(xt),real(yt),zt,cpt
           call get_settling(itime,real(xt),real(yt),zt,nsp,settling)  !bugfix
!       print*,'settle 2, aft',real(xt),real(yt),zt,cpt,w
!        print*,'settle'
      endif
      w=w+settling
    endif

! Calculate position at time step itime+lsynctime
!************************************************

!       print*,'settle 2, aft1.5',zt,settling,wp,dt
      dxsave=dxsave+(u+ux)*dt
!     if (nombre.eq.103) print*,'xt-1',dxsave,u,ux,dt
      dysave=dysave+(v+vy)*dt
      zt=zt+(w+wp)*dt*real(ldirect)
!       print*,'settle 2, aft2',zt,cpt
      if (zt.lt.0.) zt=min(h-eps2,-1.*zt)    ! if particle below ground -> reflection
!       print*,'settle 2, aft3',zt,cpt

99    continue



!****************************************************************
! Add mesoscale random disturbances
! This is done only once for the whole lsynctime interval to save
! computation time
!****************************************************************


! Mesoscale wind velocity fluctuations are obtained by scaling
! with the standard deviation of the grid-scale winds surrounding
! the particle location, multiplied by a factor turbmesoscale.
! The autocorrelation time constant is taken as half the
! time interval between wind fields
!****************************************************************

      r=exp(-2.*real(abs(lsynctime))/real(lwindinterv))
      rs=sqrt(1.-r**2)
!     if (nrand+2.gt.maxrandomp) nrand=1
      if (nrand+2.gt.maxrand2) nrand=1
!     if (nombre.eq.103) print*,'usgig0',r,usigold,rannumb(nrand)
      usigold=r*usigold+rs*rannumb(nrand)*usig*turbmesoscale
!     if (nombre.eq.103) print*,'usgig1',usigold,usig,turbmesoscale
      vsigold=r*vsigold+rs*rannumb(nrand+1)*vsig*turbmesoscale
      wsigold=r*wsigold+rs*rannumb(nrand+2)*wsig*turbmesoscale

! FLEXPART_WRF - zero u/v/wsigold if turbulence is turned off
! Note: for mesoscale model applications this component should be ignored!
!     if (turb_option .eq. turb_option_none) then
!       usigold=0.0
!       vsigold=0.0
!       wsigold=0.0
!     end if

      dxsave=dxsave+usigold*real(lsynctime)
      dysave=dysave+vsigold*real(lsynctime)

      zt=zt+wsigold*real(lsynctime)
!       print*,'settle 2, aft4',zt,cpt
      if (zt.lt.0.) zt=-1.*zt    ! if particle below ground -> refletion
!       print*,'settle 2, aft5',zt,cpt

!*************************************************************
! Transform along and cross wind components to xy coordinates,
! add them to u and v, transform u,v to grid units/second
! and calculate new position
!*************************************************************

      call windalign(dxsave,dysave,dawsave,dcwsave,ux,vy)
      dxsave=dxsave+ux
!     if (nombre.eq.103) print*,'xt0',dxsave,usigold,ux
      dysave=dysave+vy
      if (ngrid.ge.0) then
! for FLEXPART_WRF, dx & dy are in meters,
! dxconst=1/dx, dyconst=1/dy, and no cos(lat) is needed
!       cosfact=dxconst/cos((yt*dy+ylat0)*pi180)
!     if (nombre.eq.103) print*,'xt1',xt,dxsave,dxconst
!       xt=xt+real(dxsave*dxconst*real(ldirect),kind=dp)
!       yt=yt+real(dysave*dyconst*real(ldirect),kind=dp)
!      xt=xt+real(dxsave/mu*dxconst*real(ldirect),kind=dp)
!      yt=yt+real(dysave/mv*dyconst*real(ldirect),kind=dp)
       xt=xt +real(dxsave*mu*dxconst*real(ldirect),kind=dp)  !IF COOMMENTED OUT TO is to ISOLate VERTCAL FORMULAITON FOR TEST REASON BY mc
       yt=yt +real(dysave*mv*dyconst*real(ldirect),kind=dp)  !IF COOMMENTED OUT TO is to ISOLate VERTCAL FORMULAITON FOR TEST REASON BY mc
! JB: needs interpolate m_w on the coordinates
!     else if (ngrid.eq.-1) then      ! around north pole
!       xlon=xlon0+xt*dx
!       ylat=ylat0+yt*dy
!       call cll2xy(northpolemap,ylat,xlon,xpol,ypol)
!       gridsize=1000.*cgszll(northpolemap,ylat,xlon)
!       dxsave=dxsave/gridsize
!       dysave=dysave/gridsize
!       xpol=xpol+dxsave*real(ldirect)
!       ypol=ypol+dysave*real(ldirect)
!       call cxy2ll(northpolemap,xpol,ypol,ylat,xlon)
!       xt=(xlon-xlon0)/dx
!       yt=(ylat-ylat0)/dy
!     else if (ngrid.eq.-2) then    ! around south pole
!       xlon=xlon0+xt*dx
!       ylat=ylat0+yt*dy
!       call cll2xy(southpolemap,ylat,xlon,xpol,ypol)
!       gridsize=1000.*cgszll(southpolemap,ylat,xlon)
!       dxsave=dxsave/gridsize
!       dysave=dysave/gridsize
!       xpol=xpol+dxsave*real(ldirect)
!       ypol=ypol+dysave*real(ldirect)
!       call cxy2ll(southpolemap,xpol,ypol,ylat,xlon)
!       xt=(xlon-xlon0)/dx
!       yt=(ylat-ylat0)/dy
      else
        write(*,*) 'advance -- bad ngrid = ', ngrid
        stop
      endif


! If global data are available, use cyclic boundary condition
!************************************************************

      if (xglobal) then
        if (xt.ge.real(nxmin1)) xt=xt-real(nxmin1)
        if (xt.lt.0.) xt=xt+real(nxmin1)
        if (xt.le.eps) xt=eps
        if (abs(xt-real(nxmin1)).le.eps) xt=real(nxmin1)-eps
      endif


! Check position: If trajectory outside model domain, terminate it
!*****************************************************************

      if ((xt.lt.0.).or.(xt.ge.real(nxmin1)).or.(yt.lt.0.).or. &
      (yt.ge.real(nymin1))) then
        nstop=3
        return
      endif

! If particle above highest model level, set it back into the domain
!*******************************************************************

      if (zt.ge.height(nz)) zt=height(nz)-100.*eps
!       print*,'settle 2, aft6',zt,cpt


!************************************************************************
! Now we could finish, as this was done in FLEXPART versions up to 4.0.
! However, truncation errors of the advection can be significantly
! reduced by doing one iteration of the Petterssen scheme, if this is
! possible.
! Note that this is applied only to the grid-scale winds, not to
! the turbulent winds.
!************************************************************************

! The Petterssen scheme can only applied with long time steps (only then u
! is the "old" wind as required by the scheme); otherwise do nothing
!*************************************************************************

      if (ldt.ne.abs(lsynctime)) return

! The Petterssen scheme can only be applied if the ending time of the time step
! (itime+ldt*ldirect) is still between the two wind fields held in memory;
! otherwise do nothing
!******************************************************************************

      if (abs(itime+ldt*ldirect).gt.abs(memtime(2))) return  

! Apply it also only if starting and ending point of current time step are on
! the same grid; otherwise do nothing
!*****************************************************************************
      if (nglobal.and.(yt.gt.switchnorthg)) then
        ngr=-1
      else if (sglobal.and.(yt.lt.switchsouthg)) then
        ngr=-2
      else
    ngr=0
    do j=numbnests,1,-1
      if ((xt.gt.xln(j)+eps).and.(xt.lt.xrn(j)-eps).and. &
           (yt.gt.yln(j)+eps).and.(yt.lt.yrn(j)-eps)) then
        ngr=j
        goto 43
      endif
    end do
43   continue
  endif

      if (ngr.ne.ngrid) return

! Determine nested grid coordinates
!**********************************

      if (ngrid.gt.0) then
        xtn=(xt-xln(ngrid))*xresoln(ngrid)
        ytn=(yt-yln(ngrid))*yresoln(ngrid)
        ix=int(xtn)
        jy=int(ytn)
      else
        ix=int(xt)
        jy=int(yt)
      endif 
      ixp=ix+1
      jyp=jy+1


! Memorize the old wind
!**********************

      uold=u
      vold=v
      wold=w

! Interpolate wind at new position and time
!******************************************

      if (ngrid.le.0) then
        xts=real(xt)
        yts=real(yt)
        call interpol_wind_short(itime+ldt*ldirect,xts,yts,zt, &
    uprof,vprof,wprof, usigprof,vsigprof,wsigprof, &
    rhoprof,rhogradprof, tkeprof,pttprof, &
    u,v,w,usig,vsig,wsig,pvi, &
    p1,p2,p3,p4,ddx,ddy,rddx,rddy,dtt,dt1,dt2, &
    ix,jy,ixp,jyp,ngrid,indz,indzp, depoindicator, &
    indzindicator)
      else
        call interpol_wind_short_nests(itime+ldt*ldirect,xtn,ytn,zt, &
    uprof,vprof,wprof, usigprof,vsigprof,wsigprof, &
    rhoprof,rhogradprof, tkeprof,pttprof, &
    u,v,w,usig,vsig,wsig,pvi, &
    p1,p2,p3,p4,ddx,ddy,rddx,rddy,dtt,dt1,dt2, &
    ix,jy,ixp,jyp,ngrid,indz,indzp, depoindicator, &
    indzindicator)

      endif
!       print*,'settle 2, aft7',zt,cpt

  if (mdomainfill.eq.0) then
    do nsp=1,nspec
      if (xmass(nrelpoint,nsp).gt.eps2) goto 889
    end do
889   nsp=min(nsp,nspec)
   if (density(nsp).gt.0.) then
!       print*,'settle 3, bef',real(xt),real(yt),zt,cpt
        call get_settling(itime+ldt,real(xt),real(yt),zt,nsp,settling)  !bugfix
!       print*,'settle 3, aft',real(xt),real(yt),zt,cpt
!        print*,'settle'
      endif
    w=w+settling
  endif

! Determine the difference vector between new and old wind
! (use half of it to correct position according to Petterssen)
!*************************************************************

      u=(u-uold)/2.
      v=(v-vold)/2.
      w=(w-wold)/2.


! Finally, correct the old position
!**********************************

      zt=zt+w*real(ldt*ldirect)
  if (zt.lt.0.) zt=min(h-eps2,-1.*zt)    ! if particle below ground -> reflection

      if (ngrid.ge.0) then
! for FLEXPART_WRF, dx & dy are in meters,
! dxconst=1/dx, dyconst=1/dy, and no cos(lat) is needed
!       cosfact=dxconst/cos((yt*dy+ylat0)*pi180)
!       if (nombre.eq.103)  print*,'xt2',xt,u,dxconst,ldt
!       xt=xt+real(u*dxconst*real(ldt*ldirect),kind=dp)
!       yt=yt+real(v*dyconst*real(ldt*ldirect),kind=dp)
!        print*,'mw',mu,mv
!        xt=xt+real(u*dxconst/mu*real(ldt*ldirect),kind=dp)
!        yt=yt+real(v*dyconst/mv*real(ldt*ldirect),kind=dp)
         xt=xt +real(u*dxconst*mu*real(ldt*ldirect),kind=dp)  !IF COOMMENTED OUT TO is to ISOLate VERTCAL FORMULAITON FOR TEST REASON BY mc
         yt=yt +real(v*dyconst*mv*real(ldt*ldirect),kind=dp)  !IF COOMMENTED OUT TO is to ISOLate VERTCAL FORMULAITON FOR TEST REASON BY mc
!     else if (ngrid.eq.-1) then      ! around north pole
!       xlon=xlon0+xt*dx
!       ylat=ylat0+yt*dy
!       call cll2xy(northpolemap,ylat,xlon,xpol,ypol)
!       gridsize=1000.*cgszll(northpolemap,ylat,xlon)
!       u=u/gridsize
!       v=v/gridsize
!       xpol=xpol+u*real(ldt*ldirect)
!       ypol=ypol+v*real(ldt*ldirect)
!       call cxy2ll(northpolemap,xpol,ypol,ylat,xlon)
!       xt=(xlon-xlon0)/dx
!       yt=(ylat-ylat0)/dy
!     else if (ngrid.eq.-2) then    ! around south pole
!       xlon=xlon0+xt*dx
!       ylat=ylat0+yt*dy
!       call cll2xy(southpolemap,ylat,xlon,xpol,ypol)
!       gridsize=1000.*cgszll(southpolemap,ylat,xlon)
!       u=u/gridsize
!       v=v/gridsize
!       xpol=xpol+u*real(ldt*ldirect)
!       ypol=ypol+v*real(ldt*ldirect)
!       call cxy2ll(southpolemap,xpol,ypol,ylat,xlon)
!       xt=(xlon-xlon0)/dx
!       yt=(ylat-ylat0)/dy
      else
        write(*,*) 'advance -- bad ngrid = ', ngrid
        stop
      endif

! If global data are available, use cyclic boundary condition
!************************************************************

      if (xglobal) then
        if (xt.ge.real(nxmin1)) xt=xt-real(nxmin1)
        if (xt.lt.0.) xt=xt+real(nxmin1)
        if (xt.le.eps) xt=eps
        if (abs(xt-real(nxmin1)).le.eps) xt=real(nxmin1)-eps
      endif

! Check position: If trajectory outside model domain, terminate it
!*****************************************************************

      if ((xt.lt.0.).or.(xt.ge.real(nxmin1)).or.(yt.lt.0.).or. &
      (yt.ge.real(nymin1))) then
        nstop=3
        return
      endif

! If particle above highest model level, set it back into the domain
!*******************************************************************

      if (zt.ge.height(nz)) zt=height(nz)-100.*eps

!       if (nombre.eq.103)  print*,'end',xt,u,dxconst,ldt

      end subroutine advance

!       logical function isnan2(a)
!       real :: a
!        if ((a.ne.a)) then !.or.((a*0.).ne.0.)) then
!       isnan2 = .true.
!       else
!       isnan2 = .false.
!        end if
!        return
!       end 
