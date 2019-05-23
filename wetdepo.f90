!***********************************************************************
!* Copyright 2012,2013                                                *
!* Jerome Brioude, Delia Arnold, Andreas Stohl, Wayne Angevine,       *
!* John Burkhart, Massimo Cassiani, Adam Dingwell, Richard C Easter, Sabine Eckhardt,*
!* Stephanie Evan, Jerome D Fast, Don Morton, Ignacio Pisso,          *
!* Petra Seibert, Gerard Wotawa, Caroline Forster, Harald Sodemann,   *
!*                                                                     *
!* This file is part of FLEXPART WRF            
!                                                                     *
! FLEXPART is free software: you can redistribute it and/or modify    *
! it under the terms of the GNU General Public License as published by*
! the Free Software Foundation, either version 3 of the License, or   *
! (at your option) any later version.                                 *
!                                                                     *
! FLEXPART is distributed in the hope that it will be useful,         *
! but WITHOUT ANY WARRANTY; without even the implied warranty of      *
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       *
! GNU General Public License for more details.                        *
!                                                                     *
! You should have received a copy of the GNU General Public License   *
! along with FLEXPART.  If not, see <http://www.gnu.org/licenses/>.   *
!**********************************************************************

subroutine wetdepo(itime,ltsample,loutnext)
  !                     i      i        i
  !*****************************************************************************
  !                                                                            *
  ! Calculation of wet deposition using the concept of scavenging coefficients.*
  ! For lack of detailed information, washout and rainout are jointly treated. *
  ! It is assumed that precipitation does not occur uniformly within the whole *
  ! grid cell, but that only a fraction of the grid cell experiences rainfall. *
  ! This fraction is parameterized from total cloud cover and rates of large   *
  ! scale and convective precipitation.                                        *
  !                                                                            *
  !    Author: A. Stohl                                                        *
  !                                                                            *
  !    1 December 1996                                                         *
  !                                                                            *
  ! Correction by Petra Seibert, Sept 2002:                                    *
  ! use centred precipitation data for integration                             *
  ! Code may not be correct for decay of deposition!                           *
  !                                                                            *
  !* D Arnold 2012 01 25 implement patch developed by P.Seibert to avoid       *
  !* problems with cloud fraction                                              *
  !* set cloud cover to 1  when there is precip & use diagnosed clouds in      *
  !* vertransform* for in-cloud below-cloud scavenging                         *
  !* Set fraction to max val from grid cells larger than 5 km.                 *
  !* Missing: new fractions of grid cell affected by rain , coming from        *
  !* modelling studies with 150 km grid-sizes!                                 *
  !* precip interpolation is under check by developers                         *
  !* modifications identified by "CDA"                                         *
!*
!* Petra Seibert, 2011/2012: Fixing some deficiencies in this modification
  !********************************************************************************
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! cc [0-1]           total cloud cover                                       *
  ! convp [mm/h]       convective precipitation rate                           *
  ! grfraction [0-1]   fraction of grid, for which precipitation occurs        *
  ! ix,jy              indices of output grid cell for each particle           *
  ! itime [s]          actual simulation time [s]                              *
  ! jpart              particle index                                          *
  ! ldeltat [s]        interval since radioactive decay was computed           *
  ! lfr, cfr           area fraction covered by precipitation for large scale  *
  !                    and convective precipitation (dependent on prec. rate)  *
  ! loutnext [s]       time for which gridded deposition is next output        *
  ! loutstep [s]       interval at which gridded deposition is output          *
  ! lsp [mm/h]         large scale precipitation rate                          *
  ! ltsample [s]       interval over which mass is deposited                   *
  ! prec [mm/h]        precipitation rate in subgrid, where precipitation occurs*
  ! wetdeposit         mass that is wet deposited                              *
  ! wetgrid            accumulated deposited mass on output grid               *
  ! wetscav            scavenging coefficient                                  *
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************

  use point_mod
  use par_mod
  use com_mod
! use ieee_arithmetic
  implicit none

  integer :: jpart,itime,ltsample,loutnext,ldeltat,i,j,ix,jy
!  integer :: ngrid,itage,nage,hz,il,interp_time, n, clouds_v
  integer :: ngrid,itage,nage,kz,il,interp_time, n
!  integer :: ks, kp
  integer :: ks, kp, icbot,ictop, indcloud

  real :: S_i, act_temp, cl, cle ! in cloud scavenging
!  real :: clouds_h ! cloud height for the specific grid point

!  real :: xtn,ytn,lsp,convp,cc,grfraction,prec,wetscav
  real :: xtn,ytn,lsp,convp,cc,fraction,prec,wetscav,precsub,f
  real :: wetscavold
!  real :: wetdeposit(maxspec),restmass
  real :: wetdeposit(maxspec),restmass


  real,parameter :: smallnum = tiny(0.0) ! smallest number that can be handled
  save lfr,cfr


  real :: lfr(5) = (/ 0.5,0.65,0.8,0.9,0.95/)
  real :: cfr(5) = (/ 0.4,0.55,0.7,0.8,0.9 /)

  ! Compute interval since radioactive decay of deposited mass was computed
  !************************************************************************

  if (itime.le.loutnext) then
    ldeltat=itime-(loutnext-loutstep)
  else                                  ! first half of next interval
    ldeltat=itime-loutnext
  endif


  ! Loop over all particles
  !************************

  do jpart=1,numpart
    if (itra1(jpart).eq.-999999999) goto 20
    if(ldirect.eq.1)then
      if (itra1(jpart).gt.itime) goto 20
    else
      if (itra1(jpart).lt.itime) goto 20
    endif
  ! Determine age class of the particle
    itage=abs(itra1(jpart)-itramem(jpart))
    do nage=1,nageclass
      if (itage.lt.lage(nage)) goto 33
    end do
33   continue


  ! Determine which nesting level to be used
  !*****************************************

    ngrid=0
    do j=numbnests,1,-1
      if ((xtra1(jpart).gt.xln(j)).and.(xtra1(jpart).lt.xrn(j)).and. &
           (ytra1(jpart).gt.yln(j)).and.(ytra1(jpart).lt.yrn(j))) then
        ngrid=j
        goto 23
      endif
    end do
23   continue


  ! Determine nested grid coordinates
  !**********************************

    if (ngrid.gt.0) then
      xtn=(xtra1(jpart)-xln(ngrid))*xresoln(ngrid)
      ytn=(ytra1(jpart)-yln(ngrid))*yresoln(ngrid)
      ix=int(xtn)
      jy=int(ytn)
    else
      ix=int(xtra1(jpart))
      jy=int(ytra1(jpart))
    endif


  ! Interpolate large scale precipitation, convective precipitation and
  ! total cloud cover
  ! Note that interpolated time refers to itime-0.5*ltsample [PS]
  !********************************************************************
    interp_time=nint(itime-0.5*ltsample)
!   if (jpart.eq.103) print*,real(xtra1(jpart)),real(ytra1(jpart)),jpart

! CDA part of the bug fix from P. Seiber includes modification of
! CDA interpol_rain routines.
!    if (ieee_is_nan(xtra1(jpart))) then
      if (xtra1(jpart).ne.xtra1(jpart)) then
     print*,jpart,xtra1(jpart),ytra1(jpart),itra1(jpart)
       endif
    if (ngrid.eq.0) then
!CDA old code:
!      call interpol_rain(lsprec,convprec,tcc,nxmax,nymax, &
!           1,nx,ny,memind,real(xtra1(jpart)),real(ytra1(jpart)),1, &
!           memtime(1),memtime(2),interp_time,lsp,convp,cc)
!CDA new code:
        
          call interpol_rain(lsprec,convprec,tcc, &
           icloudbot,icloudthck,nxmax,nymax,1,nx,ny, &
!          memind,sngl(xtra1(jpart)),sngl(ytra1(jpart)),1,memtime(1), &
           memind,real(xtra1(jpart)),real(ytra1(jpart)),1,memtime(1), &
           memtime(2),interp_time,lsp,convp,cc,icbot,ictop,icmv)

    else
! CDA old code:
!      call interpol_rain_nests(lsprecn,convprecn,tccn, &
!           nxmaxn,nymaxn,1,maxnests,ngrid,nxn,nyn,memind,xtn,ytn,1, &
!           memtime(1),memtime(2),interp_time,lsp,convp,cc)
! CDA new code:
      call interpol_rain_nests(lsprecn,convprecn,tccn, &
           icloudbotn,icloudthckn,nxmaxn,nymaxn,1,&
             maxnests,ngrid,nxn,nyn, memind,xtn,ytn, &
           1,memtime(1), &
           memtime(2),interp_time,lsp,convp,cc,icbot,ictop,icmv)

    endif

!    if ((lsp.lt.0.01).and.(convp.lt.0.01)) goto 20

!CPS 2012: subtract a small value, eg 0.01 mm/h, to remove spurious precip
        prec = lsp+convp
        precsub = 0.01
        if (prec .lt. precsub) then
          goto 20
        else
          f = (prec-precsub)/prec
          lsp = f*lsp
          convp = f*convp
        endif


  ! get the level were the actual particle is in
      do il=2,nz
        if (height(il).gt.ztra1(jpart)) then
          kz=il-1
          goto 26
        endif
      end do
26     continue

  n=memind(2)
  if (abs(memtime(1)-interp_time).lt.abs(memtime(2)-interp_time)) &
       n=memind(1)

  ! if there is no precipitation or the particle is above the clouds no
  ! scavenging is done
! CDA  PS part of fix

        if (ztra1(jpart) .le. float(ictop)) then
          if (ztra1(jpart) .gt. float(icbot)) then
            indcloud = 2 ! in-cloud
          else
            indcloud = 1 ! below-cloud
          endif
        elseif (ictop .eq. icmv) then
          indcloud = 0 ! no cloud found, use old scheme
        else
          goto 20 ! above cloud
        endif

!  if (ngrid.eq.0) then
!     clouds_v=clouds(ix,jy,hz,n)
!     clouds_h=cloudsh(ix,jy,n)
!  else
!     clouds_v=cloudsn(ix,jy,hz,n,ngrid)
!     clouds_h=cloudsnh(ix,jy,n,ngrid)
!  endif
  !write(*,*) 'there is
  !    + precipitation',(clouds(ix,jy,ihz,n),ihz=1,20),lsp,convp,hz
!  if (clouds_v.le.1) goto 20
  !write (*,*) 'there is scavenging'

  ! 1) Parameterization of the the area fraction of the grid cell where the
  !    precipitation occurs: the absolute limit is the total cloud cover, but
  !    for low precipitation rates, an even smaller fraction of the grid cell
  !    is used. Large scale precipitation occurs over larger areas than
  !    convective precipitation.
  !**************************************************************************

    if (lsp.gt.20.) then
      i=5
    else if (lsp.gt.8.) then
      i=4
    else if (lsp.gt.3.) then
      i=3
    else if (lsp.gt.1.) then
      i=2
    else
      i=1
    endif

    if (convp.gt.20.) then
      j=5
    else if (convp.gt.8.) then
      j=4
    else if (convp.gt.3.) then
      j=3
    else if (convp.gt.1.) then
      j=2
    else
      j=1
    endif

  !CDA
     !  if (dx.le.5000.) then
     !    j=5 
     !    i=5 
     !  endif
   ! for the moment i set all to 5
       i=5
       j=5
  !CDA
  !CDA cc (coming from CLDFRA works as a mask for the precipitation
  ! in vertransform* we diagnose the clouds. We consider not that
  ! when there is precip, there is cloud as diagnosed in vertransform
  ! this is the same approx used in other models HYSPLIT, for insance
  ! set cc to 1
    cc=1. 

    fraction=max(0.05,cc*(lsp*lfr(i)+convp*cfr(j))/(lsp+convp))

  ! 2) Computation of precipitation rate in sub-grid cell
  !******************************************************

    prec=(lsp+convp)/fraction

  ! 3) Computation of scavenging coefficients for all species
  !    Computation of wet deposition
  !**********************************************************

    do 10 ks=1,nspec                                  ! loop over species
      wetdeposit(ks)=0.

          if (weta(ks).gt.0.) then
            if (indcloud .eq. 1) then ! BELOW CLOUD SCAVENGING
!C               for aerosols and not highliy soluble substances weta=5E-6
              wetscav=weta(ks)*prec**wetb(ks)                ! scavenging coeff
!c             write(*,*) 'bel. wetscav: ',wetscav
            elseif (indcloud .eq. 2) then !  IN CLOUD SCAVENGING
              if (ngrid.gt.0) then
                 act_temp=ttn(ix,jy,kz,n,ngrid)
              else
                 act_temp=tt(ix,jy,kz,n)
              endif
              cl=2E-7*prec**0.36
              if (dquer(ks).gt.0) then ! is particle
                S_i=0.9/cl
              else ! is gas
                cle=(1-cl)/(henry(ks)*(r_air/3500.)*act_temp)+cl
                S_i=1/cle
              endif
              wetscav=S_i*prec/3.6E6/(ictop-icbot) ! 3.6e6 converts mm/h to m/s
              wetscavold = 2.e-5*prec**0.8 ! 1.e-4*prec**0.62
              wetscav = min(0.1*wetscav,wetscavold) ! 0.1 is ASt current setting
!PS tihb new scheme
!PS tihc for no cloud, 1.e-5*prec**0.8 ! 1.e-4*prec**0.62
!PS tihd limit to wetscavold= 1.e-4*prec**0.8
!PS tihe mulitiply wetscav with .1
!PS tihf set the wetscavold to 1.e-5prec**0.8
!PS tihg set the wetscav & wetscavold to 2.e-5*prec**0.8 + weta=2.e-6 (instead of 1.e-6)

            else ! PS: no cloud diagnosed, old scheme,
!CPS          using with fixed a,b for simplicity, one may wish to change!!
!             wetscav = 1.e-4*prec**0.62
              wetscav = 2.e-5*prec**0.8 
            endif

            wetdeposit(ks)=xmass1(jpart,ks)*  &
              (1.-exp(-wetscav*abs(ltsample)))*fraction  ! wet deposition
! CDA test
!        if (wetdeposit(ks).gt.0) then
!           write(*,*) 'wetdepo: ',wetdeposit(ks),ks
!        endif


            restmass = xmass1(jpart,ks)-wetdeposit(ks)
            if (ioutputforeachrelease.eq.1) then
              kp=npoint(jpart)
            else
              kp=1
            endif
            if (restmass .gt. smallnum) then
              xmass1(jpart,ks)=restmass
!cccccccccccccccc depostatistic
!c            wetdepo_sum(ks,kp)=wetdepo_sum(ks,kp)+wetdeposit(ks)
!cccccccccccccccc depostatistic
            else
              xmass1(jpart,ks)=0.
            endif
!C Correct deposited mass to the last time step when radioactive decay of
!C gridded deposited mass was calculated
            if (decay(ks).gt.0.) then
              wetdeposit(ks)=wetdeposit(ks)*exp(abs(ldeltat)*decay(ks))
            endif
          else  ! weta(k)<0
             wetdeposit(ks)=0.
          endif
10      continue



  ! Sabine Eckhard, June 2008 create deposition runs only for forward runs
  ! Add the wet deposition to accumulated amount on output grid and nested output grid
  !*****************************************************************************

    if (ldirect.eq.1) then
!   print*,'kp',kp,jpart,npoint(jpart)
  ! CDA added itage for the kernel not to be applied during the first hours
    call wetdepokernel(nclass(jpart),wetdeposit,real(xtra1(jpart)), &
         real(ytra1(jpart)),itage,nage,kp)
    if (nested_output.eq.1) call wetdepokernel_nest(nclass(jpart), &
         wetdeposit,real(xtra1(jpart)),real(ytra1(jpart)),itage, &
         nage,kp)
    endif

20  continue
  end do

end subroutine wetdepo
