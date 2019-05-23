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

       subroutine timemanager(mts)

!*******************************************************************************
!                                                                              *
! Handles the computation of trajectories, i.e. determines which               *
! trajectories have to be computed at what time.                               *
! Manages dry+wet deposition routines, radioactive decay and the computation   *
! of concentrations.                                                           *
!                                                                              *
!     Author: A. Stohl                                                         *
!                                                                              *
!     20 May 1996                                                              *
!                                                                              *
!     Dec 2005, J. Fast - Only call conccalc & concoutput when (iout.ge.1)     *
!     Aug 2007, W. Wang - call KFeta convection scheme (lconvection=2or3)
!                       Note, backward is unavailabe for lconvection=2
!     Mar 2012, J. Brioude: modifications to handle openmp and mpi             *
!     Jun 2015, J. Brioude: bug fix on xmassfract when ipout=1
!*******************************************************************************
!  Changes, Bernd C. Krueger, Feb. 2001:                                       *
!        Call of convmix when new windfield is read                            *
!------------------------------------                                          *
!  Changes Petra Seibert, Sept 2002                                            *
!     fix wet scavenging problem                                               *
!     Code may not be correct for decay of deposition!                         *
!  Changes Petra Seibert, Nov 2002                                             *
!     call convection BEFORE new fields are read in BWD mode                   *
!  Changes Caroline Forster, Feb 2005
!     new interface between flexpart and convection scheme
!     Emanuel's latest subroutine convect43c.f is used
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
! DEP                .true. if either wet or dry deposition is switched on     *
! decay(maxspec) [1/s] decay constant for radioactive decay                    *
! DRYDEP             .true. if dry deposition is switched on                   *
! ideltas [s]        modelling period                                          *
! jtime [s]          actual temporal position of calculation                   *
! ldeltat [s]        time since computation of radioact. decay of depositions  *
! loutaver [s]       averaging period for concentration calculations           *
! loutend [s]        end of averaging for concentration calculations           *
! loutnext [s]       next time at which output fields shall be centered        *
! loutsample [s]     sampling interval for averaging of concentrations         *
! loutstart [s]      start of averaging for concentration calculations         *
! loutstep [s]       time interval for which concentrations shall be calculated*
! npoint(maxpart)    index, which starting point the trajectory has            *
!                    starting positions of trajectories                        *
! nstop              serves as indicator for fate of particles                 *
!                    in the particle loop                                      *
! nstop1             serves as indicator for wind fields (see getfields)       *
! outnum             number of samples for each concentration calculation      *
! outnum             number of samples for each concentration calculation      *
! prob               probability of absorption at ground due to dry deposition *
! WETDEP             .true. if wet deposition is switched on                   *
! weight             weight for each concentration sample (1/2 or 1)           *
! uap(maxpart),ucp(maxpart),uzp(maxpart) = random velocities due to turbulence *
! us(maxpart),vs(maxpart),ws(maxpart) = random velocities due to interpolation *
! xtra1(maxpart), ytra1(maxpart), ztra1(maxpart) =                             *
!                    spatial positions of trajectories                         *
!                                                                              *
! Constants:                                                                   *
! maxpart            maximum number of trajectories                            *
!                                                                              *
!*******************************************************************************

!      include 'includepar'
!      include 'includecom'
  use unc_mod
  use point_mod
!  use xmass_mod
  use flux_mod
  use outg_mod
  use oh_mod
  use par_mod
  use com_mod
  use mt_stream
  use netcdf_output_mod

!  use ran_mod
!  use interpol_mod

     implicit none
  !include 'netcdf.inc'


  integer :: ix,jy,j,ks,kp,l,n,jtime,nstop,nstop1
! integer :: ksp
  integer :: loutnext,loutstart,loutend,jj,chunksize,ncret
!,chunksize2
  integer :: chunksize3,omp_get_num_threads
  integer :: ldeltat,itage,nage,i
  integer :: clck_counts_beg, clck_counts_end, clck_rate
  integer :: clck_counts_beg2, clck_counts_end2, clck_rate2
  real :: tins,tins2
  real :: outnum,weight,prob(maxspec),nrand,decfact
!  real :: uap(maxpart),ucp(maxpart),uzp(maxpart)
!  real :: us(maxpart),vs(maxpart),ws(maxpart)
!  integer(kind=2) :: cbt(maxpart)
!  real,allocatable, dimension (:) :: uap,ucp,uzp
!  real,allocatable, dimension (:) :: us,vs,ws
!  integer(kind=2),allocatable, dimension (:) :: cbt
  real :: drydeposit(maxspec),gridtotalunc,wetgridtotalunc
  real :: drygridtotalunc,xold,yold,zold,xmassfract
!      integer j,k,l,n,jtime,nstop,nstop1
!      integer loutnext,loutstart,loutend
!      integer ix,jy,ldeltat,itage,nage
!      real outnum,weight,prob(maxspec)
!     real uap(maxpart),ucp(maxpart),uzp(maxpart),decfact
!     real us(maxpart),vs(maxpart),ws(maxpart),cbt(maxpart)
!     real drydeposit(maxspec),gridtotalunc,wetgridtotalunc
!      real drygridtotalunc,xold,yold,zold
!     real xm,xm1


  integer :: jdeb,jfin,stat
  integer,save :: cpt(maxomp)=0
! integer,save :: cpt(24)=0
  real :: ran3
  integer :: from,OMP_GET_THREAD_NUM
  real :: p1,p2,p3,p4,ddx,ddy,rddx,rddy,dtt,dt1,dt2
  integer :: ixp,jyp,ngrid,indz,indzp,nbp,jj2,ii,offset
  logical :: depoindicator(maxspec)
  logical,save :: indzindicator(nzmax)
  real :: ust,wst,ol,h,zeta,sigu,sigv,tlu,tlv,tlw
  real :: sigw,dsigwdz,dsigw2dz
  real :: start, finish
  real :: start2, finish2
  real :: uprof(nzmax),vprof(nzmax),wprof(nzmax)
  real :: usigprof(nzmax),vsigprof(nzmax),wsigprof(nzmax)
  real :: rhoprof(nzmax),rhogradprof(nzmax)
  real :: tkeprof(nzmax),pttprof(nzmax)
  real :: u,v,w,usig,vsig,wsig,pvi
  integer*4 :: now(3)
  integer :: ttime,cpttra
!  integer, dimension(MPI_STATUS_SIZE) :: status
  integer :: myid,ntasks,ierr,islave,tag2,ompid,n_threads,tag3,i_omp
  type (mt_state) :: mts (0: MAX_STREAM)
!************************

!JB
!  call MPI_COMM_RANK ( MPI_COMM_WORLD, myid, ierr )
!  call MPI_COMM_SIZE ( MPI_COMM_WORLD, ntasks, ierr )
! myid gives the info on the node id
      ntasks=1
      myid=0
      loutnext=loutstep/2
      outnum=0.
      loutstart=loutnext-loutaver/2
      loutend=loutnext+loutaver/2

!   if (myid.eq.0) then
    allocate(uap(maxpart) ,stat=stat)
    allocate(ucp(maxpart) ,stat=stat)
    allocate(uzp(maxpart) ,stat=stat)
    allocate(us(maxpart) ,stat=stat)
    allocate(vs(maxpart) ,stat=stat)
    allocate(ws(maxpart) ,stat=stat)
    allocate(cbt(maxpart) ,stat=stat)
!   endif
  call system_clock (clck_counts_beg2,clck_rate2)

!**********************************************************************
! Loop over the whole modelling period in time steps of mintime seconds
!**********************************************************************

!     print*,'time',myid,ideltas,lsynctime
      do jtime=0,ideltas,lsynctime


!         print*,'jtime',jtime
! Computation of wet deposition every lsynctime seconds
! maybe wet depo frequency can be relaxed later but better be on safe side
! wetdepo must be called BEFORE new fields are read in but should not
! be called in the very beginning before any fields are loaded, or
! before particles are in the system
! Code may not be correct for decay of deposition
! changed by Petra Seibert 9/02
!********************************************************************

        if (WETDEP .and. jtime .ne. 0 .and. numpart .gt. 0) &
          call wetdepo(jtime,lsynctime,loutnext)

    if (OHREA .and. jtime .ne. 0 .and. numpart .gt. 0) &
         call ohreaction(jtime,lsynctime,loutnext)

! compute convection for backward runs
!*************************************

!          if ((ldirect.eq.-1).and.(lconvection.eq.1).and.(jtime.lt.0))
!    &    call convmix(jtime)

           if ((ldirect.eq.-1).and.(jtime.lt.0)) then 
             if (lconvection .eq. 1) call convmix(jtime)
             if (lconvection .eq. 2 .or. lconvection .eq. 3) &
                call convmix_kfeta(jtime)
           endif

! Get necessary wind fields if not available
!*******************************************

!        call itime(now)
!        ttime=now(1)*3600+now(2)*60+now(3)
!  call system_clock (clck_counts_beg,clck_rate)
        call getfields(jtime,nstop1)
!  call system_clock (clck_counts_end,clck_rate)
!  tins2=real(clck_counts_end - clck_counts_beg)/real(clck_rate)
!              print*,'time getfield',tins2

        if (nstop1.gt.1) stop 'NO METEO FIELDS AVAILABLE'
!        call itime(now)
!        ttime=now(1)*3600+now(2)*60+now(3)-ttime
!      print*,'read wind time',ttime

! Release particles
!******************

!JB
     if (myid.eq.0) then ! I let only the master thread releasing the particles and calculate the output
!        call itime(now)
        if (mdomainfill.ge.1) then
          if (jtime.eq.0) then
            call init_domainfill()
          else
            call boundcond_domainfill(jtime,loutend)
          endif
        else
        if (numpoint_option.eq.0) then
          call releaseparticles_irreg(jtime)
         elseif (numpoint_option.eq.1) then
!      print*,'avant release'
          call releaseparticles_reg(jtime) 
          endif
        endif
!           do i=1,numpart
!         print*,'ipart 2',myid,i,ztra1(i)
!            enddo
!        print*,'test rel',npoint(1),npoint(2),npoint(3)

!         print*,'test rel1',npoint(5139),npoint(6002),npoint(100003)
! Compute convective mixing for forward runs
! for backward runs it is done before next windfield is read in
!**************************************************************

!      if ((ldirect.eq.1).and.(lconvection.eq.1)) &
!           call convmix(jtime)

          if (ldirect.eq.1) then 
           if (lconvection.eq.1)call convmix(jtime)
           if (lconvection.eq.2 .or. lconvection .eq. 3) &
             call convmix_kfeta(jtime)
          endif

! If middle of averaging period of output fields is reached, accumulated
! deposited mass radioactively decays 
!***********************************************************************

    if (DEP.and.(jtime.eq.loutnext).and.(ldirect.gt.0)) then
      do ks=1,nspec
      do kp=1,maxpointspec_act
        if (decay(ks).gt.0.) then
          do nage=1,nageclass
            do l=1,nclassunc
  ! Mother output grid
              do jy=0,numygrid-1
                do ix=0,numxgrid-1
                  wetgridunc(ix,jy,ks,kp,l,nage)= &
                       wetgridunc(ix,jy,ks,kp,l,nage)* &
                       exp(-1.*outstep*decay(ks))
                  drygridunc(ix,jy,ks,kp,l,nage)= &
                       drygridunc(ix,jy,ks,kp,l,nage)* &
                       exp(-1.*outstep*decay(ks))
                end do
              end do
  ! Nested output grid
              if (nested_output.eq.1) then
                do jy=0,numygridn-1
                  do ix=0,numxgridn-1
                    wetgriduncn(ix,jy,ks,kp,l,nage)= &
                         wetgriduncn(ix,jy,ks,kp,l,nage)* &
                         exp(-1.*outstep*decay(ks))
                    drygriduncn(ix,jy,ks,kp,l,nage)= &
                         drygriduncn(ix,jy,ks,kp,l,nage)* &
                         exp(-1.*outstep*decay(ks))
                  end do
                end do
              endif
            end do
          end do
        endif
      end do
      end do
    endif

!!! CHANGE: These lines may be switched on to check the conservation
!!! of mass within FLEXPART

!       if (mod(jtime,loutsample).eq.0) then 
!          xm=0.
!          xm1=0.
!          do 247 j=1,numpart
!47          if (itra1(j).eq.jtime) xm1=xm1+xmass1(j,1)
!          xm=xm1
!          do 248 nage=1,nageclass
!            do 248 ix=0,numxgrid-1
!              do 248 jy=0,numygrid-1
!                do 248 l=1,nclassunc
!48        xm=xm+wetgridunc(ix,jy,1,l,nage)+drygridunc(ix,jy,1,l,nage)
!          write(*,'(i6,4f10.3)') jtime,xm,xm1
!       endif
!!! CHANGE

          
! Check whether concentrations are to be calculated
!**************************************************

        if ((ldirect*jtime.ge.ldirect*loutstart).and. &
        (ldirect*jtime.le.ldirect*loutend)) then ! add to grid
          if (mod(jtime-loutstart,loutsample).eq.0) then

! If we are exactly at the start or end of the concentration averaging interval,
! give only half the weight to this sample
!*******************************************************************************

            if ((jtime.eq.loutstart).or.(jtime.eq.loutend)) then
              weight=0.5
            else
              weight=1.0
            endif
            outnum=outnum+weight
            if(iout.ge.1) then 
             if (outgrid_option.eq.0) then 
             call conccalc_irreg(jtime,weight)
             elseif (outgrid_option.eq.1) then
             call conccalc_reg(jtime,weight)
             endif
            endif
          endif


!         if ((mquasilag.eq.1).and.(jtime.eq.(loutstart+loutend)/2)) &
!         call partoutput_short(jtime)    ! dump particle positions in extremely compressed format


! Output and reinitialization of grid
! If necessary, first sample of new grid is also taken
!*****************************************************

          if ((jtime.eq.loutend).and.(outnum.gt.0.)) then
!            print*,'iout',iout,ipout,outgrid_option
            if ((iout.le.3.).or.(iout.eq.5)) then 
             if(iout.ge.1) then
             if (outgrid_option.eq.0) then
             call concoutput_irreg(jtime,outnum,gridtotalunc, &
              wetgridtotalunc,drygridtotalunc)
       if (nested_output.eq.1) call concoutput_nest_irreg(jtime,outnum)
             elseif (outgrid_option.eq.1) then
             call concoutput_reg(jtime,outnum,gridtotalunc, &
              wetgridtotalunc,drygridtotalunc)
       if (nested_output.eq.1) call concoutput_nest_reg(jtime,outnum)
             endif
            endif

!             if (nested_output.eq.1.and.iout.ge.1)
!    +           call concoutput_nest(jtime,outnum)
              outnum=0.
            endif
            if ((iout.eq.4).or.(iout.eq.5)) call plumetraj(jtime)
            if (iflux.eq.1) call fluxoutput(jtime)
            write(*,45) jtime,numpart,gridtotalunc,wetgridtotalunc, &
            drygridtotalunc
45          format(i9,' SECONDS SIMULATED: ',i9, &
            ' PARTICLES:    Uncertainty: ',3f7.3)
!  call system_clock(clck_counts_end2,clck_rate2)
!  tins2=real(clck_counts_end2 - clck_counts_beg2)/real(clck_rate2)
!              print*,'time needed',tins2
!  call system_clock (clck_counts_beg2,clck_rate2)
            if (ipout.ge.1) call partoutput(jtime)    ! dump particle positions
            loutnext=loutnext+loutstep
            loutstart=loutnext-loutaver/2
            loutend=loutnext+loutaver/2
            if (jtime.eq.loutstart) then
              weight=0.5
              outnum=outnum+weight
              if(iout.ge.1) then
               if (outgrid_option.eq.0) then
               call conccalc_irreg(jtime,weight)
               elseif (outgrid_option.eq.1) then
               call conccalc_reg(jtime,weight)
               endif
             endif
            endif


! Check, whether particles are to be split:
! If so, create new particles and attribute all information from the old
! particles also to the new ones; old and new particles both get half the
! mass of the old ones
!************************************************************************

        if (ldirect*jtime.ge.ldirect*itsplit) then
          n=numpart
          do j=1,numpart
            if (ldirect*jtime.ge.ldirect*itrasplit(j)) then
              if (n.lt.maxpart) then
                n=n+1
                itrasplit(j)=2*(itrasplit(j)-itramem(j))+itramem(j)
                itrasplit(n)=itrasplit(j)
                itramem(n)=itramem(j)
                itra1(n)=itra1(j)
                idt(n)=idt(j)
                npoint(n)=npoint(j)
                nclass(n)=nclass(j)
                xtra1(n)=xtra1(j)
                ytra1(n)=ytra1(j)
                ztra1(n)=ztra1(j)
                uap(n)=uap(j)
                ucp(n)=ucp(j)
                uzp(n)=uzp(j)
                us(n)=us(j)
                vs(n)=vs(j)
                ws(n)=ws(j)
                cbt(n)=cbt(j)
                do ks=1,nspec
                  xmass1(j,ks)=xmass1(j,ks)/2.
                  xmass1(n,ks)=xmass1(j,ks)
                end do
              endif
            endif
          end do
          numpart=n
        endif
      endif
    endif
        



! Loop over all particles
!************************


!     chunksize=int(numpart/ntasks)+1  !if sent homogeneously
!        call itime(now)
!        ttime=now(1)*3600+now(2)*60+now(3)-ttime

!      print*,'processing time',ttime
   endif !over myid
!JB
! at this stage, I assume that each node has the same shared memory because they run getfields.
! now we need to split the trajectories into pieces for each node
!   if (myid.eq.0) then

        if (jtime.eq.ideltas) exit    

! Compute interval since radioactive decay of deposited mass was computed
!************************************************************************

        if (jtime.lt.loutnext) then
          ldeltat=jtime-(loutnext-loutstep)
        else                                  ! first half of next interval
          ldeltat=jtime-loutnext
        endif


!  if (myid.eq.0) then
!       call itime(now)
!        ttime=now(1)*3600+now(2)*60+now(3)
!   do ii=1,ntasks-1
!    call MPI_SEND(chunksize,1, MPI_INTEGER, ii,3001, MPI_COMM_WORLD, ierr)
!    call MPI_SEND(numpart,1, MPI_INTEGER, ii,3002, MPI_COMM_WORLD, ierr)
!   enddo 
!   else
!    call MPI_RECV(chunksize,1, MPI_INTEGER, 0,3001, MPI_COMM_WORLD,status, ierr)
!    call MPI_RECV(numpart,1, MPI_INTEGER, 0,3002, MPI_COMM_WORLD,status, ierr)
!   endif
!  print*,'numpart',numpart

!        call itime(now)
!        ttime=now(1)*3600+now(2)*60+now(3)

! initialize the temporary drydeposition grid

            if (DRYDEP.and.ldirect.gt.0) then
  do ks=1,nspec
  do kp=1,maxpointspec_act
    do nage=1,nageclass
      do jy=0,numygrid-1
        do ix=0,numxgrid-1
          do l=1,nclassunc
            drygridunc2(ix,jy,ks,kp,l,nage)=0.
          end do
        end do
      end do
    if (nested_output.eq.1) then
      do jy=0,numygridn-1
        do ix=0,numxgridn-1
          do l=1,nclassunc
            drygriduncn2(ix,jy,ks,kp,l,nage)=0.
          end do
        end do
      end do
    endif
    end do
  end do
  end do
   endif


!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(jj, &
!$OMP decfact,clck_counts_end,clck_rate,clck_counts_beg,tins,&
!$OMP  xold, yold, zold, kp, itage, prob, nstop, xmassfract, &
!$OMP chunksize3,start,finish,ngrid,ompid,depoindicator,nbp, &
!$OMP indzindicator,cpttra,drydeposit,n_threads) &
!$OMP SHARED(height,rho,tt,vsetaver,dquer,xtra1,ytra1,ztra1, &
!$OMP density,cunningham,itra1,ioutputforeachrelease,cbt,iflux, &
!$OMP uun,vvn,wwn,ustar,wstar,oli,uupol,vvpol,uu,vv,ww,drhodz,ptt,tke, &
!$OMP rhon,drhodzn,pttn,tken,vdep,vdepn,itramem,nageclass,lage, &
!$OMP jtime,ldirect,memind,nglobal,switchnorthg,m_x,m_y,m_xn,m_yn, &
!$OMP switchsouthg,numbnests,xln,xrn,yln,yrn,memtime,xresoln, &
!$OMP yresoln,hmix,hmixn,tropopause, & 
!$OMP tropopausen,lsynctime,dxconst,dyconst,mdomainfill, &
!$OMP turb_option,turbswitch,ifine,chunksize, &
!!!maxrand, &
!$OMP xmass,xmass1,DRYDEP,DRYDEPSPEC,nspec,rannumb,uniform_rannumb,cpt, &
!$OMP lwindinterv,npart,npoint,idt,uap,ucp,uzp,us,vs,ws, &
!$OMP linit_cond,decay,ldeltat,nclass,nested_output,numpart, &
!$OMP drygridunc2,drygriduncn2, & 
!$OMP  mts)

!        call itime(now)
!        ttime=now(1)*3600+now(2)*60+now(3)
  call system_clock (clck_counts_beg,clck_rate)

!       chunksize3=int(chunksize2/omp_get_num_threads())+1
        n_threads=omp_get_num_threads()
        ompid=OMP_GET_THREAD_NUM()
        chunksize3=int(real(numpart)/real(n_threads)/20.)+1 !more efficient

!       chunksize3=int(chunksize2/omp_get_num_threads())+1
!       chunksize3=int(real(chunksize2)/real(omp_get_num_threads())/20.)+1 !more efficient

!       ompid=omp_get_num_threads()
!       ompid=0

        cpttra=0
!        print*,'chunksi',chunksize2,myid
        if (numpart.gt.0 ) then
!         print*,'test rel2',npoint(5139),npoint(6002),npoint(100003)
!        do jj=1,numpart
!        do jj=numpart,1,-1
!        print*,jj
!!!$OMP DO SCHEDULE(dynamic, max(1,int(real(numpart)/real(n_threads)/20.)))
!$OMP DO SCHEDULE(dynamic, max(1,numpart/1000))
        do jj=1,numpart

! If integration step is due, do it
!**********************************
          if (itra1(jj).eq.jtime) then
           cpttra=cpttra+1
        if (ioutputforeachrelease.eq.1) then
            kp=npoint(jj)
        else
            kp=1
        endif

! Determine age class of the particle
            itage=abs(itra1(jj)-itramem(jj))
            do nage=1,nageclass
              if (itage.lt.lage(nage)) exit
         enddo

             nbp=1
            if ((itramem(jj).eq.jtime).or.(jtime.eq.0)) &
           call initialize(jtime,idt(jj),uap(jj),ucp(jj),uzp(jj), &
        us(jj),vs(jj),ws(jj),xtra1(jj),ytra1(jj),ztra1(jj),cbt(jj), &
      ngrid,depoindicator,indzindicator,cpt(nbp),ompid,myid,n_threads,mts )

! Memorize particle positions
!****************************

            xold=xtra1(jj)
            yold=ytra1(jj)
            zold=ztra1(jj)
! Integrate Lagevin equation for lsynctime seconds
!*************************************************
!              write(*,*)'numpart,jtime, particle #=',numpart,jtime,j
        call advance(jtime,npoint(jj),idt(jj),uap(jj),ucp(jj),uzp(jj),us(jj), &
         vs(jj),ws(jj),nstop,xtra1(jj),ytra1(jj),ztra1(jj),prob,cbt(jj), &
      ngrid,depoindicator,indzindicator,cpt(nbp),ompid,myid,n_threads,mts )
! Calculate the gross fluxes across layer interfaces
!***************************************************


            if (iflux.eq.1) call calcfluxes(nage,jj,xold,yold,zold)


! Determine, when next time step is due
! If trajectory is terminated, mark it
!**************************************

        if (nstop.gt.1) then
          if (linit_cond.ge.1) call initial_cond_calc(jtime,jj)
              itra1(jj)=-999999999
            else
              itra1(jj)=jtime+lsynctime


! Dry deposition and radioactive decay for each species
!******************************************************
          xmassfract=0.

              do ks=1,nspec
                if (decay(ks).gt.0.) then             ! radioactive decay
                  decfact=exp(-real(abs(lsynctime))*decay(ks))
                else
                  decfact=1.
                endif

                if (DRYDEPSPEC(ks)) then        ! dry deposition
                  drydeposit(ks)=xmass1(jj,ks)*prob(ks)*decfact
                  xmass1(jj,ks)=xmass1(jj,ks)*(1.-prob(ks))*decfact
                  if (decay(ks).gt.0.) then   ! correct for decay (see wetdepo)
                    drydeposit(ks)=drydeposit(ks)* &
                    exp(real(abs(ldeltat))*decay(ks))
                  endif
                else                           ! no dry deposition
                  xmass1(jj,ks)=xmass1(jj,ks)*decfact
                endif

!           if (mdomainfill.eq.0) then ! JB bugfix
            if (mdomainfill.eq.0 .and. ipout.eq.0) then
              if (xmass(npoint(jj),ks).gt.0.) &
                   xmassfract=max(xmassfract,real(npart(npoint(jj)))* &
                    xmass1(jj,ks)/xmass(npoint(jj),ks))
            else
              xmassfract=1.
            endif

            end do

          if (xmassfract.lt.0.000001) then   ! terminate all particles carrying less mass
            itra1(jj)=-999999999
          endif

  !        Sabine Eckhardt, June 2008
  !        don't create depofield for backward runs
          if (DRYDEP.AND.(ldirect.eq.1)) then
           call drydepokernel(nclass(jj),drydeposit,real(xtra1(jj)), &
                 real(ytra1(jj)),itage,nage,kp)
            if (nested_output.eq.1) call drydepokernel_nest( &
              nclass(jj),drydeposit,real(xtra1(jj)),real(ytra1(jj)), &
                 itage,nage,kp)
          endif

  ! Terminate trajectories that are older than maximum allowed age
  !***************************************************************

          if (abs(itra1(jj)-itramem(jj)).ge.lage(nageclass)) then
            if (linit_cond.ge.1) &
                  call initial_cond_calc(jtime+lsynctime,jj)
            itra1(jj)=-999999999
          endif
      endif

      endif

    end do !loop over particles
!$OMP END DO 
!!!$OMP END DO NOWAIT
    endif
!      print*,'time',ttime,cpttra,myid,OMP_GET_THREAD_NUM()
             if (option_verbose.eq.1) then
  call system_clock(clck_counts_end,clck_rate)
  tins=real(clck_counts_end - clck_counts_beg)/real(clck_rate)
       print*,'time',tins,cpttra,myid,ompid
       endif
!$OMP END PARALLEL



!        call itime(now)
!        ttime=now(1)*3600+now(2)*60+now(3)-ttime
!      print*,'time',ttime,cpttra,myid,OMP_GET_THREAD_NUM()

! update the drydepo
            if (DRYDEP.and.ldirect.gt.0) then
  do ks=1,nspec
  do kp=1,maxpointspec_act
    do nage=1,nageclass

      do jy=0,numygrid-1
        do ix=0,numxgrid-1
          do l=1,nclassunc
       drygridunc(ix,jy,ks,kp,l,nage)=drygridunc(ix,jy,ks,kp,l,nage) &
               +drygridunc2(ix,jy,ks,kp,l,nage)
          end do
        end do
      end do
    if (nested_output.eq.1) then
      do jy=0,numygridn-1
        do ix=0,numxgridn-1
          do l=1,nclassunc
       drygriduncn(ix,jy,ks,kp,l,nage)=drygriduncn(ix,jy,ks,kp,l,nage) &
               +drygriduncn2(ix,jy,ks,kp,l,nage)
          end do
        end do
      end do
    endif
    end do
  end do
  end do

   endif

  end do !loop over time 


  ! Complete the calculation of initial conditions for particles not yet terminated
  !*****************************************************************************

  do j=1,numpart
    if (linit_cond.ge.1) call initial_cond_calc(jtime,j)
  end do

  if (ipout.eq.2) call partoutput(jtime)     ! dump particle positions

  if (linit_cond.ge.1) call initial_cond_output(jtime)   ! dump initial cond. field

  close(104)

  ! De-allocate memory and end
  !***************************
  if (iouttype.eq.2) then
  write(*,*) "timemanager: End of simulation reached, closing output file"
    ncret=NF90_CLOSE(ncout)
   call check_ncerror(ncret)
  if (nested_output.eq.1)  ncret=NF90_CLOSE(ncoutn)
  endif

  if (iflux.eq.1) then
      deallocate(flux)
  endif
  if (OHREA.eqv..TRUE.) then
      deallocate(OH_field,OH_field_height)
  endif
  deallocate(gridunc)
  deallocate(xpoint1,xpoint2,ypoint1,ypoint2,zpoint1,zpoint2,xmass)
  deallocate(ireleasestart,ireleaseend,npart,kindz)
!  deallocate(xmasssave)
  if (myid.eq.0) then
  if (nested_output.eq.1) then
     deallocate(orooutn, arean, volumen)
     if (ldirect.gt.0) then
     deallocate(griduncn,drygriduncn,wetgriduncn,drygriduncn2)
     endif
  endif
  if (ldirect.gt.0) then
      if (allocated(drygridunc)) deallocate(drygridunc)
      if (allocated(wetgridunc)) deallocate(wetgridunc)
      if (allocated(drygridunc2)) deallocate(drygridunc2)
      if (allocated(drygriduncn2)) deallocate(drygriduncn2)
  endif
  deallocate(outheight,outheighthalf)
  deallocate(oroout, area, volume)
  endif
end subroutine timemanager



