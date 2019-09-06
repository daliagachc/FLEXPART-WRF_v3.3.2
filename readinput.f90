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
subroutine readinput
!*******************************************************************************
!                                                                              *
!     Reads the pathnames, where input/output files are expected to be.        *
!     The file pathnames must be available in the current working directory.   *
!                                                                              *
!     Author: A. Stohl                                                         *
!                                                                              *
!     1 February 1994                                                          *
!                                                                              *
!   2015-05-07 A. Dingwell: Changes as suggested by H.M.J. Barbosa:            *
!                 - Removed unused variables                                   *
!                 - Renamed beg,end to dbeg,deng (since end is a keyword)      *
!                 - Always initialize numpart_in                               *
!  2015-05-07 A.Dingwell:                                                      *
!                 - Fixed indentation to match FLEXPART standard               *
!                 - Corrected misleading error messages                        *
!                 - Added additional debug info for errors                     *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
! len(numpath)       lengths of the path names                                 *
! path(numpath)      pathnames of input/output files                           *
!                                                                              *
! Constants:                                                                   *
! numpath            number of pathnames to be read in                         *
!                                                                              *
!*******************************************************************************
  use par_mod
  use com_mod
  use outg_mod
  use point_mod
  !use xmass_mod

  implicit none
!      include 'includepar'
!      include 'includecom'

  integer :: i, j,numtable,numpoint2,itimein,numpart_in
  integer :: hhh,mi,ss,pos_spec,iomode_xycoord_in
  !hmjb not used:      character(len=50) :: line
  !hmjb not used:      logical :: old
  integer :: idiff,ldat,ltim,wftime1(maxwf),numbwfn(maxnests),k
  integer :: wftime1n(maxnests,maxwf),wftimen(maxnests,maxwf)
  !hmjb      real(kind=dp) :: jul,juldate,beg,end,jul1,jul2
  !hmjb - end should not be used as a variable name.
  !  new names: beg=> dbeg , end=> dend
  real(kind=dp) :: jul,juldate,dbeg,dend,jul1,jul2
!      double precision juldate,jul,beg,end,jul1,jul2
  character(len=80) :: fname,wfname1(maxwf),wfname1n(maxnests,maxwf)
  character(len=10) :: spec, wfspec1(maxwf),wfspec1n(maxnests,maxwf)
  real :: outhelp,xr,xr1,yr,yr1
  real :: xtmp, xtmp1, xtmp2, ytmp, ytmp1, ytmp2
! 10-mar-2006 rce - flexpart_wrf - eps should be a small dx/dy value in meters
!     parameter(eps=1.e-4)
  real,parameter :: eps=10.0
  real :: x,y,xm,ym
  character(len=16) ::  receptor
  integer :: numpartmax,id1,it1,id2,it2,idow,ihour!hmjb not used - ,nparttoread
  integer :: emitvar,stat,nparttot,nparttot2,outgriddef,outgriddefn
  real :: vsh(ni),fracth(ni),schmih(ni),releaserate,releaserate2
  !hmjb not used  character(len=3) :: aspecnumb
  character(len=10) :: specname(maxtable)
  real :: cun
  logical :: spec_found
! Read the pathname information stored in unitpath
!*************************************************
  nparttot=0
  numpart_in=0
!     print*,'path2'
!     print*,'path',unitpath
!     open(unitpath,file='flexwrf.input',status='old',err=801)
  write(*,*) "Opening '",inputname,"' for reading"
  open(unitpath,file=inputname,status='old',err=801)

! jdf start read pathnames
!*************************************************
  write(*,*) "Reading pathnames"
  call skplin(1,unitpath)
  do i=1,numpath
    read(unitpath,'(a)',err=800) path(i) 
    length(i)=index(path(i),' ')-1
  enddo

! Check whether any nested subdomains are to be used
!***************************************************

  do i=1,maxnests
    read(unitpath,'(a)') path(numpath+2*(i-1)+1) 
    read(unitpath,'(a)') path(numpath+2*(i-1)+2) 
    if (path(numpath+2*(i-1)+1)(1:5).eq.'=====') goto 30
    length(numpath+2*(i-1)+1)=index(path(numpath+2*(i-1)+1),' ')-1
    length(numpath+2*(i-1)+2)=index(path(numpath+2*(i-1)+2),' ')-1
  enddo

  if (path(numpath+2*(i-1)+1)(1:5).ne.'=====') then
    write(*,*) ' #### FLEXPART MODEL ERROR! NUMBER OF PATHS   #### '
    write(*,*) ' #### IN FORMER PATHNAMES FILE GREATER THAN   #### '
    write(*,*) ' #### MAXNESTS IN PAR_MOD.F90, INCREASE THIS  #### '
    write(*,*) ' #### OR REMOVE EXTRA PATHS FROM PATHNAMES    #### '
    stop
  endif

! Determine number of available nested domains
!*********************************************

30    numbnests=i-1

! jdf end read pathnames
!*************************************************
! jdf start read command
!*************************************************
!
!*******************************************************************************
!                                                                              *
!     Note:  This is the FLEXPART_WRF version of subroutine readcommand.       *
!                                                                              *
!     This routine reads the user specifications for the current model run.    *
!                                                                              *
!     Author: A. Stohl                                                         *
!                                                                              *
!     18 May 1996                                                              *
!                                                                              *
!     Nov-Dec-2005, R. Easter - input turb_option, add_sfc_level, iouttype     *
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
! bdate                beginning date as Julian date                           *
! ctl                  factor by which time step must be smaller than          *
!                      Lagrangian time scale                                   *
! edate                ending date as Julian date                              *
! hhh                  hour                                                    *
! ibdate,ibtime        beginnning date and time (YYYYMMDD, HHMISS)             *
! ideltas [s]          modelling period                                        *
! iedate,ietime        ending date and time (YYYYMMDD, HHMISS)                 *
! ifine                reduction factor for vertical wind time step            *
! iflux                switch to turn on (1)/off (0) flux calculations         *
! iout                 1 for conc. (residence time for backward runs) output,  *
!                      2 for mixing ratio output, 3 both, 4 for plume          *
!                      trajectory output, 5 = options 1 and 4                  *
! ipin                 1 continue simulation with dumped particle data, 0 no   *
! ipout                0 no particle dump, 1 every output time, 3 only at end  *
! itsplit [s]          time constant for particle splitting                    *
! loutaver [s]         concentration output is an average over loutaver seconds*
! loutsample [s]       average is computed from samples taken every [s] seconds*
! loutstep [s]         time interval of concentration output                   *
! lsynctime [s]        synchronisation time interval for all particles         *
! lagespectra          switch to turn on (1)/off (0) calculation of age spectra*
! lsubgrid             switch to turn on (1)/off (0) subgrid topography        *
!                      parameterization                                        *
! lconvection          value of either 0 and 1 indicating mixing by convection *
!                      = 0 .. no convection                                    *
!                      + 1 .. parameterisation of mixing by subgrid-scale      *
!                              convection = on                                 *
! method               method used to compute the particle pseudovelocities    *
! mdomainfill          1 use domain-filling option, 0 not, 2 use strat. O3     *
! lu_option            Land-use category source: 0 for IGBP_int1.dat; 1 for WRF*
! mi                   minute                                                  *
! ss                   second                                                  *
!                                                                              *
! Constants:                                                                   *
! unitcommand          unit connected to file COMMAND                          *
!                                                                              *
! 9-10-2007, W Wang                                                            *
!   add turb_option_tke and turb_option_mytke                                  *
! 1 Oct, 2007                                                                  *
!   add dt_conv input,  time intervals to call convection, seconds             *
! 2015-03-25, Adam Dingwell                                                    *
!   Added lu_option, for toggling land-use source (WRF or IGBP_int1.dat)       *
!*******************************************************************************


! Open the command file and read user options
!********************************************

  read(unitcommand,*) ldirect
  read(unitcommand,*) ibdate,ibtime
  read(unitcommand,*) iedate,ietime
  read(unitcommand,*) loutstep
  read(unitcommand,*) loutaver
  read(unitcommand,*) loutsample
  read(unitcommand,*) itsplit
  read(unitcommand,*) lsynctime
  read(unitcommand,*) ctl
  read(unitcommand,*) ifine
  read(unitcommand,*) iout
  read(unitcommand,*) ipout
  read(unitcommand,*) lsubgrid
  read(unitcommand,*) lconvection
  read(unitcommand,*) dt_conv
  read(unitcommand,*) lagespectra
  read(unitcommand,*) ipin
  read(unitcommand,*) iflux
  read(unitcommand,*) ioutputforeachrelease
  read(unitcommand,*) mdomainfill
  read(unitcommand,*) ind_source
  read(unitcommand,*) ind_receptor
  read(unitcommand,*) nested_output
  read(unitcommand,*) linit_cond
! FLEXPART_WRF - read turb_option, add_sfc_level
  read(unitcommand,*) turb_option
  read(unitcommand,*) lu_option
  read(unitcommand,*) cblflag  ! added by mc for cbl option
!     read(unitcommand,*) add_sfc_level
  add_sfc_level=1
  read(unitcommand,*) sfc_option
  read(unitcommand,*) wind_option
  read(unitcommand,*) time_option
  read(unitcommand,*) outgrid_option
  read(unitcommand,*) numpoint_option
  read(unitcommand,*) iouttype
  read(unitcommand,*) ncnumrec
  read(unitcommand,*) option_verbose

  write(*,*) "option_verbose=",option_verbose

  ifine=max(ifine,1)

! Determine how Markov chain is formulated (for w or for w/sigw)
!***************************************************************

  if (ctl.ge.0.1) then
    turbswitch=.true.
  else
    turbswitch=.false.
    ifine=1
  endif
  fine=1./real(ifine)
  ctl=1./ctl

! Set the switches required for the various options for input/output units
!*************************************************************************
!AF Set the switches IND_REL and IND_SAMP for the release and sampling
!Af switches for the releasefile:
!Af IND_REL =  1 : xmass * rho
!Af IND_REL =  0 : xmass * 1

!Af switches for the conccalcfile:
!AF IND_SAMP =  0 : xmass * 1
!Af IND_SAMP = -1 : xmass / rho

!AF IND_SOURCE switches between different units for concentrations at the source
!Af   NOTE that in backward simulations the release of computational particles 
!Af   takes place at the "receptor" and the sampling of p[articles at the "source".
!Af          1 = mass units 
!Af          2 = mass mixing ratio units 
!Af IND_RECEPTOR switches between different units for concentrations at the receptor
!Af          1 = mass units 
!Af          2 = mass mixing ratio units 

  if ( ldirect .eq. 1 ) then  ! FWD-Run
!Af set release-switch
     if ( IND_SOURCE .eq. 1 ) then !mass
        ind_rel = 0
     else ! mass mix
        ind_rel = 1
     endif
!Af set sampling switch
     if ( IND_RECEPTOR .eq. 1) then !mass
        ind_samp = 0
     else ! mass mix
        ind_samp = -1
     endif
  elseif (ldirect .eq. -1 ) then !BWD-Run
!Af set sampling switch
     if ( IND_SOURCE .eq. 1 ) then !mass
        ind_samp = -1
     else ! mass mix
        ind_samp = 0
     endif
!Af set release-switch
     if ( IND_RECEPTOR .eq. 1) then !mass
        ind_rel = 1
     else ! mass mix
        ind_rel = 0
     endif
  else
    write(*,*) '#### FLEXPART MODEL ERROR! LDIRECT=',ldirect
    write(*,*) '#### IS NOT A VALID OPTION'
  endif



!*************************************************************
! Check whether valid options have been chosen in file COMMAND
!*************************************************************

! Check input dates
!******************

  if (iedate.lt.ibdate) then
    write(*,*) ' #### FLEXPART MODEL ERROR! BEGINNING DATE    #### '
    write(*,*) ' #### IS LARGER THAN ENDING DATE. CHANGE      #### '
    write(*,*) ' #### EITHER POINT 2 OR POINT 3               #### '
    stop
  else if (iedate.eq.ibdate) then
    if (ietime.lt.ibtime) then
    write(*,*) ' #### FLEXPART MODEL ERROR! BEGINNING TIME    #### '
    write(*,*) ' #### IS LARGER THAN ENDING TIME. CHANGE      #### '
    write(*,*) ' #### EITHER POINT 2 OR POINT 3               #### '
    stop
    endif
  endif


! Determine kind of dispersion method
!************************************

  if (ctl.gt.0.) then    
    method=1
    mintime=minstep
  else
    method=0
    mintime=lsynctime
  endif

! Check whether a valid option for gridded model output has been chosen
!**********************************************************************

!     if ((iout.lt.0).or.(iout.gt.5)) then
!       write(*,*) ' #### FLEXPART MODEL ERROR! FILE COMMAND:     #### '
!       write(*,*) ' #### IOUT MUST BE 0, 1, 2, 3, 4, OR 5!       #### '
!       stop
!     endif

!AF check consistency between units and volume mixing ratio
  if ( ((iout.eq.2).or.(iout.eq.3)).and. &
       (ind_source.gt.1 .or.ind_receptor.gt.1) ) then
    write(*,*) ' #### FLEXPART MODEL ERROR!             :     #### '
    write(*,*) ' #### VOLUME MIXING RATIO ONLY SUPPORTED      #### '
    write(*,*) ' #### FOR MASS UNITS (at the moment)          #### '
    stop
  endif

  ! For quasilag output for each release is forbidden
  !*****************************************************************************

  if ((ioutputforeachrelease.eq.1).and.(mquasilag.eq.1)) then
      write(*,*) '#### FLEXPART MODEL ERROR!             :     ####'
      write(*,*) '#### OUTPUTFOREACHRELEASE AND QUASILAGRANGIAN####'
      write(*,*) '#### MODE IS NOT POSSIBLE   !                ####'
      stop
  endif

  ! For quasilag backward is forbidden
  !*****************************************************************************

  if ((ldirect.lt.0).and.(mquasilag.eq.1)) then
      write(*,*) '#### FLEXPART MODEL ERROR!             :     ####'
      write(*,*) '#### FOR BACKWARD RUNS, QUASILAGRANGIAN MODE ####'
      write(*,*) '#### IS NOT POSSIBLE   !                     ####'
      stop
  endif
  ! For backward runs one releasefield for all releases makes no sense,
  ! For quasilag and domainfill ioutputforechrelease is forbidden
  !*****************************************************************************

!    print*,'ioutput',ioutputforeachrelease,ldirect
  if ((ldirect.lt.0).and.(ioutputforeachrelease.eq.0)) then
      write(*,*) '#### FLEXPART MODEL ERROR!             :     ####'
      write(*,*) '#### FOR BACKWARD RUNS, IOUTPUTFOREACHRLEASE ####'
      write(*,*) '#### MUST BE SET TO ONE!                     ####'
      stop
  endif

  ! For backward runs one releasefield for all releases makes no sense,
  ! and is "forbidden"
  !*****************************************************************************

  if ((mdomainfill.eq.1).and.(ioutputforeachrelease.eq.1)) then
      write(*,*) '#### FLEXPART MODEL ERROR!             :     ####'
      write(*,*) '#### FOR DOMAIN FILLING RUNS OUTPUT FOR      ####'
      write(*,*) '#### EACH RELEASE IS FORBIDDEN !             ####'
      stop
  endif


  ! For domain-filling trajectories, a plume centroid trajectory makes no sense,
  ! For backward runs, only residence time output (iout=1) or plume trajectories (iout=4),
  ! or both (iout=5) makes sense; other output options are "forbidden"
  !*****************************************************************************

  if (ldirect.lt.0) then
    if ((iout.eq.2).or.(iout.eq.3)) then
      write(*,*) '#### FLEXPART MODEL ERROR!             :     ####'
      write(*,*) '#### FOR BACKWARD RUNS, IOUT MUST BE 1,4,OR 5####'
      stop
    endif
  endif
  ! For domain-filling trajectories, a plume centroid trajectory makes no sense,
  ! and is "forbidden"
  !*****************************************************************************

  if (mdomainfill.ge.1) then
    if ((iout.eq.4).or.(iout.eq.5)) then
      write(*,*) '#### FLEXPART MODEL ERROR!             :     ####'
      write(*,*) '#### FOR DOMAIN-FILLING TRAJECTORY OPTION,   ####'
      write(*,*) '#### IOUT MUST NOT BE SET TO 4 OR 5.         ####'
      stop
    endif
  endif


! Check whether a valid options for particle dump has been chosen
!****************************************************************

  if ((ipout.ne.0).and.(ipout.ne.1).and.(ipout.ne.2)) then
    write(*,*) ' #### FLEXPART MODEL ERROR!             :     #### '
    write(*,*) ' #### IPOUT MUST BE 1, 2 OR 3!                #### '
    stop
  endif

  if(lsubgrid.ne.1) then
    write(*,*) '             ----------------               '
    write(*,*) ' INFORMATION: SUBGRIDSCALE TERRAIN EFFECT IS'
    write(*,*) ' NOT PARAMETERIZED DURING THIS SIMULATION.  '
    write(*,*) '             ----------------               '
  endif
   
  if ((ipout.ne.0).and.(ipout.ne.1).and.(ipout.ne.2)) then
    write(*,*) ' #### FLEXPART MODEL ERROR!             :     #### '
    write(*,*) ' #### IPOUT MUST BE 1, 2 OR 3!                #### '
    stop
  endif


! Check whether convection scheme is either turned on or off
!***********************************************************

!      if ((lconvection.ne.0).and.(lconvection.ne.1)) then
  if ((lconvection.lt.0).or. (lconvection.gt.3)) then
    write(*,*) ' #### FLEXPART MODEL ERROR!             :     #### '
    !        write(*,*) ' #### LCONVECTION MUST BE SET TO EITHER 1 OR 0#### '
    write(*,*) ' #### LCONVECTION MUST BE SET TO 0 or  3  #### '

    stop
  endif


! Check whether synchronisation interval is sufficiently short
!*************************************************************

  if (lsynctime.gt.(idiffnorm/2)) then
    write(*,*) ' #### FLEXPART MODEL ERROR! SYNCHRONISATION   #### '
    write(*,*) ' #### TIME IS TOO LONG. MAKE IT SHORTER.      #### '
    stop
  endif


! Check consistency of the intervals, sampling periods, etc., for model output
!*****************************************************************************

  if (loutaver.eq.0) then
    write(*,*) ' #### FLEXPART MODEL ERROR! TIME AVERAGE OF   #### '
    write(*,*) ' #### CONCENTRATION FIELD OUTPUT MUST NOT BE  #### '
    write(*,*) ' #### ZERO.                                   #### '
    write(*,*) ' #### CHANGE INPUT                            #### '
    stop
  endif

  if (loutaver.gt.loutstep) then
    write(*,*) ' #### FLEXPART MODEL ERROR! TIME AVERAGE OF   #### '
    write(*,*) ' #### CONCENTRATION FIELD OUTPUT MUST NOT BE  #### '
    write(*,*) ' #### GREATER THAN INTERVAL OF OUTPUT.        #### '
    write(*,*) ' #### CHANGE INPUT                            #### '
    stop
  endif

  if (loutsample.gt.loutaver) then
    write(*,*) ' #### FLEXPART MODEL ERROR! SAMPLING TIME OF  #### '
    write(*,*) ' #### CONCENTRATION FIELD OUTPUT MUST NOT BE  #### '
    write(*,*) ' #### GREATER THAN TIME AVERAGE OF OUTPUT.    #### '
    write(*,*) ' #### CHANGE INPUT                            #### '
    stop
  endif

  if (mod(loutaver,lsynctime).ne.0) then
    write(*,*) ' #### FLEXPART MODEL ERROR! AVERAGING TIME OF #### '
    write(*,*) ' #### CONCENTRATION FIELD MUST BE A MULTIPLE  #### '
    write(*,*) ' #### OF THE SYNCHRONISATION INTERVAL         #### '
    stop
  endif

  if ((loutaver/lsynctime).lt.2) then
    write(*,*) ' #### FLEXPART MODEL ERROR! AVERAGING TIME OF #### '
    write(*,*) ' #### CONCENTRATION FIELD MUST BE AT LEAST    #### '
    write(*,*) ' #### TWICE THE SYNCHRONISATION INTERVAL      #### '
    stop
  endif

  if (mod(loutstep,lsynctime).ne.0) then
    write(*,*) ' #### FLEXPART MODEL ERROR! INTERVAL BETWEEN  #### '
    write(*,*) ' #### CONCENTRATION FIELDS MUST BE A MULTIPLE #### '
    write(*,*) ' #### OF THE SYNCHRONISATION INTERVAL         #### '
    stop
  endif

  if ((loutstep/lsynctime).lt.2) then
    write(*,*) ' #### FLEXPART MODEL ERROR! INTERVAL BETWEEN  #### '
    write(*,*) ' #### CONCENTRATION FIELDS MUST BE AT LEAST   #### '
    write(*,*) ' #### TWICE THE SYNCHRONISATION INTERVAL      #### '
    stop
  endif

  if (mod(loutsample,lsynctime).ne.0) then
    write(*,*) ' #### FLEXPART MODEL ERROR! SAMPLING TIME OF  #### '
    write(*,*) ' #### CONCENTRATION FIELD MUST BE A MULTIPLE  #### '
    write(*,*) ' #### OF THE SYNCHRONISATION INTERVAL         #### '
    stop
  endif

  if (itsplit.lt.loutaver) then
    write(*,*) ' #### FLEXPART MODEL ERROR! SPLITTING TIME FOR#### '
    write(*,*) ' #### PARTICLES IS TOO SHORT. PLEASE INCREASE #### '
    write(*,*) ' #### SPLITTING TIME CONSTANT.                #### '
    stop
  endif

  if ((mquasilag.eq.1).and.(iout.ge.4)) then
    write(*,*) ' #### FLEXPART MODEL ERROR! CONFLICTING       #### '
    write(*,*) ' #### OPTIONS: IF MQUASILAG=1, PLUME          #### '
    write(*,*) ' #### TRAJECTORY OUTPUT IS IMPOSSIBLE.        #### '
    stop
  endif

  if (ncnumrec.le.0) then
    write(*,*) ' #### FLEXPART MODEL ERROR! NUMREC MUST       #### '
    write(*,*) ' #### BE GREATER THAN 0                       #### '
    stop
  endif

! FLEXPART_WRF - check turb_option, add_sfc_level
  if ((turb_option.ne.turb_option_none     ) .and. &
      (turb_option.ne.turb_option_diagnosed) .and. &
      (turb_option.ne.turb_option_tke      ) .and. &
      (turb_option.ne.turb_option_mytke)) then
    write(*,*) ' #### FLEXPART MODEL ERROR!                   #### '
    write(*,*) ' #### TURB_OPTION MUST BE ONE OF:             #### '
    write(*,'(5x,5i5)') turb_option_none, turb_option_diagnosed, &
               turb_option_tke,turb_option_mytke
    write(*,*) ' #### ---------------------------------       #### '
    stop
  endif

  if ((add_sfc_level.ne.0) .and. (add_sfc_level.ne.1)) then
    write(*,*) ' #### FLEXPART MODEL ERROR!                   #### '
    write(*,*) ' #### ADD_SFC_LAYER MUST BE 0 or 1            #### '
    stop
  endif

  if (sfc_option.eq.sfc_option_diagnosed) then
    if (option_verbose.ge.1) write(*,*) 'SFC_OPTION =',sfc_option
  elseif (sfc_option.eq.sfc_option_wrf) then
    write(*,*) ' #### FLEXPART MODEL ERROR! SFC_OPTION =',sfc_option
    write(*,*) ' #### Reading from WRF no longer supported.'
    write(*,*) ' #### ---------------------------------'
    stop
  else
    write(*,*) ' #### FLEXPART MODEL ERROR!                   #### '
    write(*,*) ' #### SFC_OPTION MUST BE ONE OF:              #### '
    write(*,'(5x,5i5)') sfc_option_diagnosed, sfc_option_wrf
    write(*,*) ' #### ---------------------------------'
    stop
  endif

! iouttype -- convert negative values to 0; positive out of range values to 2
  if (iouttype .lt. 0) iouttype = 0
  if (iouttype .gt. 2) iouttype = 2

! Conversion of format HHHMISS to seconds
!****************************************

  hhh=ideltas/10000
  mi=(ideltas-10000*hhh)/100
  ss=ideltas-10000*hhh-100*mi
  ideltas=hhh*3600+60*mi+ss


! Compute modeling time in seconds and beginning date in Julian date
!*******************************************************************

  outstep=real(abs(loutstep))
  if (ldirect.eq.1) then
    bdate=juldate(ibdate,ibtime)
    edate=juldate(iedate,ietime)
    ideltas=nint((edate-bdate)*86400.)
  else if (ldirect.eq.-1) then
    loutaver=-1*loutaver
    loutstep=-1*loutstep
    loutsample=-1*loutsample
    lsynctime=-1*lsynctime
    bdate=juldate(iedate,ietime)
    edate=juldate(ibdate,ibtime)
    ideltas=nint((edate-bdate)*86400.)
  else
    write(*,*) ' #### FLEXPART MODEL ERROR! DIRECTION IN      #### '
    write(*,*) ' #### INPUT FILE     MUST BE EITHER -1 OR 1.  #### '
    stop
  endif

! jdf end read command
!*************************************************
! jdf start read ageclasees
!*************************************************
!
!*******************************************************************************
!                                                                              *
!     This routine reads the age classes to be used for the current model run. *
!                                                                              *
!     Author: A. Stohl                                                         *
!                                                                              *
!     20 March 2000                                                            *
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
!                                                                              *
! Constants:                                                                   *
!                                                                              *
!*******************************************************************************

! If age spectra claculation is switched on,
! open the AGECLASSSES file and read user options
!************************************************

  call skplin(1,unitpath)
  read(unitpath,*) nageclass
  read(unitpath,*) lage(1)
  if (lagespectra.eq.1 .and. nageclass.gt.maxageclass) then
    write(*,*) ' #### FLEXPART MODEL ERROR! NUMBER OF AGE     #### '
    write(*,*) ' #### CLASSES GREATER THAN MAXIMUM ALLOWED.   #### '
    write(*,*) ' #### CHANGE SETTINGS IN      AGECLASSES OR   #### '
    write(*,*) ' #### RECOMPILE WITH LARGER MAXAGECLASS IN    #### '
    write(*,*) ' #### FILE PAR_MOD.F90                        #### '
    stop
  endif

  if (lage(1).le.0) then
    write(*,*) ' #### FLEXPART MODEL ERROR! AGE OF FIRST      #### '
    write(*,*) ' #### CLASS MUST BE GREATER THAN ZERO. CHANGE #### '
    write(*,*) ' #### SETTINGS IN      AGECLASSES.            #### '
    stop
  endif

! If age spectra calculation is switched off, set number of age classes
! to 1 and maximum age to a large number
!**********************************************************************

  if(option_verbose .ge. 1) print*,'nageclass',nageclass 
  do i=2,nageclass
    read(unitpath,*) lage(i)
!        print*,'age',lage(i),i
    if (lage(i).le.lage(i-1)) then
      write(*,*) ' #### FLEXPART MODEL ERROR! AGE CLASSES     #### '
      write(*,*) ' #### MUST BE GIVEN IN TEMPORAL ORDER.      #### '
      write(*,*) ' #### CHANGE SETTINGS IN      AGECLASSES.   #### '
      stop
    endif
  enddo
  if (lagespectra.ne.1) then
    nageclass=1
    lage(nageclass)=999999999
  endif

! jdf end read ageclasses
!*************************************************
! jdf start read available
!*************************************************
!
!*******************************************************************************
!                                                                              *
!   This routine reads the dates and times for which windfields are available. *
!                                                                              *
!     Authors: A. Stohl                                                        *
!                                                                              *
!     6 February 1994                                                          *
!     8 February 1999, Use of nested fields, A. Stohl                          *
!                                                                              *
!    12 October  2005, R. Easter -                                             *
!                      fname,wfname1,wfname1n changed from char*10 to char*80; *
!                      reads from unitavailab changed to free format           *
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
! bdate                beginning date as Julian date                           *
! beg                  beginning date for windfields                           *
! end                  ending date for windfields                              *
! fname                filename of wind field, help variable                   *
! ideltas [s]          duration of modelling period                            *
! idiff                time difference between 2 wind fields                   *
! idiffnorm            normal time difference between 2 wind fields            *
! idiffmax [s]         maximum allowable time between 2 wind fields            *
! jul                  julian date, help variable                              *
! numbwf               actual number of wind fields                            *
! wfname(maxwf)        file names of needed wind fields                        *
! wfspec(maxwf)        file specifications of wind fields (e.g., if on disc)   *
! wftime(maxwf) [s]times of wind fields relative to beginning time             *
! wfname1,wfspec1,wftime1 = same as above, but only local (help variables)     *
!                                                                              *
! Constants:                                                                   *
! maxwf                maximum number of wind fields                           *
! unitavailab          unit connected to file AVAILABLE                        *
!                                                                              *
!*******************************************************************************

! Windfields are only used, if they are within the modelling period.
! However, 1 additional day at the beginning and at the end is used for
! interpolation. -> Compute beginning and ending date for the windfields.
!************************************************************************
!hmjb beg=>dbeg    end=>dend
  if (ideltas.gt.0) then       ! forward trajectories
     dbeg=bdate-1.                  
     dend=bdate+dble(real(ideltas)/86400.)+dble(real(idiffmax)/ &
          86400.)
  else                         ! backward trajectories
     dbeg=bdate+dble(real(ideltas)/86400.)-dble(real(idiffmax)/ &
          86400.)
     dend=bdate+1.
  endif

! Open the wind field availability file and read available wind fields
! within the modelling period.
!*********************************************************************

  write(*,*) 'Opening file: ',path(3)(1:length(3)),' for reading'
  open(unitavailab,file=path(3)(1:length(3)),status='old', err=804)

  do i=1,3
    read(unitavailab,*)
  enddo
      
  numbwf=0
100     read(unitavailab,*,end=99) ldat,ltim,fname,spec
  jul=juldate(ldat,ltim)
!hmjb
  if ((jul.ge.dbeg).and.(jul.le.dend)) then
    numbwf=numbwf+1
    if (numbwf.gt.maxwf) then      ! check exceedance of dimension
      write(*,*) 'Number of wind fields needed is too great.'
      write(*,*) 'Reduce modelling period (file "COMMAND") or'
      write(*,*) 'reduce number of wind fields (file "AVAILABLE").'
      stop
    endif

    wfname1(numbwf)=fname
    wfspec1(numbwf)=spec
    wftime1(numbwf)=nint((jul-bdate)*86400.)
  endif
  goto 100       ! next wind field

99    continue

  close(unitavailab)

! Open the wind field availability file and read available wind fields
! within the modelling period (nested grids)
!*********************************************************************

  do k=1,numbnests
    write(*,*) 'Opening file: ',path(numpath+2*k)(1:length(numpath+2*k)),' for reading'
    open(unitavailab,file=path(numpath+2*(k-1)+2) &
    (1:length(numpath+2*(k-1)+2)),status='old',err=803)

    do i=1,3
      read(unitavailab,*)
    enddo

    numbwfn(k)=0
700 read(unitavailab,*,end=699) ldat,ltim,fname,spec
    jul=juldate(ldat,ltim)
!hmjb
    if ((jul.ge.dbeg).and.(jul.le.dend)) then
      numbwfn(k)=numbwfn(k)+1
      if (numbwfn(k).gt.maxwf) then      ! check exceedance of dimension
        write(*,*) 'Number of nested wind fields is too great.'
        write(*,*) 'Reduce modelling period (file "COMMAND") or'
        write(*,*) 'reduce number of wind fields (file "AVAILABLE").'
        stop
      endif

      wfname1n(k,numbwfn(k))=fname
      wfspec1n(k,numbwfn(k))=spec
      wftime1n(k,numbwfn(k))=nint((jul-bdate)*86400.)
    endif
    goto 700       ! next wind field

699 continue

  enddo
  close(unitavailab)


! Check wind field times of file AVAILABLE (expected to be in temporal order)
!****************************************************************************

  if (numbwf.eq.0) then
    write(*,*) ' #### FLEXPART MODEL ERROR! NO WIND FIELDS    #### '
    write(*,*) ' #### AVAILABLE FOR SELECTED TIME PERIOD.     #### '
    stop
  endif

  do i=2,numbwf
    if (wftime1(i).le.wftime1(i-1)) then
      write(*,*) 'FLEXPART ERROR: FILE AVAILABLE IS CORRUPT.'
      write(*,*) 'THE WIND FIELDS ARE NOT IN TEMPORAL ORDER.'
      write(*,*) 'PLEASE CHECK FIELD ',wfname1(i)
      stop
    endif
  enddo

! Check wind field times of file AVAILABLE for the nested fields
! (expected to be in temporal order)
!***************************************************************

  do k=1,numbnests
    if (numbwfn(k).eq.0) then
      write(*,*) '#### FLEXPART MODEL ERROR! NO WIND FIELDS  ####'
      write(*,*) '#### AVAILABLE FOR SELECTED TIME PERIOD.   ####'
      stop
    endif

    do i=2,numbwfn(k)
      if (wftime1n(k,i).le.wftime1n(k,i-1)) then
      write(*,*) 'FLEXTRA ERROR: FILE AVAILABLE IS CORRUPT. '
      write(*,*) 'THE NESTED WIND FIELDS ARE NOT IN TEMPORAL ORDER.'
      write(*,*) 'PLEASE CHECK FIELD ',wfname1n(k,i)
      write(*,*) 'AT NESTING LEVEL ',k
      stop
      endif
    enddo
  enddo


! For backward trajectories, reverse the order of the windfields
!***************************************************************

  do i=1,numbwf
   wfdt(i)=-999999
  enddo
  if (ideltas.ge.0) then
    do i=1,numbwf
      wfname(i)=wfname1(i)
      wfspec(i)=wfspec1(i)
      wftime(i)=wftime1(i)
    if(i.gt.1)  wfdt(i)=wftime1(i)-wftime1(i-1)  
    enddo
    do k=1,numbnests
    do i=1,numbwfn(k)
      wfnamen(k,i)=wfname1n(k,i)
      wfspecn(k,i)=wfspec1n(k,i)
      wftimen(k,i)=wftime1n(k,i)
    enddo
    enddo
  else
    do i=1,numbwf
      wfname(numbwf-i+1)=wfname1(i)
      wfspec(numbwf-i+1)=wfspec1(i)
      wftime(numbwf-i+1)=wftime1(i)
!       if(i.lt.numbwf) wfdt(numbwf-i+1)=wftime1(i+1)-wftime1(i)  
    if(i.gt.1) wfdt(numbwf-i+1)=wftime1(i)-wftime1(i-1)  
    enddo
    do k=1,numbnests
    do i=1,numbwfn(k)
      wfnamen(k,numbwfn(k)-i+1)=wfname1n(k,i)
      wfspecn(k,numbwfn(k)-i+1)=wfspec1n(k,i)
      wftimen(k,numbwfn(k)-i+1)=wftime1n(k,i)
    enddo
    enddo
  endif

! Check the time difference between the wind fields. If it is big, 
! write a warning message. If it is too big, terminate the trajectory. 
!*********************************************************************

  do i=2,numbwf
    idiff=abs(wftime(i)-wftime(i-1))
    if (idiff.gt.idiffmax) then
      write(*,*) 'FLEXPART WARNING: TIME DIFFERENCE BETWEEN TWO'
      write(*,*) 'WIND FIELDS IS TOO BIG FOR TRANSPORT CALCULATION' 
      write(*,*) 'THEREFORE, TRAJECTORIES HAVE TO BE SKIPPED.'
    else if (idiff.gt.idiffnorm) then
      write(*,*) 'FLEXPART WARNING: TIME DIFFERENCE BETWEEN TWO'
      write(*,*) 'WIND FIELDS IS BIG. THIS MAY CAUSE A DEGRADATION'
      write(*,*) 'OF SIMULATION QUALITY.'
    endif
  enddo

  do k=1,numbnests
    if (numbwfn(k).ne.numbwf) then
      write(*,*) 'FLEXTRA ERROR: THE AVAILABLE FILES FOR THE'
      write(*,*) 'NESTED WIND FIELDS ARE NOT CONSISTENT WITH'
      write(*,*) 'THE AVAILABLE FILE OF THE MOTHER DOMAIN.  '
      write(*,*) 'ERROR AT NEST LEVEL: ',k
      stop
    endif
    do i=1,numbwf
      if (wftimen(k,i).ne.wftime(i)) then
        write(*,*) 'FLEXTRA ERROR: THE AVAILABLE FILES FOR THE'
        write(*,*) 'NESTED WIND FIELDS ARE NOT CONSISTENT WITH'
        write(*,*) 'THE AVAILABLE FILE OF THE MOTHER DOMAIN.  '
        write(*,*) 'ERROR AT NEST LEVEL: ',k
        stop
      endif
    enddo
  enddo

! Reset the times of the wind fields that are kept in memory to no time
!**********************************************************************

  do i=1,2
    memind(i)=i
    memtime(i)=999999999
  enddo

  if(option_verbose.ge.10) write(*,*) "readinput() calling gridcheck()"
  call gridcheck()
  if(option_verbose.ge.10) write(*,*) "readinput() calling gridcheck_nests()"
  call gridcheck_nests()

! jdf end read available
!*************************************************
! jdf start read outgrid
!*************************************************
!
!*******************************************************************************
!                                                                              *
!     Note:  This is the FLEXPART_WRF version of subroutine readoutgrid.       *
!                                                                              *
!     This routine reads the user specifications for the output grid.          *
!                                                                              *
!     Author: A. Stohl                                                         *
!                                                                              *
!     4 June 1996                                                              *
!                                                                              *
!     Dec 2005, R. Easter -                                                    *
!             The output grid is defined by specifying its southwest and       *
!             northease corners, either in degrees-latlon or grid-meters       *
!             Changes to some read formats (wider fields).                     *
!             Changed names of "*lon0*" & "*lat0*" variables                   *
!     10 Mar 2006, R. Easter -                                                 *
!             Change eps from 1.0e-4 (degrees value) to 10.0 (meters value)    *
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
! dxout,dyout          grid distance                                           *
! numxgrid,numygrid,numzgrid    grid dimensions                                *
! out_xm0,out_ym0      lower left corner of grid                               *
! outheight(maxzgrid)  height levels of output grid [m]                        *
!                                                                              *
! Constants:                                                                   *
! unitpath             unit connected to file OUTGRID                          *
!                                                                              *
!*******************************************************************************

! Open the OUTGRID file and read output grid specifications
!**********************************************************

  write(*,*) "Reading OUTGRID"
  call skplin(1,unitpath)


! 1.  Read horizontal grid specifications
!
! *** NOTE ***
! [xtmp1, ytmp1] are the coordinates of the southwest corner
!    of the first (i.e., southwest or lower-left) output grid cell
! [xtmp2, ytmp2] are the coordinates of the northeast corner 
!    of the last (i.e,, northeast or upper-right) output grid cell
!****************************************

  !read(unitpath,'(f15.8)') xtmp1
  read(unitpath,*) xtmp1
  !read(unitpath,'(f15.8)') ytmp1
  read(unitpath,*) ytmp1
  !read(unitpath,'(2x,i7)') numxgrid
  read(unitpath,*) numxgrid
  !read(unitpath,'(2x,i7)') numygrid
  read(unitpath,*) numygrid
  !read(unitpath,'(2x,i4)') outgriddef 
  read(unitpath,*) outgriddef 
  if (outgriddef.eq.1) then 
    !read(unitpath,'(f15.8)') xtmp2
    read(unitpath,*) xtmp2
    !read(unitpath,'(f15.8)') ytmp2
    read(unitpath,*) ytmp2
  else
    !read(unitpath,'(f15.8)') dxout
    !read(unitpath,'(f15.8)') dyout
    read(unitpath,*) dxout
    read(unitpath,*) dyout
    xtmp2=dxout*real(numxgrid)+xtmp1
    ytmp2=dyout*real(numygrid)+ytmp1 
  endif

  if (option_verbose.ge. 1) then
    write(*,'(/a)') 'readoutgrid diagnostics'
    write(*,'(a,1p,2e18.10)') 'xtmp1, ytmp1 (in)', xtmp1, ytmp1
    write(*,'(a,1p,2e18.10)') 'xtmp2, ytmp2 (in)', xtmp2, ytmp2
  endif
  if (outgrid_option .eq. 1) then
! In this case, the above inputs are the actual geographical lat/lon
!   of the southwest & northeast corners of the output grid
! Need to convert from lat/lon to grid-meters
    outgrid_swlon = xtmp1
    outgrid_swlat = ytmp1

    outgrid_nelon = xtmp2
    outgrid_nelat = ytmp2
    call ll_to_xymeter_wrf( outgrid_swlon, outgrid_swlat,  &
      out_xm0, out_ym0 )
    call ll_to_xymeter_wrf( outgrid_nelon, outgrid_nelat, &
      xtmp, ytmp )
! 10-mar-2006 rce
! If out_xm0 is very close to mother grid sw corner (=xmet0), set it to that
! If xtmp    is very close to mother grid ne corner (=xmet0+nxmin1*dx), set it to that
! Do similar for out_ym0 & ytmp
    if (abs(out_xm0-xmet0) .le. eps) out_xm0 = xmet0
    if (abs(out_ym0-ymet0) .le. eps) out_ym0 = ymet0
    xr1 = xmet0 + real(nxmin1)*dx
    yr1 = ymet0 + real(nymin1)*dy
    if (abs(xtmp-xr1) .le. eps) xtmp = xr1
    if (abs(ytmp-yr1) .le. eps) ytmp = yr1
    dxout = (xtmp - out_xm0)/numxgrid
    dyout = (ytmp - out_ym0)/numygrid
    if (outgrid_option.eq.1) then ! regular
      outlat0=outgrid_swlat
      outlon0=outgrid_swlon
      dyoutl=(outgrid_nelat-outgrid_swlat)/numygrid
      dxoutl= (outgrid_nelon-outgrid_swlon)/numxgrid
    endif
  else
! In this case, the above inputs are in grid-meters 
! Need to convert from grid-meters to lat/lon
    out_xm0 = xtmp1
    out_ym0 = ytmp1
    dxout = (xtmp2 - xtmp1)/numxgrid
    dyout = (ytmp2 - ytmp1)/numygrid
    call xymeter_to_ll_wrf( xtmp1, ytmp1,  &
      outgrid_swlon, outgrid_swlat )
    call xymeter_to_ll_wrf( xtmp2, ytmp2,  &
      outgrid_nelon, outgrid_nelat )
    if (option_verbose.ge. 1) then
      write(*,'(f15.10,5x,a)') outgrid_swlon, 'outgrid_swlon'
      write(*,'(f15.10,5x,a)') outgrid_swlat, 'outgrid_swlat'
      write(*,'(f15.10,5x,a)') outgrid_nelon, 'outgrid_nelon'
      write(*,'(f15.10,5x,a)') outgrid_nelat, 'outgrid_nelat'
    endif
  end if
  if (option_verbose.ge. 1) then
    write(*,'(a,1p,2e18.10)') 'out_xm0, out_ym0 ', out_xm0, out_ym0
    write(*,'(a,1p,2e18.10)') 'dxout,   dyout   ', dxout, dyout
  endif
        
! Check validity of output grid (shall be within model domain)
!*************************************************************

  xr=out_xm0+real(numxgrid)*dxout
  yr=out_ym0+real(numygrid)*dyout
  xr1=xmet0+real(nxmin1+1)*dx
  yr1=ymet0+real(nymin1+1)*dy
  if (outgrid_option.lt.1) then
    !print*,'nx',nxmin1,nymin1
    if ((out_xm0+eps.lt.xmet0).or.(out_ym0+eps.lt.ymet0) &
        .or.(xr.gt.xr1+eps).or.(yr.gt.yr1+eps)) then
      write(*,*) ' #### FLEXPART MODEL ERROR! PART OF OUTPUT    ####'
      write(*,*) ' #### GRID IS OUTSIDE MODEL DOMAIN. CHANGE    ####'
      write(*,*) ' #### OUTGRID IN INPUT FILE                   ####'
      write(*,'(a)') path(1)(1:length(1))
      write(*,*) 'xmet0 (WRF), out_xm0 (FLX) : ',xmet0, "<", out_xm0, "?"
      write(*,*) 'ymet0 (WRF), out_ym0 (FLX) : ',ymet0, "<", out_ym0, "?"
      write(*,*) 'xr1 (WRF),   xr (FLX)      : ',xr1, ">", xr, "?"
      write(*,*) 'yr1 (WRF),   yr (FLX)      : ',yr1, ">", yr, "?"
      stop
    endif
  endif
!      if ((numxgrid.gt.maxxgrid).or.(numygrid.gt.maxygrid)) then
!        write(*,*) ' #### FLEXPART MODEL ERROR! DIMENSIONS OF     ####'
!        write(*,*) ' #### OUTPUT GRID EXCEED MAXIMUM VALUES.      ####'
!        write(*,*) ' #### CHANGE FILE $FLEXPART/options/OUTGRID.  ####'
!        stop
!      endif

! 2. Vertical levels of output grid
!**********************************
      
  read(unitpath,*) numzgrid
!        if (numzgrid.gt.maxzgrid) then
!       write(*,*) ' #### FLEXPART MODEL ERROR! TOO MANY HEIGHT   #### ' 
!       write(*,*) ' #### LEVELS ARE GIVEN FOR OUTPUT GRID.       #### ' 
!       write(*,*) ' #### MAXIMUM NUMBER IS ',maxzgrid,'          #### ' 
!       write(*,*) ' #### PLEASE MAKE CHANGES IN FILE OUTGRID.    #### ' 
!        stop
!        endif

  allocate(outheight(numzgrid), stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate outheight'
  
  allocate(outheighthalf(numzgrid), stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate outheighthalf'

  allocate(oroout(0:numxgrid-1,0:numygrid-1), stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate oroout'
  
  allocate(area(0:numxgrid-1,0:numygrid-1), stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate area'
  
  allocate(volume(0:numxgrid-1,0:numygrid-1,numzgrid), stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate volume'
  
  allocate(areaeast(0:numxgrid-1,0:numygrid-1,numzgrid), stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate areaeast'
  
  allocate(areanorth(0:numxgrid-1,0:numygrid-1,numzgrid), stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate areanorth'

  do j=1,numzgrid
    read(unitpath,*) outhelp
    outheight(j)=outhelp
  enddo


! Check whether vertical levels are specified in ascending order
!***************************************************************

  do j=2,numzgrid
    if (outheight(j).le.outheight(j-1)) then
      write(*,*) ' #### FLEXPART MODEL ERROR! YOUR SPECIFICATION#### ' 
      write(*,*) ' #### OF OUTPUT LEVELS IS CORRUPT AT LEVEL    #### ' 
      write(*,*) ' #### ',j,'                              #### ' 
      write(*,*) ' #### PLEASE MAKE CHANGES IN  OUTGRID.        #### ' 
    endif
  enddo

! Determine the half levels, i.e. middle levels of the output grid
!*****************************************************************

  outheighthalf(1)=outheight(1)/2.
  do j=2,numzgrid
    outheighthalf(j)=(outheight(j-1)+outheight(j))/2.
  enddo

  xoutshift=xmet0-out_xm0
  youtshift=ymet0-out_ym0

  if(nested_output .eq. 1) then
    write(*,*) "Reading OUTGRID nest"
    call skplin(1,unitpath)
    !read(unitpath,'(f15.8)') xtmp1
    !read(unitpath,'(f15.8)') ytmp1
    !read(unitpath,'(4x,i5)') numxgridn
    !read(unitpath,'(4x,i5)') numygridn
    !read(unitpath,'(2x,i4)') outgriddefn
    read(unitpath,*) xtmp1
    read(unitpath,*) ytmp1
    read(unitpath,*) numxgridn
    read(unitpath,*) numygridn
    read(unitpath,*) outgriddefn
    if (outgriddefn.eq.1) then
      !read(unitpath,'(f15.8)') xtmp2
      !read(unitpath,'(f15.8)') ytmp2
      read(unitpath,*) xtmp2
      read(unitpath,*) ytmp2
    else
      !read(unitpath,'(f15.8)') dxoutn
      !read(unitpath,'(f15.8)') dyoutn
      read(unitpath,*) dxoutn
      read(unitpath,*) dyoutn
      xtmp2=dxoutn*real(numxgridn)+xtmp1
      ytmp2=dyoutn*real(numygridn)+ytmp1
    endif

    write(*,'(/a)') 'readoutgrid_nest diagnostics'
    write(*,'(a,1p,2e18.10)') 'xtmp1, ytmp1  (in)', xtmp1, ytmp1
    write(*,'(a,1p,2e18.10)') 'xtmp2, ytmp2  (in)', xtmp2, ytmp2
    if (outgrid_option .eq. 1) then
      ! In this case, the above inputs are the actual geographical lat/lon
      !   of the southwest & northeast corners of the output grid
      ! Need to convert from lat/lon to grid-meters
      outgridn_swlon = xtmp1
      outgridn_swlat = ytmp1
      outgridn_nelon = xtmp2
      outgridn_nelat = ytmp2
      call ll_to_xymeter_wrf( outgridn_swlon, outgridn_swlat, &
        out_xm0n, out_ym0n )
      call ll_to_xymeter_wrf( outgridn_nelon, outgridn_nelat, &
        xtmp, ytmp )
      dxoutn = (xtmp - out_xm0n)/numxgridn
      dyoutn = (ytmp - out_ym0n)/numygridn
      if (outgrid_option.eq.1) then ! regular
        outlat0n=outgridn_swlat
        outlon0n=outgridn_swlon
        dyoutln=(outgridn_nelat-outgridn_swlat)/numygridn
        dxoutln= (outgridn_nelon-outgridn_swlon)/numxgridn
      endif

    else
      ! In this case, the above inputs are in grid-meters 
      ! Need to convert from grid-meters to lat/lon
      out_xm0n = xtmp1
      out_ym0n = ytmp1
      dxoutn = (xtmp2 - xtmp1)/numxgridn
      dyoutn = (ytmp2 - ytmp1)/numygridn
      call xymeter_to_ll_wrf( xtmp1, ytmp1,  &
        outgridn_swlon, outgridn_swlat )
      call xymeter_to_ll_wrf( xtmp2, ytmp2,  &
        outgridn_nelon, outgridn_nelat )
      write(*,'(f15.10,5x,a)') outgridn_swlon, 'outgridn_swlon'
      write(*,'(f15.10,5x,a)') outgridn_swlat, 'outgridn_swlat'
      write(*,'(f15.10,5x,a)') outgridn_nelon, 'outgridn_nelon'
      write(*,'(f15.10,5x,a)') outgridn_nelat, 'outgridn_nelat'
    endif
    write(*,'(a,1p,2e18.10)') 'out_xm0n, out_ym0n',out_xm0n,out_ym0n
    write(*,'(a,1p,2e18.10)') 'dxoutn,   dyoutn  ',dxoutn,dyoutn


  ! Check validity of output grid (shall be within model domain)
  !*************************************************************

    xr=out_xm0n+real(numxgridn)*dxoutn
    yr=out_ym0n+real(numygridn)*dyoutn
    xr1=xmet0+real(nxmin1+1)*dx
    yr1=ymet0+real(nymin1+1)*dy
    if ((out_xm0n+eps.lt.xmet0).or.(out_ym0n+eps.lt.ymet0) &
      .or.(xr.gt.xr1+eps).or.(yr.gt.yr1+eps)) then
      write(*,*) ' #### FLEXPART MODEL ERROR! PART OF OUTPUT    ####'
      write(*,*) ' #### NEST IS OUTSIDE MODEL DOMAIN. CHANGE    ####'
      write(*,*) ' ####      OUTGRID                            ####'
      write(*,'(a)') path(1)(1:length(1))
      stop
    endif
  !        if ((numxgridn.gt.maxxgridn).or.(numygridn.gt.maxygridn)) then
  !        write(*,*) ' #### FLEXPART MODEL ERROR! DIMENSIONS OF     ####'
  !        write(*,*) ' #### OUTPUT NEST EXCEED MAXIMUM VALUES.      ####'
  !        write(*,*) ' #### CHANGE FILE $FLEXPART/options/OUTGRID.  ####'
  !        stop
  !        endif
    xoutshiftn=xmet0-out_xm0n
    youtshiftn=ymet0-out_ym0n
    allocate(orooutn(0:numxgridn-1,0:numygridn-1), stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate orooutn'
    
    allocate(arean(0:numxgridn-1,0:numygridn-1), stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate arean'
    
    allocate(volumen(0:numxgridn-1,0:numygridn-1,numzgrid), stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate volumen'
  endif

! jdf end read outgrid
!*************************************************
! jdf start read receptors
!*************************************************
!
!*******************************************************************************
!                                                                              *
!     Note:  This is the FLEXPART_WRF version of subroutine assignland.        *
!            The computational grid is the WRF x-y grid rather than lat-lon.   *
!                                                                              *
!     This routine reads the user specifications for the receptor points.      *
!                                                                              *
!     Author: A. Stohl                                                         *
!                                                                              *
!     1 August 1996                                                            *
!                                                                              *
!     Oct 2005, R. Easter - change calc of receptorarea()                      *
!     Dec 2005, R. Easter - x/yrecptor values may be input either as           *
!                           degrees-latlon or grid-meters                      *
!                           Changes to some read formats (wider fields).       *
!                           Changed names of "*lon0*" & "*lat0*" variables     *
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
! receptorarea(maxreceptor)  area of dx*dy at location of receptor             *
! receptorname(maxreceptor)  names of receptors                                *
! xreceptor,yreceptor  coordinates of receptor points                          *
!                                                                              *
! Constants:                                                                   *
! unitreceptor         unit connected to file RECEPTORS                        *
!                                                                              *
!*******************************************************************************

  write(*,*) "Reading RECEPTORS"
  call skplin(1,unitpath)
  read(unitpath,*) numreceptor
! For backward runs, do not allow receptor output. Thus, set number of receptors to zero
!***************************************************************************************

  if (ldirect.lt.0) then
    numreceptor=0
  endif

! Open the RECEPTORS file and read output grid specifications
!************************************************************
! Read the names and coordinates of the receptors
!************************************************

  if (numreceptor.gt.maxreceptor) then
    write(*,*) ' #### FLEXPART MODEL ERROR! TOO MANY RECEPTOR #### ' 
    write(*,*) ' #### POINTS ARE GIVEN.                       #### ' 
    write(*,*) ' #### MAXIMUM NUMBER IS ',maxreceptor,'       #### ' 
    write(*,*) ' #### PLEASE MAKE CHANGES IN      RECEPTORS   #### ' 
    stop
  endif
  do j=1,numreceptor
    read(unitpath,'(4x,a16)') receptor
    read(unitpath,*) x
    read(unitpath,*) y
    receptorname(j)=receptor

    write(*,'(/a,i5)') 'readreceptor diagnostics, j =', j
    write(*,'(a,1p,2e18.10)') 'x, y (in)   ', x, y
    if (numpoint_option .eq. 1) then
      ! In this case, the above inputs are the actual geographical lat/lon
      ! Need to convert from lat/lon to grid-index coordinates
      receptor_lon(j) = x
      receptor_lat(j) = y
      call ll_to_xyindex_wrf( receptor_lon(j), receptor_lat(j), &
        xreceptor(j), yreceptor(j) )
    else
      ! In this case, the above inputs are in grid-meters 
      ! Need to convert from grid-meters to grid-index coordinates, then to lat/lon
      xreceptor(j)=(x-xmet0)/dx
      yreceptor(j)=(y-ymet0)/dy
      call xyindex_to_ll_wrf( 0, xreceptor(j), yreceptor(j), &
        receptor_lon(j), receptor_lat(j) )
      write(*,'(f15.10,5x,a)') receptor_lon(j), 'receptor_lon'
      write(*,'(f15.10,5x,a)') receptor_lat(j), 'receptor_lat'
    end if
    write(*,'(a,1p,2e18.10)') 'x, yreceptor',  &
       xreceptor(j), yreceptor(j)
    ! for FLEXPART_WRF, dx & dy are in meters,
    ! receptorarea() appears to be area in m**2 of a mother grid cell
    ! centers on x,y
    !       xm=r_earth*cos(y*pi/180.)*dx/180.*pi
    !       ym=r_earth*dy/180.*pi
    xm=dx
    ym=dy
    receptorarea(j)=xm*ym
  enddo

! jdf end read receptors
!*************************************************
! jdf start read species
!*************************************************
!
!*******************************************************************************
!                                                                              *
!     This routine reads names and physical constants of chemical species/     *
!     radionuclides available with FLEXPART.                                   *
!                                                                              *
!     Author: A. Stohl                                                         *
!                                                                              *
!     11 July 1996                                                             *
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
! decaytime(maxtable)  half time for radiological decay                        *
! specname(maxtable)   names of chemical species, radionuclides                *
! wetscava, wetscavb   Parameters for determining scavenging coefficient       *
!                                                                              *
! Constants:                                                                   *
!                                                                              *
!*******************************************************************************

! Open the SPECIES file and read species names and properties
!************************************************************

  write(*,*) 'Reading SPECIES'
  read(unitpath,*)
  read(unitpath,*) numtable
  read(unitpath,*)
  if(numtable.gt.maxspec) then
    write(*,*) 'Number entries in SPECIES table (numtable) exceeds maxspec'
    write(*,*) 'Reduce numtable and adjust number of entries accordingly OR'
    write(*,*) 'increase maxspec in par_mod.f90 and re-compile flexpart-wrf'
    stop
  endif

  do i=1,numtable
    read(unitpath,21) specname(i),decaytime(i), &
    wetscava(i),wetscavb(i),drydiff(i),dryhenry(i),dryactiv(i), &
    partrho(i),partmean(i),partsig(i),dryvelo(i),weightmol(i), &
    ohreact(i),spec_ass(i),kao(i)
    !read(unitpath,21) specname(i),decay(i), &
    !weta(i),wetb(i),reldiff(i),henry(i),f0(i), &
    !density(i),dquer(i),dsigma(i),dryvel(i),weightmolar(i), &
    !ohreact(i),spec_ass(i),kao(i)
    pos_spec=i

    if ((wetscava(i).gt.0.).and.(dryhenry(i).le.0.)) then
      if (partmean(i).le.0.) goto 996 ! no particle, no henry set
    endif

    if (spec_ass(i).gt.0) then
      spec_found=.FALSE.
      do j=1,pos_spec-1
        if (spec_ass(pos_spec).eq.specnum(j)) then
          spec_ass(pos_spec)=j
          spec_found=.TRUE.
          ASSSPEC=.TRUE.
        endif
      end do
      if (spec_found.eqv..FALSE.) then
        goto 997
      endif
    endif

  write(*,*) "option_verbose:",option_verbose
    if (partsig(i).eq.1.) partsig(i)=1.0001   ! avoid realing exception
    if (partsig(i).eq.0.) partsig(i)=1.0001   ! avoid realing exception

    if ((drydiff(i).gt.0.).and.(partrho(i).gt.0.)) then
      write(*,*) '#### FLEXPART MODEL ERROR! SPECIES FORMAT    ####'
      write(*,*) '#### IS CORRUPT. SPECIES CANNOT BE BOTH      ####'
      write(*,*) '#### PARTICLE AND GAS.                       ####'
      stop
    endif
  enddo

21 format(4x,a10,f10.1,e11.1,f6.2,f7.1,e9.1,f5.1,e10.1,2e8.1,2f8.2, &
          e18.1,i18,f18.2)

! jdf end read species
!*************************************************
! jdf start read releases
!*************************************************
!
!*******************************************************************************
!                                                                              *
!     This routine reads the release point specifications for the current      *
!     model run. Several release points can be used at the same time.          *
!                                                                              *
!     Note:  This is the FLEXPART_WRF version of subroutine readreleases.      *
!            The computational grid is the WRF x-y grid rather than lat-lon.   *
!                                                                              *
!     Author: A. Stohl                                                         *
!     18 May 1996                                                              *
!                                                                              *
!     Update: 29 January 2001                                                  *
!     Release altitude can be either in magl or masl                           *
!                                                                              *
!     Nov 2005, R. Easter - Do not adjust xpoint1 & 2 by +/-360 degrees        *
!     Dec 2005, R. Easter - x/ypoint1/2 values may be input either as          *
!                           degrees-latlon or grid-meters                      *
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
! decay               decay constant of species                                *
! dquer [um]          mean particle diameters                                  *
! dsigma              e.g. dsigma=10 or dsigma=0.1 means that 68% of the mass  *
!                     are between 0.1*dquer and 10*dquer                       *
! ireleasestart, ireleaseend [s] starting time and ending time of each release *
! kindz               1: zpoint is in m agl, 2: zpoint is in m asl, 3: zpoint  *
!                     is in hPa                                                *
! npart               number of particles to be released                       *
! nspec               number of species to be released                         *
! density [kg/m3]     density of the particles                                 *
! rm [s/m]            Mesophyll resistance                                     *
! species             name of species                                          *
! xmass               total mass of each species                               *
! xpoint1,ypoint1     geograf. coordinates of lower left corner of release area*
! xpoint2,ypoint2     geograf. coordinates of upper right corner of release are*
! weta, wetb          parameters to determine the wet scavenging coefficient   *
! zpoint1,zpoint2     height range, over which release takes place             *
!                                                                              *
!*******************************************************************************

  DEP=.false.
  DRYDEP=.false.
  WETDEP=.false.
  do i=1,maxspec
    DRYDEPSPEC(i)=.false.
  enddo

! Open the releases file and read user options
!*********************************************

! Check the format of the RELEASES file (either in free format,
! or using a formatted mask)
! Use of formatted mask is assumed if line 10 contains the word 'DIRECTION'
!**************************************************************************
  write(*,*) 'Reading RELEASES'
  call skplin(1,unitpath)

! Read the number of species and the link to the species information table
! Assign species-specific parameters needed for physical processes
!*************************************************************************

  read(unitpath,*) nspec
  if (nspec.gt.maxspec) then 
    write(*,*) '#####################################################'
    write(*,*) '#### FLEXPART MODEL SUBROUTINE READRELEASES:     ####'
    write(*,*) '####                                             ####'
    write(*,*) '#### ERROR - MAXIMUM NUMBER OF EMITTED SPECIES IS####'
    write(*,*) '#### TOO LARGE. PLEASE REDUCE NUMBER OF SPECIES. ####'
    write(*,*) '#####################################################'
    stop
  endif
  if (option_verbose.ge.1) write(*,*) "nspec = ",nspec
  read(unitpath,*) emitvar
  if (option_verbose.ge.1) write(*,*) "emitvar = ",emitvar
  do i=1,nspec
    if (option_verbose.ge.10) write(*,*) "Reading species:",i,"of",nspec
    read(unitpath,*) link(i)
    species(i)=specname(link(i))

! For backward runs, only 1 species is allowed
!*********************************************

    if ((ldirect.lt.0).and.(nspec.gt.1)) then
      write(*,*) '#####################################################'
      write(*,*) '#### FLEXPART MODEL SUBROUTINE READRELEASES:     ####'
      write(*,*) '#### FOR BACKWARD RUNS, ONLY 1 SPECIES IS ALLOWED####'
      write(*,*) '#####################################################'
      stop
    endif 

! Molecular weight
!*****************

    weightmolar(i)=weightmol(link(i))
    if (((iout.eq.2).or.(iout.eq.3)).and. &
    (weightmolar(i).lt.0.)) then
      write(*,*) 'For mixing ratio output, valid molar weight'
      write(*,*) 'must be specified for all simulated species.'
      write(*,*) 'Check table SPECIES or choose concentration'
      write(*,*) 'output instead if molar weight is not known.'
      stop
    endif


! Radioactive decay
!******************

    decay(i)=0.693147/decaytime(link(i)) !conversion half life to decay constant

    ! Wet deposition
    !***************

    weta(i)=wetscava(link(i))
    wetb(i)=wetscavb(link(i))

    ! Dry deposition of gases
    !************************

    reldiff(i)=drydiff(link(i))             ! Diffusivity rel. to H20
    henry(i)=dryhenry(link(i))              ! Henry constant
    f0(i)=dryactiv(link(i))                 ! activity
    if (reldiff(i).gt.0.) &
      rm(i)=1./(henry(i)/3000.+100.*f0(i))    ! mesophyll resistance

    ! Dry deposition of particles
    !****************************

    vsetaver(i)=0.
    density(i)=partrho(link(i))                 ! Particle density
    dquer(i)=partmean(link(i))*1000000.         ! Conversion m to um
    dsigma(i)=partsig(link(i))
    if (density(i).gt.0.) then                  ! Additional parameters
!          call part0(dquer(i),dsigma(i),density(i),fracth,schmih,vsh)
      call part0(dquer(i),dsigma(i),density(i),fracth,schmih,cun,vsh)
      
      do j=1,ni
        fract(i,j)=fracth(j)
        schmi(i,j)=schmih(j)
        vset(i,j)=vsh(j)
        cunningham(i)=cunningham(i)+cun*fract(i,j)
        vsetaver(i)=vsetaver(i)-vset(i,j)*fract(i,j)
      enddo
    endif

! Dry deposition for constant deposition velocity
!************************************************

    dryvel(i)=dryvelo(link(i))*0.01         ! conversion to m/s

    if (weta(i).gt.0.) WETDEP=.true.
    if ((reldiff(i).gt.0.).or.(density(i).gt.0.).or. &
        (dryvel(i).gt.0.)) then
      DRYDEP=.true.
      DRYDEPSPEC(i)=.true.
    endif


! Read in daily and day-of-week variation of emissions, if available
!*******************************************************************

    do j=1,24           ! initialize everything to no variation
      area_hour(i,j)=1.
      point_hour(i,j)=1.
    enddo
    do j=1,7
      area_dow(i,j)=1.
      point_dow(i,j)=1.
    enddo

!       write(aspecnumb,'(i3.3)') link(i)
!       open(unitemissvar,file=path(1)(1:len(1))//'EMISSION_VARIATION_'
!    +  //aspecnumb//'.dat',status='old',err=35)
!       read(unitemissvar,*)
    if(emitvar.eq.1) then
      do j=1,24     ! 24 hours, starting with 0-1 local time
        read(unitpath,*) ihour,area_hour(i,j),point_hour(i,j)
      enddo
      !       read(unitemissvar,*)
      do j=1,7      ! 7 days of the week, starting with Monday
        read(unitpath,*) idow,area_dow(i,j),point_dow(i,j)
      enddo
      !close(unitemissvar)
    endif
    !35      continue
  enddo
  if (WETDEP.or.DRYDEP) DEP=.true.


! Read specifications for each release point
!*******************************************

  numpoint=0
  numpartmax=0
  releaserate=0.
  releaserate2=0.
  read(unitpath,*) numpoint

  !     numpoint2=numpoint+500
  numpoint2=numpoint+0
  allocate(ireleasestart(numpoint2), stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate RELEASPOINT'
  
  allocate(ireleaseend(numpoint2), stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate RELEASPOINT'
  
  allocate(xpoint1(numpoint2), stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate RELEASPOINT'

  allocate(xpoint12(numpoint2), stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate RELEASPOINT'

  allocate(xpoint22(numpoint2), stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate RELEASPOINT'

  allocate(ypoint12(numpoint2), stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate RELEASPOINT'

  allocate(ypoint22(numpoint2), stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate RELEASPOINT'

  allocate(releases_swlon(numpoint2), stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate RELEASPOINT'

  allocate(releases_swlat(numpoint2), stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate RELEASPOINT'

  allocate(releases_nelon(numpoint2), stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate RELEASPOINT'

  allocate(releases_nelat(numpoint2), stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate RELEASPOINT'

  allocate(xpoint2(numpoint2), stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate RELEASPOINT'

  allocate(ypoint1(numpoint2), stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate RELEASPOINT'

  allocate(ypoint2(numpoint2), stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate RELEASPOINT'

  allocate(zpoint1(numpoint2), stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate RELEASPOINT'

  allocate(zpoint2(numpoint2), stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate RELEASPOINT'

  allocate(kindz(numpoint2), stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate RELEASPOINT'

  allocate(xmass(numpoint2,maxspec), stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate RELEASPOINT'

  allocate(rho_rel(numpoint2), stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate RELEASPOINT'

  allocate(npart(numpoint2), stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate RELEASPOINT'
  !     print*,'allocate xmassa',numpoint

  allocate(xmasssave(numpoint2), stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate RELEASPOINT'

  if (option_verbose.ge. 1) then
    write (*,*) ' Releasepoints allocated: ', numpoint
  endif
  do i=1,numpoint
    xmasssave(i)=0.
  end do



!      if (numpoint.gt.maxpoint) goto 997
  if(option_verbose.ge.1) write(*,*) 'numpoint =',numpoint
  do j=1,numpoint
    read(unitpath,*) id1,it1
    read(unitpath,*) id2,it2
    read(unitpath,*) xpoint1(j)
    read(unitpath,*) ypoint1(j)
    read(unitpath,*) xpoint2(j)
    read(unitpath,*) ypoint2(j)
    read(unitpath,*) kindz(j)
    read(unitpath,*) zpoint1(j)
    read(unitpath,*) zpoint2(j)
    read(unitpath,*) npart(j)
    do i=1,nspec
      read(unitpath,*) xmass(j,i)
    enddo
    nparttot=nparttot+npart(j)
    if (j.le.2000) then
      read(unitreleases,'(a40)',err=998) compoint(j)(1:40)
    else
      read(unitreleases,'(a40)',err=998) compoint(2001)(1:40)
    endif

    !read(unitpath,'(a20)',err=998) compoint(j)(1:20)

    if(option_verbose.ge.1) print*,'release location=',j
    
    if(option_verbose.ge.1) then ! AD: If we have many points, this extra info is extremely useful!
      if (j.le.2000) then 
        write(*,*) 'COMPOINT = "',compoint(j)(1:40)
      else
        write(*,*) 'COMPOINT = "',compoint(2001)(1:40)
      endif
    endif
        
!     write(*,'(/a,i7)') 'readreleases diagnostics - numpoint = ', j
!     write(*,'(a,1p,2e18.10)') 'x, ypoint1 (in) ',  &
!        xpoint1(j), ypoint1(j)
!     write(*,'(a,1p,2e18.10)') 'x, ypoint2 (in) ',  &
!        xpoint2(j), ypoint2(j)
! JB
! In this case, the above inputs are the actual geographical lat/lon
!   of the southwest & northeast corners of the release area
! Need to convert from lat/lon to grid-meters
!        if (outgrid_option.eq.1) then !regular
    if (numpoint_option.eq.1) then !regular
      xpoint12(j)=xpoint1(j)
      xpoint22(j)=xpoint2(j)
      ypoint12(j)=ypoint1(j)
      ypoint22(j)=ypoint2(j)

      releases_swlon(j) = xpoint1(j)
      releases_swlat(j) = ypoint1(j)
      releases_nelon(j) = xpoint2(j)
      releases_nelat(j) = ypoint2(j)
      call ll_to_xymeter_wrf( releases_swlon(j), releases_swlat(j),  &
        xpoint1(j), ypoint1(j) )
      call ll_to_xymeter_wrf( releases_nelon(j), releases_nelat(j), &
        xpoint2(j), ypoint2(j) )

    else
!        write(*,'(a,1p,2e18.10)') 'x, ypoint1      ', 
!    &      xpoint1(j), ypoint1(j)
!        write(*,'(a,1p,2e18.10)') 'x, ypoint2      ', 
!    &      xpoint2(j), ypoint2(j)
!      else
! In this case, the above inputs are in grid-meters 
! Need to convert from grid-meters to lat/lon
      call xymeter_to_ll_wrf( xpoint1(j), ypoint1(j), &
        releases_swlon(j), releases_swlat(j) )
      call xymeter_to_ll_wrf( xpoint2(j), ypoint2(j), &
        releases_nelon(j), releases_nelat(j) )
    endif
!         write(*,'(f15.10,5x,a)') releases_swlon(j), 'releases_swlon'
!         write(*,'(f15.10,5x,a)') releases_swlat(j), 'releases_swlat'
!         write(*,'(f15.10,5x,a)') releases_nelon(j), 'releases_nelon'
!         write(*,'(f15.10,5x,a)') releases_nelat(j), 'releases_nelat'
!      end if

! If a release point contains no particles, stop and issue error message
!***********************************************************************

    if (npart(j).eq.0) then
      write(*,*) 'FLEXPART MODEL ERROR'
      write(*,*) 'RELEASES file is corrupt.'
      write(*,*) 'At least for one release point, there are zero'
      write(*,*) 'particles released. Make changes to RELEASES.'
      stop
    endif

! Check whether x coordinates of release point are within model domain
!*********************************************************************

! FLEXPART_WRF - x & y coords are in meters, so the following lines 
!   (which adjust longitude by +/-360 degrees) are not needed
!
!      if (xpoint1(numpoint).lt.xlon0) 
!    +       xpoint1(numpoint)=xpoint1(numpoint)+360.
!      if (xpoint1(numpoint).gt.xlon0+(nxmin1)*dx)
!    +       xpoint1(numpoint)=xpoint1(numpoint)-360.
!      if (xpoint2(numpoint).lt.xlon0) 
!    +       xpoint2(numpoint)=xpoint2(numpoint)+360.
!      if (xpoint2(numpoint).gt.xlon0+(nxmin1)*dx)
!    +       xpoint2(numpoint)=xpoint2(numpoint)-360.

! Determine relative beginning and ending times of particle release
!******************************************************************

    jul1=juldate(id1,it1)
    jul2=juldate(id2,it2)
    if (jul1.gt.jul2) then
      write(*,*) 'FLEXPART MODEL ERROR'
      write(*,*) 'Release stops before it begins.'
      write(*,*) 'Make changes to file RELEASES.'
      stop
    endif
    if (mdomainfill.eq.0) then   ! no domain filling
      if (ldirect.eq.1) then
        if ((jul1.lt.bdate).or.(jul2.gt.edate)) then
          write(*,*) 'FLEXPART MODEL ERROR'
          write(*,*) 'Release starts before simulation begins or ends'
          write(*,*) 'after simulation stops.'
          write(*,*) 'Make files COMMAND and RELEASES consistent.'
          stop
        endif
        ireleasestart(j)=int((jul1-bdate)*86400.)
        ireleaseend(j)=int((jul2-bdate)*86400.)
      else if (ldirect.eq.-1) then
        if ((jul1.lt.edate).or.(jul2.gt.bdate)) then
          write(*,*) 'FLEXPART MODEL ERROR'
          write(*,*) 'Release starts before simulation begins or ends'
          write(*,*) 'after simulation stops.'
          write(*,*) 'Make files COMMAND and RELEASES consistent.'
          stop
        endif
        ireleasestart(j)=int((jul1-bdate)*86400.)
        ireleaseend(j)=int((jul2-bdate)*86400.)
      endif
    endif


! Check, whether the total number of particles may exceed totally allowed
! number of particles at some time during the simulation
!************************************************************************

! Determine the release rate (particles per second) and total number
! of particles released during the simulation
!*******************************************************************

    if (ireleasestart(j).ne.ireleaseend(j)) then
      releaserate=releaserate+real(npart(j))/ &
      real(ireleaseend(j)-ireleasestart(j))
      !print*,'release',ireleaseend(j),ireleasestart(j)
    else
      releaserate=99999999.
    endif
    if (ireleaseend(j)-ireleasestart(j).lt.lage(nageclass)) then
      releaserate2=releaserate2+real(npart(j))
    else
      releaserate2=releaserate2+real(npart(j))*real(lage(nageclass))/ &
        real(ireleaseend(j)-ireleasestart(j))
    endif
    numpartmax=numpartmax+npart(j)

    if (ioutputforeachrelease.eq.1) then
      maxpointspec_act=numpoint
    else
      maxpointspec_act=1
    endif

  enddo

  maxpart=int(releaserate2*1.1)
!      print*,'maxpart',maxpart,releaserate2,nparttot
!     maxpart=40000000

!     print*,'numpoint',numpoint
!     print*,'NPARTTOT',nparttot,lage(4),lage(5),lage(6),nageclass
!     print*,'maxpart',maxpart,real(lage(nageclass)),releaserate,releaserate2,real(lage(nageclass))*1.1*releaserate
!      nparttot2=nparttot+500000
!      nparttot2=nparttot

  if (ipin.eq.1) then
    if (option_verbose.ge.1) print*, "Attempting to read previously dumped particles"
    if (iouttype .eq. 0 .or. iouttype .eq.2) then
      open(unitpartin,file=path(1)(1:length(1))//'partposit_end', &
        form='unformatted',err=998)
    else
      open(unitpartin,file=path(1)(1:length(1))//'partposit_end', &
        form='formatted',err=998)
    endif
    if(option_verbose .ge. 1) print*,'numpart_in',numpart_in
    if (iouttype .eq. 0 .or. iouttype .eq.2) then
      read(unitpartin,end=101) itimein,numpart_in, iomode_xycoord_in  !AD: changed label to 101
    else
      read(unitpartin,*,end=101) itimein,numpart_in, iomode_xycoord_in  !AD: changed label to 101
    endif
101 close(unitpartin) ! AD: close when end of file is reached
    if(option_verbose.ge.1) print*, "Finished reading dumped particles"
  !else
  !  numpart_in=0	# Removed in merge with JB
  endif
  if(option_verbose .ge. 1) print*,'maxpart2',maxpart,numpart_in
  nparttot2=maxpart+numpart_in
  maxpart=nparttot2
  if( option_verbose .ge. 1) print*,'alloc part',nparttot2,maxpart,numpart_in
!      print*,'alloc part',nparttot2,maxpart,numpart_in

  allocate(xmass1(nparttot2,nspec), stat=stat)
  !allocate(drydep1(nparttot2,nspec),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate xmass1'

  allocate(itra1(nparttot2), stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate itra1'

  allocate(npoint(nparttot2), stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate npoint'

  allocate(nclass(nparttot2), stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate nclass'

  allocate(idt(nparttot2), stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate idt   '

  allocate(itramem(nparttot2), stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate itrame'

  allocate(itrasplit(nparttot2), stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate itrasp'

  allocate(xtra1(nparttot2), stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate xtra1 '

  allocate(ytra1(nparttot2), stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate ytra1 '

  allocate(ztra1(nparttot2), stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate ztra1 '

  do j=1,nparttot2
    itra1(j)=-999999999
  enddo


  if ( releaserate.gt. 0.99*real(maxpart)/real(lage(nageclass)) ) then
    !if (numpartmax.gt.maxpart) then
    !  write(*,*) '#####################################################'
    !  write(*,*) '#### FLEXPART MODEL SUBROUTINE READRELEASES:     ####'
    !  write(*,*) '####                                             ####'
    !  write(*,*) '####WARNING - TOTAL NUMBER OF PARTICLES SPECIFIED####'
    !  write(*,*) '#### IN FILE "RELEASES" MAY AT SOME POINT DURING ####'
    !  write(*,*) '#### THE SIMULATION EXCEED THE MAXIMUM ALLOWED   ####'
    !  write(*,*) '#### NUMBER (MAXPART).IF RELEASES DO NOT OVERLAP,####'
    !  write(*,*) '#### FLEXPART CAN POSSIBLY COMPLETE SUCCESSFULLY.####'
    !  write(*,*) '#### HOWEVER, FLEXPART MAY HAVE TO STOP          ####'
    !  write(*,*) '#### AT SOME TIME DURING THE SIMULATION. PLEASE  ####'
    !  write(*,*) '#### MAKE SURE THAT YOUR SETTINGS ARE CORRECT.   ####'
    !  write(*,*) '#####################################################'
    !  write(*,*) 'Maximum release rate may be: ',releaserate,' particles per second'
    !  write(*,*) 'Maximum allowed release rate is: ', &
    !    real(maxpart)/real(lage(nageclass)),' particles per second'
    !  write(*,*) &
    !    'Total number of particles released during the simulation is: ', &
    !     numpartmax
    !   write(*,*) 'Maximum allowed number of particles is: ',maxpart
    !endif
  endif
! Make a consistency check, whether the forward/backward switch is correctly set
!*******************************************************************************

  !if (ldirect.eq.1) then
  !  if (maxpointspec.lt.nspec) then
  !    write(*,*) '#####################################################'
  !    write(*,*) '#### FLEXPART MODEL SUBROUTINE READRELEASES:     ####'
  !    write(*,*) '####                                             ####'
  !    write(*,*) '#### ERROR - PARAMETER MAXPOINTSPEC IS NOT       ####'
  !    write(*,*) '#### CORRECTLY SET FOR A FORWARD SIMULATION.     ####'
  !    write(*,*) '#### CHANGE APPROPRIATELY IN FILE INCLUDEPAR.    ####'
  !    write(*,*) '#####################################################'
  !  endif
  !else
  !  if (maxpointspec.lt.numpoint) then
  !    write(*,*) '#####################################################'
  !    write(*,*) '#### FLEXPART MODEL SUBROUTINE READRELEASES:     ####'
  !    write(*,*) '####                                             ####'
  !    write(*,*) '#### ERROR - PARAMETER MAXPOINTSPEC IS NOT       ####'
  !    write(*,*) '#### CORRECTLY SET FOR A BACKWARD SIMULATION.    ####'
  !    write(*,*) '#### CHANGE APPROPRIATELY IN FILE INCLUDEPAR.    ####'
  !    write(*,*) '#####################################################'
  !  endif
  !endif

  return

! jdf end read releases
!*************************************************

997   write(*,*) '#####################################################'
      write(*,*) '#### FLEXPART MODEL SUBROUTINE READRELEASES:     ####'
      write(*,*) '####                                             ####'
      write(*,*) '#### ERROR - NUMBER OF RELEASE POINTS SPECIFIED  ####'
      write(*,*) '#### IN      "RELEASES" EXCEEDS THE MAXIMUM      ####'
      write(*,*) '#### ALLOWED NUMBER.                             ####'
      write(*,*) '#####################################################'
      stop


998   write(*,*) '#####################################################'
      write(*,*) '#### FLEXPART MODEL SUBROUTINE READRELEASES:     ####'
      write(*,*) '####                                             ####'
      write(*,*) '#### FATAL ERROR -      "RELEASES" IS            ####'
      write(*,*) '#### CORRUPT. PLEASE CHECK YOUR INPUTS FOR       ####'
      write(*,*) '#### MISTAKES OR GET A NEW "RELEASES"-           ####'
      write(*,*) '#####################################################'
      stop


999   write(*,*) '#####################################################'
      write(*,*) '   FLEXPART MODEL SUBROUTINE READRELEASES: '
      write(*,*)
      write(*,*) 'FATAL ERROR - FILE CONTAINING PARTICLE RELEASE POINTS'
      write(*,*) 'POINTS IS NOT AVAILABLE OR YOU ARE NOT'
      write(*,*) 'PERMITTED FOR ANY ACCESS'
      write(*,*) '#####################################################'
      stop

800   write(*,*) ' #### TRAJECTORY MODEL ERROR! ERROR WHILE     #### ' 
      write(*,*) ' #### READING FILE PATHNAMES.                 #### ' 
      stop

801   write(*,*) '#### TRAJECTORY MODEL ERROR! FILE '// inputname 
      write(*,*) '#### CANNOT BE OPENED IN THE CURRENT WORKING #### '
      write(*,*) '#### DIRECTORY.                              #### '
      stop

803   write(*,*) ' #### FLEXPART MODEL ERROR! FILE   #### '
      write(*,'(a)') '     '//path(numpath+2*(k-1)+2) &
      (1:length(numpath+2*(k-1)+2))
      write(*,*) ' #### CANNOT BE OPENED             #### '
      stop

804   write(*,*) ' #### FLEXPART MODEL ERROR! FILE #### '
      write(*,'(a)') '     '//path(3)(1:length(3)) 
      write(*,*) ' #### CANNOT BE OPENED           #### '
      stop
      
996   write(*,*) '#####################################################'
      write(*,*) '#### FLEXPART MODEL ERROR!                      #### '
      write(*,*) '#### WET DEPOSITION SWITCHED ON, BUT NO HENRYS  #### '
      write(*,*) '#### CONSTANT IS SET                            ####'
      write(*,*) '#### PLEASE MODIFY SPECIES DESCR.               #### '
      write(*,*) '#####################################################'
      stop
      
end subroutine readinput
