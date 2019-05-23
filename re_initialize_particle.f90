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
                                    
    
    subroutine re_initialize_particle(zp,ust,wst,h,sigmaw,wp,nrand,ol)
!                                      i   i  i   i  i    i/o  i/o 
!zp: particle position
!ust: velocity scale
!wst: velocity scale
!sigmaw: vertical velcotiy standard deviation
!wp: particle velocity
!nrand: random number counter
!=============== CBL skewed vertical profiles and formulation of LHH 1996 with profile of w3 from lHB 2000   ================================================================
!=============== LHH formulation has been modified to account for variable density profiles and backward in time or forward in time simulations                             =                    
!=============== by  Massimo Cassiani ( mc ) , NILU,  2012-2013                                                                                                             =
!=============== this routine re-initialize particle velocity if a numerical instability in the cbl scheme met a specific condition                                         =
!=============== (see routine cbl.f90 and Cassiani et al. 2013                                                                                                              =
!=============== the particle velocity is extracted from the updraft and downdraft distribution as required                                                                 =
!=============== this re-initialization si not perfectly consistent with teh well-mixed condition see Cassiani et al. 2013 for details but the error introduced is small    =
!=============== but for the rpesent this is faste and simpler and shoudl be ok                                                                                             =         
!============================================================================================================================================================================   
    use par_mod, only:pi
    use com_mod, only:ldirect,rannumb
!    use ieee_arithmetic
 
    implicit none


    real :: usurad2,usurad2p,C0,costluar4,eps 
    parameter  (usurad2=0.7071067812,usurad2p=0.3989422804,C0=2,costluar4=0.66667,eps=0.000001)

    integer idum,nrand
    real :: wp,zp,ust,wst,h,dens,ddens,sigmaw,dsigmawdz,tlw,dcas,dcas1,ran3,gasdev
    real :: w3,w2,wpold
    real ::  z, &    
    skew, &
    skew2, &
    radw2, &
    fluarw,fluarw2, &
    rluarw, &
    xluarw, &
    aluarw, &
    bluarw, &
    sigmawa, &
    sigmawb, &  
    ath, &
    bth, &
    wb,wa 
    real timedir
    real ol,transition
    !---------------------------------------------------------------------------
!     print*,'INSIDE INIT',zp,ust,wst,h,sigmaw,wp,nrand
!   wpold=wp 
          nrand=nrand+1
     dcas1=rannumb(nrand)
     timedir=ldirect !direction of time forward (1) or backward(-1)
    z=zp/h
     
     transition=1.  !comment by mc: in this version added transtion fucntion see Cassiani et al. 2013 
     if (-h/ol.lt.15) transition=((sin((((-h/ol)+10.)/10.)*pi)))/2.+0.5
    
     !   w2=((1.7*(z*(1.-0.7*z)*(1.-z))**(2./3.))+1.e-2)*(wst**2)
    w2=sigmaw*sigmaw !this is correct and use hanna routine if commented it is for test reason 
    !w3=(((1.2*z*((1.-z)**(3./2.)))+eps)*wst**3) *1.5	!the 1.5 is to test with increased skeweness see also cbl.f90
     w3=(((1.2*z*((1.-z)**(3./2.)))+eps)*wst**3) *transition !note added a transition fucntion, comemnt by mc
    skew=w3/(w2**1.5)
    skew2=skew*skew
     radw2=sigmaw !sqrt(w2) 
    
     if (skew.ne.0) then   ! the  limit must be considered explicitly to avoid NaN
     fluarw=costluar4*skew**0.333333333333333
    fluarw2=fluarw*fluarw
    rluarw=(1.+fluarw2)**3.*skew2/((3.+fluarw2)**2.*fluarw2)  !-> r
     xluarw=rluarw**0.5 !(1.+fluarw2)**1.5*skew/((3.+fluarw2)*fluarw)    !----> r^1/2
     else        
     fluarw=0.
     fluarw2=0.
     rluarw=0.        
     xluarw=0.        
     end if     
     
     aluarw=0.5*(1.-xluarw/(4.+rluarw)**0.5)
    bluarw=1.-aluarw
        
     sigmawa=radw2*(bluarw/(aluarw*(1.+fluarw2)))**0.5
     sigmawb=radw2*(aluarw/(bluarw*(1.+fluarw2)))**0.5

    wa=(fluarw*sigmawa)
    wb=(fluarw*sigmawb)

    

    if ((sign(1.,wp)*timedir).gt.0) then !updraft
     
100     wp=(dcas1*sigmawa+wa)
      if (wp.lt.0)  then
          nrand=nrand+1
          dcas1=rannumb(nrand)
          goto 100
      end if
     wp=wp*timedir
    else if ((sign(1.,wp)*timedir).lt.0) then !downdraft
101     wp=(dcas1*sigmawb-wb)
      if (wp.gt.0)  then 
           nrand=nrand+1
          dcas1=rannumb(nrand)
          goto 101
      end if
      wp=wp*timedir
    end if   
!if (ieee_is_nan(wp)) print*,'PROBLEM INSIDE',dcas1,nrand,sigmawa,fluarw,sigmawb,wb,aluarw,fluarw,xluarw,zp,ust,wst,h,sigmaw,wpold,nrand
         return
         end
         
