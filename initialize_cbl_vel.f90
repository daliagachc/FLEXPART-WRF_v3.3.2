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
    
    
    subroutine initialize_cbl_vel(idum,zp,ust,wst,h,sigmaw,wp,dcas,dcas1,ol)
!                                  i/o   i  i   i  i     i  o   i    i       i
! idum: for random number but not usednot used
! zp: particle psition
! ust: velocity scale, not used
! wst: ocnvective velcotiy scale
! sigmaW: standard deviaiton of vertical velocity
! wp: particle velocity
! dcas: for random number
! dcas1: for random number
! ol: Obukhov lenght
!=============== initilization of particle velcoity based on CBL skewed vertical profiles and formulation of LHH 1996 with profile of w3 from lHB 2000         ==================
!=============== by  Massimo Cassiani ( mc ), NILU,  2012-2013, reference to Cassiani et al. 2013..                                                            ==================
!================================================================================================================================================================================   
    use par_mod, only:pi
    use com_mod, only:ldirect
!    use ieee_arithmetic
    
    implicit none


    real :: usurad2,usurad2p,C0,costluar4,eps 
    parameter  (usurad2=0.7071067812,usurad2p=0.3989422804,C0=2,costluar4=0.66667,eps=0.000001)

    integer idum
    real :: wp,zp,ust,wst,h,dens,ddens,sigmaw,dsigmawdz,tlw,dcas,dcas1,ran3,gasdev
    real :: w3,w2
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
    real ol, transition
    
    !---------------------------------------------------------------------------  
         !timedir=dble(ldirect) !direction of time forward (1) or backward(-1)
     timedir=ldirect !time direction forward (1) or backward(-1)
	 z=zp/h  !hn is the boundarylayer top
	 
     transition=1.
     if (-h/ol.lt.15) transition=((sin((((-h/ol)+10.)/10.)*pi)))/2.+0.5 !transtion from near neutral to unstable
     
     
     !w2=(1.7*(z*(1.-0.7*z)*(1.-z))**(2./3.))*(wst**2)
	 w2=sigmaw*sigmaw
	 !w3=(((1.2*z*((1.-z)**(3./2.)))+eps)*wst**3) *1.5	!the 1.5 is to test with increased skeweness see also cbl.f90
     w3=(((1.2*z*((1.-z)**(3./2.)))+eps)*wst**3)*transition 
	 skew=w3/(w2**1.5)
	 skew2=skew*skew
     radw2=sqrt(w2) !sigmaw
     
     if (skew.ne.0) then      
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

	 !dcas=ran3(idum) !pass from outside

	 if (dcas.le.aluarw) then
	  !dcas1=gasdev(idum) !pass from outside
	  wp=timedir*(dcas1*sigmawa+wa)
	 else
	  !dcas1=gasdev(idum) !pass from outside
	  wp=timedir*(dcas1*sigmawb-wb)
	 end if   
!    if (ieee_is_nan(wp)) print*,'PROBLEM INIT',wp,timedir, &
!      dcas1,sigmawa,wa,sigmawb,wb,idum,zp,ust,wst,h,sigmaw,wp,dcas,dcas1,ol    
    return
    end
         
