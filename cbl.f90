!***********************************************************************
!* Copyright 2012,2013                                                *
!* Jerome Brioude, Delia Arnold, Andreas Stohl, Wayne Angevine,       *
!* John Burkhart, Massimo Cassiani, Adam Dingwell, Richard C Easter, Sabine Eckhardt,*
!* Stephanie Evan, Jerome D Fast, Don Morton, Ignacio Pisso,          *
!* Petra Seibert, Gerard Wotawa, Caroline Forster, Harald Sodemann,   *
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
    
    subroutine cbl(wp,zp,ust,wst,h,rhoa,rhograd,sigmaw,dsigmawdz,tlw,ptot,Q,phi,ath,bth,ol,flagrein)
!                  i  i i  i   i  i    i     i       i         i   o   o   o   o  o   o   i     i/o
!=============== well mixed formulation of CBL skewed vertical profiles following  LHH 1996 with profile of w3 from lHB 2000                    ========
!=============== LHH formulation has been modified to account for variable density profiles and backward in time or forward in time simulations ========
!=============== by  Massimo Cassiani ( mc ), NILU,  2012-2013, reference to Cassiani et al. 2013 (to be submitted...)                          ========
!=======================================================================================================================================================
!======================================================================================================
! wp: particle velocity
! zp: particle position
! ust: velcotiy scale
! wst: convective velcotiy scale
! h: boundary layer top
! rhoa: air density
! rhograd: air densiy vertical gradient
! sigmaw: turbulent flutuation of vertical velocity standard deviation
! dsigmawdz: derivative of above
! tlw: local lagrangina time scale
! ptot: pdf value for the particle velocity in drift coefficient, see Cassiani et al. 2013, not used
! Q: part of drift coefficient, not used
! phi: part of drift coeffcient, not used
! ath: drift coefficient, used
! bth: diffusion coeffcient, sued
! ol: Obukhov lenght
! flagrein: set accordingly to conditon below if 1 then re-initialize particle velocity
!======================================================================================================
    use par_mod, only:pi
    use com_mod, only:ldirect
!   use ieee_arithmetic    
    implicit none


    real :: usurad2,usurad2p,C0,costluar4,eps 
    parameter  (usurad2=0.7071067812,usurad2p=0.3989422804,C0=2,costluar4=0.66667,eps=0.000001)

    integer flagrein  !re-initlization flag for the particle velocity 
    
    real :: wp,zp,ust,wst,h,dens,ddens,sigmaw,dsigmawdz,tlw,rhoa,rhograd
    real ::fluarw,fluarw2
    real ::w3,w2
    real ::dw3,dw2
    real ::wb,wa
    real ::deltawa,deltawb
    real ::wold,wold2,wold_z
    real ::pa,pb,alfa
    real ::Phi,Q,ptot  
    real :: timedir
    real ::cuberoot
    real ::z0,ol,transition !added ol & transition with respect to cbl.f90 without transition 
    

    real :: &
    erf, &
    aperfa, &
    aperfb, &
    ath, &
    bth

    real ::  &
    pow, &
    z, &    
    skew, &
    skew2, &
    radw2, &
    rluarw, &
    xluarw, &
    aluarw, &
    bluarw, &
    sigmawa, &
    sigmawb, &
    dskew, &
    dradw2, &
    dfluarw, &
    drluarw, &
    dxluarw, &
    daluarw, &
    dbluarw, &
    dsigmawa, &
    dsigmawb, &
    dwa, &
    dwb, &
    sigmawa2, &
    sigmawb2
    
    
    dens=rhoa
    ddens=rhograd
             
    
    timedir=ldirect !direction of time forward (1) or backward(-1)
    !========================= assign z ==============================
    z=(zp/h)
    
    transition=1.
    !if (ol.lt.-50) transition=((sin(((ol+100.)/100.)*pi))-1.)/2.
    if (-h/ol.lt.15) transition=((sin((((-h/ol)+10.)/10.)*pi)))/2.+0.5  !transition fucntion to smoohtly 
    !========================= secondo moment of vertical velocity  =====================
    !!!	    w2=1.4*(z**1.5*(1.-z))**(2./3.)
    !w2=(1.7*(z*(1.-0.7*z)*(1.-z))**(2./3.))*(wst**2)
    w2=(sigmaw*sigmaw)
    dw2=(2.*sigmaw*dsigmawdz)
    !dw2=(1.7*(2./3.)*(z*(1.-0.7*z)*(1.-z))**(-1./3.)* &
    !(((1.-0.7*z)*(1.-z))+z*(-0.7)*(1.-z)+z*(1.-0.7*z)*(-1.))) *(wst**2)*1/h

    !=================== dissipation fo turbulent tke  =========================
    !alfa=0.4 !(0.75-(0.5*z*z))**(3./2.) DISSIPAZIONE ADIMENSIONALE
    alfa=2.*w2/(C0*tlw)

    !========================================================================
    wold=timedir*wp !time direction enter here for backward calculualtions
    !wold_z=wp
    ! =======================================================================
    !------------------------------ momento terzo ============================
    !!  w3=0.8*(w2**(3./2.))
    !!	dw3=0.8*1.5*w2**0.5*dw2

    ! dw3=((1.2*z*((1.-z)**(3./2.)))+eps)*(wst**3)
    ! dw3=(1.2*(((1.-z)**(3./2.))+z*1.5*((1.-z)**(1./2.))*(-1.)))*(wst**3)
    ! 3=(1.2*z*((1.-z)**(3./2.)))
    ! w3=(1.2*(((1.-z)**(3./2.))+z*1.5*((1.-z)**(1./2.))*(-1.)))
    
                            !w3=((1.2*z*((1.-z)**(3./2.)))*1.5+eps)*(wst**3)  !1.5 to increase skeweness see also initalize_cbl_vel.f90 
                            !dw3=(1.2*(((1.-z)**(3./2.))+z*1.5*((1.-z)**(1./2.))*(-1.)))*(wst**3)*(1./h)*1.5
                            
                            w3=((1.2*z*((1.-z)**(3./2.)))+eps)*(wst**3)*transition
                            dw3=(1.2*(((1.-z)**(3./2.))+z*1.5*((1.-z)**(1./2.))*(-1.)))*(wst**3)*(1./h)*transition
    !============================================================================0
   
    
                            skew=w3/(w2**1.5)
                            skew2=skew*skew
                            dskew=(dw3*w2**(1.5)-w3*1.5*w2**0.5*dw2)/w2**3
                            radw2=w2**0.5
                            dradw2=0.5*w2**(-0.5)*dw2
                            !costluar4=0.66667  ! costante da LHH
                            fluarw=costluar4*(cuberoot(skew))	!
                            fluarw2=fluarw*fluarw
    
    
    if (skew.ne.0) then
      
                             dfluarw=costluar4*(1./3.)*cuberoot(skew**(-2.))*dskew
        
                            rluarw=(1.+fluarw2)**3.*skew2/((3.+fluarw2)**2.*fluarw2)  !-> r
                            xluarw=(1.+fluarw2)**1.5*skew/((3.+fluarw2)*fluarw)       !
        
                            drluarw=( ((3.*(1.+fluarw2)**2*(2.*fluarw*dfluarw)*skew2)+ &
                            (1.+fluarw2)**3*2.*skew*dskew) *(3.+fluarw2)**2.*fluarw2 - &
                            (1.+fluarw2)**3*skew2* &
                            ( (2.*(3.+fluarw2)*(2.*fluarw*dfluarw)*fluarw2) + &
                            (3.+fluarw2)**2*2.*fluarw*dfluarw) )/ &
                            (((3.+fluarw2)**2.*fluarw2)**2)
                                      
                            dxluarw=( ((1.5*(1.+fluarw2)**0.5*(2.*fluarw*dfluarw)*skew)+ &
                            (1.+fluarw2)**1.5*dskew) *(3.+fluarw2)*fluarw - &
                            (1.+fluarw2)**1.5*skew* &
                            (3.*dfluarw+3*fluarw2*dfluarw) )/ &
                            (((3.+fluarw2)*fluarw)**2)

        
    else        
        dfluarw=0.
        rluarw=0.
        drluarw=0.
        xluarw=0.
        dxluarw=0.
    end if



                            aluarw=0.5*(1.-xluarw/(4.+rluarw)**0.5)
                            bluarw=1.-aluarw

                           daluarw=-0.5*(  (dxluarw*(4.+rluarw)**0.5) - &
                            (0.5*xluarw*(4.+rluarw)**(-0.5)*drluarw) ) &
                            /(4.+rluarw) 
                            dbluarw=-daluarw
    

                            sigmawa=radw2*(bluarw/(aluarw*(1.+fluarw2)))**0.5
                            sigmawb=radw2*(aluarw/(bluarw*(1.+fluarw2)))**0.5


                            dsigmawa=dradw2*(bluarw/(aluarw*(1.+fluarw2)))**0.5+ &
                            radw2*( &
                            (0.5*(bluarw/(aluarw*(1.+fluarw2)))**(-0.5))      * &
                            ( &
                            (dbluarw*(aluarw*(1.+fluarw2))- &
                            bluarw*(daluarw*(1.+fluarw2)+aluarw*2.*fluarw*dfluarw)) &
                            /((aluarw*(1.+fluarw2))**2) &
                            ) &
                            )

                            dsigmawb=dradw2*(aluarw/(bluarw*(1.+fluarw2)))**0.5+ &
                            radw2*( &
                            (0.5*(aluarw/(bluarw*(1.+fluarw2)))**(-0.5))      * &
                            ( &
                            (daluarw*(bluarw*(1.+fluarw2))- &
                            aluarw*(dbluarw*(1.+fluarw2)+bluarw*2.*fluarw*dfluarw)) &
                            /((bluarw*(1.+fluarw2))**2) &
                            ) &
                            )
    
                            wa=(fluarw*sigmawa)
                            wb=(fluarw*sigmawb)
                            dwa=dfluarw*sigmawa+fluarw*dsigmawa
                            dwb=dfluarw*sigmawb+fluarw*dsigmawb

                            deltawa=wold-wa
                            deltawb=wold+wb
                            wold2=wold*wold
                            sigmawa2=sigmawa*sigmawa
                            sigmawb2=sigmawb*sigmawb

                            pa=(usurad2p*(1./sigmawa))*(exp(-(0.5*((deltawa/sigmawa)**2.))))
                            pb=(usurad2p*(1./sigmawb))*(exp(-(0.5*((deltawb/sigmawb)**2.))))
                            
                            if (abs(deltawa).gt.10.*sigmawa.and.abs(deltawb).gt.10.*sigmawb) flagrein=1  !added control flag for re-initialization of velocity
!                           if (abs(deltawa).gt.6.*sigmawa.and.abs(deltawb).gt.6.*sigmawb) flagrein=1  !added control flag for re-initialization of velocity
                            	    
                            ptot=dens*aluarw*pa+dens*bluarw*pb
                               
                            aperfa=deltawa*usurad2/sigmawa
                            aperfb=deltawb*usurad2/sigmawb

!       if ((ieee_is_nan(aperfa).or.ieee_is_nan(aperfb)).and.flagrein.eq.0) &
!          print*,'PROBLEM',deltawa,deltawb,sigmawa,sigmawb,wp,zp,ust,wst,h,rhoa,rhograd,sigmaw,dsigmawdz,tlw,ptot,Q,phi,ath,bth,ol,flagrein
                            Phi=-0.5* &
                            (aluarw*dens*dwa+dens*wa*daluarw+aluarw*wa*ddens)*erf(aperfa) &
                            +sigmawa*(aluarw*dens*dsigmawa*(wold2/sigmawa2+1.)+ &
                            sigmawa*dens*daluarw+sigmawa*ddens*aluarw+ &
                            aluarw*wold*dens/sigmawa2*(sigmawa*dwa-wa*dsigmawa))*pa &
                            +0.5* &
                            (bluarw*dens*dwb+wb*dens*dbluarw+wb*bluarw*ddens)*erf(aperfb) &
                            +sigmawb*(bluarw*dens*dsigmawb*(wold2/sigmawb2+1.)+ &
                            sigmawb*dens*dbluarw+sigmawb*ddens*bluarw+ &
                            bluarw*wold*dens/sigmawb2*(-sigmawb*dwb+wb*dsigmawb))*pb


                             Q=timedir*((aluarw*dens*deltawa/sigmawa2)*pa+ &
                            (bluarw*dens*deltawb/sigmawb2)*pb)

    


                            ath=(1./ptot)*(-(C0/2.)*alfa*Q+phi)  !drift coefficient
                            bth=sqrt(C0*alfa)                    !diffusion coefficient
                           

    

    return


    end
                            




    FUNCTION CUBEROOT (X) RESULT (Y)

    IMPLICIT NONE

    real, INTENT(IN) :: X
    real:: Y

    real, PARAMETER :: THIRD = 0.333333333


    Y = SIGN((ABS(X))**THIRD, X)

    RETURN

    END FUNCTION CUBEROOT
    
    
    

    FUNCTION CUBEROOTD (X) RESULT (Y)

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: X
    DOUBLE PRECISION :: Y

    DOUBLE PRECISION, PARAMETER :: THIRD = 0.33333333333333333333333333333333333333333333333333333333333333333333333333333333333D0


    Y = SIGN((ABS(X))**THIRD, X)

    RETURN

    END FUNCTION CUBEROOTD
