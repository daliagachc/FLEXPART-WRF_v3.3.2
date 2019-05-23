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

      subroutine tke_partition_my(z, &
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
! 11 Sep 2007 implement different formulations for growing and decaying turbulences
 
!      include 'includepar'
!      include 'includecom'
!      include 'includehanna'
!      include 'includeinterpol'
  use par_mod
  use com_mod
! use hanna_mod
!  use interpol_mod

      implicit none
      real :: z,zzz,zzz1,dz,dz1,dz2,tke_z,dpttdz
      real :: fu2,fv2,fw2,ftotal,ygu,ygv,ygm,ygh,yl,ylmax,ysm,ysh
      real :: siguprof(2),sigvprof(2),sigwprof(2)
      real :: ya1
      real :: ya2
      real :: yb1
      real :: yb2
      real :: yc1
      real :: yse
      real :: yap
      real :: e1,e2,e3,e4,e5,rf,rfc,rf1,rf2,ch,cm
      real :: er,smr,shr,dudz,dvdz,tke_mid
      integer :: k,indz,indzp
     real :: uprof(nzmax),vprof(nzmax),tkeprof(nzmax),pttprof(nzmax)
    real :: ust,wst,ol,h,zeta,sigu,sigv,tlu,tlv,tlw
    real :: sigw,dsigwdz,dsigw2dz

!      data ya1,ya2,yb1,yb2,yc1,yse,yap/
!     +   0.92,0.74,16.6,10.1,0.08,0.20,0.17/
       ya1=0.92
       ya2=0.74
       yb1=16.6
       yb2=10.1
       yc1=0.08
       yse=0.20
       yap=0.17 

!- interpolate tke
          dz=1./(height(indzp)-height(indz))
          dz1=(z - height(indz))*dz
          dz2=(height(indzp)-z)*dz
 
          tke_z=max(0.01,dz1*tkeprof(indzp)+dz2*tkeprof(indz))

!- compute turbulent length scale
!  yl -- turbulence length scale
!  ylmax -- max yl
          yl=0.01
          ylmax=0.01
          do k=1,nz-1
            zzz=0.5*(height(k)+height(k+1))
            zzz1=height(k+1)-height(k)
            yl=yl+zzz*sqrt(0.5*(tkeprof(k)+tkeprof(k+1)))*zzz1
            ylmax=ylmax+zzz1*sqrt(0.5*(tkeprof(k)+tkeprof(k+1)))
          enddo
            ylmax=0.1*yl/ylmax
            yl=0.4*(z+0.1)/(1.0+0.4*(z+0.1)/ylmax)
            if (pttprof(indzp) .gt. pttprof(indz)) then

              dpttdz = (pttprof(indzp)-pttprof(indz))*dz

              yl=min(yl,0.75*sqrt(2.0*tke_z/(9.8/300.0)/dpttdz))
!            endif
!- compute nondimensional vertical gradients
          dudz=(uprof(indzp)-uprof(indz))*dz
          dvdz=(vprof(indzp)-vprof(indz))*dz
         ygu=yl/sqrt(2.0*tke_z)*dudz
         ygv=yl/sqrt(2.0*tke_z)*dvdz
         ygm=ygu*ygu+ygv*ygv+0.001       ! in case of zero
         ygh=-yl*yl/2.0/tke_z*dpttdz*9.8/300.0

!--  compute SM,SH
!-   sm,sh nondimensional eddy diffusivities

         ysm=1.0-3.0*ya2*(7.0*ya1+yb2)*ygh+ &
              27.0*ya1*ya2*ya2*(4*ya1+yb2)*ygh*ygh+ &
              6.0*ya1*ya1*(1.0-3.0*ya2*(yb2-3.0*ya2)*ygh)*ygm
         ysm=ya1*(1-3*yc1-3*ya2*(yb2*(1-3*yc1)-12*ya1*yc1-3*ya2) &
                 *ygh)/ysm

         ysh=ya2*(1-6*ya1*ysm*ygm)/(1-3*ya2*ygh*(4*ya1+yb2))   


! -- compute the equlibrium TKE, er 
        e1=yb1-6.0*ya1
        e2=12*ya1+yb1+3.0*yb2
        e3=yb1*(1.0-3.0*yc1)-6.0*ya1
        e4=yb1*(1.0-3.0*yc1)+12*ya1+9.0*ya2
        e5=yb1+3.0*ya1+3.0*yb2
        ch=ya2*e2/yb1
        cm=ya1*e4/ya2/e5
        rf1=e3/e4
        rf2=e1/e5
        rfc=e1/e2
        rf=-(ysh/ysm)*( ygh/ygm )
!          write(*,*)'ysh,ysm',ysh,ysm
!          write(*,*)'ygh,ygm,dpttdz, ptt',ygh,ygm,dpttdz,pttprof(indzp)
!          if (rf .le. 0) write(*,*)'aaaaaaaaaaaaaaa'        
        shr = 0.0
        if (rf .ge. rfc) then 
            shr=ch*1.0e-4            ! Ri > critical value, then =
        elseif ((1.0-rf) .ne. 0.0) then
            shr=ch*(rfc-rf)/(1.0-rf)
        endif

        smr = 0.0
        if (rf .ge. rf1) then 
           smr=cm*shr*1.0e-4
        elseif((rf2-rf) .ne. 0.0 .and. rf .lt. rf2) then 
           smr=cm*(rf1-rf)/(rf2-rf)*shr
        endif
        
        er=0.5*yb1*yl*yl*(smr*dudz*dudz*dvdz*dvdz- &
                          shr*dpttdz*9.8/300.0)
        er=max(er,0.01)
        if (tke_z .ge. er ) then                    ! decaying turbuelnce  
          tke_mid=tke_z
        else
          tke_mid=er
          ysm=sqrt(tke_z/er)*smr
          ysh=sqrt(tke_z/er)*shr
        endif 
 
!- fractions of u2,v2,w2
!        print*,'in tke',tke_mid,ya1,ysm,ygm,ysh,ygh,ygu
         fu2=0.333-2*ya1*(ysm*ygm+ysh*ygh)+6*ya1*ysm*ygu*ygu
         fv2=0.333-2*ya1*(ysm*ygm+ysh*ygh)+6*ya1*ysm*ygv*ygv
         fw2=0.333-2*ya1*(ysm*ygm+ysh*ygh)+6*ya1*ysh*ygh
!        fw2=0.333-2*ya1*(ysm*ygm+ysh*ygh)+6*ya1*ysm*ygh*ygh
!        fw2=fw2*0.5
         if(fu2*fv2*fw2 .le. 0.0) then 
           fu2=0.333
           fv2=0.333
           fw2=0.333
         endif

        else ! test on pttprof

           fu2=0.333
           fv2=0.333
           fw2=0.333
          tke_mid=tke_z
        endif
          !ftotal=(fu2+fv2+fw2)
          !fu2=fu2/ftotal
          !fv2=fv2/ftotal
          !fw2=fw2/ftotal

! JB: make the partitionning equal
!          fu2=0.333
!          fv2=0.333
!          fw2=0.333
!
         sigu=sqrt(2.0*tke_mid*fu2)
         sigv=sqrt(2.0*tke_mid*fv2)
         sigw=sqrt(2.0*tke_mid*fw2)

         tlu=2.0*yl/sigu
         tlv=2.0*yl/sigv
         tlw=2.0*yl/sigw

         tlu=max(10.,tlu)
         tlv=max(10.,tlv)
         tlw=max(30.,tlw)

!- dsigw/dz,dsigw2/dz ; assuming fraction is not changed

         dsigwdz=(sqrt(2.0*tkeprof(indzp)*fw2)- &
                  sqrt(2.0*tkeprof(indz)*fw2))*dz
         dsigwdz=max(1.0e-10,dsigwdz)

         dsigw2dz=(2.0*tkeprof(indzp)*fw2-2.0*tkeprof(indz)*fw2)*dz
         dsigw2dz=max(1.0e-10,dsigw2dz)

!         write(*,*)'tke=',tke_z, 'er=',er
!         write(*,*)'rfc,rf1,rf2,rf=',rfc,rf1,rf2,rf
!          write(*,*)'smr=',smr,'ysm=',ysm
!         write(*,*)'shr=',shr,'ysh=',ysh
!         write(*,*)'fu2,fv2,fw2=',fu2,fv2,fw2
!         write(*,*)'ftotal=',ftotal
!         if (fu2*fv2*fw2 .le. 0) stop
end subroutine tke_partition_my
