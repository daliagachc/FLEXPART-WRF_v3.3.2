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

! 8/30/2007 Created by Weiguo Wang 
! Notes:
! THIS CODE IS TO prepare meteorology data for REDISTRIBUTing PARTICLES 
!           INVOLVED IN UPDRAFT OR/AND  DOWDRAFT
!   Estimate fraction of particles that may be well-mixed (for well-mixing option)
 
!   prob of particle involving in clouds can be pre-computed here. will optimize code later

 
! SUBROUTINE NAME: pre_redist_kf
!   INPUT: 
!           nuvzmax-- max # of layers of flux
!           nuvz -- # of layer for work array
!           umf -- updraft mass flux (kg/s ?)
!           dmf -- downdraft mass flux
!           dz   -- different height between full levels (m) 
!           p1d -- press (pa) 
!           dx  -- horizontal grid size(m) 
!           dt  -- time step (s)
!C          cu_top1  --  cloud top index, zh(cu_top1)
!C          cu_bot1  --  cloud bottom index, ,zh(cu_bot1)
!   OUTPUT:
!            zf -- height above ground level at full levels (ddz)
!            zh -- height above ground level at half levels
!       only for option of simple mixing
!
!           umfzf -- normalized updraft mass flux*distance,min=0,max=1
!           dmfzf -- normalized downdraft mass flux*distance
!           fmix  -- fraction of paricels in cloud levels is mixed
!   CALLED by convmix_kf.f
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       Subroutine pre_redist_kf(nuvzmax,nuvz,umf,dmf,dz,p1d,dx,dt, &  ! IN &
            cu_bot1,cu_top1,    &                      ! IN
            zf,zh,              &                     ! OUT
            umfzf,dmfzf,fmix)                        ! OUT

        IMPLICIT NONE
        integer :: nuvzmax,nuvz,ldirection,n1,i,j,k
        real :: cu_bot1,cu_top1,dx,dt
        real,dimension(nuvzmax) :: umf,dmf,dz,umfdz,dmfdz,zh, &
                                    p1d
   
        real,dimension(nuvzmax+1) ::  dmfzf,umfzf,zf 
        real :: t1,t2,fmix,mass

! Compute flux height, here height is above the ground level

        zf(1)=0.0  ! ground level
        do i=2,nuvz+1
          zf(i) = zf(i-1) + dz(i)
          zh(i-1) = 0.5*( zf(i-1) + zf(i) )
        enddo 

!! cululative dz in cloud level, for reposition option 'well-mixed'
        do i=1,nuvz
          umfdz(i) = 0.0
          dmfdz(i) = 0.0
          if(umf(i).ne.0.0 .and. i.le.int(cu_top1) &
             .and.i.ge.int(cu_bot1)) umfdz(i)=dz(i)    ! mixed within cloud
!          if(umf(i).ne.0.0 .and. i.le.int(cu_top1)      ! mixed between ground and cloud top
!     &       ) umfdz(i)=dz(i)

          if(dmf(i).ne.0.0) dmfdz(i)=dz(i)
        enddo
! assume zero umf or dmf means no-cloud area
! Nomalize non-zero values (cloud up/downdraft)
          t1 = 0.0
          t2 = 0.0
         do i=1,nuvz
           t1=t1+umfdz(i)
           t2=t2+dmfdz(i) 
         enddo

         if (t1 .gt. 0.0) then 
          do i=1,nuvz
            umfdz(i)=umfdz(i)/t1
          enddo
         endif
         if (t2 .gt. 0.0) then
          DO i=1,nuvz
            dmfdz(i)=dmfdz(i)/t2
          ENDDO
         endif

! Weighted distance stating 0, ending 1.0
          umfzf(1)=0.0
          dmfzf(1)=0.0
          Do i=2,nuvz+1
            umfzf(i)=0.0
            dmfzf(i)=0.0
           if (i .le. int(cu_top1)) &
            umfzf(i)=umfzf(i-1)+abs(umfdz(i-1)) 
    
            dmfzf(i)=dmfzf(i-1)+abs(dmfdz(i-1))
            if (i .eq. 2)write(*,*)'int(cu_top1)=',int(cu_top1)
!             write(*,*)i,umfzf(i),umfdz(i)
          ENDdo 
! estimate fraction of particles in the convective column will be mixed by cloud
!   fmix=updraft flux at cloud base*dt/mass below cloud
!     fmix*dt=fraction
          mass=abs(p1d(1)-p1d( int(cu_bot1) ))
          if (mass .le. 5000.0) mass=5000.0 
          mass=dx*dx*mass/9.81
          fmix=abs(umf(int(cu_bot1)))*dt/mass
          write(*,*)'PRE_redist_kf.f, mass=,fmix=',mass,fmix
     
       return
 end subroutine pre_redist_kf


