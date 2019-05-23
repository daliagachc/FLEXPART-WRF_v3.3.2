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

! THIS CODE IS TO REDISTRIBUTE PARTICLES INVOLVED IN UPDRAFT OR/AND  DOWDRAFT
! TWO OPTIONS: 2- simply well-mixed inside updraft 
!                 mixing is based on convective mass flux, downdraft is not considered
!              3- particle positions based on entrainment/detrainment rates 
!                 if they are available
!   0- no convection 
!   1- has been used by old code

! SUBROUTINE NAME: redist_kf
!  CALLED by convmix_kf.f
!  meteoroloy data are from pre_redist_kf.f that is called inside convmix_kf.f
!
!   INPUT: 
!           dx   -- horizontal grid size
!           nuvzmax   -- # of layer for mass flux
!           umf -- updraft mass flux
!           uer -- updraft entrament flux
!           udr -- updraft detrainment flux
!           dmf -- downdraft mass flux
!           der -- downdraft entrainment flux
!           ddr -- downdraft detrainment flux
!           rho1d -- air density
!           dz1d-- delt z between full levels
!           ldirect -- flag for forward (+1) or backward (-1) run
!           delt    -- time step (s)
!           umfzf -- normalized up-mass-flux weighted distance,
!           dmfzf -- normalized down-mass-flux weighted distance,(not used)
!           fmix  -- estimated fraction of particles that will be mixed (for LCONVECTION=2)
!           zf    -- distance above ground at full levels
!           zh    --  """"""""""""""""""""""""half levels

!   IN/OUTPUT:
!           zp  -- particle z position  (IN/OUT)

       Subroutine redist_kf(mix_option,ldirect,delt,dx,            & ! IN
                          dz1d,nuvzmax, nuvz,umf,uer,udr,dmf,      & ! IN
                            der,ddr,rho1d,                         & ! IN
                            zf,zh,                                 & ! IN
                            umfzf,dmfzf,fmix,                      & ! IN
                            zp)                                      ! IN/OUT


        IMPLICIT NONE
        integer :: nuvzmax,nuvz,ldirect,i,j,k,kk,mix_option
        integer,parameter :: well_mix=2 
        integer,parameter :: prob_mix=3
        real,dimension(nuvzmax):: umf,uer,udr,zh,dz1d
        real,dimension(nuvzmax):: dmf,der,ddr,rho1d
        real,dimension(nuvzmax+1):: zf,umfzf,dmfzf
        real :: zp,delt,dx,fup,fdown,ftotal,totalmass
        real :: rn,rn1,rn2,ddz,fde,w_sub,ran3,fmix
        real :: uptop,downtop          ! top of updraft and downdraft
!        data well_mix,prob_mix/2,3/

! Check if particle position is below convection cloud top
! find top height of updraft
        uptop=0.0 
        do i=nuvz,1,-1
          if (umf(i) .gt. 1e-20) then
            if (i .eq. nuvz) uptop = zh(nuvz) 
            if (i .lt. nuvz) uptop = zh(i+1)
             goto 81
            endif
        enddo
81        continue                                ! updraft top

        downtop=0.0
        do i=nuvz,1,-1
          if (abs(dmf(i)) .gt. 1e-20) then
             if (i .eq. nuvz) downtop=zh(nuvz)
             if (i .lt. nuvz) downtop=zh(i+1)
              goto 82
          endif
        enddo
82       continue                                 ! downdraft top


         if (zp .ge. uptop) goto 89   ! no convective adjustment

           DO i=1,nuvz
             if (zp .ge. zf(i) .and. zp .lt. zf(i+1)) goto 90
           ENDDO
90          kk=i  !  kk is grid # of the particle before reposition
            totalmass=rho1d(kk)*dz1d(kk)*dx*dx   ! air mass kg


!+++++++++
! NOw start to reposition particle
!  
! Simply mix all particles by assigning a random position within cloud
! In this case, backward run is not wroking since this reposition process is not
!  reversible.

         IF (mix_option .eq. well_mix) then         ! Well-mixed

!!  only consider updraft (usually >> downdraft flux) 
!! Choose a random # evenly distributed in [0,1]

!hmjb - cannot be done this way!!! idum must change between calls otherwise it is not random
!hmjb - there are many calls like this in this function, all with the same problem. 
!hmjb - i did not dare trying to fix them...

            rn=ran3(88)
     
            if (rn .le. fmix) then       ! inside cloud 
              write(*,*)'well mixed, fmix=,rn=',fmix,rn
!    --|--  umfzf(j),  zf(j)
!    --|X-  rn2 (particle position)
!    --|--  umfzf(j-1),zf(j-1)

               rn2=ran3(881)

               do j=2,nuvz+1
                 if(umfzf(j) .ge. rn2) then 
                 zp = zf(j)-(umfzf(j)-rn)/ &
                          (umfzf(j)-umfzf(j-1))*dz1d(j-1)
                 goto 92
                 endif
               enddo
92             continue
             write(*,*)'old position=',kk, &
                       'repositioned at k=,umfzf(j),rn',j,umfzf(j),rn 
            endif                         ! inside cloud
         ENDIF                                      ! Well-mixed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! reposition based on probabilities of entrainment and detrainment
!!!!!!!!!!!!!!!! FORWARD treatment
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         IF (mix_option .eq. prob_mix .and. ldirect .gt. 0)then       !! prob_mix+ forward
           fup = abs(uer(kk))*delt/totalmass
           fdown=abs(der(kk))*delt/totalmass
           ftotal=fup+fdown                               ! ignore downward flux if 0*fdown
           if (ftotal .le. 1e-10) goto 89                    ! return if no cloud
! pick a random number to see if particle is inside cloud
            rn = ran3(88)
          

           if (rn .le. ftotal) then    ! in
            rn1=ran3(881)              
             write(*,*)'inside cloud, kk,ftotal,rn=',kk,ftotal,rn
              if (rn1 .le. fup/ftotal) then   ! updraft
! parcel will move upward till all particles inside are detrained
!   umf(k+1)=umf(k)+uer(k+1)-udr(k+1)
!
!    ---|--- 
!    ---|+++ 
!    ---|--- zf(j+2)
!    ---|+++ udr(j+1), half level, partile is detrained here
!    ---|--- zf(j+1)
!    ........
!    ---|+++ uer(kk), particle is entrained here
                   write(*,*)'updraft','fup=',fup,uer(kk)*delt,totalmass
                   write(*,*)'fdown=',fdown,fdown/ftotal
                do j=kk,nuvz-1 
                 rn=ran3(882)
                 if (umf(j) .eq. 0.0) goto 98
                 fde = abs(udr(j+1)/umf(j))   ! detrainment probability at level j+1
                 if (rn .lt. fde .or. zh(j+1) .ge. uptop) goto 98 ! being detrtained from air parcel
                                                           ! or reach cloud top 
                enddo
98              zp=zf(j+1) + rn*(zf(j+2)-zf(j+1)) ! particle is detrained between zf(j+1) and zf(j+2)              
                 write(*,*)'detrain at',j+1,rn,zp,'old=',zf(kk)
              else                            !up/down
! move downward till detrained or get to ground
                do j=kk,2,-1
                 rn=ran3(883)
                 if(dmf(j) .eq. 0.0) goto 102
                 fde = abs(ddr(j-1)/dmf(j))  ! detrainment at level j-1 
                 if(rn .lt. fde .or. j .eq. 2) goto 102
                enddo
102              zp=zf(j-1)+rn*(zf(j)-zf(j-1))

                 write(*,*)'downward prob,detrain at',fdown/ftotal,j-1

              endif                           ! downdraft

            else                        ! in /OUT
!! outside cloud
!! displace particle according to compensating subsidence 
             w_sub =(abs(umf(kk))-abs(dmf(kk)))/rho1d(kk)/dx/dx
              zp = zp - w_sub*delt
             write(*,*)'subsidence== distance(m)',w_sub*delt, &
                       'rn,ftotal=',rn,ftotal
            endif                       ! OUT
         ENDIF                                                    !! prob_mix + forward
!!!!!!!!!  END OF FORWARD treatment  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!! -------------------------------------------------------------
!!!!!!!!  backward 
!!!!!!!!    reverse all processes, e.g., up --> down, down --> up, entrainment -> detrainment
!!!!!!!!           detrainment --> entrainment,,,,,,,
!!!!!!!! -------------------------------------------------------------

         IF (mix_option .eq. prob_mix .and. ldirect .lt. 0) then   !! prob_mix + backward
           fdown = abs(udr(kk))*delt/totalmass
           fup   = abs(ddr(kk))*delt/totalmass
           ftotal=fup+fdown                          !
!           ftotal=fup                                ! ignore downdraft 
           if (ftotal .le. 1e-10) goto 89

! pick a random number to see if particle is inside cloud
            rn = ran3(88)
            if (rn .le. ftotal) then    ! in
            rn1=ran3(881)
              if (rn1 .lt. fup/ftotal) then   ! updraft
! parcel will move "upward" till all particles inside being detrained
 
                do j=kk,nuvz-1
                 rn=ran3(j)
                 if (dmf(j) .eq. 0.0) goto 981
                 fde = abs(der(j+1)/dmf(j))   ! "detrainment" probability at level j+1
                 if (rn .lt. fde .or. zh(j+1) .ge. downtop) goto 981 ! being detrtained from air parcel
                                                           ! or reach cloud top
                 enddo
981              zp=zf(j+1) + rn*(zf(j+2)-zf(j+1)) ! partice is detrained between zf(j+1) and j+2
 
              else                            !up/down
! move downward till detrained or get to ground
                do j=kk,2,-1
                 rn=ran3(j)
                 if (umf(j) .eq. 0.0) goto 1021
                 fde = abs(uer(j-1)/umf(j))  ! "detrainment" at level j-1
                 if(rn .lt. fde .or. j .eq. 2) goto 1021
                enddo
1021              zp=zf(j-1)+rn*(zf(j)-zf(j-1))
 
              endif                           ! downdraft
 
            else                        ! in /OUT
!! outside cloud
!! displace particle according to compensating subsidence
             w_sub =(abs(dmf(kk))-abs(umf(kk)))/rho1d(kk)/dx/dx
              zp = zp - w_sub*delt
            endif                       ! OUT
         ENDIF                                             !! prob_mix + backward

!-------------- END OF BACKWARD treatment ------------------





!+++++++++
! END OF calculating particle position 


89       continue

       return
       end subroutine redist_kf



