!***********************************************************************
!* Copyright 2012,2013                                                 *
!* Jerome Brioude, Delia Arnold, Jerome Fast, 
!* Andreas Stohl, Petra Seibert, A. Frank, Gerhard Wotawa              *
!* Caroline Forster, Sabine Eckhardt, John Burkhart, Harald Sodemann 
!* M. Cassiani
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
      subroutine interpol_hmix(itime,xt,yt,zt,haux, &
        p1,p2,p3,p4,ddx,ddy,rddx,rddy,dtt,dt1,dt2, &
    ix,jy,ixp,jyp)

!                               i   i  i  i
!*******************************************************************************
!                                                                              *
!  This subroutine interpolates boundary layer top (h)                          *
!  dispersion.                                                                 *
!                                                                              *
!    Author:M. cassiani 2013                                                   *
!                                                                              *
!                                                                              *
!*******************************************************************************
  use par_mod
  use com_mod
!  use interpol_mod
!  use hanna_mod
  implicit none

      integer :: itime
      real :: xt,yt,zt,h

! Auxiliary variables needed for interpolation
     
      integer :: i,m,n,indexh,n2
      
      real,parameter ::eps=1.0e-30

 
  real :: h1(2),haux
  real :: p1,p2,p3,p4,ddx,ddy,rddx,rddy,dtt,dt1,dt2
  integer :: ix,jy,ixp,jyp,ngrid,indz,indzp 
  logical :: indzindicator(nzmax)
 


!********************************************
! Multilinear interpolation in time and space
!********************************************

! Determine the lower left corner and its distance to the current position
!*************************************************************************

      ddx=xt-real(ix)                     
      ddy=yt-real(jy)
      rddx=1.-ddx
      rddy=1.-ddy
      p1=rddx*rddy
      p2=ddx*rddy
      p3=rddx*ddy
      p4=ddx*ddy

! Calculate variables for time interpolation
!*******************************************

      dt1=real(itime-memtime(1))
      dt2=real(memtime(2)-itime)
      dtt=1./(dt1+dt2)



!*****************************************

! a) Bilinear horizontal interpolation

  do m=1,2
    indexh=memind(m)

    h1(m)=p1*hmix(ix ,jy ,1,indexh) &
         + p2*hmix(ixp,jy ,1,indexh) &
         + p3*hmix(ix ,jyp,1,indexh) &
         + p4*hmix(ixp,jyp,1,indexh)
   
  end do
      

! b) Temporal interpolation

      haux=(h1(1)*dt2+h1(2)*dt1)*dtt     


end subroutine interpol_hmix

