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
subroutine drydepokernel(nunc,deposit,x,y,itage,nage,kp)
  !                          i      i    i i  i
  !*****************************************************************************
  !                                                                            *
  !     Attribution of the deposition to the grid using a uniform kernel with  *
  !     bandwidths dx and dy.                                                  *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     26 December 1996                                                       *
  !
  !     D. Arnold: modification to skip the kernel the first 3 hours.
  !     then, modification to a regular lat-lon.
  !
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! nunc             uncertainty class of the respective particle              *
  ! nage             age class of the respective particle                      *
  ! deposit          amount (kg) to be deposited                               *
  !                                                                            *
  !*****************************************************************************

  use unc_mod
  use par_mod
  use com_mod

  implicit none

  real :: x,y,deposit(maxspec),ddx,ddy,xl,yl,wx,wy,w
  integer :: ix,jy,ixp,jyp,ks,nunc,nage,kp
! CDA new declarations
   real :: rhoprof(2),rhoi,xlon,ylat,xl2,yl2
   integer :: itage
! CDA

!JB
  if (outgrid_option.eq.0) then
! CDA
  xl=(x*dx+xoutshift)/dxout
  yl=(y*dy+youtshift)/dyout
  elseif (outgrid_option.eq.1) then
! CDA new code:
  xl2=x*dx+xmet0
  yl2=y*dy+ymet0
  call xymeter_to_ll_wrf(xl2,yl2,xlon,ylat)
  xl=(xlon-outlon0)/dxoutl
  yl=(ylat-outlat0)/dyoutl
  endif


  ix=int(xl)
  jy=int(yl)

! CDA skip kernel
      if (itage.lt.7200) then ! no kernel, direct attribution to grid cell
        do ks=1,nspec
          if ((abs(deposit(ks)).gt.0).and.DRYDEPSPEC(ks)) then
!$OMP CRITICAL 
            if ((ix.ge.0).and.(jy.ge.0).and.(ix.le.numxgrid-1).and. &
                                      (jy.le.numygrid-1)) &
        drygridunc2(ix,jy,ks,kp,nunc,nage)= &
          drygridunc2(ix,jy,ks,kp,nunc,nage)+deposit(ks)
!$OMP END CRITICAL
          endif
        enddo       
 else ! attribution via uniform kernel

  ddx=xl-real(ix)                   ! distance to left cell border
  ddy=yl-real(jy)                   ! distance to lower cell border

  if (ddx.gt.0.5) then
    ixp=ix+1
    wx=1.5-ddx
  else
    ixp=ix-1
    wx=0.5+ddx
  endif

  if (ddy.gt.0.5) then
    jyp=jy+1
    wy=1.5-ddy
  else
    jyp=jy-1
    wy=0.5+ddy
  endif


  ! Determine mass fractions for four grid points
  !**********************************************
    do ks=1,nspec

    if ((abs(deposit(ks)).gt.0).and.DRYDEPSPEC(ks)) then

   if ((ix.ge.0).and.(jy.ge.0).and.(ix.le.numxgrid-1).and.(jy.le.numygrid-1)) then
    w=wx*wy
!$OMP CRITICAL 
       drygridunc2(ix,jy,ks,kp,nunc,nage)= &
            drygridunc2(ix,jy,ks,kp,nunc,nage)+deposit(ks)*w
!$OMP END CRITICAL
      continue
  endif

  if ((ixp.ge.0).and.(jyp.ge.0).and.(ixp.le.numxgrid-1).and. &
       (jyp.le.numygrid-1)) then
    w=(1.-wx)*(1.-wy)
!$OMP CRITICAL 
      drygridunc2(ixp,jyp,ks,kp,nunc,nage)= &
           drygridunc2(ixp,jyp,ks,kp,nunc,nage)+deposit(ks)*w
!$OMP END CRITICAL
  endif

  if ((ixp.ge.0).and.(jy.ge.0).and.(ixp.le.numxgrid-1).and. &
       (jy.le.numygrid-1)) then
   w=(1.-wx)*wy
!$OMP CRITICAL 
      drygridunc2(ixp,jy,ks,kp,nunc,nage)= &
           drygridunc2(ixp,jy,ks,kp,nunc,nage)+deposit(ks)*w
!$OMP END CRITICAL
  endif

  if ((ix.ge.0).and.(jyp.ge.0).and.(ix.le.numxgrid-1).and. &
       (jyp.le.numygrid-1)) then
     w=wx*(1.-wy)
!$OMP CRITICAL 
      drygridunc2(ix,jyp,ks,kp,nunc,nage)= &
            drygridunc2(ix,jyp,ks,kp,nunc,nage)+deposit(ks)*w
!$OMP END CRITICAL
  endif

  endif

    end do

  endif !kernel
end subroutine drydepokernel
