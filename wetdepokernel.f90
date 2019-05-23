!***********************************************************************
!* Copyright 2012,2013                                                *
!* Jerome Brioude, Delia Arnold, Andreas Stohl, Wayne Angevine,       *
!* John Burkhart, Massimo Cassiani, Adam Dingwell, Richard C Easter, Sabine Eckhardt,*
!* Stephanie Evan, Jerome D Fast, Don Morton, Ignacio Pisso,          *
!* Petra Seibert, Gerard Wotawa, Caroline Forster, Harald Sodemann,   *
!*                                                                     *
!* This file is part of FLEXPART WRF                                   *
!*                                  
! FLEXPART is free software: you can redistribute it and/or modify    *
! it under the terms of the GNU General Public License as published by*
! the Free Software Foundation, either version 3 of the License, or   *
! (at your option) any later version.                                 *
!                                                                     *
! FLEXPART is distributed in the hope that it will be useful,         *
! but WITHOUT ANY WARRANTY; without even the implied warranty of      *
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       *
! GNU General Public License for more details.                        *
!                                                                     *
! You should have received a copy of the GNU General Public License   *
! along with FLEXPART.  If not, see <http://www.gnu.org/licenses/>.   *
!**********************************************************************

subroutine wetdepokernel(nunc,deposit,x,y,itage,nage,kp)
  !                          i      i    i i  i
  !*****************************************************************************
  !                                                                            *
  !     Attribution of the deposition from an individual particle to the       *
  !     deposition fields using a uniform kernel with bandwidths dxout and dyout.*
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     26 December 1996                                                       *
  !
  !     D. Arnold: modification to skip the kernel the first 3 hours.
  !     then, modification to a regular lat-lon.
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
  integer :: ix,jy,ixp,jyp,nunc,nage,ks,kp
! CDA new declarations
   real :: rhoprof(2),rhoi,xlon,ylat,xl2,yl2
   integer :: itage
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
!JB I don't know what this thing is suppose to do, but I put it here in case it is necessary
! see for instance conccalc_reg.f90
!      ix=int(xl)
!      if (xl.lt.0.) ix=ix-1
!      jy=int(yl)
!      if (yl.lt.0.) jy=jy-1
!

! CDA skip kernelfor some hours to prevent smoothing close to source
   if (itage.lt.7200) then ! no kernel, direct attribution to grid cell
    do ks=1,nspec
      if ((ix.ge.0).and.(jy.ge.0).and.(ix.le.numxgrid-1).and. &
                                   (jy.le.numygrid-1)) &
                    wetgridunc(ix,jy,ks,kp,nunc,nage)= &
              wetgridunc(ix,jy,ks,kp,nunc,nage)+deposit(ks)
    enddo   
   else ! attribution via  uniform kernel

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

  if ((ix.ge.0).and.(jy.ge.0).and.(ix.le.numxgrid-1).and. &
       (jy.le.numygrid-1)) then
    w=wx*wy
      wetgridunc(ix,jy,ks,kp,nunc,nage)= &
           wetgridunc(ix,jy,ks,kp,nunc,nage)+deposit(ks)*w
  endif

  if ((ixp.ge.0).and.(jyp.ge.0).and.(ixp.le.numxgrid-1).and. &
       (jyp.le.numygrid-1)) then
    w=(1.-wx)*(1.-wy)
      wetgridunc(ixp,jyp,ks,kp,nunc,nage)= &
           wetgridunc(ixp,jyp,ks,kp,nunc,nage)+deposit(ks)*w
  endif

  if ((ixp.ge.0).and.(jy.ge.0).and.(ixp.le.numxgrid-1).and. &
       (jy.le.numygrid-1)) then
    w=(1.-wx)*wy
      wetgridunc(ixp,jy,ks,kp,nunc,nage)= &
           wetgridunc(ixp,jy,ks,kp,nunc,nage)+deposit(ks)*w
  endif

  if ((ix.ge.0).and.(jyp.ge.0).and.(ix.le.numxgrid-1).and. &
       (jyp.le.numygrid-1)) then
    w=wx*(1.-wy)
      wetgridunc(ix,jyp,ks,kp,nunc,nage)= &
           wetgridunc(ix,jyp,ks,kp,nunc,nage)+deposit(ks)*w
  endif
  end do

  endif ! Kernel
end subroutine wetdepokernel
