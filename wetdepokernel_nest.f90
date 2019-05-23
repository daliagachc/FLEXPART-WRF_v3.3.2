!***********************************************************************
!* Copyright 2012,2013                                                *
!* Jerome Brioude, Delia Arnold, Andreas Stohl, Wayne Angevine,       *
!* John Burkhart, Massimo Cassiani, Adam Dingwell, Richard C Easter, Sabine Eckhardt,*
!* Stephanie Evan, Jerome D Fast, Don Morton, Ignacio Pisso,          *
!* Petra Seibert, Gerard Wotawa, Caroline Forster, Harald Sodemann,   *
!*                                                                     *
!* This file is part of FLEXPART WRF         
!                                                                     *
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

subroutine wetdepokernel_nest &
       (nunc,deposit,x,y,itage,nage,kp)
  !        i      i    i i  i    i
  !*****************************************************************************
  !                                                                            *
  !     Attribution of the deposition from an individual particle to the       *
  !     nested deposition fields using a uniform kernel with bandwidths        *
  !     dxoutn and dyoutn.                                                     *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     26 December 1996                                                       *
  !                                                                            *
  !      2 September 2004: Adaptation from wetdepokernel.                      *
  !                                                                            *
  !     Dec 2012: modifications following wetdepokernel.f90                   *
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
  integer :: ix,jy,ixp,jyp,ks,kp,nunc,nage
! CDA new declarations
   real :: rhoprof(2),rhoi,xlon,ylat,xl2,yl2
   integer :: itage
! CDA


!JB
  if (outgrid_option.eq.0) then
! CDA
  xl=(x*dx+xoutshiftn)/dxoutn
  yl=(y*dy+youtshiftn)/dyoutn
  elseif (outgrid_option.eq.1) then
! CDA new code:
  xl2=x*dx+xmet0
  yl2=y*dy+ymet0
  call xymeter_to_ll_wrf(xl2,yl2,xlon,ylat)
  xl=(xlon-outlon0n)/dxoutln
  yl=(ylat-outlat0n)/dyoutln
  endif

  ix=int(xl)
  jy=int(yl)

! CDA skip kernelfor some hours to prevent smoothing close to source
      if (itage.lt.7200) then ! no kernel, direct attribution to grid cell
        do ks=1,nspec
          if ((abs(deposit(ks)).gt.0).and.DRYDEPSPEC(ks)) then
!        print*,ix,jy,ks,kp,nunc,nage
            if ((ix.ge.0).and.(jy.ge.0).and.(ix.le.numxgridn-1).and. &
                                      (jy.le.numygridn-1)) &
        wetgriduncn(ix,jy,ks,kp,nunc,nage)= &
          wetgriduncn(ix,jy,ks,kp,nunc,nage)+deposit(ks)
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

  if ((ix.ge.0).and.(jy.ge.0).and.(ix.le.numxgridn-1).and. &
       (jy.le.numygridn-1)) then
    w=wx*wy
      wetgriduncn(ix,jy,ks,kp,nunc,nage)= &
           wetgriduncn(ix,jy,ks,kp,nunc,nage)+deposit(ks)*w
  endif

  if ((ixp.ge.0).and.(jyp.ge.0).and.(ixp.le.numxgridn-1).and. &
       (jyp.le.numygridn-1)) then
    w=(1.-wx)*(1.-wy)
      wetgriduncn(ixp,jyp,ks,kp,nunc,nage)= &
           wetgriduncn(ixp,jyp,ks,kp,nunc,nage)+deposit(ks)*w
  endif

  if ((ixp.ge.0).and.(jy.ge.0).and.(ixp.le.numxgridn-1).and. &
       (jy.le.numygridn-1)) then
    w=(1.-wx)*wy
      wetgriduncn(ixp,jy,ks,kp,nunc,nage)= &
           wetgriduncn(ixp,jy,ks,kp,nunc,nage)+deposit(ks)*w
  endif

  if ((ix.ge.0).and.(jyp.ge.0).and.(ix.le.numxgridn-1).and. &
       (jyp.le.numygridn-1)) then
    w=wx*(1.-wy)
      wetgriduncn(ix,jyp,ks,kp,nunc,nage)= &
           wetgriduncn(ix,jyp,ks,kp,nunc,nage)+deposit(ks)*w
  endif

  end do
  endif !kernel

end subroutine wetdepokernel_nest
