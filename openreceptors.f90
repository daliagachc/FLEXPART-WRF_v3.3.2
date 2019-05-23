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

subroutine openreceptors
!*******************************************************************************
!                                                                              *
!  Note:  This is the FLEXPART_WRF version of subroutine openreceptors.        *
!                                                                              *
!  This routine opens the receptor output files and writes out the receptor    *
!  names and the receptor locations. The receptor output files are not closed, *
!  but kept open throughout the simulation. Concentrations are continuously    *
!  dumped to these files.                                                      *
!                                                                              *
!     Author: A. Stohl                                                         *
!                                                                              *
!     7 August 2002                                                            *
!                                                                              *
!     Dec 2005, J. Fast - Output files can be either binary or ascii.          *
!                         Write iomode_xycoord to output files.                *
!                         Receptor positions can be lat-lon or grid-meters.    *
!     Dec 2005, R. Easter - changed names of "*lon0*" & "*lat0*" variables     *
!     2015-05-07  H.M.J. Barbosa: outgrid modes where mixed up.                *
!     2015-05-07  A. Dingwell: Fixed indentation to match standards            *
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
! numreceptor            actual number of receptor points specified            *
! receptornames          names of the receptor points                          *
! xreceptor,yreceptor    coordinates of the receptor points                    *
!                                                                              *
!*******************************************************************************

  use par_mod
  use com_mod

  implicit none

  integer :: j
  real :: xtmp(maxreceptor),ytmp(maxreceptor)
  real :: xtmpb,ytmpb


! Open output file for receptor points and write out a short header
! containing receptor names and locations
!******************************************************************

  if (numreceptor.ge.1) then           ! do it only if receptors are specified

    do j = 1, numreceptor
      xtmp(j) = xreceptor(j)*dx + xmet0
      ytmp(j) = yreceptor(j)*dy + ymet0
  !hmjb, outgrid_option==0=wrf grid(meters), 1=regular lat/lon grid
  ! I think the call below is wrong. opt=1 is output lat/lon, but the function to 
  ! convert the output is called only in case of opt=0
  !            if (outgrid_option .eq. 0) then
      if (outgrid_option .eq. 1) then
        xtmpb = xtmp(j)
        ytmpb = ytmp(j)
        call xymeter_to_ll_wrf( xtmpb, ytmpb, xtmp(j), ytmp(j) )
      endif
    enddo
    
    ! Concentration output
    !*********************
    if ((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)) then
      if (iouttype.eq.0) then 
        open(unitoutrecept,file=path(1)(1:length(1))//'receptor_conc', &
        form='unformatted',err=997)
        write(unitoutrecept) (receptorname(j),j=1,numreceptor)
        write(unitoutrecept) (xtmp(j),ytmp(j),j=1,numreceptor), outgrid_option
      endif
      if (iouttype.eq.1) then 
        open(unitoutrecept,file=path(1)(1:length(1))//'receptor_conc', &
        form='formatted',err=997)
!         do j = 1, numreceptor
!           write(unitoutrecept,*) receptorname(j)
!         enddo
!         write(unitoutrecept,*) (xtmp(j),ytmp(j),j=1,numreceptor),
        do j=1,numreceptor
          write(unitoutrecept,1001) receptorname(j), &
            xtmp(j),ytmp(j),outgrid_option 
        enddo
      endif
    endif
    
    ! Mixing ratio output
    !********************
    if ((iout.eq.2).or.(iout.eq.3)) then
      if (iouttype.eq.0) then 
        open(unitoutreceptppt,file=path(1)(1:length(1))//'receptor_pptv', &
        form='unformatted',err=998)
        write(unitoutreceptppt) (receptorname(j),j=1,numreceptor)
        write(unitoutreceptppt) (xtmp(j),ytmp(j),j=1,numreceptor), &
          outgrid_option
      endif
      if (iouttype.eq.1) then 
        open(unitoutreceptppt,file=path(1)(1:length(1))//'receptor_pptv', &
          form='formatted',err=998)
!         do j = 1, numreceptor
!           write(unitoutreceptppt,*) receptorname(j)
!         enddo
!         write(unitoutreceptppt,*) (xtmp(j),ytmp(j),j=1,numreceptor),
        do j=1,numreceptor
          write(unitoutreceptppt,1001) receptorname(j), &
            xtmp(j),ytmp(j),outgrid_option
        enddo
      endif
    endif
  endif
  return

1001  format(a16,f10.4,f10.4,i5)

997   write(*,*) ' #### FLEXPART MODEL ERROR! THE FILE           #### '
      write(*,*) ' ####              receptor_conc               #### '
      write(*,*) ' #### CANNOT BE OPENED.                        #### '
      stop

998   write(*,*) ' #### FLEXPART MODEL ERROR! THE FILE           #### '
      write(*,*) ' ####              receptor_pptv               #### '
      write(*,*) ' #### CANNOT BE OPENED.                        #### '
      stop

end subroutine openreceptors

