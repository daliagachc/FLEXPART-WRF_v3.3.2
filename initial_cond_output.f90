!**********************************************************************
! Copyright 1998,1999,2000,2001,2002,2005,2007,2008,2009,2010         *
! Andreas Stohl, Petra Seibert, A. Frank, Gerhard Wotawa,             *
! Caroline Forster, Sabine Eckhardt, John Burkhart, Harald Sodemann   *
!                                                                     *
! This file is part of FLEXPART.                                      *
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

subroutine initial_cond_output(itime)
  !                                 i
  !*****************************************************************************
  !                                                                            *
  !     Output of the initial condition sensitivity field.                     *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     24 May 1995                                                            *
  !                                                                            *
  !     13 April 1999, Major update: if output size is smaller, dump output    *
  !                    in sparse matrix format; additional output of           *
  !                    uncertainty                                             *
  !                                                                            *
  !     05 April 2000, Major update: output of age classes; output for backward*
  !                    runs is time spent in grid cell times total mass of     *
  !                    species.                                                *
  !                                                                            *
  !     17 February 2002, Appropriate dimensions for backward and forward runs *
  !                       are now specified in file par_mod                    *
  !                                                                            *
  !     June 2006, write grid in sparse matrix with a single write command     *
  !                in order to save disk space                                 *
  !                                                                            *
  !     2008 new sparse matrix format                                          *
  !                                                                            *
  !     2015-04-30  H.M.J. Barbosa: corrected file handles from 97 to          *
  !                   unitoutgrid                                              *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! ncells          number of cells with non-zero concentrations               *
  ! sparse          .true. if in sparse matrix format, else .false.            *
  !                                                                            *
  !*****************************************************************************

  use unc_mod
  use point_mod
  use outg_mod
  use par_mod
  use com_mod

  implicit none

  integer :: itime,i,ix,jy,kz,ks,kp,sp_count_i,sp_count_r
  real :: sp_fact,fact_recept
  real,parameter :: smallnum = tiny(0.0) ! smallest number that can be handled
  logical :: sp_zer
  character(len=3) :: anspec


  !*********************************************************************
  ! Determine the standard deviation of the mean concentration or mixing
  ! ratio (uncertainty of the output) and the dry and wet deposition
  !*********************************************************************

  do ks=1,nspec
    write(anspec,'(i3.3)') ks
    open(unitoutgrid,file=path(1)(1:length(1))//'grid_initial'// &
         '_'//anspec,form='unformatted')
    write(unitoutgrid) itime

    do kp=1,maxpointspec_act

      if (ind_rel.eq.1) then
        fact_recept=rho_rel(kp)
      else
        fact_recept=1.
      endif

  !*******************************************************************
  ! Generate output: may be in concentration (ng/m3) or in mixing
  ! ratio (ppt) or both
  ! Output the position and the values alternated multiplied by
  ! 1 or -1, first line is number of values, number of positions
  ! For backward simulations, the unit is seconds, stored in grid_time
  !*******************************************************************

  ! Write out dummy "wet and dry deposition" fields, to keep same format
  ! as for concentration output
      sp_count_i=0
      sp_count_r=0
      write(unitoutgrid) sp_count_i
      write(unitoutgrid) (sparse_dump_i(i),i=1,sp_count_i)
      write(unitoutgrid) sp_count_r
      write(unitoutgrid) (sparse_dump_r(i),i=1,sp_count_r)
      write(unitoutgrid) sp_count_i
      write(unitoutgrid) (sparse_dump_i(i),i=1,sp_count_i)
      write(unitoutgrid) sp_count_r
      write(unitoutgrid) (sparse_dump_r(i),i=1,sp_count_r)


  ! Write out sensitivity to initial conditions
      sp_count_i=0
      sp_count_r=0
      sp_fact=-1.
      sp_zer=.true.
      do kz=1,numzgrid
        do jy=0,numygrid-1
          do ix=0,numxgrid-1
            if (init_cond(ix,jy,kz,ks,kp).gt.smallnum) then
              if (sp_zer.eqv..true.) then ! first non zero value
                sp_count_i=sp_count_i+1
                sparse_dump_i(sp_count_i)= &
                     ix+jy*numxgrid+kz*numxgrid*numygrid
                sp_zer=.false.
                sp_fact=sp_fact*(-1.)
              endif
              sp_count_r=sp_count_r+1
              sparse_dump_r(sp_count_r)=sp_fact* &
                   init_cond(ix,jy,kz,ks,kp)/xmass(kp,ks)*fact_recept
            else ! concentration is zero
              sp_zer=.true.
            endif
          end do
        end do
      end do
      write(unitoutgrid) sp_count_i
      write(unitoutgrid) (sparse_dump_i(i),i=1,sp_count_i)
      write(unitoutgrid) sp_count_r
      write(unitoutgrid) (sparse_dump_r(i),i=1,sp_count_r)


    end do

    close(unitoutgrid)

  end do


end subroutine initial_cond_output
