!***********************************************************************
!* Copyright 2012,2013                                                *
!* Jerome Brioude, Delia Arnold, Andreas Stohl, Wayne Angevine,       *
!* John Burkhart, Massimo Cassiani, Adam Dingwell, Richard C Easter, Sabine Eckhardt,*
!* Stephanie Evan, Jerome D Fast, Don Morton, Ignacio Pisso,          *
!* Petra Seibert, Gerard Wotawa, Caroline Forster, Harald Sodemann,   *
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
    !C--- adapted from press et al. 1992 numerial recipes in Fortran by Massimo Cassiani to jse RANLUX as uniform random number generator
      FUNCTION gasdevlux()
      USE luxury
      INTEGER idum
      REAL gasdevlux     
      INTEGER iset
      REAL fac,gset,rsq,v1,v2,RTEST(2)
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
1      call RANLUX(RTEST,2)

        v1=2.*RTEST(1)-1.
        v2=2.*RTEST(2)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdevlux=v2*fac
        iset=1
      else
        gasdevlux=gset
        iset=0
      endif
      return
      END
      
      subroutine gasdevlux2R(random1,random2)
      USE luxury
      INTEGER idum
      REAL random1,random2      
      INTEGER iset
      REAL fac,gset,rsq,v1,v2,RTEST(2)
      
1      call RANLUX(RTEST,2)

        v1=2.*RTEST(1)-1.
        v2=2.*RTEST(2)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        random1=v1*fac
        random2=v2*fac
  ! Limit the random numbers to lie within the interval -4 and +4
  !**************************************************************
       if (random1.lt.-4.) random1=-4.
       if (random2.lt.-4.) random2=-4.
       if (random1.gt.4.) random1=4.
       if (random2.gt.4.) random2=4.
      return
      END
    
