!***********************************************************************
!* Copyright 2012,2013                                                 *
!* Jerome Brioude, Delia Arnold, Jerome Fast, 
!* Andreas Stohl, Petra Seibert, A. Frank, Gerhard Wotawa *
!* Caroline Forster, Sabine Eckhardt, John Burkhart, Harald Sodemann   *
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
!
! J. Brioude: the version in fortran 9 gives some error. better use f77
!
c To be used, if the non-standard Fortran function erf does not exist on
c your machine
C
C     aus:  Numerical Recipes (FORTRAN) / Chapter 6.                       
C                                                                   
C     6.1  FUNCTION GAMMLN                                         
C     6.2  FUNCTION GAMMP   <6.2:GSER/6.2:GCF/6.1:GAMMLN>         
C     6.2  FUNCTION GAMMQ   <6.2:GSER/6.2:GCF/6.1:GAMMLN>        
C     6.2  SUBROUTINE GSER    <6.1:GAMMLN>                      
C     6.2  SUBROUTINE GCF     <6.1:GAMMLN>                     
C     6.2  FUNCTION ERF     <6.2:GAMMP/6.2:GSER/6.2:GCF/6.1:GAMMLN> 
C     6.2  FUNCTION ERFC    <6.2.:GAMMP/6.2:GAMMQ/6.2:GSER/        
C                            6.2:GCF/6.1:GAMMLN>                  
C     6.2  FUNCTION ERFCC                                             
C
      FUNCTION GAMMLN(XX)
      REAL*8 COF(6),STP,HALF,ONE,FPF,X,TMP,SER
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     *    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
11    CONTINUE
      GAMMLN=TMP+LOG(STP*SER)
      RETURN
      END
C
      FUNCTION GAMMP(A,X)
      IF(X.LT.0..OR.A.LE.0.) STOP 'GAMMP'
      IF(X.LT.A+1.)THEN
        CALL GSER(GAMSER,A,X,GLN)
        GAMMP=GAMSER
      ELSE
        CALL GCF(GAMMCF,A,X,GLN)
        GAMMP=1.-GAMMCF
      ENDIF
      RETURN
      END
C
      FUNCTION GAMMQ(A,X)
      IF(X.LT.0..OR.A.LE.0.) STOP 'GAMMQ'
      IF(X.LT.A+1.)THEN
        CALL GSER(GAMSER,A,X,GLN)
        GAMMQ=1.-GAMSER
      ELSE
        CALL GCF(GAMMCF,A,X,GLN)
        GAMMQ=GAMMCF
      ENDIF
      RETURN
      END
C
      SUBROUTINE GSER(GAMSER,A,X,GLN)
      PARAMETER (ITMAX=100,EPS=3.E-7)
      GLN=GAMMLN(A)
      IF(X.LE.0.)THEN
!       IF(X.LT.0.) PAUSE 'GSER'
        GAMSER=0.
        RETURN
      ENDIF
      AP=A
      SUM=1./A
      DEL=SUM
      DO 11 N=1,ITMAX
        AP=AP+1.
        DEL=DEL*X/AP
        SUM=SUM+DEL
        IF(ABS(DEL).LT.ABS(SUM)*EPS)GO TO 1
11    CONTINUE
      STOP 'GSER: A too large, ITMAX too small'
1     GAMSER=SUM*EXP(-X+A*LOG(X)-GLN)
      RETURN
      END
C
      SUBROUTINE GCF(GAMMCF,A,X,GLN)
      PARAMETER (ITMAX=100,EPS=3.E-7)
      GLN=GAMMLN(A)
      GOLD=0.
      A0=1.
      A1=X
      B0=0.
      B1=1.
      FAC=1.
      DO 11 N=1,ITMAX
        AN=FLOAT(N)
        ANA=AN-A
        A0=(A1+A0*ANA)*FAC
        B0=(B1+B0*ANA)*FAC
        ANF=AN*FAC
        A1=X*A0+ANF*A1
        B1=X*B0+ANF*B1
        IF(A1.NE.0.)THEN
          FAC=1./A1
          G=B1*FAC
          IF(ABS((G-GOLD)/G).LT.EPS)GO TO 1
          GOLD=G
        ENDIF
11    CONTINUE
      STOP 'GCF: A too large, ITMAX too small'
1     GAMMCF=EXP(-X+A*ALOG(X)-GLN)*G
      RETURN
      END
C
      FUNCTION erf(X)
      IF(X.LT.0.)THEN
        erf=-GAMMP(.5,X**2)
      ELSE
        erf=GAMMP(.5,X**2)
      ENDIF
      RETURN
      END
C
      FUNCTION ERFC(X)
      IF(X.LT.0.)THEN
        ERFC=1.+GAMMP(.5,X**2)
      ELSE
        ERFC=GAMMQ(.5,X**2)
      ENDIF
      RETURN
      END
C
      FUNCTION ERFCC(X)
      Z=ABS(X)
      T=1./(1.+0.5*Z)
      ERFCC=T*EXP(-Z*Z-1.26551223+T*(1.00002368+T*(.37409196+
     *    T*(.09678418+T*(-.18628806+T*(.27886807+T*(-1.13520398+
     *    T*(1.48851587+T*(-.82215223+T*.17087277)))))))))
      IF (X.LT.0.) ERFCC=2.-ERFCC
      RETURN
      END
C
