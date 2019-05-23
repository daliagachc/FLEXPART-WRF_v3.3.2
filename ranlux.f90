  
    MODULE luxury

!     Subtract-and-borrow random number generator proposed by
!     Marsaglia and Zaman, implemented by F. James with the name
!     RCARRY in 1991, and later improved by Martin Luescher
!     in 1993 to produce "Luxury Pseudorandom Numbers".
!     Fortran 77 coded by F. James, 1993

!  References:
!  M. Luscher, Computer Physics Communications  79 (1994) 100
!  F. James, Computer Physics Communications 79 (1994) 111

!   LUXURY LEVELS.
!   ------ ------      The available luxury levels are:

!  level 0  (p=24): equivalent to the original RCARRY of Marsaglia
!           and Zaman, very long period, but fails many tests.
!  level 1  (p=48): considerable improvement in quality over level 0,
!           now passes the gap test, but still fails spectral test.
!  level 2  (p=97): passes all known tests, but theoretically still
!           defective.
!  level 3  (p=223): DEFAULT VALUE.  Any theoretically possible
!           correlations have very small chance of being observed.
!  level 4  (p=389): highest possible luxury, all 24 bits chaotic.

!!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!!!  Calling sequences for RANLUX:                                  ++
!!!!      CALL RANLUX (RVEC, LEN)   returns a vector RVEC of LEN     ++
!!!!                   32-bit random floating point numbers between  ++
!!!!                   zero (not included) and one (also not incl.). ++
!!!!      CALL RLUXGO(LUX,INT,K1,K2) initializes the generator from  ++
!!!!               one 32-bit integer INT and sets Luxury Level LUX  ++
!!!!               which is integer between zero and MAXLEV, or if   ++
!!!!               LUX .GT. 24, it sets p=LUX directly.  K1 and K2   ++
!!!!               should be set to zero unless restarting at a break++
!!!!               point given by output of RLUXAT (see RLUXAT).     ++
!!!!      CALL RLUXAT(LUX,INT,K1,K2) gets the values of four integers++
!!!!               which can be used to restart the RANLUX generator ++
!!!!               at the current point by calling RLUXGO.  K1 and K2++
!!!!               specify how many numbers were generated since the ++
!!!!               initialization with LUX and INT.  The restarting  ++
!!!!               skips over  K1+K2*E9   numbers, so it can be long.++
!!!!   A more efficient but less convenient way of restarting is by: ++
!!!!      CALL RLUXIN(ISVEC)    restarts the generator from vector   ++
!!!!                   ISVEC of 25 32-bit integers (see RLUXUT)      ++
!!!!      CALL RLUXUT(ISVEC)    outputs the current values of the 25 ++
!!!!                 32-bit integer seeds, to be used for restarting ++
!!!!      ISVEC must be dimensioned 25 in the calling program        ++
!!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

IMPLICIT NONE

INTEGER            :: iseeds(24), isdext(25)
INTEGER, PARAMETER :: maxlev = 4, lxdflt = 3
INTEGER            :: ndskip(0:maxlev) = (/ 0, 24, 73, 199, 365 /)
INTEGER            :: next(24), igiga = 1000000000, jsdflt = 314159265
REAL, PARAMETER    :: twop12 = 4096.
INTEGER, PARAMETER :: itwo24 = 2**24, icons = 2147483563
INTEGER            :: luxlev = lxdflt, nskip, inseed, jseed
LOGICAL            :: notyet = .true.
INTEGER            :: in24 = 0, kount = 0, mkount = 0, i24 = 24, j24 = 10
REAL               :: seeds(24), carry = 0., twom24, twom12

!                            default
!  Luxury Level     0   1   2  *3*    4
!    ndskip        /0, 24, 73, 199, 365/
! Corresponds to p=24  48  97  223  389
!     time factor   1   2   3    6   10   on slow workstation
!                   1 1.5   2    3    5   on fast mainframe

PUBLIC notyet, i24, j24, carry, seeds, twom24, twom12, luxlev
PUBLIC nskip, ndskip, in24, next, kount, mkount, inseed


CONTAINS


SUBROUTINE ranlux(rvec, lenv)

IMPLICIT NONE

INTEGER, INTENT(IN) :: lenv
REAL, INTENT(OUT)   :: rvec(lenv)

!     Local variables

INTEGER             :: i, k, lp, ivec, isk
REAL                :: uni

!  NOTYET is .TRUE. if no initialization has been performed yet.
!              Default Initialization by Multiplicative Congruential

IF (notyet) THEN
  notyet = .false.
  jseed = jsdflt
  inseed = jseed
  WRITE (6,'(A,I12)') ' RANLUX DEFAULT INITIALIZATION: ', jseed
  luxlev = lxdflt
  nskip = ndskip(luxlev)
  lp = nskip + 24
  in24 = 0
  kount = 0
  mkount = 0
  WRITE (6,'(A,I2,A,I4)') ' RANLUX DEFAULT LUXURY LEVEL =  ', luxlev,   &
                          '    p =', lp
  twom24 = 1.
  DO i = 1, 24
    twom24 = twom24 * 0.5
    k = jseed / 53668
    jseed = 40014 * (jseed-k*53668) - k * 12211
    IF (jseed.LT.0) jseed = jseed + icons
    iseeds(i) = MOD(jseed,itwo24)
  END DO
  twom12 = twom24 * 4096.
  DO i = 1, 24
    seeds(i) = REAL(iseeds(i)) * twom24
    next(i) = i - 1
  END DO
  next(1) = 24
  i24 = 24
  j24 = 10
  carry = 0.
  IF (seeds(24).EQ.0.) carry = twom24
END IF

!          The Generator proper: "Subtract-with-borrow",
!          as proposed by Marsaglia and Zaman,
!          Florida State University, March, 1989

DO ivec = 1, lenv
  uni = seeds(j24) - seeds(i24) - carry
  IF (uni.LT.0.) THEN
    uni = uni + 1.0
    carry = twom24
  ELSE
    carry = 0.
  END IF
  seeds(i24) = uni
  i24 = next(i24)
  j24 = next(j24)
  rvec(ivec) = uni
!  small numbers (with less than 12 "significant" bits) are "padded".
  IF (uni.LT.twom12) THEN
    rvec(ivec) = rvec(ivec) + twom24 * seeds(j24)
!        and zero is forbidden in case someone takes a logarithm
    IF (rvec(ivec).EQ.0.) rvec(ivec) = twom24 * twom24
  END IF
!        Skipping to luxury.  As proposed by Martin Luscher.
  in24 = in24 + 1
  IF (in24.EQ.24) THEN
    in24 = 0
    kount = kount + nskip
    DO isk = 1, nskip
      uni = seeds(j24) - seeds(i24) - carry
      IF (uni.LT.0.) THEN
        uni = uni + 1.0
        carry = twom24
      ELSE
        carry = 0.
      END IF
      seeds(i24) = uni
      i24 = next(i24)
      j24 = next(j24)
    END DO
  END IF
END DO
kount = kount + lenv
IF (kount.GE.igiga) THEN
  mkount = mkount + 1
  kount = kount - igiga
END IF
RETURN

END SUBROUTINE ranlux


!           Subroutine to input and float integer seeds from previous run
SUBROUTINE rluxin
!     the following IF BLOCK added by Phillip Helbig, based on conversation
!     with Fred James; an equivalent correction has been published by James.

IMPLICIT NONE

!     Local variables

INTEGER             :: i, isd

IF (notyet) THEN
  WRITE (6,'(A)') ' Proper results ONLY with initialisation from 25 ',  &
  'integers obtained with RLUXUT'
  notyet = .false.
END IF

twom24 = 1.
DO i = 1, 24
  next(i) = i - 1
  twom24 = twom24 * 0.5
END DO
next(1) = 24
twom12 = twom24 * 4096.
WRITE (6,'(A)') ' FULL INITIALIZATION OF RANLUX WITH 25 INTEGERS:'
WRITE (6,'(5X,5I12)') isdext
DO i = 1, 24
  seeds(i) = REAL(isdext(i)) * twom24
END DO
carry = 0.
IF (isdext(25).LT.0) carry = twom24
isd = IABS(isdext(25))
i24 = MOD(isd,100)
isd = isd / 100
j24 = MOD(isd,100)
isd = isd / 100
in24 = MOD(isd,100)
isd = isd / 100
luxlev = isd
IF (luxlev.LE.maxlev) THEN
  nskip = ndskip(luxlev)
  WRITE (6,'(A,I2)') ' RANLUX LUXURY LEVEL SET BY RLUXIN TO: ', luxlev
ELSE IF (luxlev.GE.24) THEN
  nskip = luxlev - 24
  WRITE (6,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXIN TO:', luxlev
ELSE
  nskip = ndskip(maxlev)
  WRITE (6,'(A,I5)') ' RANLUX ILLEGAL LUXURY RLUXIN: ', luxlev
  luxlev = maxlev
END IF
inseed = -1
RETURN

END SUBROUTINE rluxin


!                    Subroutine to ouput seeds as integers
SUBROUTINE rluxut

IMPLICIT NONE

!     Local variables

INTEGER             :: i

DO i = 1, 24
  isdext(i) = INT(seeds(i)*twop12*twop12)
END DO
isdext(25) = i24 + 100 * j24 + 10000 * in24 + 1000000 * luxlev
IF (carry.GT.0.) isdext(25) = -isdext(25)
RETURN

END SUBROUTINE rluxut


!                    Subroutine to output the "convenient" restart point
SUBROUTINE rluxat(lout, inout, k1, k2)

IMPLICIT NONE

INTEGER, INTENT(OUT) :: lout, inout, k1, k2

lout = luxlev
inout = inseed
k1 = kount
k2 = mkount
RETURN

END SUBROUTINE rluxat


!                    Subroutine to initialize from one or three integers
SUBROUTINE rluxgo(lux, ins, k1, k2)

IMPLICIT NONE

INTEGER, INTENT(IN) :: lux, ins, k1, k2

!     Local variables

INTEGER             :: ilx, i, iouter, isk, k, inner, izip, izip2
REAL                :: uni

IF (lux.LT.0) THEN
  luxlev = lxdflt
ELSE IF (lux.LE.maxlev) THEN
  luxlev = lux
ELSE IF (lux.LT.24.OR.lux.GT.2000) THEN
  luxlev = maxlev
  WRITE (6,'(A,I7)') ' RANLUX ILLEGAL LUXURY RLUXGO: ', lux
ELSE
  luxlev = lux
  DO ilx = 0, maxlev
    IF (lux.EQ.ndskip(ilx)+24) luxlev = ilx
  END DO
END IF
IF (luxlev.LE.maxlev) THEN
  nskip = ndskip(luxlev)
  WRITE (6,'(A,I2,A,I4)') ' RANLUX LUXURY LEVEL SET BY RLUXGO :', luxlev,  &
                          '     P=', nskip + 24
ELSE
  nskip = luxlev - 24
  WRITE (6,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXGO TO:', luxlev
END IF
in24 = 0
IF (ins.LT.0) WRITE (6,'(A)') &
              ' Illegal initialization by RLUXGO, negative input seed'
IF (ins.GT.0) THEN
  jseed = ins
  WRITE (6,'(A,3I12)') ' RANLUX INITIALIZED BY RLUXGO FROM SEEDS', jseed, k1, k2
ELSE
  jseed = jsdflt
  WRITE (6,'(A)') ' RANLUX INITIALIZED BY RLUXGO FROM DEFAULT SEED'
END IF
inseed = jseed
notyet = .false.
twom24 = 1.
DO i = 1, 24
  twom24 = twom24 * 0.5
  k = jseed / 53668
  jseed = 40014 * (jseed-k*53668) - k * 12211
  IF (jseed.LT.0) jseed = jseed + icons
  iseeds(i) = MOD(jseed,itwo24)
END DO
twom12 = twom24 * 4096.
DO i = 1, 24
  seeds(i) = REAL(iseeds(i)) * twom24
  next(i) = i - 1
END DO
next(1) = 24
i24 = 24
j24 = 10
carry = 0.
IF (seeds(24).EQ.0.) carry = twom24
!        If restarting at a break point, skip K1 + IGIGA*K2
!        Note that this is the number of numbers delivered to
!        the user PLUS the number skipped (if luxury .GT. 0).
kount = k1
mkount = k2
IF (k1+k2.NE.0) THEN
  DO iouter = 1, k2 + 1
    inner = igiga
    IF (iouter.EQ.k2+1) inner = k1
    DO isk = 1, inner
      uni = seeds(j24) - seeds(i24) - carry
      IF (uni.LT.0.) THEN
        uni = uni + 1.0
        carry = twom24
      ELSE
        carry = 0.
      END IF
      seeds(i24) = uni
      i24 = next(i24)
      j24 = next(j24)
    END DO
  END DO
!         Get the right value of IN24 by direct calculation
  in24 = MOD(kount,nskip+24)
  IF (mkount.GT.0) THEN
    izip = MOD(igiga, nskip+24)
    izip2 = mkount * izip + in24
    in24 = MOD(izip2, nskip+24)
  END IF
!       Now IN24 had better be between zero and 23 inclusive
  IF (in24.GT.23) THEN
    WRITE (6,'(A/A,3I11,A,I5)') &
               '  Error in RESTARTING with RLUXGO:', '  The values', ins, &
               k1, k2, ' cannot occur at luxury level', luxlev
    in24 = 0
  END IF
END IF
RETURN

END SUBROUTINE rluxgo


END MODULE luxury



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine luxtst
!         Exercise for the RANLUX Pseudorandom number generator.

USE luxury

IMPLICIT NONE

REAL    :: rvec(1000)
INTEGER :: i1, i2, i3, i4, li

!         check that we get the right numbers (machine-indep.)
WRITE (6,'(/A)') '  CALL RANLUX(RVEC,100)'
CALL ranlux(rvec,100)
WRITE (6,'(A/9X,5F12.8)') ' RANLUX default numbers   1-  5:', rvec(1:5)
CALL ranlux(rvec,100)
WRITE (6,'(A/9X,5F12.8)') ' RANLUX default numbers 101-105:', rvec(1:5)

WRITE (6,'(/A)') ' CALL RLUXGO(0,0,0,0)'
CALL rluxgo(0,0,0,0)
CALL ranlux(rvec,100)
WRITE (6,'(A/9X,5F12.8)') ' RANLUX luxury level 0,   1-  5:', rvec(1:5)
CALL ranlux(rvec,100)
WRITE (6,'(A/9X,5F12.8)') ' RANLUX luxury level 0, 101-105:', rvec(1:5)

WRITE (6,'(/A)') '   CALL RLUXGO(389,1,0,0)'
CALL rluxgo(389,1,0,0)
CALL ranlux(rvec,100)
WRITE (6,'(A/9X,5F12.8)') ' RANLUX luxury p=389,   1-  5:', rvec(1:5)
CALL ranlux(rvec,100)
WRITE (6,'(A/9X,5F12.8)') ' RANLUX luxury p=389, 101-105:', rvec(1:5)

WRITE (6,'(/A)') '  CALL RLUXGO(75,0,0,0)'
CALL rluxgo(75,0,0,0)
CALL ranlux(rvec,100)
WRITE (6,'(A/9X,5F12.8)') ' RANLUX luxury p= 75,   1-  5:', rvec(1:5)
CALL ranlux(rvec,100)
WRITE (6,'(A/9X,5F12.8)') ' RANLUX luxury p= 75, 101-105:', rvec(1:5)

WRITE (6,'(/A)') '  test restarting from the full vector'
CALL rluxut
WRITE (6,'(/A/(1X,5I14))') '  current RANLUX status saved:', isdext
CALL ranlux(rvec,100)
WRITE (6,'(A/9X,5F12.8)') ' RANLUX numbers 1- 5:', rvec(1:5)
CALL ranlux(rvec,100)
WRITE (6,'(A/9X,5F12.8)') ' RANLUX numbers 101-105:', rvec(1:5)

WRITE (6,'(/A)') '   previous RANLUX status will be restored'
CALL rluxin
CALL ranlux(rvec,100)
WRITE (6,'(A/9X,5F12.8)') ' RANLUX numbers 1- 5:', rvec(1:5)
CALL ranlux(rvec,100)
WRITE (6,'(A/9X,5F12.8)') ' RANLUX numbers 101-105:', rvec(1:5)

WRITE (6,'(/A)') '     test the restarting by skipping'
CALL rluxgo(4,7674985,0,0)
CALL rluxat(i1,i2,i3,i4)
WRITE (6,'(A,4I10)') '  RLUXAT values =', i1, i2, i3, i4
DO li = 1, 10
  CALL ranlux(rvec,1000)
END DO
CALL rluxat(i1,i2,i3,i4)
WRITE (6,'(A,4I10)') '  RLUXAT values =', i1, i2, i3, i4
CALL ranlux(rvec,200)
WRITE (6,'(A,2F10.6)') '  Next and 200th numbers are:', rvec(1), rvec(200)
CALL rluxgo(i1,i2,i3,i4)
CALL ranlux(rvec,200)
WRITE (6,'(A,2F10.6)') '  Next and 200th numbers are:', rvec(1), rvec(200)

WRITE (6,'(/A)') ' The following should provoke an error message'
CALL rluxgo(4,11111,31,0)
STOP

!   OUTPUT FROM THE ABOVE TEST PROGRAM SHOULD BE:
!   --------------------------------------------
!  CALL RANLUX(RVEC,100)
! RANLUX DEFAULT INITIALIZATION:    314159265
! RANLUX DEFAULT LUXURY LEVEL =   3      p = 223
! RANLUX default numbers   1-  5:
!           0.53981817  0.76155043  0.06029940  0.79600263  0.30631220
! RANLUX default numbers 101-105:
!           0.43156743  0.03774416  0.24897110  0.00147784  0.90274453

!  CALL RLUXGO(0,0,0,0)
! RANLUX LUXURY LEVEL SET BY RLUXGO : 0     P=  24
! RANLUX INITIALIZED BY RLUXGO FROM DEFAULT SEED
! RANLUX luxury level 0,   1-  5:
!           0.53981817  0.76155043  0.06029940  0.79600263  0.30631220
! RANLUX luxury level 0, 101-105:
!           0.41538775  0.05330932  0.58195311  0.91397446  0.67034441

!   CALL RLUXGO(389,1,0,0)
! RANLUX LUXURY LEVEL SET BY RLUXGO : 4     P= 389
! RANLUX INITIALIZED BY RLUXGO FROM SEEDS           1           0           0
! RANLUX luxury p=389,   1-  5:
!           0.94589490  0.47347850  0.95152789  0.42971975  0.09127384
! RANLUX luxury p=389, 101-105:
!           0.02618265  0.03775346  0.97274780  0.13302165  0.43126065

!  CALL RLUXGO(75,0,0,0)
! RANLUX P-VALUE SET BY RLUXGO TO:   75
! RANLUX INITIALIZED BY RLUXGO FROM DEFAULT SEED
! RANLUX luxury p= 75,   1-  5:
!           0.53981817  0.76155043  0.06029940  0.79600263  0.30631220
! RANLUX luxury p= 75, 101-105:
!           0.25600731  0.23443210  0.59164381  0.59035838  0.07011414

!  test restarting from the full vector

!  current RANLUX status saved:
!       16156027      16534309      15243811       2751687       6002207
!        7979506       1301976       4567313       4305996       5872599
!       12003090       2146823      12606367       4111505       5979640
!       12739666      10489318      14036909      11729352       8061448
!        7832659       6069758       3197719       1832730      75080216
! RANLUX numbers 1- 5:
!           0.22617835  0.60655993  0.86417443  0.43920082  0.23382509
! RANLUX numbers 101-105:
!           0.08107197  0.21466845  0.84856731  0.94078046  0.85626233

!   previous RANLUX status will be restored
! FULL INITIALIZATION OF RANLUX WITH 25 INTEGERS:
!         16156027    16534309    15243811     2751687     6002207
!          7979506     1301976     4567313     4305996     5872599
!         12003090     2146823    12606367     4111505     5979640
!         12739666    10489318    14036909    11729352     8061448
!          7832659     6069758     3197719     1832730    75080216
! RANLUX P-VALUE SET BY RLUXIN TO:   75
! RANLUX numbers 1- 5:
!           0.22617835  0.60655993  0.86417443  0.43920082  0.23382509
! RANLUX numbers 101-105:
!           0.08107197  0.21466845  0.84856731  0.94078046  0.85626233

!     test the restarting by skipping
! RANLUX LUXURY LEVEL SET BY RLUXGO : 4     P= 389
! RANLUX INITIALIZED BY RLUXGO FROM SEEDS     7674985           0           0
!  RLUXAT values =         4   7674985         0         0
!  RLUXAT values =         4   7674985    161840         0
!  Next and 200th numbers are:  0.019648  0.590586
! RANLUX LUXURY LEVEL SET BY RLUXGO : 4     P= 389
! RANLUX INITIALIZED BY RLUXGO FROM SEEDS     7674985      161840           0
!  Next and 200th numbers are:  0.019648  0.590586

! The following should provoke an error message
! RANLUX LUXURY LEVEL SET BY RLUXGO : 4     P= 389
! RANLUX INITIALIZED BY RLUXGO FROM SEEDS       11111          31           0
!  Error in RESTARTING with RLUXGO:
!  The values      11111         31          0 cannot occur at luxury level    4
END subroutine luxtst

