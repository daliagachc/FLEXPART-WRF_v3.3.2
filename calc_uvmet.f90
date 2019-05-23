!! Diagnostics: U & V on earth coordinates
! from ARWpost postprocessing routine from the WRF package.

  SUBROUTINE calc_uvmet(UUU,VVV,SCRa, SCRb, i3dflag)

  IMPLICIT NONE

  !Arguments
  real, allocatable, dimension(:,:,:)             :: SCRa, SCRb
  character (len=128)                             :: cname, cdesc, cunits

  !Local
  integer                                         :: i, j, k
  integer                                         :: i3dflag
  real                                            :: cone
  real, dimension(west_east_dim,south_north_dim)  :: diff, alpha

  cname    = "uvmet"
  cdesc    = "Rotated wind component"
  cunits   = "m s-1"
  
   use com_mod

  IF ( map_proj .ge. 3 ) THEN     ! No need to rotate
    IF ( i3dflag == 1 ) THEN  
      SCRa = UUU
      SCRb = VVV
    ENDIF
    IF ( i3dflag == 0 ) THEN
      SCRa(:,:,1) = U10(:,:)
      SCRb(:,:,1) = V10(:,:)
    END IF
    RETURN
  END IF


  cone = 1.                                          !  PS
  IF ( map_proj .eq. 1) THEN                         !  Lambert Conformal mapping
    IF (ABS(truelat1-truelat2) .GT. 0.1) THEN
       cone=(ALOG(COS(truelat1*RAD_PER_DEG))-            &
             ALOG(COS(truelat2*RAD_PER_DEG))) /          &
       (ALOG(TAN((90.-ABS(truelat1))*RAD_PER_DEG*0.5 ))- &
        ALOG(TAN((90.-ABS(truelat2))*RAD_PER_DEG*0.5 )) )
    ELSE
       cone = SIN(ABS(truelat1)*RAD_PER_DEG )
    ENDIF
  END IF


  diff = XLONG - stand_lon
  DO i = 1, west_east_dim
  DO j = 1, south_north_dim
    IF ( diff(i,j) .gt. 180. ) THEN
      diff(i,j) = diff(i,j) - 360.
    END IF
    IF ( diff(i,j) .lt. -180. ) THEN
      diff(i,j) = diff(i,j) + 360.
    END IF
  END DO
  END DO


  DO i = 1, west_east_dim
  DO j = 1, south_north_dim
     IF ( XLAT(i,j) .lt. 0. ) THEN
       alpha(i,j) = - diff(i,j) * cone * RAD_PER_DEG
     ELSE
       alpha(i,j) = diff(i,j) * cone * RAD_PER_DEG
     END IF
  END DO
  END DO

  

  IF ( i3dflag == 1 ) THEN
    DO k = 1,bottom_top_dim
      SCRa(:,:,k) = VVV(:,:,k)*sin(alpha) + UUU(:,:,k)*cos(alpha)
      SCRb(:,:,k) = VVV(:,:,k)*cos(alpha) - UUU(:,:,k)*sin(alpha)
    END DO
  ELSE
     SCRa(:,:,1) = V10(:,:)*sin(alpha) + U10(:,:)*cos(alpha)
     SCRb(:,:,1) = V10(:,:)*cos(alpha) - U10(:,:)*sin(alpha)
  END IF

  END SUBROUTINE calc_uvmet

