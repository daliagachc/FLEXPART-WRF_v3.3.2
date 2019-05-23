subroutine f_get_coeff(nn,mm,rr,ww,avec,nj,id,pp,np)
!===============================================================================
! Compute MT jump ahead polynomial coefficients
! uses GF(2)[x] computation
!===============================================================================
  use mt_kind_defs
  use gf2xe
  implicit none
  integer(INT32), intent(in) :: nn,mm,rr,ww,avec,nj,id
  integer(INT32), intent(inout) :: pp(0:nn-1),np
  type(gf2x_obj) :: af,bf,ff,f1,f2
  type(gf2x_prime_obj) :: fp
  integer(INT32) :: i,ib,nws

!==============================
! MT characteristic polynomial
!  ff : MT char poly.
!==============================
  call set_coef(af,nn)
  call set_coef(af,mm)   ! af = x^nn + x^mm
  call set_coef(bf,nn-1)
  call set_coef(bf,mm-1) ! bf = x^(nn-1) + x^(mm-1)

  call pow(f1,af,ww-rr)  ! f1 = af^(ww-rr)
  call pow(f2,bf,rr)     ! f2 = bf^(rr)
  call mult(ff,f1,f2)    ! ff = f1*f2
  do i=0,rr-1
    ib = mod(i,ww)
    if (BTEST(avec,ib)) then
      call pow(f2,bf,rr-1-i)
      call mult_assign(f2,f1)
      call add_assign(ff,f2)
    endif
  enddo
  do i=rr,ww-1
    ib = mod(i,ww)
    if (BTEST(avec,ib)) then
      call pow(f1,af,ww-1-i)
      call add_assign(ff,f1)
    endif
  enddo

!#ifdef _DEBUG_
!  write(*,'("@",$)')
!  call print_hex(ff)
!#endif

  call delete(af)
  call delete(bf)
  call delete(f1)
  call delete(f2)

!==============================
! set ff for Barrett reduction
!==============================
  call set_prime(fp,ff)  ! fp = ff
  call delete(ff)
  call delete(f1)
  call delete(f2)

!===============
! jump ahead
!  long jump
!===============
  call gf2x_pow_pow_2(f1,nj,fp)   ! f1 = x**(2**nj) mod fp

!#ifdef _DEBUG_
!  write(*,'("@",$)')
!  call print_hex(f1)
!#endif

!===============
! short jump
!===============
  call pow(ff,f1,id,fp) ! ff = f1**id mod fp

!#ifdef _DEBUG_
!  write(*,'("@",$)')
!  call print_hex(ff)
!#endif

  pp(:) = 0
  np = get_deg(ff)+1
  nws = CEILING(real(np,kind=KIND(1.0d0))/32)
  pp(0:nws-1) = ff%c(0:nws-1)

  call delete(f1)
  call delete(f2)
  call delete(ff)
  call delete(fp)

  return
end subroutine
