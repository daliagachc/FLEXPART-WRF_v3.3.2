module gf2xe
!===============================================================================
! Fortran 90/95 Module for GF(2)[x] computation
!===============================================================================
  use mt_kind_defs
  implicit none
  private
  public :: gf2x_obj
  public :: gf2x_prime_obj
  public :: new,delete
  public :: print_bit,print_hex
  public :: get_deg
  public :: set_coef, set_prime
  public :: assign
  public ::  add,  add_assign
  public :: mult, mult_assign
  public ::  pow, square
  public :: div, rem, divrem
  public :: mult_by_x, div_by_x, mod_by_x
  public :: shift
  public :: deg_i32, mult_i32, square_i32, shift_i32
  public :: mult_i32_old
  public :: gf2x_pow_pow_2

  integer(INT32), parameter :: MAX_KARA  = 64

  type gf2x_obj
    integer(INT32), pointer :: c(:) => NULL()
    integer(INT32) :: deg  = -1
    integer(INT32) :: size = -1
  end type

  type gf2x_prime_obj
    type(gf2x_obj) :: prime_poly
    type(gf2x_obj) :: barrett_poly
    integer(INT32) :: deg
  end type

  interface new
    module procedure gf2x_new
    module procedure gf2x_delete_prime
  end interface

  interface delete
    module procedure gf2x_delete
    module procedure gf2x_delete_prime
  end interface

  interface print_bit
    module procedure gf2x_print_bit
  end interface

  interface print_hex
    module procedure gf2x_print_hex
  end interface

  interface set_coef
    module procedure gf2x_set_coef
  end interface

  interface set_prime
    module procedure gf2x_set_prime
  end interface

  interface assign
    module procedure gf2x_assign
  end interface

  interface add
    module procedure gf2x_add
  end interface

  interface add_assign
    module procedure gf2x_add_assign
  end interface

  interface mult
    module procedure gf2x_mult_kara
  end interface

  interface mult_assign
    module procedure gf2x_mult_assign_kara
  end interface

  interface pow
    module procedure gf2x_pow
    module procedure gf2x_pow_mod
  end interface

  interface square
    module procedure gf2x_square
  end interface

  interface mult_by_x
    module procedure gf2x_mult_by_x
  end interface

  interface mod_by_x
    module procedure gf2x_mod_by_x
  end interface

  interface div_by_x
    module procedure gf2x_div_by_x
  end interface

  interface div
    module procedure gf2x_div
  end interface

  interface rem
    module procedure gf2x_rem
    module procedure gf2x_rem_barrett
  end interface

  interface divrem
    module procedure gf2x_divrem
  end interface

  interface shift
    module procedure gf2x_shift
  end interface

contains

!!DEC$ ATTRIBUTES FORCEINLINE :: get_size
function get_size(deg) result(size)
  integer(INT32) :: deg,size
  size = CEILING(real(deg+1,kind=REAL64)/32.0_REAL64)
  return
end function

subroutine gf2x_new(this,deg)
  type(gf2x_obj), intent(inout) :: this
  integer(INT32), intent(in) :: deg
  integer(INT32) :: isize
  intrinsic :: SIZE
  if (deg < 0) then
    this%deg  = -1
    this%size = -1
    return
  endif
  isize = get_size(deg)
  this%size = isize
  this%deg  = deg
  if (.not.associated(this%c)) then
    allocate(this%c(0:isize-1))
  else
    if (SIZE(this%c) < this%size) then
      deallocate(this%c)
      NULLIFY(this%c)
      allocate(this%c(0:isize-1))
    endif
  endif
  this%c(:) = 0
  return
end subroutine

subroutine gf2x_delete(this)
  type(gf2x_obj), intent(inout) :: this
  integer(INT32) :: ierr
  if (associated(this%c)) then
    deallocate(this%c,STAT=ierr)
  endif
  NULLIFY(this%c)
  this%deg  = -1
  this%size = -1
  return
end subroutine

subroutine gf2x_print_bit(this)
  type(gf2x_obj), intent(in) :: this
  integer(INT32) :: i,ib,iw,deg
  deg = get_deg(this)
  if (deg < 0) then
    write(*,'("0")')
    return
  endif
  do i=deg,0,-1
    ib = mod(i,32)
    iw = i/32
    if (BTEST(this%c(iw),ib)) then
      write(*,'("1",$)')
    else
      write(*,'("0",$)')
    endif
  enddo
  write(*,'("")')
  return
end subroutine

subroutine gf2x_print_hex(this)
  type(gf2x_obj), intent(in) :: this
  integer(INT32) :: i,ib,iw,isize
  character(9) :: str
  if (is_zero(this)) then
    write(*,'("0")')
    return
  endif
  isize = get_size(this%deg)
  i = isize-1
  write(str,'(Z8)')this%c(i)
  write(*,'(A,$)')TRIM(ADJUSTL(str))
  do i=isize-2,0,-1
    write(*,'(Z8.8,$)')this%c(i)
  enddo
  write(*,'("")')
  return
end subroutine

subroutine gf2x_assign(c,a)
  type(gf2x_obj), intent(inout) :: c  ! c := a
  type(gf2x_obj), intent(in)    :: a
  integer(INT32) :: ia,isa,i

  call delete(c)
  if (is_zero(a)) then
    return
  endif

  ia = get_deg(a)
  isa = get_size(ia)
  call new(c,ia)
  do i=0,isa-1
    c%c(i) = a%c(i)
  enddo
  
  return
end subroutine

function is_zero(a) result(is)
  type(gf2x_obj), intent(in) :: a
  logical :: is
  integer(INT32) :: deg
  deg = get_deg(a) 
  if (deg==-1) then
    is = .true.
  else
    is = .false.
  endif
  return
end function

!!DEC$ ATTRIBUTES FORCEINLINE :: get_deg
function get_deg(a) result(deg)
  type(gf2x_obj), intent(in) :: a
  integer(INT32) :: deg
  integer(INT32) :: isize,i,top_deg
  intrinsic :: SIZE
  deg=-1
  if (.not.associated(a%c)) return
  isize = SIZE(a%c)
  do i=isize-1,0,-1
    if (a%c(i) /= 0) then
      top_deg = deg_i32(a%c(i))
      deg = 32*i + top_deg
      return
    endif
  enddo
  return
end function

subroutine gf2x_set_coef(a,i)
  type(gf2x_obj), intent(inout) :: a
  integer(INT32), intent(in) :: i
  type(gf2x_obj), pointer :: w
  integer(INT32) :: ib,iw
  NULLIFY(w)
  if (is_zero(a)) then
    call new(a,i)
  endif
  allocate(w)
  call new(w,i)
  iw =     i/32
  ib = mod(i,32)
  w%c(iw) = ibset(w%c(iw),ib)
  call add_assign(w,a)  ! w := w + a
  call assign(a,w)      ! a := w
  call delete(w)
  deallocate(w)
  NULLIFY(w)
  a%deg = get_deg(a)
  return
end subroutine

subroutine gf2x_add_assign(c,a)
  type(gf2x_obj), intent(inout) :: c  ! c := c + a
  type(gf2x_obj), intent(in)    :: a
  type(gf2x_obj), pointer :: w
  integer(INT32) :: ia,ic
  integer(INT32) :: isa,isc,i
  if (is_zero(a)) then
    return
  endif
  if (is_zero(c)) then
    call assign(c,a)
    return
  endif
  ia = a%deg
  ic = c%deg
  isa = a%size
  isc = c%size
  if (isc < isa) then
    NULLIFY(w)
    allocate(w)
    call new(w,MAX(ia,ic))
    do i=0,isc-1
      w%c(i) = IEOR(c%c(i),a%c(i))
    enddo
    do i=isc,isa-1
      w%c(i) = a%c(i)
    enddo
    call assign(c,w)
    call delete(w)
    deallocate(w)
    NULLIFY(w)
  else
    do i=0,isa-1
      c%c(i) = IEOR(c%c(i),a%c(i))
    enddo
    c%deg = get_deg(c)
    c%size = get_size(c%deg)
  endif
  return
end subroutine

subroutine gf2x_add(c,a,b)
  type(gf2x_obj), intent(inout) :: c   ! c := a + b
  type(gf2x_obj), intent(in)    :: a,b
  integer(INT32) :: ia,ib,ic
  integer(INT32) :: isa,isb,isc,i
  if (is_zero(a) .and. is_zero(b)) then
    return
  endif
  if (is_zero(a)) then
    call assign(c,b)
    return
  endif
  if (is_zero(b)) then
    call assign(c,a)
    return
  endif
  ia = get_deg(a)
  ib = get_deg(b)
  isa = get_size(ia)
  isb = get_size(ib)
  if (c%deg < MAX(ia,ib)) call new(c,MAX(ia,ib))
  if (isa < isb) then
    do i=0,isa-1
      c%c(i) = IEOR(a%c(i),b%c(i))
    enddo
    do i=isa,isb-1
      c%c(i) = b%c(i)
    enddo
  else
    do i=0,isb-1
      c%c(i) = IEOR(a%c(i),b%c(i))
    enddo
    do i=isb,isa-1
      c%c(i) = a%c(i)
    enddo
  endif
  c%deg  = get_deg(c)
  c%size = get_size(c%deg)
  return
end subroutine

subroutine gf2x_pow(c,a,e)
  type(gf2x_obj), intent(inout) :: c ! c = a**e
  type(gf2x_obj), intent(in)    :: a
  integer(INT32), intent(in) :: e
  type(gf2x_obj), pointer :: w
  integer(INT32) :: ch,cl
  integer(INT32) :: i,deg
  NULLIFY(w)
  call delete(c)
  if (e==1) then
    call assign(c,a)
    return
  endif
  if (e==0) then
    call set_coef(c,0)
    return
  endif
  if (e<0) then
    write(*,*)"pow: c = a^e : exponent should be e>=0."
    stop
  endif
  if (is_zero(a)) return

  deg = deg_i32(e)

  allocate(w)
  call set_coef(c,0)
  do i=deg,0,-1
    call square(w,c)        ! w := c**2
    if (BTEST(e,i)) then
      call mult(c,w,a)      ! c := w * a
    else
      call assign(c,w)      ! c := w
    endif
  enddo
  call delete(w)
  deallocate(w)
  NULLIFY(w)

  return
end subroutine

subroutine gf2x_square(c,a)
  type(gf2x_obj), intent(inout) :: c ! c := a**2
  type(gf2x_obj), intent(in)    :: a
  integer(INT32) :: ch,cl
  integer(INT32) :: i,deg
  call delete(c)
  if (is_zero(a)) return
  deg = a%deg*2
  call new(c,deg)
  do i=0,a%size-1
    if (a%c(i) == 0) cycle
    call square_i32(a%c(i),ch,cl)
    if (cl /= 0) c%c(2*i)   = IEOR(c%c(2*i),  cl)
    if (ch /= 0) c%c(2*i+1) = IEOR(c%c(2*i+1),ch)
  enddo
  c%deg = get_deg(c)
  c%size = get_size(c%deg)
  return
end subroutine

recursive subroutine gf2x_mult_kara(c,a,b)
!
! multiply 2 polyomials using Karatsuba algorithm
!
  type(gf2x_obj), intent(inout) :: c    ! c := a * b
  type(gf2x_obj), intent(in)    :: a,b
  type(gf2x_obj), pointer :: ah,al,bh,bl,ahbh,albl,ahl,bhl,ahlbhl
  integer(INT32) :: isa,isb,isc
  integer(INT32) :: i,j,deg

  NULLIFY(ah,al,bh,bl,ahbh,albl,ahl,bhl,ahlbhl)
  call delete(c)
  if (is_zero(a)) return
  if (is_zero(b)) return

  isa = a%size
  isb = b%size
  isc = MAX(isa,isb)
  if (isc < MAX_KARA) then
    call gf2x_mult_normal(c,a,b)
    return
  endif

  if (mod(isc,2)/=0) then
    isc = isc + 1
  endif

  allocate(ah,al,bh,bl,ahbh,albl,ahl,bhl,ahlbhl)
  deg = 32*(isc/2)-1
  call new(al,deg)
  call new(bl,deg)
  call new(ah,deg)
  call new(bh,deg)

  do i=0,MIN(isc/2-1,isa-1)
    al%c(i) = a%c(i)
  enddo
  do i=0,MIN(isc/2-1,isb-1)
    bl%c(i) = b%c(i)
  enddo
  do i=isc/2,isa-1
    ah%c(i-isc/2) = a%c(i)
  enddo
  do i=isc/2,isb-1
    bh%c(i-isc/2) = b%c(i)
  enddo
  ah%deg = get_deg(ah)
  al%deg = get_deg(al)
  bh%deg = get_deg(bh)
  bl%deg = get_deg(bl)
  ah%size = get_size(ah%deg)
  al%size = get_size(al%deg)
  bh%size = get_size(bh%deg)
  bl%size = get_size(bl%deg)

!===================================

  call add(ahl,ah,al)
  call add(bhl,bh,bl)
  call gf2x_mult_kara(ahlbhl,ahl,bhl)
  call delete(ahl)
  call delete(bhl)
  deallocate(ahl,bhl)

!===================================

  call gf2x_mult_kara(ahbh,ah,bh)
  call delete(ah)
  call delete(bh)
  deallocate(ah,bh)

  call add_assign(ahlbhl,ahbh)

!===================================

  call gf2x_mult_kara(albl,al,bl)
  call delete(al)
  call delete(bl)
  deallocate(al,bl)

  call add_assign(ahlbhl,albl)

!===================================
  deg = a%deg + b%deg
  call new(c,deg)

  do i=0,MIN(c%size,albl%size)-1
    c%c(i) = albl%c(i)
  enddo
  call delete(albl)

  if (.not. is_zero(ahlbhl)) then
    do i=isc/2,MIN(c%size,isc/2+ahlbhl%size)-1
      c%c(i) = IEOR(c%c(i),ahlbhl%c(i-isc/2))
    enddo
  endif
  call delete(ahlbhl)

  if (.not. is_zero(ahbh)) then
    do i=isc,MIN(c%size,isc+ahbh%size)-1
      c%c(i) = IEOR(c%c(i),ahbh%c(i-isc))
    enddo
  endif
  call delete(ahbh)
  deallocate(ahbh,albl,ahlbhl)
  NULLIFY(ah,al,bh,bl,ahbh,albl,ahl,bhl,ahlbhl)
  c%deg  = get_deg(c)
  c%size = get_size(c%deg)
  return
end subroutine

subroutine gf2x_mult_assign_kara(a,b)
  type(gf2x_obj), intent(inout) :: a  ! a := a * b
  type(gf2x_obj), intent(in)    :: b
  type(gf2x_obj), pointer :: w
  NULLIFY(w)
  if (is_zero(a)) then
    call delete(a)
    return
  endif
  if (is_zero(b)) then
    call delete(a)
    return
  endif
  allocate(w)
  call gf2x_mult_kara(w,a,b)
  call assign(a,w)
  call delete(w)
  deallocate(w)
  NULLIFY(w)
  return
end subroutine

subroutine gf2x_mult_assign_normal(a,b)
  type(gf2x_obj), intent(inout) :: a  ! a := a * b
  type(gf2x_obj), intent(in)    :: b
  type(gf2x_obj), pointer :: w
  integer(INT32) :: ch,cl
  integer(INT32) :: i,j,deg
  NULLIFY(w)
  allocate(w)
  deg = a%deg + b%deg
  call new(w,deg)
  call gf2x_mult_normal(w,a,b)
  call assign(a,w)  ! a := w
  call delete(w)
  deallocate(w)
  NULLIFY(w)
  return
end subroutine

subroutine gf2x_mult_normal(c,a,b)
  type(gf2x_obj), intent(inout) :: c    ! c := a * b
  type(gf2x_obj), intent(in)    :: a,b
  integer(INT32) :: ch,cl
  integer(INT32) :: i,j,ij,deg,kk,mm
  integer(INT32), allocatable :: hi(:,:),lo(:,:)

  call delete(c)
  if (is_zero(a) .or. is_zero(b) ) then
    return
  endif

  deg = a%deg + b%deg
  call new(c,deg)

!#define _NEW_
!#undef _NEW_
!#ifdef _NEW_
!  do j=0,c%size-2
!    kk = MIN(j,  a%size-1)
!    mm = MAX(0,j-b%size+1)
!    do i=mm,kk
!      call mult_i32(a%c(i),b%c(j-i),ch,cl)
!      c%c(j)   = IEOR(c%c(j),  cl)
!      c%c(j+1) = IEOR(c%c(j+1),ch)
!    enddo
!  enddo
!  j=c%size-1
!  kk = a%size-1
!  mm = c%size-b%size
!  do i=mm,kk
!    call mult_i32(a%c(i),b%c(j-i),ch,cl)
!    c%c(j)   = IEOR(c%c(j),  cl)
!  enddo
!
!#else
  do j=0,b%size-1
  if (b%c(j) == 0) cycle
  do i=0,a%size-1
  if (a%c(i) == 0) cycle

    ij = i + j
    call mult_i32(a%c(i),b%c(j),ch,cl)
                     c%c(ij)   = IEOR(c%c(ij),  cl)
    if (ij+1<c%size) c%c(ij+1) = IEOR(c%c(ij+1),ch)

  enddo
  enddo
!#endif

  c%deg = get_deg(c)
  c%size = get_size(c%deg)

  return
end subroutine

subroutine gf2x_shift(c,a,i)
  type(gf2x_obj), intent(inout) :: c  ! c := shift(a,i)
  type(gf2x_obj), intent(in)    :: a
  integer(INT32), intent(in) :: i
  integer(INT32) :: j,isn,iw,ib,ida,isa,ch,cm,cl
  if (i==0) then
    call assign(c,a)
    return
  endif
  ida = get_deg(a)
  isa = get_size(ida)
  if (ida + i < 0) then
    call delete(c)
    return
  endif
  iw = abs(i)/32
  ib = mod(abs(i),32)
  call delete(c)
  call new(c,ida+i)
  if (i > 0) then
    do j=0,isa-1
      call shift_i32(a%c(j),+ib,ch,cm,cl)
      if (ch /= 0) c%c(j+iw+1) = IEOR(c%c(j+iw+1),ch)
      if (cm /= 0) c%c(j+iw)   = IEOR(c%c(j+iw)  ,cm)
    enddo
  else 
    call shift_i32(a%c(iw),-ib,ch,cm,cl)
    if (cm /= 0) c%c(0)   = IEOR(c%c(0),cm)
    do j=iw+1,isa-1
      call shift_i32(a%c(j),-ib,ch,cm,cl)
      if (cm /= 0) c%c(j-iw)   = IEOR(c%c(j-iw)  ,cm)
      if (cl /= 0) c%c(j-iw-1) = IEOR(c%c(j-iw-1),cl)
    enddo
  endif
  c%deg = get_deg(c)
  c%size = get_size(c%deg)
  return
end subroutine

subroutine gf2x_divrem(q,r,a,b)
 ! a =: q * b + r
  type(gf2x_obj), intent(inout) :: q  ! q := a div b
  type(gf2x_obj), intent(inout) :: r  ! r := a mod b
  type(gf2x_obj), intent(in)    :: a,b
  type(gf2x_obj), pointer :: w,t,s
  integer(INT32) :: ida,idb,idw
  call delete(q)
  call delete(r)
  ida = a%deg
  idb = b%deg
  if (ida < idb) then
    call assign(r,a)
    return
  endif
  NULLIFY(w,t,s)
  allocate(w,t,s)
  call assign(w,a)
  idw = w%deg
  do
    call mult_by_x(t,b,idw-idb)  ! t := b * x^(deg(w)-deg(b))
    call set_coef(s,idw-idb)     ! s := s + x^(deg(w)-deg(b))
    call add_assign(w,t)         ! w := w + t
    call delete(t)
    idw = w%deg
    if (idw < idb) exit
  enddo
  call assign(r,w)
  call delete(w)
  call assign(q,s)
  call delete(s)
  deallocate(w,t,s)
  NULLIFY(w,t,s)
  return
end subroutine

subroutine gf2x_div(q,a,b)
 ! a =: q * b + r
  type(gf2x_obj), intent(inout) :: q  ! q := a div b
  type(gf2x_obj), intent(in)    :: a,b
  type(gf2x_obj), pointer :: w,t,s
  integer(INT32) :: ida,idb,idw
  call delete(q)
  ida = a%deg
  idb = b%deg
  if (ida < idb) then
    return
  endif
  NULLIFY(w,t,s)
  allocate(w,t,s)
  call assign(w,a)
  idw = w%deg
  do
    call mult_by_x(t,b,idw-idb)  ! t := b * x^(deg(w)-deg(b))
    call set_coef(s,idw-idb)     ! s := s + x^(deg(w)-deg(b))
    call add_assign(w,t)         ! w := w + t
    call delete(t)
    idw = w%deg
    if (idw < idb) exit
  enddo
  call delete(w)
  call assign(q,s)
  call delete(s)
  deallocate(w,t,s)
  NULLIFY(w,t,s)
  return
end subroutine

subroutine gf2x_rem(r,a,b)
  type(gf2x_obj), intent(inout) :: r   ! r := a mod b
  type(gf2x_obj), intent(in)    :: a,b
  type(gf2x_obj), pointer :: w,t
  integer(INT32) :: ida,idb,idw
  call delete(r)
  ida = a%deg
  idb = b%deg
  if (ida < idb) then
    call assign(r,a)
    return
  endif
  NULLIFY(w,t)
  allocate(w,t)
  call assign(w,a)
  idw = w%deg
  do
    call mult_by_x(t,b,idw-idb)  ! t := b * x^(deg(w)-deg(b))
    call add_assign(w,t)         ! w := w + t
    call delete(t)
    idw = w%deg
    if (idw < idb) exit
  enddo
  call assign(r,w)
  call delete(w)
  deallocate(w,t)
  NULLIFY(w,t)
  return
end subroutine

subroutine gf2x_set_prime(mp,m)
!
! Set a primitive polynomial to the cotainer.
! the container contains the prime poly and precomputed polynomial for Barrett reduciont.
! This routine does not check the primitivity.
!  mp : container
!   m : primitive polynomial
!
  type(gf2x_prime_obj), intent(inout) :: mp
  type(gf2x_obj),       intent(in)    :: m
  type(gf2x_obj), pointer :: xx
  integer(INT32) :: deg
  call delete(mp)
  call assign(mp%prime_poly,m)
  deg = get_deg(m)
  mp%deg = deg
  NULLIFY(xx)
  allocate(xx)
  call set_coef(xx,2*deg)
  call div(mp%barrett_poly,xx,m)
  call delete(xx)
  deallocate(xx)
  NULLIFY(xx)
  return
end subroutine

subroutine gf2x_delete_prime(mp)
  type(gf2x_prime_obj), intent(inout) :: mp
  call delete(mp%prime_poly)
  call delete(mp%barrett_poly)
  return
end subroutine

subroutine gf2x_rem_barrett(r,a,m)
!
! compute  r := a mod m using Barrett algorithm
!
  type(gf2x_obj), intent(inout) :: r     ! r := a mod m
  type(gf2x_obj), intent(in)    :: a
  type(gf2x_prime_obj), intent(in) :: m  ! precomputed polynomial for Barrett algorithm
  type(gf2x_obj), pointer :: q,p
  integer(INT32) :: deg

  call delete(r)
  deg = m%deg
  if (a%deg < deg) then
    call assign(r,a)
    return
  endif

  NULLIFY(q,p)
  allocate(q,p)

  call div_by_x(q,a,deg)              ! q = a  /  x**deg
  call mod_by_x(r,a,deg)              ! r = a mod x**deg

  call mult_assign(q,m%barrett_poly)  ! q = q  *  mu
  call div_by_x(p,q,deg)              ! p = q  /  x**deg

  call mult_assign(p,m%prime_poly)    ! p = p  *  m
  call mod_by_x(q,p,deg)              ! q = p mod x**deg

  call add_assign(r,q)                ! r = r + q

  call delete(p)
  call delete(q)

  deallocate(p,q)
  NULLIFY(p,q)

  return
end subroutine

subroutine gf2x_mod_by_x(c,a,i)
  type(gf2x_obj), intent(inout) :: c  ! c := a mod x^i
  type(gf2x_obj), intent(in)    :: a
  integer(INT32), intent(in) :: i
  type(gf2x_obj), pointer :: w
  integer(INT32) :: iw,ib,j
  call delete(c)
  if (a%deg < i) then
    call assign(c,a)
    return
  endif
  if (i == 0) then
    call delete(c)
    return
  endif
  if (i < 0) then
    write(*,'("mod_by_x: error, negative i:",I10)')i
    stop
  endif
  iw = i/32
  ib = mod(i,32)
  NULLIFY(w)
  allocate(w)
  call new(w,i)
  do j=0,w%size-1
    w%c(j) = a%c(j)
  enddo
  w%c(w%size-1) = IAND(w%c(w%size-1),2**ib-1)
  call assign(c,w)
  call delete(w)
  deallocate(w)
  NULLIFY(w)
  return
end subroutine

subroutine gf2x_mult_by_x(c,a,i)
  type(gf2x_obj), intent(inout) :: c  ! c := a * x^i
  type(gf2x_obj), intent(in)    :: a
  integer(INT32), intent(in) :: i
  if (i < 0) then
    write(*,'("mult_by_x: error, negative i:",I10)')i
    stop
  endif
  if (i == 0) then
    call assign(c,a)
    return
  endif
  call shift(c,a,i)
  return
end subroutine

subroutine gf2x_div_by_x(c,a,i)
  type(gf2x_obj), intent(inout) :: c  ! c := a div x^i
  type(gf2x_obj), intent(in)    :: a
  integer(INT32), intent(in) :: i
  if (i < 0) then
    write(*,'("div_by_x: error, negative i:",I10)')i
    stop
  endif
  if (i == 0) then
    call assign(c,a)
    return
  endif
  call shift(c,a,-i)
  return
end subroutine


subroutine gf2x_pow_pow_2(c,e,m)
  type(gf2x_obj), intent(inout) :: c  ! c := x**(2**e) mod m
  integer(INT32), intent(in)           :: e
  type(gf2x_prime_obj), intent(in) :: m  ! precomputed polynomial for Barrett algorithm
  integer(INT32) :: i,ee
  type(gf2x_obj), pointer :: w,s

  ee = CEILING(log(REAL(m%deg))/log(2.0))
  call delete(c)
  if (ee > e) then
    call set_coef(c,2**e)
    return
  endif

  NULLIFY(w,s)
  allocate(w,s)
  call set_coef(w,2**ee)
  call rem(s,w,m)      ! s = w mod m
  do i=ee+1,e
    call square(w,s)   ! w = s**2
    call rem(s,w,m)    ! s = w mod m
  enddo
  call assign(c,s)
  call delete(w)
  call delete(s)
  deallocate(w,s)
  NULLIFY(w,s)
  return
end subroutine

subroutine gf2x_pow_mod(c,a,e,m)
  type(gf2x_obj), intent(inout) :: c  ! c := a**e mod m
  type(gf2x_obj), intent(in)    :: a
  integer(INT32), intent(in)           :: e
  type(gf2x_prime_obj), intent(in) :: m  ! precomputed polynomial for Barrett algorithm
  type(gf2x_obj), pointer :: w
  integer(INT32) :: i,deg
  NULLIFY(w)
  call delete(c)
  if (e==1) then
    if (a%deg >= m%deg) then
       call rem(c,a,m)
       return
    else
      call assign(c,a)
      return
    endif
  endif
  if (e==0) then
    call set_coef(c,0)
    return
  endif
  if (e<0) then
    write(*,*)"pow: c = a^e mod m : exponent should be e>=0."
    stop
  endif
  if (is_zero(a)) return

  deg = deg_i32(e)

  allocate(w)
  call set_coef(c,0)
  do i=deg,0,-1
    call square(w,c)        ! c := c**2 mod m
    call rem(c,w,m)
    if (BTEST(e,i)) then
      call mult(w,c,a)      ! c := c * a mod m
      call rem(c,w,m)
    endif
  enddo
  call delete(w)
  deallocate(w)
  NULLIFY(w)
  return
end subroutine

!========================================================================
function deg_i32(a) result(d)
  integer(INT32) :: a,d,i
  d=-1
  do i=31,0,-1
    if (BTEST(a,i)) then
      d=i
      exit
    endif
  enddo
  return
end function

function deg_i64(a) result(d)
  integer(INT64) :: a
  integer(INT32) :: d,i
  do i=63,0,-1
    if (BTEST(a,i)) then
      d=i
      exit
    endif
  enddo
  return
end function

subroutine square_i32(a,ch,cl)
  integer(INT32), intent(in) :: a
  integer(INT32), intent(out) :: ch,cl   ! (ch,cl) = a**2
  integer(INT32) :: ia,i
  integer(INT64) :: da,dc
  da = a
  if (da < 0) da = da + 2_8**32 ! convert to unsigned
  dc = Z'0'
  ia = deg_i32(a)
  do i = 0,ia
    if (BTEST(a,i)) then
      dc = ibset(dc,i*2)
    endif
  enddo
  ch = ISHFT(dc,-32)
  cl = dc
  return
end subroutine

!DEC$ ATTRIBUTES FORCEINLINE :: mult_i32
subroutine mult_i32(a,b,ch,cl)
  integer(INT32), intent(in) :: a,b
  integer(INT32), intent(out) :: ch,cl  ! (ch,cl) = a*b
  integer(INT32) :: tmp,u(0:3)
  integer(INT32), parameter :: ZE = Z'eeeeeeee'
  integer(INT32), parameter :: ZC = Z'cccccccc'
  integer(INT32), parameter :: Z8 = Z'88888888'

  if (a==0 .or. b ==0) then
    ch = 0
    cl = 0
    return
  endif

  u(0) = 0
  u(1) = a
  u(2) = ISHFT(u(1),+1)
  u(3) =  IEOR(u(2),a)

  cl =                  IEOR(ISHFT(u(     ISHFT(b,-30)   ),2),u(IAND(ISHFT(b,-28),3)))
  ch =                  ISHFT(cl,-28)
  cl = IEOR(ISHFT(cl,4),IEOR(ISHFT(u(IAND(ISHFT(b,-26),3)),2),u(IAND(ISHFT(b,-24),3))))
  ch =  IOR(ISHFT(ch,4),ISHFT(cl,-28))
  cl = IEOR(ISHFT(cl,4),IEOR(ISHFT(u(IAND(ISHFT(b,-22),3)),2),u(IAND(ISHFT(b,-20),3))))
  ch =  IOR(ISHFT(ch,4),ISHFT(cl,-28))
  cl = IEOR(ISHFT(cl,4),IEOR(ISHFT(u(IAND(ISHFT(b,-18),3)),2),u(IAND(ISHFT(b,-16),3))))
  ch =  IOR(ISHFT(ch,4),ISHFT(cl,-28))
  cl = IEOR(ISHFT(cl,4),IEOR(ISHFT(u(IAND(ISHFT(b,-14),3)),2),u(IAND(ISHFT(b,-12),3))))
  ch =  IOR(ISHFT(ch,4),ISHFT(cl,-28))
  cl = IEOR(ISHFT(cl,4),IEOR(ISHFT(u(IAND(ISHFT(b,-10),3)),2),u(IAND(ISHFT(b, -8),3))))
  ch =  IOR(ISHFT(ch,4),ISHFT(cl,-28))
  cl = IEOR(ISHFT(cl,4),IEOR(ISHFT(u(IAND(ISHFT(b, -6),3)),2),u(IAND(ISHFT(b, -4),3))))
  ch =  IOR(ISHFT(ch,4),ISHFT(cl,-28))
  cl = IEOR(ISHFT(cl,4),IEOR(ISHFT(u(IAND(ISHFT(b, -2),3)),2),u(IAND(      b     ,3))))

  tmp = -IAND(ISHFT(a,-31),1)
  tmp =  IAND(tmp,ISHFT(IAND(b,ZE),-1))
  ch  =  IEOR(ch,tmp)

  tmp = -IAND(ISHFT(a,-30),1)
  tmp =  IAND(tmp,ISHFT(IAND(b,ZC),-2))
  ch  =  IEOR(ch,tmp)

  tmp = -IAND(ISHFT(a,-29),1)
  tmp =  IAND(tmp,ISHFT(IAND(b,Z8),-3))
  ch  =  IEOR(ch,tmp)

  return
end subroutine


subroutine mult_i32_old(a,b,ch,cl)
  integer(INT32), intent(in) :: a,b
  integer(INT32), intent(out) :: ch,cl  ! (ch,cl) = a*b
  integer(INT32) :: ia,ib,i
  integer(INT64) :: da,db,dc
  da = a
  db = b
  if (da < 0) da = da + 2_8**32 ! convert to unsigned
  if (db < 0) db = db + 2_8**32 ! convert to unsigned
  ia = deg_i32(a)
  ib = deg_i32(b)
  dc = Z'0'
  do i = 0,ia
    if (BTEST(a,i)) then
      dc = IEOR(dc,db)
    endif
    dc = ISHFTC(dc,-1)
  enddo
  dc = ISHFTC(dc,ia+1)
  ch = ISHFT(dc,-32)
  cl = dc
!  write(*,'(B64.64)')dc
!  write(*,'(B32.32)')ch
!  write(*,'(B64.64)')cl
  return
end subroutine

subroutine shift_i32(a,i,ch,cm,cl)
  integer(INT32), intent(in) :: a
  integer(INT32), intent(in) :: i
  integer(INT32), intent(out) :: ch,cm,cl  ! (ch,cm,cl) = shift(a,i)
  integer(INT64) :: dc
  if (abs(i) >= 32) then
    write(*,*)"shift_int32: error i=",i
    stop
  endif
  select case (i)
  case (0)
    ch = 0; cm = a; cl = 0
    return
  case (1:31)
    dc = a
    if (dc < 0) dc = dc + 2_8**32 ! convert to unsigned
    dc = ISHFT(dc,i)
    ch = ISHFT(dc,-32)
    cm = dc
    cl = 0
    return
  case (-31:-1)
    dc = a
    if (dc < 0) dc = dc + 2_8**32 ! convert to unsigned
    dc = ISHFT(dc,i+32)
    ch = 0
    cm = ISHFT(dc,-32)
    cl = dc
    return
  end select
  return
end subroutine

end module
