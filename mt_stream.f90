!
! A Fortran90/95 module program for Multiple Stream Mersenne Twister.
! 2010/02/12, Ken-Ichi Ishikawa [ishikawa[at]theo.phys.sci.hiroshima-u.ac.jp]
! 2010/07/20, M. S. Briggs, [michael.s.briggs[at]nasa.gov], minor improvements.
!
! This module provides the Mersenne Twister pseudo random number generator
! (MT-PRNG) with multiple stream ability. 
! A long period MT generaotr (backbone stream) is divided into disjoint streams,
! and the each stream has continuous stream from the backbone stream.
! The distance between each stream is chosen to have long enough length.
! The status and parameters of the each stream are encapsulated in 
! the f90 type components so that we can use multiple streams simultaneously.
! The stream length is fixed to 2^nj.
!
! This code is converted from original/sample codes located at
! Mersenne Twister Home Page:
!               [http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/mt.html]
!
! For the Original MT PRNG, see also: 
!  M. Matsumoto and T. Nishimura, 
!  "Mersenne Twister: A 623-dimensionally equidistributed uniform pseudorandom
!   number generator", 
!  ACM Trans. on Modeling and Computer Simulation Vol. 8, No. 1, 
!  January pp.3-30 (1998) DOI:10.1145/272991.272995.
!
! For the jump ahead mechanism, see also: 
!  H. Haramoto, M. Matsumoto, T. Nishimura, F. Panneton, and P. L'Ecuyer, 
!  "Efficient Jump Ahead for F_2-Linear Random Number Generators", 
!  GERAD Report G-2006-62. INFORMS Journal on Computing, 20, 3 (2008), 385-390. 
!
! This routine uses;
!  Fast arithmetic in GF(2)[x], [http://wwwmaths.anu.edu.au/~brent/software.html]
!  NTL : A Library for doing Number Theory, [http://www.shoup.net/ntl/index.html]
!
!
! Copyright (c) 2010, Ken-Ichi Ishikawa [ishikawa[at]theo.phys.sci.hiroshima-u.ac.jp]
! All rights reserved.
! 
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are
! met:
! 
! * Redistributions of source code must retain the above copyright
!   notice, this list of conditions and the following disclaimer. 
!   
! * Redistributions in binary form must reproduce the above copyright
!   notice, this list of conditions and the following disclaimer listed
!   in this license in the documentation and/or other materials
!   provided with the distribution.
!   
! * Neither the name of the copyright holders nor the names of its
!   contributors may be used to endorse or promote products derived from
!   this software without specific prior written permission.
!   
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT  
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
! A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
! OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT  
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
! 

!define _TEMPERING_

!ifdef _DEBUG_
!module mt_stream_debug
!else
 module mt_stream
!endif
  use mt_kind_defs
  implicit none
  private
  public :: MTS_SUCCESS
  public :: MTS_FAIL
  public :: mt_state
  public :: new,delete
  public :: init 
  public :: print
  public :: save
  public :: read
  public :: set_mt19937
  public :: create_stream
  public :: genrand_int32
  public :: genrand_double1
  public :: genrand_double2
  public :: genrand_double3
  public :: genrand_double4

  integer(INT32), parameter :: MTS_SUCCESS = 0
  integer(INT32), parameter :: MTS_FAIL    = -1
!#ifdef _DEBUG_
!  integer(INT32), parameter :: MT_JUMP_DISTANCE_EXP = 16   ! test jump distance (2^16 steps)
!#else
  integer(INT32), parameter :: MT_JUMP_DISTANCE_EXP = 256  ! default jump distance (2^256 steps)
!#endif

!======================================
! MT19937 parameters
!======================================
  integer(INT32), parameter :: MT19937_N = 624
  integer(INT32), parameter :: MT19937_M = 397
  integer(INT32), parameter :: MT19937_W = 32
  integer(INT32), parameter :: MT19937_R = 31
  integer(INT32), parameter :: MT19937_MATA  = INT (Z'9908b0df', INT32)
  integer(INT32), parameter :: MT19937_WMASK = INT (Z'ffffffff', INT32)
  integer(INT32), parameter :: MT19937_UMASK = INT (Z'80000000', INT32)
  integer(INT32), parameter :: MT19937_LMASK = INT (Z'7fffffff', INT32)
  integer(INT32), parameter :: MT19937_SHFT0 = 11
  integer(INT32), parameter :: MT19937_SHFT1 = 18
  integer(INT32), parameter :: MT19937_SHFTB =  7
  integer(INT32), parameter :: MT19937_SHFTC = 15
  integer(INT32), parameter :: MT19937_MASKB = INT (Z'9d2c5680', INT32)
  integer(INT32), parameter :: MT19937_MASKC = INT (Z'efc60000', INT32)

!======================================
! MT state
! Period = 2^(N*W-R)-1
!======================================
  type mt_state
    private
    integer(INT32) :: i = -1         ! state vector index
    integer(INT32) :: stream_id = -1 ! stream ID
    integer(INT32) :: istatus = -1   ! initialization status
    integer(INT32) :: nn = -1        ! MT parameter N
    integer(INT32) :: mm = -1        ! MT parameter M
    integer(INT32) :: rr = -1        ! MT parameter R
    integer(INT32) :: ww = -1        ! MT parameter W (width =32)
    integer(INT32) :: aaa = 0        ! Companion matrix parameter
    integer(INT32) :: wmask = 0      ! 32-bit mask
    integer(INT32) :: umask = 0      ! Twist mask x(1)
    integer(INT32) :: lmask = 0      ! Twist mask x(0)
    integer(INT32) :: shift0 = 0     ! Temparing parameters ...
    integer(INT32) :: shift1 = 0
    integer(INT32) :: maskB  = 0
    integer(INT32) :: maskC  = 0
    integer(INT32) :: shiftB = 0
    integer(INT32) :: shiftC = 0
    integer(INT32) :: mag(0:1) = 0   ! mag(0) = 0, mag(1) = aaa
    integer(INT32), pointer :: state(:) => NULL()  ! state vector
  end type

  type(mt_state), save :: g_mt_master ! this keeps MT parameters
  integer(INT32), save :: total_stream = 0

  interface new
    module procedure mt_new
  end interface

  interface delete
    module procedure mt_delete
  end interface

  interface print
    module procedure mt_print
  end interface

  interface read
    module procedure mt_read
  end interface

  interface save
    module procedure mt_save
  end interface

  interface set_mt19937
    module procedure mt_set_mt19937
  end interface

  interface init
    module procedure mt_init_by_scalar
    module procedure mt_init_by_array
  end interface

  interface create_stream
    module procedure mt_create_stream
  end interface

  interface genrand_int32    
   ! in [0,0xFFFFFFFF]
    module procedure mt_genrand_int32
  end interface

  interface genrand_double1
   ! in [0,1]  (53-bit resolution)
    module procedure mt_genrand_double1
  end interface

  interface genrand_double2
   ! in [0,1)  (53-bit resolution)
    module procedure mt_genrand_double2
  end interface

  interface genrand_double3
   ! in (0,1)  (52-bit resolution)
    module procedure mt_genrand_double3
  end interface

  interface genrand_double4
   ! in (0,1]  (53-bit resolution)
    module procedure mt_genrand_double4
  end interface

contains

subroutine mt_new(this,ierr)
!
!= Allocate mt state 
!
  implicit none
  type(mt_state), intent(inout)  :: this
  integer(INT32), optional, intent(out) :: ierr
  character(256), parameter :: myname="mt_new"
  integer(INT32) :: jerr,kerr
  jerr = MTS_SUCCESS
  if (0 == g_mt_master%nn) then
    write(*,'(A,": MT master parameter is not initialized.")')TRIM(myname)
    write(*,'(A,": Stop!")')TRIM(myname)
    stop
  endif
  call mt_copy(g_mt_master,this)
  this%stream_id = total_stream
  total_stream = total_stream + 1
  if (.not.associated(this%state)) then
    if (this%nn <= 0) then
      jerr = MTS_FAIL
      this%istatus = jerr
      goto 100
    endif
    allocate(this%state(0:this%nn-1),STAT=kerr)
    if (kerr /= 0) then
      jerr = MTS_FAIL
      this%istatus = jerr
      goto 100
    endif
    this%istatus = jerr
    goto 100
  endif
100 if (present(ierr)) then
    ierr = jerr
    return
  else
    if (jerr /= MTS_SUCCESS) then
      write(*,'(A,": Something wrong.")')TRIM(myname)
      write(*,'(A,": Stop!")')TRIM(myname)
      stop
    endif
  endif
end subroutine

subroutine mt_delete(this)
!
!= deallocate state vector
!
  implicit none
  type(mt_state), intent(inout) :: this
  if (associated(this%state)) then
    deallocate(this%state)
  endif
  total_stream = total_stream - 1 
  return
end subroutine

subroutine mt_set_mt19937
!
!= set MT19937 parameters on g_mt_master
!
  implicit none
  g_mt_master%nn    = MT19937_N
  g_mt_master%mm    = MT19937_M
  g_mt_master%ww    = MT19937_W
  g_mt_master%rr    = MT19937_R
  g_mt_master%aaa   = MT19937_MATA
  g_mt_master%wmask = MT19937_WMASK
  g_mt_master%umask = MT19937_UMASK
  g_mt_master%lmask = MT19937_LMASK
  g_mt_master%shift0 = MT19937_SHFT0
  g_mt_master%shift1 = MT19937_SHFT1
  g_mt_master%shiftB = MT19937_SHFTB
  g_mt_master%shiftC = MT19937_SHFTC
  g_mt_master%maskB  = MT19937_MASKB
  g_mt_master%maskC  = MT19937_MASKC
  g_mt_master%mag(0) = 0
  g_mt_master%mag(1) = g_mt_master%aaa
  total_stream = 0
  return
end subroutine

subroutine mt_copy(this,that)
!
!= duplicate/copy MT paramters
!
!  that <= this
!
  implicit none
  type(mt_state), intent(in)  :: this
  type(mt_state), intent(out) :: that
  that%i         = this%i
  that%stream_id = this%stream_id
  that%istatus   = this%istatus
  that%nn     = this%nn
  that%nn     = this%nn
  that%mm     = this%mm
  that%ww     = this%ww
  that%aaa    = this%aaa
  that%rr     = this%rr
  that%wmask  = this%wmask
  that%umask  = this%umask
  that%lmask  = this%lmask
  that%shift0 = this%shift0
  that%shift1 = this%shift1
  that%shiftB = this%shiftB
  that%shiftC = this%shiftC
  that%maskB  = this%maskB
  that%maskC  = this%maskC
  that%mag(0) = 0
  that%mag(1) = that%aaa
  return
end subroutine

subroutine mt_save(this,unit)
!
!= save MT status to file with unit
!
  implicit none
  type(mt_state), intent(in) :: this
  integer(INT32), intent(in) :: unit
  write(unit) this%i
  write(unit) this%stream_id
  write(unit) this%istatus
  write(unit) this%nn
  write(unit) this%mm
  write(unit) this%ww
  write(unit) this%aaa
  write(unit) this%rr
  write(unit) this%wmask
  write(unit) this%umask
  write(unit) this%lmask
  write(unit) this%shift0
  write(unit) this%shift1
  write(unit) this%shiftB
  write(unit) this%shiftC
  write(unit) this%maskB
  write(unit) this%maskC
  write(unit) this%mag(:)
  write(unit) this%state(:)
  return
end subroutine

subroutine mt_read(this,unit)
!
!= read MT status from file with unit
!
  implicit none
  type(mt_state), intent(inout) :: this
  integer(INT32), intent(in)    :: unit
  read(unit) this%i
  read(unit) this%stream_id
  read(unit) this%istatus
  read(unit) this%nn
  read(unit) this%mm
  read(unit) this%ww
  read(unit) this%aaa
  read(unit) this%rr
  read(unit) this%wmask
  read(unit) this%umask
  read(unit) this%lmask
  read(unit) this%shift0
  read(unit) this%shift1
  read(unit) this%shiftB
  read(unit) this%shiftC
  read(unit) this%maskB
  read(unit) this%maskC
  read(unit) this%mag(:)
  if (associated(this%state)) then
    deallocate(this%state)
    allocate(this%state(0:this%nn-1))
  endif
  read(unit) this%state(:)
  return
end subroutine

subroutine mt_print(this)
!
!= print MT stream parameter and states
!
  implicit none
  type(mt_state), intent(in)  :: this
  character(256) :: str
  write(*,'("===============================================================================")')
  write(*,'("  MT PARAMETERS")')
  write(str,*)this%stream_id
  write(*,'(" STREAMID= ",A)')TRIM(ADJUSTL(str))
  write(str,*)this%nn*this%ww-this%rr
  write(*,'("   PERIOD= 2^",A," - 1")')TRIM(ADJUSTL(str))
  write(*,'("       NN=",I8," MM=",I8," WW=",I3," RR=",I3)')this%nn,this%mm,this%ww,this%rr
  write(*,'("     AVEC= 0x",Z8.8)', advance='no' )this%aaa
  write(*,'(" WMASK= 0x",Z8.8)', advance='no' )this%wmask
  write(*,'(" UMASK= 0x",Z8.8)', advance='no' )this%umask
  write(*,'(" LMASK= 0x",Z8.8)')  this%lmask
  write(*,'("      MAG= 0x",Z8.8," 0x",Z8.8)') this%mag
  write(*,'("    SHFT0=",I3)', advance='no' )this%shift0
  write(*,'(" SHFT1=",I3)', advance='no' )this%shift1
  write(*,'(" SHFTB=",I3)', advance='no' )this%shiftB
  write(*,'(" SHFTC=",I3)', advance='no' )this%shiftC
  write(*,'("  MASKB= 0x",Z8.8)', advance='no' )this%maskB
  write(*,'(" MASKC= 0x",Z8.8)')  this%maskC
  write(str,*)this%i
  write(*,'("  POINTER= ",A)')TRIM(ADJUSTL(str))
  write(*,'("===============================================================================")')
  return
end subroutine


subroutine mt_init_by_scalar(this,iseed,ierr)
!
!= initialize MT state by a scalar seed.
!
  implicit none
  type(mt_state), intent(inout) :: this
  integer(INT32), intent(in) :: iseed
  integer(INT32), optional, intent(out) :: ierr
  character(256), parameter :: myname="mt_init_by_scalar"
  integer(INT32) :: i,jseed,jerr
  if (.not.associated(this%state)) then
    jerr = MTS_FAIL
    goto 100
  endif
  jseed = iseed
  do i=0,this%nn-1
    this%state(i) = jseed
    jseed = 1812433253 * IEOR(jseed, ISHFT(jseed,-30)) + i + 1
  enddo
  this%i = this%nn
  do i=0,this%nn-1
    this%state(i) = IAND(this%state(i),this%wmask)
  enddo
  this%mag(0) = 0
  this%mag(1) = this%aaa
  jerr = MTS_SUCCESS
100 if (present(ierr)) then
    ierr = jerr
  else
    if (jerr /= MTS_SUCCESS) then
      write(*,'(A,": State vector allocation fails.")')TRIM(myname)
      write(*,'(A,": Stop!")')TRIM(myname)
      stop
    endif
  endif
  return
end subroutine

subroutine mt_init_by_array(this,iseed,ierr)
!
!= initialize MT state by array seeds.
!
  implicit none
  type(mt_state), intent(inout) :: this
  integer(INT32), intent(in)    :: iseed(0:)
  integer(INT32), optional, intent(out) :: ierr
  integer :: isize,i,j,k,n,jerr
  character(256), parameter :: myname="mt_init_by_array"
  integer(INT32), parameter :: ASEED = 19650218
  
  call mt_init_by_scalar(this,ASEED,jerr)
  if (jerr /= MTS_SUCCESS) goto 100

  isize = SIZE(iseed)
  n = this%nn
  i = 1
  j = 0
  do k = MAX(n,isize),1,-1
    this%state(i) = IEOR(this%state(i),(IEOR(this%state(i-1),ISHFT(this%state(i-1),-30))*1664525)) + iseed(j) + j
    this%state(i) = IAND(this%state(i),this%wmask)   ! for WORDSIZE > 32 machines
    i = i + 1
    j = j + 1
    if (i >= n) then
      this%state(0) = this%state(n-1)
      i = 1
    endif
    if (j >= isize) j = 0
  enddo
  do k = n-1,1,-1
    this%state(i) = IEOR(this%state(i),(IEOR(this%state(i-1),ISHFT(this%state(i-1),-30))*1566083941)) - i
    this%state(i) = IAND(this%state(i),this%wmask)   ! for WORDSIZE > 32 machines
    i = i + 1
    if (i >= n) then
      this%state(0) = this%state(n-1)
      i = 1
    endif
  enddo
  this%state(0) = INT (Z'80000000', INT32)    ! MSB is 1; assuring non-zero initial array
  jerr = MTS_SUCCESS
100 if (present(ierr)) then
    ierr = jerr
  else
    if (jerr /= MTS_SUCCESS) then
      write(*,'(A,": mt_init_by_scalar fails.")')TRIM(myname)
      write(*,'(A,": Stop!")')TRIM(myname)
      stop
    endif
  endif
  return
end subroutine

function mt_genrand_int32(this) result(ir)
!
!= return a value in [0,0xFFFFFFFF]
!
  implicit none
  type(mt_state), intent(inout) :: this
  integer(INT32) :: ir
  integer(INT32) :: umask,lmask,n,m,is
  integer(INT32) :: k,nm,n1
  if (this%i >= this%nn) then
    n = this%nn
    m = this%mm
    lmask = this%lmask
    umask = this%umask
    nm = n - m
    n1 = n - 1
    do k=0,nm-1
      is = IOR(IAND(this%state(k),umask),IAND(this%state(k+1),lmask))
      this%state(k) = IEOR(IEOR(this%state(k+m),ISHFT(is,-1)),this%mag(IAND(is,1)))
    enddo
    do k=nm,n1-1
      is = IOR(IAND(this%state(k),umask),IAND(this%state(k+1),lmask))
      this%state(k) = IEOR(IEOR(this%state(k+m-n),ISHFT(is,-1)),this%mag(IAND(is,1)))
    enddo
    is = IOR(IAND(this%state(n-1),umask),IAND(this%state(0),lmask))
    this%state(n-1) = IEOR(IEOR(this%state(m-1),ISHFT(is,-1)),this%mag(IAND(is,1)))
    this%i = 0
  endif

  is = this%state(this%i)
  this%i = this%i + 1
!#ifdef _TEMPERING_
!  is = IEOR(is,ISHFT(is,-this%shift0))
!  is = IEOR(is,IAND(ISHFT(is, this%shiftB),this%maskB))
!  is = IEOR(is,IAND(ISHFT(is, this%shiftC),this%maskC))
!  is = IEOR(is,ISHFT(is,-this%shift1))
!#endif
  ir = is
  return
end function

subroutine mt_matvec(this,v,w)
!
!= Multiply transition matrix on a state vector v
!
!  w = B v
!
!  this : MT parameters(transition matrix)
!     v : input vector
!     w : output vector
!
  implicit none
  type(mt_state), intent(in) :: this
  integer(INT32), intent(in)  :: v(0:this%nn-1)
  integer(INT32), intent(out) :: w(0:this%nn-1)
  integer(INT32) :: umask,lmask,n,m,is
  integer(INT32) :: k
  n = this%nn
  m = this%mm
  lmask = this%lmask
  umask = this%umask
  w(0) = IAND(v(1),umask)
  do k=1,n-2
    w(k) = v(k+1)
  enddo
      is = IOR(IAND(v(0),umask),IAND(v(1),lmask))
  w(n-1) = IEOR(IEOR(v(m),ISHFT(is,-1)),this%mag(IAND(is,1)))
  return
end subroutine

subroutine mt_create_stream(this,that,id)
!
!= Create New stream (that) with distance id*2^(jp) from (this)
!
  implicit none
  type(mt_state), intent(inout) :: this ! input state
  type(mt_state), intent(inout) :: that ! output state
  integer(INT32), intent(in) :: id
  character(256), parameter :: myname="mt_create_stream"
  integer(INT32), parameter :: jp = MT_JUMP_DISTANCE_EXP
  if (id < 0) then
    write(*,'(A,": Positive ID is requried.")')TRIM(myname)
    write(*,'(A,": Stop!")')TRIM(myname)
    stop
  endif
  call mt_jumpahead(this,that,jp,id)
  that%stream_id = id
  total_stream = total_stream + 1 
  this%i = this%nn
  return
end subroutine

subroutine mt_jumpahead(this,that,jp,id)
!
!= Jump ahead by (id+1) * 2^jp steps.
!
!  this : input state
!  that : output state, proceeds by id*2^jp steps.
!
  implicit none
  type(mt_state), intent(inout) :: this ! input state
  type(mt_state), intent(inout) :: that ! output state
  integer(INT32), intent(in) :: jp     ! exponent (jump step = id*2^jp)
  integer(INT32), intent(in) :: id     ! id       (jump step = id*2^jp)
  integer(INT32) :: v(0:this%nn-1)
  integer(INT32) :: w(0:this%nn-1)
  integer(INT32) :: s(0:this%nn-1)
  integer(INT32) :: i,iwp,ibp
  character(256) :: str,str2
  integer(INT32) :: p(0:this%nn-1)
  integer(INT32) :: np
  !
  ! external routine in 
  ! jump_ahead_coeff/get_coeff.cxx
  !    written in C++ with NTL and gf2x libraries.
  !
!#ifdef _NTL_
!  interface  get_coeff_interface
!     subroutine get_coeff(nn,mm,rr,ww,avec,nj,id,pp,np) bind (C, name="get_coeff")
!        use, intrinsic :: iso_c_binding, only: c_int
!        integer (c_int), intent(in) :: nn,mm,rr,ww,avec,nj,id
!        integer (c_int), intent(inout) :: pp(0:nn-1),np
!     end subroutine get_coeff
!  end interface get_coeff_interface
!#else
  interface f_get_coeff_interface
    subroutine f_get_coeff(nn,mm,rr,ww,avec,nj,id,pp,np)
      integer, intent(in) :: nn,mm,rr,ww,avec,nj,id
      integer, intent(inout) :: pp(0:nn-1),np
    end subroutine f_get_coeff 
  end interface f_get_coeff_interface
!#endif

  !==================================
  ! state copy: this => that
  !==================================
  call mt_copy(this,that)
  
  write(str,'(I12)')id
  write(str2,'(I12)')jp
  write(*,'("# ID ",I12)')id
  write(*,'("# Jump Ahead by (",A,")*2^",A)')TRIM(ADJUSTL(str)),TRIM(ADJUSTL(str2))
  if ( this%i /= this%nn) then
    write(*,'("Error in jumpahead: input state pointer should point nn.")')
    write(*,'("this%i  = ",I12)')this%i
    write(*,'("this%nn = ",I12)')this%nn
    write(str,*)this%i
    write(str2,*)this%nn-1
    write(*,'(A,"-",A," random numbers are dropped.")') &
 &                         TRIM(ADJUSTL(str)),TRIM(ADJUSTL(str2))
    write(*,'("forced to point nn")')
    this%i = this%nn
  endif

  !==================================
  ! compute jump ahead polynomial
  ! p(x) coefficients 
  !         for this MT parameter
  !==================================
!#ifdef _NTL_
!  call get_coeff(this%nn,this%mm,this%rr,this%ww,this%aaa,jp,id,p,np)
!#else
  call f_get_coeff(this%nn,this%mm,this%rr,this%ww,this%aaa,jp,id,p,np)
!#endif

  !==================================
  ! multiply p(B) on a state vector v
  !  p(x) : jump ahead polynomial
  !    B  : transition matrix
  ! with simple Horner's method
  !       w = p(B) v = B^(2^jp) v
  !==================================
  v(:) = this%state(:)
  iwp = (np-1)/32
  ibp = mod(np-1,32)
  if (BTEST(p(iwp),ibp)) then
    w(:) = v(:)
  endif
  do i=np-2,0,-1
    iwp = i/32
    ibp = mod(i,32)
    call mt_matvec(this,w,s)   ! s = B w
    if (BTEST(p(iwp),ibp)) then
      w(:) = IEOR(v(:),s(:))   ! w = 1 v + s
    else
      w(:) = s(:)              ! w = 0 v + s
    endif
  enddo

  if (.not. associated(that%state)) then
    allocate(that%state(0:that%nn-1))
  endif
  that%state(:) = w(:)
  that%i = this%nn

  return
end subroutine

function mt_genrand_double1(this) result(r)
  !
  !  r in [0,1]  (53-bit resolution)
  !
  implicit none
  type(mt_state), intent(inout) :: this
  real(REAL64)   :: r
  real(REAL64)   :: a,b
  integer(INT32) :: ia,ib
  ia = mt_genrand_int32(this)   ! ia in [0,0xFFFFFFFF]
  ib = mt_genrand_int32(this)   ! ib in [0,0xFFFFFFFF]
  ia = ISHFT(ia,-5)             ! ia in [0,2^27-1]
  ib = ISHFT(ib,-6)             ! ib in [0,2^26-1]
  a = REAL(ia,kind=KIND(r))
  b = REAL(ib,kind=KIND(r))
  !===============================
  ! ( a*2^26 + b ) in [0,2^53-1]
  ! r = ( a*2^26 + b )/(2^53-1)
  !===============================
  r = (a*67108864.0_REAL64 + b)*(1.0_REAL64/9007199254740991.0_REAL64)
  return
end function

function mt_genrand_double2(this) result(r)
  !
  !  r in [0,1)  (53-bit resolution)
  !
  implicit none
  type(mt_state), intent(inout) :: this
  real(REAL64)   :: r
  real(REAL64)   :: a,b
  integer(INT32) :: ia,ib
  ia = mt_genrand_int32(this)   ! ia in [0,0xFFFFFFFF]
  ib = mt_genrand_int32(this)   ! ib in [0,0xFFFFFFFF]
  ia = ISHFT(ia,-5)             ! ia in [0,2^27-1]
  ib = ISHFT(ib,-6)             ! ib in [0,2^26-1]
  a = REAL(ia,kind=KIND(r))
  b = REAL(ib,kind=KIND(r))
  !===============================
  ! ( a*2^26 + b ) in [0,2^53-1]
  ! r = ( a*2^26 + b )/(2^53)
  !===============================
  r = (a*67108864.0_REAL64 + b)*(1.0_REAL64/9007199254740992.0_REAL64)
  return
end function

function mt_genrand_double3(this) result(r)
  !
  !  r in (0,1)  (52-bit resolution)
  !
  implicit none
  type(mt_state), intent(inout) :: this
  real(REAL64)   :: r
  real(REAL64)   :: a,b
  integer(INT32) :: ia,ib
  ia = mt_genrand_int32(this)   ! ia in [0,0xFFFFFFFF]
  ib = mt_genrand_int32(this)   ! ib in [0,0xFFFFFFFF]
  ia = ISHFT(ia,-6)             ! ia in [0,2^26-1]
  ib = ISHFT(ib,-6)             ! ib in [0,2^26-1]
  a = REAL(ia,kind=KIND(r))
  b = REAL(ib,kind=KIND(r))
  !===============================
  ! ( a*2^26 + b ) in [0,2^52-1]
  ! r = ( a*2^26 + b + 1/2 )/(2^52)
  !===============================
  r = (a*67108864.0_REAL64 + b + 0.5_REAL64)*(1.0_REAL64/4503599627370496.0_REAL64)
  return
end function

function mt_genrand_double4(this) result(r)
  !
  !  r in (0,1]  (53-bit resolution)
  !
  implicit none
  type(mt_state), intent(inout) :: this
  real(REAL64)   :: r
  r = 1.0_REAL64 - mt_genrand_double2(this)
  return
end function

end module
