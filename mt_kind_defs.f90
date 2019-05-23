module mt_kind_defs
  implicit none
  public
  integer, parameter :: INT32  = selected_int_kind(9)
  integer, parameter :: INT64  = selected_int_kind(18)
  integer, parameter :: REAL64 = selected_real_kind(15)
end module
