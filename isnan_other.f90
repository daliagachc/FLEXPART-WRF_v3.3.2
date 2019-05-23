       logical function isnan2(a)
! use ieee_arithmetic
       real :: a
!       if (ieee_is_nan(a)) then
        if (a.ne.a) then
       isnan2 = .true.
       else
       isnan2 = .false.
        end if
        return
       end


