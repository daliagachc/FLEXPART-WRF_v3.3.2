
       logical function isnan2(a)
       real :: a
        if (isnan(a)) then
       isnan2 = .true.
       else
       isnan2 = .false.
        end if
        return
       end

