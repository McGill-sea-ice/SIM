    if ( var1 == var2 ) then
       call success_assert_action
    else
       call failed_assert_action (to_s(var1), to_s(var2), message)
    end if
