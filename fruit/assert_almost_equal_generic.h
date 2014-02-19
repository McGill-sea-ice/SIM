! We compare the log in base 10 of the absolute difference to the precision
! required from the user. Note the factor of .01 that is substracted to `odec`
! to make sure that roundoff errors do not corrupt the results.

    IF (PRESENT(dec)) odec = dec
      if ( log10(1.0 * abs(var1 - var2)) >= odec-.01) then
       call failed_assert_action (to_s(var1), to_s(var2), message)
    else
       call success_assert_action
    end if
!  write(*,*) log10(1.0*abs(var1-var2)), odec
