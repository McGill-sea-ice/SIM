    s1 = shape(var1)
    s2 = shape(var2)

    ! Check that dimensions are identical
    IF (ANY(s1 /= s2)) THEN
       CALL failed_assert_action(to_s(s1), to_s(s2), &
            message//' -- Shape of arrays do not match.')
       RETURN
    END IF
    
    ALLOCATE(mask(s1(1), s1(2)))

    ! Compute the difference between the two arrays.
    mask = (var1 /= var2)
    notequals = COUNT(mask)

    IF (notequals>0) THEN        
       CALL failed_assert_action(to_s(var1), to_s(var2), message // ' -- ' // trim(to_s(notequals))// ' elements are not equal.' )
    ELSE
       call success_assert_action
    ENDIF
