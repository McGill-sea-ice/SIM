    n = size(value)
    l = LBOUND(value,1)
    u = UBOUND(value,1)

    IF (n <= nmax) THEN
      write(form, '(A, I0, A, A, A)') "('[', 1X,",  n,  "(", trim(cid), ",',', 1X), ']')"
      write(out, form) value
    ELSE
       write(form, '(A, I0, A, A, A, I0, A, A, A)')  &
                "('[', 1X,",      &
            d,                &
               "(",               &
           trim(cid),         &
           ",',', 1X),  A, ", &
           d,                 &
           "(",               &
           trim(cid),         &
           ", ',', 1X), ']', )"

       WRITE (out, form) value(l:l+d-1), "..., ",  value(u-d+1:u)
    END IF

    out = ADJUSTL(TRIM(out))
