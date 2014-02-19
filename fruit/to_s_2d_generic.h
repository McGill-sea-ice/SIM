    n = size(value, 1)
    l = LBOUND(value,1)
    u = UBOUND(value,1)

    k = 1
    IF (n < nmax) THEN
      DO i=l,u
        write(temp(k), '(A)') trim(to_s(value(i,:)))
    k = k + 1
      enddo
      ns = maxval(len_trim(temp(1:k-1)))
      write(form, '(A, I0, A, I0, A)') "(", n, "('\n', 1X, A", ns, "))"
!      write(*,*) form
      write(out, form) temp(1:n)
    ELSE
      do i=l, l+d-1
        write(temp(k), '(A)') trim(to_s(value(i,:)))
    k = k + 1
      enddo
      ns = maxval(len_trim(temp(1:k-1)))
      dots = index(temp(k-1), '...')
      if (dots == 0) dots = ns/2-1
      write(temp(k), '(A, A, A, A, A)') "[", repeat(" ", dots-2), "...",  repeat(" ", len_trim(temp(k-1)) - dots-3), "]"
      k = k+1
      do i=u-d+1, u
        write(temp(k), '(A)') trim(to_s(value(i,:)))
    k = k + 1
      enddo
      ns = maxval(len_trim(temp(1:k-1)))

      write(form, '(A, I0, A, I0, A)')        &
              "(",                        &
          k-1,                        &
          "('\n', 1X, A",             &
          ns,                         &
          "))"
!      write(*,*) form
      write(out, form) temp(1:k-1)
    END IF

    out = ADJUSTL(TRIM(out))
