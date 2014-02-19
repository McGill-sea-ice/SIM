MODULE test_time

  USE time
  USE fruit_ext

CONTAINS



  SUBROUTINE test_time_init()
    TYPE(time_type) :: time

    call set_unit_name('test_time_init')
    time = time_init(6, 30, 5, 200)
    CALL assert_equal(time%hour, 6, message='hour')
    CALL assert_equal(time%minute, 30, message='minute')
    CALL assert_equal(time%second, 5, message='second')
    call assert_equal(time%milli, 200, message='milli')
  END SUBROUTINE test_time_init


  SUBROUTINE test_time_operators()
    TYPE(time_type) :: t1, t2
    REAL(RP) :: ds


    t1 = time_init(5, 30, 20, 500)
    t2 = t1

    call set_unit_name('test_time_operators')

    CALL assert_equal(t2, t1, message='assignment =')
    
    CALL assert_true(t1==t2, message='equivalence ==')

    t2 = time_init(5, 30, 24)
    CALL assert_true(t2>t1, message='>')
    CALL assert_true(t2>=t1, message='>= with >')
    

    t2 = time_init(1)
    CALL assert_true(t2 < t1, message='<')
    CALL assert_true(t2 <= t1, message='<= with <')

    t2 = time_init(5,30,20,500)
    CALL assert_true(t2 >= t1, message='>= with ==' )
    CALL assert_true(t2 <= t1, message='<= with ==' )

    
    t2 = time_init(5,31)
    ds = t2-t1
    CALL assert_equal(ds, 39.5_RP, message='- small')

    t2 = time_init(6,30)
    ds = t2-t1
    CALL assert_equal(ds, 3579.5_RP, message='- large')
    
    

  END SUBROUTINE test_time_operators


  SUBROUTINE test_seconds_time()
    TYPE(time_type) :: time
    
    CALL set_unit_name('test_seconds_time')
    time = time_init(1, 3, 5, 500)
    CALL assert_equal(seconds_time(time), 3600 + 3*60._RP + 5._RP + .5_RP)

  END SUBROUTINE test_seconds_time


  SUBROUTINE test_minutes_time()
    TYPE(time_type) :: time
    
    CALL set_unit_name('test_minutes_time')
    time = time_init(1, 3, 30, 500)
    CALL assert_equal(minutes_time(time), 60._RP + 3 + 1/2._RP + 1/120._RP)

  END SUBROUTINE test_minutes_time


  SUBROUTINE test_hours_time()
    TYPE(time_type) :: time
    
    CALL set_unit_name('test_hours_time')
    time = time_init(1, 3, 5)
    CALL assert_equal(hours_time(time), 1 + 3/60._RP + 5._RP/3600._RP)

  END SUBROUTINE test_hours_time
  

  SUBROUTINE test_days_time()
    TYPE(time_type) :: time

    CALL set_unit_name('test_days_time')
    time = time_init(12, 1, 1)
    CALL assert_equal(days_time(time), .5_RP + 1._RP/24._RP/60._RP + 1._RP/24._RP/3600._RP)

  END SUBROUTINE test_days_time


  SUBROUTINE time_test_suite()  
    CALL test_time_init()
    CALL test_seconds_time()
    CALL test_minutes_time()
    CALL test_hours_time()
    CALL test_days_time()
    CALL test_time_operators()
  END SUBROUTINE time_test_suite


END MODULE test_time
