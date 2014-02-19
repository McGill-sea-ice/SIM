MODULE TEST_DATETIMEDELTA
  ! Test functions for the datetime_delta derived type and operators.

  USE DATETIMEDELTA
  USE FRUIT_ext
  
  IMPLICIT NONE
  
CONTAINS


  subroutine test_delta_init()
    type(datetime_delta_type) :: d1

    call set_unit_name('test_delta_init')
    d1 = delta_init(1,1,1,seconds=1, millis=1)
    call assert_equal(d1%seconds, 1, 'check seconds')
    call assert_equal(d1%millis, 1, 'check millis')


  end subroutine test_delta_init


  SUBROUTINE test_delta_operators() 
    TYPE(datetime_delta_type):: d1, d2, d3, res

    CALL set_unit_name('test_operators')

    d1 = delta_init(days = 1)
    d2 = delta_init(hours = 3)
    d3 = delta_init(days=1, hours=3)
    res = d1+d2
    CALL assert_equal(res, d3, 'testing +')

    res = d1-d2
    d3 = delta_init(hours=21)

    CALL assert_equal(res, d3, 'testing -')

    CALL assert_true(res == d3, 'testing ==')

    CALL assert_true(d1 > d2, 'testing >')
    
    CALL assert_true(d2 < d1, 'testing <')

    CALL assert_true(d1 >= d2, 'testing >=')
    
    CALL assert_true(d2 <= d1, 'testing <=')

    d3 = delta_init(0)
    CALL assert_equal(d3-d1, -d1, 'testing unary -')

    d3 = d1
    CALL assert_equal(d1, d3, 'testing =')

    d3 = delta_init(seconds=3) / 2._RP
    CALL assert_equal(d3, delta_init(seconds=1, millis=500), 'division small')

    d3 = d2/2._RP
    CALL assert_equal(d3, delta_init(minutes=90), 'division < day')

    d3 = d1/2._RP
    CALL assert_equal(d3, delta_init(hours=12), 'division')

    d3 = d1*2._RP
    CALL assert_equal(d3, delta_init(days=2), 'multiplication')

  END SUBROUTINE test_delta_operators


  SUBROUTINE test_millis()
    TYPE(datetime_delta_type):: d1, d2, d3, res

    CALL set_unit_name('test_millis')

    d1 = delta_init(hours=1, minutes=31, seconds=1)
    d2 = d1 / 3._RP
    d3 = delta_init(minutes=30, seconds=20, millis=333)
    call assert_equal(d2,d3)


  end SUBROUTINE test_millis



  SUBROUTINE test_seconds_delta()
    TYPE(datetime_delta_type):: d

    CALL set_unit_name('test_seconds_delta')

    d = delta_init(seconds=45)
    CALL assert_equal(seconds_delta(d), 45.0_RP, 'seconds')

    d = delta_init(minutes=120)
    CALL assert_equal(seconds_delta(d), 120*60._RP, 'minutes')

    d = delta_init(hours=2)
    CALL assert_equal(seconds_delta(d), 2.*60.*60._RP, 'hours')

    d = delta_init(days=2)
    CALL assert_equal(seconds_delta(d), 24*3600*2._RP, 'days')

  END SUBROUTINE test_seconds_delta


  SUBROUTINE test_minutes_delta()
    TYPE(datetime_delta_type) :: d

    call set_unit_name('test_minutes_delta')

    d = delta_init(days=1, hours=1, minutes=1, seconds = 1)
    CALL assert_equal(minutes_delta(d), 24*60._RP + 60 + 1 + 1._RP/60)

  END SUBROUTINE test_minutes_delta


  SUBROUTINE test_hours_delta()
    TYPE(datetime_delta_type) :: d

    call set_unit_name('test_hours_delta')

    d = delta_init(days=1, hours=1, minutes=1, seconds = 1)
    CALL assert_equal(hours_delta(d), 24._RP + 1 + 1._RP/60 + 1._RP/3600)

  END SUBROUTINE test_hours_delta


  SUBROUTINE test_days_delta()
    TYPE(datetime_delta_type):: d
    INTEGER :: seconds_in_day

    seconds_in_day = 60*60*24

    CALL set_unit_name('test_days_delta')

    d = delta_init(millis=500)
    call assert_equal(days_delta(d), .5_RP / seconds_in_day, 'millis')

    d = delta_init(seconds=45)
    CALL assert_equal(days_delta(d), 45._RP/seconds_in_day, 'seconds')

    d = delta_init(minutes=120)
    CALL assert_almost_equal(days_delta(d), 120*60._RP/seconds_in_day, 8, 'minutes')

    d = delta_init(hours=2)
    CALL assert_equal(days_delta(d), 2._RP/24., 'hours')

    d = delta_init(days=2)
    CALL assert_equal(days_delta(d), 2._RP, 'days')
  END SUBROUTINE test_days_delta


  SUBROUTINE test_delta_str()
    TYPE(datetime_delta_type) :: d
    CHARACTER(len=20) :: str

    call set_unit_name('test_delta_str')

    d = delta_init(seconds=50)
    str = delta_str(d)
    CALL assert_equal(str, '50.00 seconds')

    d = delta_init(minutes=30)
    str = delta_str(d)
    CALL assert_equal(str, '30.00 minutes')

    d = delta_init(hours=5, minutes=30)
    str = delta_str(d)
    CALL assert_equal(str, '5.50 hours')

    d = delta_init(days=2)
    str = delta_str(d)
    CALL assert_equal(str, '2.00 days')


  END SUBROUTINE test_delta_str


  SUBROUTINE delta_test_suite()
    
    ! datetime_delta tests
    call test_delta_init()
    CALL test_delta_operators()
    CALL test_seconds_delta()
    CALL test_minutes_delta()
    CALL test_hours_delta()
    CALL test_days_delta()
    CALL test_delta_str()
    call test_millis()
  END SUBROUTINE delta_test_suite

  END MODULE TEST_DATETIMEDELTA





