MODULE test_date

  ! This module tests procedures in the date module.



  USE date
  USE fruit_ext, ONLY: set_unit_name, assert_equal, assert_true

  IMPLICIT NONE

CONTAINS

  SUBROUTINE test_date_init()
    TYPE(date_type) :: date

    call set_unit_name('test_date_init')
    
    date = date_init(1989, 12, 23)
    CALL assert_equal(date%year, 1989, 'year')
    CALL assert_equal(date%month, 12, 'month')
    CALL assert_equal(date%day, 23, 'day')

  END SUBROUTINE test_date_init
    

  SUBROUTINE test_date_operators()
    TYPE(date_type) :: d1, d2, d3

    CALL set_unit_name('test_operators')

    d1 = date_init(1990,6, 1)
    d2 = date_init(1990, 6, 2)
    d3 = d2

    CALL assert_true(d2 > d1, 'testing >')
    CALL assert_true(d1 < d2, 'testing <')
    CALL assert_true(d2 == d3, 'testing ==')
    CALL assert_true(d2 <= d3, 'testing <= with =')
    CALL assert_true(d1 <= d2, 'testing <= with <')
    CALL assert_true(d2 >= d3, 'testing >= with =')
    CALL assert_true(d2 >= d1, 'testing >= with >')


  END SUBROUTINE test_date_operators



  SUBROUTINE test_date_add_days() 
    TYPE(date_type) :: d1, d2

    call set_unit_name('test_date_add_days')

    d1%year = 2000  ! Leap year
    d1%month = 1
    d1%day=1

    d2 = d1
    CALL date_add_days(d2, 4)
    CALL assert_equal(d2, date_init(2000, 1, 5), 'add a few days')
    
    d2 = d1
    CALL date_add_days(d2, 31)
    CALL assert_equal(d2, date_init(2000, 2, 1), 'add one month.')

    d2 = d1
    call date_add_days(d2, 365)
    CALL assert_equal(d2, date_init(2000, 12, 31), 'add 365 days')

    
    d1%year = 2004
    d1%month = 2
    d1%day=28

    CALL date_add_days(d1, 32)
    CALL assert_equal(d1, date_init(2004, 3, 31), 'add 32 days')


    d1%year = 2004
    d1%month = 2
    d1%day=1

    d2 = d1
    CALL date_add_days(d2, 366)
    CALL assert_equal(d2, date_init(2005, 2, 1), 'add 366 days')

    d2 = d1
    CALL date_add_days(d2, 4*366)
    CALL assert_equal(d2, date_init(2008, 2, 4), 'add four years')

    d1 = date_init(1979,1,1)
    call date_add_days(d1, -1)
    call assert_equal(d1%day, 31)
    call assert_equal(d1%year, 1978)

  END SUBROUTINE test_date_add_days


  SUBROUTINE test_diff_date_delta()
    TYPE(date_type) :: d1, d2
    TYPE(datetime_delta_type) :: delta
    
    CALL set_unit_name('test_diff_date_delta')
    
    d1 = date_init(1989, 1, 1)
    d2 = date_init(1989, 12, 31)
    delta = diff_date_delta(d2, d1)
    CALL assert_equal(delta%days, 364)
    CALL assert_equal(delta%hours, 0)

  END SUBROUTINE test_diff_date_delta



  SUBROUTINE test_diff_date()
    TYPE(date_type) :: d1, d2
    INTEGER :: days

    call set_unit_name('test_diff_date')

    d1 = date_init(1989, 1, 1)
    d2 = date_init(1989, 12, 31)
    days = diff_date(d2, d1)
    CALL assert_equal(days, 364, 'same year')

    d2 = date_init(1990, 1, 1)
    days = diff_date(d2, d1)
    CALL assert_equal(days, 365, 'one year later')

    d2 = date_init(1991, 2, 1)
    days = diff_date(d2, d1)
    CALL assert_equal(days, 2*365+31, 'two year later')

    d2 = date_init(1988, 12, 31)
    days = diff_date(d2, d1)
    CALL assert_equal(days, -1, 'one day before')

    d2 = date_init(1987, 1, 1)
    days = diff_date(d2, d1)
    CALL assert_equal(days, -365-366, 'two years before')
    

  END SUBROUTINE test_diff_date




  SUBROUTINE test_julianday_date()
    TYPE(date_type) :: d

    call set_unit_name('test_julian_day_date')

    d%year = 2000
    d%month = 3
    d%day = 23
    CALL assert_equal(julianday_date(d), 83, 'test1')

    d%year = 2004
    d%month = 12
    d%day = 31
    CALL assert_equal(julianday_date(d), 366, 'test2')

    d%year = 2005
    d%month = 12
    d%day = 31
    CALL assert_equal(julianday_date(d), 365, 'test3')

  END SUBROUTINE test_julianday_date
    
  SUBROUTINE date_test_suite()

    ! date tests
    CALL test_date_init()
    CALL test_date_add_days()
    CALL test_diff_date_delta()
    CALL test_diff_date()
    CALL test_julianday_date()
    CALL test_date_operators()
  END SUBROUTINE date_test_suite

    


END MODULE test_date
