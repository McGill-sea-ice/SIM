
MODULE test_datetime

  ! This module tests procedures that cannot be called from python due to 
  ! features not yet implemented in f2py.
  ! Note though that we can pass python functions as external objects.
  ! For example, in test_cumsum, I pass the assert_array_equal function from 
  ! numpy.testing, which checks that both arrays are identical. If they are not, 
  ! then we get an error raised and we can see the problem right away. I think this
  ! approach is a lot friendlier than using fruit.


  USE datetime
  USE fruit_ext, ONLY: set_unit_name, assert_equal, assert_true, assert_almost_equal
  
  IMPLICIT NONE

CONTAINS


  SUBROUTINE test_cumsum()
    INTEGER :: i(3), ci(3)
    REAL :: x(3), cx(3)

    call set_unit_name('test_cumsum')

    i = (/1,2,3/)
    ci = (/1,3,6/)

    x = (/1.,2.,3./)
    cx = (/1.,3.,6./)

    i = cumsum(i)
    CALL assert_equal(i, ci, 'integers')

    x = cumsum(x)
    CALL assert_equal(x, cx, 'reals')

  END SUBROUTINE test_cumsum


  SUBROUTINE test_datetime_increment() 

    TYPE(datetime_type) :: dt, dt1, dt2

    call set_unit_name('test_datetime_increment')

    dt  = datetime_init(2001, 2, 28, 12, 30, 5)

    dt2 =  datetime_increment(dt, years=10)
    dt%year = 2011
    CALL assert_equal(dt2, dt, 'add year')

    dt  = datetime_init(2001, 2, 28, 12, 30, 5)
    dt2 = datetime_increment(dt, months=13)
    dt%year = 2002
    dt%month = 3
    CALL assert_equal(dt2, dt, 'add months')


    dt  = datetime_init(2001, 2, 28, 12, 30, 5)
    dt2 = datetime_increment(dt, days=56)
    dt%month = 4
    dt%day = 25
    CALL assert_equal(dt2, dt, 'add days')


    dt  = datetime_init(2001, 2, 28, 12, 30, 5)
    dt2 = datetime_increment(dt, hours=49)
    dt%month = 3
    dt%day = 2
    dt%hour = 13
    CALL assert_equal(dt2, dt, 'add hours')


    dt  = datetime_init(2001, 2, 28, 12, 30, 5)
    dt2 = datetime_increment(dt, minutes=24*60)
    dt%month = 3
    dt%day = 1
    dt%hour = 12
    CALL assert_equal(dt2, dt, 'add minutes')

    
    dt  = datetime_init(2001, 2, 28, 12, 30, 5)
    dt2 = datetime_increment(dt, seconds=3600*24)
    dt%month = 3
    dt%day = 1
    dt%hour = 12
    CALL assert_equal(dt2, dt, 'add lots of seconds')

    
    dt  = datetime_init(2001, 2, 28, 12, 30, 5, 700)
    dt2 = datetime_increment(dt, millis=300)
    dt%second = 6
    dt%milli = 0
    CALL assert_equal(dt2, dt, 'add milliseconds')

    


    ! Check month decrementation
    dt = datetime_init(1991, 1)
    dt2 = datetime_increment(dt, months=-13)
    dt1 = datetime_init(1989, 12)
    CALL assert_equal(dt2, dt1, 'decrement months')



  END SUBROUTINE test_datetime_increment


  SUBROUTINE test_datetime_operators()
    TYPE(datetime_type) :: dt1, dt2
    type(datetime_delta_type) :: delta

    CALL set_unit_name('test_datetime_operators')


    ! Check assignment operator (=)  and equivalence (==)
    dt1 = datetime_init(2000, 1, 2, 0)
    dt2 = dt1
    CALL assert_true(dt2 == dt1, 'equivalence')


    ! Check addition of a delta
    dt1 = datetime_init(1979, 6, 21, 11)
    dt2 = datetime_init(1979, 7, 22, 12, minute=1, second=8)
    delta = delta_init(days=31, hours=1, seconds=65, millis=3000)
    CALL assert_equal(dt1+delta, dt2, 'addition')


    ! Check substraction of a day
    dt1 = datetime_init(1979, 6, 21, 11)
    dt2 = datetime_init(1979, 6, 20, 11)
    delta = delta_init(days=1)
    CALL assert_equal(dt1-delta, dt2, 'substraction of a day')

    ! Check substraction of small amounts of time
    dt1 = datetime_init(1979, 6, 21, 11, 0, 0)
    delta = delta_init(days=2, hours=49, minutes=61, seconds=90, millis=1001)
    dt2 = datetime_init(1979, 6, 17, 8, 57,28, 999)
    CALL assert_equal(dt1-delta, dt2, 'substraction of multiple elements')


    ! Check substraction of 6 hours
    dt1 = datetime_init(1979,1,1,0)
    delta = delta_init(hours=6)
    dt2 = dt1 - delta
    call assert_equal(dt2%day, 31, 'substract 6 hours')

    ! Check with year reversals through days
    dt1 = datetime_init(1989)
    delta = delta_init(days = 4)
    dt2 = datetime_init(1988, 12, 28)
    CALL assert_equal(dt1-delta, dt2, 'year reversal')

    ! Idem but with a leap year. 
    dt1 = datetime_init(1991)
    delta = delta_init(days = 4)
    dt2 = datetime_init(1990, 12, 28)
    CALL assert_equal(dt1-delta, dt2, 'leap year reversal')
  
    ! Check with a month decrement in days
    dt1 = datetime_init(1991)
    delta = delta_init(days=92 )
    dt2 = datetime_init(1990, 10 , 1)
    CALL assert_equal(dt1-delta, dt2, 'month decrement in days')

    ! Check comparison operators
    dt1 = datetime_init(1991)
    dt2 = datetime_init(1992)
    CALL assert_true(dt2 > dt1, 'testing >')
    CALL assert_true(dt1 < dt2, 'testing <')
    CALL assert_true(dt2 >= dt1, 'testing >= with >')
    CALL assert_true(dt1 <= dt2, 'testing <= with <') 
    
!    dt2 = datetime_init(1991, second=1.0E-6_RP)
!    call assert_true(dt1 <= dt2, 'testing <= with ~=') 
!    call assert_true(dt2 >= dt1, 'testing >= with ~=') 

    delta = datetime_minus(dt2,dt1)
    CALL assert_equal(delta, delta_init(days=365), 'minus')

    dt2 = datetime_init(1991, 2, 1, 1)
    CALL assert_equal(dt2-dt1, delta_init(days=31, seconds=3600), '-')
    

    

  END SUBROUTINE test_datetime_operators


  SUBROUTINE test_datetime2num()
    TYPE(datetime_type) :: dt, since

    call set_unit_name('test_datetime2num')

    since = datetime_init(1979, 1, 1, 0)

    dt = datetime_init(1979, 1, 4, 6)
    CALL assert_equal(datetime2num(dt, 'days', since), 3.25_RP, 'days')
    
    dt = datetime_init(1996,1,1,0)
    CALL assert_equal(datetime2num(dt, 'days', since), 6209.0_RP, 'lots of days')

    dt = datetime_init(1996,1,1,6)
    CALL assert_equal(datetime2num(dt, 'hours', since), 149022.0_RP, 'lots of hours')

    dt = datetime_init(1996,1,1,6,1)
    CALL assert_equal(datetime2num(dt, 'minutes', since), 8941321._RP, 'lots of minutes')

    dt = datetime_init(1996,1,1,6,1,1)
    CALL assert_equal(datetime2num(dt, 'seconds', since), 536479261._RP, 'lots of seconds')


    since = datetime_init(1, 1, 1, 0, 0)
    dt = datetime_init(1948, 1, 1, 0, 0)
    CALL assert_equal(datetime2num(dt, 'hours', since), 17067072._RP, '1948')
    
    dt = datetime_init(1900, 1, 1, 0, 0)
    CALL assert_equal(datetime2num(dt, 'hours', since), 16646328._RP, '1900')


  END SUBROUTINE test_datetime2num


  SUBROUTINE test_str2dt()
    type(datetime_type) :: dt
    CHARACTER(40) :: str
    CALL set_unit_name('test_str2dt')

    str='    1987-04-23:11:13:45   '
    dt = str2dt(str)
    CALL assert_equal(dt, datetime_init(1987, 4, 23, 11, 13, 45))

  END SUBROUTINE test_str2dt


  SUBROUTINE test_julianday()
    TYPE(date_type) :: d
    type(datetime_type) :: dt

    call set_unit_name('test_julian_day')

    d%year = 2000
    d%month = 3
    d%day = 23
    CALL assert_equal(julianday(d), 83, 'test date')

    
    dt = datetime_init(2000, 3, 23, 12)
    CALL assert_almost_equal(julianday(dt), 83.5_RP, -4,  'test datetime')

  END SUBROUTINE test_julianday


  SUBROUTINE test_date_set_from_datetime()
    TYPE(date_type) :: date
    TYPE(datetime_type) :: dt
    
    call set_unit_name('test_date_set_from_datetime')
    

    dt = datetime_init(1989, 2, 3, 5, 30, 3)
    date = date_set_from_datetime(dt)
    
    CALL assert_equal(date%year, dt%year, 'checking year')
    CALL assert_equal(date%month, dt%month, 'checking month')
    CALL assert_equal(date%day, dt%day, 'checking day')
  END SUBROUTINE test_date_set_from_datetime

  SUBROUTINE test_datetime_set_from_date()
    TYPE(date_type) :: date
    type(datetime_type) :: dt

    CALL set_unit_name('test_datetime_set_from_date')
    date = date_init(1989, 12, 31)
    dt =  datetime_set_from_date(date)

    CALL assert_equal(dt%year, date%year, 'checking year')
    CALL assert_equal(dt%month, date%month, 'checking month')
    CALL assert_equal(dt%day, date%day, 'checking day')
    CALL assert_equal(dt%hour, 0, 'checking hour')
  END SUBROUTINE test_datetime_set_from_date


  SUBROUTINE test_datetime_set_from_date_and_time()
    TYPE(time_type) :: time
    TYPE(date_type) :: date
    TYPE(datetime_type) :: dt

    CALL set_unit_name('test_datetime_set_from_date_and_time')

    time = time_init(4,5,6,0)
    date = date_init(1978, 12, 3)
    
    dt = datetime_set_from_date_and_time(date, time)

    CALL assert_equal(dt%year, 1978)
    CALL assert_equal(dt%month, 12)
    CALL assert_equal(dt%day, 3)
    CALL assert_equal(dt%hour, 4)
    CALL assert_equal(dt%minute, 5)
    CALL assert_equal(dt%second, 6)


  END SUBROUTINE test_datetime_set_from_date_and_time

  subroutine test_datetime_add_days()
    type(datetime_type) :: dt

    call set_unit_name('test_datetime_add_days')
    

    dt = datetime_init(1979,1,1,0)
    call datetime_add_days(dt, -1)
    call assert_equal(dt%day, 31)
  end subroutine test_datetime_add_days


  SUBROUTINE test_datetime_str_6()
    TYPE(datetime_type) :: dt
    CHARACTER(len=6) :: datestr
    
    CALL set_unit_name('test_datetime_str_6')
    dt = datetime_init(1996, 1, 1)
    datestr = '010196'
    CALL assert_equal(datetime_str_6(dt), datestr)
  END SUBROUTINE test_datetime_str_6



  SUBROUTINE test_datetime_str()
    TYPE(datetime_type) :: dt
    CHARACTER(len=19) :: datestr
    
    CALL set_unit_name('test_datetime_str')
    dt = datetime_init(1996, 1, 1)
    datestr = '1996-01-01 00:00:00'
    CALL assert_equal(datetime_str(dt), datestr)
  END SUBROUTINE test_datetime_str

  SUBROUTINE test_now()
    TYPE(datetime_type) :: dt1, dt2
    CALL set_unit_name('test_now')

    dt1 = now()
    CALL assert_true(dt1%year >= 2008)
    dt2 = now()
    CALL assert_equal(dt1, dt2)
  END SUBROUTINE test_now

  SUBROUTINE test_stepper()
    TYPE(datetime_type) :: dt, enddt
    TYPE(datetime_delta_type) :: delta
    INTEGER :: i
    
    CALL set_unit_name('test_stepper')

    ! This is an example of a stepping algorithm.
    
    dt = datetime_init(1900, 1, 1, 0)
    delta = delta_init(seconds=60*30)
    

    ! A do loop
    DO i=1,1234567   
       dt = dt+delta
    END DO

    ! Expected result
    enddt = datetime_init(1970,6,3,3,30, 0)
        
    CALL assert_equal(dt, enddt)
    
  END SUBROUTINE test_stepper
  
  SUBROUTINE datetime_test_suite()

    CALL test_cumsum()
    CALL test_datetime_increment()
    CALL test_datetime_operators()
    CALL test_datetime2num()
    CALL test_julianday()
    CALL test_date_set_from_datetime()
    CALL test_datetime_set_from_date()
    CALL test_datetime_set_from_date_and_time()
    call test_datetime_add_days()
    CALL test_stepper()
    CALL test_datetime_str()
    CALL test_datetime_str_6()
    CALL test_str2dt()
    CALL test_now()
  END SUBROUTINE datetime_test_suite


END MODULE test_datetime



