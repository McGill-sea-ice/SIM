

MODULE DATETIME  
  !
  ! =======================  
  !     MODULE DATETIME
  ! =======================
  !
  ! Provides derived types related to date and time keeping using the Gregorian calendar.
  ! 
  ! Derived types:
  !  #. time_type
  !  #. date_type
  !  #. datetime_type
  !
  ! Operators
  !   +, -, >, <, >=, <=
  ! 
  ! :Date: May 2008
  ! :Author: David Huard
  ! :Intitution: McGill University
  !

  USE dateutils
  USE datetimedelta
  USE date
  USE time

  IMPLICIT NONE

  TYPE datetime_type
     INTEGER :: year=1900, month=1, day=1, hour=0, minute=0, second=0, milli=0
  END TYPE datetime_type

  ! Datetime operators
  INTERFACE OPERATOR(+);  MODULE PROCEDURE datetime_add;  END INTERFACE;
  INTERFACE OPERATOR(-);  MODULE PROCEDURE datetime_minus_delta; END INTERFACE;
  INTERFACE OPERATOR(-); MODULE PROCEDURE datetime_minus; END INTERFACE
  INTERFACE OPERATOR(==); MODULE PROCEDURE datetime_eq; END INTERFACE
  INTERFACE OPERATOR(>); MODULE PROCEDURE datetime_gt; END INTERFACE
  INTERFACE OPERATOR(<); MODULE PROCEDURE datetime_lt; END INTERFACE
  INTERFACE OPERATOR(<=); MODULE PROCEDURE datetime_leq; END INTERFACE
  INTERFACE OPERATOR(>=); MODULE PROCEDURE datetime_geq; END INTERFACE



  INTERFACE julianday
     MODULE PROCEDURE julianday_date, julianday_datetime
  END INTERFACE

  INTERFACE seconds
     MODULE PROCEDURE  seconds_time, seconds_delta
  END INTERFACE

  INTERFACE minutes
     MODULE PROCEDURE minutes_time, minutes_delta
  END INTERFACE

  INTERFACE hours
     MODULE PROCEDURE hours_time, hours_delta
  END INTERFACE

  INTERFACE days
     MODULE PROCEDURE  days_time, days_delta
  END INTERFACE


CONTAINS


  ! ======================== DATETIME OPERATORS ==========================

  FUNCTION datetime_init(year, month, day, hour, minute, second, milli) RESULT (dt)
    !Return a datetype type given the year, month, day, hour, minute and second.

    INTEGER, INTENT(in) :: year
    INTEGER, OPTIONAL,  INTENT(IN) :: month, day, hour, minute, second, milli
    TYPE(datetime_type) :: dt
    LOGICAL :: isleap

    isleap = isleapyear(year)
    
    dt%year = year

    ! Check the validity of the inputs and then assign them
    IF (PRESENT(month)) THEN
       IF ((month <1) .OR. (month>12)) THEN
          CALL raiseerror('datetime_init', 'invalid month')
       ELSE
          dt%month = month
       END IF
    END IF

    IF (PRESENT(day)) THEN
       IF (day<1) CALL raiseerror('datetime_init', 'invalid day')
       
       IF (PRESENT(month)) THEN
          IF (day > days_in_month(year, month)) CALL raiseerror('datetime_init', 'invalid day')
       ELSE
          IF (day > days_in_year(year)) CALL raiseerror('datetime_init', 'invalid day')
       END IF

       dt%day = day

    END IF
    
  
    IF (PRESENT(hour)) THEN
       IF ((hour < 0) .OR. (hour > 23)) THEN
          CALL raiseerror('time_init', 'invalid hour')
       ELSE
          dt%hour = hour
       ENDIF
    ENDIF

    IF (PRESENT(minute)) THEN
       IF ((minute < 0)  .OR. (minute > 59)) THEN
          CALL raiseerror('time_init', 'invalid minute')
       ELSE
          dt%minute = minute
       END IF
    ENDIF

    IF (PRESENT(second)) THEN
       IF ((second < 0.0) .OR. (second >=60.)) THEN
          CALL raiseerror('time_init', 'invalid second')
       ELSE
          dt%second = second
       END IF
    ENDIF

    IF (present(milli)) then
       if ((milli < 0) .or. (milli >= 1000)) then
          call raiseerror('time_init', 'invalid milli')
       else
          dt%milli = milli
       end if
    end IF




  END FUNCTION datetime_init

  FUNCTION datetime_increment(dt, years, months, days, hours, minutes, seconds, millis) RESULT (out)
    ! Increment a datetime by the given years, months, days, hours, minutes, seconds.
    !
    ! years and months increments cannot be used with other increments since they are of 
    ! varying lengths and hence are order dependent. Hence, datetime_increment(years=1, days=4)
    ! is forbidden. Note that increments can be negative. 
    !
    !
    ! Input
    ! -----
    ! dt : datetime 
    !   The base date and time.
    ! years : integer
    !   The increment in years. years increments cannot be mixed with any other increment. 
    ! months : integer
    !   The increment in months. months increments cannot be mixed with any other increment.
    ! days : integer
    !   The increment in days.
    ! hours : integer
    !   The increment in hours. 
    ! minutes : integer
    !   The increment in minutes.
    ! seconds : integer
    !   The increment in seconds. 
    ! millis : integer
    !   The increment in milliseconds.
    !
    ! Output
    ! ------
    ! out : datetime
    !   The resulting datetime. 


    TYPE(datetime_type), INTENT(IN) :: dt
    INTEGER, INTENT(IN), OPTIONAL :: years, months, days, hours, minutes, seconds, millis
    TYPE(datetime_type) :: out
    INTEGER :: oyears, omonths, odays, ohours, ominutes, oseconds, omillis
    INTEGER :: rem_seconds, rem_minutes, rem_hours, rem_days, tmpmonth, rem_years

    ! Handling of optional arguments
    oyears = 0; IF (PRESENT(years)) oyears = years
    omonths = 0; IF (PRESENT(months)) omonths = months
    odays = 0; IF (PRESENT(days)) odays = days
    ohours = 0; IF (PRESENT(hours)) ohours = hours
    ominutes = 0; IF (PRESENT(minutes)) ominutes = minutes
    oseconds = 0; IF (PRESENT(seconds)) oseconds = seconds
    omillis = 0; IF (present(millis)) omillis = millis


    ! Handling of illegal operations
    ! Years and months have different lengths depending on the current date, so they 
    ! are not valid deltas unless they are specified alone, or that there is a convention 
    ! on the order of updating. But the latter solution seems to be a recipe for subtle bugs.
    IF ((oyears > 0) .AND. (omonths>0 .OR. odays>0 .OR. ohours>0 & 
         .OR. ominutes>0 .OR. oseconds>0)) & 
         CALL raiseerror('dateutils.datetime_increment', & 
         'if years is > 0, no other deltas can also be > 0.')

    IF ((omonths>0) .AND. (odays>0 .OR. ohours>0 .OR. ominutes>0 .OR. oseconds>0)) &
         CALL raiseerror('dateutils.datetime_increment', &
         'if months is > 0, no other deltas can also be > 0.')

    out = dt

    ! Add milliseconds
    call modulomilli(out%milli + omillis, out%milli, rem_seconds)

    ! Add seconds
    CALL modulosecond(out%second+oseconds+rem_seconds, out%second, rem_minutes)

    ! Add minutes
    CALL modulominute(out%minute + ominutes + rem_minutes, out%minute, rem_hours)

    ! Add hours
    CALL modulohour(out%hour + ohours + rem_hours, out%hour, rem_days)

    ! Add days
    CALL datetime_add_days(out, odays+rem_days)


    ! Add months
    CALL modulomonth(out%month - 1 + omonths, tmpmonth, rem_years)
    out%month = tmpmonth+1

    ! Add years
    out%year = out%year + oyears + rem_years

  END FUNCTION datetime_increment

  FUNCTION datetime_str_6(dt)
    ! Return a string representation (DDMMYY) of the datetime instance.

    TYPE(datetime_type), INTENT(in) :: dt
    CHARACTER(LEN=6) :: datetime_str_6
    CHARACTER(len=4) :: day, month, year

    WRITE(day, '(i0.2)') dt%day
    WRITE(month, '(i0.2)') dt%month
    WRITE(year, '(i4)') dt%year
    datetime_str_6 =  trim(day) // trim(month) // year(3:4)
  END FUNCTION datetime_str_6

  FUNCTION datetime_str(dt)
    ! Return a string representation (YYYY-MM-DD 00:00:00) of the datetime instance.

    TYPE(datetime_type), INTENT(in) :: dt
    CHARACTER(LEN=19) :: datetime_str
    WRITE(datetime_str, "(i0.4, '-', i0.2, '-', i0.2, ' ', i0.2, ':', i0.2, ':', i0.2)" ) &
         dt%year, dt%month, dt%day, dt%hour, dt%minute, dt%second
  END FUNCTION datetime_str


  FUNCTION now()
    ! Return a datetime instance of the current date and time.
    TYPE(datetime_type) :: now
    CHARACTER(len=12) :: ds, ts
    INTEGER :: year, month, day, hour, minute, milli
    REAL(RP) :: seconds
    
    CALL DATE_AND_TIME(ds, ts)
    READ(ds, fmt="(i4,i2,i2)") year, month, day
    READ(ts, fmt="(i2,i2,f5.3)") hour, minute, seconds
    milli = int(mod(seconds, 1.))*1000
    now = datetime_init(year, month, day, hour, minute, floor(seconds), milli)

  END FUNCTION now
    

  FUNCTION datetime2num(dt, unit, since) RESULT (out)
    ! Conversion from datetime to real using the base unit. 
    !
    ! :Input:
    ! 
    ! dt : datetime_type
    !   The date at which to count the number of time units. 
    ! unit : character string
    !   The definition of the counting unit and the starting point. 
    !   Example : 'days', 'hours', 'minutes', 'seconds'
    ! since : datetime_type
    !   The reference time. That is, the result will be the number of units between 
    !   `dt` and `since`. 
    !
    ! :Output:
    ! 
    ! out : real
    !   The number of units between `dt` and `since`.
    !
    ! :Notes:
    !
    ! In the UDUNITS package, the calendar used for dates before 1582-10-15 is the
    ! Julian Calendar. To avoid this hassle, we will raise an error if dt is before
    ! that date, and add an offset from that date. A Julian year has exactly 365.25 days, 
    ! while a Gregorian year has 365.2425 days.  
    
    ! Lilius originally proposed that the 10-day correction should be
    ! implemented by deleting the Julian leap day on each of its ten
    ! occurrences during a period of 40 years, thereby providing for a
    ! gradual return of the equinox to 21 March. However, Clavius's
    ! opinion was that the correction should take place in one move
    ! and it was this advice which prevailed with
    ! Gregory. Accordingly, when the new calendar was put in use, the
    ! error accumulated in the 13 centuries since the Council of
    ! Nicaea was corrected by a deletion of ten days. The last day of
    ! the Julian calendar was Thursday, 4 October 1582 and this was
    ! followed by the first day of the Gregorian calendar, Friday, 15
    ! October 1582 (the cycle of weekdays was not affected).
    

    TYPE(datetime_type), INTENT(in) :: dt, since
    CHARACTER(len=*), INTENT(in) :: unit
    REAL(RP) :: out
    TYPE(datetime_type) :: changeover
    TYPE(datetime_delta_type) :: delta

    changeover = datetime_init(1582, 10, 15)

    IF (dt < changeover) THEN
       WRITE(*,*) 'datetime2num not implemented for dates before 1582-10-15.'
       STOP
    END IF

   delta = dt-since

   IF (since < changeover) THEN
      IF (since == datetime_init(1,1,1,0)) THEN
          delta%days = delta%days + 2
       ELSE
          WRITE(*,*) "Not implemented." 
          STOP
       ENDIF
    ENDIF

    SELECT CASE(unit)
    CASE ('days')
       out = days(delta)
    CASE ('hours')
       out = hours(delta)
    CASE ('minutes')
       out = minutes(delta)
    CASE ('seconds')
       out = seconds(delta)
    END SELECT

  END FUNCTION datetime2num
    

  FUNCTION str2dt(str)
    ! Given a string of the type yyyy-mm-dd:hh:mm:ss, return a datetime_type object.
    CHARACTER(len=*), INTENT(in) :: str
    TYPE(datetime_type) :: str2dt
    character(len=19) :: ds
    integer :: year, month, day, hour, minute, second
    ds = TRIM(ADJUSTL(str))
    READ(unit=ds, fmt="(i4, tr1, i2, tr1, i2, tr1, i2, tr1, i2, tr1, i2)" )  &
         year, month, day, hour, minute, second
    str2dt = datetime_init(year, month, day, hour, minute, second)
  END FUNCTION str2dt


  FUNCTION julianday_datetime(datetime)
    ! Return the julian day of a datetime object. 
    !
    ! :Input:
    ! 
    ! datetime : datetime_type object
    !   The date and time at which the julian day is computed using the gregorian calendar.
    ! 
    ! :Output:
    ! 
    ! julian_day_datetime : real
    !   The julian day at datetime with fraction corresponding to the hours, minutes, seconds.
    
    TYPE(datetime_type), INTENT(in) :: datetime
    REAL(RP) :: julianday_datetime, fraction
    TYPE(date_type) :: date
    type(time_type) :: time
    INTEGER :: days

    date = date_set_from_datetime(datetime)

    time = time_set_from_datetime(datetime)
    days = julianday_date(date)
    fraction = seconds_time(time) / 24. / 3600.
    
    julianday_datetime = days + fraction

  END FUNCTION julianday_datetime


  FUNCTION diff_time(t1, t2) RESULT (delta)
    ! Return t2-t1 as a delta

    TYPE(time_type), INTENT(in) :: t1, t2
    type(datetime_delta_type) :: delta
    

    delta = delta_init( millis=int(seconds_time(t2) - seconds_time(t1)*1000) )
    
  END FUNCTION diff_time


    
  SUBROUTINE datetime_add_days(dt, days)
    ! Add days to datetime instance dt. 
    
    TYPE(datetime_type), intent(inout) :: dt
    INTEGER, INTENT(in) :: days
    TYPE(date_type) :: date
    
    date = date_set_from_datetime(dt)
    CALL date_add_days(date, days)
    
    dt%year = date%year
    dt%month = date%month
    dt%day = date%day
    
  END SUBROUTINE datetime_add_days



  FUNCTION datetime_set_from_date(date) RESULT (dt)
    ! Set a datetime from a date.

    TYPE(date_type), INTENT(in) :: date
    TYPE(datetime_type) :: dt

    dt = datetime_init(date%year, date%month, date%day)
  END FUNCTION datetime_set_from_date


  FUNCTION date_set_from_datetime(dt) RESULT (date)
    ! Set a date from a datetime
    TYPE(datetime_type), INTENT(IN) :: dt
    TYPE(date_type) :: date

    date = date_init(dt%year, dt%month, dt%day)
  END FUNCTION date_set_from_datetime
    

  FUNCTION time_set_from_datetime(dt) RESULT (t)
    ! Return a time_type instance by taking only the hour, minute, seconds portion
    ! of datetime.

    TYPE(datetime_type), INTENT(in) :: dt
    TYPE(time_type) :: t
    
    t = time_init(dt%hour, dt%minute, dt%second, dt%milli)
  END FUNCTION time_set_from_datetime


  FUNCTION datetime_set_from_date_and_time(date, time) RESULT (dt)
    ! Return a datetime instance created from a date and a time instance.

    TYPE(date_type), INTENT(in) :: date
    TYPE(time_type), INTENT(in) :: time
    TYPE(datetime_type) :: dt

    dt = datetime_init(date%year, date%month, date%day, time%hour, time%minute, time%second, time%milli)

  END FUNCTION datetime_set_from_date_and_time



  ! ==================== DATETIME OPERATORS =========================


  FUNCTION datetime_eq(dt1, dt2)
    ! Return True if both datetime objects are equal.

    TYPE(datetime_type), INTENT(IN) :: dt1, dt2
    TYPE(date_type) :: d1, d2
    TYPE(time_type) :: t1, t2
    LOGICAL :: datetime_eq

    d1 = date_set_from_datetime(dt1)
    d2 = date_set_from_datetime(dt2)
    
    t1 = time_set_from_datetime(dt1)
    t2 = time_set_from_datetime(dt2)

    datetime_eq = ((d1 == d2) .AND. (t1==t2))

  END FUNCTION datetime_eq

  FUNCTION datetime_gt(dt1, dt2)
    ! Return true if dt1 > dt2

    TYPE(datetime_type), INTENT(in) :: dt1, dt2
    TYPE(date_type) :: d1, d2
    TYPE(time_type) :: t1, t2
    LOGICAL :: datetime_gt


    d1 = date_set_from_datetime(dt1)
    d2 = date_set_from_datetime(dt2)
    
    t1 = time_set_from_datetime(dt1)
    t2 = time_set_from_datetime(dt2)

    datetime_gt = ((d1 > d2) .OR. ((d1==d2) .AND. (t1 > t2)))

  END FUNCTION datetime_gt
       
      

  FUNCTION datetime_lt(dt1, dt2)
    ! Return true if dt1 < dt2

    TYPE(datetime_type), INTENT(in) :: dt1, dt2
    TYPE(date_type) :: d1, d2
    TYPE(time_type) :: t1, t2
    LOGICAL :: datetime_lt


    d1 = date_set_from_datetime(dt1)
    d2 = date_set_from_datetime(dt2)
    
    t1 = time_set_from_datetime(dt1)
    t2 = time_set_from_datetime(dt2)

    datetime_lt = ((d1 < d2) .OR. ((d1==d2) .AND. (t1<t2)))

  END FUNCTION datetime_lt


  FUNCTION datetime_geq(dt1, dt2)
    ! Return true if dt1 >= dt2.

    TYPE(datetime_type), INTENT(in) :: dt1, dt2
    TYPE(date_type) :: d1, d2
    TYPE(time_type) :: t1, t2
    LOGICAL :: datetime_geq


    d1 = date_set_from_datetime(dt1)
    d2 = date_set_from_datetime(dt2)
    
    t1 = time_set_from_datetime(dt1)
    t2 = time_set_from_datetime(dt2)

    datetime_geq = ((d1 > d2) .OR. ((d1==d2) .AND. (t1>=t2)))

  END FUNCTION datetime_geq


  FUNCTION datetime_leq(dt1, dt2)
    ! Return true if dt1 <= dt2

    TYPE(datetime_type), INTENT(in) :: dt1, dt2
    TYPE(date_type) :: d1, d2
    TYPE(time_type) :: t1, t2
    LOGICAL :: datetime_leq


    d1 = date_set_from_datetime(dt1)
    d2 = date_set_from_datetime(dt2)
    
    t1 = time_set_from_datetime(dt1)
    t2 = time_set_from_datetime(dt2)

    datetime_leq = ((d1 < d2) .OR. ((d1==d2) .AND. (t1<=t2)))
    
  END FUNCTION datetime_leq




  FUNCTION datetime_add(dt, delta) RESULT (out)
    ! Add a datetime_delta to datetime dt.

    TYPE(datetime_type), INTENT(IN) :: dt
    TYPE(datetime_delta_type), INTENT(IN) :: delta
    TYPE(datetime_type) :: out

    out = datetime_increment(dt, days=delta%days, hours=delta%hours, & 
         minutes=delta%minutes, seconds=delta%seconds, millis=delta%millis)

  END FUNCTION datetime_add


  FUNCTION datetime_minus_delta(dt, delta) RESULT (out)
    ! Substract a datetime_delta from datetime dt.

    TYPE(datetime_type), INTENT(IN) :: dt
    TYPE(datetime_delta_type), INTENT(IN) :: delta
    TYPE(datetime_type) :: out

    out = datetime_add(dt, -delta)

  END FUNCTION datetime_minus_delta


  FUNCTION datetime_minus(dt1, dt2)
    ! Return the time difference between dt1 and dt2 (dt1-dt2) as a delta instance.
    ! Note that the returned delta is defined only in terms of days and seconds.

    TYPE(datetime_type), INTENT(IN) :: dt1, dt2
    TYPE(datetime_delta_type) :: datetime_minus    
    TYPE(date_type) :: d1, d2    
    TYPE(time_type) :: t1, t2
    INTEGER :: days, seconds

    
    d1 = date_set_from_datetime(dt1)
    d2 = date_set_from_datetime(dt2)

    days = diff_date(d1,d2)

    t1 = time_set_from_datetime(dt1)
    t2 = time_set_from_datetime(dt2)
    
    seconds = t1 - t2
    
    datetime_minus = delta_init(days=days, seconds=seconds)
   
  END FUNCTION datetime_minus




  ! =================== END OF DATETIME OPERATORS ===========================



END MODULE DATETIME

