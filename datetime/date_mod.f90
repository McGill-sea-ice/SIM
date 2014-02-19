MODULE DATE
  !
  ! ===================
  !     MODULE DATE
  ! ===================
  !
  ! Provide a derived type `date` based on the Gregorian calendar. 
  !
  ! A date type has three attributes:
  !   * year : integer
  !   * month : integer [1-12]
  !   * day : integer [1-31]
  !
  ! has some overloaded operators: +, -, =, ==, >, >=, <, <=, 
  !
  ! and these additional procedures:
  !   * date_init: Initialize a date instance.
  !   * julianday_date: Return the juliann day of date.
  !   * diff_date, diff_date_delta: Return the difference between dates.
  !   * date_add_days: In-place addition of days.  
  !
  !
  ! Gregorian Calendar
  ! ------------------
  ! 
  ! The Gregorian calendar has 365 days per normal year and 366 per leap year. A year
  ! is said to be a leap year if it is a multiple of 4 and is not a century (eg. 1900),
  ! unless it is a multiple of 400 (eg. 2000):
  !       
  !      MODULO(year, 4) == 0 and ((not MODULO(year, 100)) or MODULO(year, 400))
  !
  ! :Date: May 2008
  ! :Author: David Huard
  ! :Intitution: McGill University




  USE dateutils, ONLY : gregoriandate, cumsum, days_in_year, isleapyear
  USE dateutils, ONLY : days_per_month, days_per_month_leap, raiseerror
  USE datetimedelta

  IMPLICIT NONE
  
  TYPE date_type
     INTEGER :: year=0, month=1, day=1
  END TYPE date_type

  ! Date operators interface

  INTERFACE OPERATOR(+);  MODULE PROCEDURE date_add;  END INTERFACE;
  INTERFACE OPERATOR(-);  MODULE PROCEDURE date_subs;  END INTERFACE;
  INTERFACE OPERATOR(==); MODULE PROCEDURE date_eq; END INTERFACE
  INTERFACE OPERATOR(>); MODULE PROCEDURE date_gt; END INTERFACE
  INTERFACE OPERATOR(<); MODULE PROCEDURE date_lt; END INTERFACE
  INTERFACE OPERATOR(>=); MODULE PROCEDURE date_geq; END INTERFACE
  INTERFACE OPERATOR(<=); MODULE PROCEDURE date_leq; END INTERFACE


  INTERFACE OPERATOR(-);  MODULE PROCEDURE diff_date_delta;  END INTERFACE

   PRIVATE :: date_add, date_subs, date_eq
   PRIVATE :: date_gt, date_lt, date_geq, date_leq
     
CONTAINS


  
  FUNCTION date_init(year, month, day) RESULT (date)
    ! Return a date type given the year, month and day.

    INTEGER, intent(in) :: year
    INTEGER, OPTIONAL, INTENT(in) :: month, day
    TYPE(date_type) :: date

    date%year = year

    IF (PRESENT(month)) THEN
       IF ((month <1) .OR. (month>12)) THEN
          CALL raiseerror('date_init', 'invalid month')
       ELSE
          date%month = month
       END IF
    END IF

    IF (PRESENT(day)) THEN
       IF ((day < 1) .OR. (day > days_in_year(year))) THEN
          CALL raiseerror('date_init', 'invalid day')
       ELSE
          date%day = day
       END IF
    END IF
  END FUNCTION date_init




  SUBROUTINE date_add_days(date, days)
    ! Add days to the given date. 
    !
    ! Input
    ! -----
    ! date : date_type
    !   Current date.
    ! days : integer
    !  Days to add (or substract if negative)

    TYPE(date_type), INTENT(inout) :: date
    INTEGER, INTENT(in) :: days
    INTEGER :: ndays

    ndays = julianday_date(date) + days

    ! Increment or decrement years

    DO WHILE (ndays > days_in_year(date%year)) 
       ndays = ndays - days_in_year(date%year)
       date%year = date%year + 1
    ENDDO

    DO WHILE (ndays <= 0) 
       ndays = ndays + days_in_year(date%year-1)
       date%year = date%year - 1 
    ENDDO


    ! Compute month and day
    CALL gregoriandate(ndays, date%year, date%month, date%day)

  END SUBROUTINE date_add_days


  FUNCTION diff_date_delta(d1, d2) RESULT (delta)
    ! Return the difference between dates d2 and d1 as a delta object.

    TYPE(date_type), INTENT(in) :: d1, d2
    TYPE(datetime_delta_type) :: delta
    INTEGER :: days


    days = diff_date(d1,d2)
    delta = delta_init(days=days)
  END FUNCTION diff_date_delta




  FUNCTION diff_date(d1, d2) RESULT (days)
    ! Return the difference in days between dates d1 and d2 (d1-d2).

    TYPE(date_type), INTENT(in) :: d1, d2
    INTEGER ::days, jd1, jd2, year


    jd1 = julianday_date(d1)
    jd2 = julianday_date(d2)

    days = jd1 - jd2

    IF (d1 > d2)  THEN
       DO year = d2%year, d1%year-1, 1
          days = days + days_in_year(year)
       END DO
    ELSEIF (d1 < d2) THEN
       DO year = d1%year, d2%year-1, 1
          days = days - days_in_year(year)
       ENDDO
    ENDIF

  END FUNCTION diff_date



  PURE FUNCTION julianday_date(date)
    ! Return the julian day of the given date.

    TYPE(date_type), INTENT(in) :: date
    INTEGER :: julianday_date
    INTEGER :: ndays(0:12)

    ndays(0) = 0
    IF (isleapyear(date%year)) THEN
       ndays(1:12) = cumsum(days_per_month_leap)
    ELSE
       ndays(1:12) = cumsum(days_per_month)
    ENDIF

    julianday_date = ndays(date%month-1)+date%day

  END FUNCTION julianday_date



  ! ===================== DATE OPERATORS ===========================

  FUNCTION date_add(date, delta) RESULT (out)
    ! Add a datetime_delta to date. 

    TYPE(date_type), INTENT(in) :: date
    TYPE(datetime_delta_type), INTENT(in) :: delta
    TYPE(date_type) :: out

    out = date

    call date_add_days(out, INT(days_delta(delta)))
  END FUNCTION date_add

  FUNCTION date_subs(date, delta) RESULT (out)
    ! Substract a datetime_delta to date. 

    TYPE(date_type), INTENT(in) :: date
    TYPE(datetime_delta_type), INTENT(in) :: delta
    TYPE(date_type) :: out

    out = date

    call date_add_days(out, -INT(days_delta(delta)))
  END FUNCTION date_subs


  FUNCTION date_eq(d1, d2)
    ! Return True if both dates are equal.

    TYPE(date_type), INTENT(in) :: d1, d2
    LOGICAL :: date_eq

    date_eq = ((d1%year == d2%year) &
         .AND. (d1%month == d2%month) &
         .AND. (d1%day == d2%day))

  END FUNCTION date_eq

  FUNCTION date_gt(d1, d2)
    ! Return True if d1 > d2.

    TYPE(date_type), INTENT(in) :: d1, d2
    LOGICAL :: date_gt

    date_gt = ((d1%year > d2%year) &
         .OR. ((d1%year == d2%year) &
         .AND. (julianday_date(d1) > julianday_date(d2) )))
  END FUNCTION date_gt


  FUNCTION date_lt(d1, d2)
    ! Return True if d1 < d2.

    TYPE(date_type), INTENT(in) :: d1, d2
    LOGICAL :: date_lt

    date_lt = ((d1%year < d2%year) &
         .OR. ((d1%year == d2%year) &
         .AND. (julianday_date(d1) < julianday_date(d2) )))
  END FUNCTION date_lt

  FUNCTION date_geq(d1, d2)
    ! Return True if d1 >= d2.

    TYPE(date_type), INTENT(in) :: d1, d2
    LOGICAL :: date_geq

    date_geq = ((d1%year > d2%year) &
         .OR. ((d1%year == d2%year) &
         .AND. (julianday_date(d1) >= julianday_date(d2) )))
  END FUNCTION date_geq


  FUNCTION date_leq(d1, d2)
    ! Return True if d1 < d2.

    TYPE(date_type), INTENT(in) :: d1, d2
    LOGICAL :: date_leq

    date_leq = ((d1%year < d2%year) &
         .OR. ((d1%year == d2%year) &
         .AND. (julianday_date(d1) <= julianday_date(d2) )))
  END FUNCTION date_leq

  ! ===================== END DATE OPERATORS ======================


END MODULE DATE
