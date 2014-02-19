
MODULE dateutils
  ! 
  ! ===========================
  !      MODULE  DATEUTILS
  ! ===========================
  ! 
  ! This module contains generic procedures to facilitate the computation and
  ! handling of dates. 
  !
  ! Procedures: 
  !  * modulosecond
  !  * modulominute
  !  * modulohour
  !  * modulomonth
  !  * gregoriandate
  !  * days_in_year
  !  * isleapyear
  !  * cumsum
  !
  ! Parameters:
  !  * days_per_month
  !  * days_per_month_leap
  !
  ! :Date: May 2008
  ! :Author: David Huard
  ! :Intitution: McGill University
 
  USE WP, only: RP ! Working real precision


  IMPLICIT NONE



  INTEGER, PARAMETER :: days_per_month(12) =  &
       (/31,28,31,30,31,30,31,31,30,31,30,31/)
  INTEGER, PARAMETER :: days_per_month_leap(12) = & 
       (/31,29,31,30,31,30,31,31,30,31,30,31/) 


  INTERFACE cumsum
     MODULE PROCEDURE cumsum_int, cumsum_real
  END INTERFACE

  PRIVATE :: cumsum_int, cumsum_real

CONTAINS
  
  SUBROUTINE raiseerror(PROCEDURE, message, halt)
    ! Write a message to the standard output and stops execution if halt is True.

    CHARACTER(LEN=*), INTENT(IN) :: PROCEDURE, message
    LOGICAL, optional :: halt
    LOGICAL :: ohalt

    WRITE(*,*) '\nProcedure ', PROCEDURE, ' raised the following error: \n\t',&
         message
   
    ohalt = .True.
    IF (PRESENT(halt)) ohalt = halt
    IF (halt) STOP
  END SUBROUTINE raiseerror



  ! MODULO subroutines
  ! ------------------
  ! Given a number of time elements (seconds, minutes, hours) return the number of equivalent
  ! time units (minutes, hours, days)


  pure subroutine modulomilli(milli, remainder, second)
    ! Given the number of milliseconds, return the equivalent number of seconds
    ! and the remaining milliseconds.
    !
    ! Input
    ! -----
    ! milli : integer
    !   Number of milliseconds.
    ! 
    ! Output
    ! ------
    ! remainder : integer
    !   Number of remaining milliseconds.
    ! seconds : integer
    !   Equivalent number of seconds.

    integer, intent(in) :: milli
    integer, intent(out) :: remainder, second

    remainder = modulo(milli, 1000)
    second = (milli - remainder)/1000
  end subroutine modulomilli


  PURE  SUBROUTINE modulosecond(second, remainder, minute)
    ! Given a number of seconds, return the equivalent number of minutes
    ! and the remaining seconds.
    !
    ! Input
    ! -----
    ! second : integer
    !   Number of seconds.
    !
    ! Output
    ! ------
    ! remainder : real
    !   Number of seconds [0,59].
    ! minute : integer
    !   Equivalent number of minutes.

    integer, INTENT(in) :: second
    INTEGER, INTENT(out) :: remainder, minute

    remainder = MODULO(second, 60)
    minute = (second - remainder)/60
  END SUBROUTINE modulosecond


  PURE SUBROUTINE modulominute(minute, remainder, hour)
    ! Given a number of minutes, return the equivalent number of hours
    ! and the remainder in minutes.
    !
    ! Input
    ! -----
    ! minute : integer
    !   Number of minutes.
    !
    ! Output
    ! ------
    ! remainder : integer
    !   Number of minutes [0,59].
    ! hour : integer
    !   Equivalent number of hours.

    INTEGER, INTENT(in) :: minute
    INTEGER, INTENT(out) :: remainder, hour

    remainder = MODULO(minute, 60)
    hour = (minute-remainder)/60
  END SUBROUTINE modulominute


  PURE SUBROUTINE modulohour(hour, remainder, day)
    ! Given a number of hours, return the equivalent number of days
    ! and the remainder in hours.
    !
    ! Input
    ! -----
    ! hour : integer
    !   Number of hours.
    !
    ! Output
    ! ------
    ! remainder : integer
    !   Number of hours [0,23].
    ! day : integer
    !   Equivalent number of days.

    INTEGER, INTENT(in) :: hour
    INTEGER, INTENT(out) :: remainder, day

    remainder = MODULO(hour, 24)
    day = (hour-remainder)/24
  END SUBROUTINE modulohour


  PURE SUBROUTINE modulomonth(month, remainder, year)
    ! Given the number of months, return the equivalent number of 
    ! years and the remainder in months. 
    !
    ! Caution: This is not the same as the month index [1,12].
    !
    !
    ! Input
    ! -----
    ! month : integer
    !   Number of months.
    !
    ! Output
    ! ------
    ! remainder : integer
    !   Number of month [0,11]
    ! year : integer
    !   Equivalent number of years
    INTEGER, INTENT(in) :: month
    INTEGER, INTENT(out) :: remainder, year

    remainder = MODULO(month, 12)
    year = (month-remainder)/12
  END SUBROUTINE modulomonth



  PURE SUBROUTINE gregoriandate(jday, year, month, day) 
    ! Return the month index [1,12] and the day of the month given the year and 
    ! the julian day using the Gregorian calendar.
    !
    ! Input
    ! -----
    ! jday : integer
    ! The julian day [1,366] for leap years and [1,365] for non leap year.
    ! year : integer
    !   The year of the corresponding julian day.
    !
    ! Output
    ! ------
    ! month : integer
    !   The month index [1,12].
    ! day : integer
    !   The day of the month.

    INTEGER, INTENT(in) :: year, jday
    INTEGER, INTENT(out) :: month, day
    INTEGER :: ndays(0:12)
    LOGICAL :: months(13)
    ndays(0) = 0
    IF (isleapyear(year)) THEN
       ndays(1:12) = cumsum_int(days_per_month_leap)
       IF (jday > 366) RETURN
    ELSE
       ndays(1:12) = cumsum_int(days_per_month)
       IF (jday > 365) RETURN
    ENDIF
    months = (jday > ndays)
    month = COUNT(months)
    day = jday - ndays(month-1) 

  END SUBROUTINE gregoriandate



  PURE FUNCTION cumsum_real(x) RESULT (cx)
    ! Return the cumulative sum of x.

    REAL, DIMENSION(:), INTENT(in) :: x
    REAL, DIMENSION(SIZE(x)) :: cx
    INTEGER :: i, n

    n = SIZE(x)
    IF (n==0) RETURN

    cx(1) = x(1)
    DO i=2,n
       cx(i) = cx(i-1) + x(i)
    ENDDO
  END FUNCTION cumsum_real

  PURE FUNCTION cumsum_int(x) RESULT (cx)
    ! Return the cumulative sum of x.

    INTEGER, DIMENSION(:), INTENT(in) :: x
    INTEGER, DIMENSION(SIZE(x)) :: cx
    INTEGER :: i, n

    n = SIZE(x)
    IF (n==0) RETURN

    cx(1) = x(1)
    DO i=2,n
       cx(i) = cx(i-1) + x(i)
    ENDDO
  END FUNCTION cumsum_int


  PURE FUNCTION days_in_year(year)
    ! Return the number of days in the given year: 365 for a normal year and 366 for 
    ! a leap year.

    INTEGER, INTENT(in) :: year
    INTEGER :: days_in_year

    IF (isleapyear(year)) THEN
       days_in_year = 366
    ELSE
       days_in_year = 365
    ENDIF
  END FUNCTION days_in_year


  PURE FUNCTION days_in_month(year, month)
    ! Return the number of days in the given month for the given year.
    !
    ! INPUT
    ! -----
    ! year : integer
    !   The year (YYYY).
    ! month : integer [1,12]
    !   The month index.

    INTEGER, INTENT(in) :: year, month
    INTEGER :: days_in_month

    IF (isleapyear(year)) THEN
       days_in_month = days_per_month_leap(month)
    ELSE
       days_in_month = days_per_month(month)
    END IF
       
  END FUNCTION days_in_month


  PURE FUNCTION isleapyear(year)
    ! Return .True. if year is a leap year, ie, if it has 366 days.

    INTEGER, INTENT(in) :: year
    LOGICAL :: isleapyear

    isleapyear = .FALSE.
    IF (MOD(year,4)==0) THEN
       IF  ((MOD(year, 100) /= 0) .OR. (MOD(year, 400) == 0)) THEN
          isleapyear = .TRUE.
       ENDIF
    ENDIF
  END FUNCTION isleapyear


END MODULE dateutils



