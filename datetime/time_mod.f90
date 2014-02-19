MODULE TIME

  ! The module time provides a time type with operators and related procedures.
  !
  ! Operators: =, ==, >, <.
  ! 
  ! Procedures: 
  !   * time_init: Initialize a time instance.
  !   * seconds_time: Return the equivalent number of seconds.
  !   * days_time: Return the equivalent number of days. 
  !
  ! :Date: May 2008
  ! :Author: David Huard
  ! :Intitution: McGill University
  !

  USE WP, only: RP    ! Working precision.
  USE dateutils

  IMPLICIT NONE
  
  TYPE time_type
     INTEGER :: hour=0, minute=0, second=0, milli=0
  END TYPE time_type


  ! Time operators
  INTERFACE OPERATOR(==); MODULE PROCEDURE time_eq; END INTERFACE;
  INTERFACE OPERATOR(>); MODULE PROCEDURE time_gt; END INTERFACE;
  INTERFACE OPERATOR(<); MODULE PROCEDURE time_lt; END INTERFACE;
  INTERFACE OPERATOR(>=); MODULE PROCEDURE time_geq; END INTERFACE;
  INTERFACE OPERATOR(<=); MODULE PROCEDURE time_leq; END INTERFACE;
  INTERFACE OPERATOR(-); MODULE PROCEDURE time_minus; END INTERFACE;

 
  PRIVATE :: time_eq, time_gt, time_lt, time_geq, time_leq, time_minus


CONTAINS

  ! =================== TIME PROCEDURES ======================


  FUNCTION time_init(hour, minute, second, milli) RESULT (t)
    ! Return a time_type instance given hour, minute and second.
    INTEGER, OPTIONAL, INTENT(in) :: hour, minute, second, milli
    TYPE(time_type) :: t

    IF (PRESENT(hour)) THEN
       IF ((hour < 0) .OR. (hour > 23)) THEN
          CALL raiseerror('time_init', 'invalid hour')
       ELSE
          t%hour = hour
       ENDIF
    ENDIF

    IF (PRESENT(minute)) THEN
       IF ((minute < 0)  .OR. (minute > 59)) THEN
          CALL raiseerror('time_init', 'invalid minute')
       ELSE
          t%minute = minute
       END IF
    ENDIF

    IF (PRESENT(second)) THEN
       IF ((second < 0.0) .OR. (second >=60.)) THEN
          CALL raiseerror('time_init', 'invalid second')
       ELSE
          t%second = second
       END IF
    ENDIF

    IF (present(milli)) then
       if ((milli < 0) .or. (milli >= 1000)) then
          call raiseerror('time_init', 'invalid milli')
       else
          t%milli = milli
       end if
    end IF

  END FUNCTION time_init


  FUNCTION seconds_time(time)
    ! Return the equivalent number of seconds of time instance.

    TYPE(time_type), INTENT(in) :: time
    REAL(RP) :: seconds_time

    seconds_time = time%hour*3600 + time%minute*60 + time%second + time%milli/1000._RP
  END FUNCTION seconds_time


  FUNCTION minutes_time(time)
    ! Return the equivalent number of minutes of time instance.

    TYPE(time_type), INTENT(in) :: time
    REAL(RP) :: minutes_time

    minutes_time = time%hour*60 + time%minute + time%second/60._RP + time%milli/60000._RP
  END FUNCTION minutes_time


  FUNCTION hours_time(time)
    ! Return the equivalent number of hours of time instance.

    TYPE(time_type), INTENT(in) :: time
    REAL(RP) :: hours_time

    hours_time = time%hour + time%minute/60._RP + time%second/3600._RP + time%milli/3600000._RP
  END FUNCTION hours_time
  

  
  FUNCTION days_time(time)
    ! Return the equivalent number of days of time instance.

    TYPE(time_type), INTENT(in) :: time
    REAL(RP) :: days_time

    days_time = time%hour/24._RP + time%minute/(24 * 60._RP) + &
         time%second/(24 * 3600._RP) + &
         time%milli/(24*3600*1000._RP)
  END FUNCTION days_time



  ! =================== TIME OPERATORS =======================


  FUNCTION time_eq(t1, t2)
    ! Return true if t1 == t2 with millisecond accuracy.
    
    TYPE(time_type), INTENT(in) :: t1, t2
    LOGICAL :: time_eq
   
    time_eq = abs(seconds_time(t1) - seconds_time(t2)) < .0001
    

  END FUNCTION time_eq



  FUNCTION time_gt(t1, t2)
    ! Return true if t1 > t2.

    TYPE(time_type), INTENT(in) :: t1, t2
    LOGICAL :: time_gt


    time_gt = (seconds_time(t1) >  seconds_time(t2))
      
  END FUNCTION time_gt


  FUNCTION time_lt(t1, t2)
    ! Return true if t1 < t2.

    TYPE(time_type), INTENT(in) :: t1, t2
    LOGICAL :: time_lt

    time_lt = (seconds_time(t1) < seconds_time(t2))
      
  END FUNCTION time_lt
    
  FUNCTION time_geq(t1,t2)
    ! Return true if t1 >= t2. 
    !
    ! The times are considered equal if within a millionth of a second.


    TYPE(time_type), INTENT(in) :: t1, t2
    LOGICAL :: time_geq

    time_geq = (seconds_time(t1) >= seconds_time(t2))

  END FUNCTION time_geq

    
  FUNCTION time_leq(t1,t2)
    ! Return true if t1 <= t2. 
    !
    ! The times are considered equal if within a millionth of a second.


    TYPE(time_type), INTENT(in) :: t1, t2
    LOGICAL :: time_leq

    time_leq = (seconds_time(t1) <= seconds_time(t2))
    
  END FUNCTION time_leq

  FUNCTION time_minus(t1,t2)
    ! Return t1-t2 in seconds.

    TYPE(time_type), INTENT(in) :: t1, t2
    REAL(RP):: time_minus
    
    time_minus = seconds_time(t1) - seconds_time(t2)
    
  END FUNCTION time_minus



  ! ================== END OF TIME OPERATORS ==================

END MODULE time
