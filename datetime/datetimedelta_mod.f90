MODULE DATETIMEDELTA
  
  ! =============
  ! DATETIMEDELTA
  ! =============
  !
  ! An implementation of a time difference derived type. 
  !
  ! A datetime_delta type has the following attributes:
  !  * days : integer
  !  * hours : integer
  !  * minutes : integer
  !  * seconds : integer
  !  * millis  : milliseconds
  !
  ! A datetime_delta has no month or year attribute since those don't have a fixed length.
  !
  ! The following operators are defined for datetime_delta: >, <, >=, <=, ==, + and -, as
  ! well as these specific functions: days_delta and seconds_delta.
  !
  ! The datetime_delta type is not interesting in its own sake but is rather
  ! designed for operations on datetime types.
  !
  ! Examples
  ! --------
  !  TYPE(datetime_delta_type) :: d1,d2,d3
  !  d1 = delta_init(days=1, hours=6)
  !  d2 = delta_init(minutes=30)
  !  IF (d1 > d2) d3 = d2+d1
  !
  ! :Author: David Huard <david.huard@gmail.com>
  ! :Date: May 15, 2008
  ! :Institution: McGill University
  ! :License: BSD
  

  USE DATEUTILS, ONLY : modulosecond, modulohour, modulominute, modulomilli, RP


  IMPLICIT NONE

  ! =================================================================
  ! This is the fundamental type and its operators defined 
  ! in this module.
  TYPE datetime_delta_type
     INTEGER :: days=0, hours=0, minutes=0, seconds=0, millis=0
  END TYPE datetime_delta_type

  INTERFACE OPERATOR(+); MODULE PROCEDURE delta_add; END INTERFACE
  INTERFACE OPERATOR(-); MODULE PROCEDURE delta_subs; END INTERFACE
  INTERFACE OPERATOR(*); MODULE PROCEDURE delta_mul; END INTERFACE
  INTERFACE OPERATOR(/); MODULE PROCEDURE delta_div; END INTERFACE
  INTERFACE OPERATOR(-); MODULE PROCEDURE delta_minus; END INTERFACE
  INTERFACE OPERATOR(>=); MODULE PROCEDURE delta_geq; END INTERFACE
  INTERFACE OPERATOR(<=); MODULE PROCEDURE delta_leq; END INTERFACE
  INTERFACE OPERATOR(>); MODULE PROCEDURE delta_gt; END INTERFACE
  INTERFACE OPERATOR(<); MODULE PROCEDURE delta_lt; END INTERFACE
  INTERFACE OPERATOR(==); MODULE PROCEDURE delta_eq; END INTERFACE
!  INTERFACE ASSIGNMENT(=); MODULE PROCEDURE delta_assign; END INTERFACE
  ! ==================================================================

  ! Hide the low level details
  PRIVATE :: delta_geq, delta_leq, delta_gt, delta_lt, delta_eq, delta_add, delta_subs
  PRIVATE :: delta_minus


  INTEGER, PARAMETER :: dp_kind = SELECTED_REAL_KIND(16)

CONTAINS

  
  FUNCTION  delta_init(days, hours, minutes, seconds, millis) RESULT (delta)
    ! Return a datetime delta given the number of days, hours, minutes, 
    ! seconds and milliseconds.

    INTEGER, OPTIONAL, INTENT(in) :: days, hours, minutes, seconds, millis
    TYPE(datetime_delta_type) :: delta

    IF (PRESENT(days)) delta%days = days
    IF (PRESENT(hours)) delta%hours = hours
    IF (PRESENT(minutes)) delta%minutes = minutes
    IF (PRESENT(seconds)) delta%seconds = seconds
    IF (PRESENT(millis)) delta%millis = millis

  END FUNCTION delta_init


  FUNCTION delta_add(d1, d2) RESULT (out)
    ! Return a datetime delta resulting from the addition of two datetime deltas. 

    TYPE(datetime_delta_type), INTENT(IN) :: d1, d2
    TYPE(datetime_delta_type) :: out
    INTEGER :: extra_days
    
    ! Add milliseconds
    CALL modulomilli(d1%millis + d2%millis, out%millis, out%seconds)

    ! Add seconds
    CALL modulosecond(d1%seconds+d2%seconds, out%seconds, out%minutes)

    ! Add minutes
    CALL modulominute(out%minutes + d1%minutes + d2%minutes, out%minutes, out%hours)

    ! Add hours
    CALL modulohour(out%hours + d1%hours + d2%hours, out%hours, extra_days)

    ! Add days
    out%days = d1%days + d2%days + extra_days 
  END FUNCTION delta_add


  FUNCTION delta_subs(d1, d2) RESULT (out)
    ! Return a datetime delta resulting from the substraction of d1 by d2.
    
    TYPE(datetime_delta_type), INTENT(IN) :: d1, d2
    TYPE(datetime_delta_type) :: out

    out = delta_add(d1, -d2)
  END FUNCTION delta_subs


  FUNCTION delta_minus(d) RESULT (out)
    ! Return the opposite of d1, ie. all attributes are negative. 

    TYPE(datetime_delta_type), INTENT(IN) :: d
    TYPE(datetime_delta_type) :: out

    out%days = -d%days
    out%hours = -d%hours
    out%minutes = -d%minutes
    out%seconds = -d%seconds
    out%millis = -d%millis
  END FUNCTION delta_minus
    

  FUNCTION delta_mul(d, x) RESULT (out)
    ! Multiply a delta by x.
    
    TYPE(datetime_delta_type), INTENT(IN) :: d
    REAL(RP), INTENT(in) :: x
    TYPE(datetime_delta_type) :: out
    REAL(dp_kind) :: millis, days

    days = days_delta(d) * x
    millis = 86400000. * MOD(days, 1.0)
    out = delta_init(days=FLOOR(days), millis=int(nint(millis)))

  END FUNCTION delta_mul
    
  FUNCTION delta_div(d, x) RESULT (out)
    ! Divide a delta by x.

    TYPE(datetime_delta_type), INTENT(IN) :: d
    REAL(RP), INTENT(in) :: x
    TYPE(datetime_delta_type) :: out

    out = delta_mul(d, 1./x)

  END FUNCTION delta_div

    
  FUNCTION seconds_delta(d)
    ! Return the equivalent number of seconds in datetime delta d. 

    TYPE(datetime_delta_type), INTENT(IN) :: d
    REAL(RP) :: seconds_delta

    seconds_delta = d%millis/1000_RP + d%seconds + d%minutes*60 + d%hours*3600 + d%days * 3600 * 24
  END FUNCTION seconds_delta


  FUNCTION minutes_delta(d)
    ! Return the equivalent number of minutes in datetime delta d. 
    
    TYPE(datetime_delta_type), INTENT(in) :: d
    REAL(RP) :: minutes_delta

    minutes_delta = d%millis/60000_RP + d%seconds/60._RP + d%minutes + d%hours*60 + d%days*24*60

  END FUNCTION minutes_delta


  FUNCTION hours_delta(d)
    ! Return the equivalent number of hours in datetime delta d. 
    
    TYPE(datetime_delta_type), INTENT(in) :: d
    REAL(RP) :: hours_delta

    hours_delta = d%millis/ (3600*1000.0_RP) +  d%seconds/3600._RP + d%minutes/60._RP + d%hours + d%days*24

  END FUNCTION hours_delta


  FUNCTION days_delta(d)
    ! Return the equivalent number of days in datetime delta d. 
    TYPE(datetime_delta_type), INTENT(IN) :: d
    REAL(RP) :: days_delta

    days_delta = d%millis/(1000*3600*24.0_RP) + &
         d%seconds/(3600 * 24.0_RP) + &
         d%minutes/(60 *24._RP) + &
         d%hours/24._RP  + &
         d%days
  END FUNCTION days_delta
    
  ! Comparison operators
  FUNCTION delta_geq(d1, d2)
    ! Test if d1 >= d2
    LOGICAL :: delta_geq
    TYPE(datetime_delta_type), INTENT(IN) :: d1, d2
    delta_geq = (seconds_delta(d1) >= seconds_delta(d2))
  END FUNCTION delta_geq

  FUNCTION delta_leq(d1, d2)
    ! Test if d1 <= d2
    LOGICAL :: delta_leq
    TYPE(datetime_delta_type), INTENT(IN) :: d1, d2
    delta_leq = (seconds_delta(d1) <= seconds_delta(d2))
  END FUNCTION delta_leq
       
  FUNCTION delta_gt(d1, d2)
    ! Test if d1 > d2
    LOGICAL :: delta_gt
    TYPE(datetime_delta_type), INTENT(IN) :: d1, d2
    delta_gt = (seconds_delta(d1) > seconds_delta(d2))
  END FUNCTION delta_gt

  FUNCTION delta_lt(d1, d2)
    ! Test if d1 < d2
    LOGICAL :: delta_lt
    TYPE(datetime_delta_type), INTENT(IN) :: d1, d2
    delta_lt = (seconds_delta(d1) < seconds_delta(d2))
  END FUNCTION delta_lt

  FUNCTION delta_eq(d1, d2)
    ! Test if d1 and d2 are equal up to a precision of one thousandth of a second.
    LOGICAL :: delta_eq
    TYPE(datetime_delta_type), INTENT(IN) :: d1, d2
    delta_eq = abs(seconds_delta(d1) - seconds_delta(d2)) < .001
  END FUNCTION delta_eq

  FUNCTION delta_str(d, units)
    ! Return the string representation of a datetime delta.
    !
    ! Input
    ! -----
    ! d : datetime_delta_type
    !   Input delta
    ! units : {'second', 'minute', 'hour', 'day'}, optional
    !   Units in which to print the result. If not present, 
    !   the longest unit such that the value > 1 will be 
    !   chosen.

    TYPE(datetime_delta_type), INTENT(in) :: d
    CHARACTER(len=*), INTENT(in), OPTIONAL :: units
    CHARACTER(len=20) :: delta_str, fmt
    CHARACTER(len=6) :: o_units
    REAL(RP) :: x
    INTEGER :: i

    IF (PRESENT(units)) THEN 
       o_units = units
    ELSE

       i = INT(seconds_delta(d))

       SELECT CASE (i)

       CASE (:60)
          o_units = 'second'

       CASE (61:3600)
          o_units = 'minute'

       CASE (3601:24*3600)
          o_units = 'hour'

       CASE default
          o_units = 'day'
       END SELECT

    END IF
    
    SELECT CASE (o_units)

    CASE ('day')
       x = days_delta(d)
       fmt = "(f0.2,' days')"
       
    CASE ('hour')
       x = hours_delta(d)
       fmt = "(f0.2,' hours')"
       
    CASE ('minute')
       x = minutes_delta(d)
       fmt = "(f0.2,' minutes')"
          
    CASE ('second')
       x = seconds_delta(d)
       fmt = "(f0.2,' seconds')"
       
    CASE default
       WRITE(*,*) 'unit not recognized: ', o_units
       STOP
    END SELECT

    WRITE(delta_str, fmt=TRIM(fmt)) x
        
  END FUNCTION delta_str
    


END MODULE DATETIMEDELTA


