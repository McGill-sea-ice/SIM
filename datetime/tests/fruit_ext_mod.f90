MODULE fruit_ext
! This is an extension of the FRUIT module supporting datetime objects.

USE datetime, ONLY: datetime_type, datetime_delta_type, date_type, time_type, RP
USE datetime, only: operator(==)
USE fruit

! adding custom procedures to the general interfaces.
INTERFACE assert_equal
   MODULE PROCEDURE assert_equal_datetime_delta_
   MODULE PROCEDURE assert_equal_datetime_
   MODULE PROCEDURE assert_equal_date_
   MODULE PROCEDURE assert_equal_time_
  END INTERFACE

INTERFACE to_s
   MODULE PROCEDURE  to_s_datetime_delta_
   MODULE PROCEDURE  to_s_datetime_
   MODULE PROCEDURE  to_s_date_
   MODULE PROCEDURE  to_s_time_
END INTERFACE


CONTAINS

! ================================================================
!             ASSERTIONS
! ================================================================

  SUBROUTINE assert_equal_datetime_delta_(var1, var2, message)
    TYPE(datetime_delta_type), INTENT(IN) :: var1, var2
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: message

    IF (var1 == var2) THEN
       call success_assert_action
    ELSE
       CALL failed_assert_action(to_s(var1), to_s(var2), message)
    ENDIF

  END SUBROUTINE assert_equal_datetime_delta_

SUBROUTINE assert_equal_datetime_(var1, var2, message)
    TYPE(datetime_type), INTENT(IN) :: var1, var2
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: message

    IF (var1 == var2) THEN
       call success_assert_action
    ELSE
       CALL failed_assert_action(to_s(var1), to_s(var2), message)
    ENDIF

  END SUBROUTINE assert_equal_datetime_


  SUBROUTINE assert_equal_date_(var1, var2, message)
    TYPE(date_type), INTENT(in) :: var1, var2
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: message

    IF (var1 == var2) THEN
       CALL success_assert_action
    ELSE
       CALL failed_assert_action(to_s(var1), to_s(var2), message)
    ENDIF

  END SUBROUTINE assert_equal_date_
    

  SUBROUTINE assert_equal_time_(var1, var2, message)
    TYPE(time_type), INTENT(in) :: var1, var2
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: message

    IF (var1 == var2) THEN
       CALL success_assert_action
    ELSE
       CALL failed_assert_action(to_s(var1), to_s(var2), message)
    ENDIF

  END SUBROUTINE assert_equal_time_
    

! ================================================================
!            TO STRING
! ================================================================


  FUNCTION to_s_datetime_delta_(value)
    TYPE(datetime_delta_type), intent(in) :: value
    CHARACTER(len=500) :: result, to_s_datetime_delta_
    write (result, *) value
    to_s_datetime_delta_ = adjustl(trim(result))
  END FUNCTION to_s_datetime_delta_

  FUNCTION to_s_datetime_(value)
    TYPE(datetime_type), intent(in) :: value
    CHARACTER(len=500) :: result, to_s_datetime_
    write (result, *) value
    to_s_datetime_ = adjustl(trim(result))
  END FUNCTION to_s_datetime_

  FUNCTION to_s_date_(value)
    TYPE(date_type), intent(in) :: value
    CHARACTER(len=500) :: result, to_s_date_
    write (result, *) value
    to_s_date_ = adjustl(trim(result))
  END FUNCTION to_s_date_

  FUNCTION to_s_time_(value)
    TYPE(time_type), intent(in) :: value
    CHARACTER(len=500) :: result, to_s_time_
    write (result, *) value
    to_s_time_ = adjustl(trim(result))
  END FUNCTION to_s_time_


END MODULE fruit_ext
