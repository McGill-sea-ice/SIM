module fruit_util

  IMPLICIT NONE
  ! ==================== Added by Huard =======================
!  USE datetime, ONLY : datetime_type, date_type, time_type, datetime_delta_type

  ! ===========================================================

  interface equals
     module procedure equalEpsilon
     module procedure floatEqual
     module procedure integerEqual
     module procedure doublePrecisionEqual
     module procedure stringEqual
     module procedure logicalEqual
  end interface

  interface to_s
     module procedure to_s_int_
     module procedure to_s_real_
     module procedure to_s_logical_
     module procedure to_s_double_
     module procedure to_s_complex_
     module procedure to_s_double_complex_
     module procedure to_s_string_
     module procedure to_s_1d_int_
     module procedure to_s_1d_real_
     module procedure to_s_1d_logical_
     module procedure to_s_1d_double_
     module procedure to_s_1d_complex_
     module procedure to_s_2d_int_
     module procedure to_s_2d_real_
     module procedure to_s_2d_logical_
     module procedure to_s_2d_double_
  end interface

  private :: to_s_int_, to_s_real_, to_s_logical_, to_s_double_, to_s_complex_
  private :: to_s_double_complex_, to_s_string_
  private :: to_s_1d_int_, to_s_1d_real_, to_s_1d_double_, to_s_1d_complex_
  PRIVATE :: to_s_2d_int_, to_s_2d_real_

contains

  function to_s_int_ (value)
    implicit none
    character(len=500):: to_s_int_
    integer, intent(in) :: value
    character(len=500) :: result
    WRITE (RESULT, '(I0,1X)') value
    to_s_int_ = adjustl(trim(result))
  end function to_s_int_

  function to_s_real_ (value)
    implicit none
    character(len=500):: to_s_real_
    real, intent(in) :: value
    character(len=500) :: result
    write (result, *) value
    to_s_real_ = adjustl(trim(result))
  end function to_s_real_

  function to_s_double_ (value)
    implicit none
    character(len=500):: to_s_double_
    double precision, intent(in) :: value
    character(len=500) :: result
    write (result, *) value
    to_s_double_ = adjustl(trim(result))
  end function to_s_double_

  function to_s_complex_ (value)
    implicit none
    character(len=500):: to_s_complex_
    complex, intent(in) :: value
    character(len=500) :: result
    write (result, *) value
    to_s_complex_ = adjustl(trim(result))
  end function to_s_complex_

  function to_s_double_complex_ (value)
    implicit none
    character(len=500):: to_s_double_complex_
    double complex, intent(in) :: value
    character(len=500) :: result
    write (result, *) value
    to_s_double_complex_ = adjustl(trim(result))
  end function to_s_double_complex_

  function to_s_logical_ (value)
    implicit none
    character(len=500):: to_s_logical_
    logical, intent(in) :: value
    character(len=500) :: result
    write (result, *) value
    to_s_logical_ = adjustl(trim(result))
  end function to_s_logical_

  function to_s_string_ (value)
    implicit none
    character(len=500):: to_s_string_
    character(len=*), intent(in) :: value
    to_s_string_ = value
  end function to_s_string_

  FUNCTION to_s_1d_int_ (value) RESULT (out)
    INTEGER :: value(:)
    CHARACTER(len=500) :: out
    CHARACTER(len=6) :: cid
    CHARACTER(len=100) :: form
    INTEGER :: l, u, d, n, nmax !,i
    PARAMETER (d=4, nmax=10)

    cid = 'I0'

    INCLUDE 'to_s_1d_generic.h'
  END FUNCTION to_s_1d_int_

  FUNCTION to_s_1d_real_ (value) RESULT (out)
    REAL :: value(:)
    CHARACTER(len=500) :: out
    CHARACTER(len=10) :: cid
    CHARACTER(len=100) :: form
    INTEGER :: l, u, d, n, nmax !,i
    PARAMETER (d=4, nmax=10)

    cid = 'F10.5'

    INCLUDE 'to_s_1d_generic.h'
  END FUNCTION to_s_1d_real_

  FUNCTION to_s_1d_double_ (value) RESULT (out)
    DOUBLEPRECISION :: value(:)
    CHARACTER(len=500) :: out
    CHARACTER(len=10) :: cid
    CHARACTER(len=100) :: form
    INTEGER :: l, u, d, n, nmax !,i
    PARAMETER (d=3, nmax=10)

    cid = 'G9.3'

    INCLUDE 'to_s_1d_generic.h'
  END FUNCTION to_s_1d_double_

  FUNCTION to_s_1d_complex_ (value) RESULT (out)
    COMPLEX :: value(:)
    CHARACTER(len=500) :: out
    CHARACTER(len=10) :: cid
    CHARACTER(len=100) :: form
    INTEGER :: l, u, d,n,i, nmax
    PARAMETER (d=3, nmax=10)

    cid = 'I0'

    INCLUDE 'to_s_1d_generic.h'
  END FUNCTION to_s_1d_complex_

  FUNCTION to_s_1d_logical_ (value) RESULT (out)
    LOGICAL :: value(:)
    CHARACTER(len=500) :: out
    CHARACTER(len=10) :: cid
    CHARACTER(len=100) :: form
    INTEGER :: l, u, d, n, nmax !,i
    CHARACTER(len=500) :: out
    PARAMETER (d=14, nmax=30)

    cid = 'L1'
    INCLUDE 'to_s_1d_generic.h'

  END FUNCTION to_s_1d_logical_


  FUNCTION to_s_2d_int_ (value) RESULT (out)
    INTEGER :: value(:,:)
    CHARACTER(len=1500) :: out
    CHARACTER(len=10) :: cid
    CHARACTER(len=100) :: form
    CHARACTER(len=200) :: temp(12)
    INTEGER :: l, u, d, n, i, k, ns, dots, nmax
    PARAMETER (d=3, nmax=10)

    cid = 'I0'

    INCLUDE 'to_s_2d_generic.h'

  END FUNCTION to_s_2d_int_


  FUNCTION to_s_2d_real_ (value) RESULT (out)
    REAL :: value(:,:)
    CHARACTER(len=1500) :: out
    CHARACTER(len=100) :: form
    CHARACTER(len=200) :: temp(12)
    INTEGER :: l, u, d, n, i, k, ns, dots, nmax
    PARAMETER (d=3, nmax=10)

    INCLUDE 'to_s_2d_generic.h'

  END FUNCTION to_s_2d_real_


  FUNCTION to_s_2d_double_ (value) RESULT (out)
    DOUBLE PRECISION :: value(:,:)
    CHARACTER(len=1500) :: out
    CHARACTER(len=100) :: form
    CHARACTER(len=200) :: temp(12)
    INTEGER :: l, u, d, n, i, k, ns, dots, nmax
    PARAMETER (d=3, nmax=10)

    INCLUDE 'to_s_2d_generic.h'

  END FUNCTION to_s_2d_double_


  FUNCTION to_s_2d_logical_ (value) RESULT (out)
    LOGICAL :: value(:,:)
    CHARACTER(len=1000) :: out
    CHARACTER(len=100) :: form
    CHARACTER(len=200) :: temp(12)
    INTEGER :: l, u, d, n, i, k, ns, dots, nmax=10
    PARAMETER (d=4)

    INCLUDE 'to_s_2d_generic.h'

  END FUNCTION to_s_2d_logical_


  FUNCTION strip(value)
    implicit none
    CHARACTER(len=2000):: strip
    CHARACTER(len=*), INTENT(in) :: value
    strip = TRIM(ADJUSTL(value))
  END FUNCTION strip

  !------------------------
  ! test if 2 values are close
  !------------------------
  !logical function equals (number1, number2)
  !  real,  intent (in) :: number1, number2
  !
  !  return equalEpsilon (number1, number2, epsilon(number1))
  !
  !end function equals


  function equalEpsilon (number1, number2, epsilon ) result (resultValue)
    real , intent (in) :: number1, number2, epsilon
    logical :: resultValue

    resultValue = .false.

    ! test very small number1
    if ( abs(number1) < epsilon .and.  abs(number1 - number2) < epsilon ) then
       resultValue = .true.
    else
       if ((abs(( number1 - number2)) / number1) < epsilon ) then
          resultValue = .true.
       else
          resultValue = .false.
       end if
    end if

  end function equalEpsilon

  function floatEqual (number1, number2 ) result (resultValue)
    real , intent (in) :: number1, number2
    real :: epsilon
    logical :: resultValue

    resultValue = .false.
    epsilon = 1E-6

    ! test very small number1
    if ( abs(number1) < epsilon .and.  abs(number1 - number2) < epsilon ) then
       resultValue = .true.
    else
       if ((abs(( number1 - number2)) / number1) < epsilon ) then
          resultValue = .true.
       else
          resultValue = .false.
       end if
    end if
  end function floatEqual

  function doublePrecisionEqual (number1, number2 ) result (resultValue)
    double precision , intent (in) :: number1, number2
    real :: epsilon
    logical :: resultValue

    resultValue = .false.
    epsilon = 1E-6
    !epsilon = epsilon (number1)

    ! test very small number1
    if ( abs(number1) < epsilon .and.  abs(number1 - number2) < epsilon ) then
       resultValue = .true.
    else
       if ((abs(( number1 - number2)) / number1) < epsilon ) then
          resultValue = .true.
       else
          resultValue = .false.
       end if
    end if
  end function doublePrecisionEqual

  function integerEqual (number1, number2 ) result (resultValue)
    integer , intent (in) :: number1, number2
    logical :: resultValue

    resultValue = .false.

    if ( number1 .eq. number2 ) then
       resultValue = .true.
    else
       resultValue = .false.
    end if
  end function integerEqual

  function stringEqual (str1, str2 ) result (resultValue)
    character(*) , intent (in) :: str1, str2
    logical :: resultValue

    resultValue = .false.

    if ( str1 .eq. str2 ) then
       resultValue = .true.
    end if
  end function stringEqual

  function logicalEqual (l1, l2 ) result (resultValue)
    logical, intent (in) :: l1, l2
    logical              :: resultValue

    resultValue = .false.

    if ( l1 .eqv. l2 ) then
       resultValue = .true.
    end if
  end function logicalEqual

end module fruit_util
