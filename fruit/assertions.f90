MODULE assertions

! This module contains all the assertions on intrinsic types. All procedures are public since this module
! is intended to be use from the main module where the interface is defined. 
!
! Two categories of assertions are implemented: `assert_equal(var1, var2, message)` and 
! `assert_almost_equal(var1, var2, dec, message)`, where `dec` is the decimal at which the comparison 
! takes place. For real and complex types, the default value is -6. For the double precision type, the default
! is set to -8. For integers, the default is 0, which means that a difference of 1 between two integers will 
! return a failed assertion. 
!
! Assertions are implemented for scalars, 1D and 2D arrays. 

  USE management, ONLY: success_assert_action, failed_assert_action
  USE fruit_util, ONLY: to_s, strip

CONTAINS


! =====================================
!    Logical .TRUE. assertions
! =====================================

  SUBROUTINE assert_true_logical_ (var1, message)
    LOGICAL, INTENT (in) :: var1
    CHARACTER (*), INTENT (in), OPTIONAL :: message
    
    IF ( var1 .EQV. .TRUE.) THEN
       CALL success_assert_action
    ELSE
       CALL failed_assert_action(to_s(var1), to_s(.TRUE.), message)
    END IF
  END SUBROUTINE assert_true_logical_
  
  SUBROUTINE assert_true_1d_logical_ (var1, message)
    LOGICAL, INTENT(in) :: var1(:)
    CHARACTER (*), INTENT (in), OPTIONAL :: message

    IF (ALL (var1 .EQV. .TRUE.)) THEN
       CALL success_assert_action()
    ELSE
       CALL failed_assert_action(to_s(var1), to_s(.TRUE.), message)
    END IF
  END SUBROUTINE assert_true_1d_logical_

  SUBROUTINE assert_true_2d_logical_ (var1, message)
    LOGICAL, INTENT(in) :: var1(:, :)
    CHARACTER (*), INTENT (in), OPTIONAL :: message

    IF (ALL (var1 .EQV. .TRUE.)) THEN
       CALL success_assert_action()
    ELSE
       CALL failed_assert_action(to_s(var1), to_s(.TRUE.), message)
    END IF
  END SUBROUTINE assert_true_2d_logical_



! =====================================
!    Logical equality assertions
! =====================================

  SUBROUTINE assert_equal_logical_ (var1, var2, message)
    LOGICAL, INTENT (in)  :: var1, var2
    CHARACTER (*), INTENT (in), OPTIONAL :: message

    IF ( var1 .EQV. var2 ) THEN
       CALL success_assert_action
    ELSE
       CALL failed_assert_action(to_s(var1), to_s(var2), message)
    END IF
  END SUBROUTINE assert_equal_logical_


  SUBROUTINE assert_equal_1d_logical_ (var1, var2, message)
    LOGICAL, INTENT(in) :: var1(:), var2(:)     ! The two variables to compare
    CHARACTER (*), INTENT(in), OPTIONAL :: message ! The message to print in case of failure

    INTEGER :: s1, s2, notequals
    LOGICAL, ALLOCATABLE :: mask(:)

    s1 = SIZE(var1)
    s2 = size(var2)

    ! Check that dimensions are identical
    IF (s1 /= s2) THEN
       CALL failed_assert_action(to_s(s1), to_s(s2), &
            message//' -- Shape of arrays do not match.')
       RETURN
    END IF
    
    ALLOCATE(mask(s1))
    
    ! Compute the difference between the two arrays.
    mask = (var1 .NEQV.  var2)
    notequals = COUNT(mask)
    
    IF (notequals>0) THEN        
       CALL failed_assert_action(to_s(var1), to_s(var2), message // ' -- ' // TRIM(to_s(notequals)) // ' elements are not equal.' )
    ELSE
       CALL success_assert_action
    ENDIF
    
  END SUBROUTINE assert_equal_1d_logical_


  SUBROUTINE assert_equal_2d_logical_ (var1, var2, message)
    LOGICAL, INTENT(in) :: var1(:,:), var2(:,:)     ! The two variables to compare
    CHARACTER (*), INTENT(in), OPTIONAL :: message ! The message to print in case of failure

    INTEGER :: s1(2), s2(2), notequals
    LOGICAL, ALLOCATABLE :: mask(:,:)

    s1 = SHAPE(var1)
    s2 = SHAPE(var2)

    ! Check that dimensions are identical
    IF (ANY(s1 /= s2)) THEN
       CALL failed_assert_action(to_s(s1), to_s(s2), &
            message//' -- Shape of arrays do not match.')
       RETURN
    END IF
    
    ALLOCATE(mask(s1(1), s1(2)))
    
    ! Compute the difference between the two arrays.
    mask = (var1 .NEQV. var2)
    notequals = COUNT(mask)
    
    IF (notequals>0) THEN        
       CALL failed_assert_action(to_s(var1), to_s(var2), message // ' -- ' // TRIM(to_s(notequals))// ' elements are not equal.' )
    ELSE
       CALL success_assert_action
    ENDIF


  END SUBROUTINE assert_equal_2d_logical_



! =====================================
!    Integer equality assertions
! =====================================

  SUBROUTINE assert_equal_int_ (var1, var2, message)
    INTEGER, INTENT(in) :: var1, var2     ! The two variables to compare
    CHARACTER (*), INTENT(in), OPTIONAL :: message ! The message to print in case of failure

    INCLUDE 'assert_equal_generic.h'
  END SUBROUTINE assert_equal_int_


  SUBROUTINE assert_almost_equal_int_ (var1, var2, dec, message)
    INTEGER, INTENT(in) :: var1, var2     ! The two variables to compare
    INTEGER, INTENT(in), optional :: dec  ! The acceptable precision in decimal places.
    CHARACTER (*), INTENT(in), OPTIONAL :: message ! The message to print in case of failure
    INTEGER :: odec                       ! The default required precision. 

    ! Define the default value for dec
    odec = 0  

    INCLUDE 'assert_almost_equal_generic.h'
  END SUBROUTINE assert_almost_equal_int_
  

  SUBROUTINE assert_equal_1d_int_ (var1, var2, message)
    INTEGER, INTENT(in) :: var1(:), var2(:)     ! The two variables to compare
    CHARACTER (*), INTENT(in), OPTIONAL :: message ! The message to print in case of failure

    INTEGER :: s1, s2, notequals   
    LOGICAL, ALLOCATABLE :: mask(:)

    INCLUDE 'assert_equal_1d_generic.h'
  END SUBROUTINE assert_equal_1d_int_


  SUBROUTINE assert_almost_equal_1d_int_ (var1, var2, dec, message)
    INTEGER, INTENT(in) :: var1(:), var2(:)     ! The two variables to compare
    INTEGER, INTENT(in), optional :: dec  ! The acceptable precision in decimal places.
    CHARACTER (*), INTENT(in), OPTIONAL :: message ! The message to print in case of failure
    INTEGER :: odec                       ! The default required precision. 

    INTEGER :: s1, s2, notequals
    REAL, ALLOCATABLE :: diff(:)
    LOGICAL, ALLOCATABLE :: mask(:)

    ! Define the default value for dec
    odec = 0

    INCLUDE 'assert_almost_equal_1d_generic.h'
  END SUBROUTINE assert_almost_equal_1d_int_
  

  SUBROUTINE assert_equal_2d_int_ (var1, var2, message)
    INTEGER, INTENT(in) :: var1(:,:), var2(:,:)     ! The two variables to compare
    CHARACTER (*), INTENT(in), OPTIONAL :: message ! The message to print in case of failure

    INTEGER :: s1(2), s2(2), notequals
    LOGICAL, ALLOCATABLE :: mask(:,:)

    INCLUDE 'assert_equal_2d_generic.h'
  END SUBROUTINE assert_equal_2d_int_
  

  SUBROUTINE assert_almost_equal_2d_int_ (var1, var2, dec, message)
    INTEGER, INTENT(in) :: var1(:,:), var2(:,:)     ! The two variables to compare
    INTEGER, INTENT(in), optional :: dec  ! The acceptable precision in decimal places.
    CHARACTER (*), INTENT(in), OPTIONAL :: message ! The message to print in case of failure
    INTEGER :: odec                       ! The default required precision. 

    INTEGER :: s1(2), s2(2), notequals
    REAL, ALLOCATABLE :: diff(:,:)
    LOGICAL, ALLOCATABLE :: mask(:,:)

    ! Define the default value for dec
    odec = 0

    INCLUDE 'assert_almost_equal_2d_generic.h'
  END SUBROUTINE assert_almost_equal_2d_int_
  
  
! ===================================== 
!    Real assertions
! =====================================

  SUBROUTINE assert_equal_real_ (var1, var2, message)
    REAL, INTENT (in) :: var1, var2
    CHARACTER (*), INTENT(in), OPTIONAL :: message ! The message to print in case of failure

    INCLUDE 'assert_equal_generic.h'
  END SUBROUTINE assert_equal_real_


  SUBROUTINE assert_almost_equal_real_ (var1, var2, dec, message)
    REAL, INTENT (in) :: var1, var2
    INTEGER, INTENT(in), optional :: dec  ! The acceptable precision in decimal places.
    CHARACTER (*), INTENT(in), OPTIONAL :: message ! The message to print in case of failure
    INTEGER :: odec                       ! The default required precision. 

    ! Define the default value for dec
    odec = -6

    INCLUDE 'assert_almost_equal_generic.h'
  END SUBROUTINE assert_almost_equal_real_
  

  SUBROUTINE assert_equal_1d_real_ (var1, var2, message)
    REAL, INTENT (in) :: var1(:), var2(:)
    CHARACTER (*), INTENT(in), OPTIONAL :: message ! The message to print in case of failure

    INTEGER :: s1, s2, notequals
    LOGICAL, ALLOCATABLE :: mask(:)

    INCLUDE 'assert_equal_1d_generic.h'
  END SUBROUTINE assert_equal_1d_real_


  SUBROUTINE assert_almost_equal_1d_real_ (var1, var2, dec, message)
    REAL, INTENT (in) :: var1(:), var2(:)
    INTEGER, INTENT(in), optional :: dec  ! The acceptable precision in decimal places.
    CHARACTER (*), INTENT(in), OPTIONAL :: message ! The message to print in case of failure
    INTEGER :: odec                       ! The default required precision. 

    INTEGER :: s1, s2, notequals
    REAL, ALLOCATABLE :: diff(:)
    LOGICAL, ALLOCATABLE :: mask(:)

    ! Define the default value for dec
    odec = -6

    INCLUDE 'assert_almost_equal_1d_generic.h'
  END SUBROUTINE assert_almost_equal_1d_real_
  

  SUBROUTINE assert_equal_2d_real_ (var1, var2, message)
    REAL, INTENT (in) :: var1(:,:), var2(:,:)
    CHARACTER (*), INTENT(in), OPTIONAL :: message ! The message to print in case of failure

    INTEGER :: s1(2), s2(2), notequals
    LOGICAL, ALLOCATABLE :: mask(:,:)

    INCLUDE 'assert_equal_2d_generic.h'
  END SUBROUTINE assert_equal_2d_real_

  
  SUBROUTINE assert_almost_equal_2d_real_ (var1, var2, dec, message)
    REAL, INTENT (in) :: var1(:,:), var2(:,:)
    INTEGER, INTENT(in), optional :: dec  ! The acceptable precision in decimal places.
    CHARACTER (*), INTENT(in), OPTIONAL :: message ! The message to print in case of failure
    INTEGER :: odec                       ! The default required precision. 

    INTEGER :: s1(2), s2(2), notequals
    REAL, ALLOCATABLE :: diff(:,:)
    LOGICAL, ALLOCATABLE :: mask(:,:)

    ! Define the default value for dec
    odec = -6

    INCLUDE 'assert_almost_equal_2d_generic.h' 
  END SUBROUTINE assert_almost_equal_2d_real_
   
    
! =====================================
!    Double precision assertions
! =====================================

  SUBROUTINE assert_equal_double_ (var1, var2, message)
    DOUBLE PRECISION, INTENT (in) :: var1, var2
    CHARACTER (*), INTENT(in), OPTIONAL :: message ! The message to print in case of failure

    INCLUDE 'assert_equal_generic.h'
  END SUBROUTINE assert_equal_double_
 

  SUBROUTINE assert_almost_equal_double_ (var1, var2, dec, message)
    DOUBLEPRECISION, INTENT (in) :: var1, var2
    INTEGER, INTENT(in), optional :: dec  ! The acceptable precision in decimal places.
    CHARACTER (*), INTENT(in), OPTIONAL :: message ! The message to print in case of failure
    INTEGER :: odec                       ! The default required precision. 

    ! Define the default value for dec
    odec = -9

    INCLUDE 'assert_almost_equal_generic.h'
  END SUBROUTINE assert_almost_equal_double_


  SUBROUTINE assert_equal_1d_double_ (var1, var2, message)
    DOUBLEPRECISION, INTENT (in) :: var1(:), var2(:)
    CHARACTER (*), INTENT(in), OPTIONAL :: message ! The message to print in case of failure

    INTEGER :: s1, s2, notequals
    LOGICAL, ALLOCATABLE :: mask(:)

    INCLUDE 'assert_equal_1d_generic.h'
  END SUBROUTINE assert_equal_1d_double_


  SUBROUTINE assert_almost_equal_1d_double_ (var1, var2, dec, message)
    DOUBLEPRECISION, INTENT (in) :: var1(:), var2(:)
    INTEGER, INTENT(in), optional :: dec  ! The acceptable precision in decimal places.
    CHARACTER (*), INTENT(in), OPTIONAL :: message ! The message to print in case of failure
    INTEGER :: odec                       ! The default required precision. 

    INTEGER :: s1, s2, notequals
    DOUBLE PRECISION, ALLOCATABLE :: diff(:)
    LOGICAL, ALLOCATABLE :: mask(:)
 
    ! Define the default value for dec
    odec = -9

    INCLUDE 'assert_almost_equal_1d_generic.h'
  END SUBROUTINE assert_almost_equal_1d_double_


  SUBROUTINE assert_equal_2d_double_ (var1, var2, message)
    DOUBLEPRECISION, INTENT (in) :: var1(:,:), var2(:,:)
    CHARACTER (*), INTENT(in), OPTIONAL :: message ! The message to print in case of failure

    INTEGER :: s1(2), s2(2), notequals
    LOGICAL, ALLOCATABLE :: mask(:,:)

    INCLUDE 'assert_equal_2d_generic.h'
  END SUBROUTINE assert_equal_2d_double_

  SUBROUTINE assert_almost_equal_2d_double_ (var1, var2, dec, message)
    DOUBLEPRECISION, INTENT (in) :: var1(:,:), var2(:,:)
    INTEGER, INTENT(in), optional :: dec  ! The acceptable precision in decimal places.
    CHARACTER (*), INTENT(in), OPTIONAL :: message ! The message to print in case of failure
    INTEGER :: odec                       ! The default required precision. 

    INTEGER :: s1(2), s2(2), notequals
    REAL, ALLOCATABLE :: diff(:,:)
    LOGICAL, ALLOCATABLE :: mask(:,:)

    ! Define the default value for dec
    odec = -9

    INCLUDE 'assert_almost_equal_2d_generic.h'
  END SUBROUTINE assert_almost_equal_2d_double_


! =====================================
!    Complex assertions
! =====================================

  SUBROUTINE assert_equal_complex_ (var1, var2, message)
    COMPLEX, INTENT (in) :: var1, var2
    CHARACTER (*), INTENT(in), OPTIONAL :: message ! The message to print in case of failure

    INCLUDE 'assert_equal_generic.h'
  END SUBROUTINE assert_equal_complex_


!!$  SUBROUTINE assert_almost_equal_complex_ (var1, var2, dec, message)
!!$    COMPLEX, INTENT (in) :: var1, var2
!!$    INTEGER, INTENT(in), optional :: dec  ! The acceptable precision in decimal places.
!!$    CHARACTER (*), INTENT(in), OPTIONAL :: message ! The message to print in case of failure
!!$    INTEGER :: odec                       ! The default required precision. 
!!$
!!$    ! Define the default value for dec
!!$    odec = -6
!!$
!!$    INCLUDE 'assert_almost_equal_generic.h'
!!$  END SUBROUTINE assert_almost_equal_complex_


  SUBROUTINE assert_equal_1d_complex_ (var1, var2, message)
    COMPLEX, INTENT (in) :: var1(:), var2(:)
    CHARACTER (*), INTENT(in), OPTIONAL :: message ! The message to print in case of failure

    INTEGER :: s1, s2, notequals
    LOGICAL, ALLOCATABLE :: mask(:)

    INCLUDE 'assert_equal_1d_generic.h'
  END SUBROUTINE assert_equal_1d_complex_


!!$  SUBROUTINE assert_almost_equal_1d_complex_ (var1, var2, dec, message)
!!$    COMPLEX, INTENT (in) :: var1(:), var2(:)
!!$    INTEGER, INTENT(in), optional :: dec  ! The acceptable precision in decimal places.
!!$    CHARACTER (*), INTENT(in), OPTIONAL :: message ! The message to print in case of failure
!!$    INTEGER :: odec                       ! The default required precision. 
!!$
!!$    INTEGER :: s1, s2, notequals
!!$    REAL, ALLOCATABLE :: diff(:)
!!$    LOGICAL, ALLOCATABLE :: mask(:)
!!$
!!$    ! Define the default value for dec
!!$    odec = -6
!!$
!!$    INCLUDE 'assert_almost_equal_1d_generic.h'
!!$  END SUBROUTINE assert_almost_equal_1d_complex_



! =====================================
!    String assertions
! =====================================

  SUBROUTINE assert_equal_string_ (var1, var2, message)
    CHARACTER(*), INTENT (in)  :: var1, var2
    CHARACTER (*), INTENT (in), OPTIONAL :: message
    
    IF ( TRIM(strip(var1)) == TRIM(strip(var2))) THEN
       CALL success_assert_action
    ELSE
       CALL failed_assert_action(var1, var2, message)
    END IF
  END SUBROUTINE assert_equal_string_

!!$
!!$  SUBROUTINE assert_equal_1d_string_ (var1, var2, message)
!!$    CHARACTER(*), INTENT (in) :: var1(:), var2(:)
!!$    CHARACTER (*), INTENT (in), OPTIONAL :: message
!!$    INTEGER :: n1, n2
!!$    
!!$    n1 = size(var1)
!!$    n2 = size(var2)
!!$    ! Check that dimensions are identical
!!$    IF (n1 /= n2) THEN
!!$       CALL failed_assert_action(to_s(n1), to_s(n2), &
!!$            message//' -- Length of strings do not match.')
!!$       RETURN
!!$    END IF
!!$    
!!$    IF (trim(var1)  == trim(var2)) THEN 
!!$       CALL success_assert_action
!!$    ELSE
!!$       CALL failed_assert_action(var1, var2, message)
!!$    ENDIF
!!$  END SUBROUTINE assert_equal_1d_string_


! =====================================
!    2D array assertions
! =====================================

! Test equality
! -------------

!!$
!!$
!!$
!!$
!!$  SUBROUTINE assert_equal_2d_complex_ (var1, var2, dec, message)
!!$    COMPLEX, INTENT (in) :: var1(:), var2(:)
!!$    INTEGER, INTENT(in), optional :: dec  ! The acceptable precision in decimal places.
!!$    CHARACTER (*), INTENT(in), OPTIONAL :: message ! The message to print in case of failure
!!$    INTEGER :: odec                       ! The default required precision. 
!!$
!!$    INTEGER :: s1, s2, notequals
!!$    REAL, ALLOCATABLE :: diff(:)
!!$    LOGICAL, ALLOCATABLE :: mask(:)
!!$
!!$    ! Define the default value for dec
!!$    odec = -6
!!$
!!$    INCLUDE 'assert_equal_2d_generic.txt'
!!$  END SUBROUTINE assert_equal_2d_complex_


! Test near equality
! ------------------

!!$
!!$

!!$
!!$  SUBROUTINE assert_equal_2d_complex_ (var1, var2, dec, message)
!!$    COMPLEX, INTENT (in) :: var1(:), var2(:)
!!$    INTEGER, INTENT(in), optional :: dec  ! The acceptable precision in decimal places.
!!$    CHARACTER (*), INTENT(in), OPTIONAL :: message ! The message to print in case of failure
!!$    INTEGER :: odec                       ! The default required precision. 
!!$
!!$    INTEGER :: s1, s2, notequals
!!$    REAL, ALLOCATABLE :: diff(:)
!!$    LOGICAL, ALLOCATABLE :: mask(:)
!!$
!!$    ! Define the default value for dec
!!$    odec = -6
!!$
!!$    INCLUDE 'assert_equal_2d_generic.txt'
!!$  END SUBROUTINE assert_equal_2d_complex_

  
!=====================
!     UTILITIES
!=====================

END MODULE assertions
