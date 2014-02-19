MODULE test_assert

  USE fruit, ONLY: set_unit_name
  USE assertions

CONTAINS

  ! Test logical true assertions
  ! ----------------------------

  SUBROUTINE test_assert_true_logical_()
    CALL set_unit_name('test_assert_true_logical_')
    CALL assert_true_logical_(.TRUE., 'true')
    CALL assert_true_logical_(.FALSE., 'false')
  END SUBROUTINE test_assert_true_logical_


  SUBROUTINE test_assert_true_1d_logical_()
    CALL set_unit_name('test_assert_true_1d_logical_')
    CALL assert_true_1d_logical_((/.TRUE., .TRUE./), 'true')
    CALL assert_true_1d_logical_((/.FALSE., .TRUE./), 'false')
  END SUBROUTINE test_assert_true_1d_logical_


  SUBROUTINE test_assert_true_2d_logical_()
    LOGICAL :: a(2,2), b(2,2)
    a = .TRUE.
    b = RESHAPE((/.TRUE., .FALSE., .TRUE., .TRUE.  /), (/2,2/))
    CALL set_unit_name('test_assert_true_2d_logical_')
    CALL assert_true_2d_logical_(a, 'true')
    CALL assert_true_2d_logical_(b, 'false')
  END SUBROUTINE test_assert_true_2d_logical_


! Test logical equality assertions
! --------------------------------

  SUBROUTINE test_assert_equal_logical_()
    CALL set_unit_name('test_assert_equal_logical_')
    CALL assert_equal_logical_(.TRUE., .TRUE., 'equal')
    CALL assert_equal_logical_(.TRUE., .FALSE., 'not equal')
  END SUBROUTINE test_assert_equal_logical_


  SUBROUTINE test_assert_equal_1d_logical_()
    CALL set_unit_name('test_assert_equal_1d_logical_')
    CALL assert_equal_1d_logical_((/.TRUE., .FALSE./), (/.TRUE., .FALSE./), 'equal')
    CALL assert_equal_1d_logical_((/.TRUE., .FALSE./), (/.TRUE., .TRUE. /), 'not equal')
  END SUBROUTINE test_assert_equal_1d_logical_


  SUBROUTINE test_assert_equal_2d_logical_()
    LOGICAL :: a(2,2), b(2,2)
    a = RESHAPE((/.TRUE., .FALSE., .TRUE., .FALSE. /), (/2,2/))
    b = RESHAPE((/.TRUE., .FALSE., .TRUE., .TRUE.  /), (/2,2/))
    CALL set_unit_name('test_assert_equal_2d_logical_')
    CALL assert_equal_2d_logical_(a, a, 'equal')
    CALL assert_equal_2d_logical_(a, b, 'not_equal')
  END SUBROUTINE test_assert_equal_2d_logical_


  ! Test integer equality assertions
  ! --------------------------------

  SUBROUTINE test_assert_equal_int_()
    CALL set_unit_name('test_assert_equal_int_')
    CALL assert_equal_int_(5, 5, message='equal')
    CALL assert_equal_int_(4, 5, message='not equal')
  END SUBROUTINE test_assert_equal_int_


  SUBROUTINE test_assert_almost_equal_int_()
    CALL set_unit_name('test_assert_almost_equal_int_')
    CALL assert_almost_equal_int_(4, 5, 1, message='almost equal')
    CALL assert_almost_equal_int_(4, 5, 0, message='not almost equal')
  END SUBROUTINE test_assert_almost_equal_int_


  SUBROUTINE test_assert_equal_1d_int_()
    CALL set_unit_name('test_assert_equal_1d_int_')
    CALL assert_equal_1d_int_((/5, 5/), (/5, 5/), message='all equal')
    CALL assert_equal_1d_int_((/4, 5/), (/5, 5/), message='not all equal')
  END SUBROUTINE test_assert_equal_1d_int_


  SUBROUTINE test_assert_almost_equal_1d_int_()
    CALL set_unit_name('test_assert_almost_equal_1d_int_')
    CALL assert_almost_equal_1d_int_((/5, 5/), (/5, 4/), 1, message='all almost equal')
    CALL assert_almost_equal_1d_int_((/4, 5/), (/5, 5/), 0, message='not all almost equal')
  END SUBROUTINE test_assert_almost_equal_1d_int_
    

  SUBROUTINE test_assert_equal_2d_int_()
    INTEGER :: a(2,2), b(2,2)
    a = RESHAPE((/1,1,1,1/), (/2,2/))
    b = RESHAPE((/1,1,1,2/), (/2,2/))
    CALL set_unit_name('test_assert_equal_2d_int_')
    CALL assert_equal_2d_int_(a, a, message='all equal')
    CALL assert_equal_2d_int_(a, b, message='not all equal')
  END SUBROUTINE test_assert_equal_2d_int_

  SUBROUTINE test_assert_almost_equal_2d_int_()
    INTEGER :: a(2,2), b(2,2)
    a = RESHAPE((/1,1,1,1/), (/2,2/))
    b = RESHAPE((/1,1,1,2/), (/2,2/))
    CALL set_unit_name('test_assert_almost_equal_2d_int_')
    CALL assert_almost_equal_2d_int_(a, b, 1, message='all almost equal')
    CALL assert_almost_equal_2d_int_(a, b, 0, message='not all almost equal')
  END SUBROUTINE test_assert_almost_equal_2d_int_


  ! Test real equality assertions
  ! ----------------------------- 

  SUBROUTINE test_assert_equal_real_()
    CALL set_unit_name('test_assert_equal_real_')
    CALL assert_equal_real_(5., 5., message='equal')
    CALL assert_equal_real_(4.9, 5., message='not equal')
  END SUBROUTINE test_assert_equal_real_


  SUBROUTINE test_assert_almost_equal_real_()
    CALL set_unit_name('test_assert_almost_equal_real_')
    CALL assert_almost_equal_real_(4.9, 5., 0,  message='equal')
    CALL assert_almost_equal_real_(4.9, 5., -1,  message='not equal')
    CALL assert_almost_equal_real_(5.1, 5., 0,  message='equal')
    CALL assert_almost_equal_real_(5.01, 5., -2,  message='not equal')

  END SUBROUTINE test_assert_almost_equal_real_

    
  SUBROUTINE test_assert_equal_1d_real_()
    CALL set_unit_name('test_assert_equal_1d_real_')
    CALL assert_equal_1d_real_((/5., 5./), (/5., 5./), message='all equal')
    CALL assert_equal_1d_real_((/4.9, 5./), (/5., 5./), message='not all equal')
  END SUBROUTINE test_assert_equal_1d_real_
    

  SUBROUTINE test_assert_almost_equal_1d_real_()
    CALL set_unit_name('test_assert_almsot_equal_1d_real_')
    CALL assert_almost_equal_1d_real_((/4.9, 5./), (/5., 5./),0, message='all equal')
    CALL assert_almost_equal_1d_real_((/4.9, 5./), (/5., 5./),-1, message='not all equal')
  END SUBROUTINE test_assert_almost_equal_1d_real_
        
    
  SUBROUTINE test_assert_equal_2d_real_()
    REAL :: a(2,2), b(2,2)
    a = RESHAPE((/1.,1.,1.,1./), (/2,2/))
    b = RESHAPE((/1.,1.,1.,1.9/), (/2,2/))
    CALL set_unit_name('test_assert_equal_2d_real_')
    CALL assert_equal_2d_real_(a, a, message='all equal')
    CALL assert_equal_2d_real_(a, b, message='not all equal')
  END SUBROUTINE test_assert_equal_2d_real_


  SUBROUTINE test_assert_almost_equal_2d_real_()
    REAL :: a(2,2), b(2,2)
    a = RESHAPE((/1.,1.,1.,1./), (/2,2/))
    b = RESHAPE((/1.,1.,1.,1.1/), (/2,2/))
    CALL set_unit_name('test_assert_almost_equal_2d_real_')
    CALL assert_almost_equal_2d_real_(a, b, 0, message='all equal')
    CALL assert_almost_equal_2d_real_(a, b, -1, message='not all equal')
  END SUBROUTINE test_assert_almost_equal_2d_real_

  ! Test double equality assertions
  ! ----------------------------- 

  SUBROUTINE test_assert_equal_double_()
    DOUBLE PRECISION :: a, b
    a = 5.
    b = 5.1
    CALL set_unit_name('test_assert_equal_double_')
    CALL assert_equal_double_(a, a, message='equal')
    CALL assert_equal_double_(a, b, message='not equal')
  END SUBROUTINE test_assert_equal_double_


  SUBROUTINE test_assert_almost_equal_double_()
    DOUBLE PRECISION :: a, b
    a = 5.000
    b = 5.001
    CALL set_unit_name('test_assert_almost_equal_double_')
    CALL assert_almost_equal_double_(a, b, -2,  message='equal')
    CALL assert_almost_equal_double_(a, b, -3,  message='not equal')
  END SUBROUTINE test_assert_almost_equal_double_
  

  SUBROUTINE test_assert_equal_1d_double_()
    DOUBLE PRECISION :: a(2), b(2)
    a = (/5.0d0, 4.0d0/)
    b= (/5.0d0, 5.0d0/)
    CALL set_unit_name('test_assert_equal_1d_double_')
    CALL assert_equal_1d_double_(a, a, message='all equal')
    CALL assert_equal_1d_double_(a, b,  message='not all equal')
  END SUBROUTINE test_assert_equal_1d_double_
  

  SUBROUTINE test_assert_almost_equal_1d_double_()
    DOUBLE PRECISION :: a(2), b(2)
    a = (/5.0d0, 5.0d0/)
    b= (/5.0d0, 5.001d0/)
    CALL set_unit_name('test_assert_almost_equal_1d_double_')
    CALL assert_almost_equal_1d_double_(a, b, -2, message='all equal')
    CALL assert_almost_equal_1d_double_(a, b, -3,  message='not all equal')
  END SUBROUTINE test_assert_almost_equal_1d_double_
    

  SUBROUTINE test_assert_equal_2d_double_()
    DOUBLE PRECISION :: a(2,2), b(2,2)
    a = RESHAPE((/1.,1.,1.,1./), (/2,2/))
    b = RESHAPE((/1.,1.,1.,1.1/), (/2,2/))
    CALL set_unit_name('test_assert_equal_2d_double_')
    CALL assert_equal_2d_double_(a, a, message='all equal')
    CALL assert_equal_2d_double_(a, b, message='not all equal')
  END SUBROUTINE test_assert_equal_2d_double_


  SUBROUTINE test_assert_almost_equal_2d_double_()
    DOUBLE PRECISION :: a(2,2), b(2,2)
    a = RESHAPE((/1.,1.,1.,1./), (/2,2/))
    b = RESHAPE((/1.,1.,1.,1.1/), (/2,2/))
    CALL set_unit_name('test_assert_almost_equal_2d_double_')
    CALL assert_almost_equal_2d_double_(a, b, 0, message='all equal')
    CALL assert_almost_equal_2d_double_(a, b, -1, message='not all equal')
  END SUBROUTINE test_assert_almost_equal_2d_double_


  ! Test complex equality assertions
  ! ----------------------------- 

  SUBROUTINE test_assert_equal_complex_()
    COMPLEX :: a,b
    a = (1,-1)
    b = (-1,1)
    CALL set_unit_name('test_assert_equal_complex_')
    CALL assert_equal_complex_(a, a, message='equal')
    CALL assert_equal_complex_(a, b, message='not equal')
  END SUBROUTINE test_assert_equal_complex_


  ! Test string equality assertions
  ! ----------------------------- 

  SUBROUTINE test_assert_equal_string_()
    CHARACTER(6) :: a,b
    a = 'tata'
    b = 'toto'
    call set_unit_name('test_assert_equal_string_')
    CALL assert_equal_string_(a, a, message='equal')
    CALL assert_equal_string_(a, b, message='not equal')
  END SUBROUTINE test_assert_equal_string_
    

END MODULE test_assert




PROGRAM test_assertions

! This tests the success and failures of assertions. 
! There should be exactly one success for each failure. 
! That is, the results should look like
! .F.F.F.F.F.F  etc.


  USE test_assert
  USE fruit
  
  LOGICAL :: test_truth, test_logical, test_integer, test_real, test_double
  LOGICAL :: test_complex, test_string, test_all

  test_truth = .FALSE.
  test_logical = .FALSE.
  test_integer = .FALSE.
  test_real = .FALSE.
  test_double = .false.
  test_complex = .FALSE.
  test_string = .FALSE.
  test_all = .true.

  WRITE(*,*) 'The results from this test should be a sequence of one success and one failure.' 


  CALL init_fruit()
  
! Logical truth assertions
! ------------------------
  IF (test_truth .or. test_all) THEN
     CALL test_assert_true_logical_()
     CALL test_assert_true_1d_logical_()
     CALL test_assert_true_2d_logical_()
  END IF

! Logical equality assertions
! ---------------------------
  IF (test_logical .or. test_all) THEN
     CALL test_assert_equal_logical_()
     CALL test_assert_equal_1d_logical_()
     CALL test_assert_equal_2d_logical_()
  END IF

! Integer equality assertions
! ---------------------------
  IF (test_integer .or. test_all) THEN
     CALL test_assert_equal_int_()
     CALL test_assert_almost_equal_int_()
     CALL test_assert_equal_1d_int_()
     CALL test_assert_almost_equal_1d_int_()
     CALL test_assert_equal_2d_int_()
     CALL test_assert_almost_equal_2d_int_()
  END IF

! Real equality assertions
! ------------------------
  IF (test_real .or. test_all) THEN
     CALL test_assert_equal_real_()
     CALL test_assert_almost_equal_real_()
     CALL test_assert_equal_1d_real_()
     CALL test_assert_almost_equal_1d_real_()
     CALL test_assert_equal_2d_real_()
     CALL test_assert_almost_equal_2d_real_()
  END IF

! Double equality assertions
! ------------------------
  IF (test_double .or. test_all) THEN
     CALL test_assert_equal_double_()
     CALL test_assert_almost_equal_double_()
     CALL test_assert_equal_1d_double_()
     CALL test_assert_almost_equal_1d_double_()
     CALL test_assert_equal_2d_double_()
     CALL test_assert_almost_equal_2d_double_()
  END IF

! Complex equality assertions
! ------------------------
  IF (test_complex .or. test_all) THEN
     CALL test_assert_equal_complex_()
  END IF

! String equality assertions
! ------------------------
  IF (test_string .or. test_all) THEN
     CALL test_assert_equal_string_()
  END IF

  CALL fruit_summary()


END PROGRAM test_assertions
