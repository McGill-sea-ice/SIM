MODULE my_functions

CONTAINS

  FUNCTION my_sum(a,b)
    REAL, INTENT(in) :: a, b
    REAL :: my_sum
    my_sum = a+b
  END FUNCTION my_sum
END MODULE my_functions


MODULE test_suite

  USE my_functions
  USE fruit

CONTAINS

  SUBROUTINE test_my_sum()
    call set_unit_name('test_my_sum')
    call assert_equal(my_sum(1., 2.), 3., message='1+2=3')
    call assert_equal(my_sum(2., 2.), 3., message='2+2/=3 !')
  END SUBROUTINE test_my_sum

END MODULE test_suite

PROGRAM TEST

  use fruit
  use test_suite

  call init_fruit()
  call test_my_sum()
  call fruit_summary()

END PROGRAM TEST

