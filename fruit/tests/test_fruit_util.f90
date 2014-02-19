module test_fruit_util

  use fruit
  implicit none

contains

  subroutine test_to_s_should_convert_int_to_string
    call set_unit_name("int to string")
    call assert_equal (to_s(1), '1')
    call assert_equal (to_s(-1), '-1')
    call assert_Equal (to_s(0), '0')
  end subroutine test_to_s_should_convert_int_to_string

  subroutine test_to_s_should_convert_real_into_string
    call set_unit_name("real to string")
    call assert_Equal (to_s(1.0), '1.00000000')
    call assert_Equal (to_s(123456789.123), '1.23456792E+08')
    call assert_Equal (to_s(0.0), '0.0000000')
  end subroutine test_to_s_should_convert_real_into_string

  subroutine test_to_s_should_convert_logical_to_string
    call set_unit_name("logical to string")
    call assert_Equal (to_s(.true.), 'T')
    call assert_Equal (to_s(.false.), 'F')
  end subroutine test_to_s_should_convert_logical_to_string

  subroutine test_string_strip_should_remove_starting_and_ending_spaces
    call set_unit_name('test_strip')
    call assert_Equal (strip('  1.0   '), '1.0')
    call assert_Equal (strip('abc  '), 'abc')
    call assert_Equal (strip('   def'), 'def')
    call assert_Equal (strip(to_s(2)), '2')
  end subroutine test_string_strip_should_remove_starting_and_ending_spaces

  subroutine fruit_util_test_suite
    call test_to_s_should_convert_int_to_string
    call test_to_s_should_convert_real_into_string
    call test_to_s_should_convert_logical_to_string
    call test_string_strip_should_remove_starting_and_ending_spaces
  end subroutine fruit_util_test_suite

end module test_fruit_util


program test

  use test_fruit_util
  use fruit
  
  call init_fruit()
  call fruit_util_test_suite()
  call fruit_summary()

end program test
