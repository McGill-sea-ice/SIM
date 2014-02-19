MODULE fruit
  
! ===========================
! FORTRAN unit test framework
! ===========================
!
! :Author: David Huard <david.huard@gmail.com>
! :Date: May 2008
! :Version: 0.1
! :Institution: McGill University
! :License: BSD
!
!
! The fruit module provides a set of utilities to write unit tests in fortran. 
! A unit test is a small unit of code that exercises a procedure and checks the 
! output against the known answer. To facilitate this boring but necessary step in 
! software engineering, the fruit module provides three basic interfaces to 
! assertion functions: 
!   * `assert_true(var, message)`, 
!   * `assert_equal(var1, var2, message)` and
!   * `assert_almost_equal(var1, var2, dec, message)`.
!  
! along with a framework to collect messages and print the result of the test suite. 
!
! The fruit module is a derived product of the FRUIT package written by Andrew H. Chen 
! (meihome @at@ gmail.com). I modified some of the routines and adapted the interface to 
! look more like the Python UnitTest that was my reference. The main differences are the 
! introduction of a separate `assert_almost_equal`, fancier output (large vectors and matrices
! are not printed out entirely), and hopefully a more understandable summary of the results and 
! a reviewed structure making it easier for users to extend the interface for their custom derived 
! types. 
!
! See `Example` for a simple example of fruit usage. 



  USE fruit_util
  USE assertions
  USE management

  INTERFACE assert_true
     MODULE PROCEDURE assert_true_logical_, assert_true_1d_logical_
  END INTERFACE

  INTERFACE assert_equal
     MODULE PROCEDURE assert_equal_int_
     MODULE PROCEDURE assert_equal_real_
     MODULE PROCEDURE assert_equal_double_
     MODULE PROCEDURE assert_equal_logical_
     MODULE PROCEDURE assert_equal_string_
     MODULE PROCEDURE assert_equal_complex_

     MODULE PROCEDURE assert_equal_1d_int_
     MODULE PROCEDURE assert_equal_1d_real_
     MODULE PROCEDURE assert_equal_1d_double_
!     module procedure assert_equal_1d_string_
!     module procedure assert_equal_1d_complex_

     MODULE PROCEDURE assert_equal_2d_int_
     MODULE PROCEDURE assert_equal_2d_real_
     MODULE PROCEDURE assert_equal_2d_double_
  !   module procedure assert_equal_2d_complex_
  END INTERFACE

INTERFACE assert_almost_equal
   ! assert_almost_equal(var1, var2, dec, message)
   !
   ! Check that both variables are almost equal.
   !
   ! Parameters
   ! ----------
   !   var1 : Value to check.
   !   var2 : Expected value. 
   !   dec : The required precision
   !     The test will pass only if the difference between the two variables
   !     is smaller than 10**dec.
   !   message : Output message. 
   MODULE PROCEDURE assert_almost_equal_int_
   MODULE PROCEDURE assert_almost_equal_real_
   MODULE PROCEDURE assert_almost_equal_1d_int_
   MODULE PROCEDURE assert_almost_equal_1d_real_
   MODULE PROCEDURE assert_almost_equal_2d_int_
   MODULE PROCEDURE assert_almost_equal_2d_real_
   MODULE PROCEDURE assert_almost_equal_double_
   MODULE PROCEDURE assert_almost_equal_1d_double_
   MODULE PROCEDURE assert_almost_equal_2d_double_

END INTERFACE



END MODULE fruit
