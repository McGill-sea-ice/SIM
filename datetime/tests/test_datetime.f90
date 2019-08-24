 PROGRAM datetime_test

  USE fruit_ext
  USE test_datetimedelta
  USE test_date
  USE test_datetime
  USE test_time
  use test_dateutils

  CALL init_fruit

  CALL delta_test_suite()
  CALL date_test_suite()
  CALL datetime_test_suite()
  CALL time_test_suite()
  call dateutils_test_suite()



  CALL fruit_summary

END PROGRAM datetime_test

