MODULE test_dateutils


  USE datetime
  USE fruit_ext

  IMPLICIT NONE

CONTAINS


  subroutine test_modulomilli()
    integer :: sec, rem

    call set_unit_name('test_modulomilli')

    call modulomilli(0, rem, sec)
    call assert_equal(rem, 0)
    call assert_equal(sec, 0)

    call modulomilli(1000, rem, sec)
    call assert_equal(rem, 0)
    call assert_equal(sec, 1)

    call modulomilli(1001, rem, sec)
    call assert_equal(rem, 1)
    call assert_equal(sec, 1)

    call modulomilli(-1001, rem, sec)
    call assert_equal(rem, 999)
    call assert_equal(sec, -2)

  end subroutine test_modulomilli


  SUBROUTINE test_modulosecond()
    INTEGER :: min, rem


    CALL set_unit_name('test_modulosecond')


    CALL modulosecond(0, rem, min)
    CALL assert_equal(rem, 0)
    CALL assert_equal(min, 0)

    CALL modulosecond(30, rem, min)
    CALL assert_equal(rem, 30)
    CALL assert_equal(min, 0)

    CALL modulosecond(60, rem, min)
    CALL assert_equal(rem, 0)
    CALL assert_equal(min, 1)

    CALL modulosecond(123, rem, min)
    CALL assert_equal(rem, 3)
    CALL assert_equal(min, 2)

    CALL modulosecond(-5, rem, min)
    CALL assert_equal(rem, 55)
    CALL assert_equal(min, -1)

    CALL modulosecond(-65, rem, min)
    CALL assert_equal(rem, 55)
    CALL assert_equal(min, -2)

    CALL modulosecond(60*60*24, rem, min)
    CALL assert_equal(rem, 0)
    CALL assert_equal(min, 60*24)

  END SUBROUTINE test_modulosecond



  SUBROUTINE test_modulominute()
    INTEGER :: hour, min

    CALL modulominute(0, min, hour)
    CALL assert_equal((/min, hour/), (/0,0/))

    CALL modulominute(30, min, hour)
    CALL assert_equal((/min, hour/), (/30,0/))

    CALL modulominute(60, min, hour)
    CALL assert_equal((/min, hour/), (/0,1/))

    CALL modulominute(123, min, hour)
    CALL assert_equal((/min, hour/), (/3,2/))

    CALL modulominute(-30, min, hour)
    CALL assert_equal((/min, hour/), (/30,-1/))

    CALL modulominute(-90, min, hour)
    CALL assert_equal((/min, hour/), (/30,-2/))

  END SUBROUTINE test_modulominute



  SUBROUTINE test_modulohour()
    INTEGER :: hour, day

    CALL modulohour(0, hour, day)
    CALL assert_equal((/hour, day/), (/0,0/))

    CALL modulohour(3, hour, day)
    CALL assert_equal((/hour, day/), (/3,0/))

    CALL modulohour(50, hour, day)
    CALL assert_equal((/hour, day/), (/2,2/))

    CALL modulohour(-1, hour, day)
    CALL assert_equal((/hour, day/), (/23,-1/))

    CALL modulohour(-25, hour, day)
    CALL assert_equal((/hour, day/), (/23,-2/))


  END SUBROUTINE test_modulohour

  SUBROUTINE test_modulomonth()
    INTEGER :: month, year

    CALL modulomonth(1, month, year)
    CALL assert_equal((/month, year/), (/1,0/))

    CALL modulomonth(3, month, year)
    CALL assert_equal((/month, year/), (/3,0/))

    CALL modulomonth(24, month, year)
    CALL assert_equal((/month, year/), (/0,2/))

    CALL modulomonth(25, month, year)
    CALL assert_equal((/month, year/), (/1,2/))

    CALL modulomonth(-1, month, year)
    CALL assert_equal((/month, year/), (/11,-1/))

    CALL modulomonth(-12, month, year)
    CALL assert_equal((/month, year/), (/0,-1/))

    CALL modulomonth(-13, month, year)
    CALL assert_equal((/month, year/), (/11,-2/))

  END SUBROUTINE test_modulomonth

  SUBROUTINE dateutils_test_suite()
    call test_modulomilli()
    CALL test_modulosecond()
    CALL test_modulominute()
    CALL test_modulohour()
    CALL test_modulomonth()

  END SUBROUTINE dateutils_test_suite

END MODULE test_dateutils
