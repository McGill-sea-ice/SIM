
MODULE management

  USE fruit_util, ONLY : strip, to_s

  implicit none


  PRIVATE

  integer, parameter :: MSG_LENGTH = 2000
  integer, parameter :: MSG_STACK_SIZE = 300000

  integer, save :: successful_assert_count = 0
  integer, save :: failed_assert_count = 0
  character (len = MSG_LENGTH), private, dimension (MSG_STACK_SIZE), save :: messageArray
  character (len = MSG_LENGTH), private, save :: msg = '[unit name not set from set_name]: '
  character (len = MSG_LENGTH), private, save :: unit_name  = '_not_set_'
  character (len = MSG_LENGTH), private, save :: bench_name  = '_not_set_'
  character (len = MSG_LENGTH), private, dimension (1000), save :: benchmessages
  integer, save :: messageIndex = 1


  integer, save :: successful_case_count = 0
  integer, save :: failed_case_count = 0
  integer, save :: testCaseIndex = 1
  INTEGER, SAVE :: bench_count = 0
  LOGICAL, SAVE :: last_passed = .FALSE.
  
  REAL, SAVE :: start_time, end_time, bench_start_time, bench_end_time

  INTERFACE init_fruit
     module procedure init_fruit_
  end interface

  INTERFACE fruit_summary
     module procedure fruit_summary_
  end interface

  interface get_last_message
     module procedure get_last_message_
  end interface

  interface is_last_passed
     module procedure is_last_passed_
  end interface

  interface add_success
     module procedure add_success_
  end interface

  interface add_fail
     module procedure add_fail_
     module procedure add_fail_unit_
  end interface

  INTERFACE set_unit_name
     module procedure set_unit_name_
  end interface

  interface get_unit_name
     module procedure get_unit_name_
  end interface

  interface get_total_count
     module procedure get_total_count_
  end interface

  interface get_failed_count
     module procedure get_failed_count_
  end interface

  INTERFACE success_assert_action
     module procedure success_assert_action_
  end interface

  INTERFACE failed_assert_action
     module procedure failed_assert_action_
  end interface

  interface is_all_successful
     module procedure is_all_successful_
  end interface

  PUBLIC :: set_unit_name, success_assert_action, failed_assert_action, init_fruit
  PUBLIC :: init_bench, end_bench
  PUBLIC :: fruit_summary


CONTAINS


  subroutine init_fruit_
    successful_assert_count = 0
    failed_assert_count = 0
    messageIndex = 1
    CALL CPU_TIME(start_time)
  end subroutine init_fruit_

  SUBROUTINE init_bench(name)
    CHARACTER(*), intent(in) :: name
    bench_name = name
    call CPU_TIME(bench_start_time)
    bench_count = bench_count + 1
  END SUBROUTINE init_bench

  SUBROUTINE end_bench
    call CPU_TIME(bench_end_time)
    WRITE(benchmessages(bench_count), "(A, A, F0.4, A)") TRIM(bench_name), ': ', bench_end_time - bench_start_time, 's.'
  END SUBROUTINE end_bench

  subroutine fruit_summary_
    INTEGER :: i, total_asserts
    call cpu_time(end_time)
    WRITE (*,*) '\n', REPEAT("-", 80)
    WRITE(*,*) 'Ran ', TRIM(to_s(total_count())), &
         ' tests in ', TRIM(to_s(end_time-start_time)), ' seconds.'

    if ( messageIndex > 1) then
       WRITE(*,*)
       write (*,*) 'Failed tests:'

       do i = 1, messageIndex - 1
          write (*,"(A)") trim(strip(messageArray(i)))
       end do

!       write (*,*) '  -- end of failed assertion messages.'
       write (*,*)
    else
       write (*,*) 'OK\n'
    end if

    if (total_count() /= 0) then

       WRITE (*,*) 'Successful    :   ', successful_assert_count, '/ ', trim(to_s(total_count()))
       write (*,*) 'Failed        :   ', failed_assert_count, '/ ', trim(to_s(total_count()))
       write (*,'(" Success rate:   ",f6.2,"%")')  real(successful_assert_count) * 100.0 / &
            real (total_count())
       write (*, *)
    end if

    IF (bench_count > 0) THEN
       WRITE(*,*)
       WRITE(*,*) "Benchmark results"
       WRITE(*,*) "-----------------"
       DO i=1, bench_count
          WRITE(*,*) trim(strip(benchmessages(i)))
       ENDDO
    ENDIF
    
  end subroutine fruit_summary_

  subroutine add_success_
    call success_assert_action_
  end subroutine add_success_

  subroutine add_fail_ (message)
    character (*), intent (in), optional :: message
    call failed_assert_action_('none', 'none', message)
  end subroutine add_fail_

  subroutine add_fail_unit_ (unitName, message)
    character (*), intent (in) :: unitName
    character (*), intent (in) :: message

    call add_fail_ ("[in " //  unitName // "(fail)]: " //  message)
  end subroutine add_fail_unit_

  subroutine is_all_successful_(result)
    logical, intent(out) :: result
    result= (failed_assert_count .eq. 0 )
  end subroutine is_all_successful_

  subroutine success_mark_
    write(*,"(A1)",ADVANCE='NO') '.'
    IF (MOD(total_count(), 80)==0) WRITE(*, *) 

  end subroutine success_mark_

  subroutine failed_mark_
    write(*,"(A1)",ADVANCE='NO') 'F'
    IF (MOD(total_count(), 80)==0) WRITE(*, *) 

  end subroutine failed_mark_

  subroutine increase_message_stack_
    if (messageIndex > MSG_STACK_SIZE) then
       write (*, *) "Too many errors to put into stack"
       call fruit_summary ()
       stop 1
    end if

    messageArray (messageIndex) = msg
    messageIndex = messageIndex + 1
  end subroutine increase_message_stack_

  function get_last_message_()
    character(len=MSG_LENGTH) :: get_last_message_
    if (messageIndex > 1) then
       get_last_message_ = strip(messageArray(messageIndex-1))
    else
       get_last_message_ = ''
    end if
  end function get_last_message_

  FUNCTION total_count()
    INTEGER :: total_count
    total_count = successful_assert_count + failed_assert_count
  END FUNCTION total_count
  
  SUBROUTINE get_total_count_(count) 
    integer, intent(out) :: count

    count = successful_assert_count + failed_assert_count
  end subroutine get_total_count_

  subroutine get_failed_count_ (count)
    integer, intent(out) :: count
    count = failed_assert_count
  end subroutine get_failed_count_

  subroutine success_assert_action_
    successful_assert_count = successful_assert_count + 1
    last_passed = .true.
    call success_mark_  
  end subroutine success_assert_action_

  subroutine failed_assert_action_ (expected, got, message)
    character(*), intent(in) :: expected, got
    character(*), intent(in), optional :: message

    call make_error_msg_ (expected, got, message)
    call increase_message_stack_
    failed_assert_count = failed_assert_count + 1
    last_passed = .false.
    call failed_mark_
  end subroutine failed_assert_action_

  subroutine set_unit_name_(value)
    character(*), intent(in) :: value
    unit_name = strip(value)
  end subroutine set_unit_name_

  subroutine get_unit_name_(value)
    character(*), intent(out) :: value
    value = strip(unit_name)
  end subroutine get_unit_name_

  subroutine make_error_msg_ (var2, var1, message)
    character(*), intent(in) :: var1, var2
    character(*), intent(in), optional :: message
    if (present(message)) then
       msg = '\n ### ' // trim(strip(unit_name)) // ' ###\n\tExpected: ' // trim(strip(var1)) &
            // '\n\tGot: ' // TRIM(strip(var2)) // '\n\t' // message
    else
       msg = '[' // trim(strip(unit_name)) // ']: Expected [' // trim(strip(var1)) // '], Got [' // trim(strip(var2)) // ']' 
    endif
  end subroutine make_error_msg_

  SUBROUTINE make_bench_msg(unit, time)
    CHARACTER(*), INTENT(in) :: unit
    INTEGER, intent(in) :: time

  END SUBROUTINE make_bench_msg

  function is_last_passed_()
    logical:: is_last_passed_
    is_last_passed_ = last_passed 
  end function is_last_passed_



END MODULE management
