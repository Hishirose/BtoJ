!
! Copyright (C) Hishiro T. Hirose
! This file is distributed under the terms of the GNU general public license.
! See http://www.gnu.org/copyleft/gpl.txt .
!

program BtoJ

  use omp_lib
  use constants,   only : dp, i4b, maxlen, stdout, datetime_now, name_of_program, version
  use utils,       only : report_error, secs2str, execute_shell_command, int2char
  use io_data,     only : parse_cmd_args, read_config, load_Bmap_file, calc_x_interval, write_data, x_interval
  use calc,        only : interpolate_B_along_y, obtain_Jx_Jy_by_inverse, calc_J_tot, calc_invG

  implicit none

  integer(i4b)              :: io_error, is, calc_step, state_cmd_args, stdin_type
  real(dp)                  :: cpu_time0, cpu_time1, cpu_time2
  character(:), allocatable :: pwd, stdin_check

  ! Check command line arguments
  call parse_cmd_args(state_cmd_args)
  if (state_cmd_args == 1) call print_help()
  if (state_cmd_args == 2) call print_version()
  if (state_cmd_args == 3) call test()

  write(stdout, '(a)') '# ' // trim(name_of_program) // ' version ' // trim(version) // ' - standard output'
  write(stdout, '(4x, a)') 'Executed at ' // datetime_now()
  call execute_shell_command('pwd', pwd)
  write(stdout, '(4x, a)') 'Working dir = ' // '''' // pwd // ''''

  write(stdout, *)
  write(stdout, '(a)') '# Initialization'
  write(stdout, *)
  stdin_type = 0
  ! Check if stdin is pipe(0) or not(1)
  call execute_shell_command('[ -p /dev/stdin ] ; echo $?', stdin_check)
  if (trim(stdin_check) == '0') stdin_type = stdin_type + 1
  ! Check if stdin is file(0) or not(1)
  call execute_shell_command('[ -f /dev/stdin ] ; echo $?', stdin_check)
  if (trim(stdin_check) == '0') stdin_type = stdin_type + 2
  ! Check if stdin is terminal(0) or not(1)
  call execute_shell_command('[ -t 0 ] ; echo $?', stdin_check)
  if (trim(stdin_check) == '0') stdin_type = stdin_type + 4

  if (stdin_type == 1) then
    write(stdout, '("# Loading configuration from stdin(pipe)")')
  elseif (stdin_type == 2) then
    write(stdout, '("# Loading configuration from stdin(redirect)")')
  elseif (stdin_type == 4) then
    write(stdout, '("# No input is detected; set default settings")')
  else
    call report_error('main', &
         & 'Configuration data seems to be passed in more than one way (' // int2char(stdin_type) // ').')
  endif
  call cpu_time(cpu_time0)
  call read_config()

  call load_Bmap_file()

  call calc_x_interval()
  write(stdout, '(4x, "x_interval", 10x, "= ", g0)') x_interval

  call cpu_time(cpu_time1)
  call interpolate_B_along_y()
  call cpu_time(cpu_time2)
  write(stdout, '(4x, "interpolation is done (", a, ")")') trim(secs2str(cpu_time2 - cpu_time1))

  call cpu_time(cpu_time1)
  call calc_invG()
  call cpu_time(cpu_time2)
  write(stdout, '(4x, "calculation of invG is done (", a, ")")') trim(secs2str(cpu_time2 - cpu_time1))

  call cpu_time(cpu_time1)
  call obtain_Jx_Jy_by_inverse()
  call cpu_time(cpu_time2)
  write(stdout, '(4x, "calculation of Jx, Jy is done (", a, ")")') trim(secs2str(cpu_time2 - cpu_time1))

  call calc_J_tot()

  call write_data()

  call cpu_time(cpu_time2)
  write(stdout, '(4x, "done (", a, ")")') trim(secs2str(cpu_time2 - cpu_time0))



  stop

contains

  subroutine print_help()
    ! This subroutine is called when --help option is specified
    ! Output usage to stdout
    use constants, only : stdout, name_of_program

    implicit none

    write(stdout, '(a)') 'Usage: ' // name_of_program // &
      & ' [--help] [--version] [--test] [--<keyword>=<value> | --<keyword> <value>]... [<file>]'

    call exit(0)
    stop 0
  end subroutine print_help

  subroutine print_version()
    ! This subroutine is called when --version option is specified
    ! Output usage to stdout
    use constants, only : stdout, name_of_program, version

    implicit none

    write(stdout, '(a)') name_of_program // ' ' // version

    call exit(0)
    stop 0
  end subroutine print_version

  subroutine test()
    ! This subroutine is called when --test option is specified
    use constants, only : i4b, dp

    write(stdout, '(a)') 'test is called'
    call exit(0)
    stop 0
  end subroutine test

end program BtoJ









