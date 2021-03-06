!
! Copyright (C) Hishiro T. Hirose
! This file is distributed under the terms of the GNU general public license.
! See http://www.gnu.org/copyleft/gpl.txt .
!
! A part of this program uses the algorithm described in
! W. Xing, B. Heinrich, H. Zhou, A. A. Fife, and A. R. Cragg, Journal of Applied Physics 76, 4244 (1994)
! https://doi.org/10.1063/1.357308
!

program BtoJ

  use omp_lib
  use constants,   only : dp, i4b, maxlen, stdout, datetime_now, name_of_program, version
  use utils,       only : report_error, secs2str, execute_shell_command, int2char
  use io_data,     only : parse_cmd_args, read_config, load_Bmap_file, calc_x_interval, write_data, x_interval, &
                        & y_interval_interp, y_interval, Bmap, Bmap_interp, factor, offset, Bmap_sim, J_tot, &
                        & output_Bsim
  use calc,        only : interpolate_B_along_y, obtain_Jx_Jy_by_inverse, calc_J_tot, calc_invG, calc_B_from_Jx_Jy, &
                        & calc_correction_factor

  implicit none

  integer(i4b)              :: io_error, is, calc_step, state_cmd_args, stdin_type
  real(dp)                  :: cpu_time0, cpu_time1, cpu_time2, mean_difference_factor
  character(:), allocatable :: pwd, stdin_check

  ! Check command line arguments
  call parse_cmd_args(state_cmd_args)
  if (state_cmd_args == 1) call print_help()
  if (state_cmd_args == 2) call print_version()
  if (state_cmd_args == 3) call test()

  write(stdout, '(a)') '# ' // trim(name_of_program) // ' version ' // trim(version) // ' - standard output'
  write(stdout, '("#", 3x, a)') 'Executed at ' // datetime_now()
  call execute_shell_command('pwd', pwd)
  write(stdout, '("#", 3x, a)') 'Working dir = ' // '''' // pwd // ''''

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
  flush(stdout)
  call read_config()
  write(stdout, '(4x, "End of configuration")')
  write(stdout, *)
  flush(stdout)

  call load_Bmap_file()

  write(stdout, '("# Calculation")')
  flush(stdout)
  call calc_x_interval()
  write(stdout, '(4x, "x_interval", 10x, "= ", g0, " ! mm")') x_interval * 1.0e3
  flush(stdout)

  call cpu_time(cpu_time1)
  call interpolate_B_along_y()
  call cpu_time(cpu_time2)
  write(stdout, '(4x, "interpolation is done (", a, ")")') trim(secs2str(cpu_time2 - cpu_time1))
  flush(stdout)

  call cpu_time(cpu_time1)
  call calc_invG(x_interval, y_interval_interp, Bmap_interp)
  call cpu_time(cpu_time2)
  write(stdout, '(4x, "calculation of invG is done (", a, ")")') trim(secs2str(cpu_time2 - cpu_time1))
  flush(stdout)

  call cpu_time(cpu_time1)
  call obtain_Jx_Jy_by_inverse(x_interval, y_interval_interp, Bmap_interp)
  call cpu_time(cpu_time2)
  write(stdout, '(4x, "calculation of Jx, Jy is done (", a, ")")') trim(secs2str(cpu_time2 - cpu_time1))
  flush(stdout)

  call cpu_time(cpu_time1)
  call calc_J_tot()
  call cpu_time(cpu_time2)
  write(stdout, '(4x, "calculation of Jtot is done (", a, ")")') trim(secs2str(cpu_time2 - cpu_time1))
  flush(stdout)

  call cpu_time(cpu_time1)
  call calc_B_from_Jx_Jy(x_interval, y_interval_interp)
  call cpu_time(cpu_time2)
  write(stdout, '(4x, "calculation of Bz_simulated is done (", a, ")")') trim(secs2str(cpu_time2 - cpu_time1))
  flush(stdout)

  call cpu_time(cpu_time1)
  call calc_correction_factor(x_interval, y_interval_interp, factor, offset, mean_difference_factor)
  call cpu_time(cpu_time2)
  write(stdout, '(4x, "calculation of correction factor is done (", a, ")")') trim(secs2str(cpu_time2 - cpu_time1))
  write(stdout, '(4x, "estimated_ext_field", 1x, "= ", g0, " ! T")') offset
  write(stdout, '(4x, "Bz/Bz_simulated", 5x, "= ", g0)') mean_difference_factor
  write(stdout, '(4x, "correction_factor", 3x, "= ", g0)') factor
  flush(stdout)

  call cpu_time(cpu_time1)
  call write_data()
  call cpu_time(cpu_time2)
  write(stdout, '(4x, "output to file(s) is done (", a, ")")') trim(secs2str(cpu_time2 - cpu_time1))
  flush(stdout)

  call cpu_time(cpu_time2)
  write(stdout, '(4x, "done (", a, ")")') trim(secs2str(cpu_time2 - cpu_time0))
  flush(stdout)



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









