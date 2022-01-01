!
! Copyright (C) Hishiro T. Hirose
! This file is distributed under the terms of the GNU general public license.
! See http://www.gnu.org/copyleft/gpl.txt .
!

module utils

  implicit none

  logical,                 save, private     :: open_unit_list(10:99) = .false.

  public  :: new_unit
  public  :: close_unit
  public  :: execute_shell_command
  public  :: report_error
  public  :: replace_string
  public  :: rid_comment
  public  :: z2h
  public  :: lower_case
  public  :: simplify_spaces
  public  :: expand_path
  public  :: int2char
  public  :: cmd_arg_str
  public  :: secs2str
  public  :: tailless
  public  :: unique_file_name


  ! https://gcc.gnu.org/onlinedocs/gfortran/FLUSH.html
  ! Declare the interface for POSIX fsync function
  interface
    function fsync(fd) bind(c,name="fsync")
    use iso_c_binding, only: c_int
      integer(c_int), value :: fd
      integer(c_int)        :: fsync
    end function fsync
  end interface

contains

  function new_unit()
    ! Returns unused unit number
    ! Unit starts from 10 to 99 (up to 90)
    use constants, only : i4b

    implicit none

    integer(i4b) :: new_unit, i

    do i = 10, 99
      if (.not. open_unit_list(i)) then
        open_unit_list(i) = .true.
        new_unit = i
        return
      end if
    end do
    
    call report_error(routine = 'new_unit', &
         & error_msg = "All units are opened. Something is wrong!")
    return
  end function new_unit

  subroutine close_unit(u)
    ! Close unit, and set .false. in open_unit_list
    use constants, only : i4b

    implicit none

    integer(i4b) :: u
    integer(i4b) :: ret
    
    ! Flush dirty data to storage device
    flush(u)
    ret = fsync(fnum(u))
    if (ret /= 0) call report_error(routine = 'close_unit', &
                  & error_msg = 'Failed to call fsync')
    
    close(u)

    if (9 < u .and. u < 100) open_unit_list(u) = .false.
  end subroutine close_unit

  subroutine execute_shell_command(command, result_str, status)
    ! Execute command in shell
    ! Output of command is stored in result_str
    ! Success or failure is stored in status
    use constants, only : i4b, maxlen, name_of_program

    implicit none

    character(len = *), intent(in)               :: command
    character(:),       intent(out), allocatable :: result_str
    integer(i4b),       intent(out), optional    :: status
    integer(i4b)                                 :: unit_tmp, io_error, status_
    character(:),                    allocatable :: command_, path_tmp_file
    character(len = maxlen)                      :: str_line
 
    result_str = ''
    path_tmp_file = unique_file_name('.', trim(name_of_program), '.tmp')
    command_ = trim(command) // ' > ' // path_tmp_file

    call system(trim(command_), status_)
    if (present(status)) status = status_

    if (status_ /= 0 .and. .not. present(status)) &
       & call report_error('execute_shell_command', 'Execution of "' // trim(command_) // '" failed.')

    unit_tmp = new_unit()
    open(unit = unit_tmp, file = trim(path_tmp_file), &
      & form = 'formatted', action = 'read', status = 'old', &
      & iostat = io_error)
    if (io_error /= 0) call report_error('execute_shell_command', &
                                       & 'Failed to open temporary file "' // path_tmp_file // '".')

    do
      read(unit_tmp, '(a)', iostat = io_error) str_line
      if (io_error < 0) exit
      if (io_error > 0) call report_error('execute_shell_command', &
                                        & 'Temporary file cannot be read properly.')
      result_str = result_str // trim(str_line) // char(10)
    end do

    result_str = tailless(result_str, 1)

    call close_unit(unit_tmp)
   
    command_ = 'rm ' // path_tmp_file
    call system(command_, status_)
    if (status_ /= 0) call report_error('execute_shell_command', &
                                     & 'Failed to delete temporary file ".' // name_of_program // '".')

  end subroutine execute_shell_command

  subroutine report_error(routine, error_msg)
    ! Output error message in stdout
    use constants, only : i4b, stdout, name_of_program

    implicit none

    character(len = *), intent(in) :: routine
    character(len = *), intent(in) :: error_msg

    write(stdout, '(/, a, ": Error in ", a, ": ", a)') &
      & trim(name_of_program), trim(routine), trim(error_msg)
   
    call exit(1)
    stop 1
  end subroutine report_error

  function replace_string(original_str, before_str, after_str)
    ! Replace 'before_str' in 'original_str' with 'after_str'
    
    use constants, only : i4b, maxlen

    implicit none

    character(len = *), intent(in) :: original_str, before_str, after_str
    character(len = maxlen)        :: replace_string

    integer(i4b) :: pos, shift

    replace_string = original_str
    if (len(before_str) < 1) return

    shift = 1
    pos = index(replace_string(shift:), before_str)
    do while(pos > 0)
      shift = shift + pos - 1
      replace_string(shift:) = after_str//replace_string(shift + len(before_str):)    
      shift = shift + len(after_str)
      pos = index(replace_string(shift:), before_str)
    end do

    return
  end function replace_string

  subroutine rid_comment(str_line)
    ! Remove comment from text
    use constants, only : i4b

    implicit none

    character(*), intent(inout) :: str_line
    integer(i4b)                :: index_comment, ic
    integer(i4b)                :: ie, is, iq1, iq2, id1, id2

    ie = index(str_line, '!')
    is = index(str_line, '#')
    if (ie == 0 .and. is == 0) return ! No comment

    iq1 = index(str_line, "'")
    id1 = index(str_line, '"')
    if (min(z2h(ie), z2h(is)) < min(z2h(iq1), z2h(id1))) then
      ! No quotation mark appears before the start of comment
      str_line = str_line(1 : min(z2h(is), z2h(ie)) - 1)
      return
    end if

    ! There is at least 1 quotation mark before comment text
    ic = min(z2h(iq1), z2h(id1))
    do
      if (z2h(iq1) < z2h(id1)) then
        ! "'" comes before '"'
        iq2 = index(str_line(iq1 + 1:), "'")
        if (iq2 == 0) return ! No closing "'": Assume entire line as input
        ic = iq1 + iq2
      else ! z2h(iq1) > z2h(id1)
        ! '"' comes before "'"
        id2 = index(str_line(id1 + 1:), '"')
        if (id2 == 0) return ! No closing '"': Assume entire line as input
        ic = id1 + id2
      end if
      ! Redefine index
      iq1 = ic + index(str_line(ic + 1:), "'")
      id1 = ic + index(str_line(ic + 1:), '"')
      ie  = ic + index(str_line(ic + 1:), '!')
      is  = ic + index(str_line(ic + 1:), '#')

      if (ie == ic .and. is == ic) return ! No comment
      if (min(z2h(ie  - ic), z2h(is  - ic)) < &
        & min(z2h(iq1 - ic), z2h(id1 - ic))) then
        ! No more quotation marks, or
        ! next quotation mark comes after the start of comment
        str_line = str_line(1 : ic +  min(z2h(is - ic), z2h(ie - ic)) - 1)
        return
      end if
      ! Next quotation mark comes before comment text: Go to next cycle
    end do
  end subroutine rid_comment

  function z2h(input)
    ! If input integer is 0, return huge value
    ! This is useful when comparing which character comes first
    ! ex: index('test', 't') -> 1
    !     index('test', 'e') -> 2
    !     index('test', 'a') -> 0
    !     Now, which character comes first?
    !     We expect 't' < 'e' < 'a', but index returns 1 < 2 < 0 for these characters
    !     So, we convert 0 to huge to make comparison simple
    use constants, only : i4b

    implicit none

    integer(i4b), intent(in) :: input
    integer(i4b)             :: z2h

    if (input == 0) then
      z2h = huge(input)
    else
      z2h = input
    end if
  end function z2h

  function lower_case(original_str)
    ! Convert to lowercase characters
    use constants, only : i4b

    implicit none

    character(len = *), intent(in)     :: original_str
    character(len = len(original_str)) :: lower_case

    integer(i4b)                         :: i

    lower_case = original_str
    do i = 0, 25
      lower_case = replace_string(lower_case, char(ichar('A') + i), char(ichar('a') + i))
    end do
    
    return
  end function lower_case


  function simplify_spaces(original_str)
    ! Convert tab & multiple spaces to single space
    use constants, only : i4b

    implicit none

    character(len = *), intent(in)     :: original_str
    character(len = len(original_str)) :: simplify_spaces

    integer(i4b)                         :: pos

    ! Convert Tab to space
    simplify_spaces = replace_string(original_str, char(9), ' ')

    ! Remove multiplied spaces
    simplify_spaces = adjustl(simplify_spaces)
    pos = index(trim(simplify_spaces), '  ')
    do while(pos > 0)
      simplify_spaces = replace_string(trim(simplify_spaces), '  ', ' ')
      pos = index(trim(simplify_spaces), '  ')
    end do
    
    return
  end function simplify_spaces

  subroutine expand_path(path)
    ! Replace '~' in path with $HOME
    use constants, only : maxlen

    implicit none

    character(len = maxlen), intent(inout) :: path
    character(len = maxlen)                :: home

    if (path(1 : 2) /= '~/') return

    call get_environment_variable('HOME', home)
    if (len(trim(home)) < 1) &
      & call report_error('expand_path', &
             & 'Environment variable $HOME is empty. ' // &
             & 'Cannot convert "' // path // '" to real path.')
    path = trim(home) // trim(path(2:))
    return
  end subroutine

  function int2char(ival) result(cval)
    use constants, only : i4b, maxlen

    implicit none

    integer(i4b), intent(in)  :: ival
    character(len = maxlen)   :: cval_dummy
    character(:), allocatable :: cval
   
    write(cval_dummy, '(i0)') ival
    cval = trim(cval_dummy)
  end function int2char

  function cmd_arg_str(index_arg) result(arg_str)
    ! Returns string of 'index_arg'-th command line argument
    use constants, only : i4b, maxlen

    implicit none

    integer(i4b), intent(in)  :: index_arg
    character(:), allocatable :: arg_str
    integer(i4b)              :: arg_len, io_error
    character(:), allocatable :: error_msg

    if (index_arg > command_argument_count()) then
      arg_str = ''
      return
    end if

    call get_command_argument(index_arg, length = arg_len, status = io_error)
    if (io_error /= 0) then
      write(error_msg, '("Cannot get the length of ", i0, "-th argument.")') index_arg
      call report_error('cmd_arg_str', error_msg)
    end if

    allocate(character(arg_len) :: arg_str)
    call get_command_argument(index_arg, arg_str, status = io_error)
    if (io_error /= 0) then
      write(error_msg, '("Cannot read ", i0, "-th argument.")') index_arg
      call report_error('cmd_arg_str', error_msg)
    end if
    return
  end function cmd_arg_str

  function secs2str(sec) result(str)
    ! Convert sec to appropriate string

    use constants, only : dp, i4b, maxlen

    implicit none

    real(dp)           :: sec
    integer(i4b)       :: minute, hour, day
    character(len = maxlen) :: str

    if (sec < 60.0_8) then
      write(str, '(f6.3, "s")') sec
    else if (sec < 3.6e+3_8) then
      minute = floor(sec / 60.0_8)
      sec = mod(sec, 60.0_8)
      write(str, '(i2, "m", i2, "s")') minute, int(sec)
    else if (sec < 8.64e+4_8) then
      hour = floor(sec / 3.6e+3_8)
      minute = floor(sec / 60.0_8) - hour * 60
      sec = mod(sec, 60.0_8)
      write(str, '(i2, "h", i2, "m", i2, "s")') hour, minute, int(sec)
    else
      day = floor(sec / 8.64e+8_8)
      hour = floor(sec / 3.6e+3_8) - day * 24
      minute = floor(sec / 60.0_8) - hour * 60 - day * 1440
      sec = mod(sec, 60.0_8)
      write(str, '(i0, "d", i2, "h", i2, "m", i2, "s")') day, hour, minute, int(sec)
    end if
  end function secs2str

  function tailless(str, rm)
    ! Remove last characters for length specified in 'rm'
    ! This is used to preserve the intended spaces at the end from being trimmed
    ! ex: tailless('     ;', 1) -> '     '
    use constants, only : i4b, maxlen

    implicit none

    character(len = *), intent(in)           :: str
    integer(i4b)      , intent(in), optional :: rm
    character(:)      , allocatable          :: tailless
    integer(i4b)                             :: rm_, whole_len

    if (present(rm)) then
      rm_ = rm
    else
      rm_ = 1
    end if

    whole_len = len(trim(str))
    rm_ = min(rm_, whole_len - 1)
    tailless = str(1 : whole_len - rm_)
    return
  end function tailless

  function unique_file_name(prefix, seed, postfix) result(file_name)
    ! Check existence of the file whose name is supplied by prefix, seed and postfix.
    ! If it exists, append (or change) number at the end of seed until it becomes unique name.
    use constants, only : i4b, maxlen

    character(len = *), intent(in)  :: prefix, seed, postfix
    character(len = :), allocatable :: file_name
    integer(i4b)                    :: i, status

    file_name = prefix // seed // postfix
    status = access(file_name, ' ')
    if (status > 0) return

    do i = 0, 9999
      file_name = prefix // seed // int2char(i) // postfix
      status = access(file_name, ' ')
      if (status > 0) return
    end do

    call report_error('unique_file_name', &
         & 'Failed to produce unique file name.')

  end function unique_file_name

end module utils











