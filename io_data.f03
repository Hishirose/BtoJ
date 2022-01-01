!
! Copyright (C) Hishiro T. Hirose
! This file is distributed under the terms of the GNU general public license.
! See http://www.gnu.org/copyleft/gpl.txt .
!
module io_data
  use constants, only : i4b, sp, dp, maxlen
  implicit none

  integer(i4b),            save, private     :: num_lines
  character(len = maxlen), save, allocatable :: str_stdin(:)
  logical,                 save, allocatable :: cmd_arg_checked(:)

  ! ---------------------------------
  !  System configuration
  ! ---------------------------------
  integer(i4b),            public, save :: verbose
  character(len = maxlen), public, save :: b_map_file_path
  
  ! -----------------------------------------
  !  Calculation configuration
  ! -----------------------------------------
  real(dp),                public, save :: sc_thickness
  real(dp),                public, save :: dz
  real(dp),                public, save :: y_interval
  real(dp),                public, save :: y_interval_interp
  integer(i4b),            public, save :: n_overlap

  ! -----------------------------------------
  !  Output configuration
  ! -----------------------------------------
  character(len = maxlen), public, save :: prefix
  logical,                 public, save :: output_Jx_Jy
  logical,                 public, save :: output_J
  character(len = maxlen), public, save :: output_dir
  character(len = maxlen), public, save :: postfix_Jx
  character(len = maxlen), public, save :: postfix_Jy
  character(len = maxlen), public, save :: postfix_J
  character(len = maxlen), public, save :: separator

  ! -----------------------------------------
  !  Data
  ! -----------------------------------------
  real(dp), allocatable,   public, save :: Bmap(:,:)
  real(dp), allocatable,   public, save :: xpos(:)
  real(dp),                public, save :: x_interval

  real(dp), allocatable,   public, save :: Bmap_interp(:,:)
  real(dp), allocatable,   public, save :: Jx(:,:)
  real(dp), allocatable,   public, save :: Jy(:,:)
  real(dp), allocatable,   public, save :: J_tot(:,:)
  real(dp), allocatable,   public, save :: invG(:, :)

  public  :: read_config
  public  :: load_config_text_into_array
  public  :: parse_cmd_args
  private :: read_keyword_var
  private :: write_keyword_var
  public  :: load_Bmap_file
  public  :: calc_x_interval
  public  :: write_data

contains

  subroutine read_config()
    use constants, only : i4b, dp, pi, maxlen
    use utils,     only : lower_case, report_error, replace_string, expand_path, int2char

    implicit none

    logical                 :: found
    character(len = maxlen) :: keyword
    integer(i4b)            :: i, j
    real(dp)                :: stock_mat(3, 3)

    call load_config_text_into_array()

  ! -----------------------------------------
  !  System configuration
  ! -----------------------------------------
    verbose = 1
    call read_keyword_var('verbose', found, ivar = verbose)
    call write_keyword_var('verbose', found, ivar = verbose)

    b_map_file_path = 'input.dat'
    call read_keyword_var('b_map_file_path', found, cvar = b_map_file_path)
    call expand_path(b_map_file_path)
    call write_keyword_var('b_map_file_path', found, cvar = b_map_file_path)

  ! -----------------------------------------
  !  Calculation configuration
  ! -----------------------------------------
    sc_thickness = 2
    call read_keyword_var('sc_thickness', found, rvar = sc_thickness)
    call write_keyword_var('sc_thickness', found, rvar = sc_thickness)
    sc_thickness = sc_thickness * 0.000001

    dz = 2
    call read_keyword_var('dz', found, rvar = dz)
    call write_keyword_var('dz', found, rvar = dz)
    dz = dz * 0.001d0 ! in units of m

    y_interval = 0.6
    call read_keyword_var('y_interval', found, rvar = y_interval)
    call write_keyword_var('y_interval', found, rvar = y_interval)
    y_interval = y_interval * 0.001d0 ! in units of m

    y_interval_interp = 0.2
    call read_keyword_var('y_interval_interp', found, rvar = y_interval_interp)
    call write_keyword_var('y_interval_interp', found, rvar = y_interval_interp)
    y_interval_interp = y_interval_interp * 0.001d0 ! in units of m

    n_overlap = 2
    call read_keyword_var('n_overlap', found, ivar = n_overlap)
    call write_keyword_var('n_overlap', found, ivar = n_overlap)

  ! -----------------------------------------
  !  Output configuration
  ! -----------------------------------------

    prefix = 'output'
    call read_keyword_var('prefix', found, cvar = prefix)
    call write_keyword_var('prefix', found, cvar = prefix)

    output_Jx_Jy = .false.
    call read_keyword_var('output_jx_jy', found, lvar = output_Jx_Jy)
    call write_keyword_var('output_Jx_Jy', found, lvar = output_Jx_Jy)

    output_J = .false.
    call read_keyword_var('output_j', found, lvar = output_J)
    call write_keyword_var('output_J', found, lvar = output_J)

    output_dir = './'
    call read_keyword_var('output_dir', found, cvar = output_dir)
    call write_keyword_var('output_dir', found, cvar = output_dir)
    call prepare_dir(trim(output_dir))

    postfix_Jx = '_Jx.dat'
    call read_keyword_var('output_postfix_jx', found, cvar = postfix_Jx)
    call write_keyword_var('output_postfix_Jx', found, cvar = postfix_Jx)

    postfix_Jy = '_Jy.dat'
    call read_keyword_var('output_postfix_jy', found, cvar = postfix_Jy)
    call write_keyword_var('output_postfix_Jy', found, cvar = postfix_Jy)

    postfix_J = '_J.dat'
    call read_keyword_var('output_postfix_j', found, cvar = postfix_J)
    call write_keyword_var('output_postfix_J', found, cvar = postfix_J)

    separator = ',\t'
    call read_keyword_var('separator', found, cvar = separator)
    call write_keyword_var('separator', found, cvar = separator)
    separator = replace_string(trim(separator) // ';', '\t', char(9))
    separator = replace_string(trim(separator), '\s', char(32))

  end subroutine read_config
  
  subroutine load_config_text_into_array()
    ! Load configuration data into text array
    ! 1. Check how configuration data comes from (pipe, file or neither)
    ! 2. Check number of lines
    ! 3. Allocate text array with the sufficient size
    ! 4. Load text data and remove comment text
    use constants, only : i4b, maxlen, stdout, stdin
    use utils,     only : report_error, rid_comment, new_unit, close_unit, &
                        & replace_string, execute_shell_command, lower_case, z2h

    implicit none

    integer(i4b)              :: i, tot_num_lines, io_error
    character(len = maxlen)   :: str_line, error_msg
    character(:), allocatable :: stdin_check

    call execute_shell_command('[ -t 0 ] ; echo $?', stdin_check)
    if (trim(stdin_check) == '0') then
      ! There is neither of pipe nor redirect
      allocate(str_stdin(1))
      str_stdin = ' '
      return
    end if
    
    ! Count up number of non-blank lines
    tot_num_lines = 0
    do
      read(stdin, '(a)', iostat = io_error) str_line

      if (io_error > 0) then
        write(error_msg, '(i0, "th line cannot be read properly (1).")') tot_num_lines
        call report_error('load_config_text_into_array', error_msg)
      end if
      if (io_error < 0) exit
      tot_num_lines = tot_num_lines + 1
      if (index(lower_case(str_line), 'end of configuration') > 0) then
        if (index(lower_case(str_line), 'end of configuration') < &
          & min(z2h(index(str_line, '!')) , z2h(index(str_line, '#')))) exit
      end if
    end do

    ! Allocate str_stdin
    if (tot_num_lines == 0) then
      allocate(str_stdin(1))
      str_stdin = ' '
      return
    end if
    allocate(str_stdin(tot_num_lines), stat = io_error)
    if (io_error /= 0) &
      & call report_error('load_config_text_into_array', 'Allocation failed.')
    
    ! Store contents of configuration data (pipe, redirect or file) into str_stdin
    rewind(stdin)
    num_lines = 0
    do i = 1, tot_num_lines
      read(stdin, '(a)', iostat = io_error) str_line
      if (io_error < 0) exit
      if (io_error > 0) then
        write(error_msg, '(i0, "-th line cannot be read properly (2).")') i
        call report_error('load_config_text_into_array', error_msg)
      end if

      num_lines = num_lines + 1
      call rid_comment(str_line)

      ! Convert Tab to space
      str_line = replace_string(str_line, char(9), char(32))

      str_stdin(num_lines) = str_line
    end do

  end subroutine load_config_text_into_array

  subroutine parse_cmd_args(state)
    ! Parse list of command line arguments
    ! Meanings of returned states are as follows:
    !   0 : no special arguments (possibly contains arguments for configuration)
    !   1 : '--help' is specified
    !   2 : '--version' is specified
    !   3 : '--test' is specified
    use constants, only : i4b
    use utils,     only : cmd_arg_str, z2h, report_error

    implicit none

    integer(i4b), intent(out) :: state
    character(:), allocatable :: arg_str
    integer(i4b)              :: i
    logical                   :: skip

    state = 0
    if (command_argument_count() < 1) return ! No command line argument

    allocate(cmd_arg_checked(command_argument_count()))
    cmd_arg_checked = .false.

    skip = .false.
    do i = 1, command_argument_count()
     
      if (skip) then
        ! This argument is a value of the flag in the previous argument
        skip = .false.
        cycle
      end if

      arg_str = cmd_arg_str(i)
      if (arg_str == '-h' .or. arg_str == '--help') then
        state = 1
        return
      end if
      if (arg_str == '-v' .or. arg_str == '--version') then
        state = 2
        return
      end if
      if (arg_str == '-t' .or. arg_str == '--test') then
        state = 3
        return
      end if

      if (len(trim(arg_str)) > 2) then
        if (arg_str(1 : 2) == '--') then
          if (index(arg_str, '=') == 0) then
            ! This argument contains keyword, but do not have value
            if (i == command_argument_count()) call report_error('parse_cmd_args', &
              & 'Keyword "' // trim(arg_str) // '" must has a value.')

            skip = .true.
            cycle
          end if

          if ((index(arg_str, '=') < z2h(index(arg_str, '"'))) .and. &
            & (index(arg_str, '=') < z2h(index(arg_str, "'")))) then
            ! This argument contains both keyword and value
            cycle
          end if
        end if
      end if

    end do
  end subroutine parse_cmd_args

  subroutine read_keyword_var(keyword, found, rvar, ivar, lvar, cvar)
    ! Get variable named 'keyword', and store it into
    ! one of rvar, ivar, lvar, cvar according to its kind
   
    use constants, only : i4b, dp, maxlen
    use utils,     only : report_error, cmd_arg_str, lower_case, simplify_spaces

    implicit none

    character(len = *),           intent(in)    :: keyword
    logical,                      intent(out)   :: found
    real(dp),           optional, intent(inout) :: rvar
    integer(i4b),       optional, intent(inout) :: ivar
    logical,            optional, intent(inout) :: lvar
    character(len = *), optional, intent(inout) :: cvar
    integer(i4b)                                :: i, pos, behind_keyword, io_error
    character(len = maxlen)                     :: param_str, error_msg
    character(:),       allocatable             :: arg_str

    found = .false.

    ! Search keyword in command line arguments
    if (command_argument_count() > 0) then
      do i = 1, command_argument_count()
        if (cmd_arg_checked(i)) cycle

        arg_str = cmd_arg_str(i)

        if (arg_str(1 : 2) /= '--') cycle
        if (len(trim(arg_str)) < 2 + len(trim(keyword))) cycle
        if (lower_case(arg_str(3 : 2 + len(trim(keyword)))) /= trim(keyword)) cycle

        if (len(trim(arg_str)) == 2 + len(trim(keyword))) then
          ! Argument contains only the keyword; value is in the next argument
          ! --<keyword> <value>
          if (i >= command_argument_count()) call report_error('read_keyword_var', &
            & 'Argument --' // trim(keyword) // ' need value.')
          if (cmd_arg_checked(i + 1)) call report_error('read_keyword_var', &
            & 'Argument --' // trim(keyword) // ' need value. If "' // trim(cmd_arg_str(i + 1)) // '" ' // &
            & 'is the value, probably it should be quoted.')
          param_str = cmd_arg_str(i + 1)
          cmd_arg_checked(i : i + 1) = .true.
          found = .true.
          exit
        elseif (arg_str(3 + len(trim(keyword)) : 3 + len(trim(keyword))) == '=') then
          ! Argument contains both keyword and value
          ! --<keyword>=<value>
          if (len(trim(arg_str)) >= 4 + len(trim(keyword))) then
            param_str = arg_str(4 + len(trim(keyword)):)
          else
            ! Value is empty
            param_str = ' '
          end if
          cmd_arg_checked(i) = .true.
          found = .true.
          exit
        end if

        ! keyword in arg_str is longer than the target keyword -> skip
      end do
    end if

    ! Find keyword line in config text
    if (.not. found) then
      do i = 1, num_lines
        pos = index(lower_case(str_stdin(i)), trim(keyword))
        if (pos < 1) cycle
        if (len_trim(str_stdin(i)(: pos - 1)) > 1) cycle ! Contains non-blank character before <keyword>

        behind_keyword = pos + len_trim(keyword)
        if (str_stdin(i)(behind_keyword:behind_keyword) /= ' ' .and. &
          & str_stdin(i)(behind_keyword:behind_keyword) /= '=' .and. &
          & str_stdin(i)(behind_keyword:behind_keyword) /= ':') cycle

        if (found) then
          write(error_msg, '("Second same <keyword> ''", a ,"'' found at ", i0 ,"-th line.")') &
            & trim(keyword), i
          call report_error('read_keyword_var', trim(error_msg))
        end if

        found = .true.
        param_str = str_stdin(i)(behind_keyword:)
      end do
    end if

    if (.not. found) return

    ! Remove '=' or ':'
    param_str = adjustl(param_str)
    if (param_str(1:1) == '=' .or. param_str(1:1) == ':') param_str = adjustl(param_str(2:))

    ! Convert to the appropriate type
    if (present(rvar)) then
      read(param_str, *, iostat = io_error) rvar
      if (io_error /= 0) call report_error('read_keyword_var', &
        & 'Failed to convert "' // trim(param_str) // '" into real.')
    end if
  
    if (present(ivar)) then
      read(param_str, *, iostat = io_error) ivar
      if (io_error /= 0) call report_error('read_keyword_var', &
        & 'Failed to convert "' // trim(param_str) // '" into integer.')
    end if
  
    if (present(lvar)) then
      param_str = lower_case(param_str)
      param_str = param_str(1 : index(param_str, ' ') - 1)
      if (index(param_str, 't') > 0) lvar = .true.
      if (index(param_str, 'f') > 0) lvar = .false.
      if ((index(param_str, 't') == 0 .and. index(param_str, 'f') == 0) .or. &
        & (index(param_str, 't')  > 0 .and. index(param_str, 'f')  > 0)) &
        & call report_error('read_keyword_var', &
                          & 'Failed to convert "' // trim(param_str) // '" into logical.')
    end if
  
    if (present(cvar)) then
      if (param_str(1:1) == '''') then
        cvar = param_str(2:index(param_str(2:), ''''))
      else if(param_str(1:1) == '"') then
        cvar = param_str(2:index(param_str(2:), '"'))
      else
        cvar = trim(param_str)
      end if
    end if
    
  end subroutine read_keyword_var

  subroutine write_keyword_var(keyword, found, rvar, ivar, lvar, cvar, comment)
    ! Output variable named <keyword> depending on 'found' and <verbose>
   
    use constants, only : i4b, dp, maxlen, stdout

    implicit none

    character(len = *),           intent(in) :: keyword
    logical,            optional, intent(in) :: found
    real(dp),           optional, intent(in) :: rvar
    integer(i4b),       optional, intent(in) :: ivar
    logical,            optional, intent(in) :: lvar
    character(len = *), optional, intent(in) :: cvar
    character(len = *), optional, intent(in) :: comment
    character(len = maxlen)                  :: format_str, comment_part

    if (verbose == 0) return
    if (present(found)) then
      if (.not. found .and. verbose == 1) return
    end if

    if (present(comment)) then
      comment_part = ' ' // trim(adjustl(comment))
    else
      comment_part = ' '
    end if

    write(format_str, '("(4x, a, ", i0, "x, ""= "", g0, a)")') max(1, 20 - len(trim(keyword)))
    if (present(rvar)) write(stdout, format_str) trim(keyword), rvar, trim(comment_part)
    if (present(ivar)) write(stdout, format_str) trim(keyword), ivar, trim(comment_part)
    if (present(lvar)) write(stdout, format_str) trim(keyword), lvar, trim(comment_part)
    if (present(cvar)) write(stdout, format_str) trim(keyword), '''' // trim(cvar) // '''', trim(comment_part)
    flush(stdout)
  end subroutine write_keyword_var

  subroutine prepare_dir(dir)
    ! Check existence of dir.
    ! If it doesn't exist, call mkdir.
    use constants, only : i4b
    use utils,     only : report_error, execute_shell_command

    implicit none

    character(len = *), intent(in)  :: dir
    character(:),       allocatable :: dir_check

    call execute_shell_command('[ -e ' // dir // ' ] ; echo $?', dir_check)
    if (trim(dir_check) == '0') return ! Exist

    call execute_shell_command('mkdir ' // dir // ' ; echo $?', dir_check)
    if (trim(dir_check) == '0') return ! Succeeded to create dir

    call report_error('prepare_dir', "Failed to create '" // dir // "'. Output is '" // dir_check // "'.")

  end subroutine prepare_dir

  subroutine load_Bmap_file()
    ! Load TapeStar B(x) file
    
    use constants, only : i4b, dp, maxlen, stdout
    use utils,     only : new_unit, close_unit, report_error, lower_case, int2char, simplify_spaces

    implicit none

    character(len = maxlen)   :: str_line, error_msg
    integer(i4b)              :: unit_Bmap, io_error, ip, ncol, nrow
    real(dp), allocatable     :: buffer(:,:)

    unit_Bmap = new_unit()

    open(unit = unit_Bmap, file = trim(b_map_file_path), &
      & form = 'formatted', action = 'read', status = 'old', &
      & iostat = io_error)
    if (io_error /= 0) call report_error('load_Bmap_file', &
                                       & 'Failed to open ''' // trim(b_map_file_path) // '''.')
    
    read(unit_Bmap, '(a)') str_line
    str_line = simplify_spaces(str_line)

    ! Examine the number of column
    ip = 0
    ncol = 0
    do while(ip < len(str_line))
      ip = ip + 1
      if (str_line(ip:ip) == ' ') cycle
      ncol = ncol + 1
      do while(str_line(ip:ip) /= ' ')
        ip = ip + 1
      end do
    end do

    ! Examine the number of row
    nrow = 1
    do
      read(unit_Bmap, *, iostat = io_error)
      if (io_error < 0) exit
      nrow = nrow + 1
    end do

    allocate(buffer(ncol, nrow))
    rewind(unit_Bmap)

    read(unit_Bmap, *, iostat = io_error) buffer
    if (io_error /= 0) then
      call close_unit(unit_Bmap)
      call report_error('load_Bmap_file', 'Failed to read data file.')
    end if

    call close_unit(unit_Bmap)
    
    allocate(Bmap(ncol - 1, nrow))
    allocate(xpos(nrow))
    Bmap = buffer(2:ncol,:)
    xpos = buffer(1,:)

  end subroutine load_Bmap_file

  subroutine calc_x_interval()

    implicit none

    integer(i4b)          :: i, nx
    real(dp), allocatable :: diff_xpos(:)

    nx = ubound(xpos, 1)
    allocate(diff_xpos(nx - 1))
    forall(i = 1:nx) diff_xpos(i) = xpos(i + 1) - xpos(i)
    x_interval = sum(diff_xpos) / (nx - 1) * 1.0e-3 ! in units of m

  end subroutine calc_x_interval

  subroutine write_data()

    use utils,     only : new_unit, close_unit, report_error, tailless
    implicit none

    character(:), allocatable :: file_path_Jx, file_path_Jy, file_path_J
    integer(i4b)              :: nx, ny
    integer(i4b)              :: unit_Jx, unit_Jy, unit_J, io_error, i, j

    nx = ubound(J_tot, 2); ny = ubound(J_tot, 1)
    file_path_Jx = trim(output_dir) // trim(prefix) // trim(postfix_Jx)
    file_path_Jy = trim(output_dir) // trim(prefix) // trim(postfix_Jy)
    file_path_J  = trim(output_dir) // trim(prefix) // trim(postfix_J)

    if (output_Jx_Jy) then
      unit_Jx = new_unit()
      open(unit = unit_Jx, file = trim(file_path_Jx), &
        & form = 'formatted', status = 'replace', iostat = io_error)
      if (io_error /= 0) call report_error('write_data', &
                         & 'Failed to open ''' // trim(file_path_Jx) // '''.')
      do i = 1, nx
        do j = 1, ny - 1
          write(unit_Jx, '(es21.13, a)', advance="no") Jx(j, i), tailless(separator)
        end do
        write(unit_Jx, '(es21.13)') Jx(ny, i)
      end do
      call close_unit(unit_Jx)

      unit_Jy = new_unit()
      open(unit = unit_Jy, file = trim(file_path_Jy), &
        & form = 'formatted', status = 'replace', iostat = io_error)
      if (io_error /= 0) call report_error('write_data', &
                         & 'Failed to open ''' // trim(file_path_Jy) // '''.')
      do i = 1, nx
        do j = 1, ny - 1
          write(unit_Jy, '(es21.13, a)', advance="no") Jy(j, i), tailless(separator)
        end do
        write(unit_Jy, '(es21.13)') Jy(ny, i)
      end do
      call close_unit(unit_Jy)
    end if

    if (output_J) then
      unit_J = new_unit()
      open(unit = unit_J, file = trim(file_path_J), &
        & form = 'formatted', status = 'replace', iostat = io_error)
      if (io_error /= 0) call report_error('write_data', &
                         & 'Failed to open ''' // trim(file_path_J) // '''.')
      do i = 1, nx
        do j = 1, ny - 1
          write(unit_J, '(es21.13, a)', advance="no") J_tot(j, i), tailless(separator)
        end do
        write(unit_J, '(es21.13)') J_tot(ny, i)
      end do
      call close_unit(unit_J)
    end if

  end subroutine write_data

end module io_data


  
