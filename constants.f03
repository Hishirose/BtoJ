!
! Copyright (C) Hishiro T. Hirose
! This file is distributed under the terms of the GNU general public license.
! See http://www.gnu.org/copyleft/gpl.txt .
!

module constants

  implicit none

  character(4), public, parameter :: name_of_program = 'BtoJ.x'
  character(4), public, parameter :: version = '1.0'
  integer,      public, parameter :: sp = selected_real_kind(6, 30)
  integer,      public, parameter :: dp = selected_real_kind(14, 200)
  integer,      public, parameter :: i1b = selected_int_kind(2)
  integer,      public, parameter :: i4b = selected_int_kind(9)
  integer,      public, parameter :: i8b = selected_int_kind(18)
  integer,      public, parameter :: stdin  = 5   ! unit for standard input
  integer,      public, parameter :: stdout = 6   ! unit for standard output
  integer,      public, parameter :: maxlen = 1024

  real(dp),     public, parameter :: pi = 4 * atan(1.0_8)

  public :: datetime_now

contains

  function datetime_now()
    ! Returns a text representing date and time
    ! 'Wed Sep 12 16:32:58 (UTC+0900) 2018'
    ! Weekday is obtained from Zeller's congruence
    ! https://en.wikipedia.org/wiki/Zeller%27s_congruence
    !
    ! It is easy to use UNIX command via `call system('date')`
    implicit none

    integer(i8b)  :: sys_time(8), Zeller, month_shift
    character(35) :: datetime_now
    character(3)  :: weekdays(0 : 6) = &
      & (/'Sat', 'Sun', 'Mon', 'Tue', 'Wed', 'Thu', 'Fri'/)
    character(7)  :: month(12) = &
      & (/'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', &
      &   'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'/)

    call date_and_time(values = sys_time)
    ! 1: year
    ! 2: month
    ! 3: day
    ! 4: shift from UTC in units of minute
    ! 5: hour
    ! 6: minute
    ! 7: second
    ! 8: milisecond

    if (sys_time(2) < 3) then
      month_shift = 12
    else
      month_shift = 0
    end if

    Zeller = sys_time(1) + sys_time(1) / 4 - sys_time(1) / 100 + &
           & sys_time(1) / 400 + sys_time(3) +                   &
           & (13 * (sys_time(2) + month_shift) + 8) / 5

    write(datetime_now, '(a3, x, a3, x, i2.2,                      &
                        & x, i2.2, ":", i2.2, ":", i2.2,           &
                        & x, "(UTC", sp, i3.2, ss, i2.2, ")",      &
                        & x, i4.4)')                               &
      & weekdays(mod(Zeller, 7)), month(sys_time(2)), sys_time(3), & ! Date
      & sys_time(5), sys_time(6), sys_time(7),                     & ! Time
      & sys_time(4) / 60, mod(sys_time(4) + 1440, 60),             & ! Shift
      & sys_time(1)                                                  ! Year

    return
  end function datetime_now

end module constants










