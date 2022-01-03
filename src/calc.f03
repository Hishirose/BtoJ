!
! Copyright (C) Hishiro T. Hirose
! This file is distributed under the terms of the GNU general public license.
! See http://www.gnu.org/copyleft/gpl.txt .
!
module calc
  use constants, only : i4b, sp, dp, maxlen
  implicit none

  public  :: interpolate_B_along_y
  public  :: obtain_Jx_Jy_by_inverse
  public  :: calc_J_tot
  public  :: calc_invG
  private :: obtain_Mz_by_inverse_mini
  private :: diff_2d
  public  :: calc_B_from_Jx_Jy
  public  :: calc_correction_factor
  private :: evaluate_difference
  private :: renew_factor
  private :: renew_offset

  contains

  subroutine interpolate_B_along_y()
    use io_data, only : Bmap, y_interval, y_interval_interp, Bmap_interp, Jx, Jy, J_tot, outside_points
    use algebra, only : interpolate1D

    integer(i4b)             :: i, nx, ny, ny_interp
    real(dp),    allocatable :: y(:), y_interp(:), B1d(:), B1d_interp(:)

    ny = ubound(Bmap, 1)
    nx = ubound(Bmap, 2)
    allocate(y(ny))
    forall(i = 1 : ny) y(i) = (i - 1) * y_interval
    ny_interp = nint(y(ny) / y_interval_interp) + 1 + outside_points * 2
    allocate(y_interp(ny_interp))
    forall(i = 1 : ny_interp) y_interp(i) = (i - 1 - outside_points) * y_interval_interp

    allocate(Bmap_interp(ny_interp, nx))
    allocate(B1d(ny), B1d_interp(ny_interp))

    !$omp parallel do private(i, B1d, B1d_interp), shared(nx, Bmap, Bmap_interp, y, y_interp)
    do i = 1, nx
      B1d = Bmap(:, i)
      call interpolate1D(y, y_interp, B1d, B1d_interp)
      Bmap_interp(:, i) = B1d_interp(:)
    end do
    !$omp end parallel do

  end subroutine interpolate_B_along_y

  subroutine obtain_Jx_Jy_by_inverse(dx, dy, Bmap)
    use constants, only : stdout
    use io_data, only   : n_overlap, Jx, Jy, sc_thickness, verbose, f_d_factor
    use algebra, only   : fermi_dirac
    implicit none

    real(dp), intent(in)  :: Bmap(:, :)
    real(dp), intent(in)  :: dx, dy
    integer(i4b)          :: ny, nx, n
    real(dp), allocatable :: Mz_block(:, :, :), Mz(:, :)
    integer(i4b)          :: i, j, step
    real(dp)              :: inv_sc_thickness

    ny = ubound(Bmap, 1)
    nx = ubound(Bmap, 2)
    step = ny - n_overlap
    n = ceiling(dble(nx - n_overlap) / step)
    allocate(Mz_block(ny, ny, n), Mz(ny, nx), Jx(ny, nx), Jy(ny, nx))

    if (verbose > 1) write(stdout, '(6x, "allocation is done")'); flush(stdout)
    !$omp parallel do private(i), shared(ny, step, Bmap, Mz_block)
    do i = 1, n - 1
      call obtain_Mz_by_inverse_mini(Bmap(1:ny, (i - 1) * step + 1: (i - 1) * step + ny), Mz_block(:, :, i))
    end do
    !$omp end parallel do
    call obtain_Mz_by_inverse_mini(Bmap(1:ny, nx - ny + 1: nx), Mz_block(:, :, n))
    if (verbose > 1) write(stdout, '(6x, "calculation of Mz is done")'); flush(stdout)

    ! Connect calculated Jx,Jy blocks
    Mz(:, 1:n_overlap) = Mz_block(:, 1:n_overlap, 1)
    Mz(:, nx - ny + 1:nx) = Mz_block(:, :, n)
    if (verbose > 1) write(stdout, '(6x, "connection of blocks (1st step) is done")'); flush(stdout)
    do i = 1, n - 1
      Mz(:, n_overlap + (i - 1) * step + 1 : i * step) = Mz_block(:, n_overlap + 1: ny - n_overlap, i)
    end do
    if (verbose > 1) write(stdout, '(6x, "connection of blocks (2nd step) is done")'); flush(stdout)
    do i = 1, n - 2
      forall(j = 1:n_overlap) Mz(:, i * step + j) = &
          & (Mz_block(:, ny - n_overlap + j, i) * fermi_dirac((j - 1) / dble(n_overlap - 1), beta = f_d_factor) + &
          &  Mz_block(:, j, i + 1) * fermi_dirac((n_overlap - j) / dble(n_overlap - 1), beta = f_d_factor))
    end do
    if (verbose > 1) write(stdout, '(6x, "connection of blocks (3rd step) is done")'); flush(stdout)
    do i = 1, n_overlap
      Mz(:, (n - 1) * step + i) = &
        & (Mz_block(:, ny - n_overlap + i, n - 1) * fermi_dirac((i - 1) / dble(n_overlap - 1), beta = f_d_factor) + &
        &  Mz_block(:, ((n - 1) * step + i) - (nx - ny), n) &
        &* fermi_dirac((n_overlap - i) / dble(n_overlap - 1), beta = f_d_factor))
    end do
    if (verbose > 1) write(stdout, '(6x, "connection of blocks (4th step) is done")'); flush(stdout)

    ! Calculate Jx, Jy from Mz
    forall(i = 1:ny, j = 1:nx) Jx(i, j) =  diff_2d(Mz, 1, dx, dy, j, i)
    forall(i = 1:ny, j = 1:nx) Jy(i, j) = -diff_2d(Mz, 2, dx, dy, j, i)

  end subroutine obtain_Jx_Jy_by_inverse

  subroutine calc_J_tot()
    use io_data, only : J_tot, Jx, Jy
    implicit none
    integer(i4b)          :: ny, nx

    ny = ubound(Jx, 1)
    nx = ubound(Jx, 2)
    allocate(J_tot(ny, nx))

    J_tot(:, :) = sqrt(Jx(:, :)**2 + Jy(:, :)**2)
  end subroutine calc_J_tot

  subroutine calc_invG(dx, dy, Bmap)
    use io_data, only : dz, invG
    use algebra, only : mat_inv
    implicit none

    real(dp), intent(in)  :: Bmap(:, :)
    real(dp), intent(in)  :: dx, dy
    integer(i4b)          :: n
    integer(i4b)          :: i, j, k, l
    real(dp), allocatable :: G(:, :)
    real(dp)              :: r, dxy

    n = ubound(Bmap, 1)
    allocate(G(n * n, n * n), invG(n * n, n * n))

    ! B(y, x)
    dxy = dx * dy
    !$omp parallel do private(i, j, k, l, r), shared(n, dz, dxy, dx, dy, G)
    do i = 1, n
      do j = 0, n - 1
        do k = 1, n
          do l = 0, n - 1
            r = sqrt(((i - k) * dx)**2 + ((j - l) * dy)**2 + dz**2)
            G(i + j * n, k + l * n) = 1.0e-7 * (3.0d0 * dz**2 - r**2) / (r**5) * dxy
          end do
        end do
      end do
    end do
    !$omp end parallel do

    invG = mat_inv(G)
  end subroutine calc_invG

  subroutine obtain_Mz_by_inverse_mini(Bmap_mini, Mz_mini)
    use io_data, only : invG
    implicit none

    real(dp), intent(in)  :: Bmap_mini(:, :)
    real(dp), intent(out) :: Mz_mini(:, :)
    integer(i4b)          :: n, i, j, k, l
    real(dp), allocatable :: Bmap_linear(:), Mz_linear(:)

    n = ubound(Bmap_mini, 1)
    allocate(Bmap_linear(n * n), Mz_linear(n * n))
    Bmap_linear = Reshape(Transpose(Bmap_mini), (/n * n/))

    Mz_linear = matmul(invG, Bmap_linear)
    Mz_mini = Transpose(Reshape(Mz_linear, (/n, n/)))
  end subroutine obtain_Mz_by_inverse_mini

pure function diff_2d(map, d, dx, dy, px, py) result(slope)
    implicit none

    real(dp),     intent(in) :: map(:,:)
    integer(i4b), intent(in) :: d, px, py
    real(dp),     intent(in) :: dx, dy
    real(dp)                 :: slope
    integer(i4b)             :: n


    n = ubound(map, d)
    if (d == 1) then ! y
      if (py > 1 .and. py < n) then
        slope = (map(py + 1, px) - map(py - 1, px)) / (2.0d0 * dy)
      else if (py > 1) then
        slope = (map(py, px) - map(py - 1, px)) / dy
      else
        slope = (map(py + 1, px) - map(py, px)) / dy
      end if
    else if (d == 2) then ! x
      if (px > 1 .and. px < n) then
        slope = (map(py, px + 1) - map(py, px - 1)) / (2.0d0 * dx)
      else if (px > 1) then
        slope = (map(py, px) - map(py, px - 1)) / dx
      else
        slope = (map(py, px + 1) - map(py, px)) / dx
      end if
    end if
  end function diff_2d

  subroutine calc_B_from_Jx_Jy(dx, dy)
    use io_data, only : Jx, Jy, dz, Bmap_sim
    implicit none

    real(dp), intent(in) :: dx, dy
    integer(i4b)         :: nx, ny, i, j, k, l
    real(dp)             :: dxy, dz2
    
    ny = ubound(Jx, 1)
    nx = ubound(Jx, 2)
    allocate(Bmap_sim(ny, nx))
    Bmap_sim = 0.0d0
    dxy = dx * dy
    dz2 = dz * dz

    !$omp parallel do private(i, j, k, l), shared(nx, ny, dx, dy, dxy, dz2, Jx, Jy, Bmap_sim)
    do i = 1, ny
      do j = 1, nx
        do k = 1, ny
          do l = max(1, j - ny), min(nx, j + ny)
            Bmap_sim(i, j) = Bmap_sim(i, j) + &
                           & 1.0e-7 * (Jx(k, l) * ((i - k) * dy) - Jy(k, l) * ((j - l) * dx)) / &
                           & sqrt(((j - l) * dx)**2 + ((i - k) * dy)**2 + dz2)**3 * dxy
          end do
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine calc_B_from_Jx_Jy

  subroutine calc_correction_factor(dx, dy, factor, offset, diff)
    use io_data, only : Jx, Jy, J_tot, Bmap_sim, invG, Bmap_interp, use_correction_factor
    implicit none

    real(dp), intent(in)  :: dx, dy
    real(dp), intent(out) :: factor, offset, diff
    integer(i4b)          :: nx, ny, i
    real(dp), allocatable :: Bmap_diff(:, :)

    ny = ubound(Bmap_interp, 1)
    nx = ubound(Bmap_interp, 2)
    allocate(Bmap_diff(ny, nx))

    ! Give initial values
    Bmap_diff = Bmap_interp(:, :) - Bmap_sim(:, :)
    factor = 1.0d0
    offset = sum(Bmap_diff) / (nx * ny)

    ! Get better values
    if (use_correction_factor) then
      do i = 1, 20
        call renew_factor(Bmap_interp, Bmap_sim, factor, offset)
        call renew_offset(Bmap_interp, Bmap_sim, factor, offset)
      end do
    end if

    Bmap_diff = abs(Bmap_interp(:, :) / (Bmap_sim(:, :) * factor + offset))
    diff = sum(Bmap_diff) / dble(nx * ny)

    ! correct magnitude by using the obtained factor and offset
    Bmap_sim = Bmap_sim * factor + offset
    Jx = Jx * factor
    Jy = Jy * factor
    J_tot = J_tot * factor

  end subroutine calc_correction_factor

  function evaluate_difference(Bmap, Bmap_sim, factor, offset) result(diff)
    implicit none
    
    real(dp), intent(in)  :: Bmap(:, :), Bmap_sim(:, :)
    real(dp), intent(in)  :: factor, offset
    real(dp)              :: diff
    integer(i4b)          :: nx, ny
    real(dp), allocatable :: Bmap_diff(:, :)

    ny = ubound(Bmap, 1)
    nx = ubound(Bmap, 2)
    allocate(Bmap_diff(ny, nx))

    Bmap_diff = abs(Bmap(:, :) - (Bmap_sim(:, :) * factor + offset))
    diff = sum(Bmap_diff)
  end function evaluate_difference

  subroutine renew_factor(Bmap, Bmap_sim, factor, offset)
    implicit none

    real(dp), intent(in)    :: Bmap(:, :), Bmap_sim(:, :)
    real(dp), intent(inout) :: factor
    real(dp), intent(in)    :: offset
    real(dp)                :: diff0, diff_new, step
    integer(i4b)            :: i

    diff0 = evaluate_difference(Bmap, Bmap_sim, factor, offset)
    step = 1.0d0
    do i = 1, 10
      plus: do
        diff_new = evaluate_difference(Bmap, Bmap_sim, factor + step, offset)
        if (diff_new < diff0) then
          diff0 = diff_new
          factor = factor + step
          cycle plus
        end if
        exit plus
      end do plus

      minus: do
        diff_new = evaluate_difference(Bmap, Bmap_sim, factor - step, offset)
        if (diff_new < diff0) then
          diff0 = diff_new
          factor = factor - step
          cycle minus
        end if
        exit minus
      end do minus

      step = step / 10.0d0
    end do

  end subroutine renew_factor

  subroutine renew_offset(Bmap, Bmap_sim, factor, offset)
    implicit none

    real(dp), intent(in)    :: Bmap(:, :), Bmap_sim(:, :)
    real(dp), intent(inout) :: offset
    real(dp), intent(in)    :: factor
    integer(i4b)            :: nx, ny
    real(dp), allocatable   :: Bmap_diff(:, :)

    ny = ubound(Bmap, 1)
    nx = ubound(Bmap, 2)
    allocate(Bmap_diff(ny, nx))

    Bmap_diff = Bmap(:, :) - (Bmap_sim(:, :) * factor + offset)
    offset = offset + sum(Bmap_diff) / dble(nx * ny)

  end subroutine renew_offset

end module calc














