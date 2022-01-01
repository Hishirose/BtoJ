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

  contains

  subroutine interpolate_B_along_y()
    use io_data, only : Bmap, y_interval, y_interval_interp, Bmap_interp, Jx, Jy, J_tot
    use algebra, only : interpolate1D

    integer(i4b)             :: i, nx, ny, ny_interp
    real(dp),    allocatable :: y(:), y_interp(:), B1d(:), B1d_interp(:)

    ny = ubound(Bmap, 1)
    nx = ubound(Bmap, 2)
    allocate(y(ny))
    forall(i = 1 : ny) y(i) = (i - 0) * y_interval
    ny_interp = nint(y(ny) / y_interval_interp) + 1
    allocate(y_interp(ny_interp))
    forall(i = 1 : ny_interp) y_interp(i) = (i - 0) * y_interval_interp

    allocate(Bmap_interp(ny_interp, nx), Jx(ny_interp, nx), Jy(ny_interp, nx), J_tot(ny_interp, nx))
    allocate(B1d(ny), B1d_interp(ny_interp))

    !$omp parallel do private(i, B1d, B1d_interp), shared(nx, Bmap, Bmap_interp, y, y_interp)
    do i = 1, nx
      B1d = Bmap(:, i)
      call interpolate1D(y, y_interp, B1d, B1d_interp)
      Bmap_interp(:, i) = B1d_interp(:)
    end do
    !$omp end parallel do

  end subroutine interpolate_B_along_y

  subroutine obtain_Jx_Jy_by_inverse()
    use constants, only : stdout
    use io_data, only   : n_overlap, Bmap_interp, Jx, Jy, sc_thickness, y_interval_interp, &
                        & x_interval, verbose, external_field, f_d_factor
    use algebra, only   : fermi_dirac
    implicit none

    integer(i4b)          :: ny, nx, n
    real(dp), allocatable :: Mz_block(:, :, :), Mz(:, :)
    integer(i4b)          :: i, j, step
    real(dp)              :: inv_sc_thickness

    ny = ubound(Bmap_interp, 1)
    nx = ubound(Bmap_interp, 2)
    step = ny - n_overlap
    n = ceiling(dble(nx - n_overlap) / step)
    allocate(Mz_block(ny, ny, n), Mz(ny, nx))

    Bmap_interp = Bmap_interp - external_field

    if (verbose > 1) write(stdout, '(4x, "allocation is done")'); flush(stdout)
    !$omp parallel do private(i), shared(ny, step, Bmap_interp, Mz_block)
    do i = 1, n - 1
      call obtain_Mz_by_inverse_mini(Bmap_interp(1:ny, (i - 1) * step + 1: (i - 1) * step + ny), Mz_block(:, :, i))
    end do
    !$omp end parallel do
    call obtain_Mz_by_inverse_mini(Bmap_interp(1:ny, nx - ny + 1: nx), Mz_block(:, :, n))
    if (verbose > 1) write(stdout, '(4x, "calculation of Mz is done")'); flush(stdout)

    ! Connect calculated Jx,Jy blocks
    Mz(:, 1:n_overlap) = Mz_block(:, 1:n_overlap, 1)
    Mz(:, nx - ny + 1:nx) = Mz_block(:, :, n)
    if (verbose > 1) write(stdout, '(4x, "connection of blocks (1st step) is done")'); flush(stdout)
    do i = 1, n - 1
      Mz(:, n_overlap + (i - 1) * step + 1 : i * step) = Mz_block(:, n_overlap + 1: ny - n_overlap, i)
    end do
    if (verbose > 1) write(stdout, '(4x, "connection of blocks (2nd step) is done")'); flush(stdout)
    do i = 1, n - 2
      forall(j = 1:n_overlap) Mz(:, i * step + j) = &
          & (Mz_block(:, ny - n_overlap + j, i) * fermi_dirac((j - 1) / dble(n_overlap - 1), beta = f_d_factor) + &
          &  Mz_block(:, j, i + 1) * fermi_dirac((n_overlap - j) / dble(n_overlap - 1), beta = f_d_factor))
    end do
    if (verbose > 1) write(stdout, '(4x, "connection of blocks (3rd step) is done")'); flush(stdout)
    do i = 1, n_overlap
      Mz(:, (n - 1) * step + i) = &
        & (Mz_block(:, ny - n_overlap + i, n - 1) * fermi_dirac((i - 1) / dble(n_overlap - 1), beta = f_d_factor) + &
        &  Mz_block(:, ((n - 1) * step + i) - (nx - ny), n) &
        &* fermi_dirac((n_overlap - i) / dble(n_overlap - 1), beta = f_d_factor))
    end do
    if (verbose > 1) write(stdout, '(4x, "connection of blocks (4th step) is done")'); flush(stdout)

    ! Calculate Jx, Jy from Mz
    inv_sc_thickness = 1 / sc_thickness * 1e3
    forall(i = 1:ny, j = 1:nx) Jx(i, j) = diff_2d(Mz, 1, x_interval, y_interval_interp, j, i) * inv_sc_thickness
    forall(i = 1:ny, j = 1:nx) Jy(i, j) = -diff_2d(Mz, 2, x_interval, y_interval_interp, j, i) * inv_sc_thickness

! Jx = Mz
  end subroutine obtain_Jx_Jy_by_inverse

  subroutine calc_J_tot()
    use io_data, only : J_tot, Jx, Jy
    implicit none

    J_tot(:, :) = sqrt(Jx(:, :)**2 + Jy(:, :)**2)
  end subroutine calc_J_tot

  subroutine calc_invG()
    use io_data, only : sc_thickness, dz, y_interval_interp, x_interval, invG, Bmap_interp
    use algebra, only : mat_inv
    implicit none

    integer(i4b)          :: n
    integer(i4b)          :: i, j, k, l
    real(dp), allocatable :: G(:, :)
    real(dp)              :: r, dxy, inv_sc_thickness

    n = ubound(Bmap_interp, 1)
    allocate(G(n * n, n * n), invG(n * n, n * n))

    ! B(y, x)
    dxy = x_interval * y_interval_interp
    do i = 1, n
      do j = 0, n - 1
        do k = 1, n
          do l = 0, n - 1
            r = sqrt(((i - k) * x_interval)**2 + ((j - l) * y_interval_interp)**2 + dz**2)
            G(i + j * n, k + l * n) = 1.0e-7 * (3.0d0 * dz**2 - r**2) / (r**5) * dxy
          end do
        end do
      end do
    end do

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


    if (d == 1) then ! y
      n = ubound(map, 1)
      if (py > 1 .and. py < n) then
        slope = (map(py + 1, px) - map(py - 1, px)) / (2.0d0 * dy)
      else if (py > 1) then
        slope = (map(py, px) - map(py - 1, px)) / dy
      else
        slope = (map(py + 1, px) - map(py, px)) / dy
      end if
    else if (d == 2) then ! x
      n = ubound(map, 2)
      if (px > 1 .and. px < n) then
        slope = (map(py, px + 1) - map(py, px - 1)) / (2.0d0 * dx)
      else if (px > 1) then
        slope = (map(py, px) - map(py, px - 1)) / dx
      else
        slope = (map(py, px + 1) - map(py, px)) / dx
      end if
    end if
  end function diff_2d

end module calc














