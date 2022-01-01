!! Copyright (C) Hishiro T. Hirose
! This file is distributed under the terms of the GNU general public license.
! See http://www.gnu.org/copyleft/gpl.txt .
!
module algebra

  implicit none

!  public  :: interpolate_y
  public  :: interpolate1D
  public  :: mat_inv
  public  :: fermi_dirac

  contains

  subroutine interpolate1D(y, y_interp, B, B_interp)
    use constants, only : dp, i4b, pi

    implicit none

    real(dp), intent(in)            :: y(:)
    real(dp), intent(in)            :: y_interp(:)
    real(dp), intent(in)            :: B(:)
    real(dp), intent(out)           :: B_interp(:)
    integer(i4b)                    :: n, m, i, j
    real(dp)                        :: dy, t
    real(dp), allocatable           :: an(:), bn(:), cn(:), dn(:), B_prime(:), delta(:)
    real(dp), allocatable           :: matA(:,:)

    n = ubound(y, 1)
    m = ubound(y_interp, 1)
    allocate(an(n - 1), bn(n - 1), cn(n - 1), dn(n - 1), B_prime(n), delta(n))
    allocate(matA(n, n))
    i = 0
    matA = 0.0d0 ; matA(1,1) = 2.0d0 ; matA(n, n) = 2.0d0
    forall(i = 2:n - 1) matA(i, i) = 4.0d0
    forall(i = 2:n) matA(i - 1, i) = 1.0d0
    forall(i = 2:n) matA(i, i - 1) = 1.0d0
    forall(i = 2:n - 1) delta(i) = 3.0d0 * (B(i + 1) - B(i - 1))
    delta(1) = 3.0d0 * (B(2) - B(1))
    delta(n) = 3.0d0 * (B(n) - B(n - 1))
    
    B_prime = matmul(mat_inv(matA), delta)
    dn = B(1:n - 1)
    cn = B_prime(1:n - 1)
    forall(i = 1:n - 1) bn(i) = 3.0d0 * (B(i + 1) - B(i)) - (2.0d0 * B_prime(i) + B_prime(i + 1))
    forall(i = 1:n - 1) an(i) = 2.0d0 * (B(i) - B(i + 1)) + (B_prime(i) + B_prime(i + 1))

    dy = y(2) - y(1)
    do i = 1, m
      j = 1 + floor((y_interp(i) - y(1)) / dy)
      if (j < 1) then
        j = 1
      else if (j > n - 1) then
        j = n - 1
      end if
      t = (y_interp(i) - y(j))/dy
      B_interp(i) = an(j) * t**3 + bn(j) * t**2 + cn(j) * t + dn(j)
    end do
  end subroutine interpolate1D

  function mat_inv(A) result(C)
    ! Calculate the pseudo-inverse matrix C(n, m) of A(m, n) using LAPACK
    ! Ref: http://www.netlib.org/lapack/explore-html/d7/d3b/group__double_g_esolve_gaa6ed601d0622edcecb90de08d7a218ec.html#gaa6ed601d0622edcecb90de08d7a218ec
    use constants, only : i4b, dp

    implicit none

    real(dp),    intent(in)   :: A(:, :)        ! input matrix
    real(dp)                  :: C(ubound(A, 2), ubound(A, 1))
                                                ! The pseudo-inverse matrix
    real(dp)                  :: A_(ubound(A, 1), ubound(A, 2))
    integer(i4b)              :: m              ! Number of rows of the matrix A
    integer(i4b)              :: n              ! Number of columns of the matrix A
    integer(i4b)              :: NRHS           ! Number of right hand sides
    integer(i4b)              :: LDA            ! Leading dimension of A
    integer(i4b)              :: LDB            ! Leading dimension of B
    real(dp)                  :: B(max(ubound(A, 1), ubound(A, 2)), max(ubound(A, 1), ubound(A, 2)))
                                                ! Right hand side matrix
    real(dp)                  :: S(min(ubound(A, 1), ubound(A, 2)))
                                                ! Singular values of A in decreasing order
    real(dp)                  :: RCond          ! The effective rank of A is determined in machine precision
    integer(i4b)              :: Rank           ! The effective rank of A
    integer(i4b)              :: LWork          ! The dimension of Work
    real(dp),    allocatable  :: Work(:)        ! Array for working area
    real(dp)                  :: tmp(1)         ! Temporary array to receive the optimal LWork
    integer(i4b)              :: Info           ! = 0: successful exit
                                                ! < 0: '-Info'-th argument has an illegal value
                                                ! > 0: SVD computing failed to converge
    integer(i4b)              :: i, j
   
    A_ = A
    m = ubound(A, 1)
    n = ubound(A, 2)
    NRHS = max(m, n)
    LDA = m
    LDB = NRHS
    RCond = -1.0_8
    
    ! Receive the optimal LWork
    call DGELSS(m, n, NRHS, A_, LDA, B, LDB, S, RCond, Rank, tmp, -1, Info)
    LWork = int(tmp(1))
    allocate(Work(LWork))
    
    ! Initialize B as identity matrix
    B = 0.0_8
    do i = 1, LDB
      B(i, i) = 1.0_8
    end do
    
    ! Compute pseudo-inverse
    call DGELSS(m, n, NRHS, A_, LDA, B, LDB, S, RCond, Rank, Work, LWork, Info)
    
    do i = 1, n
      do j = 1, m
        C(i, j) = B(i, j)
      end do
    end do
  end function mat_inv

pure function fermi_dirac(val, beta) result(fd)
    use constants, only : dp, i4b, pi

    implicit none

    real(dp), intent(in)           :: val
    real(dp), intent(in), optional :: beta
    real(dp)                       :: fd

    if (.not. present(beta)) then
      fd = 1.0d0 / (dexp(1.0d0 * (val - 0.5d0)) + 1.0d0)
    else
      fd = 1.0d0 / (dexp(beta * (val - 0.5d0)) + 1.0d0)
    end if
  end function fermi_dirac

end module algebra

  
