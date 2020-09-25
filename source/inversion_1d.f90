module inversion_1d
use utils, only: dp, pi, imagunit
use integration, only: defint_fn_trap_real, defint_fn_trap_cmplx, indefint_fn_trap_real
use fem, only: fem_vec_cmplx_1d, fem_vec_evalpts_1d
use mkl_spblas
use mkl_dss
use blas95, only: gemv
use iso_c_binding
implicit none
private
public invert1_1d, invert2_1d, invert3_1d

interface

! Function interface for the initial condition
real(dp) function ic(x)
    use utils, only: dp
    implicit none
    real(dp), intent(in) :: x
end function ic

! Function interface for the source term
real(dp) function f(x,t)
    use utils, only: dp
    implicit none
    real(dp), intent(in) :: x
    real(dp), intent(in) :: t
end function f

! Function interface for the Laplace-transformed source
complex(dp) function fhat(x,z)
    use utils, only: dp
    implicit none
    real(dp), intent(in) :: x
    complex(dp), intent(in) :: z
end function fhat



! Solves a parabolic problem using a Laplace transform method. First of three
! methods described in McLean and Thomee (2009).
module subroutine invert1_1d(Amat_rowptrb, Amat_rowptre, Amat_cols, Amat_vals,  &
    Imat_vals, order, homogeneous, &
    fhat_fn, ic_fn, x_pts, dx, focus, incline, t1, t2, nz, theta, times, soln, stat)

    ! --- Inputs --- !

    ! Laplacian matrix in MKL CSR format
    integer, intent(in) :: Amat_rowptrb(:)
    integer, intent(in) :: Amat_rowptre(:)
    integer, intent(in) :: Amat_cols(:)
    real(dp), intent(in) :: Amat_vals(:)

    ! Identity operator matrix in MKL CSR format
    real(dp), intent(in) :: Imat_vals(:)

    ! Fractional differential order
    real(dp), intent(in) :: order

    ! Set to true if the problem is homogeneous, false if not
    logical :: homogeneous

    ! Laplace transform of the source term. Should be a complex-valued
    ! function of one complex-valued variable. See the interface above.
    ! Will be ignored if homogeneous = .true.
    procedure(fhat) :: fhat_fn

    ! Initial condition function
    procedure(ic) :: ic_fn

    ! Spatial locations to be considered
    real(dp), intent(in) :: x_pts(:)

    ! Integration interval parameter for computing the FEM vector
    real(dp), intent(in) :: dx

    ! Focus of the integration contour
    real(dp), intent(in) :: focus

    ! Limiting incline of the integration contour
    real(dp), intent(in) :: incline

    ! Time interval in which accurate results are desired
    real(dp), intent(in) :: t1, t2

    ! Number of contour points to use minus one
    integer, intent(in) :: nz

    ! Parameter related to contour spacing. Should be between 0 and 1.
    real(dp), intent(in) :: theta

    ! Times to solve at
    real(dp), intent(in) :: times(:)

    
    ! --- Outputs --- !

    ! Upon successful exit, an nx-times-nt array containing the solution,
    ! where nx is the number of interior points (NOT including endpoints).
    ! It is not necessary to allocate before calling the subroutine.
    real(dp), allocatable, intent(out) :: soln(:,:)

    ! Error code. Will be zero upon a successful exit, nonzero otherwise.
    integer, intent(out) :: stat

end subroutine invert1_1d



! Solves a parabolic problem using a Laplace transform method. Second of three
! methods described in McLean and Thomee (2009).
module subroutine invert2_1d(Amat_rowptrb, Amat_rowptre, Amat_cols, Amat_vals,  &
    Imat_rowptrb, Imat_rowptre, Imat_cols, Imat_vals, order, homogeneous, f_fn, &
    fhat_fn, ic_fn, x_pts, dx, focus, incline, t2, nz, sigma, times, dt, soln, stat)

    ! --- Inputs --- !

    ! Laplacian matrix in MKL CSR format
    integer, intent(in) :: Amat_rowptrb(:)
    integer, intent(in) :: Amat_rowptre(:)
    integer, intent(in) :: Amat_cols(:)
    real(dp), intent(in) :: Amat_vals(:)

    ! Identity operator matrix in MKL CSR format
    integer, intent(in) :: Imat_rowptrb(:)
    integer, intent(in) :: Imat_rowptre(:)
    integer, intent(in) :: Imat_cols(:)
    real(dp), intent(in) :: Imat_vals(:)

    ! Fractional differential order
    real(dp), intent(in) :: order

    ! Set to true if the problem is homogeneous, false if not
    logical :: homogeneous

    ! The source term. Should be a real-valued function of one real-valued
    ! variable. See the interface above. Will be ignored if
    ! homogeneous = .true.
    procedure(f) :: f_fn

    ! Laplace transform of the source term. Should be a complex-valued
    ! function of one complex-valued variable. See the interface above.
    ! Will be ignored if homogeneous = .true.
    procedure(fhat) :: fhat_fn

    ! Initial condition function
    procedure(ic) :: ic_fn

    ! Spatial locations to be considered
    real(dp), intent(in) :: x_pts(:)

    ! Integration interval parameter for computing the FEM vector
    real(dp), intent(in) :: dx

    ! Focus of the integration contour
    real(dp), intent(in) :: focus

    ! Limiting incline of the integration contour
    real(dp), intent(in) :: incline

    ! Latest time at which accurate results are desired
    real(dp), intent(in) :: t2

    ! Number of contour points to use minus one
    integer, intent(in) :: nz

    ! Parameter related to problem regularity
    real(dp), intent(in) :: sigma

    ! Times to solve at
    real(dp), intent(in) :: times(:)

    ! Integration interval parameter for time integration
    real(dp), intent(in) :: dt

    
    ! --- Outputs --- !

    ! Upon successful exit, an nx-times-nt array containing the solution,
    ! where nx is the number of interior points (NOT including endpoints).
    ! It is not necessary to allocate before calling the subroutine.
    real(dp), allocatable, intent(out) :: soln(:,:)

    ! Error code. Will be zero upon a successful exit, nonzero otherwise.
    integer, intent(out) :: stat

end subroutine invert2_1d



! Solves a parabolic problem using a Laplace transform method. Second of three
! methods described in McLean and Thomee (2009).
module subroutine invert3_1d(Amat_rowptrb, Amat_rowptre, Amat_cols, Amat_vals,  &
    Imat_rowptrb, Imat_rowptre, Imat_cols, Imat_vals, order, homogeneous, f_fn, &
    ic_fn, x_pts, dx, t2, nz, sigma, times, dt, soln, stat)

    ! --- Inputs --- !

    ! Laplacian matrix in MKL CSR format
    integer, intent(in) :: Amat_rowptrb(:)
    integer, intent(in) :: Amat_rowptre(:)
    integer, intent(in) :: Amat_cols(:)
    real(dp), intent(in) :: Amat_vals(:)

    ! Identity operator matrix in MKL CSR format
    integer, intent(in) :: Imat_rowptrb(:)
    integer, intent(in) :: Imat_rowptre(:)
    integer, intent(in) :: Imat_cols(:)
    real(dp), intent(in) :: Imat_vals(:)

    ! Fractional differential order
    real(dp), intent(in) :: order

    ! Set to true if the problem is homogeneous, false if not
    logical :: homogeneous

    ! The source term. Should be a real-valued function of one real-valued
    ! variable. See the interface above. Will be ignored if
    ! homogeneous = .true.
    procedure(f) :: f_fn

    ! Initial condition function
    procedure(ic) :: ic_fn

    ! Spatial locations to be considered
    real(dp), intent(in) :: x_pts(:)

    ! Integration interval parameter for computing the FEM vector
    real(dp), intent(in) :: dx

    ! Latest time at which accurate results are desired
    real(dp), intent(in) :: t2

    ! Number of contour points to use minus one
    integer, intent(in) :: nz

    ! Parameter related to problem regularity
    real(dp), intent(in) :: sigma

    ! Times to solve at
    real(dp), intent(in) :: times(:)

    ! Integration interval parameter for time integration
    real(dp), intent(in) :: dt

    
    ! --- Outputs --- !

    ! Upon successful exit, an nx-times-nt array containing the solution,
    ! where nx is the number of interior points (NOT including endpoints).
    ! It is not necessary to allocate before calling the subroutine.
    real(dp), allocatable, intent(out) :: soln(:,:)

    ! Error code. Will be zero upon a successful exit, nonzero otherwise.
    integer, intent(out) :: stat

end subroutine invert3_1d

end interface

end module inversion_1d
