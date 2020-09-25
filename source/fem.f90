module fem
use utils, only: dp
use sparse_utils, only: csr_kron_sym
use integration, only: defint_fn_trap_real, defint_trap_cmplx
implicit none
private
public fem_mass_1d, fem_stiff_1d, fem_vec_cmplx_1d, fem_vec_evalpts_1d, &
    fem_mats_2d, fem_vec_cmplx_2d

! Module containing subroutines related to the finite element method.
! Submodules: fem_mat_1d.f90, fem_vec_1d.f90, fem_mat_2d.f90

interface

! Function interface for the problem diffusivity.
real(dp) function diffusivity(x)
    use utils, only: dp
    implicit none
    real(dp), intent(in) :: x
end function



! ------------------------------------------------------------------------------
! One-dimensional subroutines
! ------------------------------------------------------------------------------

! Computes the Galerkin linear spline FEM mass matrix in MKL CSR format.
! This matrix is a representation of the identity operator. NOTE: since the
! matrix is symmetric, this subroutine computes only the upper triangle.
! Source code: fem_mat_1d.f90
module subroutine fem_mass_1d(pts, rowptrb, rowptre, cols, vals, stat)

    ! --- Inputs --- !

    ! Grid points
    real(dp), intent(in) :: pts(:)


    ! --- Outputs --- !

    ! Mass matrix in MKL CSR format. There is no need to allocate these arrays
    ! before calling the subroutine.
    integer, allocatable, intent(out) :: rowptrb(:)
    integer, allocatable, intent(out) :: rowptre(:)
    integer, allocatable, intent(out) :: cols(:)
    real(dp), allocatable, intent(out) :: vals(:)
    
    ! Error code
    integer, intent(out) :: stat

end subroutine fem_mass_1d



! Computes the Galerkin linear spline FEM stiffness matrix in MKL CSR format.
! This matrix is a representation of  d2/dx2 (sigma(x) u(x)),  where u is the
! unknown function and sigma is the diffusivity function.
! Source code: fem_mat_1d.f90
module subroutine fem_stiff_1d(pts, diffusivity_fn, dx, rowptrb, rowptre, &
    cols, vals, stat)

    ! --- Inputs --- !

    ! Grid points
    real(dp), intent(in) :: pts(:)

    ! Diffusivity function
    procedure(diffusivity) :: diffusivity_fn

    ! Maximum allowable integration interval width
    real(dp), intent(in) :: dx


    ! --- Outputs --- !

    ! Stiffness matrix in MKL CSR format. There is no need to allocate these
    ! arrays before calling the subroutine.
    integer, allocatable, intent(out) :: rowptrb(:)
    integer, allocatable, intent(out) :: rowptre(:)
    integer, allocatable, intent(out) :: cols(:)
    real(dp), allocatable, intent(out) :: vals(:)
    
    ! Error code
    integer, intent(out) :: stat

end subroutine fem_stiff_1d



! Computes the evaluation points needed to compute the load vector using the 
! desired dx, so that the source term evaluations may be performed by the 
! calling subroutine.
! Source code: fem_vec_1d.f90
module subroutine fem_vec_evalpts_1d(gridpts, dx, n_evals_total, n_evals, evalpts, stat)

    ! --- Inputs --- !

    ! Grid points
    real(dp), intent(in) :: gridpts(:)

    ! Maximum allowable integration interval width
    real(dp), intent(in) :: dx


    ! --- Outputs --- !

    ! The total number of evaluation points to be used
    integer, intent(out) :: n_evals_total

    ! Array containing the number of interior evaluation points per interval.
    ! It is not necessary to allocate this array beforehand. Upon exit, its
    ! size is  size(gridpts) - 1.
    integer, allocatable, intent(out) :: n_evals(:)

    ! The evaluation points. It is not necessary to allocate this array
    ! beforehand. Upon exit, its size is n_evals_total.
    real(dp), allocatable, intent(out) :: evalpts(:)

    ! Error code
    integer :: stat

end subroutine fem_vec_evalpts_1d



! Computes the FEM load vector from the given source term evaluations and
! grid information computed by fem_vec_evalpts_1d.
! Source code: fem_vec_1d.f90
module subroutine fem_vec_cmplx_1d(n_evals, evalpts, evals, femvec, stat)

    ! --- Inputs --- !

    ! Array containing the number of interior evaluation points per interval,
    ! as computed by fem_vec_evalpts_1d. It should be of size n+1, where n is
    ! the number of interior grid points (not including the endpoints).
    integer, intent(in) :: n_evals(:)

    ! The evaluation points, as computed by fem_vec_evalpts_1d. It should be of
    ! size sum(n_evals) + n + 2, where n is the number of interior grid points.
    real(dp), intent(in) :: evalpts(:)

    ! The function evaluations at evalpts. It should be the same size as
    ! evalpts.
    complex(dp), intent(in) :: evals(:)


    ! --- Outputs --- !

    ! The FEM load vector. NOTE: this array must be of size n, allocated
    ! before calling the subroutine. 
    complex(dp), intent(out) :: femvec(:)

    ! Error code
    integer, intent(out) :: stat

end subroutine fem_vec_cmplx_1d



! ------------------------------------------------------------------------------
! Two-dimensional subroutines
! ------------------------------------------------------------------------------

! Computes both FEM matrices in 2d, using bilinear basis functions. NOTE: since
! both matrices are symmetric, this subroutine computes only the upper triangle.
! Source code: fem_mat_2d.f90
module subroutine fem_mats_2d(xpts, ypts, diffusivity_x_fn, diffusivity_y_fn, &
    dx, dy, mass_rowptrb, mass_rowptre, mass_cols, mass_vals, stiff_rowptrb, &
    stiff_rowptre, stiff_cols, stiff_vals, stat)

    ! --- Inputs --- !

    ! Grid points
    real(dp), intent(in) :: xpts(:), ypts(:)

    ! Diffusivity functions
    procedure(diffusivity) :: diffusivity_x_fn, diffusivity_y_fn

    ! Maximum allowable integration interval widths
    real(dp), intent(in) :: dx, dy


    ! --- Outputs --- !

    ! Matrices in MKL CSR format. There is no need to allocate these
    ! arrays before calling the subroutine.
    integer, allocatable, intent(out) :: mass_rowptrb(:)
    integer, allocatable, intent(out) :: mass_rowptre(:)
    integer, allocatable, intent(out) :: mass_cols(:)
    real(dp), allocatable, intent(out) :: mass_vals(:)

    integer, allocatable, intent(out) :: stiff_rowptrb(:)
    integer, allocatable, intent(out) :: stiff_rowptre(:)
    integer, allocatable, intent(out) :: stiff_cols(:)
    real(dp), allocatable, intent(out) :: stiff_vals(:)
    
    ! Error code
    integer, intent(out) :: stat

end subroutine fem_mats_2d



! Computes the FEM load vector from the given source term evaluations and
! grid information, which should be computed as two separate calls to
! fem_vec_evalpts_1d (one for x, one for y).
! Source code: fem_vec_2d.f90
module subroutine fem_vec_cmplx_2d(n_evals_x, n_evals_y, evalpts_x, evalpts_y, &
    evals, femvec, stat)

    ! --- Inputs --- !

    ! Grid information - see the descriptions in fem_vec_cmplx_1d above.
    integer, intent(in) :: n_evals_x(:), n_evals_y(:)
    real(dp), intent(in) :: evalpts_x(:), evalpts_y(:)

    ! The function evaluations - x index first, y second.
    complex(dp), intent(in) :: evals(:,:)


    ! --- Outputs --- !

    ! The load vector. NOTE: this array must be of size (nx*ny), allocated
    ! before calling the subroutine.
    complex(dp), intent(out) :: femvec(:)

    ! Error code
    integer, intent(out) :: stat

end subroutine fem_vec_cmplx_2d



end interface

end module fem