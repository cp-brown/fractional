module inversion_0d
use utils, only: dp, pi, imagunit
use integration, only: defint_fn_trap_real, defint_fn_trap_cmplx, indefint_fn_trap_real
implicit none
private
public invert1_0d, invert2_0d, invert3_0d

interface

! Function prototype for source term
real(dp) function f(t)
    use utils, only: dp
    implicit none
    real(dp), intent(in) :: t
end function f

! Function prototype for Laplace-transformed source
complex(dp) function fhat(z)
    use utils, only: dp
    implicit none
    complex(dp), intent(in) :: z
end function fhat



module subroutine invert1_0d(a, order, fhat_fn, ic, focus, incline, t1, t2, &
    nz, theta, times, soln, stat)

    ! --- Inputs --- !

    ! Scalar in the ODE
    real(dp), intent(in) :: a

    ! Fractional differential order
    real(dp), intent(in) :: order

    ! Laplace transform of the source term, see prototype above
    procedure(fhat) :: fhat_fn

    ! Initial value
    real(dp), intent(in) :: ic

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

    ! Times at which a solution is requested
    real(dp), intent(in) :: times(:)


    ! --- Outputs --- !

    ! Solution at requested times. No need to allocate beforehand.
    real(dp), allocatable, intent(out) :: soln(:)

    ! Error code
    integer, intent(out) :: stat

end subroutine invert1_0d



module subroutine invert2_0d(a, order, f_fn, fhat_fn, ic, focus, incline, t2, &
    nz, sigma, dt, times, soln, stat)

    ! --- Inputs --- !

    ! Scalar in the ODE
    real(dp), intent(in) :: a

    ! Fractional differential order
    real(dp), intent(in) :: order

    ! The source term. See prototype above
    procedure(f) :: f_fn

    ! Laplace transform of the source term, see prototype above
    procedure(fhat) :: fhat_fn

    ! Initial value
    real(dp), intent(in) :: ic

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

    ! Timestep for numerical integration
    real(dp), intent(in) :: dt

    ! Times at which a solution is requested
    real(dp), intent(in) :: times(:)


    ! --- Outputs --- !

    ! Solution at requested times. No need to allocate beforehand.
    real(dp), allocatable, intent(out) :: soln(:)

    ! Error code
    integer, intent(out) :: stat

end subroutine invert2_0d



module subroutine invert3_0d(a, order, f_fn, ic, t2, nz, sigma, dt, times, soln, stat)

    ! --- Inputs --- !

    ! Scalar in the ODE
    real(dp), intent(in) :: a

    ! Fractional differential order
    real(dp), intent(in) :: order

    ! The source term. See prototype above
    procedure(f) :: f_fn

    ! Initial value
    real(dp), intent(in) :: ic

    ! Latest time at which accurate results are desired
    real(dp), intent(in) :: t2

    ! Number of contour points to use minus one
    integer, intent(in) :: nz

    ! Parameter related to problem regularity
    real(dp), intent(in) :: sigma

    ! Timestep for numerical integration
    real(dp), intent(in) :: dt

    ! Times at which a solution is requested
    real(dp), intent(in) :: times(:)


    ! --- Outputs --- !

    ! Solution at requested times. No need to allocate beforehand.
    real(dp), allocatable, intent(out) :: soln(:)

    ! Error code
    integer, intent(out) :: stat

end subroutine invert3_0d

end interface

end module inversion_0d
