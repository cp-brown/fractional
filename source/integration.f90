module integration
use utils, only: dp
implicit none
private
public defint_fn_trap_real, defint_fn_trap_cmplx, indefint_fn_trap_real, defint_trap_cmplx

! Module containing subroutines for numerical integration.
! Submodules: integration_trap_real.f90, integration_trap_cmplx.f90

interface

! Interface for a real-valued integrand of a real variable
real(dp) function real_integrand(x)
    use utils, only: dp
    implicit none
    real(dp), intent(in) :: x
end function

! Interface for a complex-valued integrand of a real variable
complex(dp) function cmplx_integrand(x)
    use utils, only: dp
    implicit none
    real(dp), intent(in) :: x
end function



! Computes the integral of a real integrand over the interval [a,b], using the
! trapezoid rule, where the integrand is a given function.
! Source code: integration_trap_real.f90
real(dp) module function defint_fn_trap_real(fn, a, b, dt) result(integral)

    ! --- Inputs --- !

    ! The integrand. See the interface above.
    procedure(real_integrand) :: fn
    
    ! The interval of integration
    real(dp), intent(in) :: a, b

    ! The maximum allowable interval width for numerical integration.
    real(dp), intent(in) :: dt

end function defint_fn_trap_real



! Computes the numerical "indefinite integral" of a real integrand using the 
! trapezoid rule at the locations delims(2:), where the value at delims(1)
! is given as ic, and the integrand is a given function.
! Source code: integration_trap_real.f90
module subroutine indefint_fn_trap_real(fn, ic, delims, dt, vals, stat)

    ! --- Inputs --- !

    ! Integrand function. See the interface above.
    procedure(real_integrand) :: fn

    ! Value of the integral at delims(1)
    real(dp), intent(in) :: ic

    ! Integration interval delimiters.
    real(dp), intent(in) :: delims(:)

    ! Maximum allowable interval width for numerical integration
    real(dp), intent(in) :: dt


    ! --- Outputs --- !

    ! Values of the integral at delims(2:). Should be allocated beforehand.
    real(dp), intent(out) :: vals(:)

    ! Error code
    integer, intent(out) :: stat

end subroutine indefint_fn_trap_real



! Computes the integral of a complex integrand over the interval [a,b], using
! the trapezoid rule, where the integrand is a given function.
! Source code: integration_trap_cmplx.f90
complex(dp) module function defint_fn_trap_cmplx(fn, a, b, dt) result(integral)

    ! --- Inputs --- !

    ! The integrand. See the interface above.
    procedure(cmplx_integrand) :: fn

    ! The interval
    real(dp), intent(in) :: a, b
    
    ! Maximum allowable interval width for numerical integration.
    real(dp), intent(in) :: dt

end function defint_fn_trap_cmplx



! Computes the integral of a complex integrand over the interval [a,b], using
! the trapezoid rule, where evaluations of the integrand are given. It is
! assumed that the evaluations are given at uniformly-spaced points from a to
! b, inclusive.
! Source code: integration_trap_cmplx.f90
complex(dp) module function defint_trap_cmplx(evals, a, b) result(integral)

    ! --- Inputs --- !

    ! Integrand evaluations, assumed to be at uniformly-spaced points.
    complex(dp), intent(in) :: evals(0:)

    ! The interval of integration.
    real(dp), intent(in) :: a, b

end function defint_trap_cmplx

end interface

end module integration