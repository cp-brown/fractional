submodule (integration) integration_trap_cmplx
implicit none

! Submodule containing subroutines for the integration of complex-valued
! integrands.

contains

module procedure defint_fn_trap_cmplx

    ! --- Declarations of local variables --- !

    integer :: npts     ! Number of interior points
    real(dp) :: h       ! Spacing
    integer :: i        ! Iterator

    ! --- Begin program --- !

    npts = ceiling((b-a)/dt) - 1

    ! Return 0 if npts is negative (a = b)
    if (npts < 0) then
        integral = (0.0_dp, 0.0_dp)
        return
    end if

    h = (b - a) / (npts + 1.0_dp)
    integral = 0.5_dp * (fn(a) + fn(b))
    do i = 1, npts
        integral = integral + fn( a + i*h )
    end do
    integral = integral * h
    return

end procedure defint_fn_trap_cmplx



module procedure defint_trap_cmplx

    ! --- Declarations of local variables --- !

    ! Number of interior points
    integer :: npts

    ! Iterator
    integer :: i
    

    ! --- Begin program --- !

    npts = size(evals) - 2
    integral = 0.5_dp * (evals(0) + evals(npts+1))
    do i = 1, npts
        integral = integral + evals(i)
    end do
    integral = integral * (b - a) / (npts + 1.0_dp)
    return

end procedure defint_trap_cmplx

end submodule integration_trap_cmplx
