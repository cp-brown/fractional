submodule (integration) integration_trap_real
implicit none

! Submodule containing subroutines to integrate real-valued functions
! using the trapezoid rule.

contains

module procedure defint_fn_trap_real

    ! --- Declarations of local variables --- !

    ! Number of interior points
    integer :: npts

    ! Interval widths
    real(dp) :: h

    ! Iterator
    integer :: i


    ! --- Begin program --- !

    npts = ceiling((b-a)/dt) - 1

    ! Return 0 if npts is negative (a = b)
    if (npts < 0) then
        integral = 0.0_dp
        return
    end if

    h = (b - a) / (npts + 1.0_dp)
    integral = 0.5_dp * (fn(a) + fn(b))
    do i = 1, npts
        integral = integral + fn( a + i*h )
    end do
    integral = integral * h
    return

end procedure defint_fn_trap_real



module procedure indefint_fn_trap_real

    ! --- Declarations --- !

    ! Number of intervals to integrate over
    integer :: n

    ! Iterator
    integer :: j

    ! For error messages
    character(100) :: msg


    ! --- Begin program --- !

    n = size(delims) - 1

    ! Check vals is the right size
    if (size(vals) /= n) then
        write(*,*) 'indefint_fn_trap_real: unexpected size for input vals.'
        stat = -1
        return
    end if

    ! Compute integrals over intervals
    do j = 1, n
        vals(j) = defint_fn_trap_real(fn, delims(j), delims(j+1), dt)
    end do

    ! Compute partial sums
    vals(1) = vals(1) + ic
    do j = 2, n
        vals(j) = vals(j) + vals(j-1)
    end do

    ! Successful exit
    stat = 0
    return

end procedure indefint_fn_trap_real

end submodule integration_trap_real
