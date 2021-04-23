program main
use utils, only: dp, pi
use inversion_0d, only: invert1_0d, invert2_0d, invert3_0d
implicit none

! ------------------------------------------------------------------------------
! Declarations
! ------------------------------------------------------------------------------

    ! Inputs
    integer :: method
    real(dp) :: a
    real(dp) :: order
    real(dp) :: ic
    real(dp) :: focus
    real(dp) :: t1
    real(dp) :: t2
    integer :: nz
    real(dp) :: theta
    real(dp) :: dt
    integer :: nt
    real(dp), allocatable :: times(:)

    ! Output
    real(dp), allocatable :: soln(:)

    ! For file I/O
    character(100) :: filein
    character(100) :: ignore
    integer :: u

    ! Local variables
    real(dp) :: incline
    integer :: j

    ! Error catching
    integer :: stat
    character(100) :: msg

! ------------------------------------------------------------------------------
! Get input
! ------------------------------------------------------------------------------

    ! Make sure an input has been specified
    if (command_argument_count() /= 1) then
        write(*,*) 'main: expected 1 input argument.'
        stop
    end if

    ! Get input file name
    call get_command_argument(1, filein)

    ! Open input file
    open(newunit=u, file=filein, status='old', action='read', iostat=stat, &
        iomsg=msg)
    if (stat /= 0) then
        write(*,*) 'main: open: ', msg
        stop
    end if

    ! Read input
    read(u,*) method, ignore
    if (method < 1 .or. method > 3) then
        write(*,*) 'main: method must be 1, 2, or 3.'
        stop
    end if
    read(u,*) a, ignore
    read(u,*) order, ignore
    read(u,*) ic, ignore
    read(u,*) focus, ignore
    read(u,*) t1, ignore
    read(u,*) t2, ignore
    read(u,*) nz, ignore
    read(u,*) theta, ignore
    read(u,*) dt, ignore
    read(u,*) nt, ignore
    allocate(times(nt))
    read(u,*) times

    close(u)

! ------------------------------------------------------------------------------
! Get solution
! ------------------------------------------------------------------------------

    ! This is problem-specific
    incline = 0.5_dp * atan((1.0_dp + focus) / pi)

    if (method == 1) then
        call invert1_0d(a, order, fhat, ic, focus, incline, t1, t2, nz, theta, &
            times, soln, stat)
    else if (method == 2) then
        call invert2_0d(a, order, f, fhat, ic, focus, incline, t2, nz, theta, &
            dt, times, soln, stat)
    else if (method == 3) then
        call invert3_0d(a, order, f, ic, t2, nz, theta, dt, times, &
            soln, stat)
    end if
    if (stat /= 0) stop

! ------------------------------------------------------------------------------
! Output, cleaning up
! ------------------------------------------------------------------------------

    write(*,*) times
    write(*,*) soln
    deallocate(soln)
    deallocate(times)

! ------------------------------------------------------------------------------
! Problem functions
! ------------------------------------------------------------------------------

contains

    ! Source term
    real(dp) function f(t) result(val)
        real(dp), intent(in) :: t
        val = exp(-t) * cos(pi*t)
    end function f

    ! Laplace-transformed source term
    complex(dp) function fhat(z) result(val)
        complex(dp), intent(in) :: z
        val = (z + 1.0_dp) / ((z + 1.0_dp)**2.0_dp + pi**2.0_dp)
    end function fhat

end program main
