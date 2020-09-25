program main
use utils, only: dp, pi
use inversion_0d, only: invert1_0d, invert2_0d, invert3_0d
use mpi_f08
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

    ! For MPI
    integer, parameter :: master = 0
    integer :: rank
    real(dp) :: tstart, tstop

    ! Local variables
    real(dp) :: incline
    integer :: j

    ! Error catching
    integer :: stat
    character(100) :: msg

! ------------------------------------------------------------------------------
! Initialize MPI
! ------------------------------------------------------------------------------

    call MPI_Init(stat)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, stat)

! ------------------------------------------------------------------------------
! Get input
! ------------------------------------------------------------------------------

    ! Master core does I/O
    if (rank == master) then

        ! Make sure an input has been specified
        if (command_argument_count() /= 1) then
            write(*,*) 'main: expected 1 input argument.'
            call MPI_Abort(MPI_COMM_WORLD, 1, stat)
        end if

        ! Get input file name
        call get_command_argument(1, filein)

        ! Open input file
        open(newunit=u, file=filein, status='old', action='read', iostat=stat, &
            iomsg=msg)
        if (stat /= 0) then
            write(*,*) 'main: open: ', msg
            call MPI_Abort(MPI_COMM_WORLD, 1, stat)
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

    end if

    ! Distribute scalar-valued data to all cores
    call MPI_Bcast(method, 1, MPI_INTEGER, master, MPI_COMM_WORLD, stat)
    call MPI_Bcast(a, 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, stat)
    call MPI_Bcast(order, 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, stat)
    call MPI_Bcast(ic, 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, stat)
    call MPI_Bcast(focus, 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, stat)
    call MPI_Bcast(t1, 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, stat)
    call MPI_Bcast(t2, 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, stat)
    call MPI_Bcast(nz, 1, MPI_INTEGER, master, MPI_COMM_WORLD, stat)
    call MPI_Bcast(theta, 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, stat)
    call MPI_Bcast(dt, 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, stat)
    call MPI_Bcast(nt, 1, MPI_INTEGER, master, MPI_COMM_WORLD, stat)

    ! Non-master cores need to allocate array-valued data first
    if (rank /= master) then
        allocate(times(nt))
    end if

    ! Distribute array-valued data
    call MPI_Bcast(times, nt, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, stat)

! ------------------------------------------------------------------------------
! Get solution
! ------------------------------------------------------------------------------

    if (rank == master) tstart = MPI_Wtime()

    ! This is problem-specific
    incline = 0.5_dp * atan((1.0_dp + focus) / pi)

    if (method == 1) then
        call invert1_0d(a, order, fhat, ic, focus, incline, t1, t2, nz, theta, &
            times, master, MPI_COMM_WORLD, soln, stat)
    else if (method == 2) then
        call invert2_0d(a, order, f, fhat, ic, focus, incline, t2, nz, theta, &
            dt, times, master, MPI_COMM_WORLD, soln, stat)
    else if (method == 3) then
        call invert3_0d(a, order, f, ic, t2, nz, theta, dt, times, master, &
            MPI_COMM_WORLD, soln, stat)
    end if
    if (stat /= 0) call MPI_Abort(MPI_COMM_WORLD, 1, stat)

    if (rank == master) tstop = MPI_Wtime()

! ------------------------------------------------------------------------------
! Output, cleaning up
! ------------------------------------------------------------------------------

    if (rank == master) then
        write(*,*) tstop - tstart
        write(*,*) times
        write(*,*) soln
        deallocate(soln)
    end if

    deallocate(times)
    call MPI_Finalize(stat)

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
