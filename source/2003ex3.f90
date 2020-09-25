program main
use utils, only: dp
use fem, only: fem_mass_1d, fem_stiff_1d
use inversion_1d, only: invert1_1d, invert2_1d, invert3_1d
implicit none

! ------------------------------------------------------------------------------
! Declarations
! ------------------------------------------------------------------------------

! Command-line inputs
character(100) :: filein

! Inputs from filein
integer :: method
real(dp) :: diffusionCoeff
real(dp) :: order
real(dp) :: focus
real(dp) :: t1
real(dp) :: t2
integer :: nz
real(dp) :: theta
real(dp) :: dx
real(dp) :: dt
integer :: nx
integer :: nt
real(dp), allocatable :: xpts(:)
real(dp), allocatable :: times(:)

! Output to fileout
real(dp), allocatable :: soln(:,:)

! For file I/O
character(100) :: ignore
integer :: u

! Local variables
real(dp) :: incline
real(dp), allocatable :: Imat_vals(:)
integer, allocatable :: Imat_cols(:)
integer, allocatable :: Imat_rowptrb(:), Imat_rowptre(:)
real(dp), allocatable :: Amat_vals(:)
integer, allocatable :: Amat_cols(:)
integer, allocatable :: Amat_rowptrb(:), Amat_rowptre(:)

! Error catching
integer :: stat
character(100) :: msg

! ------------------------------------------------------------------------------
! Read parameters
! ------------------------------------------------------------------------------

! Make sure an input and output file have been specified
if (command_argument_count() /= 1) then
    write(*,*) 'main: expected 1 input argument.'
    stop
end if

! Get input filename
call get_command_argument(1, filein)

! Read data from filein
open(newunit=u, file=filein, status='old', action='read', iostat=stat, iomsg=msg)
if (stat /= 0) then
    write(*,*) 'main: open: ', msg
    stop
end if
read(u,*) method, ignore
if (method < 1 .or. method > 3) then
    write(*,*) 'main: method must be 1, 2, or 3.'
    stop
end if
read(u,*) diffusionCoeff, ignore
read(u,*) order, ignore
read(u,*) focus, ignore
read(u,*) t1, ignore
read(u,*) t2, ignore
read(u,*) nz, ignore
read(u,*) theta, ignore
read(u,*) dx, ignore
read(u,*) dt, ignore
read(u,*) nx, ignore
read(u,*) nt, ignore
allocate(xpts(nx))
read(u,*) xpts
allocate(times(nt))
read(u,*) times
close(u)

! ------------------------------------------------------------------------------
! Computation
! ------------------------------------------------------------------------------

incline = 0.5_dp * atan(2.0_dp + focus)

! Get matrices
call fem_mass_1d(xpts, Imat_rowptrb, Imat_rowptre, Imat_cols, Imat_vals, stat)
if (stat /= 0) stop
call fem_stiff_1d(xpts, diffusivity_fn, dx, Amat_rowptrb, Amat_rowptre, &
    Amat_cols, Amat_vals, stat)
if (stat /= 0) stop
Amat_vals = -1.0_dp * Amat_vals

! Get solution
if (method == 1) then
    call invert1_1d(Amat_rowptrb, Amat_rowptre, Amat_cols, Amat_vals, &
        Imat_vals, &
        order, .false., fhat_fn, ic_fn, xpts, dx, focus, incline, t1, t2, nz, &
        theta, times, soln, stat)
else if (method == 2) then
    call invert2_1d(Amat_rowptrb, Amat_rowptre, Amat_cols, Amat_vals, &
        Imat_rowptrb, Imat_rowptre, Imat_cols, Imat_vals, &
        order, .false., f_fn, fhat_fn, ic_fn, xpts, dx, focus, incline, &
        t2, nz, theta, times, dt, soln, stat)
else if (method == 3) then
    call invert3_1d(Amat_rowptrb, Amat_rowptre, Amat_cols, Amat_vals, &
        Imat_rowptrb, Imat_rowptre, Imat_cols, Imat_vals, &
        order, .false., f_fn, ic_fn, xpts, dx, t2, nz, theta, &
        times, dt, soln, stat)
end if
if (stat /= 0) stop

! Clean up
deallocate(Amat_rowptrb, Amat_rowptre, Amat_cols, Amat_vals, Imat_rowptrb, &
    Imat_rowptre, Imat_cols, Imat_vals)

! ------------------------------------------------------------------------------
! Output
! ------------------------------------------------------------------------------

write(*,*) xpts(2:nx-1)
write(*,*) times
write(*,*) soln

deallocate(xpts, times, soln)

! --- End of program --- !

! ------------------------------------------------------------------------------
! Problem functions
! ------------------------------------------------------------------------------

contains

    ! Initial condition
    real(dp) function ic_fn(x) result(val)
        real(dp), intent(in) :: x
        val = sin(x) + sin(2.0_dp*x)
    end function ic_fn

    ! Uniform diffusivity
    real(dp) function diffusivity_fn(x) result(val)
        real(dp), intent(in) :: x
        val = diffusionCoeff
    end function diffusivity_fn

    ! A constant function
    real(dp) function one_fn(x) result(val)
        real(dp), intent(in) :: x
        val = 1.0_dp
    end function one_fn

    ! Source term
    real(dp) function f_fn(x,t) result(val)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: t
        val = sin(x) * exp(-t) + sin(2.0_dp*x) * exp(-2.0_dp*t) * (2.0_dp*cos(t) - sin(t))
    end function f_fn

    ! Laplace-transformed source term
    complex(dp) function fhat_fn(x,z) result(val)
        real(dp), intent(in) :: x
        complex(dp), intent(in) :: z
        val = sin(x) * (1.0_dp + z)**(-1.0_dp) + sin(2.0_dp*x) * (2.0_dp*z + 3.0_dp) / ((z + 2.0_dp)**(2.0_dp) + 1.0_dp)
    end function fhat_fn

end program main
