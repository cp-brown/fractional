program main
use utils, only: dp, pi
use fem, only: fem_mats_2d
use inversion_2d, only: invert1_2d, invert2_2d, invert3_2d
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
real(dp) :: dy
real(dp) :: dt
integer :: nx
integer :: ny
integer :: nt
real(dp), allocatable :: xpts(:)
real(dp), allocatable :: ypts(:)
real(dp), allocatable :: times(:)

! Output to fileout
real(dp), allocatable :: soln(:,:,:)

! For file I/O
character(100) :: ignore
integer :: u

! Local variables
real(dp), allocatable :: Imat_vals(:)
integer, allocatable :: Imat_cols(:)
integer, allocatable :: Imat_rowptrb(:), Imat_rowptre(:)
real(dp), allocatable :: Amat_vals(:)
integer, allocatable :: Amat_cols(:)
integer, allocatable :: Amat_rowptrb(:), Amat_rowptre(:)
real(dp) :: incline

! Error catching
integer :: stat
character(100) :: msg

! ------------------------------------------------------------------------------
! Read parameters
! ------------------------------------------------------------------------------

! Make sure an input and output file have been specified
if (command_argument_count() /= 1) then
    write(*,*) 'main: expected 3 input arguments.'
    stop
end if

! Get input filenames
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
read(u,*) dy, ignore
read(u,*) dt, ignore
read(u,*) nx, ignore
read(u,*) ny, ignore
read(u,*) nt, ignore
allocate(xpts(nx))
read(u,*) xpts
allocate(ypts(ny))
read(u,*) ypts
allocate(times(nt))
read(u,*) times
close(u)

! ------------------------------------------------------------------------------
! Computation
! ------------------------------------------------------------------------------

incline = 0.5_dp * atan((1.0_dp + focus) / pi)

! Get matrices
call fem_mats_2d(xpts, ypts, diffusivity_x_fn, diffusivity_y_fn, dx, dy, &
    Imat_rowptrb, Imat_rowptre, Imat_cols, Imat_vals, Amat_rowptrb, Amat_rowptre, &
    Amat_cols, Amat_vals, stat)
if (stat /= 0) stop
Amat_vals = -1.0_dp * Amat_vals

! Get solution
if (method == 1) then
    call invert1_2d(Amat_rowptrb, Amat_rowptre, Amat_cols, Amat_vals, &
        Imat_rowptrb, Imat_rowptre, Imat_cols, Imat_vals, &
        order, .false., fhat_fn, ic_fn, xpts, ypts, dx, dy, focus, incline, &
        t1, t2, nz, theta, times, soln, stat)
else if (method == 2) then
    call invert2_2d(Amat_rowptrb, Amat_rowptre, Amat_cols, Amat_vals, &
        Imat_rowptrb, Imat_rowptre, Imat_cols, Imat_vals, &
        order, .false., f_fn, fhat_fn, ic_fn, xpts, ypts, dx, dy, focus, incline, &
        t2, nz, theta, times, dt, soln, stat)
else if (method == 3) then
    call invert3_2d(Amat_rowptrb, Amat_rowptre, Amat_cols, Amat_vals, &
        Imat_rowptrb, Imat_rowptre, Imat_cols, Imat_vals, &
        order, .false., f_fn, ic_fn, xpts, ypts, dx, dy, t2, nz, theta, &
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
write(*,*) ypts(2:ny-1)
write(*,*) times
write(*,*) soln

deallocate(xpts, ypts, times, soln)

! --- End of program --- !

! ------------------------------------------------------------------------------
! Problem functions
! ------------------------------------------------------------------------------

contains

    ! Useful for simplifying the other functions
    real(dp) function phi(x, y, j, k) result(val)
        real(dp), intent(in) :: x, y
        integer, intent(in) :: j, k
        val = sin(j*pi*x/4.0_dp) * sin(k*pi*y/4.0_dp)
    end function phi

    ! Initial condition
    real(dp) function ic_fn(x, y) result(val)
        real(dp), intent(in) :: x, y
        val = phi(x,y,1,1) - phi(x,y,2,1)
    end function ic_fn

    ! Uniform diffusivity in x
    real(dp) function diffusivity_x_fn(x) result(val)
        real(dp), intent(in) :: x
        val = diffusionCoeff
    end function diffusivity_x_fn

    ! Uniform diffusivity in y
    real(dp) function diffusivity_y_fn(x) result(val)
        real(dp), intent(in) :: x
        val = 1.0_dp
    end function diffusivity_y_fn

    ! ! Source term
    ! real(dp) function f_fn(x, y, t) result(val)
    !     real(dp), intent(in) :: x, y
    !     real(dp), intent(in) :: t
    !     val = exp(-t/2.0_dp) * cos(t) * phi(x,y,1,1) &
    !         + 0.5_dp * exp(-t) * cos(pi*t) * phi(x,y,2,1)
    ! end function f_fn

    ! ! Laplace-transformed source term
    ! complex(dp) function fhat_fn(x, y, z) result(val)
    !     real(dp), intent(in) :: x, y
    !     complex(dp), intent(in) :: z
    !     val = phi(x,y,1,1) * (z+0.5_dp) / ((z+0.5_dp)**2 + 1.0_dp) &
    !         + 0.5_dp * phi(x,y,2,1) * (z+1.0_dp) / ((z+1.0_dp)**2 + pi**2)
    ! end function fhat_fn

    ! Source term
    real(dp) function f_fn(x, y, t) result(val)
        real(dp), intent(in) :: x, y
        real(dp), intent(in) :: t
        val = 0.0_dp
    end function f_fn

    ! Laplace-transformed source term
    complex(dp) function fhat_fn(x, y, z) result(val)
        real(dp), intent(in) :: x, y
        complex(dp), intent(in) :: z
        val = 0.0_dp
    end function fhat_fn

end program main
