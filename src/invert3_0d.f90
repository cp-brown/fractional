submodule (inversion_0d) inversion_invert3_0d
implicit none

contains

module procedure invert3_0d
! ------------------------------------------------------------------------------
! Declarations for local variables
! ------------------------------------------------------------------------------

    ! For computing the contour
    real(dp) :: incline, r, gam, kappa, lambda, z_spacing
    complex(dp), allocatable :: z_pts(:), z_derivs(:)

    integer :: nt
    complex(dp) :: lhs
    real(dp) :: t, tm1
    complex(dp) :: z
    complex(dp), allocatable :: weights(:)
    complex(dp), allocatable :: g(:)
    complex(dp), allocatable :: rhs(:)

    ! Iterators
    integer :: j, k

    ! For error messages
    character(100) :: msg

! ------------------------------------------------------------------------------
! Input processing
! ------------------------------------------------------------------------------
    
    if (a <= 0.0_dp) then
        write(*,*) 'invert3: a must be positive.'
        stat = 1
        return
    end if
    if (order <= -1.0_dp .or. order >= 1.0_dp) then
        write(*,*) 'invert3: order must satisfy -1 < order < 1.'
        stat = 2
        return
    end if
    nt = size(times)

! ------------------------------------------------------------------------------
! Generate contour
! ------------------------------------------------------------------------------

    ! Derived contour parameters
    incline = 0.5 * (min(pi / (1.0_dp + order), pi) - pi/2.0_dp)
    r = 0.9_dp * incline
    gam = (1.0_dp + order) * sigma
    kappa = 1.0_dp - sin(incline - r)
    lambda = gam / (kappa * t2)
    z_spacing = sqrt(2.0_dp * pi * r / (gam * nz))

    ! Allocate pts and derivs
    allocate(z_pts(0:nz), z_derivs(0:nz), stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'contour: allocate: ', msg
        return
    end if

    ! Compute pts and derivs
    do j = 0, nz
        z_pts(j) = lambda*(1.0_dp - sin(incline - imagunit*z_spacing*j))
        z_derivs(j) = imagunit*lambda*cos(incline - imagunit*z_spacing*j)
    end do

    ! Multiply derivs by the constants
    z_derivs = z_spacing/pi * z_derivs
    z_derivs(0) = z_derivs(0) / 2

! ------------------------------------------------------------------------------
! Solve
! ------------------------------------------------------------------------------

    ! Allocate solution
    if (allocated(soln)) deallocate(soln, stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'invert3: deallocate: ', msg
        return
    end if
    allocate(soln(nt), stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'invert3: allocate: ', msg
        return
    end if

    ! Initialize to the integral term
    call indefint_fn_trap_real(f_fn, ic, [0.0_dp, (times(j), j=1,nt)], dt, soln, stat)
    if (stat /= 0) return

    ! Allocate weights, g, rhs
    allocate(g(nt), rhs(nt), stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'invert3: allocate: ', msg
        return
    end if

    ! Compute the Laplace transform part
    do j = 0, nz

        z = z_pts(j)

        lhs = z**(order+1.0_dp) + a

        ! Compute g(z,t) = e^(zt) u0 + int_0^t e^(z(t-s)) f(s) ds.
        ! Note g(z,t2) = e^(z(t2-t1)) * g(z,t1) + int_t1^t2 e^(z(t2-s)) f(s) dx
        t = times(1)
        g(1) = exp(z*times(1)) * ic + defint_fn_trap_cmplx(expf, 0.0_dp, t, dt)
        do k = 2, nt
            t = times(k)
            tm1 = times(k-1)
            g(k) = exp(z*(t-tm1)) * g(k-1) + defint_fn_trap_cmplx(expf, tm1, t, dt)
        end do

        rhs = z**order * g

        soln = soln + aimag( z_derivs(j) * (rhs / lhs - g / z) )

    end do
    
    ! Clean up
    deallocate(z_pts, z_derivs, g, rhs, stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'invert3: deallocate: ', msg
        return
    end if

    ! Successful exit
    stat = 0
    return    

    contains
        complex(dp) function expf(s) result(val)
            real(dp), intent(in) :: s
            val = exp(z*(t-s)) * f_fn(s)
        end function expf

end procedure invert3_0d

end submodule inversion_invert3_0d
