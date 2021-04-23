submodule (inversion_0d) inversion_invert1_0d
implicit none

contains

module procedure invert1_0d
! ------------------------------------------------------------------------------
! Declarations for local variables
! ------------------------------------------------------------------------------

    ! For computing the contour
    real(dp) :: r, b, lambda, z_spacing
    complex(dp), allocatable :: z_pts(:), z_derivs(:)

    ! 
    integer :: nt
    complex(dp) :: uhat
    complex(dp) :: lhs
    complex(dp) :: rhs

    ! Iterator
    integer :: j

    ! For error messages
    character(100) :: msg

! ------------------------------------------------------------------------------
! Input processing
! ------------------------------------------------------------------------------
    
    if (a <= 0.0_dp) then
        write(*,*) 'invert1: a must be positive.'
        stat = 1
        return
    end if
    if (order <= -1.0_dp .or. order >= 1.0_dp) then
        write(*,*) 'invert1: order must satisfy -1 < order < 1.'
        stat = 1
        return
    end if
    if (theta <= 0.0_dp .or. theta >= 1.0_dp) then
        write(*,*) 'invert1: order must satisfy 0 < theta < 1.'
        stat = 1
        return
    end if
    nt = size(times)

! ------------------------------------------------------------------------------
! Generate contour
! ------------------------------------------------------------------------------

    ! Derived contour parameters
    r = 0.9_dp * incline
    b = acosh((theta * t1 / t2 * sin(incline))**(-1.0_dp))
    lambda = 2.0_dp * pi * theta * r * nz / (b * t2)
    z_spacing = b / nz

    ! Allocate pts and derivs
    allocate(z_pts(0:nz), z_derivs(0:nz), stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'contour: allocate: ', msg
        return
    end if

    ! Compute pts and derivs
    do j = 0, nz
        z_pts(j) = focus + lambda*(1.0_dp - sin(incline - imagunit*z_spacing*j))
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
        write(*,*) 'invert1: deallocate: ', msg
        return
    end if
    allocate(soln(nt), stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'invert1: allocate: ', msg
        return
    end if
    soln = 0.0_dp

    ! Add respective Laplace space terms
    do j = 0, nz

        lhs = z_pts(j)**(order+1.0_dp) + a
        rhs = z_pts(j)**order * (ic + fhat_fn(z_pts(j)))
        uhat = z_derivs(j) * rhs / lhs

        soln = soln + aimag( uhat * exp(times * z_pts(j)) )

    end do
    
    ! Clean up
    deallocate(z_pts, z_derivs, stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'invert1: deallocate: ', msg
        return
    end if

    ! Successful exit
    stat = 0
    return    

end procedure invert1_0d

end submodule inversion_invert1_0d
