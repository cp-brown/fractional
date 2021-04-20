submodule (inversion_0d) inversion_invert2_0d
implicit none

contains

module procedure invert2_0d
! ------------------------------------------------------------------------------
! Declarations for local variables
! ------------------------------------------------------------------------------

    ! For computing the contour
    real(dp) :: r, gam, kappa, lambda, z_spacing
    complex(dp), allocatable :: z_pts(:), z_derivs(:)

    integer :: nt
    complex(dp) :: uhat
    complex(dp) :: lhs
    complex(dp) :: g
    complex(dp) :: rhs
    real(dp), allocatable :: solnpt(:)

    ! For MPI and work distribution
    integer :: rank, num_cores
    integer :: ntasks, rem
    integer :: zfirst, zlast

    ! Iterator
    integer :: j
    
    ! For error messages
    character(100) :: msg

! ------------------------------------------------------------------------------
! Input processing
! ------------------------------------------------------------------------------

    if (a <= 0.0_dp) then
        write(*,*) 'invert2: a must be positive.'
        stat = 1
        return
    end if
    if (order <= -1.0_dp .or. order >= 1.0_dp) then
        write(*,*) 'invert2: order must satisfy -1 < order < 1.'
        stat = 2
        return
    end if
    nt = size(times)

! ------------------------------------------------------------------------------
! Distribute work
! ------------------------------------------------------------------------------

    call MPI_Comm_rank(comm, rank, stat)
    call MPI_Comm_size(comm, num_cores, stat)

    ntasks = (nz+2) / num_cores
    rem = mod((nz+2), num_cores)
    
    if (rank <= rem) then
        zfirst = rank * (ntasks+1)
    else
        zfirst = rem * (ntasks+1) + (rank-rem) * ntasks
    end if
    
    if (rank < rem) then
        zlast = zfirst + ntasks
    else
        zlast = zfirst + ntasks - 1
    end if

    if (rank == num_cores-1) zlast = zlast - 1

! ------------------------------------------------------------------------------
! Generate contour
! ------------------------------------------------------------------------------

    ! Derived contour parameters
    r = 0.9_dp * incline
    gam = (1.0_dp + order) * sigma
    kappa = 1.0_dp - sin(incline - r)
    lambda = gam / (kappa * t2)
    z_spacing = sqrt(2.0_dp * pi * r / (gam * nz))

    ! Allocate pts and derivs
    allocate(z_pts(zfirst:zlast), z_derivs(zfirst:zlast), stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'contour: allocate: ', msg
        return
    end if

    ! Compute pts and derivs
    do j = zfirst, zlast
        z_pts(j) = focus + lambda*(1.0_dp - sin(incline - imagunit*z_spacing*j))
        z_derivs(j) = imagunit*lambda*cos(incline - imagunit*z_spacing*j)
    end do

    ! Multiply derivs by the constants
    z_derivs = z_spacing/pi * z_derivs
    if (zfirst == 0) then
        z_derivs(0) = z_derivs(0) / 2
    end if

! ------------------------------------------------------------------------------
! Solve
! ------------------------------------------------------------------------------

    ! Allocate partial solutions
    allocate(solnpt(nt), stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'invert2: allocate: ', msg
        return
    end if

    ! Last core gets to compute integral term, others initialize to zero
    if (rank == num_cores-1) then
        call indefint_fn_trap_real(f_fn, ic, [0.0_dp, (times(j), j=1,nt)], dt, solnpt, stat)
        if (stat /= 0) return
    else
        solnpt = 0.0_dp
    end if

    ! Add respective Laplace space terms
    do j = zfirst, zlast

        lhs = z_pts(j)**(order+1.0_dp) + a
        g = ic + fhat_fn(z_pts(j))
        rhs = z_pts(j)**order * g
        uhat = z_derivs(j) * ( rhs / lhs - g / z_pts(j) )

        solnpt = solnpt + aimag( uhat * exp(times * z_pts(j)) )

    end do

    ! User-specified solncore computes final solution by adding all partial solutions
    if (rank == solncore) then

        ! Allocate solution
        if (allocated(soln)) deallocate(soln, stat=stat, errmsg=msg)
        if (stat /= 0) then
            write(*,*) 'invert2: deallocate: ', msg
            return
        end if
        allocate(soln(nt), stat=stat, errmsg=msg)
        if (stat /= 0) then
            write(*,*) 'invert2: allocate: ', msg
            return
        end if

    end if
    call MPI_Reduce(solnpt, soln, nt, MPI_DOUBLE_PRECISION, MPI_SUM, solncore, comm, stat)
    
    ! Clean up
    deallocate(z_pts, z_derivs, solnpt, stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'invert2: deallocate: ', msg
        return
    end if

    ! Successful exit
    stat = 0
    return    

end procedure invert2_0d

end submodule inversion_invert2_0d
