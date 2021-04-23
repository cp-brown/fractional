submodule (inversion_1d) inversion_invert1_1d
implicit none

contains

module procedure invert1_1d
! ------------------------------------------------------------------------------
! Declarations for local variables
! ------------------------------------------------------------------------------

    ! --- Problem size --- !

    ! Number of interior spatial points
    integer :: nx

    ! Number of times to solve at
    integer :: nt


    ! --- For computing the contour --- !

    real(dp) :: r, b, lambda, z_spacing
    complex(dp), allocatable :: z_pts(:), z_derivs(:)


    ! --- For the lhs matrix --- !

    complex(dp), allocatable :: lhs_vals(:)

    ! Array used by mkl_dss, to be computed from ptrb and ptre
    integer, allocatable :: lhs_rowidx(:)


    ! --- For the rhs vector --- !

    ! Array to hold the number of interior integration points in each element
    integer, allocatable :: n_rhspts(:)

    ! Total number of evaluation points needed to compute the rhs vector
    integer :: n_rhspts_total

    ! Locations of the evaluation points needed
    real(dp), allocatable :: rhspts(:)

    ! Array to hold the values of the evaluations
    complex(dp), allocatable :: rhsevals(:)

    ! Array to hold the final rhs vector
    complex(dp), allocatable :: rhs(:)


    ! --- For the solver --- !

    ! Handle for the MKL direct sparse solver
    type(MKL_DSS_HANDLE) :: dss_handle

    ! To temporarily store the linear system solution
    complex(dp), allocatable :: solnvec(:)

    ! 
    integer, dimension(1) :: dss_reorder_perm = [0]


    ! --- Miscellaneous --- !

    ! For the constant and z' parts of the summation
    complex(dp), allocatable :: weights(:)

    ! For all parts of the summation except the exponential term
    complex(dp), allocatable :: uhat(:,:)

    ! Iterators
    integer :: j, k

    ! For error messages
    character(100) :: msg

! ------------------------------------------------------------------------------
! Input processing
! ------------------------------------------------------------------------------

    ! Deallocate soln if already allocated
    if (allocated(soln)) deallocate(soln, stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'invert1: deallocate: ', msg
        return
    end if

    ! Get problem sizes
    nx = size(x_pts) - 2
    nt = size(times)

    ! Verify Amat is the correct size
    if (size(Amat_rowptrb) /= nx .or. size(Amat_rowptre) /= nx) then
        write(*,*) 'invert1: unexpected size for Amat.'
        return
    end if
    
    ! Verify -1 < order < 1
    if (order <= -1.0_dp .or. order >= 1.0_dp) then
        write(*,*) 'invert1: order must satisfy -1 < order < 1.'
        stat = -1
        return
    end if

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
    
! ------------------------------------------------------------------------------
! Matrix and solver setup
! ------------------------------------------------------------------------------

    ! --- Generate lhs matrix structure for the solver setup --- !

    ! Allocate lhs_rowidx, lhs_vals
    allocate(lhs_rowidx(nx+1), lhs_vals(Amat_rowptre(nx)-1), stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'invert1: allocate: ', msg
        return
    end if

    ! Compute lhs_rowidx
    lhs_rowidx(1:nx) = Amat_rowptrb
    lhs_rowidx(nx+1) = Amat_rowptre(nx)


    ! --- Solver setup and re-ordering --- !

    ! Initialize solver
    stat = dss_create(dss_handle, MKL_DSS_DEFAULTS)
    if (stat /= 0) then
        write(*,*) 'invert1: dss_create: exited with code ', stat
        return
    end if

    ! Define structure
    stat = dss_define_structure(dss_handle, MKL_DSS_SYMMETRIC_COMPLEX, &
        lhs_rowidx, nx, nx, Amat_cols, lhs_rowidx(nx+1)-1)
    if (stat /= 0) then
        write(*,*) 'invert1: dss_define_structure: exited with code ', stat
        return
    end if

    ! Perform re-ordering
    stat = dss_reorder(dss_handle, MKL_DSS_DEFAULTS, dss_reorder_perm)
    if (stat /= 0) then
        write(*,*) 'invert1: dss_reorder: exited with code ', stat
        return
    end if


    ! --- Allocate and compute arrays needed for the rhs --- !

    call fem_vec_evalpts_1d(x_pts, dx, n_rhspts_total, n_rhspts, rhspts, stat)
    if (stat /= 0) return

    allocate(rhsevals(n_rhspts_total), rhs(nx), stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'invert1: allocate: ', msg
        return
    end if


    ! --- Allocate and compute weights --- !

    allocate(weights(0:nz), stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'invert1: allocate: ', msg
        return
    end if

    weights = z_spacing/pi * z_derivs
    weights(0) = weights(0) / 2


    ! --- Other allocation --- !

    allocate(uhat(nx, 0:nz), solnvec(nx), stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'invert1: allocate: ', msg
        return
    end if

! ------------------------------------------------------------------------------
! Compute Laplace-space terms
! ------------------------------------------------------------------------------

    do j = 0, nz

        ! --- Generate the lhs matrix entries --- !

        ! Compute the sum
        lhs_vals = z_pts(j)**(order+1.0_dp) * Imat_vals + Amat_vals


        ! --- Generate the rhs vector --- !

        ! Compute the evals needed for the FEM vector
        if (homogeneous) then
            do k = 1, n_rhspts_total
                rhsevals(k) = ic_fn(rhspts(k))
            end do 
        else
            do k = 1, n_rhspts_total
                rhsevals(k) = ic_fn(rhspts(k)) + fhat_fn(rhspts(k), z_pts(j))
            end do 
        end if

        ! Compute the FEM vector
        call fem_vec_cmplx_1d(n_rhspts, rhspts, rhsevals, rhs, stat)
        if (stat /= 0) return
        rhs = z_pts(j)**order * rhs


        ! --- Solve and save --- !
    
        ! Factor the system
        stat = dss_factor_complex(dss_handle, MKL_DSS_INDEFINITE, lhs_vals)
        if (stat /= 0) then
            write(*,*) 'invert1: dss_factor_complex: exited with code ', stat
            return
        end if

        ! Solve the system
        stat = dss_solve_complex(dss_handle, MKL_DSS_DEFAULTS, rhs, 1, solnvec)
        if (stat /= 0) then
            write(*,*) 'invert1: dss_solve_complex: exited with code ', stat
            return
        end if

        ! Save to uhat, multiply by weight
        uhat(:,j) = weights(j) * solnvec

    end do

! ------------------------------------------------------------------------------
! Solver clean-up
! ------------------------------------------------------------------------------

    ! Uninitialize solver
    stat = dss_delete(dss_handle, MKL_DSS_DEFAULTS)
    if (stat /= 0) then
        write(*,*) 'invert1: dss_delete: exited with code ', stat
        return
    end if

    ! rhs deallocations
    deallocate(n_rhspts, rhspts, rhsevals, rhs, stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'invert1: deallocate: ', msg
        return
    end if

    ! Other deallocations
    deallocate(z_derivs, lhs_rowidx, lhs_vals, &
        weights, stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'invert1: deallocate: ', msg
        return
    end if

! ------------------------------------------------------------------------------
! Compute solution
! ------------------------------------------------------------------------------

    ! Allocate solution
    allocate(soln(nx, nt), stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'invert1: allocate: ', msg
        return
    end if

    ! Compute solution
    do j = 1, nt
        call gemv(uhat, exp(times(j)*z_pts), solnvec)
        soln(:,j) = aimag(solnvec)
    end do

! ------------------------------------------------------------------------------
! Final clean-up, exit
! ------------------------------------------------------------------------------
    
    deallocate(z_pts, uhat, solnvec, stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'invert1: deallocate: ', msg
        return
    end if

    ! Successful exit
    stat = 0
    return    

end procedure invert1_1d

end submodule inversion_invert1_1d
