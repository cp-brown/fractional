submodule (inversion_2d) inversion_invert3_2d
implicit none

contains

module procedure invert3_2d
! ------------------------------------------------------------------------------
! Declarations for local variables
! ------------------------------------------------------------------------------

    ! --- Problem size --- !

    ! Number of interior spatial points
    integer :: nx, ny

    ! Number of times to solve at
    integer :: nt


    ! --- For computing the contour --- !

    real(dp) :: incline, r, gam, kappa, lambda, z_spacing
    complex(dp), allocatable :: z_pts(:), z_derivs(:)


    ! --- For the lhs matrix --- !

    complex(dp), allocatable :: lhs_vals(:)

    ! Array used by mkl_dss, to be computed from ptrb and ptre
    integer, allocatable :: lhs_rowidx(:)


    ! --- For the rhs vector --- !

    ! Array to hold the number of interior integration points in each element
    integer, allocatable :: n_rhspts_x(:), n_rhspts_y(:)

    ! Total number of evaluation points needed to compute the rhs vector
    integer :: n_rhspts_total_x, n_rhspts_total_y

    ! Locations of the evaluation points needed
    real(dp), allocatable :: rhspts_x(:), rhspts_y(:)

    ! Array to hold the values of the evaluations
    complex(dp), allocatable :: rhsevals(:,:)

    ! Array to hold the final rhs vector
    complex(dp), allocatable :: rhs(:)

    ! Array to hold the correction vector gvec
    complex(dp), allocatable :: gvec(:)


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

    ! To store the current locations
    real(dp) :: x, y, t, tm1
    complex(dp) :: z

    ! Iterators
    integer :: j, k, l, p

    ! For error messages
    character(100) :: msg

! ------------------------------------------------------------------------------
! Input processing
! ------------------------------------------------------------------------------
    
    ! Deallocate soln if already allocated
    if (allocated(soln)) deallocate(soln, stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'invert3: deallocate: ', msg
        return
    end if

    ! Get problem sizes
    nx = size(x_pts) - 2
    ny = size(y_pts) - 2
    nt = size(times)

    ! Verify Amat and Imat are the correct size
    if (size(Amat_rowptrb) /= nx*ny .or. size(Amat_rowptre) /= nx*ny) then
        write(*,*) 'invert3: unexpected size for Amat.'
        return
    end if
    if (size(Imat_rowptrb) /= nx*ny .or. size(Imat_rowptre) /= nx*ny) then
        write(*,*) 'invert3: unexpected size for Imat.'
        return
    end if
    
    ! Verify -1 < order < 1
    if (order <= -1.0_dp .or. order >= 1.0_dp) then
        write(*,*) 'invert3: order must satisfy -1 < order < 1.'
        stat = -1
        return
    end if

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

! ------------------------------------------------------------------------------
! Compute the integral part of the solution
! ------------------------------------------------------------------------------

    ! Allocate solution
    allocate(soln(nx, ny, nt), stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'invert3: allocate: ', msg
        return
    end if

    ! Initialize solution as the integral part
    do j = 1, nx
        do l = 1, ny

        ! This is needed so that the local function fx_fn knows the correct
        ! spatial location
        x = x_pts(j+1)
        y = y_pts(l+1)

        call indefint_fn_trap_real(fx_fn, ic_fn(x,y), [0.0_dp, (times(k), &
            k=1,nt)], dt, soln(j,l,:), stat)
        if (stat /= 0) return

        end do

    end do
    
! ------------------------------------------------------------------------------
! Matrix and solver setup
! ------------------------------------------------------------------------------

    ! --- Generate lhs matrix structure for the solver setup --- !

    ! Allocate lhs_rowidx
    allocate(lhs_rowidx(nx*ny+1), lhs_vals(Amat_rowptre(nx*ny)-1), stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'invert3: allocate: ', msg
        return
    end if

    ! Compute lhs_rowidx
    lhs_rowidx(1:nx*ny) = Amat_rowptrb
    lhs_rowidx(nx*ny+1) = Amat_rowptre(nx*ny)


    ! --- Solver setup and re-ordering --- !

    ! Initialize solver
    stat = dss_create(dss_handle, MKL_DSS_DEFAULTS)
    if (stat /= 0) then
        write(*,*) 'invert3: dss_create: exited with code ', stat
        return
    end if

    ! Define structure
    stat = dss_define_structure(dss_handle, MKL_DSS_SYMMETRIC_COMPLEX, &
        lhs_rowidx, nx*ny, nx*ny, Amat_cols, lhs_rowidx(nx*ny+1)-1)
    if (stat /= 0) then
        write(*,*) 'invert3: dss_define_structure: exited with code ', stat
        return
    end if

    ! Perform re-ordering
    stat = dss_reorder(dss_handle, MKL_DSS_DEFAULTS, dss_reorder_perm)
    if (stat /= 0) then
        write(*,*) 'invert3: dss_reorder: exited with code ', stat
        return
    end if


    ! --- Allocate and compute arrays needed for the rhs --- !

    call fem_vec_evalpts_1d(x_pts, dx, n_rhspts_total_x, n_rhspts_x, &
        rhspts_x, stat)
    if (stat /= 0) return
    call fem_vec_evalpts_1d(y_pts, dy, n_rhspts_total_y, n_rhspts_y, &
        rhspts_y, stat)
    if (stat /= 0) return

    allocate(rhsevals(n_rhspts_total_x, n_rhspts_total_y), rhs(nx*ny), gvec(nx*ny), stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'invert3: allocate: ', msg
        return
    end if


    ! --- Allocate and compute weights --- !

    allocate(weights(0:nz), stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'invert3: allocate: ', msg
        return
    end if

    weights = z_spacing/pi * z_derivs
    weights(0) = weights(0) / 2


    ! --- Other allocation --- !

    allocate(solnvec(nx*ny), stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'invert3: allocate: ', msg
        return
    end if

! ------------------------------------------------------------------------------
! Add Laplace-space terms
! ------------------------------------------------------------------------------

    do j = 0, nz

        z = z_pts(j)

        ! --- Generate the lhs matrix entries, factor --- !

        ! Compute the sum
        lhs_vals = z_pts(j)**(order+1.0_dp) * Imat_vals + Amat_vals

        ! Factor the system
        stat = dss_factor_complex(dss_handle, MKL_DSS_INDEFINITE, lhs_vals)
        if (stat /= 0) then
            write(*,*) 'invert1: dss_factor_complex: exited with code ', stat
            return
        end if


        ! --- Generate first rhs, gvec vectors --- !

        t = times(1)

        ! Compute the evals needed for the FEM vector
        if (homogeneous) then
            do k = 1, n_rhspts_total_y
                do l = 1, n_rhspts_total_x
                    x = rhspts_x(l)
                    y = rhspts_y(k)
                    rhsevals(l,k) = exp(z*t) * ic_fn(x,y)
                end do
            end do 
        else
            do k = 1, n_rhspts_total_y
                do l = 1, n_rhspts_total_x
                    x = rhspts_x(l)
                    y = rhspts_y(k)
                    rhsevals(l,k) = exp(z*t) * ic_fn(x,y) &
                    + defint_fn_trap_cmplx(expf_fn, 0.0_dp, t, dt)
                end do
            end do 
        end if

        ! Compute the FEM vector
        call fem_vec_cmplx_2d(n_rhspts_x, n_rhspts_y, rhspts_x, rhspts_y, &
            rhsevals, rhs, stat)
        if (stat /= 0) return
        rhs = z**order * rhs

        ! Compute gvec
        if (homogeneous) then
            do k = 1, nx
                do l = 1, ny
                    x = x_pts(k+1)
                    y = y_pts(l+1)
                    gvec((l-1)*nx + k) = exp(z*t) * ic_fn(x,y)
                end do
            end do
        else
            do k = 1, nx
                do l = 1, ny
                    x = x_pts(k+1)
                    y = y_pts(l+1)
                    gvec((l-1)*nx + k) = exp(z*t) * ic_fn(x,y) &
                    + defint_fn_trap_cmplx(expf_fn, 0.0_dp, t, dt)
                end do
            end do
        end if


        ! --- Solve first system, save --- !

        ! Solve the system
        stat = dss_solve_complex(dss_handle, MKL_DSS_DEFAULTS, rhs, 1, solnvec)
        if (stat /= 0) then
            write(*,*) 'invert1: dss_solve_complex: exited with code ', stat
            return
        end if

        ! Apply correction, save to uhat, multiply by weight
        soln(:,:,1) = soln(:,:,1) + reshape(aimag( weights(j) * (solnvec - gvec/z) ), &
            [nx, ny])


        ! --- Subsequent times --- !

        do l = 2, nt

            tm1 = times(l-1)
            t = times(l)

            ! Compute the evals needed for the FEM vector
            if (homogeneous) then
                do k = 1, n_rhspts_total_y
                    do p = 1, n_rhspts_total_x
                        x = rhspts_x(p)
                        y = rhspts_y(k)
                        rhsevals(p,k) = exp(z*(t-tm1)) * rhsevals(p,k)
                    end do
                end do
            else
                do k = 1, n_rhspts_total_y
                    do p = 1, n_rhspts_total_x
                        x = rhspts_x(p)
                        y = rhspts_y(k)
                        rhsevals(p,k) = exp(z*(t-tm1)) * rhsevals(p,k) &
                        + defint_fn_trap_cmplx(expf_fn, tm1, t, dt)
                    end do
                end do
            end if

            ! Compute the FEM vector
            call fem_vec_cmplx_2d(n_rhspts_x, n_rhspts_y, rhspts_x, rhspts_y, &
                rhsevals, rhs, stat)
            if (stat /= 0) return
            rhs = z**order * rhs

            ! Compute gvec
            if (homogeneous) then
                do k = 1, nx
                    do p = 1, ny
                        x = x_pts(k+1)
                        y = y_pts(p+1)
                        gvec((p-1)*nx+k) = exp(z*(t-tm1)) * gvec((p-1)*nx+k)
                    end do
                end do
            else
                do k = 1, nx
                    do p = 1, ny
                        x = x_pts(k+1)
                        y = y_pts(p+1)
                        gvec((p-1)*nx+k) = exp(z*(t-tm1)) * gvec((p-1)*nx+k) &
                        + defint_fn_trap_cmplx(expf_fn, tm1, t, dt)
                    end do
                end do
            end if

            ! Solve the system
            stat = dss_solve_complex(dss_handle, MKL_DSS_DEFAULTS, rhs, 1, solnvec)
            if (stat /= 0) then
                write(*,*) 'invert1: dss_solve_complex: exited with code ', stat
                return
            end if

            ! Apply correction, save to uhat, multiply by weight
            soln(:,:,l) = soln(:,:,l) + reshape(aimag( weights(j) * (solnvec &
                - gvec/z) ), [nx,ny])

        end do

    end do

! ------------------------------------------------------------------------------
! Clean-up, exit
! ------------------------------------------------------------------------------

    ! Uninitialize solver
    stat = dss_delete(dss_handle, MKL_DSS_DEFAULTS)
    if (stat /= 0) then
        write(*,*) 'invert3: dss_delete: exited with code ', stat
        return
    end if

    ! rhs deallocations
    deallocate(n_rhspts_x, n_rhspts_y, rhspts_x, rhspts_y, rhsevals, &
        rhs, gvec, stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'invert3: deallocate: ', msg
        return
    end if

    ! Other deallocations
    deallocate(z_derivs, lhs_rowidx, lhs_vals, &
        weights, stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'invert3: deallocate: ', msg
        return
    end if
    
    deallocate(solnvec, stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'invert3: deallocate: ', msg
        return
    end if

    ! Successful exit
    stat = 0
    return

! ------------------------------------------------------------------------------
! Local functions
! ------------------------------------------------------------------------------

    contains

    ! For the integral part of the solution
    real(dp) function fx_fn(t) result(val)
        real(dp), intent(in) :: t
        val = f_fn(x,y,t)
    end function fx_fn

    !
    complex(dp) function expf_fn(s) result(val)
        real(dp), intent(in) :: s
        val = exp(z*(t-s)) * f_fn(x,y,s)
    end function expf_fn

end procedure invert3_2d

end submodule inversion_invert3_2d
