submodule (inversion_1d) inversion_invert3_1d
implicit none

contains

module procedure invert3_1d
! ------------------------------------------------------------------------------
! Declarations for local variables
! ------------------------------------------------------------------------------

    ! --- Problem size --- !

    ! Number of interior spatial points
    integer :: nx

    ! Number of times to solve at
    integer :: nt


    ! --- For computing the contour --- !

    real(dp) :: incline, r, gam, kappa, lambda, z_spacing
    complex(dp), allocatable :: z_pts(:), z_derivs(:)


    ! --- For the Laplacian matrix Amat and identity matrix Imat --- !

    ! For converting to complex format
    complex(C_DOUBLE_COMPLEX), allocatable :: Amat_vals_cmplx(:)
    complex(C_DOUBLE_COMPLEX), allocatable :: Imat_vals_cmplx(:)

    ! MKL sparse BLAS handles
    type(SPARSE_MATRIX_T) :: Amat
    type(SPARSE_MATRIX_T) :: Imat


    ! --- For the lhs matrix --- !

    ! MKL sparse BLAS handle
    type(SPARSE_MATRIX_T) :: lhs

    ! Variables needed for mkl_sparse_z_export_csr, otherwise unused
    integer(C_INT) :: lhs_idxing
    integer(C_INT) :: lhs_nrows, lhs_ncols

    ! C pointers for exporting
    type(c_ptr) :: lhs_rowptrb_c, lhs_rowptre_c
    type(c_ptr) :: lhs_cols_c
    type(c_ptr) :: lhs_vals_c

    ! Fortran pointers for exporting
    integer(C_INT), pointer :: lhs_rowptrb(:), lhs_rowptre(:)
    integer(C_INT), pointer :: lhs_cols(:)
    complex(C_DOUBLE_COMPLEX), pointer :: lhs_vals(:)

    ! Number of non-zero elements
    integer :: lhs_nnz

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

    ! Array to hold the correction vector gvec
    complex(dp), allocatable :: gvec(:)


    ! --- For the solver --- !

    ! Handle for the MKL direct sparse solver
    type(MKL_DSS_HANDLE) :: dss_handle

    ! To temporarily store the linear system solution
    complex(dp), allocatable :: solnvec(:)


    ! --- Miscellaneous --- !

    ! For the constant and z' parts of the summation
    complex(dp), allocatable :: weights(:)

    ! To store the current locations
    real(dp) :: x, t, tm1
    complex(dp) :: z

    ! Iterators
    integer :: j, k, l

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

    ! Verify Amat and Imat are the correct size
    if (size(Amat_rowptrb) /= nx .or. size(Amat_rowptre) /= nx) then
        write(*,*) 'invert1: unexpected size for Amat.'
        return
    end if
    if (size(Imat_rowptrb) /= nx .or. size(Imat_rowptre) /= nx) then
        write(*,*) 'invert1: unexpected size for Imat.'
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
    allocate(soln(nx, nt), stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'invert1: allocate: ', msg
        return
    end if

    ! Initialize solution as the integral part
    do j = 1, nx

        ! This is needed so that the local function fx_fn knows the correct
        ! spatial location
        x = x_pts(j+1)

        call indefint_fn_trap_real(fx_fn, ic_fn(x), [0.0_dp, (times(k), &
            k=1,nt)], dt, soln(j,:), stat)
        if (stat /= 0) return

    end do
    
! ------------------------------------------------------------------------------
! Matrix and solver setup
! ------------------------------------------------------------------------------

    ! --- Given matrix processing --- !

    ! Convert Amat, Imat values to complex format
    allocate(Amat_vals_cmplx(size(Amat_vals)), Imat_vals_cmplx(size(Imat_vals)), &
        stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'invert1: allocate: ', msg
        return
    end if
    Amat_vals_cmplx = cmplx(Amat_vals, kind=C_DOUBLE_COMPLEX)
    Imat_vals_cmplx = cmplx(Imat_vals, kind=C_DOUBLE_COMPLEX)

    ! Export Amat to MKL sparse format
    stat = mkl_sparse_z_create_csr(Amat, SPARSE_INDEX_BASE_ONE, nx, nx, &
        Amat_rowptrb, Amat_rowptre, Amat_cols, Amat_vals_cmplx)
    if (stat /= 0) then
        write(*,*) 'invert1: mkl_sparse_d_create_csr: exited with code ', stat
        return
    end if
    
    ! Export Imat to MKL sparse format
    stat = mkl_sparse_z_create_csr(Imat, SPARSE_INDEX_BASE_ONE, nx, nx, &
        Imat_rowptrb, Imat_rowptre, Imat_cols, Imat_vals_cmplx)
    if (stat /= 0) then
        write(*,*) 'invert1: mkl_sparse_d_create_csr: exited with code ', stat
        return
    end if


    ! --- Generate lhs matrix structure for the solver setup --- !

    ! Compute the sum - the values will be ignored
    stat = mkl_sparse_z_add(SPARSE_OPERATION_NON_TRANSPOSE, Imat, &
        (1.0_dp,0.0_dp), Amat, lhs)
    if (stat /= 0) then
        write(*,*) 'invert1: mkl_sparse_z_add: exited with code ', stat
        return
    end if

    ! Export to c_ptr arrays
    stat = mkl_sparse_z_export_csr(lhs, lhs_idxing, lhs_nrows, lhs_ncols, &
        lhs_rowptrb_c, lhs_rowptre_c, lhs_cols_c, lhs_vals_c)
    if (stat /= 0) then
        write(*,*) 'invert1: mkl_sparse_z_export_csr: exited with code ', stat
        return
    end if

    ! Convert c pointers to fortran pointers - structure arrays only
    call c_f_pointer(lhs_rowptrb_c, lhs_rowptrb, [nx])
    call c_f_pointer(lhs_rowptre_c, lhs_rowptre, [nx])
    lhs_nnz = lhs_rowptre(nx) - 1
    call c_f_pointer(lhs_cols_c, lhs_cols, [lhs_nnz])

    ! Allocate lhs_rowidx
    allocate(lhs_rowidx(nx+1), stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'invert1: allocate: ', msg
        return
    end if

    ! Compute lhs_rowidx
    lhs_rowidx(1:nx) = lhs_rowptrb
    lhs_rowidx(nx+1) = lhs_rowptre(nx)


    ! --- Solver setup and re-ordering --- !

    ! Initialize solver
    stat = dss_create(dss_handle, MKL_DSS_DEFAULTS)
    if (stat /= 0) then
        write(*,*) 'invert1: dss_create: exited with code ', stat
        return
    end if

    ! Define structure
    stat = dss_define_structure(dss_handle, MKL_DSS_SYMMETRIC_COMPLEX, &
        lhs_rowidx, nx, nx, lhs_cols, lhs_nnz)
    if (stat /= 0) then
        write(*,*) 'invert1: dss_define_structure: exited with code ', stat
        return
    end if

    ! Perform re-ordering
    stat = dss_reorder(dss_handle, MKL_DSS_DEFAULTS, [0])
    if (stat /= 0) then
        write(*,*) 'invert1: dss_reorder: exited with code ', stat
        return
    end if


    ! --- Allocate and compute arrays needed for the rhs --- !

    call fem_vec_evalpts_1d(x_pts, dx, n_rhspts_total, n_rhspts, rhspts, stat)
    if (stat /= 0) return

    allocate(rhsevals(n_rhspts_total), rhs(nx), gvec(nx), stat=stat, errmsg=msg)
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

    allocate(solnvec(nx), stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'invert1: allocate: ', msg
        return
    end if

! ------------------------------------------------------------------------------
! Add Laplace-space terms
! ------------------------------------------------------------------------------

    do j = 0, nz

        z = z_pts(j)

        ! --- Generate the lhs matrix entries, factor --- !

        ! Compute the sum
        stat = mkl_sparse_z_add(SPARSE_OPERATION_NON_TRANSPOSE, Imat, &
            z**(order+1.0_dp), Amat, lhs)
        if (stat /= 0) then
            write(*,*) 'invert1: mkl_sparse_z_add: exited with code ', stat
            return
        end if

        ! Export to c_ptr arrays
        stat = mkl_sparse_z_export_csr(lhs, lhs_idxing, lhs_nrows, lhs_ncols, &
            lhs_rowptrb_c, lhs_rowptre_c, lhs_cols_c, lhs_vals_c)
        if (stat /= 0) then
            write(*,*) 'invert1: mkl_sparse_z_export_csr: exited with code ', stat
            return
        end if

        ! Convert c pointer to fortran pointer - value array only
        call c_f_pointer(lhs_vals_c, lhs_vals, [lhs_nnz])

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
            do k = 1, n_rhspts_total
                x = rhspts(k)
                rhsevals(k) = exp(z*t) * ic_fn(x)
            end do 
        else
            do k = 1, n_rhspts_total
                x = rhspts(k)
                rhsevals(k) = exp(z*t) * ic_fn(x) &
                + defint_fn_trap_cmplx(expf_fn, 0.0_dp, t, dt)
            end do 
        end if

        ! Compute the FEM vector
        call fem_vec_cmplx_1d(n_rhspts, rhspts, rhsevals, rhs, stat)
        if (stat /= 0) return
        rhs = z**order * rhs

        ! Compute gvec
        if (homogeneous) then
            do k = 1, nx
                x = x_pts(k+1)
                gvec(k) = exp(z*t) * ic_fn(x)
            end do
        else
            do k = 1, nx
                x = x_pts(k+1)
                gvec(k) = exp(z*t) * ic_fn(x) &
                + defint_fn_trap_cmplx(expf_fn, 0.0_dp, t, dt)
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
        soln(:,1) = soln(:,1) + aimag( weights(j) * (solnvec - gvec/z) )


        ! --- Subsequent times --- !

        do l = 2, nt

            tm1 = times(l-1)
            t = times(l)

            ! Compute the evals needed for the FEM vector
            if (homogeneous) then
                do k = 1, n_rhspts_total
                    x = rhspts(k)
                    rhsevals(k) = exp(z*(t-tm1)) * rhsevals(k)
                end do
            else
                do k = 1, n_rhspts_total
                    x = rhspts(k)
                    rhsevals(k) = exp(z*(t-tm1)) * rhsevals(k) &
                    + defint_fn_trap_cmplx(expf_fn, tm1, t, dt)
                end do
            end if

            ! Compute the FEM vector
            call fem_vec_cmplx_1d(n_rhspts, rhspts, rhsevals, rhs, stat)
            if (stat /= 0) return
            rhs = z**order * rhs

            ! Compute gvec
            if (homogeneous) then
                do k = 1, nx
                    x = x_pts(k+1)
                    gvec(k) = exp(z*(t-tm1)) * gvec(k)
                end do
            else
                do k = 1, nx
                    x = x_pts(k+1)
                    gvec(k) = exp(z*(t-tm1)) * gvec(k) &
                    + defint_fn_trap_cmplx(expf_fn, tm1, t, dt)
                end do
            end if

            ! Solve the system
            stat = dss_solve_complex(dss_handle, MKL_DSS_DEFAULTS, rhs, 1, solnvec)
            if (stat /= 0) then
                write(*,*) 'invert1: dss_solve_complex: exited with code ', stat
                return
            end if

            ! Apply correction, save to uhat, multiply by weight
            soln(:,l) = soln(:,l) + aimag( weights(j) * (solnvec - gvec/z) )

        end do

    end do

! ------------------------------------------------------------------------------
! Clean-up, exit
! ------------------------------------------------------------------------------

    ! Uninitialize solver
    stat = dss_delete(dss_handle, MKL_DSS_DEFAULTS)
    if (stat /= 0) then
        write(*,*) 'invert1: dss_delete: exited with code ', stat
        return
    end if

    ! Delete lhs MKL object
    stat = mkl_sparse_destroy(lhs)
    if (stat /= 0) then
        write(*,*) 'invert1: mkl_sparse_destroy: exited with code ', stat
        return
    end if

    ! rhs deallocations
    deallocate(n_rhspts, rhspts, rhsevals, rhs, gvec, stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'invert1: deallocate: ', msg
        return
    end if

    ! Other deallocations
    deallocate(z_derivs, Amat_vals_cmplx, Imat_vals_cmplx, lhs_rowidx, &
        weights, stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'invert1: deallocate: ', msg
        return
    end if
    
    deallocate(solnvec, stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'invert1: deallocate: ', msg
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
        val = f_fn(x,t)
    end function fx_fn

    !
    complex(dp) function expf_fn(s) result(val)
        real(dp), intent(in) :: s
        val = exp(z*(t-s)) * f_fn(x,s)
    end function expf_fn

end procedure invert3_1d

end submodule inversion_invert3_1d
