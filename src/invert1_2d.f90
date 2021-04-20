submodule (inversion_2d) inversion_invert1_2d
implicit none

contains

module procedure invert1_2d
! ------------------------------------------------------------------------------
! Declarations for local variables
! ------------------------------------------------------------------------------

    ! --- Problem size --- !

    ! Number of interior spatial points
    integer :: nx, ny

    ! Number of times to solve at
    integer :: nt


    ! --- For computing the contour --- !

    real(dp) :: r, b, lambda, z_spacing
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
    integer, allocatable :: n_rhspts_x(:), n_rhspts_y(:)

    ! Total number of evaluation points needed to compute the rhs vector
    integer :: n_rhspts_total_x, n_rhspts_total_y

    ! Locations of the evaluation points needed
    real(dp), allocatable :: rhspts_x(:), rhspts_y(:)

    ! Array to hold the values of the evaluations
    complex(dp), allocatable :: rhsevals(:,:)

    ! Array to hold the final rhs vector
    complex(dp), allocatable :: rhs(:)


    ! --- For the solver --- !

    ! Handle for the MKL direct sparse solver
    type(MKL_DSS_HANDLE) :: dss_handle

    ! To temporarily store the linear system solution
    complex(dp), allocatable :: solnvec(:)


    ! --- Miscellaneous --- !

    ! For the constant and z' parts of the summation
    complex(dp), allocatable :: weights(:)

    ! For all parts of the summation except the exponential term
    complex(dp), allocatable :: uhat(:,:)

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
    ny = size(y_pts) - 2
    nt = size(times)

    ! Verify Amat and Imat are the correct size
    if (size(Amat_rowptrb) /= nx*ny .or. size(Amat_rowptre) /= nx*ny) then
        write(*,*) 'invert1: unexpected size for Amat.'
        return
    end if
    if (size(Imat_rowptrb) /= nx*ny .or. size(Imat_rowptre) /= nx*ny) then
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

    ! --- Given matrix processing --- !

    ! Convert Amat, Imat values to complex format
    allocate(Amat_vals_cmplx(size(Amat_vals)), Imat_vals_cmplx(size(Imat_vals)), stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'invert1: allocate: ', msg
        return
    end if
    Amat_vals_cmplx = cmplx(Amat_vals, kind=C_DOUBLE_COMPLEX)
    Imat_vals_cmplx = cmplx(Imat_vals, kind=C_DOUBLE_COMPLEX)

    ! Export Amat to MKL sparse format
    stat = mkl_sparse_z_create_csr(Amat, SPARSE_INDEX_BASE_ONE, nx*ny, nx*ny, &
        Amat_rowptrb, Amat_rowptre, Amat_cols, Amat_vals_cmplx)
    if (stat /= 0) then
        write(*,*) 'invert1: mkl_sparse_d_create_csr: exited with code ', stat
        return
    end if
    
    ! Export Imat to MKL sparse format
    stat = mkl_sparse_z_create_csr(Imat, SPARSE_INDEX_BASE_ONE, nx*ny, nx*ny, &
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
    call c_f_pointer(lhs_rowptrb_c, lhs_rowptrb, [nx*ny])
    call c_f_pointer(lhs_rowptre_c, lhs_rowptre, [nx*ny])
    lhs_nnz = lhs_rowptre(nx*ny) - 1
    call c_f_pointer(lhs_cols_c, lhs_cols, [lhs_nnz])

    ! Allocate lhs_rowidx
    allocate(lhs_rowidx(nx*ny+1), stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'invert1: allocate: ', msg
        return
    end if

    ! Compute lhs_rowidx
    lhs_rowidx(1:nx*ny) = lhs_rowptrb
    lhs_rowidx(nx*ny+1) = lhs_rowptre(nx*ny)


    ! --- Solver setup and re-ordering --- !

    ! Initialize solver
    stat = dss_create(dss_handle, MKL_DSS_DEFAULTS)
    if (stat /= 0) then
        write(*,*) 'invert1: dss_create: exited with code ', stat
        return
    end if

    ! Define structure
    stat = dss_define_structure(dss_handle, MKL_DSS_SYMMETRIC_COMPLEX, &
        lhs_rowidx, nx*ny, nx*ny, lhs_cols, lhs_nnz)
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

    call fem_vec_evalpts_1d(x_pts, dx, n_rhspts_total_x, n_rhspts_x, &
        rhspts_x, stat)
    if (stat /= 0) return
    call fem_vec_evalpts_1d(y_pts, dy, n_rhspts_total_y, n_rhspts_y, &
        rhspts_y, stat)
    if (stat /= 0) return

    allocate(rhsevals(n_rhspts_total_x, n_rhspts_total_y), rhs(nx*ny), &
        stat=stat, errmsg=msg)
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

    allocate(uhat(nx*ny, 0:nz), solnvec(nx*ny), stat=stat, errmsg=msg)
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
        stat = mkl_sparse_z_add(SPARSE_OPERATION_NON_TRANSPOSE, Imat, &
            z_pts(j)**(order+1.0_dp), Amat, lhs)
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


        ! --- Generate the rhs vector --- !

        ! Compute the evals needed for the FEM vector
        if (homogeneous) then
            do k = 1, n_rhspts_total_y
                do l = 1, n_rhspts_total_x
                    rhsevals(l,k) = ic_fn(rhspts_x(l), rhspts_y(k))
                end do
            end do 
        else
            do k = 1, n_rhspts_total_y
                do l = 1, n_rhspts_total_x
                    rhsevals(l,k) = ic_fn(rhspts_x(l), rhspts_y(k)) &
                        + fhat_fn(rhspts_x(l), rhspts_y(k), z_pts(j))
                end do
            end do 
        end if

        ! Compute the FEM vector
        call fem_vec_cmplx_2d(n_rhspts_x, n_rhspts_y, rhspts_x, rhspts_y, &
            rhsevals, rhs, stat)
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

    ! Delete lhs MKL object
    stat = mkl_sparse_destroy(lhs)
    if (stat /= 0) then
        write(*,*) 'invert1: mkl_sparse_destroy: exited with code ', stat
        return
    end if

    ! rhs deallocations
    deallocate(n_rhspts_x, n_rhspts_y, rhspts_x, rhspts_y, rhsevals, rhs, &
        stat=stat, errmsg=msg)
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

! ------------------------------------------------------------------------------
! Compute solution
! ------------------------------------------------------------------------------

    ! Allocate solution
    allocate(soln(nx, ny, nt), stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'invert1: allocate: ', msg
        return
    end if

    ! Compute solution
    do j = 1, nt
        call gemv(uhat, exp(times(j)*z_pts), solnvec)
        soln(:,:,j) = reshape(aimag(solnvec), [nx, ny])
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

end procedure invert1_2d

end submodule inversion_invert1_2d
