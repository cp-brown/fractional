submodule (fem) fem_mat_2d
implicit none

! Submodule containing a subroutine generating the FEM matrices in 2d.

contains

module procedure fem_mats_2d
! ------------------------------------------------------------------------------
! Declarations
! ------------------------------------------------------------------------------

    ! 1d matrices
    integer, allocatable :: xmass_rowptrb(:)
    integer, allocatable :: xmass_rowptre(:)
    integer, allocatable :: xmass_cols(:)
    real(dp), allocatable :: xmass_vals(:)

    integer, allocatable :: ymass_rowptrb(:)
    integer, allocatable :: ymass_rowptre(:)
    integer, allocatable :: ymass_cols(:)
    real(dp), allocatable :: ymass_vals(:)

    integer, allocatable :: xstiff_rowptrb(:)
    integer, allocatable :: xstiff_rowptre(:)
    integer, allocatable :: xstiff_cols(:)
    real(dp), allocatable :: xstiff_vals(:)

    integer, allocatable :: ystiff_rowptrb(:)
    integer, allocatable :: ystiff_rowptre(:)
    integer, allocatable :: ystiff_cols(:)
    real(dp), allocatable :: ystiff_vals(:)

    ! Kronecker products
    integer, allocatable :: ysxm_rowptrb(:)
    integer, allocatable :: ysxm_rowptre(:)
    integer, allocatable :: ysxm_cols(:)
    real(dp), allocatable :: ysxm_vals(:)

    ! For error messages
    character(100) :: msg

! ------------------------------------------------------------------------------
! Computation
! ------------------------------------------------------------------------------

    ! Get 1d mass matrices
    call fem_mass_1d(xpts, xmass_rowptrb, xmass_rowptre, xmass_cols, &
        xmass_vals, stat)
    if (stat /= 0) return
    call fem_mass_1d(ypts, ymass_rowptrb, ymass_rowptre, ymass_cols, &
        ymass_vals, stat)
    if (stat /= 0) return

    ! Get 1d stiffness matrices
    call fem_stiff_1d(xpts, diffusivity_x_fn, dx, xstiff_rowptrb, &
        xstiff_rowptre, xstiff_cols, xstiff_vals, stat)
    if (stat /= 0) return
    call fem_stiff_1d(ypts, diffusivity_y_fn, dy, ystiff_rowptrb, &
        ystiff_rowptre, ystiff_cols, ystiff_vals, stat)
    if (stat /= 0) return

    ! Compute (2d) mass matrix as kron(ymass, xmass)
    call csr_kron_sym(ymass_rowptrb, ymass_rowptre, ymass_cols, ymass_vals, &
        xmass_rowptrb, xmass_rowptre, xmass_cols, xmass_vals, mass_rowptrb, &
        mass_rowptre, mass_cols, mass_vals, stat)
    if (stat /= 0) return


    ! --- Compute tensor products needed for stiffness matrix --- !

    ! xmys = kron(ymass, xstiff)
    call csr_kron_sym(ymass_rowptrb, ymass_rowptre, ymass_cols, ymass_vals, &
        xstiff_rowptrb, xstiff_rowptre, xstiff_cols, xstiff_vals, &
        stiff_rowptrb, stiff_rowptre, stiff_cols, stiff_vals, stat)
    if (stat /= 0) return

    ! No longer need ymass, xstiff
    deallocate(ymass_rowptrb, ymass_rowptre, ymass_cols, ymass_vals, &
        xstiff_rowptrb, xstiff_rowptre, xstiff_cols, xstiff_vals, &
        stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'fem_mats_2d: deallocate: ', msg
        return
    end if

    ! xsym = kron(ystiff, xmass)
    call csr_kron_sym(ystiff_rowptrb, ystiff_rowptre, ystiff_cols, ystiff_vals, &
        xmass_rowptrb, xmass_rowptre, xmass_cols, xmass_vals, &
        ysxm_rowptrb, ysxm_rowptre, ysxm_cols, ysxm_vals, stat)
    if (stat /= 0) return

    ! No longer need ystiff, xmass
    deallocate(xmass_rowptrb, xmass_rowptre, xmass_cols, xmass_vals, &
        ystiff_rowptrb, ystiff_rowptre, ystiff_cols, ystiff_vals, &
        stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'fem_mats_2d: deallocate: ', msg
        return
    end if


    ! --- Add xmys, xsym to get stiffness matrix --- !

    ! It would probably be good practice to use the MKL sparse BLAS routines
    ! to compute the sum, however since all the 1d matrices are tridiagonal,
    ! xmys and xsym have the same structure, so we can simply add the values.
    stiff_vals = stiff_vals + ysxm_vals

    ! No longer need ysxm
    deallocate(ysxm_rowptrb, ysxm_rowptre, ysxm_cols, ysxm_vals, &
        stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'fem_mats_2d: deallocate: ', msg
        return
    end if

end procedure fem_mats_2d

end submodule fem_mat_2d