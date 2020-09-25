submodule (fem) fem_mat_1d
implicit none

! Submodule containing subroutines generating the FEM matrices in 1d.

contains

module procedure fem_mass_1d
! ------------------------------------------------------------------------------
! Declarations of local variables
! ------------------------------------------------------------------------------

    ! Size of matrix
    integer :: n

    ! Number of non-zero entries
    integer :: nnz

    ! For keeping track of location within arrays
    integer :: k

    ! Iterator
    integer :: j

    ! For error messages
    character(100) :: msg

! ------------------------------------------------------------------------------
! Setup
! ------------------------------------------------------------------------------

    ! Size of matrix
    n = size(pts) - 2

    ! Number of nonzero elements: full diagonal and first superdiagonal
    nnz = n + (n - 1)

    ! De-allocate the CSR representation arrays, if needed
    if (allocated(rowptrb)) deallocate(rowptrb, stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'fem_mass_1d: deallocate: ', msg
        return
    end if
    if (allocated(rowptre)) deallocate(rowptre, stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'fem_mass_1d: deallocate: ', msg
        return
    end if
    if (allocated(cols)) deallocate(cols, stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'fem_mass_1d: deallocate: ', msg
        return
    end if
    if (allocated(vals)) deallocate(vals, stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'fem_mass_1d: deallocate: ', msg
        return
    end if
    
    ! Allocate the CSR representation arrays
    allocate(vals(nnz), cols(nnz), rowptrb(n), rowptre(n), stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'fem_mass_1d: allocate: ', msg
        return
    end if

! ------------------------------------------------------------------------------
! Computation
! ------------------------------------------------------------------------------

    ! Populate the arrays - first n-1 rows
    k = 1
    do j = 1, n-1

        rowptrb(j) = k
        rowptre(j) = k+2

        ! Diagonal value
        vals(k) = 1.0_dp/3.0_dp * (pts(j+2) - pts(j))
        cols(k) = j

        ! Superdiagonal value
        vals(k+1) = 1.0_dp/6.0_dp * (pts(j+2) - pts(j+1))
        cols(k+1) = j+1

        k = k + 2

    end do

    ! Last row
    rowptrb(n) = k
    rowptre(n) = k+1
    vals(k) = 1.0_dp/3.0_dp * (pts(n+2) - pts(n))
    cols(k) = n

    ! Successful exit
    stat = 0
    return

end procedure fem_mass_1d



module procedure fem_stiff_1d
! ------------------------------------------------------------------------------
! Declarations of local variables
! ------------------------------------------------------------------------------

    ! Size of matrix
    integer :: n

    ! Number of non-zero entries
    integer :: nnz

    ! For holding integrals, to avoid re-computation
    real(dp) :: intprev, intcur

    ! To keep track of location within the arrays
    integer :: k

    ! Iterator
    integer :: j

    ! For error messages
    character(100) :: msg

! ------------------------------------------------------------------------------
! Setup
! ------------------------------------------------------------------------------

    ! Size of matrix
    n = size(pts) - 2

    ! Number of nonzero elements: full diagonal and first superdiagonal
    nnz = n + (n - 1)

    ! De-allocate the CSR representation arrays, if needed
    if (allocated(rowptrb)) deallocate(rowptrb, stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'csr_identity: deallocate: ', msg
        return
    end if
    if (allocated(rowptre)) deallocate(rowptre, stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'csr_identity: deallocate: ', msg
        return
    end if
    if (allocated(cols)) deallocate(cols, stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'csr_identity: deallocate: ', msg
        return
    end if
    if (allocated(vals)) deallocate(vals, stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'csr_identity: deallocate: ', msg
        return
    end if
    
    ! Allocate the CSR representation arrays
    allocate(vals(nnz), cols(nnz), rowptrb(n), rowptre(n), stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'fdm_id_1d: allocate: ', msg
        return
    end if

! ------------------------------------------------------------------------------
! Computation
! ------------------------------------------------------------------------------

    ! We need the integrals  I(i) = hi^-2 * int_xi-1^xi sigma(x) dx:
    ! Diagonal entries are  -(I(i) + I(i+1)), where i = row index
    ! superdiagonal entries are  I(i+1).

    ! I(1)
    intprev = defint_fn_trap_real(diffusivity_fn, pts(1), pts(2), dx) &
        / (pts(2) - pts(1))**2.0_dp
    
    ! First n-1 rows
    k = 1
    do j = 1, n-1

        rowptrb(j) = k
        rowptre(j) = k+2

        ! Compute I(i+1)
        intcur = defint_fn_trap_real(diffusivity_fn, pts(j), pts(j+1), dx) &
            / (pts(j+1) - pts(j))**2.0_dp

        ! Diagonal value
        vals(k) = -1.0_dp * (intprev + intcur)
        cols(k) = j

        ! Superdiagonal value
        vals(k+1) = intcur
        cols(k+1) = j+1

        ! So that intprev holds I(i) in the next iteration
        intprev = intcur

        k = k + 2

    end do
    
    ! Last row
    rowptrb(n) = k
    rowptre(n) = k+1
    intcur = defint_fn_trap_real(diffusivity_fn, pts(n+1), pts(n+2), dx) &
        / (pts(n+2) - pts(n+1))**2.0_dp
    vals(k) = -1.0_dp * (intprev + intcur)
    cols(k) = n

    ! Successful exit
    stat = 0
    return

end procedure fem_stiff_1d

end submodule fem_mat_1d