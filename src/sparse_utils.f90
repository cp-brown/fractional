module sparse_utils
use utils, only: dp
use mkl_spblas
use iso_c_binding
implicit none
private
public csr_identity, csr_kron_sym

contains

! Computes the n x n identity matrix in MKL CSR format.
subroutine csr_identity(n, rowptrb, rowptre, cols, vals, stat)
! ------------------------------------------------------------------------------
! Declarations
! ------------------------------------------------------------------------------

    ! --- Inputs --- !

    ! Size of matrix
    integer, intent(in) :: n


    ! --- Outputs --- !

    ! Identity matrix in MKL CSR format. There is no need to allocate these
    ! arrays before calling the subroutine.
    integer, allocatable, intent(out) :: rowptrb(:)
    integer, allocatable, intent(out) :: rowptre(:)
    integer, allocatable, intent(out) :: cols(:)
    real(dp), allocatable, intent(out) :: vals(:)
    
    ! Error code
    integer, intent(out) :: stat


    ! --- Local --- !

    ! Iterator
    integer :: j

    ! For error messages
    character(100) :: msg

! ------------------------------------------------------------------------------
! Execution
! ------------------------------------------------------------------------------

    ! De-allocate output arrays, if needed.
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

    ! Allocate output arrays
    allocate(vals(n), cols(n), rowptrb(n), rowptre(n), stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'fdm_id_1d: allocate: ', msg
        return
    end if

    ! Populate the arrays
    do j = 1, n
        rowptrb(j) = j
        rowptre(j) = j+1
        cols(j) = j
        vals(j) = 1.0_dp
    end do

    ! Successful exit
    stat = 0
    return

end subroutine csr_identity



! Computes the Kronecker product of two MKL CSR symmetric matrices. Assumes both
! are square and that all diagonal entries are represented explicitly.
subroutine csr_kron_sym(A_rowptrb, A_rowptre, A_cols, A_vals, B_rowptrb, &
    B_rowptre, B_cols, B_vals, C_rowptrb, C_rowptre, C_cols, C_vals, stat)
! ------------------------------------------------------------------------------
! Declarations
! ------------------------------------------------------------------------------

    ! --- Inputs --- !

    ! Matrix A in MKL CSR format
    integer, intent(in) :: A_rowptrb(:), A_rowptre(:), A_cols(:)
    real(dp), intent(in) :: A_vals(:)

    ! Matrix B in MKL CSR format
    integer, intent(in) :: B_rowptrb(:), B_rowptre(:), B_cols(:)
    real(dp), intent(in) :: B_vals(:)


    ! --- Outputs --- !

    ! Matrix C = kron(A,B) in MKL CSR format. There is no need to allocate these
    ! arrays beforehand.
    integer, allocatable, intent(out) :: C_rowptrb(:), C_rowptre(:), C_cols(:)
    real(dp), allocatable, intent(out) :: C_vals(:)

    ! Error code
    integer, intent(out) :: stat


    ! --- Local --- !

    ! Numbers of rows
    integer :: na, nb, nc

    ! Numbers of non-zero elements
    integer :: nnza, nnzb, nnzc

    ! To hold indexes of the current rows of A and B
    integer :: astart, astop, bstart, bstop, btstart, btstop

    ! To hold current column index of A and B
    integer :: acol, bcol
    
    ! To hold current entries of A and B
    real(dp) :: aval, bval

    ! To keep track of our location within C_cols and C_vals
    integer :: loc

    ! For the transpose of B
    integer, allocatable :: Z_rowptrb(:), Z_rowptre(:), Z_cols(:)
    real(dp), allocatable :: Z_vals(:)
    type(SPARSE_MATRIX_T) :: Bmat, Zmat, Bt
    integer(C_INT) :: Bt_idxing
    integer :: Bt_nrows, Bt_ncols, Bt_nnz
    type(c_ptr) :: Bt_rowptrb_c, Bt_rowptre_c, Bt_cols_c, Bt_vals_c
    integer, pointer :: Bt_rowptrb(:), Bt_rowptre(:), Bt_cols(:)
    real(dp), pointer :: Bt_vals(:)

    ! Iterators
    integer :: i, j, k, l

    ! For error messages
    character(100) :: msg

! ------------------------------------------------------------------------------
! Setup
! ------------------------------------------------------------------------------

    ! Compute sizes
    na = size(A_rowptrb)
    nb = size(B_rowptrb)
    nc = na * nb

    ! Check A and B rowptrs are the same length
    if (size(A_rowptre) /= na) then
        write(*,*) 'csr_kron: matrix A rowptr''s must be the same length.'
        stat = -1
        return
    end if
    if (size(B_rowptre) /= nb) then
        write(*,*) 'csr_kron: matrix B rowptr''s must be the same length.'
        stat = -1
        return
    end if

    ! Compute number of non-zero entries
    nnza = size(A_cols)
    nnzb = size(B_cols)
    nnzc = na*nnzb + (nnza-na)*(2*nnzb-nb)

    ! Check A and B cols and vals are the same length
    if (size(A_vals) /= nnza) then
        write(*,*) 'csr_kron: matrix A cols and vals must be the same length.'
        stat = -1
        return
    end if
    if (size(B_vals) /= nnzb) then
        write(*,*) 'csr_kron: matrix B cols and vals must be the same length.'
        stat = -1
        return
    end if

    ! De-allocate output arrays, if needed
    if (allocated(C_rowptrb)) deallocate(C_rowptrb, stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'csr_kron: deallocate: ', msg
        return
    end if
    if (allocated(C_rowptre)) deallocate(C_rowptre, stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'csr_kron: deallocate: ', msg
        return
    end if
    if (allocated(C_cols)) deallocate(C_cols, stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'csr_kron: deallocate: ', msg
        return
    end if
    if (allocated(C_vals)) deallocate(C_vals, stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'csr_kron: deallocate: ', msg
        return
    end if

    ! Allocate output arrays
    allocate(C_rowptrb(nc), C_rowptre(nc), C_cols(nnzc), C_vals(nnzc), &
        stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'csr_kron: allocate: ', msg
        return
    end if

! ------------------------------------------------------------------------------
! Get transpose of B
! ------------------------------------------------------------------------------

    ! Export Bmat to MKL sparse format
    stat = mkl_sparse_d_create_csr(Bmat, SPARSE_INDEX_BASE_ONE, nb, nb, &
        B_rowptrb, B_rowptre, B_cols, B_vals)
    if (stat /= 0) then
        write(*,*) 'csr_kron_sym: there was an error'
        return
    end if

    ! Export a zero matrix to MKL sparse format
    call csr_identity(nb, Z_rowptrb, Z_rowptre, Z_cols, Z_vals, stat)
    if (stat /= 0) return
    stat = mkl_sparse_d_create_csr(Zmat, SPARSE_INDEX_BASE_ONE, nb, nb, Z_rowptrb, &
        Z_rowptre, Z_cols, Z_vals)
    if (stat /= 0) then
        write(*,*) 'csr_kron_sym: there was an error'
        return
    end if

    ! Transpose
    stat = mkl_sparse_d_add(SPARSE_OPERATION_TRANSPOSE, Bmat, 1.0_dp, Zmat, Bt)
    if (stat /= 0) then
        write(*,*) 'csr_kron_sym: there was an error'
        return
    end if

    ! Destroy Zmat
    stat = mkl_sparse_destroy(Zmat)
    if (stat /= 0) then
        write(*,*) 'csr_kron_sym: there was an error'
        return
    end if

    ! Export to c_ptr arrays
    stat = mkl_sparse_z_export_csr(Bt, Bt_idxing, Bt_nrows, Bt_ncols, &
        Bt_rowptrb_c, Bt_rowptre_c, Bt_cols_c, Bt_vals_c)
    if (stat /= 0) then
        write(*,*) 'csr_kron_sym: there was an error'
        return
    end if

    ! Convert c pointers to fortran pointers
    call c_f_pointer(Bt_rowptrb_c, Bt_rowptrb, [nb])
    call c_f_pointer(Bt_rowptre_c, Bt_rowptre, [nb])
    Bt_nnz = Bt_rowptre(nb) - 1
    call c_f_pointer(Bt_cols_c, Bt_cols, [Bt_nnz])
    call c_f_pointer(Bt_vals_c, Bt_vals, [Bt_nnz])

! ------------------------------------------------------------------------------
! Computation
! ------------------------------------------------------------------------------

    loc = 1
    ! For each row of blocks
    do i = 1, na

        astart = A_rowptrb(i)
        astop = A_rowptre(i)

        ! For each row within the block
        do j = 1, nb

            C_rowptrb( (i-1)*nb + j ) = loc

            bstart = B_rowptrb(j)
            bstop = B_rowptre(j)

            btstart = Bt_rowptrb(j)
            btstop = Bt_rowptre(j)

            ! Diagonal entry of A -- Can approach directly
            acol = A_cols(astart)
            aval = A_vals(astart)
            do l = 0, bstop - bstart - 1

                bcol = B_cols(bstart+l)
                bval = B_vals(bstart+l)

                C_cols(loc) = (acol-1)*nb + bcol
                C_vals(loc) = aval * bval
                loc = loc + 1
                
            end do

            ! Nondiagonal entries of A -- need to include subdiagonal entries of B
            do k = 1, astop - astart - 1

                acol = A_cols(astart+k)
                aval = A_vals(astart+k)

                ! Subdiagonal entries of B
                do l = 0, btstop - btstart - 2

                    bcol = Bt_cols(btstart+l)
                    bval = Bt_vals(btstart+l)

                    C_cols(loc) = (acol - 1) * nb + bcol
                    C_vals(loc) = aval * bval
                    loc = loc + 1

                end do

                ! Diagonal and superdiagonal entries
                do l = 0, bstop - bstart - 1

                    bcol = B_cols(bstart+l)
                    bval = B_vals(bstart+l)

                    C_cols(loc) = (acol - 1) * nb + bcol
                    C_vals(loc) = aval * bval
                    loc = loc + 1

                end do

            end do

            C_rowptre( (i-1)*nb + j ) = loc

        end do

    end do

    ! Destroy Bt
    stat = mkl_sparse_destroy(Bt)
    if (stat /= 0) then
        write(*,*) 'csr_kron_sym: there was an error'
        return
    end if

    ! Successful exit
    stat = 0
    return

end subroutine csr_kron_sym

end module sparse_utils