submodule (fem) fem_vec_1d
implicit none

! Submodule containing subroutines generating the FEM load vector in 1d.

contains

module procedure fem_vec_evalpts_1d

    ! --- Declarations --- !

    ! Number of sub-intervals
    integer :: n

    ! To hold sizes of intervals
    real(dp) :: h

    ! To keep track of where we are in the evalpts array
    integer :: l

    ! Iterators
    integer :: j, k
    
    ! For error messages
    character(100) :: msg


    ! --- Begin program --- !

    ! De-allocate output arrays, if needed
    if (allocated(n_evals)) deallocate(n_evals, stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'fem_vec_evalpts_1d: deallocate: ', msg
        return
    end if
    if (allocated(evalpts)) deallocate(evalpts, stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'fem_vec_evalpts_1d: deallocate: ', msg
        return
    end if

    ! Number of sub-intervals
    n = size(gridpts) - 1

    ! Allocate n_evals
    allocate(n_evals(n), stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'fem_vec_evalpts_1d: allocate: ', msg
        return
    end if

    ! Compute n_evals
    do j = 1, n
        n_evals(j) = ceiling((gridpts(j+1) - gridpts(j)) / dx) - 1
    end do

    ! Total number of evals needed
    n_evals_total = sum(n_evals) + n + 1

    ! Allocate evalpts to the size just computed
    allocate(evalpts(n_evals_total), stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'fem_vec_evalpts_1d: allocate: ', msg
        return
    end if

    ! Compute the evalpts
    l = 1
    do j = 1, n

        ! Interval width
        h = (gridpts(j+1) - gridpts(j)) / (n_evals(j) + 1.0_dp)

        ! Eval point at left of interval
        evalpts(l) = gridpts(j)
        l = l + 1

        ! Eval points inside the interval
        do k = 1, n_evals(j)
            evalpts(l) = evalpts(l-1) + h
            l = l + 1
        end do

        ! Eval point at right of interval will be added later
        
    end do

    ! Don't forget the last grid point!
    evalpts(l) = gridpts(n+1)

    ! Successful exit
    stat = 0
    return

end procedure fem_vec_evalpts_1d



module procedure fem_vec_cmplx_1d

    ! --- Declarations of local variables --- !

    ! Number of interior grid points
    integer :: n

    ! To hold the indexes of grid points
    integer :: xm1idx, xidx, xp1idx

    ! To hold the grid points
    real(dp) :: xm1, x, xp1

    ! Iterator
    integer :: j


    ! --- Begin program --- !

    ! Number of interior points
    n = size(n_evals) - 1

    ! Inputs checking
    if (size(evalpts) /= sum(n_evals) + n + 2) then
        write(*,*) 'fem_vec_cmplx_1d: Unexpected number of evalpts'
        stat = -1
        return
    end if
    if (size(evals) /= size(evalpts)) then
        write(*,*) 'fem_vec_cmplx_1d: Unexpected number of evals'
        stat = -1
        return
    end if
    if (size(femvec) /= n) then
        write(*,*) 'fem_vec_cmplx_1d: Unexpected size for femvec'
        stat = -1
        return
    end if

    ! First and second grid point indexes
    xm1idx = 1
    xidx = 2 + n_evals(1)

    ! Compute entries
    do j = 1, n

        ! Index of x(j+1)
        xp1idx = xidx + n_evals(j+1) + 1

        ! Relevant grid points, for brevity
        xm1 = evalpts(xm1idx)
        x = evalpts(xidx)
        xp1 = evalpts(xp1idx)

        ! Entries are:
        ! 1/h(i) * int_x(j-1)^x(j) f(x) (x - x(i-1)) dx
        ! - 1/h(i+1) * int_x(j)^x(j+1) f(x) (x - x(i+1)) dx
        femvec(j) = defint_trap_cmplx( evals(xm1idx:xidx) &
            * (evalpts(xm1idx:xidx) - xm1), xm1, x) / (x - xm1) &
            + defint_trap_cmplx( evals(xidx:xp1idx) * (evalpts(xidx:xp1idx) &
            - xp1), x, xp1) / (x - xp1)

        xm1idx = xidx
        xidx = xp1idx

    end do

    ! Successful exit
    stat = 0
    return

end procedure fem_vec_cmplx_1d

end submodule fem_vec_1d