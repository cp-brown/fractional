submodule (fem) fem_vec_2d
implicit none

! Submodule containing subroutines generating the FEM load vector in 2d.

contains

module procedure fem_vec_cmplx_2d

    ! --- Declarations of local variables --- !

    ! Sizes of things
    integer :: nx, ny, nxtot, nytot

    ! For the temporary result
    complex(dp), allocatable :: xreduc(:,:)

    ! Iterator
    integer :: j

    ! For error messages
    character(100) :: msg


    ! --- Begin program --- !

    ! Number of interior points
    nx = size(n_evals_x) - 1
    ny = size(n_evals_y) - 1

    ! Total number of points
    nxtot = size(evalpts_x)
    nytot = size(evalpts_y)

    ! Allocate xreduc to size Nx x nytot
    allocate(xreduc(nx, nytot), stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'fem_vec_cmplx_2d: allocate: ', msg
        return
    end if

    ! Integrate in the x-direction
    do j = 1, nytot
        call fem_vec_cmplx_1d(n_evals_x, evalpts_x, evals(:,j), xreduc(:,j), stat)
    end do

    ! Integrate in the y-direction
    do j = 1, nx
        call fem_vec_cmplx_1d(n_evals_y, evalpts_y, xreduc(j,:), &
            femvec( j : (ny-1)*nx + j : nx ), stat)
    end do

    ! Deallocate xreduc
    deallocate(xreduc, stat=stat, errmsg=msg)
    if (stat /= 0) then
        write(*,*) 'fem_vec_cmplx_2d: deallocate: ', msg
        return
    end if

    ! Successful exit
    stat = 0
    return

end procedure fem_vec_cmplx_2d

end submodule fem_vec_2d