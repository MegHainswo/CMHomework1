include 'toolbox.f90'
program Matrices

    use toolbox

    ! define variables
    implicit none
    real*8 :: A(4, 4), L(4, 4), U(4, 4), LxU(4,4), b(4)
    integer :: i, j

    ! define matrix A
    A(1, :) = (/1d0, 5d0, 2d0, 3d0/)
    A(2, :) = (/1d0, 6d0, 8d0, 6d0/)
    A(3, :) = (/1d0, 6d0, 11d0, 2d0/)
    A(4, :) = (/1d0, 7d0, 17d0, 4d0/)
    b       = (/1d0, 2d0, 1d0, 1d0/)
    
    ! print A
    write(*,'(a)')'A = '
    write(*,'(4f7.1)')((A(j, i),i = 1, 4), j = 1, 4)

    ! decompose A to L and U, solve for b
    call lu_dec(A, L, U)
    call lu_solve(A,b)

    ! check results
    LxU = matmul(L, U)

    ! print results
    write(*,'(a)')'L = '
    write(*,'(4f7.1)')((L(j, i), i = 1,4), j = 1, 4)
    write(*,'(/a)')'U = '
    write(*,'(4f7.1)')((U(j, i), i = 1, 4), j = 1, 4)
    write(*,'(/a)')'LxU = '
    write(*,'(4f7.1)')((LxU(j, i), i = 1, 4), j = 1, 4)
    write(*,'(30x,a)')'= A'
    write(*,'(a,4f7.1/)')' x=', (b(j),j=1,4)


end program
