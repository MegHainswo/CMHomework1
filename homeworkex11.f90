include 'toolbox.f90'
program Cost_Minimization

    use toolbox 
    implicit none
    
    ! parameters :)
    integer, parameter :: n=3   ! # gravel pits
    integer, parameter :: m=4   ! # building sites
    
    ! variables
    real*8 :: c(n*m), x(n*m)
    real*8 :: A(n+m, n*m), b(n+m)
    
    ! cost table
    c(:) = (/ 10d0, 790d0, 10d0, 80d0, &
            130d0, 90d0, 120d0, 11d0, &
            50d0, 30d0, 80d0, 10d0 /)

    
    ! supply constraints
    A(1,1:4) = 1d0
    A(2,5:8) = 1d0
    A(3,9:12) = 1d0
    
    ! demand constraints
    A(4,(/1,5,9/)) = 1d0
    A(5,(/2,6,10/)) = 1d0
    A(6,(/3,7,11/)) = 1d0
    A(7,(/4,8,12/)) = 1d0
    
    ! right-hand side vector
    b(:) = (/ 11d0, 13d0, 10d0, 5d0, 7d0, 13d0, 6d0 /)
    
    call solve_lin(x,c,A,b,3,0,4)
    
    ! results
    write(*,'(a)') 'Building site      :     B1     B2     B3     B4    Total  Total Cost'
    write(*,'(a,4f8.1,2f10.1)') 'Gravel-pit A1     :', x(1:4), sum(x(1:4)), sum(x(1:4)*c(1:4))
    write(*,'(a,4f8.1,2f10.1)') 'Gravel-pit A2     :', x(5:8), sum(x(5:8)), sum(x(5:8)*c(5:8))
    write(*,'(a,4f8.1,2f10.1)') 'Gravel-pit A3     :', x(9:12), sum(x(9:12)), sum(x(9:12)*c(9:12))
    write(*,'(a,4f8.1,2f10.1)') 'Total delivered   :', &
     x(1)+x(5)+x(9), x(2)+x(6)+x(10), x(3)+x(7)+x(11), x(4)+x(8)+x(12), sum(x), sum(x*c)


end program
    



