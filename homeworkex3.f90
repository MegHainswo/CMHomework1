include 'toolbox.f90'
program Golden_method

    use toolbox

    ! parameters :)
    integer, parameter :: n = 5            ! indicates # of subintervals
    integer, parameter :: n_plot = 200      ! indicates # of plot points
    real*8, parameter :: x_u = 5d0          ! creates the upper bound
    real*8, parameter :: x_l = 0d0          ! creates the lower bound
    real*8, parameter :: tol = 1d-6         ! sets the tolerance (as accuracy requirement)
    
    ! variables
    real*8 :: x(0:n), xplot(0:n_plot)
    real*8 :: minimum_x(n), fmin(n), yplot(0:n_plot)
    real*8 :: min_global, fmin_global
    integer :: i, i_global
    
    ! generate a uniform grid
    call grid_Cons_Equi(x, x_l, x_u)
    
    ! find minimums at each interval
    do i=1, n
        minimum_x(i) = minimize(x(i-1), x(i))
        fmin(i) = f(minimum_x(i))
    end do
    
    ! find global minimum
    i_global = minloc(fmin,1)
    min_global = minimum_x(i_global)
    fmin_global = fmin(i_global)
    
    ! display result
    write(*,'(a)') 'The global min. is found at:'
    write(*,'(2(a,f10.6))') 'x=', min_global, 'f(x)=', fmin_global
    
    ! create a visual for the results
    call grid_Cons_Equi(xplot, x_l, x_u)
    do i=0, n_plot
        yplot(i)=f(xplot(i))
    enddo
    
    call plot(xplot,yplot)
    call execplot(xlabel='x', ylabel='f(x)=(x-2.3)^2+sin(5x)')

contains

    ! establishing the function given to us in the problem
    function f(x)
        real*8, intent(in) :: x
        real*8 :: f
        f = x*cos(x**2)
    endfunction 
    
    ! implementing the golden section search algorithm
    function minimize(a,b)
        real*8, intent(in) :: a,b
        real*8 :: minimize
        real*8 :: a1, b1, x1, x2, f1, f2
        integer :: iter
        
        a1 = a
        b1 = b
        
        do iter = 1, 200
            x1 = a1+(3d0-sqrt(5d0))/2d0*(b1-a1)
            x2 = a1+ (sqrt(5d0)-1d0)/2d0*(b1-a1)
            f1 = f(x1)
            f2 = f(x2)
            
            if (f1<f2) then
                b1=x2
            else
                a1=x1
            endif
            
            if (abs(b1-a1)<tol) exit
        enddo
        
        if (f1 < f2) then
            minimize = x1
        else
            minimize = x2
        endif
    endfunction

endprogram

