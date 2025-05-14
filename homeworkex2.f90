include "toolbox.f90"

program Household

    use toolbox
    implicit none
    
    ! parameters :)
    real*8, parameter :: beta = 1.0d0
    real*8, parameter :: r = 0.0d0
    real*8, parameter :: w = 1.0d0
    real*8, parameter :: gamma = 0.5d0
    real*8, parameter :: ex = 1d0-1d0/gamma 

    ! declaration of variables
    real*8 :: x_root(3), x_minimize, a, b, fret
    logical :: check

    ! initialize x_root
    x_root = -50

    ! call subroutine fzero, check
    call fzero(x_root, foc, check)
    if(check)stop 'Looks like something is wrong'

    ! display fzero result 
    write(*,'(a)') 'Result with fzero:'
    write(*,'(a,f10.6)') 'c_1 = ', x_root(1)
    write(*,'(a,f10.6)') 'c_2 = ', x_root(2)
    write(*,'(a,f10.6)') 'lambda = ', x_root(3)

    ! initialize interval for fminsearch and x_minimize
    a = 0d0
    b = w
    x_minimize = w/2d0

    ! call subroutine fminsearch
    call fminsearch(x_minimize, fret, a, b, utility)

    ! display fminsearch result
    write(*,'(/a)') 'Result with fminsearch:'
    write(*,'(a,f10.6)') 'c_1 = ', x_minimize
    write(*,'(a,f10.6)') 'c_2 = ', (w - x_minimize)*(1d0+r)

contains
    
    function foc(x_in) result(fval)
        implicit none
        real*8, intent(in) :: x_in(:)
        real*8 :: fval(size(x_in))
        
        fval(1)=x_in(1)**(-1d0/gamma)-x_in(3)
        fval(2)=beta*x_in(2)**(-1d0/gamma)-x_in(3)/(1+r)
        fval(3)=w-x_in(1)-x_in(2)/(1+r)
    end function foc
    
    function utility(x_in) result(uval)
        implicit none
        real*8, intent(in) :: x_in
        real*8 :: uval
        
        uval = -((x_in**ex)/ex+beta *((w -x_in)**ex)/ex)
    end function utility

end program

