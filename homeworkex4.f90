include 'toolbox.f90'
program Household_2
    use toolbox
    implicit none

    ! parameters :)
    real*8, parameter :: gamma = 0.5d0
    real*8, parameter :: ex = 1d0 - 1d0 / gamma
    real*8, parameter :: r = 0.0d0
    real*8, parameter :: beta = 1d0
    real*8, parameter :: w = 1d0

    ! variables
    real*8 :: x_in(2), a(2), b(2), fret, c1

    ! initialization
    a = (/ 0.01d0, 0.01d0 /)                  ! Lower bounds for c2, c3
    b = (/ 2d0, 2d0 /)                        ! Upper bounds for c2, c3
    x_in = (/ w / 2d0, w / 2d0 /)             ! Initial guess: equal split

    ! optimization using fminsearch 
    call fminsearch(x_in, fret, a, b, utility)

    ! recover c1 from the budget constraint
    c1 = w + w / (1d0 + r) - x_in(1) / (1d0 + r) - x_in(2) / (1d0 + r)**2

    ! output 
    write(*, '(/a)') 'Result with fminsearch:'
    write(*, '(a,f10.6)') 'c_1 = ', c1
    write(*, '(a,f10.6)') 'c_2 = ', x_in(1)
    write(*, '(a,f10.6)') 'c_3 = ', x_in(2)

contains

    function utility(x) result(uval)
        implicit none
        real*8, intent(in) :: x(:)
        real*8 :: uval, c1, c2, c3

        ! local variables use outer scope values (module constants)

        c2 = x(1)
        c3 = x(2)

        ! derive c1 using the budget constraint
        c1 = w + w / (1d0 + r) - c2 / (1d0 + r) - c3 / (1d0 + r)**2

        ! penalize infeasible (negative) consumption
        if (c1 <= 0d0 .or. c2 <= 0d0 .or. c3 <= 0d0) then
            uval = 1d6
            return
        end if

        ! compute utility (negated for minimization)
        uval = - (c1**ex + beta * c2**ex + beta**2d0 * c3**ex) / ex
    end function

end program
    



