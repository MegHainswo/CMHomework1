include 'toolbox.f90'
program Oligolopy
    use toolbox
    implicit none

    ! parameters :)
    integer, parameter :: N = 10, NP = 1000
    real*8, parameter :: p_l = 0.1d0, p_u = 3.0d0
    real*8, parameter :: alpha = 1, eta = 1.5d0
    integer, parameter :: m = 3

    ! variables and arrays
    real*8 :: P(0:N), q(0:N), coeff_q(N+3)
    real*8 :: p_plot(0:NP), q_s_plot(0:NP), q_d_plot(0:NP)
    real*8 :: q_in, p_star, q_star
    integer :: ip, ip_com
    logical :: check

    ! task a) compute individual firm output across price grid
    call grid_Cons_Equi(P, p_l, p_u)
    coeff_q = 0d0

    do ip = 0, N
        q_in = 1d0
        ip_com = ip
        call fzero(q_in, foc, check)
        if (check) stop 'Error: fzero did not converge'
        q(ip) = q_in
    end do

    ! print results task a)
    write(*,'(/a)') 'Computed individual quantities at each price grid point:'
    write(*,'(a)') '  Price      Quantity'
    write(*,'(a)') '---------------------------'
    do ip = 0, N
        write(*,'(f7.4,2x,f10.6)') P(ip), q(ip)
    end do

    ! task b) interpolate and plot supply and demand curves
    call spline_interp(q, coeff_q)
    call grid_Cons_Equi(p_plot, p_l, p_u)

    do ip = 0, NP
        q_s_plot(ip) = m * spline_eval(p_plot(ip), coeff_q, p_l, p_u)
        q_d_plot(ip) = p_plot(ip)**(-eta)
    end do

    ! print results task b)
    write(*,'(/a)') 'Sampled market supply and demand curves:'
    write(*,'(a)') '  Price    Demand     Supply'
    write(*,'(a)') '-------------------------------'
    do ip = 0, NP, 100
        write(*,'(f7.4,2x,f10.6,2x,f10.6)') p_plot(ip), q_d_plot(ip), q_s_plot(ip)
    end do

    ! plot supply and demand curves
    call plot(p_plot, q_d_plot, legend='Market Demand')
    call plot(p_plot, q_s_plot, legend='Aggregate Supply')
    call execplot(xlabel='Price', ylabel='Quantity', &
                  xlim=(/p_l, p_u/), ylim=(/0d0, 2d0/))

    ! task c) find equilibrium p such that: D(p)=m*q(p)
    p_star = 1d0  ! Initial guess
    call fzero(p_star, market_eq, check)
    if (check) stop 'Error: fzero did not converge for equilibrium'

    q_star = p_star**(-eta)

    ! print results task c)
    write(*,'(/a,f10.6)') 'Equilibrium Price:    ', p_star
    write(*,'(a,f10.6)')  'Equilibrium Quantity: ', q_star

contains

    ! foc
    function foc(q_in)
        implicit none
        real*8, intent(in) :: q_in
        real*8 :: foc, p_now, dpdq
        integer :: ip_com_local
        save ip_com_local
        ip_com_local = ip_com

        p_now = P(ip_com_local)
        dpdq = -p_now**(1d0 + eta) / eta
        foc = p_now + q_in * dpdq - (alpha * sqrt(q_in) + q_in**2)
    end function

    ! this function finds residual in market-clearing condition:
    function market_eq(p_val)
        implicit none
        real*8, intent(in) :: p_val
        real*8 :: market_eq, q_interp

        q_interp = spline_eval(p_val, coeff_q, p_l, p_u)
        market_eq = p_val**(-eta) - m * q_interp
    end function

end program
    


