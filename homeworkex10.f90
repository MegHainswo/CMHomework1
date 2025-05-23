include 'toolbox.f90'
program NewspaperMonop

  use toolbox
  implicit none

  ! parameters :)
  integer, parameter :: NP       = 3           ! number of grid points
  real*8,   parameter :: p_min   = 0.5d0       ! lower price bound
  real*8,   parameter :: p_max   = 12.5d0      ! upper price bound
  real*8,   parameter :: c       = 0.1d0       ! marginal cost
  integer, parameter :: Mplots   = 3           
  integer, parameter :: Nplot_list(Mplots) = (/100, 1000, 10000/)

  ! coarse grid and sample profits
  real*8 :: pr(0:NP), pa(0:NP)
  real*8 :: G_sample(0:NP,0:NP)
  real*8 :: coeffG(NP+3,NP+3)

  real*8, allocatable :: pr_plot(:), pa_plot(:), G_plot(:,:)

  integer :: ip, ir, ia, ix(2), Nplot
  real*8 :: pr_star_a, pa_star_a, profit_star_a

  ! for part (b)
  real*8 :: x0(2), fret
  real*8 :: pr_star_b, pa_star_b, profit_star_b

  ! task (a): 2D spline interpolation to find otimum
  call grid_Cons_Equi(pr, p_min, p_max)
  call grid_Cons_Equi(pa, p_min, p_max)

  ! manually specified profit values at coarse grid points
  G_sample = 0d0
  G_sample(0,:) = (/ 11.5d0,  70.9d0,  98.3d0,  93.7d0 /)
  G_sample(1,:) = (/ 31.1d0,  82.5d0, 101.9d0,  89.3d0 /)
  G_sample(2,:) = (/ 18.7d0,  62.1d0,  73.5d0,  52.9d0 /)
  G_sample(3,:) = (/-25.7d0,   9.7d0,  13.1d0, -15.5d0 /)

  coeffG = 0d0
  call spline_interp(G_sample, coeffG)

  write(*,*) 'Part (a): spline-based optimum'
  do ip = 1, Mplots
    Nplot = Nplot_list(ip)
    allocate(pr_plot(0:Nplot), pa_plot(0:Nplot), G_plot(0:Nplot,0:Nplot))

    call grid_Cons_Equi(pr_plot, p_min, p_max)
    call grid_Cons_Equi(pa_plot, p_min, p_max)

    ! evaluate interpolated profit surface
    do ir = 0, Nplot
      do ia = 0, Nplot
        G_plot(ir,ia) = spline_eval( (/ pr_plot(ir), pa_plot(ia) /), &
                                     coeffG, (/p_min,p_min/), (/p_max,p_max/) )
      end do
    end do

    ! find grid location of maximum profit
    ix = maxloc(G_plot) - 1
    pr_star_a     = pr_plot(ix(1))
    pa_star_a     = pa_plot(ix(2))
    profit_star_a = G_plot(ix(1),ix(2))

    write(*,'(a,i6,3(2x,a,f10.6))') ' Nplot=', Nplot,       &
         ' pR*=', pr_star_a, ' pA*=', pa_star_a, ' G*=', profit_star_a

    deallocate(pr_plot, pa_plot, G_plot)
  end do

  ! task (b): true optimum via fminsearch
  x0 = (/1.0d0, 1.0d0/)
  call fminsearch(x0, fret, (/p_min,p_min/), (/p_max,p_max/), profit_fn)
  pr_star_b     = x0(1)
  pa_star_b     = x0(2)
  profit_star_b = -fret
  write(*,'(a,2x,a,f10.6,2x,a,f10.6,2x,a,f10.6)')                  &
       'Part (b) fminsearch optimum:', ' pR*=', pr_star_b,         &
       ' pA*=', pa_star_b, ' G*=', profit_star_b

  ! task (c): interpolation error at true optimum
  write(*,*)
  write(*,*) 'Part (c): |G_true - G_interp|'
  do ip = 1, Mplots
    Nplot = Nplot_list(ip)
    profit_star_a = spline_eval( (/ pr_star_b, pa_star_b /), &
                                 coeffG, (/p_min,p_min/), (/p_max,p_max/) )
    write(*,'(a,i6,2x,a,f10.6,2x,a,f10.6,2x,a,f10.6)')            &
      ' Nplot=', Nplot, ' G_true=', profit_star_b,               &
      ' G_interp=', profit_star_a,                                &
      ' |err|=', abs(profit_star_b - profit_star_a)
  end do

contains

  function profit_fn(x)
    implicit none
    real*8, intent(in) :: x(:)
    real*8 :: profit_fn, pR, pA, qR, qA
    pR = x(1);  pA = x(2)
    qR = 10d0 - pR
    qA = 20d0 - pA - 0.5d0*pR
    profit_fn = - (pR*qR + pA*qA - c*(qR + qA))
  end function

  function foc_fn(x)
    implicit none
    real*8, intent(in) :: x(:)
    real*8 :: foc_fn(size(x)), pR, pA
    pR = x(1);  pA = x(2)
    foc_fn(1) = 10d0 - 2d0*pR - 0.5d0*pA + 1.5d0*c
    foc_fn(2) = 20d0 - 2d0*pA - 0.5d0*pR +    1.0d0*c
  end function

end program

            
        
    
    


