!=======================================================================
!  Program: bound_states_polar
!
!  Purpose:
!       This program computes polar gravitational bound states inside
!       a Schwarzschild black hole by solving the corresponding
!       perturbation equation for the polar sector.
!
!       The code searches for purely imaginary eigenfrequencies
!       omega = i omegaI that satisfy physical boundary conditions at
!       the black-hole horizon (r → 1) and the origin (r → 0).
!       A shooting method combined with a random-walk minimization
!       is used to determine the correct omegaI.
!
!  Physics and Mathematical Background:
!       • The master equation for polar perturbations is solved for a
!         given multipole number l.
!       • A Frobenius-type power-series expansion is used near r = 0
!         to obtain initial boundary conditions (c_series).
!       • The wavefunction is defined as
!         Z(r) = r g(r)/[e^{i omega r} (r - 1)^{i omega} (1.5 + lambda r)^2].
!       • The solutions Zm(r) (integrated outward) and Zp(r)
!         (integrated inward) are matched via the Wronskian
!         W = Zm Zp' − Zm' Zp.
!       • The eigenvalue omegaI corresponds to a minimum of |W|.
!
!  Numerical Methods:
!       • High-order Runge–Kutta (RK8) integrator for solving the ODEs.
!       • Series expansion up to n_series terms around r = 0.
!       • Uniform grid of n_array points between r ≈ 0 and r ≈ 1.
!       • Trapezoidal-rule normalization of Zm and Zp.
!       • Random-walk search in the interval [omegaI_min, omegaI_max] with
!         dynamically shrinking step size until delta_omegaI < delta_omegaI_tol.
!
!  Parallelization Strategy (MPI):
!       • Each MPI process evaluates the Wronskian for a different trial
!        value of omegaI.
!       • Values of (omegaI, f(omegaI)) are gathered on rank 0.
!       • Rank 0 selects the best candidate (smallest |W|), updates the
!        interval, and broadcasts the new omegaI to all processes.
!       • Iteration continues until convergence is achieved.
!
!  Input Parameters:
!       • l: multipole number (e.g. l = 2, 3, 4, …)
!       • M: black hole mass (here fixed to M = 0.5)
!       • omegaI_start: initial guess for the imaginary frequency
!       • (omegaI_min, omegaI_max): convergence interval
!
!  Output:
!       • The converged eigenfrequency omegaI.
!       • Arrays r, Zm(r), Zp(r), W(r), and |W(r)|.
!
!  Dependencies:
!       • Requires an MPI compiler.
!       • All computations use REAL(16) and COMPLEX(16) quadruple precision.
!
!  USAGE:
!       • Compile with MPI Fortran compiler, e.g.:
!         mpiifort -O3 bound_states_polar.f90 -o bound_states_polar
!       • Run with desired number of parallel chains:
!         mpirun -np 8 ./bound_states_polar
!
!  REPOSITORY:
!       https://github.com/krezazadeh/bound_states_polar
!
!  REFERENCE:
!       https://arxiv.org/abs/2511.15140
!
!  AUTHOR:
!       Kazem Rezazadeh
!       School of Astronomy,
!       Institute for Research in Fundamental Sciences (IPM)
!       Email: kazem.rezazadeh@ipm.ir
!       GitHub: https://github.com/krezazadeh
!
!  DATE:
!       14 November 2025
!=======================================================================

program bound_states_polar

    use mpi

    implicit none

    integer :: ierr, rank, nprocess
    character(len=30) :: filename
    integer :: unit

    real(16), parameter :: M = 0.5q0
    
    integer, parameter :: nr = 100000
    integer :: i, j
    real(16) :: i_real

    real(16) :: l
    real(16) :: omegaI, LL, lambda
    complex(16) :: omega
    real(16) :: rbc1, rbc2, dr
    real(16) :: C1
    real(16) :: rinitial
    complex(16) :: ginitial, gpinitial

    integer, parameter :: n_array = 100000
    real(16) :: array_rm(n_array), array_rp(n_array)
    complex(16) :: array_Zm(n_array)
    complex(16) :: array_Zmp(n_array) 
    complex(16) :: array_Zp(n_array)
    complex(16) :: array_Zpp(n_array)
    complex(16) :: array_W(n_array)
    real(16) :: array_absW(n_array)

    integer, parameter :: n_series = 100
    complex(16) :: c_series(n_series)

    real(16) :: ri
    complex(16) :: gi, gpi
    complex(16) :: Zi

    integer, parameter :: n_try_max = 30
    integer :: n_try

    real(16) :: rand

    real(16) :: omegaI_start
    real(16) :: omegaI_result, omegaI_min, omegaI_max
    real(16) :: delta_omegaI, delta_omegaI_start
    real(16) :: max_absW, max_absW_result
    ! real(16) :: f_omegaI, f_omegaI_result
    real(16) :: f_omegaI1, f_omegaI1_result
    real(16) :: f_omegaI2, f_omegaI2_result

    ! real(16), parameter :: delta_omegaI_tol = 1.0q-10
    real(16) :: omegaI_result_old
    
    ! real(16), allocatable :: omegaI_process(:), fomegaI_process(:)
    real(16), allocatable :: omegaI_process(:), fomegaI1_process(:), fomegaI2_process(:)

    logical :: improved

    ! exact results

    ! l = 2.0q0
    ! omegaI = 4.000000000000000000000000000000000000000q0

    ! l = 3.0q0
    ! omegaI = 1.705649855867412886401121336882458301887q0
    ! omegaI = 20.00000000000000000000000000000000000000q0

    ! l = 4.0q0
    ! omegaI = 1.394461407336213459264897395377440300424q0
    ! omegaI = 5.745972088109590244621221888097850268509q0
    ! omegaI = 60.00000000000000000000000000000000000000q0

    ! l = 5.0q0
    ! omegaI = 1.274578720545848220578308097601014984390q0
    ! omegaI = 3.968693503443007241452978377145027652482q0
    ! omegaI = 13.76767449911914180521941737341018867631q0
    ! omegaI = 140.0000000000000000000000000000000000000q0

    ! l = 6.0q0
    ! omegaI = 1.211013016622633341214320861522351447501q0
    ! omegaI = 3.315604822709691677976475072452344351991q0
    ! omegaI = 8.389356733797011655824954515669729486190q0
    ! omegaI = 27.78537685419365464281140450876736315183q0
    ! omegaI = 280.0000000000000000000000000000000000000q0

    ! l = 7.0q0
    ! omegaI = 1.171550415416408937839788644242667740414q0
    ! omegaI = 2.982338126696478388906388236995354668964q0
    ! omegaI = 6.487859859176734474087274490565748183308q0
    ! omegaI = 15.42616838129598421337714051813568435383q0
    ! omegaI = 50.20490177786582672678667336176488879974q0
    ! omegaI = 504.0000000000000000000000000000000000000q0

    ! end exact results

    ! mpi
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nprocess, ierr)

    allocate(omegaI_process(nprocess))
    ! allocate(fomegaI_process(nprocess))
    allocate(fomegaI1_process(nprocess))
    allocate(fomegaI2_process(nprocess))

    call random_seed()

    ! l = 2.0q0
    ! l = 3.0q0
    ! l = 4.0q0
    l = 5.0q0
    ! l = 6.0q0
    ! l = 7.0q0

    omegaI_start = 1.274q0

    LL = l*(l + 1.0q0)
    lambda = ((-1.0q0 + l)*(2.0q0 + l))/2.0q0

    ! Important note: It is necessary to choose the values ​​of the parameters
    ! omegaI_min and omegaI_max in such a way that they ensure that omegaI 
    ! varies always within the radius of convergence of the root. 
    ! Otherwise, the code may diverge, and the sign of a divergence
    ! of the code is that it tends to one of the endpoints of the interval,
    ! namely omegaI_min and omegaI_max.

    omegaI_result = omegaI_start
    delta_omegaI_start = 1.0q-3
    omegaI_min = omegaI_start - delta_omegaI_start
    omegaI_max = omegaI_start + delta_omegaI_start
    delta_omegaI = delta_omegaI_start
    omegaI_result_old = omegaI_result

    omegaI = omegaI_result
    call process()
    ! f_omegaI_result = fomegaI()
    f_omegaI1_result = fomegaI1()
    f_omegaI2_result = fomegaI2()

    ! print *, omegaI_result, f_omegaI_result, delta_omegaI
    ! write(*, "(3(1x,es42.34))") omegaI_result, f_omegaI_result, delta_omegaI
    write(*, "(4(1x,es42.34))") omegaI_result, f_omegaI1_result, f_omegaI2_result, delta_omegaI

    ! examine only one omegaI

    ! call store_omegaI()
    ! call store_rZmZpW()
    ! stop

    ! end examine only one omegaI

    ! Random-walk minimization

    n_try = 1

    do while(n_try < n_try_max)

    ! do while(delta_omegaI > delta_omegaI_tol)

        ! Generate random number for omegaI update
        call random_number(rand)
        omegaI = omegaI_result + delta_omegaI*(2.0q0*rand - 1.0q0)
        omegaI = max(omegaI_min, min(omegaI, omegaI_max))

        ! Perform your custom process (function)
        call process()
        ! f_omegaI = fomegaI()
        f_omegaI1 = fomegaI1()
        f_omegaI2 = fomegaI2()

        ! print *, omegaI_result, f_omegaI_result, delta_omegaI

        call MPI_Barrier(MPI_COMM_WORLD, ierr)

        ! Gather omegaI and f_omegaI values from all processes
        call MPI_GATHER(omegaI, 1, MPI_REAL16, omegaI_process, 1, MPI_REAL16, 0, MPI_COMM_WORLD, ierr)
        ! call MPI_GATHER(f_omegaI, 1, MPI_REAL16, fomegaI_process, 1, MPI_REAL16, 0, MPI_COMM_WORLD, ierr)
        call MPI_GATHER(f_omegaI1, 1, MPI_REAL16, fomegaI1_process, 1, MPI_REAL16, 0, MPI_COMM_WORLD, ierr)
        call MPI_GATHER(f_omegaI2, 1, MPI_REAL16, fomegaI2_process, 1, MPI_REAL16, 0, MPI_COMM_WORLD, ierr)

        ! Root process compares the gathered results
        if (rank == 0) then
            improved = .false.
            do j = 1, nprocess
                if ((fomegaI1_process(j) < f_omegaI1_result) .and. (fomegaI2_process(j) < f_omegaI2_result)) then
                    improved = .true.
                    omegaI_result = omegaI_process(j)
                    ! f_omegaI_result = fomegaI_process(j)
                    f_omegaI1_result = fomegaI1_process(j)
                    f_omegaI2_result = fomegaI2_process(j)
                    delta_omegaI = min(delta_omegaI_start, &
                    10.0q0*abs(omegaI_result - omegaI_result_old))
                    omegaI_result_old = omegaI_result
                end if
            end do

            ! If improvement found, reset n_try; otherwise, increment n_try
            if (improved) then
                n_try = 1
            else
                n_try = n_try + 1
            end if

            ! delta_omegaI = abs(maxval(omegaI_process(1:nprocess)) - &
            ! minval(omegaI_process(1:nprocess)))

            ! Print current try and best results
            ! print *, omegaI_result, f_omegaI_result, delta_omegaI
            ! write(*, "(3(1x,es42.34))") omegaI_result, f_omegaI_result, delta_omegaI
            write(*, "(4(1x,es42.34))") omegaI_result, f_omegaI1_result, f_omegaI2_result, delta_omegaI

        end if ! (rank == 0)

        call MPI_Barrier(MPI_COMM_WORLD, ierr)

        ! Broadcast the best omegaI_result to all processes for the next iteration
        call MPI_BCAST(omegaI_result, 1, MPI_REAL16, 0, MPI_COMM_WORLD, ierr)

        call MPI_BCAST(n_try, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

        call MPI_BCAST(delta_omegaI, 1, MPI_REAL16, 0, MPI_COMM_WORLD, ierr)

    end do

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    if (rank == 0) then

        omegaI = omegaI_result

        call process()

        ! f_omegaI = fomegaI()
        f_omegaI1 = fomegaI1()
        f_omegaI2 = fomegaI2()

        call store_omegaI()

        call store_rZmZpW()

    end if ! (rank == 0)

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    deallocate(omegaI_process)
    ! deallocate(fomegaI_process)
    deallocate(fomegaI1_process)
    deallocate(fomegaI2_process)

    ! stop

    ! mpi
    call MPI_Finalize(ierr)

contains

subroutine process()

    c_series = 0.0q0
    array_rm = 0.0q0
    array_rp = 0.0q0
    array_Zm = 0.0q0
    array_Zmp = 0.0q0
    array_Zp = 0.0q0
    array_Zpp = 0.0q0
    array_W = 0.0q0
    array_absW = 0.0q0

    omega = (0.0q0, 1.0q0)*omegaI

    call compute_c_series()

    call compute_Zm()

    call compute_Zmp()

    call compute_Zp()
    
    call compute_Zpp()

    call compute_W()

    call compute_absW()

end subroutine process

function dg(r, g, gp)

    complex(16) :: dg
    real(16) :: r
    complex(16) :: g, gp

    dg = gp

end function dg

function dgp(r, g, gp)

    complex(16) :: dgp
    real(16) :: r
    complex(16) :: g, gp

    dgp = -(((-6.0q0 + 8.0q0*r**2 - 3.0q0*l**5*r**2 - l**6*r**2 + &
    (0.0q0, 2.0q0)*omega*r*(-9.0q0 + 4.0q0*r**2) + &
    l**4*r*(-6.0q0 + 3.0q0*r + (0.0q0, 2.0q0)*omega*r**2) + &
    l**3*r*(-12.0q0 + 11.0q0*r + (0.0q0, 4.0q0)*omega*r**2) + &
    l**2*(3.0q0 + 6.0q0*r - 6.0q0*r**2 - &
    (0.0q0, 6.0q0)*omega*r**3) + &
    l*(3.0q0 + 12.0q0*r - 12.0q0*r**2 - (0.0q0, 8.0q0)*omega*r**3) &
    )*g + (3.0q0 + (-2.0q0 + l + l**2)*r)* &
    (-3.0q0 + (4.0q0 - (0.0q0, 6.0q0)*omega)*r**2 + &
    (0.0q0, 4.0q0)*omega*r**3 + &
    l*r*(3.0q0 - 2.0q0*r - (0.0q0, 2.0q0)*omega*r**2) + &
    l**2*r*(3.0q0 - 2.0q0*r - (0.0q0, 2.0q0)*omega*r**2))* &
    gp)/((-1.0q0 + r)*r*(3.0q0 + (-2.0q0 + l + l**2)*r)**2))

end function dgp

function gbc(r)

    complex(16) :: gbc
    real(16) :: r

    integer :: i

    gbc = 0.0q0

    do i = 0, n_series
        gbc = gbc + c_series(i)*r**i
    end do

end function gbc

function gpbc(r)

    complex(16) :: gpbc
    real(16) :: r

    integer :: i

    gpbc = 0.0q0

    do i = 1, n_series
        gpbc = gpbc + real(i, 16)*c_series(i)*r**(i-1)
    end do

end function gpbc

function Z(r, g)

    complex(16) :: Z
    real(16) :: r
    complex(16) :: g

    Z = (r*g)/(exp((0.0q0, 1.0q0)*omega*r)*(-1.0q0 + r)**((0.0q0, 1.0q0)*omega)* &
    (1.5q0 + r*lambda)**2)

end function Z

subroutine compute_c_series()

    c_series(0) = 1.0q0
    c_series(1) = -0.3333333333333333q0*((2.0q0 - l - l**2)*c_series(0))
    c_series(2) = -0.08333333333333333q0*((-2.0q0*l - l**2 + 2.0q0*l**3 + &
    l**4 + (0.0q0, 6.0q0)*omega)*c_series(0))
    c_series(3) = -0.027777777777777776q0* &
    (l*(1.0q0 + l)*(-2.0q0*l - l**2 + 2.0q0*l**3 + l**4 + &
    (0.0q0, 6.0q0)*omega)*c_series(0))

    do i = 4, n_series
        i_real = real(i, 16)
        c_series(i) = ((0.0q0,-0.2222222222222222q0)*(-2.0q0 + l + l**2)**2* &
        (-5.0q0 + i_real)*omega*c_series(-4 + i))/i_real**2 - &
        ((-2.0q0 + l + l**2)* &
        (2.0q0*l**3 + l**4 - l**2*(21.0q0 - 9.0q0*i_real + i_real**2) - &
        l*(22.0q0 - 9.0q0*i_real + i_real**2) + &
        2.0q0*(20.0q0 + i_real**2 + i_real*(-9.0q0 + (0.0q0, 6.0q0)*omega) - &
        (0.0q0, 18.0q0)*omega))*c_series(-3 + i))/(9.0q0*i_real**2) - &
        ((l**2*(-78.0q0 + 54.0q0*i_real - 9.0q0*i_real**2) + &
        2.0q0*l**3*(18.0q0 - 8.0q0*i_real + i_real**2) + &
        l**4*(18.0q0 - 8.0q0*i_real + i_real**2) - &
        2.0q0*l*(48.0q0 - 31.0q0*i_real + 5.0q0*i_real**2) + &
        2.0q0*(60.0q0 + 8.0q0*i_real**2 + i_real*(-46.0q0 + &
        (0.0q0, 9.0q0)*omega) - &
        (0.0q0, 9.0q0)*omega))*c_series(-2 + i))/(9.0q0*i_real**2) - &
        ((-10.0q0 + 19.0q0*i_real - 7.0q0*i_real**2 + l*(5.0q0 - &
        8.0q0*i_real + 2.0q0*i_real**2) + &
        l**2*(5.0q0 - 8.0q0*i_real + 2.0q0*i_real**2))*c_series(-1 + i))/ &
        (3.0q0*i_real**2)
    end do

end subroutine compute_c_series

subroutine compute_Zm()

    integer :: i
    real(16) :: CC

    rbc1 = 1.0q-6
    rbc2 = 1.0q0 - 1.0q-6
    dr = (rbc2 - rbc1)/real(n_array - 1, 16)

    rinitial = rbc1
    ginitial = gbc(rinitial)
    gpinitial = gpbc(rinitial)

    i = 1

    ri = rinitial
    gi = ginitial
    gpi = gpinitial

    Zi = Z(ri, gi)

    array_rm(i) = ri
    array_Zm(i) = Zi

    do i = 2, n_array

        call RK8()

        array_rm(i) = ri

        if ((isnan(real(Zi)) .or. isnan(aimag(Zi)))) then
            array_Zm(i) = 0.0q0
        else
            array_Zm(i) = Zi
        end if

    end do

    ! normalization

    ! trapezoidal rule integraration

    ! CC = (abs(array_Zm(1))**2 + abs(array_Zm(n_array))**2)/2.0q0

    ! do i = 2, n_array - 1
    !     CC = CC + abs(array_Zm(i))**2
    ! end do

    ! CC = CC*dr

    ! do i = 1, n_array
    !     array_Zm(i) = array_Zm(i)/sqrt(CC)
    ! end do

    ! test normalization
    
    ! CC = 0.0q0

    ! do i = 1, n_array
    !     CC = CC + abs(array_Zm(i))**2
    ! end do

    ! CC = CC*dr

    ! print *, "CC = ", CC

    ! end test normalization

    ! d1cubicspline = (array_Zm(2) - array_Zm(1))/(array_rm(2) - array_rm(1))
    ! d2cubicspline = (array_Zm(nr) - array_Zm(nr - 1))/ &
    ! (array_rm(nr) - array_rm(nr - 1))

    ! call spline(array_rm(1:nr), array_Zm(1:nr), nr, &
    ! d1cubicspline, d2cubicspline, array_Zm_buffer(1:nr))

end subroutine compute_Zm

! function Zm(r)

!     complex(16) :: Zm
!     real(16) :: r
!     complex(16) :: Zm_cubicspline

!     call spline_out(array_rm(1:nr), array_Zm(1:nr), array_Zm_buffer(1:r), &
!     nr, r, Zm_cubicspline)

!     Zm = Zm_cubicspline

! end function Zm

subroutine compute_Zmp()

    integer :: i

    do i = 2, n_array - 1

        array_Zmp(i) = (array_Zm(i + 1) - array_Zm(i - 1))/ &
        (array_rm(i + 1) - array_rm(i - 1))

    end do

    array_Zmp(1) = (array_Zm(2) - array_Zm(1))/ &
    (array_rm(2) - array_rm(1))

    array_Zpp(n_array) = (array_Zp(n_array) - array_Zp(n_array - 1))/ &
    (array_rm(n_array) - array_rm(n_array - 1))

    ! d1cubicspline = (array_Zmp(2) - array_Zmp(1))/(array_rm(2) - array_rm(1))
    ! d2cubicspline = (array_Zmp(nr) - array_Zmp(nr - 1))/ &
    ! (array_rm(nr) - array_rm(nr - 1))

    ! call spline(array_rm(1:nr), array_Zmp(1:nr), nr, &
    ! d1cubicspline, d2cubicspline, array_Zmp_buffer(1:nr))

end subroutine compute_Zmp

! function Zmp(r)

!     complex(16) :: Zmp
!     real(16) :: r
!     complex(16) :: Zmp_cubicspline

!     call spline_out(array_rm(1:nr), array_Zmp(1:nr), array_Zmp_buffer(1:nr), &
!     nr, r, Zmp_cubicspline)

!     Zmp = Zmp_cubicspline

! end function Zmp

subroutine compute_Zp()

    integer :: i
    real(16) :: CC

    rbc1 = 1.0q-6
    rbc2 = 1.0q0 - 1.0q-20
    dr = -abs((rbc2 - rbc1)/real(n_array - 1, 16))

    rinitial = rbc2
    ginitial = gbc(rinitial)
    gpinitial = gpbc(rinitial)

    i = n_array

    ri = rinitial
    gi = ginitial
    gpi = gpinitial

    Zi = Z(ri, gi)

    array_rp(i) = ri 
    array_Zp(i) = Zi

    do while(i > 1)

        i = i - 1

        call RK8()

        array_rp(i) = ri

        if ((isnan(real(Zi)) .or. isnan(aimag(Zi)))) then
            array_Zp(i) = 0.0q0
        else
            array_Zp(i) = Zi
        end if

    end do

    ! normalization

    ! trapezoidal rule integraration

    ! CC = (abs(array_Zp(1))**2 + abs(array_Zp(n_array))**2)/2.0q0

    ! do i = 2, n_array - 1
    ! CC = CC + abs(array_Zp(i))**2
    ! end do

    ! CC = CC*abs(dr)

    ! do i = 1, n_array
    ! array_Zp(i) = array_Zp(i)/sqrt(CC)
    ! end do

    ! test normalization
    
    ! CC = 0.0q0

    ! do i = 1, n_array
    !     if (isnan(abs(array_Zp(i)))) array_Zp(i) = (0.0q0, 0.0q0)
    !     CC = CC + abs(array_Zp(i))**2
    ! end do

    ! CC = CC*abs(dr)

    ! print *, "CC = ", CC

    ! ! end test normalization

    ! d1cubicspline = (array_Zp(2) - array_Zp(1))/(array_rp(2) - array_rp(1))
    ! d2cubicspline = (array_Zp(nr) - array_Zp(nr - 1))/ &
    ! (array_rp(nr) - array_rp(nr - 1))

    ! call spline(array_rp(1:nr), array_Zp(1:nr), nr, &
    ! d1cubicspline, d2cubicspline, array_Zp_buffer(1:nr))

end subroutine compute_Zp

! function Zp(r)

!     complex(16) :: Z
!     real(16) :: r
!     complex(16) :: Zp_cubicspline

!     call spline_out(array_rp(1:nr), array_Zp(1:nr), array_Zp_buffer(1:nr), &
!     nr, r, Zp_cubicspline)

!     Zp = Zp_cubicspline

! end function Zp

subroutine compute_Zpp()

    integer :: i

    do i = 2, n_array - 1

        array_Zpp(i) = (array_Zp(i + 1) - array_Zp(i - 1))/ &
        (array_rp(i + 1) - array_rp(i - 1))

    end do

    array_Zpp(1) = (array_Zp(2) - array_Zp(1))/ &
    (array_rp(2) - array_rp(1))

    array_Zpp(n_array) = (array_Zp(n_array) - array_Zp(n_array - 1))/ &
    (array_rp(n_array) - array_rp(n_array - 1))

    ! d1cubicspline = (array_Zpp(2) - array_Zpp(1))/(array_rp(2) - array_rp(1))
    ! d2cubicspline = (array_Zpp(nr) - array_Zpp(nr - 1))/ &
    ! (array_rp(nr) - array_rp(nr - 1))

    ! call spline(array_rp(1:nr), array_Zpp(1:nr), nr, &
    ! d1cubicspline, d2cubicspline, array_Zpp_buffer(1:nr))

end subroutine compute_Zpp

! function Zpp(r)

!     complex(16) :: Zpp
!     real(16) :: r
!     complex(16) :: Zpp_cubicspline

!     call spline_out(array_rp(1:nr), array_Zpp(1:nr), array_Zpp_buffer(1:nr), &
!     nr, r, Zpp_cubicspline)

!     Zpp = Zpp_cubicspline

! end function Zpp

subroutine compute_W()

    integer :: i

    do i = 1, n_array

        array_W(i) = array_Zm(i)*array_Zpp(i) - array_Zmp(i)*array_Zp(i)

    end do

end subroutine compute_W

subroutine compute_absW()

    integer :: i

    do i = 1, n_array

        array_absW(i) = abs(array_W(i))

    end do

end subroutine compute_absW

subroutine store_omegaI()

    open(unit=11, file="omegaI.txt", status='replace')

    ! write(11, "(3e42.32)") omegaI, f_omegaI, delta_omegaI
    write(11, "(4e42.32)") omegaI, f_omegaI1, f_omegaI2, delta_omegaI

    close(11)

end subroutine store_omegaI

subroutine store_rZmZpW()

    integer :: i

    open(unit=11, file="rZmZpW.txt", status='replace')

    do i = 1, n_array
        write(11, "(8e25.16)") array_rm(i), &
        real(array_Zm(i)), aimag(array_Zm(i)), &
        real(array_Zp(i)), aimag(array_Zp(i)), &
        real(array_W(i)), aimag(array_W(i)), &
        array_absW(i)
    end do

    close(11)

end subroutine store_rZmZpW

function fomegaI()

    real(16) :: fomegaI

    fomegaI = minval(array_absW(1:n_array))

    ! fomegaI = sum(array_absW(1:n_array))/real(n_array, 16)

end function fomegaI

function fomegaI1()

    real(16) :: fomegaI1

    fomegaI1 = minval(array_absW(1:n_array))

end function fomegaI1

function fomegaI2()

    real(16) :: fomegaI2

    fomegaI2 = sum(array_absW(1:n_array))/real(n_array, 16)

    if (isnan(fomegaI2)) fomegaI2 = fomegaI1()

end function fomegaI2

subroutine RK8()

    real(16) :: rim1
    complex(16) :: gim1, gpim1
    complex(16) :: kg01, kg02, kg03, kg04, kg05, kg06, kg07, kg08, kg09, kg10
    complex(16) :: kgp01, kgp02, kgp03, kgp04, kgp05, kgp06, kgp07, kgp08, kgp09, kgp10
    real(16) :: rtemp
    complex(16) :: gtemp, gptemp

    rim1 = ri
    gim1 = gi
    gpim1 = gpi

    kg01 = dg(rim1, gim1, gpim1)

    kgp01 = dgp(rim1, gim1, gpim1)

    rtemp = rim1 + (4.0q0/27.0q0)*dr

    gtemp = gim1 + (4.0q0/27.0q0)*dr*kg01

    gptemp = gpim1 + (4.0q0/27.0q0)*dr*kgp01

    kg02 = dg(rtemp, gtemp, gptemp)

    kgp02 = dgp(rtemp, gtemp, gptemp)

    rtemp = rim1 + (2.0q0/9.0q0)*dr

    gtemp = gim1 + (1.0q0/18.0q0)*dr*(kg01 + 3.0q0*kg02)

    gptemp = gpim1 + (1.0q0/18.0q0)*dr*(kgp01 + 3.0q0*kgp02)

    kg03 = dg(rtemp, gtemp, gptemp)

    kgp03 = dgp(rtemp, gtemp, gptemp)

    rtemp = rim1 + (1.0q0/3.0q0)*dr

    gtemp = gim1 + (1.0q0/12.0q0)*dr*(kg01 + 3.0q0*kg03)

    gptemp = gpim1 + (1.0q0/12.0q0)*dr*(kgp01 + 3.0q0*kgp03)

    kg04 = dg(rtemp, gtemp, gptemp)

    kgp04 = dgp(rtemp, gtemp, gptemp)

    rtemp = rim1 + (1.0q0/2.0q0)*dr

    gtemp = gim1 + (1.0q0/8.0q0)*dr*(kg01 + 3.0q0*kg04)

    gptemp = gpim1 + (1.0q0/8.0q0)*dr*(kgp01 + 3.0q0*kgp04)

    kg05 = dg(rtemp, gtemp, gptemp)

    kgp05 = dgp(rtemp, gtemp, gptemp)

    rtemp = rim1 + (2.0q0/3.0q0)*dr

    gtemp = gim1 + (1.0q0/54.0q0)*dr*(13.0q0*kg01 - 27.0q0*kg03 + 42.0q0*kg04 + 8.0q0*kg05)

    gptemp = gpim1 + (1.0q0/54.0q0)*dr*(13.0q0*kgp01 - 27.0q0*kgp03 + 42.0q0*kgp04 + 8.0q0*kgp05)

    kg06 = dg(rtemp, gtemp, gptemp)

    kgp06 = dgp(rtemp, gtemp, gptemp)

    rtemp = rim1 + (1.0q0/6.0q0)*dr

    gtemp = gim1 + (1.0q0/4320.0q0)*dr*(389.0q0*kg01 - 54.0q0*kg03 + 966.0q0*kg04 - 824.0q0*kg05 + 243.0q0*kg06)

    gptemp = gpim1 + (1.0q0/4320.0q0)*dr*(389.0q0*kgp01 - 54.0q0*kgp03 + 966.0q0*kgp04 - 824.0q0*kgp05 + 243.0q0*kgp06)

    kg07 = dg(rtemp, gtemp, gptemp)

    kgp07 = dgp(rtemp, gtemp, gptemp)

    rtemp = rim1 + dr

    gtemp = gim1 + (1.0q0/20.0q0)*dr*(-231.0q0*kg01 + 81.0q0*kg03 - 1164.0q0*kg04 + 656.0q0*kg05 - 122.0q0*kg06 + 800.0q0*kg07)

    gptemp = gpim1 + (1.0q0/20.0q0)*dr*(-231.0q0*kgp01 + 81.0q0*kgp03 - 1164.0q0*kgp04 + 656.0q0*kgp05 - 122.0q0*kgp06 + 800.0q0*kgp07)

    kg08 = dg(rtemp, gtemp, gptemp)

    kgp08 = dgp(rtemp, gtemp, gptemp)

    rtemp = rim1 + (5.0q0/6.0q0)*dr

    gtemp = gim1 + (1.0q0/288.0q0)*dr*(-127.0q0*kg01 + 18.0q0*kg03 - 678.0q0*kg04 + 456.0q0*kg05 - 9.0q0*kg06 + 576.0q0*kg07 + 4.0q0*kg08)

    gptemp = gpim1 + (1.0q0/288.0q0)*dr*(-127.0q0*kgp01 + 18.0q0*kgp03 - 678.0q0*kgp04 + 456.0q0*kgp05 - 9.0q0*kgp06 + 576.0q0*kgp07 + 4.0q0*kgp08)

    kg09 = dg(rtemp, gtemp, gptemp)

    kgp09 = dgp(rtemp, gtemp, gptemp)

    rtemp = rim1 + dr

    gtemp = gim1 + (1.0q0/820.0q0)*dr*(1481.0q0*kg01 - 81.0q0*kg03 + 7104.0q0*kg04 - 3376.0q0*kg05 + 72.0q0*kg06 - 5040.0q0*kg07 - 60.0q0*kg08 + 720.0q0*kg09)

    gptemp = gpim1 + (1.0q0/820.0q0)*dr*(1481.0q0*kgp01 - 81.0q0*kgp03 + 7104.0q0*kgp04 - 3376.0q0*kgp05 + 72.0q0*kgp06 - 5040.0q0*kgp07 - 60.0q0*kgp08 + 720.0q0*kgp09)

    kg10 = dg(rtemp, gtemp, gptemp)

    kgp10 = dgp(rtemp, gtemp, gptemp)

    ri = rim1 + dr

    gi = gim1 + (1.0q0/840.0q0)*dr*(41.0q0*kg01 + 27.0q0*kg04 + 272.0q0*kg05 + 27.0q0*kg06 + 216.0q0*kg07 + 216.0q0*kg09 + 41.0q0*kg10)

    gpi = gpim1 + (1.0q0/840.0q0)*dr*(41.0q0*kgp01 + 27.0q0*kgp04 + 272.0q0*kgp05 + 27.0q0*kgp06 + 216.0q0*kgp07 + 216.0q0*kgp09 + 41.0q0*kgp10)

    Zi = Z(ri, gi)

end subroutine RK8

end program bound_states_polar
