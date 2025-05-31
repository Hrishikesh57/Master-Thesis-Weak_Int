program many
    use bspqed_discretization, only: nodwgt, storewf, storefunc, matelems
    use bspqed_dkb, only: solve
    use bspqed_function, only: function_type
    use bspqed_physconsts, only: alpha, electron_reduced_compton
    use bspqed_potentials, only: sphere_potential, free_potential
    use bspqed_quadrature, only: nodes10, weights10, nodes20, weights20
    use bspqed_utilities, only: kntseq1
    use bspse_manypotential, only: manypot_low, manypot_high
    implicit none

    double precision, dimension(:, :), allocatable :: coefmat, flow, fhgh
    double precision, dimension(:), allocatable :: knots, speca, x1, w1, x2, w2, v1, ga1, fa1, ga2, fa2, xtmp,&
    & knots_b, knots_f
    double precision :: rms, epsabs, epsrel, zalpha, rnucl, rcavity, eps1s, leading, fln, flp, fhn, fhp
    integer :: znucl, kapa, nmain, npoints, modkapilow, modkapiupp, ngl, order, nknots, nintervals, nbsplines,&
    & nspectrum, order_b, npoints_b, nknots_b, order_f, npoints_f, nknots_f, la, na, modkapi
    type(sphere_potential) :: nuclpot
    type(free_potential) :: freepot

    read (*, *) znucl
    read (*, *) rms
    read (*, *) kapa
    read (*, *) nmain
    read (*, *) npoints
    read (*, *) epsabs
    read (*, *) epsrel
    read (*, *) modkapilow
    read (*, *) modkapiupp

    write (*, "('# Nuclear charge number: ', I3)") znucl
    write (*, "('# Nuclear RMS in fm: ', F7.4)") rms
    write (*, "('# Dirac quantum number of reference state: ', I2)") kapa
    write (*, "('# Principal quantum number of reference state: ', I2)") nmain
    write (*, "('# Number of breakpoints: ', I3)") npoints
    write (*, "('# Expected absoulute error: ', ES7.1)") epsabs
    write (*, "('# Expected relative error: ', ES7.1)") epsrel
    write (*, "('# Lower Dirac quantum number of sequence: ', I2)") modkapilow
    write (*, "('# Upper Dirac quantum number of sequence: ', I2)") modkapiupp

    ngl = size(nodes10)
    order = 9
    nknots = npoints + 2 * (order - 1)
    nintervals = npoints - 1
    nbsplines = nknots - order
    nspectrum = 2 * (nbsplines - 2)

    order_b = order
    npoints_b = npoints
    nknots_b = npoints_b + 2 * (order_b - 1)

    order_f = order
    npoints_f = npoints_b - 20
    nknots_f = npoints_f + 2 * (order_f - 1)

    allocate(knots(nknots))
    knots = 0.0d0

    allocate (speca(nspectrum))
    speca = 0.0d0

    allocate (coefmat(2*nbsplines, nspectrum))
    coefmat = 0.0d0

    allocate (x1(ngl * nintervals))
    x1 = 0.0d0

    allocate (w1(ngl * nintervals))
    w1 = 0.0d0

    allocate (x2(ngl**2 * nintervals))
    x2 = 0.0d0

    allocate (w2(ngl**2 * nintervals))
    w2 = 0.0d0

    allocate (v1(ngl * nintervals))
    v1 = 0.0d0

    allocate (ga1(ngl * nintervals))
    ga1 = 0.0d0

    allocate (fa1(ngl * nintervals))
    fa1 = 0.0d0

    allocate (ga2(ngl**2 * nintervals))
    ga2 = 0.0d0

    allocate (fa2(ngl**2 * nintervals))
    fa2 = 0.0d0

    allocate (xtmp(ngl * nintervals + 1))
    xtmp = 0.0d0

    allocate (knots_b(nknots_b))
    knots_b = 0.0d0

    allocate (knots_f(nknots_f))
    knots_f = 0.0d0

    allocate (flow(modkapilow:modkapiupp, 2))
    flow = 0.0d0

    allocate (fhgh(modkapilow:modkapiupp, 2))
    fhgh = 0.0d0

    zalpha = znucl * alpha

    rnucl = sqrt(5.0d0/3.0d0) * rms / electron_reduced_compton
    rcavity = 2.5d0 / alpha

    call kntseq1(rnucl, rcavity, npoints/10, order, knots)

    call nodwgt(knots(order:(nknots-order+1)), nodes10, weights10, x1, w1)
    xtmp(1) = 0.0d0
    xtmp(2:(ngl * nintervals + 1)) = x1
    call nodwgt(xtmp, nodes10, weights10, x2, w2)

    nuclpot%zalpha = zalpha
    nuclpot%rad = rnucl
    call storefunc(nuclpot, x1, v1)

    call solve(knots, order, -1, nuclpot, nodes20, weights20, speca, coefmat)
    na = nspectrum / 2 + 1
    eps1s = speca(na)

    la = (abs(2 * kapa + 1) - 1) / 2
    na = nspectrum / 2 + nmain - la
    call solve(knots, order, kapa, nuclpot, nodes20, weights20, speca, coefmat)
    call storewf(x1, knots, order, kapa, coefmat(:, na), ga1, fa1)
    call storewf(x2, knots, order, kapa, coefmat(:, na), ga2, fa2)

    call kntseq1(rnucl, rcavity, npoints/10, order_b, knots_b)
    call kntseq1(rnucl, rcavity, npoints/10, order_f, knots_f)

    leading = zalpha**4 / nmain**3
    epsabs = epsabs * leading

    !$OMP PARALLEL DO
    do modkapi = modkapilow, modkapiupp
        call comp_low(modkapi, knots_b, order_b, knots_f, order_f, nuclpot, freepot, x1, w1, v1, na, kapa, speca,&
        & ga1, fa1, eps1s, epsabs, epsrel, flow(modkapi, 1), flow(modkapi, 2))
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO
    do modkapi = modkapilow, modkapiupp
        call comp_high(modkapi, knots_b, order_b, knots_f, order_f, nuclpot, freepot, x1, w1, x2, w2, v1, na, kapa&
        &, speca, ga1, fa1, ga2, fa2, eps1s, epsabs, epsrel, fhgh(modkapi, 1), fhgh(modkapi, 2))
    end do
    !$OMP END PARALLEL DO

    do modkapi = modkapilow, modkapiupp
        fln = flow(modkapi, 1) / leading
        flp = flow(modkapi, 2) / leading
        fhn = fhgh(modkapi, 1) / leading
        fhp = fhgh(modkapi, 2) / leading
        write (*, *) modkapi, fln, flp, fhn, fhp
    end do
    return

contains

    subroutine comp_low(modkapi, knots_b, order_b, knots_f, order_f, nuclpot, freepot, x1, w1, v1, na, kapa, speca&
    &, ga1, fa1, eps1s, epsabs, epsrel, empn, empp)
        integer, intent(in) :: modkapi
        double precision, dimension(:), intent(in) :: knots_b
        integer, intent(in) :: order_b
        double precision, dimension(:), intent(in) :: knots_f
        integer, intent(in) :: order_f
        class(function_type), intent(in) :: nuclpot
        class(function_type), intent(in) :: freepot
        double precision, dimension(:), intent(in) :: x1
        double precision, dimension(:), intent(in) :: w1
        double precision, dimension(:), intent(in) :: v1
        integer, intent(in) :: na
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: speca
        double precision, dimension(:), intent(in) :: ga1
        double precision, dimension(:), intent(in) :: fa1
        double precision, intent(in) :: eps1s
        double precision, intent(in) :: epsabs
        double precision, intent(in) :: epsrel

        double precision, intent(out) :: empn
        double precision, intent(out) :: empp

        double precision, dimension(2*(size(knots_b)-order_b), 2*(size(knots_b)-order_b-2)) :: coefmat_b
        double precision, dimension(2*(size(knots_f)-order_f), 2*(size(knots_f)-order_f-2)) :: coefmat_f
        double precision, dimension(size(x1), 2*(size(knots_b)-order_b-2)) :: gb1, fb1
        double precision, dimension(size(x1), 2*(size(knots_f)-order_f-2)) :: gf1, ff1
        double precision, dimension(2*(size(knots_f)-order_f-2), 2*(size(knots_b)-order_b-2)) :: vmat
        double precision, dimension(2*(size(knots_b)-order_b-2)) :: specb
        double precision, dimension(2*(size(knots_f)-order_f-2)) :: specf
        double precision, dimension(-1:1) :: emps
        double precision :: emp
        integer :: nspectrum_b, nspectrum_f, sgnkapi, kapi, n

        nspectrum_b = 2 * (size(knots_b) - order_b - 2)
        nspectrum_f = 2 * (size(knots_f) - order_f - 2)
        do sgnkapi = -1, 1, 2
            kapi = sgnkapi * modkapi
            call solve(knots_b, order_b, kapi, nuclpot, nodes20, weights20, specb, coefmat_b)
            do n = 1, nspectrum_b
                call storewf(x1, knots_b, order_b, kapi, coefmat_b(:, n), gb1(:, n), fb1(:, n))
            end do
            call solve(knots_f, order_f, kapi, freepot, nodes20, weights20, specf, coefmat_f)
            do n = 1, nspectrum_f
                call storewf(x1, knots_f, order_f, kapi, coefmat_f(:, n), gf1(:, n), ff1(:, n))
            end do
            call matelems(w1, v1, gf1, ff1, gb1, fb1, vmat)
            emp = manypot_low(na, kapa, speca, ga1, fa1, kapi, specb, specf, gf1, ff1, vmat, eps1s, x1, w1, epsabs&
            &, epsrel)
            emps(sgnkapi) = emp
        end do
        empn = emps(-1)
        empp = emps(+1)
        return
    end subroutine comp_low

    subroutine comp_high(modkapi, knots_b, order_b, knots_f, order_f, nuclpot, freepot, x1, w1, x2, w2, v1, na,&
    & kapa, speca, ga1, fa1, ga2, fa2, eps1s, epsabs, epsrel, empn, empp)
        integer, intent(in) :: modkapi
        double precision, dimension(:), intent(in) :: knots_b
        integer, intent(in) :: order_b
        double precision, dimension(:), intent(in) :: knots_f
        integer, intent(in) :: order_f
        class(function_type), intent(in) :: nuclpot
        class(function_type), intent(in) :: freepot
        double precision, dimension(:), intent(in) :: x1
        double precision, dimension(:), intent(in) :: w1
        double precision, dimension(:), intent(in) :: x2
        double precision, dimension(:), intent(in) :: w2
        double precision, dimension(:), intent(in) :: v1
        integer, intent(in) :: na
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: speca
        double precision, dimension(:), intent(in) :: ga1
        double precision, dimension(:), intent(in) :: fa1
        double precision, dimension(:), intent(in) :: ga2
        double precision, dimension(:), intent(in) :: fa2
        double precision, intent(in) :: eps1s
        double precision, intent(in) :: epsabs
        double precision, intent(in) :: epsrel

        double precision, intent(out) :: empn
        double precision, intent(out) :: empp

        double precision, dimension(2*(size(knots_b)-order_b), 2*(size(knots_b)-order_b-2)) :: coefmat_b
        double precision, dimension(2*(size(knots_f)-order_f), 2*(size(knots_f)-order_f-2)) :: coefmat_f
        double precision, dimension(size(x1), 2*(size(knots_b)-order_b-2)) :: gb1, fb1
        double precision, dimension(size(x1), 2*(size(knots_f)-order_f-2)) :: gf1, ff1
        double precision, dimension(size(x2), 2*(size(knots_b)-order_b-2)) :: gb2, fb2
        double precision, dimension(size(x2), 2*(size(knots_f)-order_f-2)) :: gf2, ff2
        double precision, dimension(2*(size(knots_f)-order_f-2), 2*(size(knots_b)-order_b-2)) :: vmat
        double precision, dimension(2*(size(knots_b)-order_b-2)) :: specb
        double precision, dimension(2*(size(knots_f)-order_f-2)) :: specf
        double precision, dimension(-1:1) :: emps
        double precision :: emp
        integer :: nspectrum_b, nspectrum_f, sgnkapi, kapi, n

        nspectrum_b = 2 * (size(knots_b) - order_b - 2)
        nspectrum_f = 2 * (size(knots_f) - order_f - 2)
        do sgnkapi = -1, 1, 2
            kapi = sgnkapi * modkapi
            call solve(knots_b, order_b, kapi, nuclpot, nodes20, weights20, specb, coefmat_b)
            do n = 1, nspectrum_b
                call storewf(x1, knots_b, order_b, kapi, coefmat_b(:, n), gb1(:, n), fb1(:, n))
                call storewf(x2, knots_b, order_b, kapi, coefmat_b(:, n), gb2(:, n), fb2(:, n))
            end do
            call solve(knots_f, order_f, kapi, freepot, nodes20, weights20, specf, coefmat_f)
            do n = 1, nspectrum_f
                call storewf(x1, knots_f, order_f, kapi, coefmat_f(:, n), gf1(:, n), ff1(:, n))
                call storewf(x2, knots_f, order_f, kapi, coefmat_f(:, n), gf2(:, n), ff2(:, n))
            end do
            call matelems(w1, v1, gf1, ff1, gb1, fb1, vmat)
            emp = manypot_high(na, kapa, speca, ga1, fa1, ga2, fa2, kapi, specb, specf, gf1, ff1, gf2, ff2, vmat,&
            & eps1s, x1, w1, x2, w2, epsabs, epsrel)
            emps(sgnkapi) = emp
        end do
        empn = emps(-1)
        empp = emps(+1)
        return
    end subroutine comp_high

end program many
