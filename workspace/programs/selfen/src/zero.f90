program zero
    use bspqed_coefficients, only: bsppbn
    use bspqed_dkb, only: solve
    use bspqed_physconsts, only: alpha, electron_reduced_compton
    use bspqed_potentials, only: sphere_potential
    use bspqed_quadrature, only: nodes20, weights20
    use bspqed_utilities, only: kntseq1
    use bspse_zeropotential, only: zeropot
    implicit none

    double precision, dimension(:, :), allocatable :: coefmat, d1, d2
    double precision, dimension(:), allocatable :: knots, spectrum
    double precision :: rms, epsabs, epsrel, zalpha, rnucl, rcavity, epsa, leading, f0, abserr
    integer :: znucl, kapa, nmain, npoints, order, nintervals, nknots, nbsplines, nspectrum, la, na
    type(sphere_potential) :: nuclpot

    read (*, *) znucl
    read (*, *) rms
    read (*, *) kapa
    read (*, *) nmain
    read (*, *) npoints
    read (*, *) epsabs
    read (*, *) epsrel

    write (*, "('# Nuclear charge number: ', I3)") znucl
    write (*, "('# Nuclear RMS in fm: ', F7.4)") rms
    write (*, "('# Dirac quantum number of reference state: ', I2)") kapa
    write (*, "('# Principal quantum number of reference state: ', I2)") nmain
    write (*, "('# Number of breakpoints: ', I3)") npoints
    write (*, "('# Expected absoulute error: ', ES7.1)") epsabs
    write (*, "('# Expected relative error: ', ES7.1)") epsrel

    order = 9
    nintervals = npoints - 1
    nknots = npoints + 2 * (order - 1)
    nbsplines = nknots - order
    nspectrum = 2 * (nbsplines - 2)

    allocate(knots(nknots))
    knots = 0.0d0

    allocate (spectrum(nspectrum))
    spectrum = 0.0d0

    allocate (coefmat(2*nbsplines, nspectrum))
    coefmat = 0.0d0

    allocate (d1(order, nintervals))
    d1 = 0.0d0

    allocate (d2(order, nintervals))
    d2 = 0.0d0

    zalpha = znucl * alpha

    rnucl = sqrt(5.0d0/3.0d0) * rms / electron_reduced_compton
    rcavity = 2.5d0 / alpha

    call kntseq1(rnucl, rcavity, npoints/10, order, knots)

    nuclpot%zalpha = zalpha
    nuclpot%rad = rnucl

    call solve(knots, order, kapa, nuclpot, nodes20, weights20, spectrum, coefmat)

    la = (abs(2 * kapa + 1) - 1) / 2
    na = nspectrum / 2 + nmain - la
    epsa = spectrum(na)
    call bsppbn(knots, order, coefmat(:, na), d1, d2)

    leading = zalpha**4 / nmain**3
    epsabs = epsabs * leading
    call zeropot(knots, order, kapa, epsa, d1, d2, f0, abserr, epsabs, epsrel)
    f0 = f0 / leading
    abserr = abserr / leading
    write (*, *) f0, abserr
    return
end program zero
