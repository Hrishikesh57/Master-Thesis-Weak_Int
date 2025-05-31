program g_wichmann_kroll
    use bspqed_dkb, only: solve
    use bspqed_physconsts, only: alpha, electron_reduced_compton
    use bspqed_potentials, only: sphere_potential, linear_potential
    use bspqed_quadrature, only: nodes20, weights20
    use bspqed_utilities, only: meanr, kntseq1, perturbation
    use bspqed_wichmann_kroll, only: wichmann_kroll_potential
    use bspgf_hydrogen, only: gfactor
    implicit none

    double precision, parameter :: pi = 4.0d0 * atan(1.0d0)

    double precision :: rms, zalpha, rnucl, rcavity
    integer :: znucl, kappa, nmain, order

    read (*, *) znucl
    read (*, *) rms
    read (*, *) kappa
    read (*, *) nmain

    write (*, "('# Nuclear charge number: ', I3)") znucl
    write (*, "('# Nuclear RMS in fm: ', F7.4)") rms
    write (*, "('# Dirac quantum number of reference state: ', I2)") kappa
    write (*, "('# Principal quantum number of reference state: ', I2)") nmain

    order = 9
    zalpha = znucl * alpha
    rnucl = sqrt(5.0d0 / 3.0d0) * rms / electron_reduced_compton
    rcavity = meanr(zalpha, kappa, nmain) * 100.0d0

    call compute_gfactor(order, kappa, nmain, zalpha, rnucl, rcavity)
    return

contains

    subroutine compute_gfactor(order, kappa, nmain, zalpha, rnucl, rcavity)
        integer, intent(in) :: order
        integer, intent(in) :: kappa
        integer, intent(in) :: nmain
        double precision, intent(in) :: zalpha
        double precision, intent(in) :: rnucl
        double precision, intent(in) :: rcavity

        double precision, dimension(:), allocatable :: knots, spectrum, coefa, wichcoef
        double precision, dimension(:, :), allocatable :: coefmat
        double precision :: g
        integer :: nmin, npoints, nknots, nbsplines, nspectrum, na
        type(sphere_potential) :: nuclpot
        type(linear_potential) :: magnpot
        type(wichmann_kroll_potential) :: wichpot

        nuclpot%zalpha = zalpha
        nuclpot%rad = rnucl

        wichpot%zalpha = zalpha

        nmin = (abs(2 * kappa + 1) - 1) / 2
        do npoints = 60, 90
            nknots = npoints + 2 * (order - 1)
            nbsplines = nknots - order
            nspectrum = 2 * (nbsplines - 2)

            allocate(knots(nknots))
            knots = 0.0d0

            allocate (spectrum(nspectrum))
            spectrum = 0.0d0

            allocate (coefa(2*nbsplines))
            coefa = 0.0d0

            allocate (wichcoef(2*nbsplines))
            wichcoef = 0.0d0

            allocate (coefmat(2*nbsplines, nspectrum))
            coefmat = 0.0d0

            call kntseq1(rnucl, rcavity, npoints/10, order, knots)

            call solve(knots, order, kappa, nuclpot, nodes20, weights20, spectrum, coefmat)

            na = nspectrum / 2 - nmin + nmain
            coefa = coefmat(:, na)

            call perturbation(knots, order, na, kappa, spectrum, coefmat, wichpot, nodes20, weights20, wichcoef)

            g = gfactor(knots, order, kappa, coefa, wichcoef, magnpot, nodes20, weights20) * 2.0d0
            g = g * (alpha / pi) * zalpha**3

            write (*, *) npoints, g

            deallocate (knots, spectrum, coefa, wichcoef, coefmat)
        end do
        return
    end subroutine compute_gfactor

end program g_wichmann_kroll
