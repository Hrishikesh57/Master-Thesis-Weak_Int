program g_dirac
    use bspqed_dkb, only: solve
    use bspqed_physconsts, only: alpha, electron_reduced_compton
    use bspqed_potentials, only: sphere_potential, linear_potential
    use bspqed_quadrature, only: nodes20, weights20
    use bspqed_utilities, only: meanr, kntseq1
    use bspgf_hydrogen, only: gfactor
    implicit none

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

        double precision, dimension(:), allocatable :: knots, spectrum, coefa
        double precision, dimension(:, :), allocatable :: coefmat
        double precision :: g
        integer :: nmin, npoints, nknots, nbsplines, nspectrum, na
        type(sphere_potential) :: nuclpot
        type(linear_potential) :: magnpot

        nuclpot%zalpha = zalpha
        nuclpot%rad = rnucl

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

            allocate (coefmat(2*nbsplines, nspectrum))
            coefmat = 0.0d0

            call kntseq1(rnucl, rcavity, npoints/10, order, knots)

            call solve(knots, order, kappa, nuclpot, nodes20, weights20, spectrum, coefmat)

            na = nspectrum / 2 - nmin + nmain
            coefa = coefmat(:, na)

            g = gfactor(knots, order, kappa, coefa, coefa, magnpot, nodes20, weights20)

            write (*, *) npoints, g

            deallocate (knots, spectrum, coefa, coefmat)
        end do
        return
    end subroutine compute_gfactor

end program g_dirac
