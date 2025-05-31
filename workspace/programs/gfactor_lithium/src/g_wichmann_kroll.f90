program g_wichmann_kroll
    use bspqed_physconsts, only: alpha, electron_reduced_compton
    use bspqed_potentials, only: sphere_potential, linear_potential
    use bspqed_quadrature, only: nodes20, weights20
    use bspqed_utilities, only: gam, kntseq1
    use bspqed_wichmann_kroll, only: wichmann_kroll_potential
    use bspgf_lithium, only: fbza_el, cbza_el
    implicit none

    double precision, parameter :: pi = 4.0d0 * atan(1.0d0)

    double precision :: rms, zalpha, rnucl, rcavity
    integer :: znucl, order

    read (*, *) znucl
    read (*, *) rms

    write (*, "('# Nuclear charge number: ', I3)") znucl
    write (*, "('# Nuclear RMS in fm: ', F7.4)") rms

    order = 9
    zalpha = znucl * alpha
    rnucl = sqrt(5.0d0 / 3.0d0) * rms / electron_reduced_compton
    rcavity = gam(zalpha) / zalpha * 100.0d0

    write (*, "('# Results in Feynman gauge')")
    call compute_gfactor_feynman(order, zalpha, rnucl, rcavity)

    write (*, "('# Results in Coulomb gauge')")
    call compute_gfactor_coulomb(order, zalpha, rnucl, rcavity)

    return

contains

    subroutine compute_gfactor_feynman(order, zalpha, rnucl, rcavity)
        integer, intent(in) :: order
        double precision, intent(in) :: zalpha
        double precision, intent(in) :: rnucl
        double precision, intent(in) :: rcavity

        double precision, dimension(:), allocatable :: knots
        double precision :: a, b, e, g1, g2, g3, tot
        integer :: npoints, nknots
        type(sphere_potential) :: nuclpot
        type(linear_potential) :: magnpot
        type(wichmann_kroll_potential) :: elecpot

        nuclpot%zalpha = zalpha
        nuclpot%rad = rnucl

        elecpot%zalpha = zalpha

        do npoints = 60, 90
            nknots = npoints + 2 * (order - 1)

            allocate(knots(nknots))
            knots = 0.0d0

            call kntseq1(rnucl, rcavity, npoints/10, order, knots)

            call fbza_el(knots, order, nuclpot, magnpot, elecpot, nodes20, weights20, a, b, e, g1, g2, g3)
            a  = a  * 2.0d0 * alpha * (alpha / pi) * zalpha**3
            b  = b  * 2.0d0 * alpha * (alpha / pi) * zalpha**3
            e  = e  * 2.0d0 * alpha * (alpha / pi) * zalpha**3
            g1 = g1 * 2.0d0 * alpha * (alpha / pi) * zalpha**3
            g2 = g2 * 2.0d0 * alpha * (alpha / pi) * zalpha**3
            g3 = g3 * 2.0d0 * alpha * (alpha / pi) * zalpha**3
            tot = a + b + e + g1 + g2 + g3

            write (*, *) npoints, a, b, e, g1, g2, g3, tot

            deallocate (knots)
        end do
        return
    end subroutine compute_gfactor_feynman

    subroutine compute_gfactor_coulomb(order, zalpha, rnucl, rcavity)
        integer, intent(in) :: order
        double precision, intent(in) :: zalpha
        double precision, intent(in) :: rnucl
        double precision, intent(in) :: rcavity

        double precision, dimension(:), allocatable :: knots
        double precision :: a, b, e, g1, g2, g3, tot
        integer :: npoints, nknots
        type(sphere_potential) :: nuclpot
        type(linear_potential) :: magnpot
        type(wichmann_kroll_potential) :: elecpot

        nuclpot%zalpha = zalpha
        nuclpot%rad = rnucl

        elecpot%zalpha = zalpha

        do npoints = 60, 90
            nknots = npoints + 2 * (order - 1)

            allocate(knots(nknots))
            knots = 0.0d0

            call kntseq1(rnucl, rcavity, npoints/10, order, knots)

            call cbza_el(knots, order, nuclpot, magnpot, elecpot, nodes20, weights20, a, b, e, g1, g2, g3)
            a  = a  * 2.0d0 * alpha * (alpha / pi) * zalpha**3
            b  = b  * 2.0d0 * alpha * (alpha / pi) * zalpha**3
            e  = e  * 2.0d0 * alpha * (alpha / pi) * zalpha**3
            g1 = g1 * 2.0d0 * alpha * (alpha / pi) * zalpha**3
            g2 = g2 * 2.0d0 * alpha * (alpha / pi) * zalpha**3
            g3 = g3 * 2.0d0 * alpha * (alpha / pi) * zalpha**3
            tot = a + b + e + g1 + g2 + g3

            write (*, *) npoints, a, b, e, g1, g2, g3, tot

            deallocate (knots)
        end do
        return
    end subroutine compute_gfactor_coulomb

end program g_wichmann_kroll
