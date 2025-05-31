program g_dirac
    use bspqed_physconsts, only: alpha, electron_reduced_compton
    use bspqed_potentials, only: sphere_potential, linear_potential
    use bspqed_quadrature, only: nodes20, weights20
    use bspqed_utilities, only: gam, kntseq1
    use bspgf_boron, only: fbza, cbza
    implicit none

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
        double precision :: irr, red, tot
        integer :: npoints, nknots
        type(sphere_potential) :: nuclpot
        type(linear_potential) :: magnpot

        nuclpot%zalpha = zalpha
        nuclpot%rad = rnucl

        do npoints = 60, 90
            nknots = npoints + 2 * (order - 1)

            allocate(knots(nknots))
            knots = 0.0d0

            call kntseq1(rnucl, rcavity, npoints/10, order, knots)

            call fbza(knots, order, nuclpot, magnpot, nodes20, weights20, irr, red)
            irr = irr * 2.0d0 * alpha
            red = red * 2.0d0 * alpha
            tot = irr + red

            write (*, *) npoints, irr, red, tot

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
        double precision :: irr, red, tot
        integer :: npoints, nknots
        type(sphere_potential) :: nuclpot
        type(linear_potential) :: magnpot

        nuclpot%zalpha = zalpha
        nuclpot%rad = rnucl

        do npoints = 60, 90
            nknots = npoints + 2 * (order - 1)

            allocate(knots(nknots))
            knots = 0.0d0

            call kntseq1(rnucl, rcavity, npoints/10, order, knots)

            call cbza(knots, order, nuclpot, magnpot, nodes20, weights20, irr, red)
            irr = irr * 2.0d0 * alpha
            red = red * 2.0d0 * alpha
            tot = irr + red

            write (*, *) npoints, irr, red, tot

            deallocate (knots)
        end do
        return
    end subroutine compute_gfactor_coulomb

end program g_dirac
