program dhfcb
    use bspwf_constants, only: bohr_radius, ialpha
    use bspwf_dirac, only: dirac
    use bspwf_function, only: function_type
    use bspwf_hartree_fock, only: dhfc_dir_matrix, dhfcb_matrix
    use bspwf_knots, only: knots_01
    use bspwf_potentials, only: fermi_potential
    use bspwf_quadrature, only: nodes20, weights20
    implicit none

    double precision, parameter :: pi = 4.0d0 * atan(1.0d0)

    double precision, dimension(:, :), allocatable :: core_coefficients, valence_coefficients
    double precision, dimension(:), allocatable :: knots, core_energies, valence_energies
    double precision :: rms, afermi, cfermi, rcavity, epsb, epsv
    real :: start, finish
    integer, dimension(:), allocatable :: core_kappas, core_nmains, valence_kappas, valence_nmains
    integer :: znucl, ncores, nvalences, i, order, npoints, nknots, nbsplines, kapb, nmnb, kapv, nmnv
    type(fermi_potential) :: nuclpot

    read (*, *) znucl
    read (*, *) rms

    read (*, *) ncores
    allocate (core_kappas(ncores))
    core_kappas = 0
    allocate (core_nmains(ncores))
    core_nmains = 0
    do i = 1, ncores
        read (*, *) core_nmains(i), core_kappas(i)
    end do

    read (*, *) nvalences
    allocate (valence_kappas(nvalences))
    valence_kappas = 0
    allocate (valence_nmains(nvalences))
    valence_nmains = 0
    do i = 1, nvalences
        read (*, *) valence_nmains(i), valence_kappas(i)
    end do

    write (*, *) "# ", repeat("=", 86)
    write (*, *) znucl
    write (*, *) rms
    do i = 1, ncores
        write (*, *) i, core_nmains(i), core_kappas(i)
    end do
    do i = 1, nvalences
        write (*, *) i, valence_nmains(i), valence_kappas(i)
    end do
    write (*, *) "# ", repeat("-", 86)

    order = 9
    npoints = 100
    nknots = npoints + 2 * (order - 1)
    nbsplines = nknots - order

    allocate (knots(nknots))
    knots = 0.0d0

    afermi = 2.3d0 / (4.0d0 * log(3.0d0)) / bohr_radius
    cfermi = sqrt((5.0d0 * rms**2 - 7.0d0 * (pi * afermi)**2) / 3.0d0) / bohr_radius

    rcavity = 50.0d0
    call knots_01(cfermi, rcavity, npoints/10, order, knots)

    nuclpot%znucl = znucl
    nuclpot%a = afermi
    nuclpot%c = cfermi

    allocate (core_energies(ncores))
    core_energies = 0.0d0

    allocate (core_coefficients(2*nbsplines, ncores))
    core_coefficients = 0.0d0

    call cpu_time(start)
    call starting_solutions(knots, order, nuclpot, core_kappas, core_nmains, nodes20, weights20, core_energies,&
    & core_coefficients, print_output=.true.)
    call cpu_time(finish)

    write (*, *) "# ", repeat("-", 86)
    write (*, *) "# Calculation of starting approximations done!"
    write (*, *) "# Total duration in seconds: ", finish - start
    write (*, *) "# Energies of the core states are ..."
    write (*, *) "# ", repeat("-", 86)
    do i = 1, ncores
        kapb = core_kappas(i)
        nmnb = core_nmains(i)
        epsb = core_energies(i)
        write (*, *) nmnb, kapb, epsb
    end do
    write (*, *) "# ", repeat("-", 86)

    call cpu_time(start)
    call refining_solutions(knots, order, nuclpot, core_kappas, core_nmains, nodes20, weights20, core_energies,&
    & core_coefficients, print_output=.true.)
    call cpu_time(finish)

    write (*, *) "# ", repeat("-", 86)
    write (*, *) "# Refining the solutions done!"
    write (*, *) "# Total duration in seconds: ", finish - start
    write (*, *) "# Energies of the core states are ..."
    write (*, *) "# ", repeat("-", 86)
    do i = 1, ncores
        kapb = core_kappas(i)
        nmnb = core_nmains(i)
        epsb = core_energies(i)
        write (*, *) nmnb, kapb, epsb
    end do
    write (*, *) "# ", repeat("-", 86)

    allocate (valence_energies(nvalences))
    valence_energies = 0.0d0

    allocate (valence_coefficients(2*nbsplines, nvalences))
    valence_coefficients = 0.0d0

    call valence_solutions(knots, order, nuclpot, core_kappas, core_coefficients, valence_kappas, valence_nmains,&
    & nodes20, weights20, valence_energies, valence_coefficients)

    write (*, *) "# ", repeat("-", 86)
    write (*, *) "# Energies of the valence states are ..."
    write (*, *) "# ", repeat("-", 86)
    do i = 1, nvalences
        kapv = valence_kappas(i)
        nmnv = valence_nmains(i)
        epsv = valence_energies(i)
        write (*, *) nmnv, kapv, epsv
    end do
    write (*, *) "# ", repeat("=", 86)

    return

contains

    subroutine starting_solutions(knots, order, potential, core_kappas, core_nmains, nodes, weights, core_energies&
    &, core_coefficients, epsabs, epsrel, maxiter, print_output)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        class(function_type), intent(in) :: potential
        integer, dimension(:), intent(in) :: core_kappas
        integer, dimension(:), intent(in) :: core_nmains
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision, dimension(:), intent(out) :: core_energies
        double precision, dimension(:, :), intent(out) :: core_coefficients

        double precision, intent(in), optional :: epsabs
        double precision, intent(in), optional :: epsrel
        integer, intent(in), optional :: maxiter
        logical, intent(in), optional :: print_output

        double precision, dimension(2*(size(knots)-order-2), 2*(size(knots)-order-2), size(core_kappas)) :: hfdir,&
        & hfdir_curr
        double precision, dimension(2*(size(knots)-order-2), 2*(size(knots)-order-2)) :: ahf
        double precision, dimension(2*(size(knots)-order), 2*(size(knots)-order-2)) :: coefmat
        double precision, dimension(size(core_coefficients, 1), size(core_coefficients, 2)) :: prev_core_coefs
        double precision, dimension(size(core_kappas), 2) :: epsbs
        double precision, dimension(2*(size(knots)-order-2)) :: spectrum
        double precision, dimension(2*(size(knots)-order)) :: coefb, coefc
        double precision, dimension(size(core_kappas)) :: relerrs
        double precision :: tolabs, tolrel, epsb, abserr, tolerance, relerr
        integer, dimension(1) :: temp
        integer :: maxstep, ncores, nspectrum, step, i, nelectrons, j, kapb, nmnb, lb, nb, kapc, imax
        logical :: show

        tolabs = 0.0d0
        if (present(epsabs)) tolabs = epsabs

        tolrel = 5.0d-04
        if (present(epsrel)) tolrel = epsrel

        maxstep = 50
        if (present(maxiter)) maxstep = maxiter

        show = .false.
        if (present(print_output)) show = print_output

        ncores = size(core_kappas)
        nspectrum = 2 * (size(knots) - order - 2)

        step = 1
        do i = 1, ncores
            kapb = core_kappas(i)
            nmnb = core_nmains(i)

            hfdir(:, :, i) = 0.0d0

            call dirac(knots, order, kapb, potential, nodes, weights, spectrum, coefmat, hfmat=hfdir(:, :, i))

            lb = (abs(2 * kapb + 1) - 1) / 2
            nb = nspectrum / 2 + nmnb - lb
            epsb = spectrum(nb) - ialpha**2
            coefb = coefmat(:, nb)
            core_coefficients(:, i) = coefb

            epsbs(i, 1) = epsb

            if (show) write (*, *) step, nmnb, kapb, epsb
        end do
        if (show) write (*, *) ""

        step = 2
        prev_core_coefs = core_coefficients
        do i = 1, ncores
            kapb = core_kappas(i)
            nmnb = core_nmains(i)

            hfdir_curr(:, :, i) = 0.0d0
            nelectrons = 0
            do j = 1, ncores
                kapc = core_kappas(j)
                coefc = prev_core_coefs(:, j)
                call dhfc_dir_matrix(knots, order, kapb, kapc, coefc, nodes, weights, ahf)
                hfdir_curr(:, :, i) = hfdir_curr(:, :, i) + 2 * abs(kapc) * ahf
                nelectrons = nelectrons + 2 * abs(kapc)
            end do
            hfdir_curr(:, :, i) = (nelectrons - 1) * hfdir_curr(:, :, i) / nelectrons

            hfdir(:, :, i) = 0.5d0 * hfdir_curr(:, :, i) + 0.5d0 * hfdir(:, :, i)
            call dirac(knots, order, kapb, potential, nodes, weights, spectrum, coefmat, hfmat=hfdir(:, :, i))

            lb = (abs(2 * kapb + 1) - 1) / 2
            nb = nspectrum / 2 + nmnb - lb
            epsb = spectrum(nb) - ialpha**2
            coefb = coefmat(:, nb)
            core_coefficients(:, i) = coefb

            epsbs(i, 2) = epsb

            if (show) write (*, *) step, nmnb, kapb, epsb
        end do
        if (show) write (*, *) ""

        do step = 3, maxstep
            do i = 1, ncores
                relerrs(i) = abs((epsbs(i, 1) - epsbs(i, 2)) / epsbs(i, 2))
            end do
            temp = maxloc(relerrs)
            imax = temp(1)

            abserr = abs(epsbs(imax, 2) - epsbs(imax, 1))
            tolerance = max(tolabs, tolrel * abs(epsbs(imax, 2)))
            if (abserr <= tolerance) exit

            prev_core_coefs = core_coefficients
            do i = 1, ncores
                kapb = core_kappas(i)
                nmnb = core_nmains(i)

                hfdir_curr(:, :, i) = 0.0d0
                nelectrons = 0
                do j = 1, ncores
                    kapc = core_kappas(j)
                    coefc = prev_core_coefs(:, j)
                    call dhfc_dir_matrix(knots, order, kapb, kapc, coefc, nodes, weights, ahf)
                    hfdir_curr(:, :, i) = hfdir_curr(:, :, i) + 2 * abs(kapc) * ahf
                    nelectrons = nelectrons + 2 * abs(kapc)
                end do
                hfdir_curr(:, :, i) = (nelectrons - 1) * hfdir_curr(:, :, i) / nelectrons

                hfdir(:, :, i) = 0.5d0 * hfdir_curr(:, :, i) + 0.5d0 * hfdir(:, :, i)
                call dirac(knots, order, kapb, potential, nodes, weights, spectrum, coefmat, hfmat=hfdir(:, :, i))

                lb = (abs(2 * kapb + 1) - 1) / 2
                nb = nspectrum / 2 + nmnb - lb
                epsb = spectrum(nb) - ialpha**2
                coefb = coefmat(:, nb)
                core_coefficients(:, i) = coefb

                epsbs(i, 1) = epsbs(i, 2)
                epsbs(i, 2) = epsb

                relerr = abs((epsbs(i, 1) - epsbs(i, 2)) / epsbs(i, 2))

                if (show) write (*, *) step, nmnb, kapb, epsb, relerr
            end do
            if (show) write (*, *) ""
        end do

        core_energies = epsbs(:, 2)

        return
    end subroutine starting_solutions

    subroutine refining_solutions(knots, order, potential, core_kappas, core_nmains, nodes, weights, core_energies&
    &, core_coefficients, epsabs, epsrel, maxiter, print_output)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        class(function_type), intent(in) :: potential
        integer, dimension(:), intent(in) :: core_kappas
        integer, dimension(:), intent(in) :: core_nmains
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision, dimension(:), intent(in out) :: core_energies
        double precision, dimension(:, :), intent(in out) :: core_coefficients

        double precision, intent(in), optional :: epsabs
        double precision, intent(in), optional :: epsrel
        integer, intent(in), optional :: maxiter
        logical, intent(in), optional :: print_output

        double precision, dimension(2*(size(knots)-order-2), 2*(size(knots)-order-2)) :: hfmat, ahf
        double precision, dimension(2*(size(knots)-order), 2*(size(knots)-order-2)) :: coefmat
        double precision, dimension(size(core_coefficients, 1), size(core_coefficients, 2)) :: prev_core_coefs
        double precision, dimension(size(core_kappas), 2) :: epsbs
        double precision, dimension(2*(size(knots)-order-2)) :: spectrum
        double precision, dimension(2*(size(knots)-order)) :: coefb, coefc
        double precision, dimension(size(core_kappas)) :: relerrs
        double precision :: tolabs, tolrel, epsb, abserr, tolerance, relerr
        integer, dimension(1) :: temp
        integer :: maxstep, ncores, nspectrum, step, i, kapb, nmnb, j, kapc, lb, nb, imax
        logical :: show

        tolabs = 0.0d0
        if (present(epsabs)) tolabs = epsabs

        tolrel = 5.0d-08
        if (present(epsrel)) tolrel = epsrel

        maxstep = 100
        if (present(maxiter)) maxstep = maxiter

        show = .false.
        if (present(print_output)) show = print_output

        ncores = size(core_kappas)
        nspectrum = 2 * (size(knots) - order - 2)

        do step = 1, 2
            prev_core_coefs = core_coefficients
            do i = 1, ncores
                kapb = core_kappas(i)
                nmnb = core_nmains(i)

                hfmat = 0.0d0
                do j = 1, ncores
                    kapc = core_kappas(j)
                    coefc = prev_core_coefs(:, j)
                    call dhfcb_matrix(knots, order, kapb, kapc, coefc, nodes, weights, ahf)
                    hfmat = hfmat + 2 * abs(kapc) * ahf
                end do

                call dirac(knots, order, kapb, potential, nodes, weights, spectrum, coefmat, hfmat=hfmat)

                lb = (abs(2 * kapb + 1) - 1) / 2
                nb = nspectrum / 2 + nmnb - lb
                epsb = spectrum(nb) - ialpha**2
                coefb = coefmat(:, nb)
                core_coefficients(:, i) = coefb

                epsbs(i, step) = epsb

                if (show) write (*, *) step, nmnb, kapb, epsb
            end do
            if (show) write (*, *) ""
        end do

        do step = 3, maxstep
            do i = 1, ncores
                relerrs(i) = abs((epsbs(i, 1) - epsbs(i, 2)) / epsbs(i, 2))
            end do
            temp = maxloc(relerrs)
            imax = temp(1)

            abserr = abs(epsbs(imax, 2) - epsbs(imax, 1))
            tolerance = max(tolabs, tolrel * abs(epsbs(imax, 2)))
            if (abserr <= tolerance) exit

            prev_core_coefs = core_coefficients

            kapb = core_kappas(imax)
            nmnb = core_nmains(imax)

            hfmat = 0.0d0
            do j = 1, ncores
                kapc = core_kappas(j)
                coefc = prev_core_coefs(:, j)
                call dhfcb_matrix(knots, order, kapb, kapc, coefc, nodes, weights, ahf)
                hfmat = hfmat + 2 * abs(kapc) * ahf
            end do

            call dirac(knots, order, kapb, potential, nodes, weights, spectrum, coefmat, hfmat=hfmat)

            lb = (abs(2 * kapb + 1) - 1) / 2
            nb = nspectrum / 2 + nmnb - lb
            epsb = spectrum(nb) - ialpha**2
            coefb = coefmat(:, nb)
            core_coefficients(:, imax) = coefb

            epsbs(imax, 1) = epsbs(imax, 2)
            epsbs(imax, 2) = epsb

            relerr = abs((epsbs(imax, 1) - epsbs(imax, 2)) / epsbs(imax, 2))

            if (show) write (*, *) step, nmnb, kapb, epsb, relerr
        end do
        if (show) write (*, *) ""

        core_energies = epsbs(:, 2)

        return
    end subroutine refining_solutions

    subroutine valence_solutions(knots, order, potential, core_kappas, core_coefficients, valence_kappas,&
    & valence_nmains, nodes, weights, valence_energies, valence_coefficients)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        class(function_type), intent(in) :: potential
        integer, dimension(:), intent(in) :: core_kappas
        double precision, dimension(:, :), intent(in) :: core_coefficients
        integer, dimension(:), intent(in) :: valence_kappas
        integer, dimension(:), intent(in) :: valence_nmains
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision, dimension(:), intent(out) :: valence_energies
        double precision, dimension(:, :), intent(out) :: valence_coefficients

        double precision, dimension(2*(size(knots)-order-2), 2*(size(knots)-order-2)) :: hfmat, ahf
        double precision, dimension(2*(size(knots)-order), 2*(size(knots)-order-2)) :: coefmat
        double precision, dimension(2*(size(knots)-order-2)) :: spectrum
        double precision, dimension(2*(size(knots)-order)) :: coefb
        integer :: ncores, nvalences, nspectrum, i, kapa, nmna, j, kapb, la, na

        ncores = size(core_kappas)
        nvalences = size(valence_kappas)
        nspectrum = 2 * (size(knots) - order - 2)

        do i = 1, nvalences
            kapa = valence_kappas(i)
            nmna = valence_nmains(i)

            hfmat = 0.0d0
            do j = 1, ncores
                kapb = core_kappas(j)
                coefb = core_coefficients(:, j)
                call dhfcb_matrix(knots, order, kapa, kapb, coefb, nodes, weights, ahf)
                hfmat = hfmat + 2 * abs(kapb) * ahf
            end do

            call dirac(knots, order, kapa, potential, nodes, weights, spectrum, coefmat, hfmat=hfmat)

            la = (abs(2 * kapa + 1) - 1) / 2
            na = nspectrum / 2 + nmna - la
            valence_energies(i) = spectrum(na) - ialpha**2
            valence_coefficients(:, i) = coefmat(:, na)
        end do

        return
    end subroutine valence_solutions

end program dhfcb
