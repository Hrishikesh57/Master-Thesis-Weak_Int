program shf
    use iso_fortran_env, only: error_unit
    use bspwf_function, only: function_type
    use bspwf_hartree_fock, only: shf_dir_matrix, shf_matrix
    use bspwf_knots, only: knots_exp
    use bspwf_potentials, only: coulomb_potential
    use bspwf_quadrature, only: ngl10, nodes10, weights10
    use bspwf_schroedinger, only: schroedinger
    implicit none

    double precision, dimension(:, :), allocatable :: core_coefficients, valence_coefficients
    double precision, dimension(:), allocatable :: knots, core_energies, valence_energies
    double precision :: rcavity, epsb, ecore_total, epsv
    real :: start, finish
    integer, dimension(:), allocatable :: core_orbitals, core_nmains, valence_orbitals, valence_nmains
    integer :: znucl, ncores, nvalences, i, order, npoints, nknots, nbsplines, nspectrum, lb, nmnb, lv, nmnv
    type(coulomb_potential) :: nuclpot

    read (*, *) znucl

    read (*, *) ncores
    allocate (core_orbitals(ncores))
    core_orbitals = 0
    allocate (core_nmains(ncores))
    core_nmains = 0
    do i = 1, ncores
        read (*, *) core_nmains(i), core_orbitals(i)
    end do

    read (*, *) nvalences
    allocate (valence_orbitals(nvalences))
    valence_orbitals = 0
    allocate (valence_nmains(nvalences))
    valence_nmains = 0
    do i = 1, nvalences
        read (*, *) valence_nmains(i), valence_orbitals(i)
    end do

    write (*, *) "# ", repeat("=", 86)
    write (*, *) znucl
    do i = 1, ncores
        write (*, *) i, core_nmains(i), core_orbitals(i)
    end do
    do i = 1, nvalences
        write (*, *) i, valence_nmains(i), valence_orbitals(i)
    end do
    write (*, *) "# ", repeat("-", 86)

    order = 6
    npoints = 36
    nknots = npoints + 2 * (order - 1)
    nbsplines = nknots - order
    nspectrum = nbsplines - 2

    allocate (knots(nknots))
    knots = 0.0d0

    rcavity = 40.0d0
    call knots_exp(5.0d0, rcavity, order, knots)

    nuclpot%znucl = znucl

    allocate (core_energies(ncores))
    core_energies = 0.0d0

    allocate (core_coefficients(nbsplines, ncores))
    core_coefficients = 0.0d0

    call cpu_time(start)
    call starting_solutions(knots, order, nuclpot, core_orbitals, core_nmains, nodes10, weights10, core_energies,&
    & core_coefficients, print_output=.true., epsrel=5.0D-10)
    call cpu_time(finish)

    write (*, *) "# ", repeat("-", 86)
    write (*, *) "# Calculation of starting approximations done!"
    write (*, *) "# Total duration in seconds: ", finish - start
    write (*, *) "# Energies of the core states are ..."
    write (*, *) "# ", repeat("-", 86)
    do i = 1, ncores
        lb = core_orbitals(i)
        nmnb = core_nmains(i)
        epsb = core_energies(i)
        write (*, *) nmnb, lb, epsb
    end do
    write (*, *) "# ", repeat("-", 86)

    call cpu_time(start)
    call refining_solutions(knots, order, nuclpot, core_orbitals, core_nmains, nodes10, weights10, core_energies,&
    & core_coefficients, print_output=.true., epsrel=5.0D-14)
    call cpu_time(finish)

    write (*, *) "# ", repeat("-", 86)
    write (*, *) "# Refining the solutions done!"
    write (*, *) "# Total duration in seconds: ", finish - start
    write (*, *) "# Energies of the core states are ..."
    write (*, *) "# ", repeat("-", 86)
    do i = 1, ncores
        lb = core_orbitals(i)
        nmnb = core_nmains(i)
        epsb = core_energies(i)
        write (*, *) nmnb, lb, epsb
    end do
    write (*, *) "# ", repeat("-", 86)

    ecore_total = total_core_energy(knots, order, core_orbitals, core_energies, core_coefficients, nodes10,&
    & weights10)

    write (*, *) "# ", repeat("-", 86)
    write (*, *) "# Total core energy: ", ecore_total
    write (*, *) "# ", repeat("-", 86)

    allocate (valence_energies(nvalences))
    valence_energies = 0.0d0

    allocate (valence_coefficients(nbsplines, nvalences))
    valence_coefficients = 0.0d0

    call valence_solutions(knots, order, nuclpot, core_orbitals, core_coefficients, valence_orbitals,&
    & valence_nmains, nodes10, weights10, valence_energies, valence_coefficients)

    write (*, *) "# ", repeat("-", 86)
    write (*, *) "# Energies of the valence states are ..."
    write (*, *) "# ", repeat("-", 86)
    do i = 1, nvalences
        lv = valence_orbitals(i)
        nmnv = valence_nmains(i)
        epsv = valence_energies(i)
        write (*, *) nmnv, lv, epsv
    end do
    write (*, *) "# ", repeat("=", 86)

    return

contains

    subroutine starting_solutions(knots, order, potential, core_orbitals, core_nmains, nodes, weights,&
    & core_energies, core_coefficients, epsabs, epsrel, maxiter, print_output)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        class(function_type), intent(in) :: potential
        integer, dimension(:), intent(in) :: core_orbitals
        integer, dimension(:), intent(in) :: core_nmains
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision, dimension(:), intent(out) :: core_energies
        double precision, dimension(:, :), intent(out) :: core_coefficients

        double precision, intent(in), optional :: epsabs
        double precision, intent(in), optional :: epsrel
        integer, intent(in), optional :: maxiter
        logical, intent(in), optional :: print_output

        double precision, dimension(size(knots)-order-2, size(knots)-order-2) :: hfmat, hfmat_curr, ahf
        double precision, dimension(size(knots)-order, size(knots)-order-2) :: coefmat
        double precision, dimension(size(core_orbitals), 2) :: epsbs
        double precision, dimension(size(knots)-order-2) :: spectrum
        double precision, dimension(size(knots)-order) :: coefb
        double precision, dimension(size(core_orbitals)) :: relerrs
        double precision :: tolabs, tolrel, epsb, abserr, tolerance, relerr
        integer, dimension(1) :: temp
        integer :: maxstep, ncores, step, nelectrons, i, lb, nmnb, nb, imax
        logical :: show

        tolabs = 0.0d0
        if (present(epsabs)) tolabs = epsabs

        tolrel = 5.0d-04
        if (present(epsrel)) tolrel = epsrel

        maxstep = 50
        if (present(maxiter)) maxstep = maxiter

        show = .false.
        if (present(print_output)) show = print_output

        ncores = size(core_orbitals)

        hfmat = 0.0d0
        do step = 1, 2
            hfmat_curr = 0.0d0
            nelectrons = 0
            do i = 1, ncores
                lb = core_orbitals(i)
                nmnb = core_nmains(i)

                call schroedinger(knots, order, lb, potential, nodes, weights, spectrum, coefmat, hfmat=hfmat)

                nb = nmnb - lb
                epsb = spectrum(nb)
                coefb = coefmat(:, nb)
                core_coefficients(:, i) = coefb

                call shf_dir_matrix(knots, order, coefb, nodes, weights, ahf)
                hfmat_curr = hfmat_curr + 2 * (2 * lb + 1) * ahf
                nelectrons = nelectrons + 2 * (2 * lb + 1)

                epsbs(i, step) = epsb

                if (show) write (*, *) step, nmnb, lb, epsb
            end do
            hfmat_curr = (nelectrons - 1) * hfmat_curr / nelectrons
            hfmat = 0.5d0 * hfmat_curr + 0.5d0 * hfmat

            if (show) write (*, *) ""
        end do

        do step = 3, maxstep
            do i = 1, ncores
                relerrs(i) = abs((epsbs(i, 2) - epsbs(i, 1)) / epsbs(i, 2))
            end do
            temp = maxloc(relerrs)
            imax = temp(1)

            abserr = abs(epsbs(imax, 2) - epsbs(imax, 1))
            tolerance = max(tolabs, tolrel * abs(epsbs(imax, 2)))
            if (abserr <= tolerance) exit

            hfmat_curr = 0.0d0
            nelectrons = 0
            do i = 1, ncores
                lb = core_orbitals(i)
                nmnb = core_nmains(i)

                call schroedinger(knots, order, lb, potential, nodes, weights, spectrum, coefmat, hfmat=hfmat)

                nb = nmnb - lb
                epsb = spectrum(nb)
                coefb = coefmat(:, nb)
                core_coefficients(:, i) = coefb

                call shf_dir_matrix(knots, order, coefb, nodes, weights, ahf)
                hfmat_curr = hfmat_curr + 2 * (2 * lb + 1) * ahf
                nelectrons = nelectrons + 2 * (2 * lb + 1)

                epsbs(i, 1) = epsbs(i, 2)
                epsbs(i, 2) = epsb

                relerr = abs((epsbs(i, 1) - epsbs(i, 2)) / epsbs(i, 2))

                if (show) write (*, *) step, nmnb, lb, epsb, relerr
            end do
            hfmat_curr = (nelectrons - 1) * hfmat_curr / nelectrons
            hfmat = 0.5d0 * hfmat_curr + 0.5d0 * hfmat

            if (show) write (*, *) ""
        end do

        core_energies = epsbs(:, 2)

        return
    end subroutine starting_solutions

    subroutine refining_solutions(knots, order, potential, core_orbitals, core_nmains, nodes, weights,&
    & core_energies, core_coefficients, epsabs, epsrel, maxiter, print_output)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        class(function_type), intent(in) :: potential
        integer, dimension(:), intent(in) :: core_orbitals
        integer, dimension(:), intent(in) :: core_nmains
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision, dimension(:), intent(in out) :: core_energies
        double precision, dimension(:, :), intent(in out) :: core_coefficients

        double precision, intent(in), optional :: epsabs
        double precision, intent(in), optional :: epsrel
        integer, intent(in), optional :: maxiter
        logical, intent(in), optional :: print_output

        double precision, dimension(size(knots)-order-2, size(knots)-order-2) :: hfmat
        double precision, dimension(size(knots)-order, size(knots)-order-2) :: coefmat
        double precision, dimension(size(core_coefficients, 1), size(core_coefficients, 2)) :: prev_core_coefs
        double precision, dimension(size(core_orbitals), 2) :: epsbs
        double precision, dimension(size(knots)-order-2) :: spectrum
        double precision, dimension(size(knots)-order) :: coefb
        double precision, dimension(size(core_orbitals)) :: relerrs
        double precision :: tolabs, tolrel, epsb, abserr, tolerance, relerr
        integer, dimension(1) :: temp
        integer :: maxstep, ncores, step, i, lb, nmnb, nb, imax
        logical :: show

        tolabs = 0.0d0
        if (present(epsabs)) tolabs = epsabs

        tolrel = 5.0d-10
        if (present(epsrel)) tolrel = epsrel

        maxstep = 100
        if (present(maxiter)) maxstep = maxiter

        show = .false.
        if (present(print_output)) show = print_output

        ncores = size(core_orbitals)

        do step = 1, 2
            prev_core_coefs = core_coefficients
            do i = 1, ncores
                lb = core_orbitals(i)
                nmnb = core_nmains(i)

                call hartree_fock_matrix(lb, core_orbitals, prev_core_coefs, nodes, weights, hfmat)
                call schroedinger(knots, order, lb, potential, nodes, weights, spectrum, coefmat, hfmat=hfmat)

                nb = nmnb - lb
                epsb = spectrum(nb)
                coefb = coefmat(:, nb)
                core_coefficients(:, i) = coefb

                epsbs(i, step) = epsb

                if (show) write (*, *) step, nmnb, lb, epsb
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

            lb = core_orbitals(imax)
            nmnb = core_nmains(imax)

            call hartree_fock_matrix(lb, core_orbitals, prev_core_coefs, nodes, weights, hfmat)
            call schroedinger(knots, order, lb, potential, nodes, weights, spectrum, coefmat, hfmat=hfmat)

            nb = nmnb - lb
            epsb = spectrum(nb)
            coefb = coefmat(:, nb)
            core_coefficients(:, imax) = coefb

            epsbs(imax, 1) = epsbs(imax, 2)
            epsbs(imax, 2) = epsb

            relerr = abs((epsbs(imax, 1) - epsbs(imax, 2)) / epsbs(imax, 2))

            if (show) write (*, *) step, nmnb, lb, epsb, relerr
        end do

        core_energies = epsbs(:, 2)

        return
    end subroutine refining_solutions

    double precision function total_core_energy(knots, order, core_orbitals, core_energies, core_coefficients,&
    & nodes, weights)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, dimension(:), intent(in) :: core_orbitals
        double precision, dimension(:), intent(in) :: core_energies
        double precision, dimension(:, :), intent(in) :: core_coefficients
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision, dimension(size(knots)-order-2, size(knots)-order-2) :: hfmat
        double precision, dimension(size(knots)-order) :: coefa
        double precision :: etot, epsa, vhf_aa
        integer :: ncores, nspectrum, i, la, k, j

        ncores = size(core_orbitals)
        nspectrum = size(knots) - order - 2
        etot = 0.0d0
        do i = 1, ncores
            la = core_orbitals(i)
            epsa = core_energies(i)
            coefa = core_coefficients(:, i)
            call hartree_fock_matrix(la, core_orbitals, core_coefficients, nodes, weights, hfmat)
            vhf_aa = 0.0d0
            do k = 1, nspectrum
                do j = 1, nspectrum
                    vhf_aa = vhf_aa + coefa(j+1) * hfmat(j, k) * coefa(k+1)
                end do
            enddo
            etot = etot + 2 * (2 * la + 1) * (epsa - 0.5d0 * vhf_aa)
        end do
        total_core_energy = etot
        return
    end function total_core_energy

    subroutine valence_solutions(knots, order, potential, core_orbitals, core_coefficients, valence_orbitals,&
    & valence_nmains, nodes, weights, valence_energies, valence_coefficients)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        class(function_type), intent(in) :: potential
        integer, dimension(:), intent(in) :: core_orbitals
        double precision, dimension(:, :), intent(in) :: core_coefficients
        integer, dimension(:), intent(in) :: valence_orbitals
        integer, dimension(:), intent(in) :: valence_nmains
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision, dimension(:), intent(out) :: valence_energies
        double precision, dimension(:, :), intent(out) :: valence_coefficients

        double precision, dimension(size(knots)-order-2, size(knots)-order-2) :: hfmat
        double precision, dimension(size(knots)-order, size(knots)-order-2) :: coefmat
        double precision, dimension(size(knots)-order-2) :: spectrum
        integer :: ncores, nvalences, i, la, nmna, na

        ncores = size(core_orbitals)
        nvalences = size(valence_orbitals)

        do i = 1, nvalences
            la = valence_orbitals(i)
            nmna = valence_nmains(i)

            call hartree_fock_matrix(la, core_orbitals, core_coefficients, nodes, weights, hfmat)
            call schroedinger(knots, order, la, potential, nodes, weights, spectrum, coefmat, hfmat=hfmat)

            na = nmna - la
            valence_energies(i) = spectrum(na)
            valence_coefficients(:, i) = coefmat(:, na)
        end do

        return
    end subroutine valence_solutions

    subroutine hartree_fock_matrix(l, core_orbitals, core_coefficients, nodes, weights, hfmat)
        integer, intent(in) :: l
        integer, dimension(:), intent(in) :: core_orbitals
        double precision, dimension(:, :), intent(in) :: core_coefficients
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision, dimension(:, :), intent(out) :: hfmat

        double precision, dimension(size(hfmat, 1), size(hfmat, 2)) :: ahf
        double precision, dimension(size(core_coefficients, 1)) :: coefb
        integer :: ncores, i, lb

        ncores = size(core_orbitals)
        hfmat = 0.0d0
        do i = 1, ncores
            lb = core_orbitals(i)
            coefb = core_coefficients(:, i)
            call shf_matrix(knots, order, l, lb, coefb, nodes, weights, ahf)
            hfmat = hfmat + 2 * (2 * lb + 1) * ahf
        end do
        return
    end subroutine hartree_fock_matrix

end program shf
