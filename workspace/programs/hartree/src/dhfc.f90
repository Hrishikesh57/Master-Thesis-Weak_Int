program dhfc
    use iso_fortran_env, only: error_unit
    use bspwf_constants, only: bohr_radius, ialpha
    use bspwf_dirac, only: dirac
    use bspwf_function, only: function_type
    use bspwf_hartree_fock, only: dhfc_dir_matrix, dhfc_matrix
    use bspwf_knots, only: knots_01
    use bspwf_potentials, only: fermi_potential
    use bspqed_quadrature, only: nodes20, weights20
    use bspqed_dkb, only: radialg, radialf
    use bspqed_physconsts, only: alpha, electron_reduced_compton
    use bspqed_utilities, only: kntseq1
    use bspqed_potentials, only: sphere_potential
    implicit none

    double precision, parameter :: pi = 4.0d0 * atan(1.0d0)

    double precision, dimension(:, :), allocatable :: core_coefficients, valence_coefficients
    double precision, dimension(:), allocatable :: knots, core_energies, valence_energies
    double precision :: rms, zalpha, rnucl, rcavity, epsb, ecore_total, epsv
    real :: start, finish
    integer, dimension(:), allocatable :: core_kappas, core_nmains, valence_kappas, valence_nmains
    integer :: znucl, ncores, nvalences, i, order, npoints, nknots, nbsplines, kapb, nmnb, kapv, nmnv
    type(sphere_potential) :: nuclpot

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

    ! Use same setup as practise.f90
    order = 9
    zalpha = znucl * alpha
    
    ! Define the nuclear radius for a uniform sphere
    rnucl = sqrt(5.0d0/3.0d0) * rms / electron_reduced_compton
    
    ! Choose a cavity (box) radius in natural units
    rcavity = 2.5d0 / alpha
    
    ! B-spline grid
    npoints = 1000
    nknots = npoints + 2 * (order - 1)
    nbsplines = nknots - order

    allocate (knots(nknots))
    knots = 0.0d0
    
    call kntseq1(rnucl, rcavity, npoints/10, order, knots)
    if (knots(1) <= 0.0d0) knots(1) = 1.0d-6

    ! Define the nuclear potential
    nuclpot%zalpha = zalpha
    nuclpot%rad = rnucl

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

    ecore_total = total_core_energy(knots, order, core_kappas, core_energies, core_coefficients, nodes20,&
    & weights20)

    write (*, *) "# ", repeat("-", 86)
    write (*, *) "# Total core energy: ", ecore_total
    write (*, *) "# ", repeat("-", 86)

    ! Print core wavefunctions (multi-electron)
    write (*, '(A)') "# *** ABOUT TO CALL print_wavefunctions FOR CORE ***"
    write (*, '(A,I0,A,I0)') "# Core coefficients size: ", size(core_coefficients, 1), "x", size(core_coefficients, 2)
    write (*, '(A,999I5)') "# Core kappas: ", core_kappas
    write (*, '(A,999I5)') "# Core nmains: ", core_nmains
    write (*, '(A,999ES12.5)') "# Core energies: ", core_energies
    
    ! Force output before calling subroutine
    call flush(6)
    
    call print_wavefunctions(knots, order, core_coefficients, core_kappas, core_nmains, core_energies, "core")
    
    write (*, '(A)') "# *** FINISHED CALLING print_wavefunctions FOR CORE ***"
    call flush(6)

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
    write (*, *) "# ", repeat("-", 86)

    ! Print valence wavefunctions (in multi-electron field)
    write (*, '(A)') "# *** ABOUT TO CALL print_wavefunctions FOR VALENCE ***"
    write (*, '(A,I0,A,I0)') "# Valence coefficients size: ", size(valence_coefficients, 1), "x", size(valence_coefficients, 2)
    write (*, '(A,999I5)') "# Valence kappas: ", valence_kappas
    write (*, '(A,999I5)') "# Valence nmains: ", valence_nmains
    write (*, '(A,999ES12.5)') "# Valence energies: ", valence_energies
    
    ! Force output before calling subroutine
    call flush(6)
    
    call print_wavefunctions(knots, order, valence_coefficients, valence_kappas, valence_nmains, valence_energies, "valence")
    
    write (*, '(A)') "# *** FINISHED CALLING print_wavefunctions FOR VALENCE ***"
    call flush(6)

    write (*, *) "# ", repeat("=", 86)
    write (*, *) "# PROGRAM COMPLETED SUCCESSFULLY"
    call flush(6)

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

        double precision, dimension(2*(size(knots)-order-2), 2*(size(knots)-order-2)) :: hfmat
        double precision, dimension(2*(size(knots)-order), 2*(size(knots)-order-2)) :: coefmat
        double precision, dimension(size(core_coefficients, 1), size(core_coefficients, 2)) :: prev_core_coefs
        double precision, dimension(size(core_kappas), 2) :: epsbs
        double precision, dimension(2*(size(knots)-order-2)) :: spectrum
        double precision, dimension(2*(size(knots)-order)) :: coefb
        double precision, dimension(size(core_kappas)) :: relerrs
        double precision :: tolabs, tolrel, epsb, abserr, tolerance, relerr
        integer, dimension(1) :: temp
        integer :: maxstep, ncores, nspectrum, step, i, kapb, nmnb, lb, nb, imax
        logical :: show

        tolabs = 0.0d0
        if (present(epsabs)) tolabs = epsabs

        tolrel = 5.0d-10
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

                call hartree_fock_matrix(kapb, core_kappas, prev_core_coefs, nodes, weights, hfmat)
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

            call hartree_fock_matrix(kapb, core_kappas, prev_core_coefs, nodes, weights, hfmat)
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

    double precision function total_core_energy(knots, order, core_kappas, core_energies, core_coefficients, nodes&
    &, weights)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, dimension(:), intent(in) :: core_kappas
        double precision, dimension(:), intent(in) :: core_energies
        double precision, dimension(:, :), intent(in) :: core_coefficients
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision, dimension(2*(size(knots)-order-2), 2*(size(knots)-order-2)) :: hfmat
        double precision, dimension(2*(size(knots)-order)) :: coefa
        double precision, dimension(2*(size(knots)-order-2)) :: ca
        double precision :: etot, epsa, vhf_aa
        integer :: ncores, nbsplines, nspectrum, i, kapa, k, j

        ncores = size(core_kappas)
        nbsplines = size(knots) - order
        nspectrum = 2 * (size(knots) - order - 2)
        etot = 0.0d0
        do i = 1, ncores
            kapa = core_kappas(i)
            epsa = core_energies(i)
            coefa = core_coefficients(:, i)
            ca(:(nbsplines-2)) = coefa(2:(nbsplines-1))
            ca((nbsplines-1):) = coefa((nbsplines+2):(2*nbsplines-1))
            call hartree_fock_matrix(kapa, core_kappas, core_coefficients, nodes, weights, hfmat)
            vhf_aa = 0.0d0
            do k = 1, nspectrum
                do j = 1, nspectrum
                    vhf_aa = vhf_aa + ca(j) * hfmat(j, k) * ca(k)
                end do
            enddo
            etot = etot + 2 * abs(kapa) * (epsa - 0.5d0 * vhf_aa)
        end do
        total_core_energy = etot
        return
    end function total_core_energy

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

        double precision, dimension(2*(size(knots)-order-2), 2*(size(knots)-order-2)) :: hfmat
        double precision, dimension(2*(size(knots)-order), 2*(size(knots)-order-2)) :: coefmat
        double precision, dimension(2*(size(knots)-order-2)) :: spectrum
        integer :: nvalences, nspectrum, i, kapa, nmna, la, na

        nvalences = size(valence_kappas)
        nspectrum = 2 * (size(knots) - order - 2)

        do i = 1, nvalences
            kapa = valence_kappas(i)
            nmna = valence_nmains(i)

            call hartree_fock_matrix(kapa, core_kappas, core_coefficients, nodes, weights, hfmat)
            call dirac(knots, order, kapa, potential, nodes, weights, spectrum, coefmat, hfmat=hfmat)

            la = (abs(2 * kapa + 1) - 1) / 2
            na = nspectrum / 2 + nmna - la
            valence_energies(i) = spectrum(na) - ialpha**2
            valence_coefficients(:, i) = coefmat(:, na)
        end do

        return
    end subroutine valence_solutions

    subroutine print_wavefunctions(knots, order, coefficients, kappas, nmains, energies, label)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, dimension(:), intent(in) :: kappas, nmains
        double precision, dimension(:, :), intent(in) :: coefficients
        double precision, dimension(:), intent(in) :: energies
        character(len=*), intent(in) :: label
        
        double precision :: x, gx, fx, x_start, x_end, x_step
        integer :: i, j, npoints_output, inbv, nstates, unit_num, iostat
        character(len=100) :: filename
        
        ! Use explicit format for better output control
        write(*, '(A)') "# ========================================="
        write(*, '(A)') "# ENTERING print_wavefunctions subroutine"
        write(*, '(A)') "# ========================================="
        write(*, '(A,A)') "# Label: ", trim(label)
        write(*, '(A,I0)') "# znucl = ", znucl
        call flush(6)
        
        nstates = size(kappas)
        write(*, '(A,I0)') "# Number of states to process: ", nstates
        
        if (nstates <= 0) then
            write(*, '(A)') "# ERROR: No states to process!"
            call flush(6)
            return
        endif
        
        npoints_output = 1000
        x_start = 1.0d-6
        x_end = 20.0d0
        x_step = (x_end - x_start) / (npoints_output - 1)
        
        write(*, '(A,A,A)') "# Writing ", trim(label), " wavefunctions to files..."
        write(*, '(A,I0,A)') "# Will create ", nstates, " files"
        
        ! Check if arrays are properly allocated
        if (size(coefficients, 1) == 0) then
            write(*, '(A)') "# ERROR: coefficients array is empty!"
            call flush(6)
            return
        endif
        
        write(*, '(A,I0,A,I0)') "# Coefficient matrix size: ", size(coefficients, 1), "x", size(coefficients, 2)
        call flush(6)
        
        ! Print wavefunctions for each state
        do j = 1, nstates
            write(*, '(A)') "# ==============================================="
            write(*, '(A,I0,A,I0)') "# Processing state ", j, " of ", nstates
            write(*, '(A,I0,A,I0,A,ES12.5)') "# n=", nmains(j), " kappa=", kappas(j), " E=", energies(j)
            call flush(6)
            
            ! Create filename for each state
            write(filename, '(I0,A,A,I0,A,I0,A)') znucl, '_', trim(label), '_n', nmains(j), '_kappa', kappas(j), '.txt'
            
            write(*, '(A,A,A)') "# Creating file: ", trim(filename), ""
            call flush(6)
            
            unit_num = 20 + j
            
            ! Try to open the file
            open(unit=unit_num, file=filename, status='replace', action='write', iostat=iostat)
            
            if (iostat /= 0) then
                write(*, '(A,A)') '# ERROR: Could not create file: ', trim(filename)
                write(*, '(A,I0)') '# IOSTAT = ', iostat
                call flush(6)
                cycle
            endif
            
            write(*, '(A,A)') '# File opened successfully: ', trim(filename)
            call flush(6)
            
            ! Write file header
            write(unit_num, '(A)') '# Multi-electron Dirac-Hartree-Fock wavefunction'
            write(unit_num, '(A,I0,A,I0,A,ES15.8)') '# n=', nmains(j), ' kappa=', kappas(j), ' E=', energies(j)
            write(unit_num, '(A)') '# r (natural units, 1/MeV)    g_nk(r)           f_nk(r)'
            
            ! Test the first point
            x = x_start
            inbv = 1
            
            write(*, '(A,ES12.5)') '# Testing functions at r=', x
            call flush(6)
            
            gx = radialg(knots, order, kappas(j), coefficients(:,j), x, inbv)
            fx = radialf(knots, order, kappas(j), coefficients(:,j), x, inbv)
            
            write(*, '(A,ES12.5,A,ES12.5,A,ES12.5)') '# Test: r=', x, ' g=', gx, ' f=', fx
            call flush(6)
            
            ! Write all data points
            write(*, '(A,I0,A)') '# Writing ', npoints_output, ' points to file...'
            call flush(6)
            
            do i = 1, npoints_output
                x = x_start + (i - 1) * x_step
                inbv = 1
                gx = radialg(knots, order, kappas(j), coefficients(:,j), x, inbv)
                fx = radialf(knots, order, kappas(j), coefficients(:,j), x, inbv)
                
                write(unit_num, '(F12.6, 3X, ES18.10E3, 3X, ES18.10E3)') x, gx, fx
            end do
            
            close(unit_num)
            write(*, '(A,A)') '# File completed: ', trim(filename)
            write(*, '(A,I0,A)') '# Contains ', npoints_output, ' data points'
            
            ! Print sample data to terminal
            write(*, '(A)') '# Sample data (first 3 points):'
            write(*, '(A)') '#    r            g(r)             f(r)'
            do i = 1, 3
                x = x_start + (i - 1) * x_step
                inbv = 1
                gx = radialg(knots, order, kappas(j), coefficients(:,j), x, inbv)
                fx = radialf(knots, order, kappas(j), coefficients(:,j), x, inbv)
                write(*, '(A,F8.6,2X,ES12.5E2,2X,ES12.5E2)') '#  ', x, gx, fx
            end do
            write(*, '(A)') "# -----------------------------------------------"
            call flush(6)
        end do
        
        write(*, '(A)') "# ========================================="
        write(*, '(A,A,A)') "# All ", trim(label), " files created successfully!"
        write(*, '(A)') "# Files are in current working directory"
        write(*, '(A)') "# EXITING print_wavefunctions subroutine"
        write(*, '(A)') "# ========================================="
        call flush(6)
        return
    end subroutine print_wavefunctions

    subroutine hartree_fock_matrix(kappa, core_kappas, core_coefficients, nodes, weights, hfmat)
        integer, intent(in) :: kappa
        integer, dimension(:), intent(in) :: core_kappas
        double precision, dimension(:, :), intent(in) :: core_coefficients
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision, dimension(:, :), intent(out) :: hfmat

        double precision, dimension(size(hfmat, 1), size(hfmat, 2)) :: ahf
        double precision, dimension(size(core_coefficients, 1)) :: coefb
        integer :: ncores, i, kapb

        ncores = size(core_kappas)
        hfmat = 0.0d0
        do i = 1, ncores
            kapb = core_kappas(i)
            coefb = core_coefficients(:, i)
            call dhfc_matrix(knots, order, kappa, kapb, coefb, nodes, weights, ahf)
            hfmat = hfmat + 2 * abs(kapb) * ahf
        end do
        return
    end subroutine hartree_fock_matrix

end program dhfc
