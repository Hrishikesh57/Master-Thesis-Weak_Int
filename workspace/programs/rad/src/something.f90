program compute_radial_integrals
    use bspqed_potentials
    use bspqed_physconsts
    use bspqed_quadrature       ! <--- We bring in ngl40, nodes40, weights40, etc.
    implicit none

    double precision, parameter :: e_charge = 1.0d0
    character(len=*), parameter :: rad_base = '.'
    character(len=20) :: main_choice

    call prompt_main_choice(main_choice)

    select case (trim(main_choice))
    case ('energy')
        call run_energy_mode()
    case ('amplitudes')
        call run_amplitude_mode()
    case default
        print *, 'Invalid main choice. Please enter energy or amplitudes. Exiting.'
        stop
    end select

contains
    subroutine prompt_main_choice(main_choice)
        implicit none
        character(len=*), intent(out) :: main_choice

        print *, 'Do you want to calculate energy corrections or amplitudes? (enter energy or amplitudes)'
        read *, main_choice
    end subroutine prompt_main_choice

    subroutine run_energy_mode()
        implicit none
        character(len=20) :: order_choice
        character(len=20) :: state_initial, state_intermediate
        character(len=256) :: filename_n, filename_s, reference_filename
        integer :: npoints, npoints_n, npoints_s
        double precision :: z, n, r_fm, r_rel
        double precision :: e_pot, m_n_s, delta_e
        double precision :: e_initial, e_intermediate
        double precision, allocatable :: r(:), g(:), f(:)
        double precision, allocatable :: r_n(:), g_n(:), f_n(:)
        double precision, allocatable :: r_s(:), g_s(:), f_s(:)
        type(weak_potential) :: potential

        print *, 'Do you want 1st-order or 2nd-order energy corrections? (enter 1 or 2)'
        read *, order_choice

        call prompt_nuclear_inputs(z, n, r_fm)
        call initialize_weak_potential(z, n, r_fm, potential, r_rel)
        call print_potential_summary(potential, r_rel)

        reference_filename = resolve_reference_radial_filename()

        if (trim(order_choice) == '1') then
            call read_radial_functions(reference_filename, npoints, r, g, f)
            e_pot = integrate_single_channel(0d0, r_rel, npoints, r, g, f, potential, .true.)
            print *, 'E_pot (0 to R) using Gauss–Legendre = ', e_pot, ' in relativistic units (m_e = 1)'
            deallocate(r, g, f)
        else if (trim(order_choice) == '2') then
            print *, 'Enter the state for the initial energy (e.g., 1s, 2s):'
            read *, state_initial
            print *, 'Enter the state for the intermediate energy (e.g., 2p, 3p, 4p, ..., 10p):'
            read *, state_intermediate

            call prompt_state_energy('initial', state_initial, e_initial)
            call prompt_state_energy('intermediate', state_intermediate, e_intermediate)

            call load_state_wavefunction(nint(z), trim(state_intermediate), .false., filename_n, npoints_n, r_n, g_n, f_n)
            call load_state_wavefunction(nint(z), trim(state_initial), .false., filename_s, npoints_s, r_s, g_s, f_s)

            print *, 'Reading wave functions from:'
            print *, 'Intermediate state (', trim(state_intermediate), '): ', trim(filename_n)
            print *, 'Initial state (', trim(state_initial), '): ', trim(filename_s)

            m_n_s = -integrate_two_channel(0d0, r_rel, npoints_n, r_n, g_n, f_n, &
                                           npoints_s, r_s, g_s, f_s, potential, .true.)

            print *, 'M_n,s (0 to R) using Gauss–Legendre = ', m_n_s

            delta_e = (m_n_s**2) / (e_initial - e_intermediate)
            print *, 'For initial state ', trim(state_initial), ' (E = ', e_initial, &
                      ') and intermediate state ', trim(state_intermediate), ' (E = ', e_intermediate, '):'
            print *, 'Second-order energy correction ΔE^(2) = ', delta_e

            call update_second_order_table(trim(state_initial), trim(state_intermediate), delta_e)
            call deallocate_wavefunction(r_n, g_n, f_n)
            call deallocate_wavefunction(r_s, g_s, f_s)
        else
            print *, 'Invalid order choice. Exiting.'
            stop
        end if
    end subroutine run_energy_mode

    subroutine run_amplitude_mode()
        implicit none
        integer :: num_intermediate_states, i_state
        logical :: use_scf
        character(len=20) :: state_initial, state_final, system_label
        character(len=20), allocatable :: intermediate_states(:)
        character(len=256) :: filename_s, filename_f
        integer :: npoints_s, npoints_f
        double precision :: z, n, r_fm, r_rel
        double precision :: e_initial_stark, e_final_stark, e0_field
        double precision :: a_s, a_w, eta_asymmetry
        double precision, allocatable :: e_intermediate_array(:)
        double precision, allocatable :: r_s(:), g_s(:), f_s(:)
        double precision, allocatable :: r_f(:), g_f(:), f_f(:)
        type(weak_potential) :: potential

        call prompt_nuclear_inputs(z, n, r_fm)
        call initialize_weak_potential(z, n, r_fm, potential, r_rel)

        print *, 'Calculating both stark and weak amplitudes...'
        call print_potential_summary(potential, r_rel)

        call prompt_system_type(use_scf, system_label)

        print *, 'Enter the initial state (e.g., 1s, 2s):'
        read *, state_initial
        print *, 'Enter the final state (e.g., 1s, 2s):'
        read *, state_final

        print *, 'How many intermediate states do you want to sum over?'
        read *, num_intermediate_states

        allocate(intermediate_states(num_intermediate_states))
        allocate(e_intermediate_array(num_intermediate_states))

        print *, 'Enter the intermediate states to sum over:'
        print *, 'For p states, you can specify both j=1/2 and j=3/2 components:'
        print *, 'Examples: 2p, 2p3 (for 2p3/2), 3p, 3p3, etc.'
        do i_state = 1, num_intermediate_states
            print *, 'Enter intermediate state ', i_state, ':'
            read *, intermediate_states(i_state)
        end do

        call load_state_wavefunction(nint(z), trim(state_initial), use_scf, filename_s, npoints_s, r_s, g_s, f_s)
        call load_state_wavefunction(nint(z), trim(state_final), use_scf, filename_f, npoints_f, r_f, g_f, f_f)

        print *, 'Reading wave functions from:'
        print *, 'Initial state (', trim(state_initial), '): ', trim(filename_s)
        print *, 'Final state (', trim(state_final), '): ', trim(filename_f)

        call prompt_state_energy('initial', state_initial, e_initial_stark)
        call prompt_state_energy('final', state_final, e_final_stark)
        print *, 'Enter the electric field strength E0 (in relativistic units, m_e = 1):'
        read *, e0_field

        call prompt_intermediate_energies(intermediate_states, num_intermediate_states, e_intermediate_array)

        call accumulate_amplitudes(r_rel, potential, e_initial_stark, e_final_stark, &
                                   nint(z), use_scf, intermediate_states, e_intermediate_array, num_intermediate_states, &
                                   npoints_s, r_s, g_s, f_s, npoints_f, r_f, g_f, f_f, a_s, a_w)

        a_s = -e_charge**2 * e0_field * (1.0d0 / 27.0d0) * a_s
        a_w = a_w / (3.0d0 * sqrt(3.0d0))
        eta_asymmetry = 2.0d0 * (a_s * a_w) / (a_s * a_s)

        print *, ''
        print *, '========================================='
        print *, 'AMPLITUDE CALCULATION RESULTS'
        print *, '========================================='
        print *, 'Transition: ', trim(state_initial), ' -> ', trim(state_final)
        print *, 'System type: ', trim(system_label)
        print *, 'Electric field E0: ', e0_field
        print *, ''
        print *, 'Stark Amplitude A_S = ', a_s
        print *, 'Weak Amplitude A_W = ', a_w
        print *, 'Asymmetry Factor η = 2*(A_S*A_W)/|A_S|² = ', eta_asymmetry
        print *, '========================================='

        call write_amplitude_results_to_file(z, n, r_fm, trim(state_initial), trim(state_final), &
                                             e_initial_stark, e_final_stark, e0_field, &
                                             intermediate_states, e_intermediate_array, &
                                             num_intermediate_states, a_s, a_w, eta_asymmetry, &
                                             trim(system_label))

        call deallocate_wavefunction(r_s, g_s, f_s)
        call deallocate_wavefunction(r_f, g_f, f_f)
        deallocate(intermediate_states, e_intermediate_array)
    end subroutine run_amplitude_mode

    subroutine prompt_nuclear_inputs(z, n, r_fm)
        implicit none
        double precision, intent(out) :: z, n, r_fm

        print *, 'Enter atomic number Z (e.g., 92 for Uranium):'
        read *, z
        print *, 'Enter number of neutrons N (e.g., 146 for Uranium):'
        read *, n
        print *, 'Enter nuclear radius R in fm (e.g., 5.8571 for Uranium):'
        read *, r_fm
    end subroutine prompt_nuclear_inputs

    subroutine initialize_weak_potential(z, n, r_fm, potential, r_rel)
        implicit none
        double precision, intent(in) :: z, n, r_fm
        double precision, intent(out) :: r_rel
        type(weak_potential), intent(out) :: potential

        r_rel = convert_R_fm_to_relativistic(r_fm)
        call potential%init(z, n, r_rel, 1.0d0)
    end subroutine initialize_weak_potential

    subroutine print_potential_summary(potential, r_rel)
        implicit none
        type(weak_potential), intent(in) :: potential
        double precision, intent(in) :: r_rel

        print *, 'rho0 is :', potential%rho0
        print *, 'Nuclear radius R in relativistic units (m_e = 1) :', r_rel
        print *, 'Gfermi rel is : ', Gfermi_relativistic
        print *, 'Cv_pZ_minus_N is : ', potential%Cv_pZ_minus_N
    end subroutine print_potential_summary

    subroutine prompt_system_type(use_scf, system_label)
        implicit none
        logical, intent(out) :: use_scf
        character(len=*), intent(out) :: system_label
        integer :: system_type

        print *, 'Choose the type of wavefunction:'
        print *, '1 - Single-electron system'
        print *, '2 - Multi-electron system (SCF)'
        print *, 'Enter your choice (1 or 2):'
        read *, system_type

        select case (system_type)
        case (1)
            use_scf = .false.
            system_label = 'single'
            print *, 'Using single-electron wavefunctions'
        case (2)
            use_scf = .true.
            system_label = 'scf'
            print *, 'Using multi-electron SCF wavefunctions'
        case default
            use_scf = .false.
            system_label = 'single'
            print *, 'Invalid choice. Defaulting to single-electron system.'
        end select
    end subroutine prompt_system_type

    subroutine load_state_wavefunction(z_int, state_name, use_scf, filename, npoints, r, g, f)
        implicit none
        integer, intent(in) :: z_int
        character(len=*), intent(in) :: state_name
        logical, intent(in) :: use_scf
        character(len=*), intent(out) :: filename
        integer, intent(out) :: npoints
        double precision, allocatable, intent(out) :: r(:), g(:), f(:)

        call build_wavefunction_filename(z_int, state_name, use_scf, filename)
        call read_radial_functions(filename, npoints, r, g, f)
    end subroutine load_state_wavefunction

    subroutine deallocate_wavefunction(r, g, f)
        implicit none
        double precision, allocatable, intent(inout) :: r(:), g(:), f(:)

        if (allocated(r)) deallocate(r)
        if (allocated(g)) deallocate(g)
        if (allocated(f)) deallocate(f)
    end subroutine deallocate_wavefunction

    subroutine prompt_state_energy(role_label, state_name, energy_value)
        implicit none
        character(len=*), intent(in) :: role_label, state_name
        double precision, intent(out) :: energy_value

        print *, 'Enter the energy of the ', trim(role_label), ' state (', trim(state_name), ') in relativistic units (m_e = 1):'
        read *, energy_value
    end subroutine prompt_state_energy

    subroutine prompt_intermediate_energies(intermediate_states, num_intermediate_states, e_intermediate_array)
        implicit none
        integer, intent(in) :: num_intermediate_states
        character(len=20), dimension(num_intermediate_states), intent(in) :: intermediate_states
        double precision, dimension(num_intermediate_states), intent(out) :: e_intermediate_array

        integer :: i_state

        do i_state = 1, num_intermediate_states
            print *, 'Enter the energy of intermediate state ', trim(intermediate_states(i_state)), ' in relativistic units (m_e = 1):'
            read *, e_intermediate_array(i_state)
        end do
    end subroutine prompt_intermediate_energies

subroutine accumulate_amplitudes(r_rel, potential, e_initial, e_final, z_int, use_scf, &
                                 intermediate_states, e_intermediate_array, num_intermediate_states, &
                                 npoints_s, r_s, g_s, f_s, npoints_f, r_f, g_f, f_f, a_s, a_w)
    implicit none
    integer, intent(in) :: z_int, num_intermediate_states, npoints_s, npoints_f
    logical, intent(in) :: use_scf
    character(len=20), dimension(num_intermediate_states), intent(in) :: intermediate_states
    double precision, intent(in) :: r_rel, e_initial, e_final
    double precision, dimension(num_intermediate_states), intent(in) :: e_intermediate_array
    double precision, dimension(npoints_s), intent(in) :: r_s, g_s, f_s
    double precision, dimension(npoints_f), intent(in) :: r_f, g_f, f_f
    type(weak_potential), intent(in) :: potential
    double precision, intent(out) :: a_s, a_w

    character(len=256) :: filename_n
    character(len=20)  :: k_state
    integer :: i_state, npoints_n
    double precision :: r_temp_1, r_temp_2
    double precision :: r_weak_temp_2, r_weak_temp_3
    double precision :: weak_contribution_1, weak_contribution_2
    double precision :: state_contribution, denominator_sum, contribution_1, state_factor
    double precision, allocatable :: r_n(:), g_n(:), f_n(:)

    a_s = 0.0d0
    a_w = 0.0d0

    print *, ''
    print *, '========================================='
    print *, 'PROCESSING INTERMEDIATE STATES'
    print *, '========================================='

    do i_state = 1, num_intermediate_states
        k_state = trim(intermediate_states(i_state))
        print *, 'Processing intermediate state: ', trim(k_state)

        call load_state_wavefunction(z_int, k_state, use_scf, filename_n, npoints_n, r_n, g_n, f_n)
        print *, 'Reading intermediate state from: ', trim(filename_n)

        !----------------------------------------------------------
        ! Dipole radial integrals R_fk and R_ki
        !----------------------------------------------------------
        r_temp_1 = integrate_two_channel(0d0, 1000.0d0, npoints_f, r_f, g_f, f_f, &
                                         npoints_n, r_n, g_n, f_n, potential, .false.)

        r_temp_2 = integrate_two_channel(0d0, 1000.0d0, npoints_n, r_n, g_n, f_n, &
                                         npoints_s, r_s, g_s, f_s, potential, .false.)

        !----------------------------------------------------------
        ! Stark amplitude: keep all dipole-allowed p_j states
        !----------------------------------------------------------
        denominator_sum = 1.0d0 / (e_initial - e_intermediate_array(i_state)) + &
                          1.0d0 / (e_final   - e_intermediate_array(i_state))

        state_factor = j_state_factor(k_state)

        contribution_1 = state_factor**2 * (r_temp_1 * r_temp_2) * denominator_sum
        a_s = a_s + contribution_1

        ! Only p_{1/2} states contribute to the spin-independent weak sum.
        if (is_p12_state(k_state)) then

            r_weak_temp_2 = integrate_two_channel(0d0, r_rel, npoints_s, r_s, g_s, f_s, &
                                                  npoints_n, r_n, g_n, f_n, potential, .true.)

            r_weak_temp_3 = integrate_two_channel(0d0, r_rel, npoints_n, r_n, g_n, f_n, &
                                                  npoints_f, r_f, g_f, f_f, potential, .true.)

            weak_contribution_1 = r_temp_1 * r_weak_temp_2 / (e_initial - e_intermediate_array(i_state))
            weak_contribution_2 = r_weak_temp_3 * r_temp_2 / (e_final   - e_intermediate_array(i_state))

            state_contribution = state_factor * (weak_contribution_1 + weak_contribution_2)
            a_w = a_w + state_contribution

        end if

        call deallocate_wavefunction(r_n, g_n, f_n)
    end do

end subroutine accumulate_amplitudes

    logical function is_p12_state(state_name)
        implicit none
        character(len=*), intent(in) :: state_name

        character(len=64) :: normalized
        character(len=64) :: suffix
        integer :: i, j, code
        integer :: ppos

        normalized = ''
        j = 0
        do i = 1, len_trim(state_name)
            code = iachar(state_name(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) code = code + 32
            if ((code >= iachar('a') .and. code <= iachar('z')) .or. &
                (code >= iachar('0') .and. code <= iachar('9'))) then
                j = j + 1
                if (j <= len(normalized)) normalized(j:j) = achar(code)
            end if
        end do

        ppos = index(normalized, 'p')
        if (ppos == 0) then
            is_p12_state = .false.
            return
        end if

        ! Accept both explicit p1/2 labels and the existing bare `np` shorthand.
        suffix = adjustl(normalized(ppos + 1:))
        select case (trim(suffix))
        case ('', '1', '12')
            is_p12_state = .true.
        case default
            is_p12_state = .false.
        end select
    end function is_p12_state

    double precision function j_state_factor(state_name)
        implicit none
        character(len=*), intent(in) :: state_name

        j_state_factor = 1.0d0
        if (index(trim(state_name), '3') > 0) j_state_factor = sqrt(2.0d0)
    end function j_state_factor

    !> Convert `R` from fm to relativistic units with `hbar = c = m_e = 1`.
    double precision function convert_R_fm_to_relativistic(r_fm)
        implicit none
        double precision, intent(in) :: r_fm

        convert_R_fm_to_relativistic = r_fm / electron_reduced_compton
    end function convert_R_fm_to_relativistic

    !----------------------------------------------------------------
    ! 2) Read radial functions from a file: r, G(r), F(r)
    !----------------------------------------------------------------
    subroutine read_radial_functions(filename, npoints, r, g, f)
        implicit none
        character(len=*), intent(in) :: filename
        integer, intent(out) :: npoints
        double precision, allocatable, intent(out) :: r(:), g(:), f(:)

        integer :: unit_num, ios, i
        character(len=256) :: line, line_trimmed
        double precision :: r_tmp, g_tmp, f_tmp
        logical :: convert_from_au

        convert_from_au = .false.
        open(newunit=unit_num, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, 'Error opening file: ', filename
            stop
        end if

        ! Count valid numeric rows and ignore comments/metadata.
        npoints = 0
        do
            read(unit_num, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line_trimmed = adjustl(line)
            if (len_trim(line_trimmed) == 0) cycle
            if (index(line_trimmed, 'a.u.') > 0 .or. index(line_trimmed, 'atomic units') > 0) convert_from_au = .true.
            if (line_trimmed(1:1) == '#') cycle
            read(line_trimmed, *, iostat=ios) r_tmp, g_tmp, f_tmp
            if (ios /= 0) then
                ios = 0
                cycle
            end if
            npoints = npoints + 1
        end do

        if (npoints <= 0) then
            print *, 'No radial data rows found in file: ', filename
            stop
        end if

        rewind(unit_num)

        allocate(r(npoints), g(npoints), f(npoints))

        i = 0
        do
            read(unit_num, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line_trimmed = adjustl(line)
            if (len_trim(line_trimmed) == 0) cycle
            if (line_trimmed(1:1) == '#') cycle
            read(line_trimmed, *, iostat=ios) r_tmp, g_tmp, f_tmp
            if (ios /= 0) then
                ios = 0
                cycle
            end if
            i = i + 1
            if (i > npoints) exit
            if (convert_from_au) then
                r(i) = r_tmp / alpha
                g(i) = sqrt(alpha) * g_tmp
                f(i) = sqrt(alpha) * f_tmp
            else
                r(i) = r_tmp
                g(i) = g_tmp
                f(i) = f_tmp
            end if
        end do

        if (i /= npoints) then
            print *, 'Mismatch while reading radial data from: ', filename
            stop
        end if

        close(unit_num)
    end subroutine read_radial_functions

    function resolve_reference_radial_filename() result(filename)
        implicit none
        character(len=256) :: filename

        filename = trim(rad_base) // '/data/reference/radial_functions/radial_function_relunits.txt'
        call try_fallback_candidate(trim(rad_base) // '/radial_function_relunits.txt', filename)
        call try_fallback_candidate(trim(rad_base) // '/data/reference/radial_functions/radial_function.dat', filename)
        call try_fallback_candidate(trim(rad_base) // '/radial_function.dat', filename)
    end function resolve_reference_radial_filename

    subroutine build_wavefunction_filename(Z_int, state_name, use_scf, out_filename)
        implicit none
        integer, intent(in) :: Z_int
        character(len=*), intent(in) :: state_name
        logical, intent(in) :: use_scf
        character(len=*), intent(out) :: out_filename

        character(len=10) :: z_local
        character(len=256) :: candidate
        logical :: exists

        exists = .false.
        write(z_local, '(I0)') Z_int

        if (use_scf) then
            candidate = trim(rad_base) // '/outputs/wavefunctions/scf/' // trim(z_local) // '_scf_polynomial_table_' // trim(state_name) // '.txt'
            out_filename = candidate
            call try_candidate(candidate, exists, out_filename)
            call try_candidate(trim(rad_base) // '/data/reference/wavefunctions/scf/' // trim(z_local) // '_scf_polynomial_table_' // trim(state_name) // '.txt', exists, out_filename)
            call try_candidate(trim(rad_base) // '/data/reference/wavefunctions/scf_legacy/' // trim(z_local) // '_scf_polynomial_table_' // trim(state_name) // '.txt', exists, out_filename)
            call try_candidate(trim(rad_base) // '/multi electron /' // trim(z_local) // '_scf_polynomial_table_' // trim(state_name) // '.txt', exists, out_filename)
        else
            candidate = trim(rad_base) // '/outputs/wavefunctions/single_electron/' // trim(z_local) // '_polynomial_table_' // trim(state_name) // '.txt'
            out_filename = candidate
            call try_candidate(candidate, exists, out_filename)
            call try_candidate(trim(rad_base) // '/data/reference/wavefunctions/single_electron/' // trim(z_local) // '_polynomial_table_' // trim(state_name) // '.txt', exists, out_filename)
            call try_candidate(trim(rad_base) // '/data/reference/wavefunctions/single_electron_legacy/' // trim(z_local) // '_polynomial_table_' // trim(state_name) // '.txt', exists, out_filename)
            call try_candidate(trim(rad_base) // '/single electron wavefunctions/' // trim(z_local) // '_polynomial_table_' // trim(state_name) // '.txt', exists, out_filename)
        end if

        call try_candidate(trim(rad_base) // '/data/grasp/rel/' // trim(z_local) // '_' // trim(state_name) // '_rel.txt', exists, out_filename)
        call try_candidate(trim(rad_base) // '/GRASP_rel/' // trim(z_local) // '_' // trim(state_name) // '_rel.txt', exists, out_filename)
    end subroutine build_wavefunction_filename

    subroutine try_candidate(candidate, exists, selected)
        implicit none
        character(len=*), intent(in) :: candidate
        logical, intent(inout) :: exists
        character(len=*), intent(inout) :: selected
        logical :: candidate_exists

        if (exists) return
        inquire(file=trim(candidate), exist=candidate_exists)
        if (candidate_exists) then
            exists = .true.
            selected = candidate
        end if
    end subroutine try_candidate

    subroutine try_fallback_candidate(candidate, selected)
        implicit none
        character(len=*), intent(in) :: candidate
        character(len=*), intent(inout) :: selected
        logical :: exists

        inquire(file=trim(selected), exist=exists)
        if (.not. exists) then
            inquire(file=trim(candidate), exist=exists)
            if (exists) selected = candidate
        end if
    end subroutine try_fallback_candidate

    subroutine interpolate_components(rval, npoints, r_array, g_array, f_array, g_val, f_val)
        implicit none
        integer, intent(in) :: npoints
        double precision, intent(in) :: rval
        double precision, dimension(npoints), intent(in) :: r_array, g_array, f_array
        double precision, intent(out) :: g_val, f_val

        integer :: i_left, i_right, i_mid
        double precision :: dr, fraction

        if (rval <= 0.0d0) then
            g_val = 0.0d0
            f_val = 0.0d0
            return
        end if

        if (rval < r_array(1)) then
            if (r_array(1) > 0.0d0) then
                fraction = rval / r_array(1)
            else
                fraction = 0.0d0
            end if
            g_val = fraction * g_array(1)
            f_val = fraction * f_array(1)
            return
        end if

        if (rval > r_array(npoints)) then
            g_val = 0.0d0
            f_val = 0.0d0
            return
        end if

        i_left = 1
        i_right = npoints
        do while (i_right - i_left > 1)
            i_mid = (i_left + i_right) / 2
            if (r_array(i_mid) <= rval) then
                i_left = i_mid
            else
                i_right = i_mid
            end if
        end do

        dr = r_array(i_right) - r_array(i_left)
        if (dr == 0d0) then
            fraction = 0d0
        else
            fraction = (rval - r_array(i_left)) / dr
        end if

        g_val = g_array(i_left) + fraction * (g_array(i_right) - g_array(i_left))
        f_val = f_array(i_left) + fraction * (f_array(i_right) - f_array(i_left))
    end subroutine interpolate_components

    !> Integrate one radial channel against the weak potential.
    function integrate_single_channel(a, b, npoints, r_array, g_array, f_array, potential, use_density_norm) result(ival)
        implicit none
        double precision, intent(in) :: a, b
        integer, intent(in) :: npoints
        double precision, dimension(npoints), intent(in) :: r_array, g_array, f_array
        type(weak_potential), intent(in) :: potential
        logical, intent(in) :: use_density_norm
        double precision :: ival

        integer :: k
        double precision :: sum_val, xm, xr, r_node, upper_bound

        upper_bound = min(b, r_array(npoints))
        if (upper_bound <= a) then
            ival = 0.0d0
            return
        end if

        xm = 0.5d0 * (upper_bound + a)
        xr = 0.5d0 * (upper_bound - a)
        sum_val = 0.0d0

        do k = 1, ngl40
            r_node = xm + xr * nodes40(k)
            sum_val = sum_val + weights40(k) * single_channel_integrand(r_node, npoints, r_array, g_array, f_array, &
                                                                        potential, use_density_norm)
        end do

        ival = xr * sum_val
    end function integrate_single_channel

    function single_channel_integrand(rval, npoints, r_array, g_array, f_array, potential, use_density_norm) result(val)
        implicit none
        integer, intent(in) :: npoints
        double precision, intent(in) :: rval
        double precision, dimension(npoints), intent(in) :: r_array, g_array, f_array
        type(weak_potential), intent(in) :: potential
        logical, intent(in) :: use_density_norm
        double precision :: val, g_val, f_val

        call interpolate_components(rval, npoints, r_array, g_array, f_array, g_val, f_val)
        if (use_density_norm) then
            val = (g_val * g_val + f_val * f_val) * potential%eval(rval)
        else
            val = (g_val * f_val) * potential%eval(rval)
        end if
    end function single_channel_integrand

    !> Integrate a two-channel matrix element for either weak or dipole form.
    function integrate_two_channel(a, b, npoints_1, r_1, g_1, f_1, npoints_2, r_2, g_2, f_2, &
                                   potential, use_weak_potential) result(r_val)
        implicit none
        double precision, intent(in) :: a, b
        integer, intent(in) :: npoints_1, npoints_2
        double precision, dimension(npoints_1), intent(in) :: r_1, g_1, f_1
        double precision, dimension(npoints_2), intent(in) :: r_2, g_2, f_2
        type(weak_potential), intent(in) :: potential
        logical, intent(in) :: use_weak_potential
        double precision :: r_val

        integer :: k
        double precision :: sum_val, xm, xr, r_node, upper_bound

        upper_bound = min(b, r_1(npoints_1), r_2(npoints_2))
        if (upper_bound <= a) then
            r_val = 0.0d0
            return
        end if

        xm = 0.5d0 * (upper_bound + a)
        xr = 0.5d0 * (upper_bound - a)
        sum_val = 0.0d0

        do k = 1, ngl40
            r_node = xm + xr * nodes40(k)
            sum_val = sum_val + weights40(k) * two_channel_integrand(r_node, npoints_1, r_1, g_1, f_1, &
                                                                     npoints_2, r_2, g_2, f_2, &
                                                                     potential, use_weak_potential)
        end do

        r_val = xr * sum_val
    end function integrate_two_channel

    function two_channel_integrand(rval, npoints_1, r_1, g_1, f_1, npoints_2, r_2, g_2, f_2, &
                                   potential, use_weak_potential) result(val)
        implicit none
        integer, intent(in) :: npoints_1, npoints_2
        double precision, intent(in) :: rval
        double precision, dimension(npoints_1), intent(in) :: r_1, g_1, f_1
        double precision, dimension(npoints_2), intent(in) :: r_2, g_2, f_2
        type(weak_potential), intent(in) :: potential
        logical, intent(in) :: use_weak_potential
        double precision :: val
        double precision :: g1_val, f1_val, g2_val, f2_val

        call interpolate_components(rval, npoints_1, r_1, g_1, f_1, g1_val, f1_val)
        call interpolate_components(rval, npoints_2, r_2, g_2, f_2, g2_val, f2_val)

        if (use_weak_potential) then
            val = potential%eval(rval) * (g1_val * f2_val + f1_val * g2_val)
        else
            val = rval * (g1_val * g2_val + f1_val * f2_val)
        end if
    end function two_channel_integrand

    !----------------------------------------------------------------
    ! Update the second order energy correction table
    !----------------------------------------------------------------
    subroutine update_second_order_table(state_init, state_int, energy_corr)
        implicit none
        character(len=*), intent(in) :: state_init, state_int
        double precision, intent(in) :: energy_corr
        character(len=256) :: filename
        character(len=20) :: file_init, file_int
        double precision :: file_corr
        logical :: found
        integer :: unit, ios

        call ensure_directory(trim(rad_base) // '/outputs/results/energy/hydrogen_like')
        filename = trim(rad_base) // '/outputs/results/energy/hydrogen_like/second_order_energy_corrections.txt'
        found = .false.

        inquire(file=filename, exist=found)
        if (found) then
            open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
            if (ios /= 0) then
                print *, 'Error opening table file for reading.'
                return
            end if
            do
                read(unit, *, iostat=ios) file_init, file_int, file_corr
                if (ios /= 0) exit
                if (trim(file_init) == trim(state_init) .and. trim(file_int) == trim(state_int)) then
                    found = .true.
                    exit
                end if
            end do
            close(unit)
        end if

        if (.not. found) then
            open(newunit=unit, file=filename, status='unknown', action='write', position='append', iostat=ios)
            if (ios /= 0) then
                print *, 'Error opening table file for writing.'
                return
            end if
            write(unit, '(A,1x,A,1x,ES24.16E3)') trim(state_init), trim(state_int), energy_corr
            close(unit)
            print *, 'Entry added to second_order_energy_corrections.txt'
        else
            print *, 'Entry already exists in second_order_energy_corrections.txt'
        end if
    end subroutine update_second_order_table

    !----------------------------------------------------------------
    ! Write amplitude calculation results to file
    !----------------------------------------------------------------
    subroutine write_amplitude_results_to_file(Z_val, N_val, R_fm_val, state_init, state_final, &
                                             E_init, E_final, E0, intermediate_states_arr, &
                                             E_intermediate_arr, num_inter, A_S_val, A_W_val, &
                                             eta_val, system_type)
        implicit none
        double precision, intent(in) :: Z_val, N_val, R_fm_val, E_init, E_final, E0
        character(len=*), intent(in) :: state_init, state_final, system_type
        integer, intent(in) :: num_inter
        character(len=20), dimension(num_inter), intent(in) :: intermediate_states_arr
        double precision, dimension(num_inter), intent(in) :: E_intermediate_arr
        double precision, intent(in) :: A_S_val, A_W_val, eta_val
        
        character(len=256) :: amplitude_filename
        integer :: unit, ios, i
        character(len=500) :: current_inter_states
        ! Create filename
        call ensure_directory(trim(rad_base) // '/outputs/results/amplitudes')
        amplitude_filename = trim(rad_base) // '/outputs/results/amplitudes/amplitude_calculation_results.txt'
        ! Build intermediate states string for current calculation
        current_inter_states = ''
        do i = 1, num_inter
            if (i == 1) then
                current_inter_states = trim(intermediate_states_arr(i))
            else
                current_inter_states = trim(current_inter_states) // ',' // trim(intermediate_states_arr(i))
            end if
        end do
        
        ! Open file for writing (append if exists, create if not)
        open(newunit=unit, file=amplitude_filename, status='unknown', action='write', &
             position='append', iostat=ios)

        if (ios /= 0) then
            print *, 'Error opening amplitude results file for writing. IOSTAT:', ios
            print *, 'Check that the current directory is writeable and that the file is not locked.'
            return
        end if

        ! Write entry separator and header
        write(unit, '(A)') '==============================================================================='
        write(unit, '(A, I0, A, F0.1, A, F0.4, A)') '# NEW ENTRY - Z=', nint(Z_val), ', N=', N_val, ', R=', R_fm_val, ' fm'
        write(unit, '(A)') '==============================================================================='
        
        ! System information
        write(unit, '(A)') '# SYSTEM PARAMETERS:'
        write(unit, '(A, F8.1)') '#   Atomic Number (Z):              ', Z_val
        write(unit, '(A, F8.1)') '#   Number of Neutrons (N):         ', N_val
        write(unit, '(A, F8.4, A)') '#   Nuclear Radius (R):             ', R_fm_val, ' fm'
        write(unit, '(A, A)') '#   Wavefunction Type:              ', trim(system_type)
        write(unit, '(A)') ''
        
        ! Transition information
        write(unit, '(A)') '# TRANSITION INFORMATION:'
        write(unit, '(A, A)') '#   Initial State:                  ', trim(state_init)
        write(unit, '(A, A)') '#   Final State:                    ', trim(state_final)
        write(unit, '(A, ES15.8E2, A)') '#   Initial State Energy:           ', E_init, ' (relativistic units, m_e = 1)'
        write(unit, '(A, ES15.8E2, A)') '#   Final State Energy:             ', E_final, ' (relativistic units, m_e = 1)'
        write(unit, '(A, ES15.8E2, A)') '#   Electric Field Strength (E0):   ', E0, ' (relativistic units, m_e = 1)'
        write(unit, '(A)') ''
        
        ! Intermediate states information
        write(unit, '(A)') '# INTERMEDIATE STATES:'
        do i = 1, num_inter
            write(unit, '(A, I0, A, A, A, ES15.8E2)') '#   State ', i, ': ', &
                  trim(intermediate_states_arr(i)), ', Energy: ', E_intermediate_arr(i)
        end do
        write(unit, '(A)') ''
        
        ! Results
        write(unit, '(A)') '# CALCULATION RESULTS:'
        write(unit, '(A, ES15.8E2)') '#   Stark Amplitude (A_S):          ', A_S_val
        write(unit, '(A, ES15.8E2)') '#   Weak Amplitude (A_W):           ', A_W_val
        write(unit, '(A, ES15.8E2)') '#   Asymmetry Factor (η):           ', eta_val
        write(unit, '(A)') ''
        
        ! Simple summary line for easy parsing
        write(unit, '(A)') '# SUMMARY:'
        write(unit, '(A, F0.1, A, F0.1, A, F0.4, A, A, A, A, A, A, A, ES15.8E2, A, ES15.8E2, A, ES15.8E2)') &
            'DATA: Z=', Z_val, ' N=', N_val, ' R=', R_fm_val, ' ', &
            trim(state_init), '->', trim(state_final), ' Inter:', trim(current_inter_states), &
            ' A_S=', A_S_val, ' A_W=', A_W_val, ' eta=', eta_val
        
        write(unit, '(A)') ''
        write(unit, '(A)') ''
        
        close(unit)
        print *, 'Amplitude calculation results written to: ', trim(amplitude_filename)
    end subroutine write_amplitude_results_to_file

    subroutine ensure_directory(dir_path)
        implicit none
        character(len=*), intent(in) :: dir_path
        character(len=512) :: command
        integer :: cmd_stat, exit_stat

        command = 'mkdir -p "' // trim(dir_path) // '"'
        call execute_command_line(trim(command), wait=.true., exitstat=exit_stat, cmdstat=cmd_stat)
        if (cmd_stat /= 0 .or. exit_stat /= 0) then
            print *, 'Error creating directory: ', trim(dir_path)
            stop
        end if
    end subroutine ensure_directory
end program compute_radial_integrals
