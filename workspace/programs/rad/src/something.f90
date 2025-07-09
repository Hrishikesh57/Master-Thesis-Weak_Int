program compute_radial_integrals
    use bspqed_potentials
    use bspqed_physconsts
    use bspqed_quadrature       ! <--- We bring in ngl40, nodes40, weights40, etc.
    implicit none

    ! ==============================
    ! Existing declarations
    ! ==============================
    integer :: npoints, i
    integer :: iL, iR
    double precision, allocatable :: r(:), g(:), f(:), V(:)
    double precision :: E_pot, Z, N, R_fm, R_rel, c = 1
    character(len=256) :: filename
    type(weak_potential) :: potential
    ! character(len=6) :: choice  ! <-- Removed unused variable

    ! ==============================
    ! NEW declarations for 2nd order and energy prompt
    ! (Placed at the top so the compiler sees them
    !  before any executable statements.)
    ! ==============================
    character(len=256) :: filename_n, filename_s
    character(len=256) :: matrix_element_filename ! New filename for M_n_s output
    integer :: npoints_n, npoints_s
    double precision, allocatable :: r_n(:), G_n(:), F_n(:)
    double precision, allocatable :: r_s(:), G_s(:), F_s(:)
    double precision :: M_n_s

    ! For the energy state prompt:
    character(len=20) :: order_choice
    character(len=20) :: state_initial, state_intermediate
    double precision :: E_initial, E_intermediate, Delta_E
    integer :: output_unit, io_stat ! For file writing
    character(len=10) :: z_str ! To format Z in filename
    logical :: file_exists ! To check if output file exists

    ! NEW: For amplitude vs energy choice
    character(len=20) :: main_choice, amplitude_choice

    ! NEW: For weak amplitude calculations
    character(len=20) :: state_final
    character(len=256) :: filename_f
    integer :: npoints_f
    double precision, allocatable :: r_f(:), G_f(:), F_f(:)
    double precision :: R_ik, R_kf
    character(len=256) :: weak_amplitude_filename

    ! --------------------------------------------------------------
    ! NEW: Prompt for energy corrections or amplitudes
    ! --------------------------------------------------------------
    print *, 'Do you want to calculate energy corrections or amplitudes? (enter energy or amplitudes)'
    read *, main_choice

    ! --------------------------------------------------------------
    ! Branch based on main choice
    ! --------------------------------------------------------------
    if (trim(main_choice) == 'energy') then
        ! --------------------------------------------------------------
        ! Prompt for 1st or 2nd order (existing code)
        ! --------------------------------------------------------------
        print *, 'Do you want 1st-order or 2nd-order energy corrections? (enter 1 or 2)'
        read *, order_choice

        ! Prompt user for atomic number Z, N, and R in fm
        print *, 'Enter atomic number Z (e.g., 92 for Uranium):'
        read *, Z
        print *, 'Enter number of neutrons N (e.g., 146 for Uranium):'
        read *, N
        print *, 'Enter nuclear radius R in fm (e.g., 5.8571 for Uranium):'
        read *, R_fm

        ! Convert R from fm to relativistic units (1/MeV)
        R_rel = convert_R_fm_to_relativistic()

        ! Initialize the weak potential with the new parameter c
        call potential%init(Z, N, R_rel, c)

        ! Debug print
        print *, 'rho0 is :', potential%rho0
        print *, 'Nuclear radius R in relativistic units :', R_rel
        print *, 'Gfermi rel is : ', Gfermi_relativistic
        print *, 'Cv_pZ_minus_N is : ', potential%Cv_pZ_minus_N

        ! Set filename to read radial functions from (for the 1st order case)
        filename = './radial_function.dat'

        ! -------------------------------------------------------------------------
        ! If user wants 1st order, do the old code exactly as is
        ! -------------------------------------------------------------------------
        if (trim(order_choice) == '1') then

           ! Read the radial functions from the file
           call read_radial_functions(filename, npoints, r, g, f)

           !------------------------------------------------------
           ! Compute potential V(r) for each data point (optional).
           !------------------------------------------------------
           allocate(V(npoints))
           V = 0.0d0
           do i = 1, npoints
               V(i) = potential%eval(r(i))
           end do

           !------------------------------------------------------
           ! E_pot = ∫[0 to R_rel] (g(r)^2 + f(r)^2) * r^2 * V(r) dr
           ! via Gauss–Legendre with N=40
           !------------------------------------------------------
           E_pot = gauss_legendre_integration(0d0, R_rel, npoints, r, g, f, potential, &
                                              .true.)   ! .true. => use (g^2+f^2)*r^2

           print *, 'E_pot (0 to R) using Gauss–Legendre = ', E_pot, ' in relativistic units'

           ! Deallocate arrays
           deallocate(r, g, f, V)

        ! -------------------------------------------------------------------------
        ! Else if user wants 2nd order, do the new integral and energy correction
        ! -------------------------------------------------------------------------
        else if (trim(order_choice) == '2') then
           !---------------------------------------------------------
           ! Prompt for the states (e.g., 1s, 2s) to compute
           ! the energy denominator in the second-order correction.
           !---------------------------------------------------------
           print *, 'Enter the state for the initial energy (e.g., 1s, 2s):'
           read *, state_initial
           print *, 'Enter the state for the intermediate energy (e.g., 2p, 3p, 4p, ..., 10p):'
           read *, state_intermediate

           !---------------------------------------------------------
           ! Prompt for the actual energy values
           !---------------------------------------------------------
           print *, 'Enter the energy for the initial state (', trim(state_initial), ') in relativistic units:'
           read *, E_initial
           print *, 'Enter the energy for the intermediate state (', trim(state_intermediate), ') in relativistic units:'
           read *, E_intermediate

        !Always use fixed values for energies (no special fork)
        !    select case (trim(state_initial))
        !    case ('1s')
        !        E_initial = 0.74152412753477781d0
        !    case ('2s')
        !        E_initial = 0.93311595815167880d0
        !    case default
        !        print *, 'Unknown initial state, defaulting to 1s'
        !        E_initial = 1.0d0
        !    end select

        !    select case (trim(state_intermediate))
        !    case ('2p')
        !        E_intermediate = 0.93305061918473864d0
        !    case ('2p3')
        !         E_intermediate = 0.94197671628375501d0
        !    case ('3p')
        !        E_intermediate = 0.97129559739544835d0
        !    case ('4p')
        !        E_intermediate = 0.98427844487291205d0
        !    case ('5p')
        !        E_intermediate = 0.99012334249205169d0
        !    case ('6p')
        !        E_intermediate = 0.99323248359537508d0
        !    case ('7p')
        !        E_intermediate = 0.99507764832056722d0
        !    case ('8p')
        !        E_intermediate = 0.99626055954639836d0
        !    case ('9p')
        !        E_intermediate = 0.99706371585576403d0
        !    case ('10p')
        !        E_intermediate = 0.99763994924249344d0
        !    case ('11p')
        !         E_intermediate = 0.99814020802280712d0
        !    case ('12p')
        !         E_intermediate = 0.99872729278951022d0
        !    case default
        !        print *, 'Unknown intermediate state, defaulting to 2p'
        !        state_intermediate = '2p'
        !        E_intermediate = 3.0d0
        !    end select

           ! Convert Z to string for filename
           write(z_str, '(I0)') nint(Z) ! nint(Z) converts Z to nearest integer

           filename_n = '/home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int/workspace/programs/rad/' // &
                        trim(z_str) // '_polynomial_table_' // trim(state_intermediate) // '.txt'
           filename_s = '/home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int/workspace/programs/rad/' // &
                        trim(z_str) // '_polynomial_table_' // trim(state_initial) // '.txt'
           print *, 'Reading wave functions from:'
           print *, 'Intermediate state (', trim(state_intermediate), '): ', trim(filename_n)
           print *, 'Initial state (', trim(state_initial), '): ', trim(filename_s)

           call read_radial_functions(filename_n, npoints_n, r_n, G_n, F_n)
           call read_radial_functions(filename_s, npoints_s, r_s, G_s, F_s)

           M_n_s = gauss_legendre_integration_2nd(0d0, R_rel, &
                        npoints_n, r_n, G_n, F_n, &
                        npoints_s, r_s, G_s, F_s, &
                        potential)
           M_n_s = - M_n_s

           print *, 'M_n,s (0 to R) using Gauss–Legendre = ', M_n_s

        !    !----------------------------------------------------------
        !    ! Write Z and M_n_s to a single file (append if exists)
        !    !----------------------------------------------------------
        !    ! Use a fixed filename
        !    matrix_element_filename = '/home/sriman/Documents/Thesis/Z_Matrix_elements.txt'

        !    ! Check if file exists
        !    inquire(file=trim(matrix_element_filename), exist=file_exists)

        !    ! Open for append (creates if not existing)
        !    open(newunit=output_unit, file=trim(matrix_element_filename), status='unknown', &
        !         action='write', position='append', iostat=io_stat)

        !    if (io_stat == 0) then
        !        ! Write header only if the file was newly created
        !        if (.not. file_exists) then
        !            ! Reverted to 2 columns with spacing
        !            write(output_unit, '(A5, A28)') '# Z', ' M_n_s'
        !        end if
        !        ! Write the data with adjusted spacing (Z and M_n_s only)
        !        write(output_unit, '(F5.1, 5x, ES24.16E3)') Z, M_n_s
        !        close(output_unit)
        !        print *, 'Matrix element M_n_s written/appended to: ', trim(matrix_element_filename)
        !    else
        !        print *, 'Error opening file for writing M_n_s: ', trim(matrix_element_filename)
        !        print *, 'IOSTAT = ', io_stat
        !    end if

           Delta_E = (M_n_s**2) / (E_initial - E_intermediate)
           print *, 'For initial state ', trim(state_initial), ' (E = ', E_initial, &
             ') and intermediate state ', trim(state_intermediate), ' (E = ', E_intermediate, '):'
           print *, 'Second-order energy correction ΔE^(2) = ', Delta_E

           !----------------------------------------------------------
           ! Write/update the results table for second order corrections
           !----------------------------------------------------------
           call update_second_order_table(trim(state_initial), trim(state_intermediate), Delta_E) ! <-- Commented out call
        else
           print *, 'Invalid order choice. Exiting.'
           stop
        end if

    else if (trim(main_choice) == 'amplitudes') then
        ! --------------------------------------------------------------
        ! NEW: Prompt for weak or stark amplitude
        ! --------------------------------------------------------------
        print *, 'Do you want weak amplitude or stark amplitude? (enter weak or stark)'
        read *, amplitude_choice

        if (trim(amplitude_choice) == 'weak') then
            ! Prompt user for atomic number Z, N, and R in fm
            print *, 'Enter atomic number Z (e.g., 92 for Uranium):'
            read *, Z
            print *, 'Enter number of neutrons N (e.g., 146 for Uranium):'
            read *, N
            print *, 'Enter nuclear radius R in fm (e.g., 5.8571 for Uranium):'
            read *, R_fm

            ! Convert R from fm to relativistic units (1/MeV)
            R_rel = convert_R_fm_to_relativistic()

            ! Initialize the weak potential with the new parameter c
            call potential%init(Z, N, R_rel, c)

            print *, 'Calculating weak amplitude...'
            print *, 'rho0 is :', potential%rho0
            print *, 'Nuclear radius R in relativistic units :', R_rel
            print *, 'Gfermi rel is : ', Gfermi_relativistic
            print *, 'Cv_pZ_minus_N is : ', potential%Cv_pZ_minus_N

            !---------------------------------------------------------
            ! Prompt for the states (initial, intermediate, final)
            !---------------------------------------------------------
            print *, 'Enter the initial state (e.g., 1s, 2s):'
            read *, state_initial
            print *, 'Enter the intermediate state (e.g., 2p, 3p, 4p, ..., 10p):'
            read *, state_intermediate
            print *, 'Enter the final state (e.g., 1s, 2s):'
            read *, state_final

            ! Convert Z to string for filename
            write(z_str, '(I0)') nint(Z)

            ! Build filenames for all three states
            filename_s = '/home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int/workspace/programs/rad/' // &
                        trim(z_str) // '_polynomial_table_' // trim(state_initial) // '.txt'
            filename_n = '/home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int/workspace/programs/rad/' // &
                        trim(z_str) // '_polynomial_table_' // trim(state_intermediate) // '.txt'
            filename_f = '/home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int/workspace/programs/rad/' // &
                        trim(z_str) // '_polynomial_table_' // trim(state_final) // '.txt'

            print *, 'Reading wave functions from:'
            print *, 'Initial state (', trim(state_initial), '): ', trim(filename_s)
            print *, 'Intermediate state (', trim(state_intermediate), '): ', trim(filename_n)
            print *, 'Final state (', trim(state_final), '): ', trim(filename_f)

            ! Read all three wave functions
            call read_radial_functions(filename_s, npoints_s, r_s, G_s, F_s)
            call read_radial_functions(filename_n, npoints_n, r_n, G_n, F_n)
            call read_radial_functions(filename_f, npoints_f, r_f, G_f, F_f)

            ! Calculate R_ik (initial to intermediate)
            R_ik = gauss_legendre_integration_weak(0d0, R_rel, &
                        npoints_s, r_s, G_s, F_s, &
                        npoints_n, r_n, G_n, F_n, &
                        potential)

            ! Calculate R_kf (intermediate to final)
            R_kf = gauss_legendre_integration_weak(0d0, R_rel, &
                        npoints_f, r_f, G_f, F_f, &
                        npoints_n, r_n, G_n, F_n, &
                        potential)

            print *, 'R_ik (', trim(state_initial), ' -> ', trim(state_intermediate), ') = ', R_ik
            print *, 'R_kf (', trim(state_final), ' -> ', trim(state_intermediate), ') = ', R_kf

            !----------------------------------------------------------
            ! Write results to weak amplitude file
            !----------------------------------------------------------
            weak_amplitude_filename = '/home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int/workspace/programs/rad/weak_amplitudes.txt'
            
            ! Check if file exists to write header
            inquire(file=trim(weak_amplitude_filename), exist=file_exists)
            
            open(newunit=output_unit, file=trim(weak_amplitude_filename), status='unknown', &
                 action='write', position='append', iostat=io_stat)
            
            if (io_stat == 0) then
                ! Write header only if the file was newly created
                if (.not. file_exists) then
                    write(output_unit, '(A5, A10, A15, A15, A15, A25, A25)') &
                        '# Z', ' R_fm', ' Initial', ' Intermediate', ' Final', ' R_ik', ' R_kf'
                end if
                ! Write the data
                write(output_unit, '(F5.1, F10.6, 3(A15), 2(ES25.16E3))') &
                    Z, R_fm, trim(state_initial), trim(state_intermediate), trim(state_final), R_ik, R_kf
                close(output_unit)
                print *, 'Weak amplitude results written to: ', trim(weak_amplitude_filename)
            else
                print *, 'Error opening file for writing weak amplitudes: ', trim(weak_amplitude_filename)
                print *, 'IOSTAT = ', io_stat
            end if

            ! Deallocate arrays
            deallocate(r_s, G_s, F_s)
            deallocate(r_n, G_n, F_n)
            deallocate(r_f, G_f, F_f)

        else if (trim(amplitude_choice) == 'stark') then
            print *, 'Stark amplitude calculation not yet implemented.'
        else
            print *, 'Invalid amplitude choice. Exiting.'
            stop
        end if

    else
        print *, 'Invalid main choice. Please enter energy or amplitudes. Exiting.'
        stop
    end if

contains
    !----------------------------------------------------------------
    ! 1) Convert R from fm to relativistic units (1/MeV)
    !----------------------------------------------------------------
    double precision function convert_R_fm_to_relativistic()
        implicit none
        ! Conversion factor: 1 fm = 1/197 MeV^{-1}
        double precision, parameter :: conversion_factor = 386.1592
        convert_R_fm_to_relativistic = R_fm / conversion_factor
    end function convert_R_fm_to_relativistic

    !----------------------------------------------------------------
    ! 2) Read radial functions from a file: r, g(r), f(r)
    !----------------------------------------------------------------
    subroutine read_radial_functions(filename, npoints, r, g, f)
        implicit none
        character(len=*), intent(in) :: filename
        integer, intent(out) :: npoints
        double precision, allocatable, intent(out) :: r(:), g(:), f(:)

        integer :: unit_num, ios, i
        character(len=256) :: line

        open(newunit=unit_num, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, 'Error opening file: ', filename
            stop
        end if

        ! Skip header line
        read(unit_num, '(A)', iostat=ios) line
        if (ios /= 0) then
            print *, 'Error reading header from file: ', filename
            stop
        end if

        ! Count lines
        npoints = 0
        do
            read(unit_num, '(A)', iostat=ios) line
            if (ios /= 0) exit
            npoints = npoints + 1
        end do

        rewind(unit_num)

        ! Skip header line again
        read(unit_num, '(A)', iostat=ios) line

        allocate(r(npoints), g(npoints), f(npoints))

        do i = 1, npoints
            read(unit_num, *, iostat=ios) r(i), g(i), f(i)
            if (ios /= 0) then
                print *, 'Error reading data at line ', i
                stop
            end if
        end do

        close(unit_num)
    end subroutine read_radial_functions

    !----------------------------------------------------------------
    ! 3) Gauss–Legendre integration from a to b, using ngl40 (1st-order).
    !----------------------------------------------------------------
    function gauss_legendre_integration(a, b, npoints, r_array, g_array, f_array, potential, use_radial_factor) result(Ival)
        use bspqed_quadrature, only: ngl40, nodes40, weights40
        implicit none

        ! Inputs:
        double precision, intent(in) :: a, b
        integer, intent(in) :: npoints
        double precision, dimension(npoints), intent(in) :: r_array, g_array, f_array
        type(weak_potential), intent(in) :: potential
        logical, intent(in) :: use_radial_factor

        ! Output:
        double precision :: Ival

        ! Locals:
        integer :: k, nGL
        double precision :: sum_val, xm, xr, t, w, r_node

        nGL = ngl40
        xm  = 0.5d0 * (b + a)
        xr  = 0.5d0 * (b - a)

        sum_val = 0d0
        do k = 1, nGL
            t = nodes40(k)   ! in [-1,1]
            w = weights40(k)
            r_node = xm + xr*t
            sum_val = sum_val + w * integrand(r_node, npoints, r_array, g_array, f_array, potential, use_radial_factor)
        end do

        Ival = xr * sum_val
    end function gauss_legendre_integration

    !----------------------------------------------------------------
    ! 4) integrand(r): used for 1st-order E_pot
    !----------------------------------------------------------------
    function integrand(rval, local_npoints, r_array, g_array, f_array, local_potential, use_radial_factor) result(val)
        implicit none
        integer, intent(in) :: local_npoints
        double precision, intent(in) :: rval
        double precision, dimension(local_npoints), intent(in) :: r_array, g_array, f_array
        type(weak_potential), intent(in) :: local_potential
        logical, intent(in) :: use_radial_factor
        double precision :: val, g_val, f_val, dr, fraction, V_val
        integer :: iL, iR ! Local loop indices

        if (rval <= r_array(1)) then
            g_val = g_array(1)
            f_val = f_array(1)
        else if (rval >= r_array(local_npoints)) then
            g_val = g_array(local_npoints)
            f_val = f_array(local_npoints)
        else
            iR = 1 ! Initialize iR to avoid potential warning
            do iL = 1, local_npoints-1
                iR = iL + 1
                if (r_array(iL) <= rval .and. rval < r_array(iR)) exit
            end do
            dr = r_array(iR) - r_array(iL)
            fraction = (rval - r_array(iL)) / dr
            g_val = g_array(iL) + fraction*(g_array(iR)-g_array(iL))
            f_val = f_array(iL) + fraction*(f_array(iR)-f_array(iL))
        end if

        V_val = local_potential%eval(rval)
        if (use_radial_factor) then
            val = (g_val*g_val + f_val*f_val) * V_val
        else
            val = (g_val * f_val) * V_val
        end if
    end function integrand

    !----------------------------------------------------------------
    ! NEW: Gauss–Legendre integration for the 2nd-order integral
    !----------------------------------------------------------------
    function gauss_legendre_integration_2nd(a, b, npoints_n, r_n, G_n, F_n, &
                                            npoints_s, r_s, G_s, F_s, &
                                            potential) result(M_val)
        use bspqed_quadrature, only: ngl40, nodes40, weights40
        implicit none
        double precision, intent(in) :: a, b
        integer, intent(in) :: npoints_n, npoints_s
        double precision, dimension(npoints_n), intent(in) :: r_n, G_n, F_n
        double precision, dimension(npoints_s), intent(in) :: r_s, G_s, F_s
        type(weak_potential), intent(in) :: potential
        double precision :: M_val
        integer :: k, nGL
        double precision :: sum_val, xm, xr, t, w, r_node

        nGL = ngl40
        xm  = 0.5d0 * (b + a)
        xr  = 0.5d0 * (b - a)
        sum_val = 0d0
        do k = 1, nGL
            t = nodes40(k)
            w = weights40(k)
            r_node = xm + xr*t
            sum_val = sum_val + w * integrand_2nd(r_node, npoints_n, r_n, G_n, F_n, &
                                                 npoints_s, r_s, G_s, F_s, &
                                                 potential)
        end do
        M_val = xr * sum_val
    end function gauss_legendre_integration_2nd

    !----------------------------------------------------------------
    ! integrand_2nd(r): returns rho(r)*[G_n(r)*F_s(r) + F_n(r)*G_s(r)]
    ! (with r^2 implicitly included if G, F are r*g, r*f).
    !----------------------------------------------------------------
    function integrand_2nd(rval, npoints_n, r_n, G_n, F_n, &
                           npoints_s, r_s, G_s, F_s, &
                           local_potential) result(val2)
        implicit none
        integer, intent(in) :: npoints_n, npoints_s
        double precision, intent(in) :: rval
        double precision, dimension(npoints_n),  intent(in) :: r_n, G_n, F_n
        double precision, dimension(npoints_s), intent(in) :: r_s, G_s, F_s
        type(weak_potential), intent(in) :: local_potential
        double precision :: val2
        double precision :: Gn_val, Fn_val, Gs_val, Fs_val
        double precision :: dr, fraction, rho_val
        integer :: iL, iR ! Local loop indices

        if (rval <= r_n(1)) then
            Gn_val = G_n(1)
            Fn_val = F_n(1)
        else if (rval >= r_n(npoints_n)) then
            Gn_val = G_n(npoints_n)
            Fn_val = F_n(npoints_n)
        else
            iR = 1 ! Initialize iR to avoid potential warning
            do iL = 1, npoints_n-1
                iR = iL + 1
                if (r_n(iL) <= rval .and. rval < r_n(iR)) exit
            end do
            dr = r_n(iR) - r_n(iL)
            fraction = (rval - r_n(iL)) / dr
            Gn_val = G_n(iL) + fraction*(G_n(iR)-G_n(iL))
            Fn_val = F_n(iL) + fraction*(F_n(iR)-F_n(iL))
        end if

        if (rval <= r_s(1)) then
            Gs_val = G_s(1)
            Fs_val = F_s(1)
        else if (rval >= r_s(npoints_s)) then
            Gs_val = G_s(npoints_s)
            Fs_val = F_s(npoints_s)
        else
            iR = 1 ! Initialize iR to avoid potential warning
            do iL = 1, npoints_s-1
                iR = iL + 1
                if (r_s(iL) <= rval .and. rval < r_s(iR)) exit
            end do
            dr = r_s(iR) - r_s(iL)
            fraction = (rval - r_s(iL)) / dr
            Gs_val = G_s(iL) + fraction*(G_s(iR)-G_s(iL))
            Fs_val = F_s(iL) + fraction*(F_s(iR)-F_s(iL))
        end if

        rho_val = local_potential%eval(rval)
        val2 = rho_val * (Gn_val*Fs_val + Fn_val*Gs_val)
    end function integrand_2nd

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
        integer :: nlines, i
        character(len=256) :: line
        character(len=20), allocatable :: init_arr(:), int_arr(:)
        double precision, allocatable :: corr_arr(:)

        filename = '/home/sriman/Documents/Thesis/workspace/programs/rad/second_order_energy_corrections.txt'
        found = .false.
        nlines = 0

        ! First, count lines and read existing entries if file exists
        inquire(file=filename, exist=found)
        if (found) then
            open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
            if (ios /= 0) then
                print *, 'Error opening table file for reading.'
                return
            end if
            ! Count lines
            do
                read(unit, '(A)', iostat=ios) line
                if (ios /= 0) exit
                nlines = nlines + 1
            end do
            rewind(unit)
            if (nlines > 0) then
                allocate(init_arr(nlines), int_arr(nlines), corr_arr(nlines))
                do i = 1, nlines
                    read(unit, *, iostat=ios) file_init, file_int, file_corr
                    if (ios /= 0) exit
                    init_arr(i) = file_init
                    int_arr(i) = file_int
                    corr_arr(i) = file_corr
                end do
            end if
            close(unit)
        else
            nlines = 0
        end if

        ! Check for redundancy
        found = .false.
        if (nlines > 0) then
            do i = 1, nlines
                if (trim(init_arr(i)) == trim(state_init) .and. trim(int_arr(i)) == trim(state_int)) then
                    found = .true.
                    exit
                end if
            end do
        end if

        if (.not. found) then
            ! Append new entry
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

        if (allocated(init_arr)) deallocate(init_arr)
        if (allocated(int_arr)) deallocate(int_arr)
        if (allocated(corr_arr)) deallocate(corr_arr)
    end subroutine update_second_order_table

    !----------------------------------------------------------------
    ! NEW: Gauss–Legendre integration for weak amplitude calculation
    !----------------------------------------------------------------
    function gauss_legendre_integration_weak(a, b, npoints_1, r_1, G_1, F_1, &
                                            npoints_2, r_2, G_2, F_2, &
                                            potential) result(R_val)
        use bspqed_quadrature, only: ngl40, nodes40, weights40
        implicit none
        double precision, intent(in) :: a, b
        integer, intent(in) :: npoints_1, npoints_2
        double precision, dimension(npoints_1), intent(in) :: r_1, G_1, F_1
        double precision, dimension(npoints_2), intent(in) :: r_2, G_2, F_2
        type(weak_potential), intent(in) :: potential
        double precision :: R_val
        integer :: k, nGL
        double precision :: sum_val, xm, xr, t, w, r_node

        nGL = ngl40
        xm  = 0.5d0 * (b + a)
        xr  = 0.5d0 * (b - a)
        sum_val = 0d0
        do k = 1, nGL
            t = nodes40(k)
            w = weights40(k)
            r_node = xm + xr*t
            sum_val = sum_val + w * integrand_weak(r_node, npoints_1, r_1, G_1, F_1, &
                                                 npoints_2, r_2, G_2, F_2, &
                                                 potential)
        end do
        R_val = xr * sum_val
    end function gauss_legendre_integration_weak

    !----------------------------------------------------------------
    ! integrand_weak(r): returns r * rho(r) * [G_1(r)*G_2(r) + F_1(r)*F_2(r)]
    !----------------------------------------------------------------
    function integrand_weak(rval, npoints_1, r_1, G_1, F_1, &
                           npoints_2, r_2, G_2, F_2, &
                           local_potential) result(val_weak)
        implicit none
        integer, intent(in) :: npoints_1, npoints_2
        double precision, intent(in) :: rval
        double precision, dimension(npoints_1), intent(in) :: r_1, G_1, F_1
        double precision, dimension(npoints_2), intent(in) :: r_2, G_2, F_2
        type(weak_potential), intent(in) :: local_potential
        double precision :: val_weak
        double precision :: G1_val, F1_val, G2_val, F2_val
        double precision :: dr, fraction, rho_val
        integer :: iL, iR ! Local loop indices

        ! Interpolate first state
        if (rval <= r_1(1)) then
            G1_val = G_1(1)
            F1_val = F_1(1)
        else if (rval >= r_1(npoints_1)) then
            G1_val = G_1(npoints_1)
            F1_val = F_1(npoints_1)
        else
            iR = 1
            do iL = 1, npoints_1-1
                iR = iL + 1
                if (r_1(iL) <= rval .and. rval < r_1(iR)) exit
            end do
            dr = r_1(iR) - r_1(iL)
            fraction = (rval - r_1(iL)) / dr
            G1_val = G_1(iL) + fraction*(G_1(iR)-G_1(iL))
            F1_val = F_1(iL) + fraction*(F_1(iR)-F_1(iL))
        end if

        ! Interpolate second state
        if (rval <= r_2(1)) then
            G2_val = G_2(1)
            F2_val = F_2(1)
        else if (rval >= r_2(npoints_2)) then
            G2_val = G_2(npoints_2)
            F2_val = F_2(npoints_2)
        else
            iR = 1
            do iL = 1, npoints_2-1
                iR = iL + 1
                if (r_2(iL) <= rval .and. rval < r_2(iR)) exit
            end do
            dr = r_2(iR) - r_2(iL)
            fraction = (rval - r_2(iL)) / dr
            G2_val = G_2(iL) + fraction*(G_2(iR)-G_2(iL))
            F2_val = F_2(iL) + fraction*(F_2(iR)-F_2(iL))
        end if

        rho_val = local_potential%eval(rval)
        val_weak = rval * rho_val * (G1_val*G2_val + F1_val*F2_val)
    end function integrand_weak
end program compute_radial_integrals
