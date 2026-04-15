program print_polynomial_table_relativistic
   use bspqed_dkb
   use bspqed_coefficients
   use bspqed_utilities, only: kntseq1
   use bspqed_potentials, only: sphere_potential
   use bspqed_quadrature, only: nodes20, weights20
   use bspqed_physconsts        ! This module should define alpha and electron_rest_energy
   use bspwf_constants, only: bohr_radius
   use bspwf_dirac, only: dirac, au_radialg => radialg, au_radialf => radialf
   use bspwf_hartree_fock, only: dhfc_dir_matrix, dhfc_matrix
   use bspwf_knots, only: knots_01
   use bspwf_potentials, only: fermi_potential
   use bspwf_quadrature, only: nodes10, weights10
   implicit none

   ! Constants
   double precision, parameter :: pi_local = 4.0d0 * atan(1.0d0)
   double precision, parameter :: output_rmax_rel = 200.0d0
   character(len=*), parameter :: single_output_dir = 'outputs/wavefunctions/single_electron'
   character(len=*), parameter :: scf_output_dir = 'outputs/wavefunctions/scf'
   integer :: wavefunction_type

   call prompt_wavefunction_type(wavefunction_type)

   select case (wavefunction_type)
   case (1)
       call run_single_electron_calculation()
   case (2)
       call run_multi_electron_calculation()
   end select

contains

   subroutine prompt_wavefunction_type(wavefunction_type)
       integer, intent(out) :: wavefunction_type

       print *, 'Choose the type of wavefunction calculation:'
       print *, '1 - Single-electron wavefunctions'
       print *, '2 - Multi-electron wavefunctions'
       print *, 'Enter your choice (1 or 2):'
       read(*,*) wavefunction_type

       if (wavefunction_type /= 1 .and. wavefunction_type /= 2) then
           print *, 'Error: Invalid choice. Please enter 1 or 2.'
           stop
       endif
   end subroutine prompt_wavefunction_type

   subroutine run_single_electron_calculation()
       integer :: znucl, order, kappa, nmain, l
       integer :: nbspl, nspectrum, na
       character(len=10) :: state_str
       double precision :: rms_fm, rcavity
       double precision, allocatable :: knots(:), coefs(:), nodes(:), weights(:)
       double precision, allocatable :: spectrum(:), coefmat(:,:)
       type(sphere_potential) :: potential

       print *, 'Single-electron calculation selected.'
       print *, 'Enter the value for znucl (atomic number Z):'
       read(*,*) znucl

       print *, 'Enter the state (choose among: 1s, 2s, 3s, 4s, 5s, 6s, 7s, 8s, 2p, 2p3, 3p, 4p, 5p, 6p, 6p3, 7p, 7p3, 8p, 9p, 10p):'
       read(*,*) state_str
       call parse_state(state_str, kappa, nmain, l)

       print *, 'Enter the value for rms radius (in fm):'
       read(*,*) rms_fm

       if (znucl <= 0) then
           print *, 'Error: znucl must be positive.'
           stop
       endif

       call initialize_single_electron_problem(znucl, rms_fm, order, rcavity, knots, nodes, weights, potential)

       nbspl = size(knots) - order
       nspectrum = 2 * (nbspl - 2)
       na = nspectrum / 2 + (nmain - l)

       allocate(coefs(2 * nbspl))
       allocate(spectrum(nspectrum))
       allocate(coefmat(2 * nbspl, nspectrum))

       print *, 'nspectrum = ', nspectrum
       call solve(knots, order, kappa, potential, nodes, weights, spectrum, coefmat)
       print *, 'Eigenvalue (energy) for chosen state = ', spectrum(na)

       coefs = coefmat(:, na)
       call enforce_phase_convention_single(knots, order, kappa, l, coefs)
       call write_single_wavefunction_output(znucl, state_str, knots, order, kappa, coefs, rcavity)
   end subroutine run_single_electron_calculation

   subroutine run_multi_electron_calculation()
       integer :: znucl, order, kappa, nmain, l
       integer :: num_electrons, ncores, nvalences, nbsplines
       integer, allocatable :: core_kappas(:), core_nmains(:), valence_kappas(:), valence_nmains(:)
       double precision :: rms_fm, ecore_total, rcavity
       double precision, allocatable :: knots(:), coefs(:), core_energies(:), valence_energies(:)
       double precision, allocatable :: core_coefficients(:,:), valence_coefficients(:,:)
       type(fermi_potential) :: nuclpot_scf
       character(len=10) :: state_str
       real :: start, finish

       print *, 'Enter the number of electrons in the system (1-90):'
       read(*,*) num_electrons

       print *, 'Enter the atomic number Z:'
       read(*,*) znucl

       if (num_electrons < 1 .or. num_electrons > znucl) then
           print *, 'Error: Number of electrons must be between 1 and Z=', znucl
           stop
       endif

       if (num_electrons > 90) then
           print *, 'Error: Number of electrons must be <= 90 for this implementation.'
           stop
       endif

       print *, 'Multi-electron calculation selected with', num_electrons, 'electrons.'
       print *, 'This creates a ', trim(ion_name(znucl, num_electrons)), ' system.'

       print *, 'Enter the nuclear RMS radius (in fm):'
       read(*,*) rms_fm

       print *, 'Enter the state for wavefunction output (1s, 2s, 2p, 2p3, 3s, 3p, 3p3, 4s, 5s, 5p, 5p3, 6s, 6p, 6p3, 7s, 7p, 7p3, 8s):'
       read(*,*) state_str
       call parse_state(state_str, kappa, nmain, l)

       call setup_electron_configuration(num_electrons, nmain, kappa, core_kappas, core_nmains, &
                                         valence_kappas, valence_nmains, ncores, nvalences)
       call print_configuration_summary(znucl, num_electrons, rms_fm, state_str, &
                                        core_nmains, core_kappas, valence_nmains, valence_kappas)

       call initialize_scf_problem(znucl, nmain, rms_fm, order, rcavity, knots, nuclpot_scf, nbsplines)

       allocate(core_energies(ncores), core_coefficients(2 * nbsplines, ncores))
       core_energies = 0.0d0
       core_coefficients = 0.0d0

       print *, "# Starting SCF calculation for core electrons..."
       call cpu_time(start)
       call starting_solutions_scf(knots, order, nuclpot_scf, core_kappas, core_nmains, nodes10, weights10, &
                                   core_energies, core_coefficients, .true.)
       call cpu_time(finish)

       write (*, *) "# ", repeat("-", 86)
       write (*, *) "# Starting approximations completed!"
       write (*, '(A,F8.3,A)') "# Duration: ", finish - start, " seconds"
       call print_energy_report('# Core state energies:', core_nmains, core_kappas, core_energies)
       write (*, *) "# ", repeat("-", 86)

       call cpu_time(start)
       call refining_solutions_scf(knots, order, nuclpot_scf, core_kappas, core_nmains, nodes10, weights10, &
                                   core_energies, core_coefficients, .true.)
       call cpu_time(finish)

       write (*, *) "# ", repeat("-", 86)
       write (*, *) "# SCF refinement completed!"
       write (*, '(A,F8.3,A)') "# Duration: ", finish - start, " seconds"
       call print_energy_report('# Final core state energies:', core_nmains, core_kappas, core_energies)

       ecore_total = total_core_energy_scf(knots, order, core_kappas, core_energies, core_coefficients, &
                                           nodes10, weights10)
       write (*, '(A,ES15.8,A,ES15.8)') "# Total core energy (au): ", ecore_total, &
                                         ", Total core energy (rel): ", alpha**2 * ecore_total
       write (*, *) "# ", repeat("-", 86)

       allocate(valence_energies(nvalences), valence_coefficients(2 * nbsplines, nvalences))
       valence_energies = 0.0d0
       valence_coefficients = 0.0d0

       if (nvalences > 0) then
           call valence_solutions_scf(knots, order, nuclpot_scf, core_kappas, core_coefficients, &
                                      valence_kappas, valence_nmains, nodes10, weights10, &
                                      valence_energies, valence_coefficients)
           call print_energy_report('# Valence state energies:', valence_nmains, valence_kappas, valence_energies)
           write (*, *) "# ", repeat("-", 86)
       endif

       allocate(coefs(2 * nbsplines))
       call extract_target_wavefunction_scf(kappa, nmain, core_kappas, core_nmains, core_coefficients, &
                                            valence_kappas, valence_nmains, valence_coefficients, &
                                            ncores, nvalences, coefs)
       call enforce_phase_convention_scf(knots, order, kappa, l, coefs)
       call write_scf_wavefunction_output(znucl, num_electrons, rms_fm, state_str, knots, order, kappa, coefs, rcavity)
   end subroutine run_multi_electron_calculation

   !> Build the single-electron B-spline problem in relativistic units.
   !> The later call to `solve()` performs the full generalized eigensystem
   !> diagonalization for the finite-size nuclear Dirac Hamiltonian.
   subroutine initialize_single_electron_problem(znucl, rms_fm, order, rcavity, knots, nodes, weights, potential)
       integer, intent(in) :: znucl
       integer, intent(out) :: order
       double precision, intent(in) :: rms_fm
       double precision, intent(out) :: rcavity
       double precision, allocatable, intent(out) :: knots(:), nodes(:), weights(:)
       type(sphere_potential), intent(out) :: potential

       integer :: npoints, nknots
       double precision :: rnucl

       order = 9
       rnucl = sqrt(5.d0 / 3.d0) * rms_fm / electron_reduced_compton
       rcavity = 2.5d0 / alpha

       npoints = 1000
       nknots = npoints + 2 * (order - 1)
       allocate(knots(nknots))
       call kntseq1(rnucl, rcavity, npoints / 10, order, knots)
       if (knots(1) <= 0.0d0) knots(1) = 1.0d-6

       allocate(nodes(size(nodes20)), weights(size(weights20)))
       nodes = nodes20
       weights = weights20

       potential%zalpha = znucl * alpha
       potential%rad = rnucl
   end subroutine initialize_single_electron_problem

   !> Build the SCF knot grid and Fermi potential in atomic units.
   !> The Fermi radius is computed in fm first and then converted to bohr so
   !> the Dirac and Hartree-Fock routines use a consistent length scale.
   subroutine initialize_scf_problem(znucl, nmain, rms_fm, order, rcavity, knots, nuclpot_scf, nbsplines)
       integer, intent(in) :: znucl, nmain
       integer, intent(out) :: order, nbsplines
       double precision, intent(in) :: rms_fm
       double precision, intent(out) :: rcavity
       double precision, allocatable, intent(out) :: knots(:)
       type(fermi_potential), intent(out) :: nuclpot_scf

       integer :: npoints, nknots
       double precision :: afermi_fm, cfermi_fm

       order = 9
       npoints = 100
       nknots = npoints + 2 * (order - 1)
       nbsplines = nknots - order

       allocate(knots(nknots))
       knots = 0.0d0

       afermi_fm = 2.3d0 / (4.0d0 * log(3.0d0))
       cfermi_fm = sqrt((5.0d0 * rms_fm**2 - 7.0d0 * (pi_local * afermi_fm)**2) / 3.0d0)
       rcavity = choose_scf_rcavity(nmain)

       call knots_01(cfermi_fm / bohr_radius, rcavity, npoints / 10, order, knots)
       nuclpot_scf%znucl = znucl
       nuclpot_scf%a = afermi_fm / bohr_radius
       nuclpot_scf%c = cfermi_fm / bohr_radius
   end subroutine initialize_scf_problem

   subroutine print_configuration_summary(znucl, num_electrons, rms_fm, state_str, core_nmains, core_kappas, &
                                          valence_nmains, valence_kappas)
       integer, intent(in) :: znucl, num_electrons
       integer, dimension(:), intent(in) :: core_nmains, core_kappas, valence_nmains, valence_kappas
       double precision, intent(in) :: rms_fm
       character(len=*), intent(in) :: state_str

       write (*, *) "# ", repeat("=", 86)
       write (*, '(A,I0)') "# Atomic number Z = ", znucl
       write (*, '(A,I0)') "# Number of electrons = ", num_electrons
       write (*, '(A,F8.4)') "# Nuclear RMS radius = ", rms_fm, " fm"
       write (*, '(A,A)') "# Target state for output: ", trim(state_str)
       call print_shells('# Core electron configuration:', core_nmains, core_kappas)
       if (size(valence_nmains) > 0) call print_shells('# Valence electron configuration:', valence_nmains, valence_kappas)
       write (*, *) "# ", repeat("-", 86)
   end subroutine print_configuration_summary

   subroutine print_shells(title, nmains, kappas)
       character(len=*), intent(in) :: title
       integer, dimension(:), intent(in) :: nmains, kappas

       integer :: i

       write (*, '(A)') trim(title)
       do i = 1, size(nmains)
           write (*, '(A,I0,A,I0,A,I0)') "#   ", i, ": n=", nmains(i), ", kappa=", kappas(i)
       end do
   end subroutine print_shells

   subroutine print_energy_report(title, nmains, kappas, energies)
       character(len=*), intent(in) :: title
       integer, dimension(:), intent(in) :: nmains, kappas
       double precision, dimension(:), intent(in) :: energies

       integer :: i

       write (*, '(A)') trim(title)
       do i = 1, size(energies)
           write (*, '(A,I0,A,I0,A,ES15.8,A,ES15.8)') "#   n=", nmains(i), ", kappa=", kappas(i), &
                                                       ", E_au=", energies(i), ", E_rel=", alpha**2 * energies(i)
       end do
   end subroutine print_energy_report

   subroutine open_wavefunction_output_file(output_dir, file_tag, znucl, state_str, unit_num, filename)
       character(len=*), intent(in) :: output_dir, file_tag, state_str
       integer, intent(in) :: znucl
       integer, intent(out) :: unit_num
       character(len=256), intent(out) :: filename

       character(len=16) :: znucl_str
       integer :: iostat

       call ensure_directory(output_dir)
       write(znucl_str, '(I0)') znucl
       filename = trim(output_dir) // '/' // trim(znucl_str) // '_' // trim(file_tag) // '_' // trim(state_str) // '.txt'
       open(newunit=unit_num, file=filename, status='replace', action='write', iostat=iostat)
       if (iostat /= 0) then
           print *, 'Error opening file: ', trim(filename)
           stop
       endif

       print *, 'File opened successfully: ', trim(filename)
   end subroutine open_wavefunction_output_file

   subroutine write_single_wavefunction_output(znucl, state_str, knots, order, kappa, coefs, rcavity)
       integer, intent(in) :: znucl, order, kappa
       character(len=*), intent(in) :: state_str
       double precision, intent(in) :: rcavity
       double precision, dimension(:), intent(in) :: knots, coefs

       integer :: i, inbv, unit_num
       integer, parameter :: npoints_output = 1000
       character(len=256) :: filename
       double precision :: gx, fx, x, x_start, x_end, x_step

       call open_wavefunction_output_file(single_output_dir, 'polynomial_table', znucl, state_str, unit_num, filename)

       write(unit_num, '(A)') 'r (relativistic units, m_e=1)    G(r)           F(r)'
       print *, '   r (relativistic units, m_e=1)    G(r)           F(r)'

       x_start = 1.0d-6
       x_end = min(output_rmax_rel, rcavity - 1.0d-6)
       x_step = (x_end - x_start) / (npoints_output - 1)

       do i = 1, npoints_output
           x = x_start + (i - 1) * x_step
           inbv = 1
           gx = radialg(knots, order, kappa, coefs, x, inbv)
           fx = radialf(knots, order, kappa, coefs, x, inbv)
           write(unit_num, '(F12.6, 3X, ES18.10E3, 3X, ES18.10E3)') x, gx, fx
           print '(F12.6, 3X, ES18.10E3, 3X, ES18.10E3)', x, gx, fx
       end do

       close(unit_num)
   end subroutine write_single_wavefunction_output

   !> Write the SCF orbital on a logarithmic atomic-unit grid.
   !> Conversion to relativistic units is deferred to the final formatting
   !> step so the internal evaluation stays in atomic units.
   subroutine write_scf_wavefunction_output(znucl, num_electrons, rms_fm, state_str, knots, order, kappa, coefs, rcavity)
       integer, intent(in) :: znucl, num_electrons, order, kappa
       character(len=*), intent(in) :: state_str
       double precision, intent(in) :: rms_fm, rcavity
       double precision, dimension(:), intent(in) :: knots, coefs

       integer :: i, inbv, unit_num
       integer, parameter :: npoints_output = 2000
       character(len=256) :: filename
       double precision :: gx, fx, x_au, x_nat, scale
       double precision :: x_start, x_end, log_x_start, log_x_end, log_x_step

       call open_wavefunction_output_file(scf_output_dir, 'scf_polynomial_table', znucl, state_str, unit_num, filename)

       write(unit_num, '(A)') 'r (relativistic units, m_e=1)    G(r)           F(r)'
       write(unit_num, '(A,A)') '# State: ', trim(state_str)
       write(unit_num, '(A,I0)') '# Z = ', znucl
       write(unit_num, '(A,I0)') '# Number of electrons = ', num_electrons
       write(unit_num, '(A,F8.4)') '# RMS radius = ', rms_fm
       write(unit_num, '(A)') '#'

       x_start = 1.0d-6
       x_end = rcavity - 1.0d-6
       log_x_start = log(x_start)
       log_x_end = log(x_end)
       log_x_step = (log_x_end - log_x_start) / (npoints_output - 1)
       scale = sqrt(alpha)

       do i = 1, npoints_output
           x_au = exp(log_x_start + (i - 1) * log_x_step)
           inbv = 1
           gx = scale * au_radialg(knots, order, kappa, coefs, x_au, inbv)
           fx = scale * au_radialf(knots, order, kappa, coefs, x_au, inbv)
           x_nat = x_au / alpha
           write(unit_num, '(F12.6, 3X, ES18.10E3, 3X, ES18.10E3)') x_nat, gx, fx
       end do

       close(unit_num)

       write (*, *) "# ", repeat("=", 86)
       write (*, *) "# Multi-electron calculation completed successfully!"
       write (*, '(A,A)') "# Wavefunction written to: ", trim(filename)
       write (*, *) "# ", repeat("=", 86)
   end subroutine write_scf_wavefunction_output

   !> Choose an SCF cavity radius large enough for diffuse high-n states.
   double precision function choose_scf_rcavity(nmain)
       integer, intent(in) :: nmain

       choose_scf_rcavity = max(50.0d0, 4.0d0 * dble(nmain * nmain))
   end function choose_scf_rcavity

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

   !> Enforce a deterministic phase convention for single-electron output.
   subroutine enforce_phase_convention_single(knots, order, kappa, l, coefs)
       double precision, dimension(:), intent(in) :: knots
       integer, intent(in) :: order, kappa, l
       double precision, dimension(:), intent(in out) :: coefs

       double precision :: desired_sign, xref, gref
       integer :: inbv, itry

       desired_sign = 1.0d0
       if (mod(l, 2) == 1) desired_sign = -1.0d0

       xref = 5.0d-2
       gref = 0.0d0
       do itry = 1, 20
           inbv = 1
           gref = radialg(knots, order, kappa, coefs, xref, inbv)
           if (abs(gref) > 1.0d-14) exit
           xref = min(2.0d0 * xref, 5.0d0)
       end do

       if (desired_sign * gref < 0.0d0) coefs = -coefs
   end subroutine enforce_phase_convention_single

   !> Enforce a deterministic phase convention for SCF output.
   subroutine enforce_phase_convention_scf(knots, order, kappa, l, coefs)
       double precision, dimension(:), intent(in) :: knots
       integer, intent(in) :: order, kappa, l
       double precision, dimension(:), intent(in out) :: coefs

       double precision :: desired_sign, xref, gref
       integer :: inbv, itry

       desired_sign = 1.0d0
       if (mod(l, 2) == 1) desired_sign = -1.0d0

       xref = 5.0d-2 * alpha
       gref = 0.0d0
       do itry = 1, 20
           inbv = 1
           gref = au_radialg(knots, order, kappa, coefs, xref, inbv)
           if (abs(gref) > 1.0d-14) exit
           xref = min(2.0d0 * xref, 5.0d0 * alpha)
       end do

       if (desired_sign * gref < 0.0d0) coefs = -coefs
   end subroutine enforce_phase_convention_scf

   !> Parse a spectroscopic state label into `kappa`, `nmain`, and `l`.
   subroutine parse_state(state_str, kappa, nmain, l)
       character(len=*), intent(in) :: state_str
       integer, intent(out) :: kappa, nmain, l

       integer :: split_idx, ios
       character(len=32) :: state, suffix

       state = to_lower_copy(adjustl(trim(state_str)))
       split_idx = scan(state, 'spdf')
       if (split_idx <= 1) then
           call stop_invalid_state(state_str)
       endif

       read(state(:split_idx - 1), *, iostat=ios) nmain
       if (ios /= 0) then
           call stop_invalid_state(state_str)
       endif

       suffix = trim(state(split_idx:))
       select case (suffix)
       case ('s')
           l = 0
           kappa = -1
       case ('p', 'p1')
           l = 1
           kappa = 1
       case ('p3')
           l = 1
           kappa = -2
       case ('d', 'd3')
           l = 2
           kappa = 2
       case ('d5')
           l = 2
           kappa = -3
       case ('f', 'f5')
           l = 3
           kappa = 3
       case ('f7')
           l = 3
           kappa = -4
       case default
           call stop_invalid_state(state_str)
       end select

       if (nmain <= l) call stop_invalid_state(state_str)
   end subroutine parse_state

   function to_lower_copy(text) result(lowered)
       character(len=*), intent(in) :: text
       character(len=len(text)) :: lowered

       integer :: i, code

       lowered = text
       do i = 1, len(text)
           code = iachar(lowered(i:i))
           if (code >= iachar('A') .and. code <= iachar('Z')) lowered(i:i) = achar(code + 32)
       end do
   end function to_lower_copy

   subroutine stop_invalid_state(state_str)
       character(len=*), intent(in) :: state_str

       print *, 'Invalid state: ', trim(state_str)
       print *, 'Use labels such as 1s, 2p, 2p3, 4d, 4d5, 5f, or 5f7.'
       stop
   end subroutine stop_invalid_state

   !> Return a compact ion label for terminal output.
   function ion_name(z, nelec) result(name)
       integer, intent(in) :: z, nelec
       character(len=20) :: name
       character(len=2) :: element_symbols(20) = ['H ', 'He', 'Li', 'Be', 'B ', 'C ', &
                                                  'N ', 'O ', 'F ', 'Ne', 'Na', 'Mg', &
                                                  'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', &
                                                  'K ', 'Ca']
       integer :: charge
       
       charge = z - nelec
       if (z <= 20) then
           if (charge == 0) then
               name = trim(element_symbols(z))
           else
               write(name, '(A,I0,"+")') trim(element_symbols(z)), charge
           endif
       else
           if (charge == 0) then
               write(name, '("Z=",I0)') z
           else
               write(name, '("Z=",I0,"(",I0,"+)")') z, charge
           endif
       endif
   end function ion_name

   !> Build the occupied shell list for the requested electron count.
   !> The shells follow a relativistic Aufbau ordering up to 8s, with
   !> capacities given by `2*abs(kappa)`. The target shell is forced into the
   !> returned configuration even when it is the first unfilled shell.
   subroutine setup_electron_configuration(num_electrons, target_nmain, target_kappa, &
                                          core_kappas, core_nmains, &
                                          valence_kappas, valence_nmains, &
                                          ncores, nvalences)
       integer, intent(in) :: num_electrons, target_nmain, target_kappa
       integer, dimension(:), allocatable, intent(out) :: core_kappas, core_nmains
       integer, dimension(:), allocatable, intent(out) :: valence_kappas, valence_nmains
       integer, intent(out) :: ncores, nvalences

       integer, parameter :: max_aufbau = 32
       integer, dimension(max_aufbau), parameter :: aufbau_n = (/ &
           1, &                                       ! 1s
           2, 2, 2, &                                 ! 2s, 2p1/2, 2p3/2
           3, 3, 3, &                                 ! 3s, 3p1/2, 3p3/2
           4, 3, 3, 4, 4, &                           ! 4s, 3d3/2, 3d5/2, 4p1/2, 4p3/2
           5, 4, 4, 5, 5, &                           ! 5s, 4d3/2, 4d5/2, 5p1/2, 5p3/2
           6, 4, 4, 5, 5, 6, 6, &                     ! 6s, 4f5/2, 4f7/2, 5d3/2, 5d5/2, 6p1/2, 6p3/2
           7, 5, 5, 6, 6, 7, 7, &                     ! 7s, 5f5/2, 5f7/2, 6d3/2, 6d5/2, 7p1/2, 7p3/2
           8 /)                                        ! 8s
       integer, dimension(max_aufbau), parameter :: aufbau_kappa = (/ &
          -1, &                                        ! 1s
          -1,  1, -2, &                                ! 2s, 2p1/2, 2p3/2
          -1,  1, -2, &                                ! 3s, 3p1/2, 3p3/2
          -1,  2, -3,  1, -2, &                        ! 4s, 3d3/2, 3d5/2, 4p1/2, 4p3/2
          -1,  2, -3,  1, -2, &                        ! 5s, 4d3/2, 4d5/2, 5p1/2, 5p3/2
          -1,  3, -4,  2, -3,  1, -2, &                ! 6s, 4f5/2, 4f7/2, 5d3/2, 5d5/2, 6p1/2, 6p3/2
          -1,  3, -4,  2, -3,  1, -2, &                ! 7s, 5f5/2, 5f7/2, 6d3/2, 6d5/2, 7p1/2, 7p3/2
          -1 /)                                         ! 8s

       integer :: electrons_left, shell_capacity, electrons_in_shell, idx, i
       logical :: target_found
       integer, dimension(max_aufbau) :: temp_core_kappas, temp_core_nmains

       ncores = 0
       electrons_left = num_electrons
       target_found = .false.

       do idx = 1, max_aufbau
           shell_capacity = 2 * abs(aufbau_kappa(idx))
           electrons_in_shell = min(electrons_left, shell_capacity)

           if (electrons_in_shell > 0 .or. &
               (aufbau_n(idx) == target_nmain .and. aufbau_kappa(idx) == target_kappa)) then
               ncores = ncores + 1
               temp_core_nmains(ncores) = aufbau_n(idx)
               temp_core_kappas(ncores) = aufbau_kappa(idx)

               if (aufbau_n(idx) == target_nmain .and. aufbau_kappa(idx) == target_kappa) then
                   target_found = .true.
               endif

               electrons_left = electrons_left - electrons_in_shell
           endif

           if (electrons_left <= 0 .and. target_found) exit
       end do

       if (.not. target_found) then
           print *, 'Error: Target state not reachable with given electron configuration'
           print *, 'Target: n=', target_nmain, ', kappa=', target_kappa
           stop
       endif

       print *, '# Debug: Configuration setup for', num_electrons, 'electrons:'
       do i = 1, ncores
           print *, '#   Shell', i, ': n=', temp_core_nmains(i), ', kappa=', temp_core_kappas(i)
       end do

       nvalences = 0

       allocate(core_kappas(ncores), core_nmains(ncores))
       core_kappas(1:ncores) = temp_core_kappas(1:ncores)
       core_nmains(1:ncores) = temp_core_nmains(1:ncores)

       allocate(valence_kappas(0), valence_nmains(0))
   end subroutine setup_electron_configuration

   !> Extract the requested orbital coefficients from the SCF result arrays.
   subroutine extract_target_wavefunction_scf(target_kappa, target_nmain, &
                                              core_kappas, core_nmains, core_coefficients, &
                                              valence_kappas, valence_nmains, valence_coefficients, &
                                              ncores, nvalences, coefs)
       integer, intent(in) :: target_kappa, target_nmain
       integer, dimension(:), intent(in) :: core_kappas, core_nmains
       integer, dimension(:), intent(in) :: valence_kappas, valence_nmains
       double precision, dimension(:,:), intent(in) :: core_coefficients, valence_coefficients
       integer, intent(in) :: ncores, nvalences
       double precision, dimension(:), intent(out) :: coefs
       
       integer :: i
       logical :: found
       
       found = .false.
       
       do i = 1, ncores
           if (core_kappas(i) == target_kappa .and. core_nmains(i) == target_nmain) then
               coefs = core_coefficients(:, i)
               found = .true.
               print *, '# Target state found in core electrons'
               exit
           endif
       end do
       
       if (.not. found .and. nvalences > 0) then
           do i = 1, nvalences
               if (valence_kappas(i) == target_kappa .and. valence_nmains(i) == target_nmain) then
                   coefs = valence_coefficients(:, i)
                   found = .true.
                   print *, '# Target state found in valence electrons'
                   exit
               endif
           end do
       endif
       
       if (.not. found) then
           print *, 'Error: Target state not found in electron configuration'
           stop
       endif
   end subroutine extract_target_wavefunction_scf

   !> Generate the initial SCF orbitals by starting from the bare nuclear
   !> potential and then iterating with damped direct Hartree-Fock updates.
   subroutine starting_solutions_scf(knots, order, potential_scf, core_kappas, core_nmains, nodes, weights, &
                                     core_energies, core_coefficients, print_output)
       double precision, dimension(:), intent(in) :: knots
       integer, intent(in) :: order
       type(fermi_potential), intent(in) :: potential_scf
       integer, dimension(:), intent(in) :: core_kappas, core_nmains
       double precision, dimension(:), intent(in) :: nodes, weights
       double precision, dimension(:), intent(out) :: core_energies
       double precision, dimension(:, :), intent(out) :: core_coefficients
       logical, intent(in) :: print_output

       double precision, dimension(2*(size(knots)-order-2), 2*(size(knots)-order-2), size(core_kappas)) :: hfdir, hfdir_curr
       double precision, dimension(2*(size(knots)-order-2), 2*(size(knots)-order-2)) :: ahf
       double precision, dimension(2*(size(knots)-order), 2*(size(knots)-order-2)) :: coefmat
       double precision, dimension(size(core_coefficients, 1), size(core_coefficients, 2)) :: prev_core_coefs
       double precision, dimension(size(core_kappas), 2) :: epsbs
       double precision, dimension(2*(size(knots)-order-2)) :: spectrum
       double precision, dimension(2*(size(knots)-order)) :: coefb, coefc
       double precision, dimension(size(core_kappas)) :: relerrs
       double precision :: tolrel, epsb, abserr, tolerance
       integer, dimension(1) :: temp
       integer :: maxstep, ncores, nspectrum, step, i, nelectrons, j, kapb, nmnb, lb, nb, kapc, imax

       tolrel = 5.0d-04
       maxstep = 50
       ncores = size(core_kappas)
       nspectrum = 2 * (size(knots) - order - 2)

       step = 1
       do i = 1, ncores
           kapb = core_kappas(i)
           nmnb = core_nmains(i)
           hfdir(:, :, i) = 0.0d0
           call dirac(knots, order, kapb, potential_scf, nodes, weights, spectrum, coefmat, hfmat=hfdir(:, :, i))
           lb = (abs(2 * kapb + 1) - 1) / 2
           nb = nspectrum / 2 + nmnb - lb
           epsb = spectrum(nb) - ialpha**2
           coefb = coefmat(:, nb)
           core_coefficients(:, i) = coefb
           epsbs(i, 1) = epsb
           if (print_output) write (*, *) step, nmnb, kapb, epsb
       end do
       if (print_output) write (*, *) ""

       do step = 2, maxstep
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
               call dirac(knots, order, kapb, potential_scf, nodes, weights, spectrum, coefmat, hfmat=hfdir(:, :, i))
               lb = (abs(2 * kapb + 1) - 1) / 2
               nb = nspectrum / 2 + nmnb - lb
               epsb = spectrum(nb) - ialpha**2
               coefb = coefmat(:, nb)
               core_coefficients(:, i) = coefb
               epsbs(i, 2) = epsb
               if (print_output) write (*, *) step, nmnb, kapb, epsb
           end do
           if (print_output) write (*, *) ""

           if (step > 2) then
               do i = 1, ncores
                   relerrs(i) = abs((epsbs(i, 1) - epsbs(i, 2)) / epsbs(i, 2))
               end do
               temp = maxloc(relerrs)
               imax = temp(1)
               abserr = abs(epsbs(imax, 2) - epsbs(imax, 1))
               tolerance = tolrel * abs(epsbs(imax, 2))
               if (abserr <= tolerance) exit
           endif

           epsbs(:, 1) = epsbs(:, 2)
       end do

       core_energies = epsbs(:, 2)
   end subroutine starting_solutions_scf

   !> Refine the SCF solution by repeatedly updating the most slowly
   !> converging orbital until the energy tolerance is met.
   subroutine refining_solutions_scf(knots, order, potential_scf, core_kappas, core_nmains, nodes, weights, &
                                     core_energies, core_coefficients, print_output)
       double precision, dimension(:), intent(in) :: knots
       integer, intent(in) :: order
       type(fermi_potential), intent(in) :: potential_scf
       integer, dimension(:), intent(in) :: core_kappas, core_nmains
       double precision, dimension(:), intent(in) :: nodes, weights
       double precision, dimension(:), intent(in out) :: core_energies
       double precision, dimension(:, :), intent(in out) :: core_coefficients
       logical, intent(in) :: print_output

       double precision, dimension(2*(size(knots)-order-2), 2*(size(knots)-order-2)) :: hfmat
       double precision, dimension(2*(size(knots)-order), 2*(size(knots)-order-2)) :: coefmat
       double precision, dimension(size(core_coefficients, 1), size(core_coefficients, 2)) :: prev_core_coefs
       double precision, dimension(size(core_kappas), 2) :: epsbs
       double precision, dimension(2*(size(knots)-order-2)) :: spectrum
       double precision, dimension(2*(size(knots)-order)) :: coefb
       double precision, dimension(size(core_kappas)) :: relerrs
       double precision :: tolrel, epsb, abserr, tolerance, relerr
       integer, dimension(1) :: temp
       integer :: maxstep, ncores, nspectrum, step, i, kapb, nmnb, lb, nb, imax

       tolrel = 5.0d-10
       maxstep = 100
       ncores = size(core_kappas)
       nspectrum = 2 * (size(knots) - order - 2)

       do step = 1, 2
           prev_core_coefs = core_coefficients
           do i = 1, ncores
               kapb = core_kappas(i)
               nmnb = core_nmains(i)
               call hartree_fock_matrix_scf(kapb, core_kappas, prev_core_coefs, knots, order, nodes, weights, hfmat)
               call dirac(knots, order, kapb, potential_scf, nodes, weights, spectrum, coefmat, hfmat=hfmat)
               lb = (abs(2 * kapb + 1) - 1) / 2
               nb = nspectrum / 2 + nmnb - lb
               epsb = spectrum(nb) - ialpha**2
               coefb = coefmat(:, nb)
               core_coefficients(:, i) = coefb
               epsbs(i, step) = epsb
               if (print_output) write (*, *) step, nmnb, kapb, epsb
           end do
           if (print_output) write (*, *) ""
       end do

       do step = 3, maxstep
           do i = 1, ncores
               relerrs(i) = abs((epsbs(i, 1) - epsbs(i, 2)) / epsbs(i, 2))
           end do
           temp = maxloc(relerrs)
           imax = temp(1)
           abserr = abs(epsbs(imax, 2) - epsbs(imax, 1))
           tolerance = tolrel * abs(epsbs(imax, 2))
           if (abserr <= tolerance) exit

           prev_core_coefs = core_coefficients
           kapb = core_kappas(imax)
           nmnb = core_nmains(imax)
           call hartree_fock_matrix_scf(kapb, core_kappas, prev_core_coefs, knots, order, nodes, weights, hfmat)
           call dirac(knots, order, kapb, potential_scf, nodes, weights, spectrum, coefmat, hfmat=hfmat)
           lb = (abs(2 * kapb + 1) - 1) / 2
           nb = nspectrum / 2 + nmnb - lb
           epsb = spectrum(nb) - ialpha**2
           coefb = coefmat(:, nb)
           core_coefficients(:, imax) = coefb
           epsbs(imax, 1) = epsbs(imax, 2)
           epsbs(imax, 2) = epsb
           relerr = abs((epsbs(imax, 1) - epsbs(imax, 2)) / epsbs(imax, 2))
           if (print_output) write (*, *) step, nmnb, kapb, epsb, relerr
       end do
       if (print_output) write (*, *) ""

       core_energies = epsbs(:, 2)
   end subroutine refining_solutions_scf

   !> Compute the total SCF core energy including the Hartree-Fock
   !> double-counting correction.
   double precision function total_core_energy_scf(knots, order, core_kappas, core_energies, core_coefficients, &
                                                   nodes, weights)
       double precision, dimension(:), intent(in) :: knots
       integer, intent(in) :: order
       integer, dimension(:), intent(in) :: core_kappas
       double precision, dimension(:), intent(in) :: core_energies
       double precision, dimension(:, :), intent(in) :: core_coefficients
       double precision, dimension(:), intent(in) :: nodes, weights

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
           call hartree_fock_matrix_scf(kapa, core_kappas, core_coefficients, knots, order, nodes, weights, hfmat)
           vhf_aa = 0.0d0
           do k = 1, nspectrum
               do j = 1, nspectrum
                   vhf_aa = vhf_aa + ca(j) * hfmat(j, k) * ca(k)
               end do
           enddo
           etot = etot + 2 * abs(kapa) * (epsa - 0.5d0 * vhf_aa)
       end do
       total_core_energy_scf = etot
   end function total_core_energy_scf

   !> Solve the valence orbitals in the frozen-core Hartree-Fock field.
   subroutine valence_solutions_scf(knots, order, potential_scf, core_kappas, core_coefficients, valence_kappas, &
                                    valence_nmains, nodes, weights, valence_energies, valence_coefficients)
       double precision, dimension(:), intent(in) :: knots
       integer, intent(in) :: order
       type(fermi_potential), intent(in) :: potential_scf
       integer, dimension(:), intent(in) :: core_kappas
       double precision, dimension(:, :), intent(in) :: core_coefficients
       integer, dimension(:), intent(in) :: valence_kappas, valence_nmains
       double precision, dimension(:), intent(in) :: nodes, weights
       double precision, dimension(:), intent(out) :: valence_energies
       double precision, dimension(:, :), intent(out) :: valence_coefficients

       double precision, dimension(2*(size(knots)-order-2), 2*(size(knots)-order-2)) :: hfmat
       double precision, dimension(2*(size(knots)-order), 2*(size(knots)-order-2)) :: coefmat
       double precision, dimension(2*(size(knots)-order-2)) :: spectrum
       integer :: nvalences, nspectrum, i, kapa, nmna, la, na

       if (size(valence_kappas) == 0) then
           return
       endif

       nvalences = size(valence_kappas)
       nspectrum = 2 * (size(knots) - order - 2)

       do i = 1, nvalences
           kapa = valence_kappas(i)
           nmna = valence_nmains(i)
           call hartree_fock_matrix_scf(kapa, core_kappas, core_coefficients, knots, order, nodes, weights, hfmat)
           call dirac(knots, order, kapa, potential_scf, nodes, weights, spectrum, coefmat, hfmat=hfmat)
           la = (abs(2 * kapa + 1) - 1) / 2
           na = nspectrum / 2 + nmna - la
           valence_energies(i) = spectrum(na) - ialpha**2
           valence_coefficients(:, i) = coefmat(:, na)
       end do
   end subroutine valence_solutions_scf

   !> Assemble the Hartree-Fock matrix for one target `kappa` channel from
   !> the current set of occupied core orbitals.
   subroutine hartree_fock_matrix_scf(kappa, core_kappas, core_coefficients, knots, order, nodes, weights, hfmat)
       integer, intent(in) :: kappa
       integer, dimension(:), intent(in) :: core_kappas
       double precision, dimension(:, :), intent(in) :: core_coefficients
       double precision, dimension(:), intent(in) :: knots
       integer, intent(in) :: order
       double precision, dimension(:), intent(in) :: nodes, weights
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
   end subroutine hartree_fock_matrix_scf

end program print_polynomial_table_relativistic
