program print_polynomial_table_relativistic
   use bspqed_dkb
   use bspqed_coefficients
   use bspqed_utilities, only: meanr, kntseq1
   use bspqed_potentials, only: sphere_potential
   use bspqed_quadrature, only: nodes20, weights20
   use bspqed_physconsts        ! This module should define alpha and electron_rest_energy
   use bspwf_constants, only: bohr_radius
   use bspwf_dirac, only: dirac
   use bspwf_function, only: function_type
   use bspwf_hartree_fock, only: dhfc_dir_matrix, dhfc_matrix
   use bspwf_knots, only: knots_01
   use bspwf_potentials, only: fermi_potential
   use bspwf_quadrature, only: nodes10, weights10
   implicit none

   ! Constants
   double precision, parameter :: pi_local = 4.0d0 * atan(1.0d0)
   double precision, parameter :: ialpha_scf = 1.0d0/137.036d0  ! Fine structure constant

   ! Variable declarations
   character(len=10) :: state_str
   character(len=5)  :: znucl_str ! Added to store Z as string

   ! Other declarations
   integer :: znucl, order, kappa, inbv, nmain, l
   integer :: i, unit_num, iostat ! Added unit for energy file
   integer :: npoints, nknots, nbspl, nintervals, ngl, nspectrum, na
   integer :: npoints_output
   double precision :: gx, fx, x
   double precision, allocatable :: knots(:), coefs(:), nodes(:), weights(:)
   double precision :: rnucl, rcavity, zalpha
   double precision :: x_start, x_end, x_step
   character(len=100) :: filename
   double precision, allocatable :: spectrum(:), coefmat(:,:)
   type(sphere_potential) :: potential
   double precision :: rms_fm

   ! New variables for multi-electron choice
   integer :: wavefunction_type, num_electrons

   ! Multi-electron SCF variables
   integer :: ncores, nvalences, nbsplines
   integer, dimension(:), allocatable :: core_kappas, core_nmains, valence_kappas, valence_nmains
   double precision, dimension(:, :), allocatable :: core_coefficients, valence_coefficients
   double precision, dimension(:), allocatable :: core_energies, valence_energies
   double precision :: rms, afermi, cfermi, ecore_total
   type(fermi_potential) :: nuclpot_scf
   real :: start, finish

   !-----------------------------------------------------------------------
   ! Ask user to choose between single-electron and multi-electron wavefunctions
   !-----------------------------------------------------------------------
   print *, 'Choose the type of wavefunction calculation:'
   print *, '1 - Single-electron wavefunctions'
   print *, '2 - Multi-electron wavefunctions'
   print *, 'Enter your choice (1 or 2):'
   read(*,*) wavefunction_type

   if (wavefunction_type /= 1 .and. wavefunction_type /= 2) then
      print *, 'Error: Invalid choice. Please enter 1 or 2.'
      stop
   endif

   if (wavefunction_type == 2) then
      print *, 'Enter the number of electrons in the system (1-90):'
      read(*,*) num_electrons
      
      if (num_electrons < 1 .or. num_electrons > 90) then
         print *, 'Error: Number of electrons must be between 1 and 90.'
         stop
      endif
      
      print *, 'Multi-electron calculation selected with', num_electrons, 'electrons.'
      
      ! Get atomic number and nuclear parameters
      print *, 'Enter the atomic number Z:'
      read(*,*) znucl
      
      print *, 'Enter the nuclear RMS radius (in fm):'
      read(*,*) rms_fm
      rms = rms_fm ! Convert for SCF calculation
      
      ! Get the state for wavefunction output
      print *, 'Enter the state for wavefunction output (1s, 2s, 2p, 2p3, 3s, 3p, 3p3, 4s):'
      read(*,*) state_str
      
      ! Parse the requested state
      call parse_state(state_str, kappa, nmain, l)
      
      ! Setup electron configuration based on number of electrons
      call setup_electron_configuration(num_electrons, nmain, kappa, core_kappas, core_nmains, &
                                        valence_kappas, valence_nmains, ncores, nvalences)
      
      ! Display configuration
      write (*, *) "# ", repeat("=", 86)
      write (*, '(A,I0)') "# Atomic number Z = ", znucl
      write (*, '(A,I0)') "# Number of electrons = ", num_electrons
      write (*, '(A,F8.4)') "# Nuclear RMS radius = ", rms, " fm"
      write (*, '(A,A)') "# Target state for output: ", trim(state_str)
      write (*, *) "# Core electron configuration:"
      do i = 1, ncores
          write (*, '(A,I0,A,I0,A,I0)') "#   ", i, ": n=", core_nmains(i), ", kappa=", core_kappas(i)
      end do
      if (nvalences > 0) then
          write (*, *) "# Valence electron configuration:"
          do i = 1, nvalences
              write (*, '(A,I0,A,I0,A,I0)') "#   ", i, ": n=", valence_nmains(i), ", kappa=", valence_kappas(i)
          end do
      endif
      write (*, *) "# ", repeat("-", 86)
      
      ! Setup B-spline grid and nuclear potential for SCF
      order = 9
      npoints = 1000
      nknots = npoints + 2 * (order - 1)
      nbsplines = nknots - order
      
      allocate(knots(nknots))
      knots = 0.0d0
      
      ! Nuclear potential parameters (Fermi distribution)
      afermi = 2.3d0 / (4.0d0 * log(3.0d0)) / bohr_radius
      cfermi = sqrt((5.0d0 * rms**2 - 7.0d0 * (pi_local * afermi)**2) / 3.0d0) / bohr_radius
      
      rcavity = 40.0d0
      call knots_01(cfermi, rcavity, npoints/10, order, knots)
      
      nuclpot_scf%znucl = znucl
      nuclpot_scf%a = afermi
      nuclpot_scf%c = cfermi
      
      ! Allocate arrays for SCF calculation
      allocate(core_energies(ncores))
      core_energies = 0.0d0
      allocate(core_coefficients(2*nbsplines, ncores))
      core_coefficients = 0.0d0
      
      ! Perform SCF calculation
      print *, "# Starting SCF calculation for core electrons..."
      call cpu_time(start)
      call starting_solutions_scf(knots, order, nuclpot_scf, core_kappas, core_nmains, nodes10, weights10, &
                                  core_energies, core_coefficients, .true.)
      call cpu_time(finish)
      
      write (*, *) "# ", repeat("-", 86)
      write (*, *) "# Starting approximations completed!"
      write (*, '(A,F8.3,A)') "# Duration: ", finish - start, " seconds"
      write (*, *) "# Core state energies:"
      do i = 1, ncores
          write (*, '(A,I0,A,I0,A,ES15.8)') "#   n=", core_nmains(i), ", kappa=", core_kappas(i), &
                                           ", E=", core_energies(i)
      end do
      write (*, *) "# ", repeat("-", 86)
      
      call cpu_time(start)
      call refining_solutions_scf(knots, order, nuclpot_scf, core_kappas, core_nmains, nodes10, weights10, &
                                  core_energies, core_coefficients, .true.)
      call cpu_time(finish)
      
      write (*, *) "# ", repeat("-", 86)
      write (*, *) "# SCF refinement completed!"
      write (*, '(A,F8.3,A)') "# Duration: ", finish - start, " seconds"
      write (*, *) "# Final core state energies:"
      do i = 1, ncores
          write (*, '(A,I0,A,I0,A,ES15.8)') "#   n=", core_nmains(i), ", kappa=", core_kappas(i), &
                                           ", E=", core_energies(i)
      end do
      
      ! Calculate total core energy
      ecore_total = total_core_energy_scf(knots, order, core_kappas, core_energies, core_coefficients, &
                                          nodes10, weights10)
      write (*, '(A,ES15.8)') "# Total core energy: ", ecore_total
      write (*, *) "# ", repeat("-", 86)
      
      ! Handle valence electrons if any
      if (nvalences > 0) then
          allocate(valence_energies(nvalences))
          valence_energies = 0.0d0
          allocate(valence_coefficients(2*nbsplines, nvalences))
          valence_coefficients = 0.0d0
          
          call valence_solutions_scf(knots, order, nuclpot_scf, core_kappas, core_coefficients, &
                                     valence_kappas, valence_nmains, nodes10, weights10, &
                                     valence_energies, valence_coefficients)
          
          write (*, *) "# Valence state energies:"
          do i = 1, nvalences
              write (*, '(A,I0,A,I0,A,ES15.8)') "#   n=", valence_nmains(i), ", kappa=", valence_kappas(i), &
                                               ", E=", valence_energies(i)
          end do
          write (*, *) "# ", repeat("-", 86)
      endif
      
      ! Extract target wavefunction and generate output
      allocate(coefs(2 * nbsplines))
      call extract_target_wavefunction_scf(kappa, nmain, core_kappas, core_nmains, core_coefficients, &
                                           valence_kappas, valence_nmains, valence_coefficients, &
                                           ncores, nvalences, coefs)
      
      ! Apply sign correction for multi-electron case only for 1s state
      if (trim(state_str) .eq. '1s') then
          coefs = -coefs
      endif
      
      ! Create output filename
      write(znucl_str, '(I0)') znucl
      filename = trim(znucl_str) // '_scf_polynomial_table_' // trim(state_str) // '.txt'
      unit_num = 10
      open(unit=unit_num, file=filename, status='replace', action='write', iostat=iostat)
      if (iostat /= 0) then
          print *, 'Error opening file: ', trim(filename)
          stop
      else
          print *, 'File opened successfully: ', trim(filename)
      endif
      
      ! Write header
      write(unit_num, '(A)') 'r (natural units, 1/MeV)    g_nk(r)           f_nk(r)'
      write(unit_num, '(A,A)') '# State: ', trim(state_str)
      write(unit_num, '(A,I0)') '# Z = ', znucl
      write(unit_num, '(A,I0)') '# Number of electrons = ', num_electrons
      write(unit_num, '(A,F8.4)') '# RMS radius = ', rms_fm
      write(unit_num, '(A)') '#'
      
      ! Generate wavefunction data
      npoints_output = 1000
      x_start = 1.0d-6
      x_end = 20.0d0
      x_step = (x_end - x_start) / (npoints_output - 1)
      
      do i = 1, npoints_output
          x = x_start + (i - 1) * x_step
          inbv = 1
          gx = radialg(knots, order, kappa, coefs, x, inbv)
          fx = radialf(knots, order, kappa, coefs, x, inbv)
          write(unit_num, '(F12.6, 3X, ES18.10E3, 3X, ES18.10E3)') x, gx, fx
      end do
      
      close(unit_num)
      
      write (*, *) "# ", repeat("=", 86)
      write (*, *) "# Multi-electron calculation completed successfully!"
      write (*, '(A,A)') "# Wavefunction written to: ", trim(filename)
      write (*, *) "# ", repeat("=", 86)
      
      ! Cleanup
      deallocate(knots, coefs, core_energies, core_coefficients)
      if (allocated(valence_energies)) deallocate(valence_energies, valence_coefficients)
      deallocate(core_kappas, core_nmains)
      if (allocated(valence_kappas)) deallocate(valence_kappas, valence_nmains)
      
      return
   endif

   !-----------------------------------------------------------------------
   ! Single-electron calculation (existing code)
   !-----------------------------------------------------------------------
   print *, 'Single-electron calculation selected.'
   
   print *, 'Enter the value for znucl (atomic number Z):'
   read(*,*) znucl

   print *, 'Enter the state (choose among: 1s, 2s, 2p, 3p, 4p, 5p, 6p, 7p, 8p, 9p, 10p):'
   read(*,*) state_str
   if (state_str .eq. '1s') then
      kappa = -1
      nmain = 1
      l = 0
   elseif (state_str .eq. '2s') then
      kappa = -1
      nmain = 2
      l = 0
   elseif (state_str .eq. '3s') then
      kappa = -1
      nmain = 3
      l =  0
   elseif (state_str .eq. '4s') then
      kappa = -1
      nmain = 4
      l =  0
   elseif (state_str .eq. '2p') then
      kappa = 1
      nmain = 2
      l = 1
   elseif (state_str .eq. '2p3') then 
      kappa = -2
      nmain = 2
      l = 1
   elseif (state_str .eq. '3p') then
      kappa = 1 
      nmain = 3
      l = 1 
   elseif (state_str .eq. '4p') then
      kappa = 1
      nmain = 4
      l = 1
   elseif (state_str .eq. '5p') then
      kappa = 1
      nmain = 5 
      l = 1 
   elseif (state_str .eq. '6p') then
      kappa = 1
      nmain = 6
      l = 1
   elseif (state_str .eq. '7p') then
      kappa = 1
      nmain = 7
      l = 1
   elseif (state_str .eq. '8p') then
      kappa = 1
      nmain = 8
      l = 1
   elseif (state_str .eq. '9p') then
      kappa = 1
      nmain = 9
      l = 1
   elseif (state_str .eq. '10p') then
      kappa = 1
      nmain = 10
      l = 1
   elseif(state_str .eq. '11p') then 
      kappa = 1
      nmain = 11
      l = 1
   elseif(state_str .eq. '12p') then 
      kappa = 1
      nmain = 12
      l = 1
   elseif(state_str .eq. '13p') then 
      kappa = 1
      nmain = 13
      l = 1
   elseif(state_str .eq. '14p') then 
      kappa = 1
      nmain = 14
      l = 1
   elseif(state_str .eq. '15p') then 
      kappa = 1
      nmain = 15
      l = 1
   elseif(state_str .eq. '16p') then 
      kappa = 1
      nmain = 16
      l = 1
   elseif(state_str .eq. '17p') then 
      kappa = 1
      nmain = 17
      l = 1
   elseif(state_str .eq. '18p') then 
      kappa = 1
      nmain = 18
      l = 1
   elseif(state_str .eq. '19p') then 
      kappa = 1
      nmain = 19
      l = 1
   elseif(state_str .eq. '20p') then 
      kappa = 1
      nmain = 20
      l = 1
   elseif(state_str .eq. '5f') then 
      kappa  = 3
      nmain  = 5
      l = 3
   else
      print *, 'Invalid state entered.'
      stop
   endif

   print *, 'Enter the value for rms radius (in fm):'
   read(*,*) rms_fm

   if (znucl <= 0) then
      print *, 'Error: znucl must be positive.'
      stop
   endif

   ! Initialize other parameters
   order  = 9
   ! zalpha = Z * alpha (dimensionless)
   zalpha = znucl * alpha

   ! Define the nuclear radius for a uniform sphere:
   !    rnucl = sqrt(5/3) * rms,
   ! expressed now in natural length units.
   rnucl = sqrt(5.d0/3.d0) * rms_fm / electron_reduced_compton

   ! Choose a cavity (box) radius in natural units.
   rcavity = 2.5d0 / alpha

   ! B-spline grid:
   npoints = 1000
   nknots  = npoints + 2*(order - 1)
   allocate(knots(nknots))
   call kntseq1(rnucl, rcavity, npoints/10, order, knots)
   nbspl      = nknots - order
   nintervals = npoints - 1
   if (knots(1) <= 0.0d0) knots(1) = 1.0d-6

   ! Allocate arrays for coefficients and quadrature.
   allocate(coefs(2 * nbspl))
   ngl = 20
   allocate(nodes(ngl), weights(ngl))
   nodes   = nodes20
   weights = weights20

   nspectrum = 2 * (nbspl - 2)
   print *, 'nspectrum = ', nspectrum
   ! ----------------------------------------------------------------------
   ! The array 'spectrum' now holds ALL positive‑energy eigenvalues in ascending order.
   ! We split them into two “sub‑blocks” (negative/positive κ) of size nspectrum/2 each.
   ! Within each block, the bound states appear in order of principal quantum number n.
   ! Index of our chosen level (nmain, l) in the 'spectrum' array:

   na = nspectrum/2 + (nmain - l)
   allocate(spectrum(nspectrum))
   allocate(coefmat(2 * nbspl, nspectrum))
   
   ! Define the nuclear potential.
   potential%zalpha = zalpha
   potential%rad    = rnucl

 ! ---------------------------------------------------------------------
   ! Solve the Dirac equation H ψ = E ψ in a B‑spline basis *including*
   ! finite‑size nuclear potential.  The subroutine “solve” does three things:
   !
   !  1) It constructs the Hamiltonian matrix
   !       H_{ij} = < B_i | α·p + β m + V_Coulomb(r;R) | B_j >
   !     and the overlap matrix
   !       S_{ij} = < B_i | B_j >.
   !
   !  2) It then calls a LAPACK generalized‑eigenvalue solver (e.g. DSYGV)
   !     to diagonalise the pair (H, S), solving
   !         H · c = E · S · c
   !     for all eigenvalues E and eigenvectors c.
   !
   !  3) On return, the array `spectrum(1:nspectrum)` contains the sorted
   !     positive‑energy eigenvalues E_n, and `coefmat(:,n)` holds the
   !     corresponding eigenvector c_n (the B‑spline coefficients).
   !
   ! Thus *this one call* performs the entire matrix diagonalisation for you.
   ! You can then pick off the unperturbed 2s and 2p₁/₂ energies simply by
   ! indexing into `spectrum(...)` as shown below.
   ! ---------------------------------------------------------------------
   call solve(knots, order, kappa, potential, nodes, weights, spectrum, coefmat)
   ! Because we use a *finite‑size* nuclear potential, the exact Dirac degeneracy
   ! between 2s and 2p1/2 is lifted: those two entries in 'spectrum' will differ slightly.
   ! -----------------------------------------------------------------------
   print *, 'Eigenvalue (energy) for chosen state = ', spectrum(na)

   ! Use the selected eigenstate's coefficients.
   coefs = coefmat(:, na)
   if (state_str .eq. '1s') then
      coefs = -coefs
   endif

   ! Write the computed wavefunctions to an output file.
   ! Convert znucl to string
   write(znucl_str, '(I0)') znucl
   !Construct the filename including znucl
   filename = trim(znucl_str) // '_polynomial_table_' // trim(state_str) // '.txt'
   unit_num = 10
   open(unit=unit_num, file=filename, status='replace', action='write', iostat=iostat)
   if (iostat /= 0) then
     print *, 'Error opening file: ', trim(filename)
     stop
   else
     print *, 'File opened successfully: ', trim(filename)
   endif

   !Comment out the table header writing
   write(unit_num, '(A)') 'r (natural units, 1/MeV)    g_nk(r)           f_nk(r)'
   print *, '   r (natural units, 1/MeV)    g_nk(r)           f_nk(r)'

   npoints_output = 1000
   x_start = 1.0d-6
   x_end   = 20.0d0
   x_step  = (x_end - x_start) / (npoints_output - 1)

   do i = 1, npoints_output
      x = x_start + (i - 1)*x_step
      inbv = 1
      gx = radialg(knots, order, kappa, coefs, x, inbv)
      fx = radialf(knots, order, kappa, coefs, x, inbv)
      !Comment out the table data writing
      write(unit_num, '(F12.6, 3X, ES18.10E3, 3X, ES18.10E3)') x, gx, fx
      print '(F12.6, 3X, ES18.10E3, 3X, ES18.10E3)', x, gx, fx
   end do

   close(unit_num)

   ! ! Write Z, RMS radius, and Energy to the state-specific file
   ! unit_num_energy = 11 ! Choose a different unit number
   ! filename_energy = 'Energy_states_' // trim(state_str) // '.txt'

   ! ! Check if the file exists
   ! inquire(file=filename_energy, exist=file_exists)

   ! if (file_exists) then
   !    ! File exists, open in append mode
   !    open(unit=unit_num_energy, file=filename_energy, position='append', action='write', iostat=iostat)
   ! else
   !    ! File does not exist, open as new and write header
   !    open(unit=unit_num_energy, file=filename_energy, status='new', action='write', iostat=iostat)
   !    if (iostat == 0) then
   !       write(unit_num_energy, '(A)') 'Z    RMS_Radius(fm)    Energy(natural units)' ! Write header only for new file
   !    endif
   ! endif

   ! if (iostat /= 0) then
   !    print *, 'Error opening energy file: ', trim(filename_energy)
   !    ! Decide whether to stop or continue without writing
   !    stop ! Stop if file cannot be opened/created
   ! else
   !    ! Write data: Z, R_fm, E_state
   !    write(unit_num_energy, '(I5, 3X, F12.6, 3X, ES18.10E3)') znucl, rms_fm, spectrum(na)
   !    close(unit_num_energy)
   !    print *, 'Energy data written to: ', trim(filename_energy)
   ! endif

   deallocate(knots, coefs, nodes, weights, spectrum, coefmat)

contains

   !-----------------------------------------------------------------------
   ! Parse State String
   !-----------------------------------------------------------------------
   subroutine parse_state(state_str, kappa, nmain, l)
       character(len=*), intent(in) :: state_str
       integer, intent(out) :: kappa, nmain, l
       
       select case (trim(state_str))
       case ('1s')
           kappa = -1; nmain = 1; l = 0
       case ('2s')
           kappa = -1; nmain = 2; l = 0
       case ('3s')
           kappa = -1; nmain = 3; l = 0
       case ('4s')
           kappa = -1; nmain = 4; l = 0
       case ('2p', '2p1')
           kappa = 1; nmain = 2; l = 1
       case ('2p3')
           kappa = -2; nmain = 2; l = 1
       case ('3p', '3p1')
           kappa = 1; nmain = 3; l = 1
       case ('3p3')
           kappa = -2; nmain = 3; l = 1
       case default
           print *, 'Invalid state: ', trim(state_str)
           print *, 'Valid states: 1s, 2s, 2p, 2p3, 3s, 3p, 3p3, 4s'
           stop
       end select
   end subroutine parse_state

   !-----------------------------------------------------------------------
   ! Setup Electron Configuration
   !-----------------------------------------------------------------------
   subroutine setup_electron_configuration(num_electrons, target_nmain, target_kappa, &
                                          core_kappas, core_nmains, &
                                          valence_kappas, valence_nmains, &
                                          ncores, nvalences)
       integer, intent(in) :: num_electrons, target_nmain, target_kappa
       integer, dimension(:), allocatable, intent(out) :: core_kappas, core_nmains
       integer, dimension(:), allocatable, intent(out) :: valence_kappas, valence_nmains
       integer, intent(out) :: ncores, nvalences
       
       integer :: electrons_placed, n, kap
       logical :: target_is_core
       
       ! Fill electron shells in order: 1s, 2s, 2p1/2, 2p3/2, 3s, 3p1/2, 3p3/2, etc.
       ! Determine configuration up to target state
       target_is_core = .false.
       electrons_placed = 0
       
       ! Count electrons up to target state
       do n = 1, target_nmain
           do kap = -n, n
               if (kap == 0) cycle
               electrons_placed = electrons_placed + 2 * abs(kap)
               if (n == target_nmain .and. kap == target_kappa) then
                   if (electrons_placed <= num_electrons) target_is_core = .true.
                   exit
               endif
               if (electrons_placed >= num_electrons) exit
           end do
           if (electrons_placed >= num_electrons) exit
       end do
       
       ! Build configuration
       ncores = 0
       nvalences = 0
       electrons_placed = 0
       
       ! Count states first
       do n = 1, max(target_nmain + 2, 4)
           do kap = -n, n
               if (kap == 0) cycle
               if (electrons_placed + 2 * abs(kap) <= num_electrons .or. &
                   (n == target_nmain .and. kap == target_kappa .and. target_is_core)) then
                   ncores = ncores + 1
                   electrons_placed = electrons_placed + 2 * abs(kap)
               else if (n == target_nmain .and. kap == target_kappa .and. .not. target_is_core) then
                   nvalences = nvalences + 1
               else if (electrons_placed < num_electrons .and. n <= target_nmain + 1) then
                   nvalences = nvalences + 1
               endif
               if (n == target_nmain .and. kap == target_kappa) exit
           end do
           if (n == target_nmain .and. kap == target_kappa) exit
       end do
       
       ! Allocate arrays
       allocate(core_kappas(ncores), core_nmains(ncores))
       if (nvalences > 0) then
           allocate(valence_kappas(nvalences), valence_nmains(nvalences))
       endif
       
       ! Fill arrays
       ncores = 0
       nvalences = 0
       electrons_placed = 0
       
       do n = 1, max(target_nmain + 2, 4)
           do kap = -n, n
               if (kap == 0) cycle
               if (electrons_placed + 2 * abs(kap) <= num_electrons .or. &
                   (n == target_nmain .and. kap == target_kappa .and. target_is_core)) then
                   ncores = ncores + 1
                   core_nmains(ncores) = n
                   core_kappas(ncores) = kap
                   electrons_placed = electrons_placed + 2 * abs(kap)
               else if (n == target_nmain .and. kap == target_kappa .and. .not. target_is_core) then
                   nvalences = nvalences + 1
                   valence_nmains(nvalences) = n
                   valence_kappas(nvalences) = kap
               else if (electrons_placed < num_electrons .and. n <= target_nmain + 1) then
                   nvalences = nvalences + 1
                   valence_nmains(nvalences) = n
                   valence_kappas(nvalences) = kap
               endif
               if (n == target_nmain .and. kap == target_kappa) exit
           end do
           if (n == target_nmain .and. kap == target_kappa) exit
       end do
   end subroutine setup_electron_configuration

   !-----------------------------------------------------------------------
   ! Extract Target Wavefunction from SCF results
   !-----------------------------------------------------------------------
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
       
       ! Check core states first
       do i = 1, ncores
           if (core_kappas(i) == target_kappa .and. core_nmains(i) == target_nmain) then
               coefs = core_coefficients(:, i)
               found = .true.
               print *, '# Target state found in core electrons'
               exit
           endif
       end do
       
       ! Check valence states if not found in core
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

   !-----------------------------------------------------------------------
   ! SCF subroutines (adapted from dhfc.f90)
   !-----------------------------------------------------------------------
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

       ! Initial step without electron-electron interactions
       step = 1
       do i = 1, ncores
           kapb = core_kappas(i)
           nmnb = core_nmains(i)
           hfdir(:, :, i) = 0.0d0
           call dirac(knots, order, kapb, potential_scf, nodes, weights, spectrum, coefmat, hfmat=hfdir(:, :, i))
           lb = (abs(2 * kapb + 1) - 1) / 2
           nb = nspectrum / 2 + nmnb - lb
           ! For the binding energy of the electron in the SCF
           epsb = spectrum(nb) - ialpha_scf**2
           coefb = coefmat(:, nb)
           core_coefficients(:, i) = coefb
           epsbs(i, 1) = epsb
           if (print_output) write (*, *) step, nmnb, kapb, epsb
       end do
       if (print_output) write (*, *) ""

       ! SCF iterations
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
               epsb = spectrum(nb) - ialpha_scf**2
               coefb = coefmat(:, nb)
               core_coefficients(:, i) = coefb
               epsbs(i, 2) = epsb
               if (print_output) write (*, *) step, nmnb, kapb, epsb
           end do
           if (print_output) write (*, *) ""

           ! Check convergence
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

       ! Initial steps
       do step = 1, 2
           prev_core_coefs = core_coefficients
           do i = 1, ncores
               kapb = core_kappas(i)
               nmnb = core_nmains(i)
               call hartree_fock_matrix_scf(kapb, core_kappas, prev_core_coefs, knots, order, nodes, weights, hfmat)
               call dirac(knots, order, kapb, potential_scf, nodes, weights, spectrum, coefmat, hfmat=hfmat)
               lb = (abs(2 * kapb + 1) - 1) / 2
               nb = nspectrum / 2 + nmnb - lb
               epsb = spectrum(nb) - ialpha_scf**2
               coefb = coefmat(:, nb)
               core_coefficients(:, i) = coefb
               epsbs(i, step) = epsb
               if (print_output) write (*, *) step, nmnb, kapb, epsb
           end do
           if (print_output) write (*, *) ""
       end do

       ! Refinement iterations
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
           epsb = spectrum(nb) - ialpha_scf**2
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

       ! Check if valence arrays have any elements
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
           valence_energies(i) = spectrum(na) - ialpha_scf**2
           valence_coefficients(:, i) = coefmat(:, na)
       end do
   end subroutine valence_solutions_scf

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