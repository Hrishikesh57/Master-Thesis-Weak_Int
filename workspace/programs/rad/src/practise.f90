program print_polynomial_table_relativistic
   use bspqed_dkb
   use bspqed_coefficients
   use bspqed_utilities, only: meanr, kntseq1
   use bspqed_potentials
   use bspqed_quadrature, only: nodes20, weights20
   use bspqed_physconsts        ! This module should define alpha and electron_rest_energy
   implicit none

   ! Variable declarations
   character(len=10) :: state_str
   character(len=5)  :: znucl_str ! Added to store Z as string
   character(len=100) :: filename_energy ! For energy output file

   ! Other declarations
   integer :: znucl, order, kappa, inbv, nmain, l
   integer :: i, unit_num, iostat, unit_num_energy ! Added unit for energy file
   integer :: npoints, nknots, nbspl, nintervals, ngl, nspectrum, na
   integer :: npoints_output
   logical :: file_exists ! Variable to check if file exists
   double precision :: gx, fx, x
   double precision, allocatable :: knots(:), coefs(:), nodes(:), weights(:)
   double precision :: rnucl, rcavity, zalpha
   double precision :: x_start, x_end, x_step
   character(len=100) :: filename
   double precision, allocatable :: spectrum(:), coefmat(:,:)
   type(sphere_potential) :: potential
   double precision :: rms_fm

   !-----------------------------------------------------------------------
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

end program print_polynomial_table_relativistic