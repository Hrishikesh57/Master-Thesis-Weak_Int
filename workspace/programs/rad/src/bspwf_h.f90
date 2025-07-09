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

   ! Other declarations
   integer :: znucl, order, kappa, inbv, nmain, l
   integer :: i, unit_num, iostat
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

   !-----------------------------------------------------------------------
   print *, 'Enter the value for znucl (atomic number Z):'
   read(*,*) znucl

   print *, 'Enter the state (choose among: 1s, 2s, 2p):'
   read(*,*) state_str
   if (state_str .eq. '1s') then
      kappa = -1
      nmain = 1
      l = 0
   elseif (state_str .eq. '2s') then
      kappa = -1
      nmain = 2
      l = 0
   elseif (state_str .eq. '2p') then
      kappa = 1
      nmain = 2
      l = 1
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
   na = nspectrum/2 + (nmain - l)
   allocate(spectrum(nspectrum))
   allocate(coefmat(2 * nbspl, nspectrum))

   ! Define the nuclear potential.
   potential%zalpha = zalpha
   potential%rad    = rnucl

   ! Solve the Dirac equation using the B-spline basis.
   call solve(knots, order, kappa, potential, nodes, weights, spectrum, coefmat)

   ! Use the selected eigenstate's coefficients.
   coefs = coefmat(:, na)
   if (state_str .eq. '1s') then
      coefs = -coefs
   endif

   ! Write the computed wavefunctions to an output file.
   filename = 'polynomial_table_' // trim(state_str) // '.txt'
   unit_num = 10
   open(unit=unit_num, file=filename, status='replace', action='write', iostat=iostat)
   if (iostat /= 0) then
      print *, 'Error opening file: ', trim(filename)
      stop
   else
      print *, 'File opened successfully: ', trim(filename)
   endif

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
      write(unit_num, '(F12.6, 3X, ES18.10E3, 3X, ES18.10E3)') x, gx, fx
      print '(F12.6, 3X, ES18.10E3, 3X, ES18.10E3)', x, gx, fx
   end do

   close(unit_num)
   deallocate(knots, coefs, nodes, weights, spectrum, coefmat)

end program print_polynomial_table_relativistic