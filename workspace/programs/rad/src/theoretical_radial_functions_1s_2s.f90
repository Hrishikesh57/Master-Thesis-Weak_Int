program compute_radial_function_1s_2s
   use bspqed_physconsts
   implicit none
   ! Variable declarations
   double precision :: Z, n, nr, kappa, alpha_val, gamma_val, E, lambda, N1
   double precision :: r, rmin, rmax, dr
   double precision :: func_g = 0.0d0, func_f = 0.0d0
   double precision :: mc2
   double precision :: prefactor, sqrt_term, power_term, exp_term, bracket_term
   double precision :: temp1, temp2
   double precision :: j  ! for the relativistic case
   integer :: i, npoints, a0
   character(len=2) :: state
   character(len=20) :: c_choice
   integer :: unit_number
   double precision, parameter :: zero_tol = 1.0d-300
   double precision, parameter :: pi_val = 3.14159265358979323846d0
   double precision :: c_val
   logical :: nonrel_flag
 
   !------------------------------------------------------------
   ! Get inputs from the user
   !------------------------------------------------------------
   print *, 'Enter atomic number Z (e.g., 1 for Hydrogen, 92 for Uranium):'
   read *, Z
 
   print *, 'Enter state (1s or 2s):'
   read *, state
 
   print *, 'Enter "rel" or "nonrel" for speed of light:'
   read *, c_choice
 
   if (trim(adjustl(c_choice)) == 'rel') then
      c_val = 1! Relativistic speed of light
      nonrel_flag = .false.
   else if (trim(adjustl(c_choice)) == 'nonrel') then
      c_val = 137036.0d0       ! Nonrelativistic option
      nonrel_flag = .true.
   else
      print *, 'Invalid choice for c. Defaulting to relativistic (c = 137).'
      c_val = 137.0d0
      nonrel_flag = .false.
   end if
 
   if (trim(state) == '1s') then
      n = 1.0d0
      kappa = -1.0d0
      nr = n - abs(kappa)
   else if (trim(state) == '2s') then
      n = 2.0d0
      kappa = -1.0d0
      nr = n - abs(kappa)
   else
      print *, 'Invalid state. Please enter 1s or 2s.'
      stop
   end if
 
   ! Compute mc^2 and the fine-structure constant
   mc2 = c_val**2
   alpha_val = alpha
 
   !------------------------------------------------------------
   ! Set parameters for nonrelativistic vs. relativistic case
   !------------------------------------------------------------
   if (nonrel_flag) then
      gamma_val = 1.0d0
      E = mc2
      ! In atomic units, a0 = 1.
      a0 = 1
      lambda = Z / (n*a0)

   else
      gamma_val = sqrt(kappa**2 - (Z * alpha_val)**2)
      if (gamma_val <= 0.0d0) then
         print *, 'Error: gamma_val is non-positive.'
         stop
      end if
      j = dabs(kappa) - 0.5d0
      E = mc2 * ( 1.d0 + (Z*alpha_val)**2 / &
      &           ((n - j - 0.5d0) + sqrt((j + 0.5d0)**2 - (Z*alpha_val)**2)) )**(-0.5d0) 
      lambda = sqrt((mc2**2-E**2))/c_val
   end if

 

 
   ! Set up the r grid
   rmin = 0.0001d0
   rmax = 0.1d0
   npoints = 1000
   dr = (rmax - rmin) / (npoints - 1)
 
   unit_number = 10
   open(unit=unit_number, file='radial_function.dat', status='unknown')
 
   write(unit=unit_number, fmt='(A)') 'r (a.u.)       r*g(r)         r*f(r)'
   print *, '   r (a.u.)       r*g(r)         r*f(r)'
 
   !------------------------------------------------------------
   ! Loop over r and compute the radial functions
   !------------------------------------------------------------
   do i = 0, npoints - 1
      r = rmin + i * dr
 
      if (trim(state) == '1s') then
         if (nonrel_flag) then
            !---------------------------------------------------------
            ! NONRELATIVISTIC 1s formula override:
            ! We now set:
            !   g_{1s}(r) = 2 * lambda^(3/2) * exp(-2 * lambda * r)
            !   f_{1s}(r) = 0.0
            !---------------------------------------------------------
            func_g = 2.0d0 * lambda**(1.5d0) * exp(-  lambda * r)
            func_f = 0.0d0
         else
            !---------------------------------------------------------
            ! RELATIVISTIC 1s formulas :
            prefactor = (2.0d0*lambda)**(1.5d0) / gamma_function(2.0d0*gamma_val + 1.0d0)
            sqrt_term = sqrt( ((E + mc2) * gamma_function(2.0d0*gamma_val + 1.0d0)) / &
                          ((gamma_val**2) * 4.0d0 * mc2**3 + 4.0d0 * gamma_val * E * mc2**2) )
            power_term = (2.0d0*lambda*r)**(gamma_val - 1.0d0)
            exp_term = exp(-lambda*r)
            bracket_term = gamma_val * mc2 + E 
            N1 = (4.0d0 * lambda**3) / gamma_function(2.0d0 * gamma_val + 1.0d0) * &
            ((E + gamma_val * mc2) / (gamma_val * mc2))
    
            func_g = prefactor * sqrt_term * power_term * exp_term * bracket_term
 
            sqrt_term = sqrt( ((-E + mc2) * gamma_function(2.0d0*gamma_val + 1.0d0)) / &
                          ((gamma_val**2) * 4.0d0 * mc2**3 + 4.0d0 * gamma_val * E * mc2**2) )
            func_f = - prefactor * sqrt_term * power_term * exp_term * bracket_term
         end if
 
      else if (trim(state) == '2s') then
         if (nonrel_flag) then
            !---------------------------------------------------------
            ! NONRELATIVISTIC 2s formulas:
            !   g2s(r) = (2λ)^(3/2)/Γ(2γ+1) *
            !            sqrt( ((mc²+E)*Γ(2γ+2))/(4mc²*((γ+1)mc²/E)*(((γ+1)mc²/E)+1)) )
            !            * (2λr)^(γ-1) e^(–λr) *
            !            { (( (γ+1)mc²/E + 1)*(1 - (2λr)/(2γ+1)) - 1) }
            !
            !   f2s(r) = - (2λ)^(3/2)/Γ(2γ+1) *
            !             sqrt( ((mc²-E)*Γ(2γ+2))/(4mc²*((γ+1)mc²/E)*(((γ+1)mc²/E)+1)) )
            !             * (2λr)^(γ-1) e^(–λr) *
            !             { (( (γ+1)mc²/E + 1)*(1 - (2λr)/(2γ+1)) + 1) }
            !---------------------------------------------------------
            prefactor = (2.0d0*lambda)**(1.5d0) / gamma_function(2.0d0*gamma_val + 1.0d0)
            sqrt_term = sqrt( (mc2 + E) * gamma_function(2.0d0*gamma_val + 2.0d0) / &
                          (4.0d0 * mc2 * (((gamma_val+1.0d0)*mc2/E)) * (((gamma_val+1.0d0)*mc2/E)+1.0d0) ) )
            power_term = (2.0d0*lambda*r)**(gamma_val - 1.0d0)
            exp_term = exp(-lambda*r)
            bracket_term = (((gamma_val+1.0d0)*mc2/E)+1.0d0) * (1.0d0 - (2.0d0*lambda*r)/(2.0d0*gamma_val+1.0d0)) - 1.0d0
            func_g = prefactor * sqrt_term * power_term * exp_term * bracket_term
 
            sqrt_term = sqrt( (mc2 - E) * gamma_function(2.0d0*gamma_val + 2.0d0) / &
                          (4.0d0 * mc2 * (((gamma_val+1.0d0)*mc2/E)) * (((gamma_val+1.0d0)*mc2/E)+1.0d0) ) )
            bracket_term = (((gamma_val+1.0d0)*mc2/E)+1.0d0) * (1.0d0 - (2.0d0*lambda*r)/(2.0d0*gamma_val+1.0d0)) + 1.0d0
            func_f = - prefactor * sqrt_term * power_term * exp_term * bracket_term
         else
            !---------------------------------------------------------
            ! RELATIVISTIC 2s formulas:
            prefactor = (2.0d0*lambda)**(1.5d0) / gamma_function(2.0d0*gamma_val + 1.0d0)
            sqrt_term = sqrt( ((mc2 + E) * gamma_function(2.0d0*gamma_val + 2.0d0)) / &
                          (4.0d0 * mc2 * (gamma_val+1.0d0) * (((gamma_val+1.0d0)*mc2/E + 1.0d0)) ) )
            power_term = (2.0d0*lambda*r)**(gamma_val - 1.0d0)
            exp_term = exp(-lambda*r)
            temp1 = ((gamma_val+1.0d0)*mc2/E) + 1.0d0
            temp2 = 1.0d0 - (2.0d0*lambda*r)/(2.0d0*gamma_val+1.0d0)
            bracket_term = temp1 * temp2 - 1.0d0
            func_g = prefactor * sqrt_term * power_term * exp_term * bracket_term
 
            sqrt_term = sqrt( ((mc2 - E) * gamma_function(2.0d0*gamma_val + 2.0d0)) / &
                          (4.0d0 * mc2 * (gamma_val+1.0d0) * (((gamma_val+1.0d0)*mc2/E + 1.0d0)) ) )
            bracket_term = temp1 * temp2 + 1.0d0
            func_f = - prefactor * sqrt_term * power_term * exp_term * bracket_term
         end if
      end if
 
      ! Write out the values: file and screen show r, r*g(r), r*f(r)
      write(unit_number, '(F10.5, 3X, ES18.10E3, 3X, ES18.10E3)') r, r*func_g, r*func_f
      print '(F10.5, 3X, ES18.10E3, 3X, ES18.10E3)', r, r*func_g, r*func_f
   end do
 
   close(unit_number)
   print *, "Calculation complete. Data written to radial_function.dat"
   print *, "E = " , E
   print *, "Lamda = ", lambda
   print * , "N1 = ", N1
 
 contains
 
   !------------------------------------------------------------
   ! Gamma function
   !------------------------------------------------------------
   double precision function gamma_function(x)
     implicit none
     double precision, intent(in) :: x
     gamma_function = gamma(x)
   end function gamma_function

 
 end program compute_radial_function_1s_2s
 