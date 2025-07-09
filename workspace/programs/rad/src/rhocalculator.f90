program calculate_rho0
    implicit none
    
    ! Physical constants
    double precision, parameter :: pi = 4.0d0 * atan(1.0d0)
    double precision, parameter :: hbar_c = 197.3269804d0  ! MeV·fm (or natural units)
    double precision, parameter :: sqrt2 = sqrt(2.0d0)
    
    ! Input variables
    double precision :: R_bohr, Cv_e, Cv_p, Z, N, G_fermi
    double precision :: rho0, R_fm
    
    ! Intermediate calculations
    double precision :: hbar_c_over_R, term1, term2, term3
    
    ! Conversion factor from Bohr to fm
    double precision, parameter :: bohr_to_fm = 0.529177249d0  ! 1 Bohr = 0.529 fm
    
    print *, '==============================================='
    print *, '        Calculate rho0 using the formula:'
    print *, '  rho0 = (3/4π) * (ℏc/R)³ * (Cv^p*Z - N) * (G/(2√2))'
    print *, '==============================================='
    print *, ''
    
    ! Get input values from user
    print *, 'Enter the nuclear radius R in Bohr radii:'
    print *, '(Default example: R = 1.1069e-4 Bohr)'
    read *, R_bohr
    
    print *, 'Enter Cv^e (electron coupling constant):'
    print *, '(Default example: Cv^e = 0.08)'
    read *, Cv_e
    
    print *, 'Enter Cv^p (proton coupling constant):'
    print *, '(Default example: Cv^p = 0.04 or similar)'
    read *, Cv_p
    
    print *, 'Enter atomic number Z:'
    print *, '(Example: Z = 6 for Carbon, Z = 92 for Uranium)'
    read *, Z
    
    print *, 'Enter number of neutrons N:'
    print *, '(Example: N = 6 for Carbon-12, N = 146 for Uranium-238)'
    read *, N
    
    print *, 'Enter Fermi coupling constant G in appropriate units:'
    print *, '(Example: G = 1.166e-5 GeV^-2 = 1.166e-11 MeV^-2)'
    read *, G_fermi
    
    print *, ''
    print *, '==============================================='
    print *, '                CALCULATION'
    print *, '==============================================='
    
    ! Convert R from Bohr to fm for display
    R_fm = R_bohr * bohr_to_fm
    
    ! Calculate step by step
    print *, 'Input values:'
    print *, '  R =', R_bohr, 'Bohr =', R_fm, 'fm'
    print *, '  Cv^e =', Cv_e
    print *, '  Cv^p =', Cv_p
    print *, '  Z =', Z
    print *, '  N =', N
    print *, '  G =', G_fermi, 'MeV^-2'
    print *, ''
    
    ! Step 1: Calculate (ℏc/R)³
    hbar_c_over_R = hbar_c / R_fm  ! Using R in fm
    term1 = (3.0d0 / (4.0d0 * pi)) * (hbar_c_over_R**3)
    
    print *, 'Step 1: (3/4π) * (ℏc/R)³'
    print *, '  ℏc =', hbar_c, 'MeV·fm'
    print *, '  ℏc/R =', hbar_c_over_R, 'MeV'
    print *, '  (ℏc/R)³ =', hbar_c_over_R**3, 'MeV³'
    print *, '  (3/4π) * (ℏc/R)³ =', term1, 'MeV³'
    print *, ''
    
    ! Step 2: Calculate (Cv^p*Z - N)
    term2 = Cv_p * Z - N
    
    print *, 'Step 2: (Cv^p*Z - N)'
    print *, '  Cv^p * Z =', Cv_p * Z
    print *, '  Cv^p*Z - N =', term2
    print *, ''
    
    ! Step 3: Calculate G/(2√2)
    term3 = G_fermi / (2.0d0 * sqrt2)
    
    print *, 'Step 3: G/(2√2)'
    print *, '  2√2 =', 2.0d0 * sqrt2
    print *, '  G/(2√2) =', term3, 'MeV^-2'
    print *, ''
    
    ! Final calculation
    rho0 = term1 * term2 * term3
    
    print *, '==============================================='
    print *, '                 RESULT'
    print *, '==============================================='
    print *, 'rho0 =', rho0, 'MeV⁵'
    print *, ''
    
    ! Additional useful information
    print *, 'Additional information:'
    print *, '  Units: [rho0] = MeV⁵'
    print *, '  Sign: rho0 is', merge('positive', 'negative', rho0 > 0.0d0)
    print *, '  Magnitude: |rho0| =', abs(rho0), 'MeV⁵'
    print *, '  Scientific notation: rho0 =', rho0, 'MeV⁵'
    
    ! Show the complete formula with substituted values
    print *, ''
    print *, 'Complete formula with your values:'
    print *, 'rho0 = (3/4π) × (', hbar_c, '/', R_fm, ')³ × (', Cv_p, '×', Z, ' - ', N, ') × (', G_fermi, '/(2√2))'
    print *, '     = ', term1, ' × ', term2, ' × ', term3
    print *, '     = ', rho0, ' MeV⁵'
    
end program calculate_rho0