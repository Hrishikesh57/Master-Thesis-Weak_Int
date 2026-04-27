module bspqed_physconsts
    implicit none
    save

    !> The value of the inverse fine-structure constant taken from NIST 2014
    !> CODATA: 137.035 999 139(31)
    double precision, parameter :: ialpha = 137.035999139d0

    !> The value of the fine-structure constant given as : 1.0d0 / ialpha
    double precision, parameter :: alpha = 1.0d0 / ialpha

    !> The value of the reduced Planck constant times the speed of light in
    !> vacuum in MeV fm taken from NIST 2014 CODATA: 197.326 9788(12) MeV fm
    double precision, parameter :: hbarc = 197.3269788d0

    !> The value of the electron's rest energy in MeV taken from NIST 2014
    !> CODATA: 0.510 998 9461(31) MeV
    double precision, parameter :: electron_rest_energy = 0.5109989461d0

    !> The reduced Compton wavelength of the electron
    double precision, parameter :: electron_reduced_compton = hbarc / electron_rest_energy



    ! ! Add definitions of additional constants
        !> Conversion factors
    double precision, parameter :: eV_to_Hartree = 1.0d0 / 27.21138602d0  ! 1 eV = 1/27.21138602 Hartree
    double precision, parameter :: MeV_to_Hartree = eV_to_Hartree * 1.0d6  ! 1 MeV = 1e6 eV
    double precision, parameter :: GeV_to_Hartree = eV_to_Hartree * 1.0d9  ! 1 GeV = 1e9 eV
    double precision, parameter :: fm_to_Bohr = 1.0d-15 / 0.52917721067d-10  ! Convert fm to Bohr
    double precision, parameter :: fm_to_Bohr_factor = 1.8897259886d-5  ! More precise value

    !> Electron rest energy in atomic units 
    double precision, parameter :: electron_rest_energy_au = 137.0d0**2
    double precision, parameter :: electron_rest_energy_GeV = 0.0005109989461d0

    !> hbarc in atomic units (Hartree·Bohr)
    !> hbarc_au = hbarc (MeV·fm) * MeV_to_Hartree (Hartree/MeV) * fm_to_Bohr (Bohr/fm)
    double precision, parameter :: hbarc_au = 137.035999139d0  ! Hartree·Bohr
    double precision, parameter :: hbarc_rel = 1 


 !> The reduced Compton wavelength of the electron in atomic units 
    double precision, parameter :: electron_reduced_compton_au = hbarc_au / electron_rest_energy_au

    

    !> Fermi coupling constant G_F in GeV^{-2}
    double precision, parameter :: Gfermi = 1.1663787d-5  ! GeV^{-2}
    double precision, parameter :: Gfermi_MeV = 1.1663787d-2  ! MeV^{-2}
    double precision , parameter :: Gfermi_relativistic = Gfermi * electron_rest_energy_GeV *electron_rest_energy_GeV
    !> G_Fermi in atomic units (Hartree^{-2})
    !> Gfermi_au = Gfermi / (GeV_to_Hartree)^2
    double precision, parameter :: Gfermi_au = Gfermi / (GeV_to_Hartree**2 ) ! Hartree^{-2}

    !> sin^2(theta_W) (Weinberg angle squared)
    double precision, parameter :: sin2thetaW = 0.2312d0

    !> Internal common-radius weak-source convention:
    !> C_V^e = 1 - 4 sin^2(theta_W).
    double precision, parameter :: Cv_e = 1.0d0 - 4.0d0 * sin2thetaW  ! Approximately 

    !> Internal common-radius weak-source convention:
    !> C_V^p = 1 - 4 sin^2(theta_W).
    double precision, parameter :: Cv_p = 1.0d0 - 4.0d0 * sin2thetaW   ! Approximately 

end module bspqed_physconsts
