module bspwf_constants
    implicit none
    save

    !> The value of the inverse fine-structure constant taken from NIST 2018
    !> CODATA: 137.035 999 084(21)
    double precision, parameter :: ialpha = 137.035999084d0

    !> The value of the fine-structure constant given as : 1.0d0 / ialpha
    double precision, parameter :: alpha = 1.0d0 / ialpha

    !> The value of the reduced planck constant times the speed of light in
    !> vacuum in MeV fm taken from NIST 2018 CODATA: 197.326 980 4 MeV fm
    double precision, parameter :: hbarc = 197.3269804d0

    !> The value of the electron's rest energy in MeV taken from NIST 2018
    !> CODATA: 0.510 998 950 00(15) MeV
    double precision, parameter :: electron_rest_energy = 0.51099895000d0

    !> The reduced compton wavelength of the electron
    double precision, parameter :: reduced_compton_wavelength = hbarc / electron_rest_energy

    !> The Bohr radius of the electron
    double precision, parameter :: bohr_radius = reduced_compton_wavelength / alpha

    !> The Rydberg energy of the electron in eV
    double precision, parameter :: rydberg_energy = 0.5d0 * alpha**2 * electron_rest_energy * 1.0d+06

end module bspwf_constants
