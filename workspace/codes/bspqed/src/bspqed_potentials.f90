module bspqed_potentials
    use iso_fortran_env, only: error_unit  ! Use ISO Fortran environment for error handling
    use bspqed_function, only: function_type ! Use the abstract function type
    use bspqed_physconsts,  ! Use physical constants
    use bspqed_uehling, only: uehling_sphere_potential
    implicit none
    public

    public :: free_potential
    public :: sphere_potential
    public :: shell_potential
    public :: gauss_potential
    public :: fermi_potential
    public :: linear_potential
    public :: sphere_uehlsph_potential
    public :: weak_potential 

    type, extends(function_type) :: free_potential
    contains
        procedure, pass(self) :: eval => free_potential_eval
    end type free_potential

    type, extends(function_type) :: sphere_potential
        double precision :: zalpha
        double precision :: rad
    contains
        procedure, pass(self) :: eval => sphere_potential_eval
    end type sphere_potential

    type, extends(function_type) :: shell_potential
        double precision :: zalpha
        double precision :: rad
    contains
        procedure, pass(self) :: eval => shell_potential_eval
    end type shell_potential

    type, extends(function_type) :: gauss_potential
        double precision :: zalpha
        double precision :: rms
    contains
        procedure, pass(self) :: eval => gauss_potential_eval
    end type gauss_potential

    type, extends(function_type) :: fermi_potential
        double precision :: zalpha
        double precision :: a
        double precision :: c
    contains
        procedure, pass(self) :: eval => fermi_potential_eval
    end type fermi_potential

    type, extends(function_type) :: linear_potential
    contains
        procedure, pass(self) :: eval => linear_potential_eval
    end type linear_potential

    type, extends(function_type) :: sphere_uehlsph_potential
        double precision :: alpha
        double precision :: zalpha
        double precision :: rad
    contains
        procedure, pass(self) :: eval => sphere_uehlsph_potential_eval
    end type sphere_uehlsph_potential

    double precision, parameter :: pi = 4.0d0 * atan(1.0d0)


    type, extends(function_type) :: weak_potential
    double precision :: Z            ! Atomic number
    double precision :: N            ! Number of neutrons
    double precision :: R_au         ! Nuclear radius R in atomic units
    double precision :: rho0         ! Computed constant rho0 in atomic units
    double precision :: Cv_pZ_minus_N
    double precision :: c  ! New parameter
contains
    procedure, pass(self) :: init => weak_potential_init
    procedure, pass(self) :: eval => weak_potential_eval
end type weak_potential

contains

subroutine weak_potential_init(self, Z_in, N_in, R_au_in, c_in)
    class(weak_potential), intent(inout) :: self
    double precision, intent(in) :: Z_in, N_in, R_au_in, c_in

    double precision :: Cv_pZ_minus_N
    double precision :: rho0_au

    ! Assign Z, N, R_au, and c
    self%Z = Z_in
    self%N = N_in
    self%R_au = R_au_in
    self%c = c_in

    ! Compute (Cv_p * Z - N)
    Cv_pZ_minus_N = (Cv_p * self%Z) - self%N
    self%Cv_pZ_minus_N = Cv_pZ_minus_N

    ! Compute rho0 in atomic units
    rho0_au = (3.0d0 / (4.0d0 * pi)) * ((self%c / self%R_au)**3) * Cv_pZ_minus_N * (Gfermi_relativistic / (2.0d0 * sqrt(2.0d0)))

    ! Assign rho0
    self%rho0 = rho0_au

end subroutine weak_potential_init

    double precision function weak_potential_eval(x, self)
        double precision, intent(in) :: x
        class(weak_potential), intent(in) :: self

        if (x < 0.0d0) then
            write (error_unit, *) "* ERROR in eval of weak_potential"
            write (error_unit, *) "* x must be nonnegative."
            stop
        end if

        if (self%R_au <= 0.0d0) then
            write (error_unit, *) "* ERROR in eval of weak_potential"
            write (error_unit, *) "* rad must be positive."
            stop
        end if

        if (x <= self%R_au) then
            ! Inside the nuclear radius:
            weak_potential_eval = Cv_e * self%rho0
        else
            ! Outside the nuclear radius:
            weak_potential_eval = 0.0d0
        end if

        return
    end function weak_potential_eval


    double precision function free_potential_eval(x, self)
        double precision, intent(in) :: x
        class(free_potential), intent(in) :: self

        if (x < 0.0d0) then
            write (error_unit, *) "* ERROR in eval of free_potential:"
            write (error_unit, *) "* x must be nonnegative."
            stop
        end if

        free_potential_eval = 0.0d0
        return
    end function free_potential_eval

    double precision function sphere_potential_eval(x, self)
        double precision, intent(in) :: x
        class(sphere_potential), intent(in) :: self

        double precision :: zalpha, rad

        zalpha = self%zalpha
        rad = self%rad

        if (x < 0.0d0) then
            write (error_unit, *) "* ERROR in eval of sphere_potential"
            write (error_unit, *) "* x must be nonnegative."
            stop
        end if

        if (rad <= 0.0d0) then
            write (error_unit, *) "* ERROR in eval of sphere_potential"
            write (error_unit, *) "* rad must be positive."
            stop
        end if

        if (x < rad) then
            sphere_potential_eval = -zalpha * (3.0d0 - (x / rad)**2) / (2.0d0 * rad)
            return
        else
            sphere_potential_eval = -zalpha / x
            return
        end if
    end function sphere_potential_eval

    double precision function shell_potential_eval(x, self)
        double precision, intent(in) :: x
        class(shell_potential), intent(in) :: self

        double precision :: zalpha, rad

        zalpha = self%zalpha
        rad = self%rad

        if (x < 0.0d0) then
            write (error_unit, *) "* ERROR in eval of sphere_potential"
            write (error_unit, *) "* x must be nonnegative."
            stop
        end if

        if (rad <= 0.0d0) then
            write (error_unit, *) "* ERROR in eval of sphere_potential"
            write (error_unit, *) "* rad must be positive."
            stop
        end if

        if (x < rad) then
            shell_potential_eval = -zalpha / rad
            return
        else
            shell_potential_eval = -zalpha / x
            return
        end if
    end function shell_potential_eval

    double precision function gauss_potential_eval(x, self)
        double precision, intent(in) :: x
        class(gauss_potential), intent(in) :: self

        double precision :: zalpha, rms

        zalpha = self%zalpha
        rms = self%rms

        if (x < 0.0d0) then
            write (error_unit, *) "* ERROR in eval of gauss_potential"
            write (error_unit, *) "* x must be nonnegative."
            stop
        end if

        if (rms <= 0.0d0) then
            write (error_unit, *) "* ERROR in eval of gauss_potential"
            write (error_unit, *) "* rms must be positive."
            stop
        end if

        if (x > 0.0) then
            gauss_potential_eval = -zalpha * erf(sqrt(3.0d0 / 2.0d0) * (x / rms)) / x
            return
        else  ! i.e. x == 0.0d0
            gauss_potential_eval = -zalpha * sqrt(6.0d0 / pi) / rms
            return
        end if
    end function gauss_potential_eval

    double precision function linear_potential_eval(x, self)
        double precision, intent(in) :: x
        class(linear_potential), intent(in) :: self

        if (x < 0.0d0) then
            write (error_unit, *) "* ERROR in eval of linear_potential:"
            write (error_unit, *) "* x must be nonnegative."
            stop
        end if

        linear_potential_eval = x
        return
    end function linear_potential_eval

    double precision function fermi_potential_eval(x, self)
        double precision, intent(in) :: x
        class(fermi_potential), intent(in) :: self

        double precision :: zalpha, a, c, r, s3a, n, s2r, s3r, factor

        zalpha = self%zalpha
        a = self%a
        c = self%c

        if (x < 0.0d0) then
            write (error_unit, *) "* ERROR in eval of fermi_potential"
            write (error_unit, *) "* x must be nonnegative."
            stop
        end if

        if (a <= 0.0d0) then
            write (error_unit, *) "* ERROR in eval of fermi_potential"
            write (error_unit, *) "* a must be positive."
            stop
        end if

        if (c <= 0.0d0) then
            write (error_unit, *) "* ERROR in eval of fermi_potential"
            write (error_unit, *) "* c must be positive."
            stop
        end if

        if (x == 0.0d0) then
            fermi_potential_eval = 0.0d0
            return
        end if

        a = a / c
        r = x / c
        s3a = sfunc(3, -1.0d0 / a)
        n = 1.0d0 + a**2 * (pi**2 - 6.0d0 * a * s3a)
        if (r < 1.0d0) then
            s2r = sfunc(2, (r - 1.0d0) / a)
            s3r = sfunc(3, (r - 1.0d0) / a)
            factor = 6.0d0 * a**3 * (s3r - s3a)
            factor = factor + 0.5d0 * r * a**2 * (pi**2 - 6.0d0 * s2r)
            factor = factor + 0.5d0 * r * (3.0d0 - r**2)
        else
            s2r = sfunc(2, (1.0d0 - r) / a)
            s3r = sfunc(3, (1.0d0 - r) / a)
            factor = n + 3.0d0 * a**2 * (r * s2r + 2.0d0 * a * s3r)
        end if
        fermi_potential_eval = -zalpha * factor / (n * x)
        return
    end function fermi_potential_eval

    double precision function sphere_uehlsph_potential_eval(x, self)
        double precision, intent(in) :: x
        class(sphere_uehlsph_potential), intent(in) :: self

        double precision :: alpha, zalpha, rad, sphpot, uehlpot
        type(sphere_potential) :: sph_pot
        type(uehling_sphere_potential) :: uehl_sph_pot

        alpha = self%alpha
        zalpha = self%zalpha
        rad = self%rad

        sph_pot%zalpha = zalpha
        sph_pot%rad = rad

        uehl_sph_pot%rad = rad

        sphpot = sph_pot%eval(x)
        uehlpot = uehl_sph_pot%eval(x)
        uehlpot = (alpha / pi) * zalpha * uehlpot

        sphere_uehlsph_potential_eval = sphpot + uehlpot
        return
    end function sphere_uehlsph_potential_eval

    double precision function sfunc(l, x)
        integer, intent(in) :: l
        double precision, intent(in) :: x

        double precision, dimension(0:24), parameter :: weights = [&
            & 9.9999999999999988898d-01, -9.9999999999999988898d-01, &
            & 9.9999999999996203037d-01, -9.9999999999683419905d-01, &
            & 9.9999999985920762224d-01, -9.9999999613411694188d-01, &
            & 9.9999992840519291715d-01, -9.9999905165055258394d-01, &
            & 9.9999063480600425269d-01, -9.9992891127931748940d-01, &
            & 9.9957546287386789263d-01, -9.9796887921273413902d-01, &
            & 9.9210135627641904410d-01, -9.7473348838492690493d-01, &
            & 9.3283006109116761717d-01, -8.5017916311864960566d-01, &
            & 7.1687126316297511419d-01, -5.4150365145658008181d-01, &
            & 3.5444486563642518817d-01, -1.9426081291418739738d-01, &
            & 8.5828531071441754863d-02, -2.9156606763386567344d-02, &
            & 7.1108687661727044971d-03, -1.1042039205550148104d-03, &
            & 8.1792883003067385142d-05 &
        &]

        double precision :: y, s, a
        integer :: k

        y = exp(x)
        s = 0.0d0
        do k = 0, 24
            a = -y**k / (k + 1)**l
            s = s + weights(k) * a
        end do
        sfunc = s * y
        return
    end function sfunc



end module bspqed_potentials
