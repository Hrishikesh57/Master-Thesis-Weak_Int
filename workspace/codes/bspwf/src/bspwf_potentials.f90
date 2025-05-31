module bspwf_potentials
    use iso_fortran_env, only: error_unit
    use bspwf_function, only: function_type
    implicit none
    private

    public :: coulomb_potential
    public :: sphere_potential
    public :: fermi_potential

    type, extends(function_type) :: coulomb_potential
        double precision :: znucl
    contains
        procedure, pass(self) :: eval => coulomb_potential_eval
    end type coulomb_potential

    type, extends(function_type) :: sphere_potential
        double precision :: znucl
        double precision :: rad
    contains
        procedure, pass(self) :: eval => sphere_potential_eval
    end type sphere_potential

    type, extends(function_type) :: fermi_potential
        double precision :: znucl
        double precision :: a
        double precision :: c
    contains
        procedure, pass(self) :: eval => fermi_potential_eval
    end type fermi_potential

    double precision, parameter :: pi = 4.0d0 * atan(1.0d0)

contains

    double precision function coulomb_potential_eval(x, self)
        double precision, intent(in) :: x
        class(coulomb_potential), intent(in) :: self

        double precision :: znucl

        znucl = self%znucl

        if (x <= 0.0d0) then
            write (error_unit, *) "* ERROR in 'coulomb_potential_eval'"
            write (error_unit, *) "* 'x' must be positive."
            stop
        end if

        coulomb_potential_eval = -znucl / x
        return
    end function coulomb_potential_eval

    double precision function sphere_potential_eval(x, self)
        double precision, intent(in) :: x
        class(sphere_potential), intent(in) :: self

        double precision :: znucl, rad

        znucl = self%znucl
        rad = self%rad

        if (x < 0.0d0) then
            write (error_unit, *) "* ERROR in 'sphere_potential_eval'"
            write (error_unit, *) "* 'x' must be nonnegative."
            stop
        end if

        if (rad <= 0.0d0) then
            write (error_unit, *) "* ERROR in 'sphere_potential_eval'"
            write (error_unit, *) "* 'rad' must be positive."
            stop
        end if

        if (x < rad) then
            sphere_potential_eval = 0.5d0 * znucl * ((x / rad)**2 - 3.0d0) / rad
            return
        else
            sphere_potential_eval = -znucl / x
            return
        end if
    end function sphere_potential_eval

    double precision function fermi_potential_eval(x, self)
        double precision, intent(in) :: x
        class(fermi_potential), intent(in) :: self

        double precision :: znucl, a, c, a0, r0, s3a, n, s2r, s3r, factor

        znucl = self%znucl
        a = self%a
        c = self%c

        if (x < 0.0d0) then
            write (error_unit, *) "* ERROR in 'fermi_potential_eval'"
            write (error_unit, *) "* 'x' must be nonnegative."
            stop
        end if

        if (a <= 0.0d0) then
            write (error_unit, *) "* ERROR in 'fermi_potential_eval'"
            write (error_unit, *) "* 'a' must be positive."
            stop
        end if

        if (c <= 0.0d0) then
            write (error_unit, *) "* ERROR in 'fermi_potential_eval'"
            write (error_unit, *) "* 'c' must be positive."
            stop
        end if

        if (x == 0.0d0) then
            fermi_potential_eval = 0.0d0
            return
        end if
        a0 = a / c
        r0 = x / c
        s3a = sfunc(3, -1.0d0 / a0)
        n = 1.0d0 + a0**2 * (pi**2 - 6.0d0 * a0 * s3a)
        if (r0 < 1.0d0) then
            s2r = sfunc(2, (r0 - 1.0d0) / a0)
            s3r = sfunc(3, (r0 - 1.0d0) / a0)
            factor = 6.0d0 * a0**3 * (s3r - s3a)
            factor = factor + 0.5d0 * r0 * a0**2 * (pi**2 - 6.0d0 * s2r)
            factor = factor + 0.5d0 * r0 * (3.0d0 - r0**2)
        else
            s2r = sfunc(2, (1.0d0 - r0) / a0)
            s3r = sfunc(3, (1.0d0 - r0) / a0)
            factor = n + 3.0d0 * a0**2 * (r0 * s2r + 2.0d0 * a0 * s3r)
        end if
        fermi_potential_eval = -znucl * factor / (n * x)
        return
    end function fermi_potential_eval

    double precision function sfunc(l, x)
        integer, intent(in) :: l
        double precision, intent(in) :: x

        double precision, dimension(0:24), parameter :: weights = [&
        & 9.9999999999999988898d-01, -9.9999999999999988898d-01,&
        & 9.9999999999996203037d-01, -9.9999999999683419905d-01,&
        & 9.9999999985920762224d-01, -9.9999999613411694188d-01,&
        & 9.9999992840519291715d-01, -9.9999905165055258394d-01,&
        & 9.9999063480600425269d-01, -9.9992891127931748940d-01,&
        & 9.9957546287386789263d-01, -9.9796887921273413902d-01,&
        & 9.9210135627641904410d-01, -9.7473348838492690493d-01,&
        & 9.3283006109116761717d-01, -8.5017916311864960566d-01,&
        & 7.1687126316297511419d-01, -5.4150365145658008181d-01,&
        & 3.5444486563642518817d-01, -1.9426081291418739738d-01,&
        & 8.5828531071441754863d-02, -2.9156606763386567344d-02,&
        & 7.1108687661727044971d-03, -1.1042039205550148104d-03,&
        & 8.1792883003067385142d-05]
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

end module bspwf_potentials
