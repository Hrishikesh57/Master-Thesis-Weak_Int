module bspqed_uehling
    use iso_fortran_env, only: error_unit
    use bspqed_function, only: function_type
    implicit none
    private

    public :: uehling_coulomb_potential
    public :: uehling_sphere_potential

    double precision, parameter :: pi = 4.0d0 * atan(1.0d0)

    !+ ========================================================================
    !> @brief Abstract type for the Uehling-Coulomb potential energy function
    !>
    !> @author Halil Cakir
    !> @date 2017/12/19
    !>
    !> #### Description
    !>
    !> Definition of an abstract type for the potential energy function
    !> resulting from the first order vacuum polarization correction and one
    !> nuclear potential insertion. The nucleus is considered to be a point.
    !> The potential energy function is given by
    !>
    !>     V(x) = (alpha / pi) * (Z * alpha) * v(x),
    !>
    !> where alpha is the fine-structure constant and Z is the nuclear charge
    !> number. This function returns the value of v for positive values of the
    !> the argument x. The function v is evaluated using expressions from
    !> [1, 2].
    !>
    !> #### References
    !>
    !> [1] L. W. Fullerton, G. A. Rinkler, Jr., Phys. Rev. A 13, 1283 (1976)
    !> [2] S. Klarsfeld, Phys. Lett. 66B, 86 (1977)
    !>
    !+ #### Parameters
    !>
    !+ ========================================================================
    type, extends(function_type) :: uehling_coulomb_potential
    contains
        procedure, pass(self) :: eval => uehling_coulomb_potential_eval
    end type uehling_coulomb_potential

    !+ ========================================================================
    !> @brief Abstract type for the Uehling-sphere potential energy function
    !>
    !> @author Halil Cakir
    !> @date 2017/04/27
    !>
    !> #### Description
    !>
    !> Definition of an abstract type for the potential energy function
    !> resulting from the first order vacuum polarization correction and one
    !> nuclear potential insertion. The nucleus is considered to be a
    !> homogeneously charged sphere. The potential energy function is given by
    !>
    !>     V(x) = (alpha / pi) * (Z * alpha) * v(x),
    !>
    !> where alpha is the fine-structure constant and Z is the nuclear charge
    !> number. This function returns the value of v for non-negative values of
    !> the argument x. The function v is evaluated using expressions from [1].
    !>
    !> #### References
    !>
    !> [1] S. Klarsfeld, Phys. Lett. 66B, 86 (1977)
    !>
    !+ #### Parameters
    !>
    !> @param rad : Nuclear radius
    !+ ========================================================================
    type, extends(function_type) :: uehling_sphere_potential
        double precision :: rad
    contains
        procedure, pass(self) :: eval => uehling_sphere_potential_eval
    end type uehling_sphere_potential

contains

    !+ ========================================================================
    !+ @brief Evaluation of the Uehling-Coulomb potential energy function
    !+
    !+ @author Halil Cakir
    !+ @date 2017/12/19
    !+ ========================================================================
    double precision function uehling_coulomb_potential_eval(x, self)
        double precision, intent(in) :: x
        class(uehling_coulomb_potential), intent(in) :: self

        if (x <= 0.0d0) then
            write(error_unit, *) "* ERROR in eval of uehling_coulomb_potential"
            write(error_unit, *) "* x must be positive."
            stop
        endif

        uehling_coulomb_potential_eval = uehlclmbfunc(x) / x
        return
    end function uehling_coulomb_potential_eval

    !+ ========================================================================
    !+ @brief Evaluation of the Uehling-sphere potential energy function
    !+
    !+ @author Halil Cakir
    !+ @date 2017/04/27
    !+ ========================================================================
    double precision function uehling_sphere_potential_eval(x, self)
        double precision, intent(in) :: x
        class(uehling_sphere_potential), intent(in) :: self

        double precision :: rad

        rad = self%rad

        if (x < 0.0d0) then
            write(error_unit, *) "* ERROR in eval of uehling_sphere_potential"
            write(error_unit, *) "* x must be nonnegative."
            stop
        endif

        if (rad <= 0.0d0) then
            write(error_unit, *) "* ERROR in eval of uehling_sphere_potential"
            write(error_unit, *) "* rad must be positive."
            stop
        endif

        uehling_sphere_potential_eval = uehlsph(x, rad)
        return
    end function uehling_sphere_potential_eval

    double precision function uehlclmbfunc(x)
        double precision, intent(in) :: x
        double precision :: chi1
        if (x < 0.0d0) then
            chi1 = 0.0d0
        else if (x < 1.0d0) then
            chi1 = chip1(2.0d0 * x)
        else if (x < 100.0d0) then
            chi1 = chi(1, 2.0d0 * x)
        else
            chi1 = 0.0d0
        end if
        uehlclmbfunc = -2.0d0 * chi1 / 3.0d0
        return
    end function uehlclmbfunc

    double precision function chi(n, x)
        integer, intent(in) :: n
        double precision, intent(in) :: x
        double precision, parameter :: epsabs = 0.0d0
        double precision, parameter :: epsrel = 5.0d-14
        integer, parameter :: maxiter = 6
        double precision :: res, abserr

        call tsquad(f, maxiter, epsabs, epsrel, res, abserr)
        chi = 0.5d0**(n-1) * res
        return

    contains

        double precision function f(t)
            double precision, intent(in) :: t
            double precision :: a, b, c, d
            a = exp(-2.0d0 * x / t)
            b = t**(n-2)
            c = 1.0d0 + 0.125d0 * t**2
            d = 1.0d0 - 0.5d0 * t
            if (d < 0.0d0) d = 0.0d0
            d = sqrt(d * (1.0d0 + 0.5d0 * t))
            f = a * b * c * d
            return
        end function f
    end function chi

    !+ ========================================================================
    !+ Tanh-Sinh quadrature on [-1, 1] for a function f which can be expressed
    !+ as f(t) = g(1-t) for some other function g on [0, 2].
    !+ ========================================================================
    subroutine tsquad(f, maxiter, epsabs, epsrel, res, abserr)
        interface
            double precision function f(x)
                double precision, intent(in) :: x
            end function f
        end interface
        integer, intent(in) :: maxiter
        double precision, intent(in) :: epsabs, epsrel
        double precision, intent(out) :: res, abserr
        double precision, parameter :: tin = 1.0d-30
        double precision :: prv, cur, h, t, w, xpos, xneg, eps, bnd, cosht, &
        & exppt, expmt
        integer :: j, k
        k = 0
        h = 0.5d0**k
        j = 0
        cur = f(1.0d0) * 0.5d0 * pi
        eps = huge(1.0d0)
        do
            j = j + 1
            t = 0.5d0 * pi * sinh(j * h)
            cosht = cosh(t)
            exppt = exp(+t)
            expmt = exp(-t)
            w = 0.5d0 * pi * cosh(j * h) / cosht**2
            xpos = 1.0d0 / (exppt * cosht)
            xneg = 1.0d0 / (expmt * cosht)
            if (w < tin) exit
            cur = cur + (f(xpos) + f(xneg)) * w
        end do
        cur = h * cur
        do k = 1, maxiter
            prv = cur
            h = 0.5d0**k
            j = 1
            cur = 0.0d0
            do
                t = 0.5d0 * pi * sinh(j * h)
                cosht = cosh(t)
                exppt = exp(+t)
                expmt = exp(-t)
                w = 0.5d0 * pi * cosh(j * h) / cosht**2
                xpos = 1.0d0 / (exppt * cosht)
                xneg = 1.0d0 / (expmt * cosht)
                if (w < tin) exit
                cur = cur + (f(xpos) + f(xneg)) * w
                j = j + 2
            end do
            cur = h * cur + 0.5d0 * prv
            eps = abs(cur - prv)
            bnd = max(epsabs, epsrel * abs(cur))
            if (eps < bnd) exit
        end do
        res = cur
        abserr = eps
        return
    end subroutine tsquad

    double precision function uehlsph(x, rad)
        double precision, intent(in) :: x
        double precision, intent(in) :: rad

        double precision, parameter :: bound = 1.0d0
        double precision :: factor, term01, term02, term03

        factor = -0.5d0 / rad**3
        if (x == 0.0d0) then
            term01 = 0.2d0
            term02 = -rad * chip2(2.0d0 * rad)
            term03 = -0.5d0 * chip3(2.0d0 * rad)
            uehlsph = 2.0d0 * factor * (term01 + term02 + term03)
            return
        end if
        if (0.0d0 < x .and. x <= rad) then
            term01 = 0.4d0 * x
            term02 = phi(2.0d0 * (rad + x), rad)
            term03 = -phi(2.0d0 * (rad - x), rad)
            uehlsph = factor * (term01 + term02 + term03) / x
            return
        end if
        if (rad < x .and. x <= bound) then
            term01 = rad * chip3(2.0d0 * (x - rad))
            term02 = phi(2.0d0 * (x + rad), rad)
            term03 = -phi(2.0d0 * (x - rad), rad)
            uehlsph = factor * (term01 + term02 + term03) / x
            return
        end if
        if (bound < x .and. x <= 100.0d0) then
            uehlsph = usphasy(x, rad)
            return
        end if
        uehlsph = 0.0d0
        return
    end function uehlsph

    double precision function usphasy(x, rad)
        double precision, intent(in) :: x
        double precision, intent(in) :: rad

        double precision :: factor, twox, mom1, mom2, term01, term02, term03

        factor = -2.0d0 / 3.0d0
        mom1 = 0.6d0 * rad**2
        mom2 = 3.0d0 * rad**4 / 7.0d0
        twox = 2.0d0 * x
        term01 = chip1(twox)
        term02 = 2.0d0 * mom1 * chim1(twox) / 3.0d0
        term03 = 2.0d0 * mom2 * chim3(twox) / 15.0d0
        usphasy = factor * (term01 + term02 + term03) / x
        return
    end function usphasy

    double precision function phi(x, rad)
        double precision, intent(in) :: x, rad
        double precision :: a, b
        a = 0.5d0 * rad * chip3(x)
        b = 0.25d0 * chip4(x)
        phi = a + b
        return
    end function phi

    double precision function chip4(x)
        double precision, intent(in) :: x
        double precision, external :: dbesk0, dbesk1
        double precision, dimension(1) :: y
        double precision :: k0, k1, ki1, p1, p2, p3, t1, t2, t3
        integer :: nz, ierr
        external :: dbskin
        k0 = dbesk0(x)
        k1 = dbesk1(x)
        call dbskin(x, 1, 1, 1, y, nz, ierr)
        ki1 = y(1)
        p1 = x * (225.0d0 - x**2 * (48.0d0 + x**2))
        p2 = -x**2 * (351.0d0 - x**2 * (46.0d0 + x**2))
        p3 = 225.0d0 + x**2 * (405.0d0 - x**2 * (45.0d0 + x**2))
        t1 = p1 * k0
        t2 = p2 * k1
        t3 = p3 * ki1
        chip4 = (t1 + t2 + t3) / 1440.0d0
        return
    end function chip4

    double precision function chip3(x)
        double precision, intent(in) :: x
        double precision, external :: dbesk0, dbesk1
        double precision, dimension(1) :: y
        double precision :: k0, k1, ki1, p1, p2, p3, t1, t2, t3
        integer :: nz, ierr
        external :: dbskin
        k0 = dbesk0(x)
        k1 = dbesk1(x)
        call dbskin(x, 1, 1, 1, y, nz, ierr)
        ki1 = y(1)
        p1 = x**2 * (33.0d0 + x**2)
        p2 = x * (96.0d0 - x**2 * (31.0d0 + x**2))
        p3 = -x * (135.0d0 - x**2 * (30.0d0 + x**2))
        t1 = p1 * k0
        t2 = p2 * k1
        t3 = p3 * ki1
        chip3 = (t1 + t2 + t3) / 240.0d0
        return
    end function chip3

    double precision function chip2(x)
        double precision, intent(in) :: x
        double precision, external :: dbesk0, dbesk1
        double precision, dimension(1) :: y
        double precision :: k0, k1, ki1, p1, p2, p3, t1, t2, t3
        integer :: nz, ierr
        external :: dbskin
        k0 = dbesk0(x)
        k1 = dbesk1(x)
        call dbskin(x, 1, 1, 1, y, nz, ierr)
        ki1 = y(1)
        p1 = -x * (21.0d0 + x**2)
        p2 = x**2 * (19.0d0 + x**2)
        p3 = 27.0d0 - x**2 * (18.0d0 + x**2)
        t1 = p1 * k0
        t2 = p2 * k1
        t3 = p3 * ki1
        chip2 = (t1 + t2 + t3) / 48.0d0
        return
    end function chip2

    double precision function chip1(x)
        double precision, intent(in) :: x
        double precision, external :: dbesk0, dbesk1
        double precision, dimension(1) :: y
        double precision :: k0, k1, ki1, p1, p2, p3, t1, t2, t3
        integer :: nz, ierr
        external :: dbskin
        k0 = dbesk0(x)
        k1 = dbesk1(x)
        call dbskin(x, 1, 1, 1, y, nz, ierr)
        ki1 = y(1)
        p1 = 1.0d0 + x**2 / 12.0d0
        p2 = -x * (5.0d0 + 0.5d0 * x**2) / 6.0d0
        p3 = x * (3.0d0 + x**2 / 3.0d0) / 4.0d0
        t1 = p1 * k0
        t2 = p2 * k1
        t3 = p3 * ki1
        chip1 = t1 + t2 + t3
        return
    end function chip1

    double precision function chim1(x)
        double precision, intent(in) :: x
        double precision, external :: dbesk0, dbesk1
        double precision, dimension(1) :: y
        double precision :: k0, k1, ki1, p1, p2, p3, t1, t2, t3
        integer :: nz, ierr
        external :: dbskin
        k0 = dbesk0(x)
        k1 = dbesk1(x)
        call dbskin(x, 1, 1, 1, y, nz, ierr)
        ki1 = y(1)
        p1 = 0.5d0
        p2 = 1.0d0 / x - 0.5d0 * x
        p3 = 0.5d0 * x
        t1 = p1 * k0
        t2 = p2 * k1
        t3 = p3 * ki1
        chim1 = t1 + t2 + t3
        return
    end function chim1

    double precision function chim3(x)
        double precision, intent(in) :: x
        double precision, external :: dbesk0, dbesk1
        double precision :: k0, k1, p1, p2, t1, t2
        k0 = dbesk0(x)
        k1 = dbesk1(x)
        p1 = 3.0d0 / x**2
        p2 = 3.0d0 * (2.0d0 / x**2 + 0.5d0) / x
        t1 = p1 * k0
        t2 = p2 * k1
        chim3 = t1 + t2
        return
    end function chim3

end module bspqed_uehling
