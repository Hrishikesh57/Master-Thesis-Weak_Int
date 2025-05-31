module bspse_onepotential
    use iso_fortran_env, only: error_unit
    use bspqed_function, only: function_type
    use bspqed_momentum, only: momentumgf
    implicit none
    private

    public :: onepot

    double precision, parameter :: pihalf = 2.0d0 * atan(1.0d0)
    double precision, parameter :: sqrtpihalf = sqrt(pihalf)

contains

    subroutine onepot(knots, order, kapa, epsa, d1, d2, densq, res, abserr, epsabs, epsrel)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, intent(in) :: kapa
        double precision, intent(in) :: epsa
        double precision, dimension(0:, :), intent(in) :: d1
        double precision, dimension(0:, :), intent(in) :: d2
        class(function_type), intent(in) :: densq

        double precision, intent(out) :: res
        double precision, intent(out) :: abserr

        double precision, optional, intent(in) :: epsabs
        double precision, optional, intent(in) :: epsrel

        integer, parameter :: limit = 100
        integer, parameter :: lenw = 4 * limit

        double precision, dimension(lenw) :: work
        integer, dimension(limit) :: iwork

        double precision :: epsab, epsrl
        integer :: neval, ier, last

        epsab = 1.0d-10
        if (present(epsabs)) epsab = epsabs
        epsrl = 5.0d-07
        if (present(epsrel)) epsrl = epsrel
        call dqagi(f, 0.0d0, 1, epsab, epsrl, res, abserr, neval, ier, limit, lenw, last, iwork, work)
        res = res
        abserr = abserr
        return

    contains

        double precision function f(x)
            double precision, intent(in) :: x
            f = integral2(knots, order, kapa, epsa, d1, d2, densq, x, epsab/100.0d0, epsrl/100.0d0) * sqrtpihalf
            return
        end function f
    end subroutine onepot

    double precision function integral2(knots, order, kapa, epsa, d1, d2, densq, x, epsabs, epsrel)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, intent(in) :: kapa
        double precision, intent(in) :: epsa
        double precision, dimension(0:, :), intent(in) :: d1
        double precision, dimension(0:, :), intent(in) :: d2
        class(function_type), intent(in) :: densq
        double precision, intent(in) :: x
        double precision, intent(in) :: epsabs
        double precision, intent(in) :: epsrel

        integer, parameter :: limit = 100
        integer, parameter :: lenw = 4 * limit

        double precision, dimension(lenw) :: work
        integer, dimension(limit) :: iwork

        double precision :: res, abserr
        integer :: neval, ier, last

        call dqags(f, 0.0d0, x, epsabs, epsrel, res, abserr, neval, ier, limit, lenw, last, iwork, work)
        integral2 = res
        return

    contains

        double precision function f(y)
            double precision, intent(in) :: y
            f = integral1(knots, order, kapa, epsa, d1, d2, densq, x, y, epsabs/100.0d0, epsrel/100.0d0)
            return
        end function f
    end function integral2

    double precision function integral1(knots, order, kapa, epsa, d1, d2, densq, x, y, epsabs, epsrel)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, intent(in) :: kapa
        double precision, intent(in) :: epsa
        double precision, dimension(0:, :), intent(in) :: d1
        double precision, dimension(0:, :), intent(in) :: d2
        class(function_type), intent(in) :: densq
        double precision, intent(in) :: x
        double precision, intent(in) :: y
        double precision, intent(in) :: epsabs
        double precision, intent(in) :: epsrel

        integer, parameter :: limit = 100
        integer, parameter :: lenw = 4 * limit

        double precision, dimension(lenw) :: work
        integer, dimension(limit) :: iwork

        double precision :: p1, p2, g1, f1, g2, f2, a, b, res, abserr
        integer :: l1, l2, neval, ier, last

        integral1 = 0.0d0
        if (x == 0.0d0) return
        p1 = 0.5d0 * (x + y)
        p2 = 0.5d0 * (x - y)
        call momentumgf(knots, order, kapa, d1, d2, p1, g1, f1)
        call momentumgf(knots, order, kapa, d1, d2, p2, g2, f2)
        l1 = (abs(2 * kapa + 1) - 1) / 2
        l2 = (abs(2 * kapa - 1) - 1) / 2
        if (y == 0.0d0) then
            a = -log(x)
            call dqagi(f, a, 1, epsabs, epsrel, res, abserr, neval, ier, limit, lenw, last, iwork, work)
        else
            a = -log(x)
            b = -log(y)
            call dqags(f, a, b, epsabs, epsrel, res, abserr, neval, ier, limit, lenw, last, iwork, work)
        end if
        integral1 = res
        return

    contains

        double precision function f(v)
            double precision, intent(in) :: v
            double precision :: q
            q = exp(-v)
            f = formfunc(l1, l2, epsa, q, p1, g1, f1, p2, g2, f2, densq)
            return
        end function f
    end function integral1

    double precision function formfunc(l1, l2, epsa, q, p1, g1, f1, p2, g2, f2, densq)
        integer, intent(in) :: l1
        integer, intent(in) :: l2
        double precision, intent(in) :: epsa
        double precision, intent(in) :: q
        double precision, intent(in) :: p1
        double precision, intent(in) :: g1
        double precision, intent(in) :: f1
        double precision, intent(in) :: p2
        double precision, intent(in) :: g2
        double precision, intent(in) :: f2
        class(function_type), intent(in) :: densq

        double precision, dimension(0:max(l1, l2)) :: legp
        double precision :: xi, rhoq, ff1, ff2

        formfunc = 0.0d0
        if (p1 * p2 == 0.0d0) return
        rhoq = densq%eval(q)
        call formfactors(epsa, q**2, p1, g1, f1, p2, g2, f2, ff1, ff2)
        xi = 0.5d0 * (p1**2 + p2**2 - q**2) / (p1 * p2)
        call legendre(xi, legp)
        formfunc = rhoq * (ff1 * legp(l1) + ff2 * legp(l2))
        return
    end function formfunc

    subroutine legendre(x, p)
        double precision, intent(in) :: x

        double precision, dimension(0:), intent(out) :: p

        double precision :: pnm2, pnm1, pn
        integer :: nmax, n

        nmax = size(p) - 1
        pnm2 = 0.0d0
        pnm1 = 1.0d0
        p(0) = 1.0d0
        do n = 1, nmax
            pn = ((2 * n - 1) * x * pnm1 - (n - 1) * pnm2) / n
            pnm2 = pnm1
            pnm1 = pn
            p(n) = pn
        end do
        return
    end subroutine legendre

    subroutine formfactors(epsa, q2, p1, g1, f1, p2, g2, f2, ff1, ff2)
        double precision, intent(in) :: epsa
        double precision, intent(in) :: q2
        double precision, intent(in) :: p1
        double precision, intent(in) :: g1
        double precision, intent(in) :: f1
        double precision, intent(in) :: p2
        double precision, intent(in) :: g2
        double precision, intent(in) :: f2

        double precision, intent(out) :: ff1
        double precision, intent(out) :: ff2

        double precision :: a, b1, b2, c1, c2, d, h1, h2, b12, c12, h12, u1, u2

        call abcdh(epsa**2, p1**2, p2**2, q2, a, b1, b2, c1, c2, d, h1, h2)
        b12 = b1 + b2
        c12 = c1 + c2
        h12 = h1 + h2
        u1 = epsa * g1 + p1 * f1
        u2 = epsa * g2 + p2 * f2
        ff1 = epsa * (b12 * u1 * g2 + c12 * g1 * u2) + d * u1 * u2 + (a + epsa * h12) * g1 * g2
        u1 = epsa * f1 + p1 * g1
        u2 = epsa * f2 + p2 * g2
        ff2 = epsa * (b12 * u1 * f2 + c12 * f1 * u2) + d * u1 * u2 + (a - epsa * h12) * f1 * f2
        return
    end subroutine formfactors

    subroutine abcdh(epsa2, p12, p22, q2, a, b1, b2, c1, c2, d, h1, h2)
        double precision, intent(in) :: epsa2
        double precision, intent(in) :: p12
        double precision, intent(in) :: p22
        double precision, intent(in) :: q2

        double precision, intent(out) :: a
        double precision, intent(out) :: b1
        double precision, intent(out) :: b2
        double precision, intent(out) :: c1
        double precision, intent(out) :: c2
        double precision, intent(out) :: d
        double precision, intent(out) :: h1
        double precision, intent(out) :: h2

        double precision :: c0, c11, c12, c21, c22, c23, c24

        call coefs(epsa2, p12, p22, q2, c0, c11, c12, c21, c22, c23, c24)
        a = c24 - 2.0d0 + c11 + c12 - 2.0d0 * c0 + epsa2 * (4.0d0 * c0 + 5.0d0 * (c11 + c12)) &
        & - p12 * (2.0d0 * (c0 + c12) + 3.0d0 * c11) - p22 * (2.0d0 * (c0 + c11) + 3.0d0 * c12) &
        & + 2.0d0 * q2 * (c0 + c11 + c12)
        b1 = -4.0d0 * (c11 + c21)
        b2 = -4.0d0 * (c0 + c11 + c12 + c23)
        c1 = b2
        c2 = -4.0d0 * (c12 + c22)
        d = 2.0d0 * (c0 + c11 + c12)
        h1 = 4.0d0 * (c0 + 2.0d0 * c11)
        h2 = 4.0d0 * (c0 + 2.0d0 * c12)
        return
    end subroutine abcdh

    subroutine coefs(epsa2, p12, p22, q2, c0, c11, c12, c21, c22, c23, c24)
        double precision, intent(in) :: epsa2
        double precision, intent(in) :: p12
        double precision, intent(in) :: p22
        double precision, intent(in) :: q2

        double precision, intent(out) :: c0
        double precision, intent(out) :: c11
        double precision, intent(out) :: c12
        double precision, intent(out) :: c21
        double precision, intent(out) :: c22
        double precision, intent(out) :: c23
        double precision, intent(out) :: c24

        integer, parameter :: limit = 100
        integer, parameter :: lenw = 4 * limit

        double precision, dimension(lenw) :: work
        integer, dimension(limit) :: iwork

        double precision, parameter :: a = 0.0d0
        double precision, parameter :: b = 1.0d0
        double precision, parameter :: epsabs = 0.0d0
        double precision, parameter :: epsrel = 5.0d-05

        double precision :: abserr
        integer :: neval, ier, last

        call dqags(f0,  a, b, epsabs, epsrel, c0,  abserr, neval, ier, limit, lenw, last, iwork, work)
        call dqags(f11, a, b, epsabs, epsrel, c11, abserr, neval, ier, limit, lenw, last, iwork, work)
        call dqags(f12, a, b, epsabs, epsrel, c12, abserr, neval, ier, limit, lenw, last, iwork, work)
        call dqags(f21, a, b, epsabs, epsrel, c21, abserr, neval, ier, limit, lenw, last, iwork, work)
        call dqags(f22, a, b, epsabs, epsrel, c22, abserr, neval, ier, limit, lenw, last, iwork, work)
        call dqags(f23, a, b, epsabs, epsrel, c23, abserr, neval, ier, limit, lenw, last, iwork, work)
        c24 = coef24(q2)
        return

    contains

        double precision function f0(y)
            double precision, intent(in) :: y
            double precision :: a1, a2, t
            a1 = funca1(epsa2, p12, p22, q2, y)
            a2 = funca2(epsa2, p12, p22, y)
            t = a1 / a2
            f0 = func0(t) / a2
            return
        end function f0

        double precision function f11(y)
            double precision, intent(in) :: y
            double precision :: a1, a2, t
            a1 = funca1(epsa2, p12, p22, q2, y)
            a2 = funca2(epsa2, p12, p22, y)
            t = a1 / a2
            f11 = y * func1(t) / a2
            return
        end function f11

        double precision function f12(y)
            double precision, intent(in) :: y
            double precision :: a1, a2, t
            a1 = funca1(epsa2, p12, p22, q2, y)
            a2 = funca2(epsa2, p12, p22, y)
            t = a1 / a2
            f12 = (1.0d0 - y) * func1(t) / a2
            return
        end function f12

        double precision function f21(y)
            double precision, intent(in) :: y
            double precision :: a1, a2, t
            a1 = funca1(epsa2, p12, p22, q2, y)
            a2 = funca2(epsa2, p12, p22, y)
            t = a1 / a2
            f21 = y**2 * func2(t) / a2
            return
        end function f21

        double precision function f22(y)
            double precision, intent(in) :: y
            double precision :: a1, a2, t
            a1 = funca1(epsa2, p12, p22, q2, y)
            a2 = funca2(epsa2, p12, p22, y)
            t = a1 / a2
            f22 = (1.0d0 - y)**2 * func2(t) / a2
            return
        end function f22

        double precision function f23(y)
            double precision, intent(in) :: y
            double precision :: a1, a2, t
            a1 = funca1(epsa2, p12, p22, q2, y)
            a2 = funca2(epsa2, p12, p22, y)
            t = a1 / a2
            f23 = y * (1.0d0 - y) * func2(t) / a2
            return
        end function f23
    end subroutine coefs

    double precision function coef24(q2)
        double precision, intent(in) :: q2
        double precision :: x2, x
        x2 = q2 / (4.0d0 + q2)
        if (x2 < 0.0025d0) then
            coef24 = -2.0d0 * x2 * (1.0d0 / 3.0d0 + x2 * (1.0d0 / 5.0d0 + x2 * (1.0d0 / 7.0d0 + x2 * (1.0d0 /&
            & 9.0d0 + x2 / 11.0d0))))
            return
        end if
        x = sqrt(x2)
        coef24 = 2.0d0 * (1.0d0 - atanh(x) / x)
        return
    end function coef24

    double precision function funca0(epsa2, p12, p22, y)
        double precision, intent(in) :: epsa2
        double precision, intent(in) :: p12
        double precision, intent(in) :: p22
        double precision, intent(in) :: y
        funca0 = epsa2 - y * p12 - (1.0d0 - y) * p22
        return
    end function funca0

    double precision function funca1(epsa2, p12, p22, q2, y)
        double precision, intent(in) :: epsa2
        double precision, intent(in) :: p12
        double precision, intent(in) :: p22
        double precision, intent(in) :: q2
        double precision, intent(in) :: y
        double precision :: a0
        a0 = funca0(epsa2, p12, p22, y)
        funca1 = y * (1.0d0 - y) * q2 + a0
        return
    end function funca1

    double precision function funca2(epsa2, p12, p22, y)
        double precision, intent(in) :: epsa2
        double precision, intent(in) :: p12
        double precision, intent(in) :: p22
        double precision, intent(in) :: y
        double precision :: a0
        a0 = funca0(epsa2, p12, p22, y)
        funca2 = 1.0d0 - a0
        return
    end function funca2

    double precision function func0(t)
        double precision, intent(in) :: t
        double precision :: y, x
        if (abs(t) < 1.0d-02) then
            func0 = expandf(0, 5, t)
            return
        end if
        y = 1.0d0 / t
        x = 1.0d0 + 1.0d0 / y
        func0 = -log(x) * y
        return
    end function func0

    double precision function func1(t)
        double precision, intent(in) :: t
        double precision :: y, x
        if (abs(t) < 1.0d-02) then
            func1 = expandf(1, 5, t)
            return
        end if
        y = 1.0d0 / t
        x = 1.0d0 + 1.0d0 / y
        func1 = (1.0d0 - y * log(x)) * y
        return
    end function func1

    double precision function func2(t)
        double precision, intent(in) :: t
        double precision :: y, x
        if (abs(t) < 1.0d-02) then
            func2 = expandf(2, 5, t)
            return
        end if
        y = 1.0d0 / t
        x = 1.0d0 + 1.0d0 / y
        func2 = (-0.5d0 + y - y**2 * log(x)) * y
        return
    end function func2

    double precision function expandf(m, n, t)
        integer, intent(in) :: m
        integer, intent(in) :: n
        double precision, intent(in) :: t

        double precision :: x, y
        integer :: k

        x = -t
        y = 0.0d0
        do k = n, 1, -1
            y = y * x + 1.0d0 / (k + m)
        end do
        expandf = (-1)**(m+1) * y
        return
    end function expandf

end module bspse_onepotential
