module bspqed_delbrueck
    use iso_fortran_env, only: error_unit
    use bspqed_function, only: function_type
    use bspqed_quadrature, only: ngl20, nodes20, weights20
    implicit none
    private

    public :: delbrueck_potential

    double precision, parameter :: pi = 4.0d0 * atan(1.0d0)
    double precision, parameter :: pihalf = 0.5d0 * pi

    type, extends(function_type) :: delbrueck_potential
    contains
        procedure, pass(self) :: eval => delbrueck_potential_eval
    end type delbrueck_potential

contains

    double precision function delbrueck_potential_eval(x, self)
        double precision, intent(in) :: x
        class(delbrueck_potential), intent(in) :: self

        if (x < 0.0d0) then
            write(error_unit, *) "* ERROR in eval of delbrueck_potential"
            write(error_unit, *) "* x must be nonnegative."
            stop
        endif

        delbrueck_potential_eval = delbrckpot(x)
        return
    end function delbrueck_potential_eval

    double precision function delbrckpot(x)
        double precision, intent(in) :: x

        if (x == 0.0d0) then
            delbrckpot = 0.0d0
            return
        end if
        delbrckpot = 4.0d0 * ktilde(x) / x
        return
    end function delbrckpot

    double precision function ktilde(x)
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

        double precision, dimension(0:24) :: seq
        double precision :: xi, s
        integer :: k

        xi = abs(x)
        s = 0.0d0
        if (xi == 0.0d0) then
            s = 0.0d0
        else if (xi < 3.0d-01) then
            call sequence_small(xi, seq)
            do k = 24, 0, -1
                s = s + weights(k) * seq(k)
            end do
        else
            call sequence_large(xi, seq)
            do k = 24, 0, -1
                s = s + weights(k) * seq(k)
            end do
        end if
        ktilde = s
        return
    end function ktilde

    subroutine sequence_small(x, seqvals)
        double precision, intent(in) :: x
        double precision, dimension(0:), intent(out) :: seqvals

        double precision :: left, rght, seqn, term
        integer :: n, k

        n = size(seqvals) - 1

        ! compute seqvals(0)
        seqn = 0.0d0
        left = 0.0d0
        rght = pi * x
        do while (rght < pi)
            seqn = seqn + quadsinc_small(left, rght, x)
            left = rght
            rght = 2.0 * rght
        end do
        rght = pi
        seqn = seqn + quadsinc_small(left, rght, x)
        seqn = abs(seqn) / x
        seqvals(0) = seqn

        ! comupute seqvals(1)
        seqn = 0.0d0
        left = 0.0d0
        rght = 0.5d0 * pi * x
        do while (rght < 0.5d0 * pi)
            seqn = seqn + quadcos_small(left, rght, x)
            left = rght
            rght = 2.0 * rght
        end do
        rght = 0.5d0 * pi
        seqn = seqn + quadcos_small(left, rght, x)
        seqn = abs(seqn) / x
        term = seqn
        left = pi
        rght = pi * 2.0d0
        seqn = quadsinc_small(left, rght, x)
        seqn = abs(seqn) / x
        term = term + seqn
        seqvals(1) = term

        ! compute seqvals(k) for k = 2, 3, ..., n
        do k = 2, n
            left = 0.5d0 * pi * (2.0d0 * k - 3.0d0)
            rght = 0.5d0 * pi * (2.0d0 * k - 1.0d0)
            seqn = quadcos_small(left, rght, x)
            seqn = abs(seqn) / x
            term = seqn
            left = pi * k
            rght = pi * (k + 1.0d0)
            seqn = quadsinc_small(left, rght, x)
            seqn = abs(seqn) / x
            term = term + seqn
            seqvals(k) = term
        end do
        return
    end subroutine sequence_small

    double precision function quadcos_small(a, b, x)
        double precision, intent(in) :: a
        double precision, intent(in) :: b
        double precision, intent(in) :: x

        double precision :: s, q, f
        integer :: i

        s = 0.0d0
        do i = 1, ngl20
            q = 0.5d0 * (b - a) * (nodes20(i) + 1.0d0) + a
            f = delbrck(q / x) * cos(q)
            s = s + f * weights20(i)
        end do
        s = 0.5d0 * (b - a) * s
        quadcos_small = s
        return
    end function quadcos_small

    double precision function quadsinc_small(a, b, x)
        double precision, intent(in) :: a
        double precision, intent(in) :: b
        double precision, intent(in) :: x

        double precision :: s, q, f
        integer :: i

        s = 0.0d0
        do i = 1, ngl20
            q = 0.5d0 * (b - a) * (nodes20(i) + 1.0d0) + a
            f = delbrck(q / x) * sinc(q)
            s = s + f * weights20(i)
        end do
        s = 0.5d0 * (b - a) * s
        quadsinc_small = s
        return
    end function quadsinc_small

    subroutine sequence_large(x, seqvals)
        double precision, intent(in) :: x
        double precision, dimension(0:), intent(out) :: seqvals

        double precision :: left, rght, seqn, term
        integer :: n, k

        n = size(seqvals) - 1

        ! compute seqvals(0)
        left = 0.0d0
        rght = pi / x
        seqn = quadsinc_large(left, rght, x)
        seqn = abs(seqn)
        seqvals(0) = seqn

        ! comupute seqvals(1)
        left = 0.0d0
        rght = 0.5d0 * pi / x
        seqn = quadcos_large(left, rght, x)
        seqn = abs(seqn)
        term = seqn
        left = pi / x
        rght = pi / x * 2.0d0
        seqn = quadsinc_large(left, rght, x)
        seqn = abs(seqn)
        term = term + seqn
        seqvals(1) = term

        ! compute seqvals(k) for k = 2, 3, ..., n
        do k = 2, n
            left = 0.5d0 * pi / x * (2.0d0 * k - 3.0d0)
            rght = 0.5d0 * pi / x * (2.0d0 * k - 1.0d0)
            seqn = quadcos_large(left, rght, x)
            seqn = abs(seqn)
            term = seqn
            left = pi / x * k
            rght = pi / x * (k + 1.0d0)
            seqn = quadsinc_large(left, rght, x)
            seqn = abs(seqn)
            term = term + seqn
            seqvals(k) = term
        end do
        return
    end subroutine sequence_large

    double precision function quadcos_large(a, b, x)
        double precision, intent(in) :: a
        double precision, intent(in) :: b
        double precision, intent(in) :: x

        double precision :: s, q, f
        integer :: i

        s = 0.0d0
        do i = 1, ngl20
            q = 0.5d0 * (b - a) * (nodes20(i) + 1.0d0) + a
            f = delbrck(q) * cos(q * x)
            s = s + f * weights20(i)
        end do
        s = 0.5d0 * (b - a) * s
        quadcos_large = s
        return
    end function quadcos_large

    double precision function quadsinc_large(a, b, x)
        double precision, intent(in) :: a
        double precision, intent(in) :: b
        double precision, intent(in) :: x

        double precision :: s, q, f
        integer :: i

        s = 0.0d0
        do i = 1, ngl20
            q = 0.5d0 * (b - a) * (nodes20(i) + 1.0d0) + a
            f = delbrck(q) * sinc(q * x)
            s = s + f * weights20(i)
        end do
        s = 0.5d0 * (b - a) * s
        quadsinc_large = s
        return
    end function quadsinc_large

    double precision function sinc(x)
        double precision, intent(in) :: x

        double precision, dimension(11), parameter :: coefs = [&
            &  1.0000000000000000000000d+00, &
            & -1.6666666666666666666667d-01, &
            &  8.3333333333333333333333d-03, &
            & -1.9841269841269841269841d-04, &
            &  2.7557319223985890652557d-06, &
            & -2.5052108385441718775052d-08, &
            &  1.6059043836821614599392d-10, &
            & -7.6471637318198164759011d-13, &
            &  2.8114572543455207631989d-15, &
            & -8.2206352466243297169560d-18, &
            &  1.9572941063391261230848d-20  &
        ]
        double precision :: xi, arg

        xi = abs(x)
        if (xi <= 1.0d-04) then
            arg = xi**2
            sinc = polyeval(coefs, arg)
            return
        else
            sinc = sin(xi) / xi
            return
        end if
    end function sinc

    double precision function delbrck(x)
        double precision, intent(in) :: x
        if (x < 4.5d-01) then
            delbrck = delbrck0(x)
            return
        else
            delbrck = delbrck1(x)
            return
        end if
    end function delbrck

    double precision function delbrck0(x)
        double precision, intent(in) :: x
        double precision, dimension(14), parameter :: coefs = [&
            &  1.0000000000000000000d+00,  2.2857142857142857143d-01, &
            & -8.1428571428571428571d-02, -1.6326530612244897959d-01, &
            &  3.2025555758017492711d-02,  4.6530612244897959184d-02, &
            & -9.6983064058956916100d-03, -1.1774891774891774892d-02, &
            &  2.6266050732825413223d-03,  2.8814722692273712682d-03, &
            & -6.7767739998722517181d-04, -6.9726192175171767009d-04, &
            &  1.7075416899550840656d-04,  1.6813252107369754429d-04  &
        &]
        delbrck0 = polyeval(coefs, x) * 7.0d0 / 1152.0d0
        return
    end function delbrck0

    double precision function delbrck1(x)
        double precision, intent(in) :: x
        double precision, external :: dspenc
        double precision :: s, a1, a2, a3, a4, a5, a6, a7, a8
        s = sqrt(x**2 + 4.0d0)
        a1 = (3.0d0 / (4.0d0 * x**3) + 4.0d0 / x**5) * asinh(0.5d0 * x)**2
        a2 = 0.25d0 * (x**2 + 8.0d0) / (x**4 * s)
        a2 = a2 + 0.125d0 * (x**2 + 8.0d0) / x**3
        a2 = a2 - 0.5d0 * log(x) / x**3
        a2 = a2 * asinh(0.5d0 * x)
        a3 = -0.5d0 * (3.0d0 * x**2 + 16.0d0) * dspenc(-0.5d0 * (s + x)) / x**5
        a4 = 0.5d0 * dspenc(0.5d0 * (s - x)) / x**3
        a5 = 0.5d0 * (x**2 + 8.0d0) * dspenc(-0.25d0 * (s + x)**2) / x**5
        a6 = (x**4 + 32.0d0 * x**2 + 128.0d0) / (16.0d0 * x**4 * s)
        a7 = 8.0d0 * pi**2 * (x**2 + 2.0d0)
        a7 = a7 + 3.0d0 * x * (x**3 + 4.0d0 * x**2 + 16.0d0 * x + 64.0d0)
        a7 = a7 / (-48.0d0 * x**5)
        a8 = 0.125d0 * (x**2 + 8.0d0) * s * log((s + 2.0d0) / s**2) / x**4
        delbrck1 = a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8
        return
    end function delbrck1

    double precision function polyeval(coefs, x)
        double precision, dimension(:), intent(in) :: coefs
        double precision, intent(in) :: x

        double precision :: y
        integer :: n, k

        n = size(coefs)
        y = coefs(n)
        do k = n-1, 1, -1
            y = y * x + coefs(k)
        end do
        polyeval = y
        return
    end function polyeval

end module bspqed_delbrueck
