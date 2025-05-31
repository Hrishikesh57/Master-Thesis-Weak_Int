module bspwf_dirac
    use iso_fortran_env, only: error_unit
    use bspwf_constants, only: ialpha
    use bspwf_function, only: function_type
    implicit none
    private

    public :: dirac
    public :: radialg
    public :: radialf
    public :: radialgf

    double precision, parameter :: speed_of_light_in_au = ialpha

contains

    subroutine dirac(knots, order, kappa, potential, nodes, weights, spectrum, coefmat, hfmat, speed_of_light)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, intent(in) :: kappa
        class(function_type), intent(in) :: potential
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision, dimension(:), intent(out) :: spectrum
        double precision, dimension(:, :), intent(out) :: coefmat

        double precision, dimension(:, :), intent(in), optional :: hfmat
        double precision, intent(in), optional :: speed_of_light

        double precision, dimension(2*(size(knots)-order-2), 2*(size(knots)&
        &-order-2)) :: a, b
        double precision, dimension(order, 0:2) :: dbsx
        double precision, dimension((order+2)*(order+1)/2) :: work1
        double precision, dimension(6*(size(knots)-order-2)) :: work2
        double precision :: c, left, right, halfd, x, w, v, b0k, b1k, b2k, b0j, b1j, b2j, fnc
        integer :: nknots, nintervals, nbsplines, nspectrum, ngl, s, ileft, imin, imax, lwork2, info, j, k, l, m

        external :: dbspvd, dsygv

        c = speed_of_light_in_au
        if (present(speed_of_light)) c = speed_of_light

        nknots = size(knots)
        nintervals = nknots - 2 * (order - 1) - 1
        nbsplines = nknots - order
        nspectrum = 2 * (nbsplines - 2)
        ngl = size(nodes)
        s = nbsplines - 2

        a = 0.0d0
        b = 0.0d0
        do l = 1, nintervals
            ileft = order - 1 + l
            left = knots(ileft)
            right = knots(ileft+1)
            halfd = 0.5d0 * (right - left)
            imin = ileft - order + 1
            if (imin == 1) imin = 2
            imax = ileft
            if (imax == nbsplines) imax = nbsplines - 1
            do m = 1, ngl
                x = halfd * (nodes(m) + 1.0d0) + left
                w = halfd * weights(m)
                v = potential%eval(x)
                call dbspvd(knots, order, 3, x, ileft, order, dbsx, work1)
                do k = imin, imax
                    b0k = dbsx(k+order-ileft, 0)
                    b1k = dbsx(k+order-ileft, 1)
                    b2k = dbsx(k+order-ileft, 2)
                    do j = imin, imax
                        b0j = dbsx(j+order-ileft, 0)
                        b1j = dbsx(j+order-ileft, 1)
                        b2j = dbsx(j+order-ileft, 2)

                        fnc = (v * (1.0d0 + (0.5d0 * kappa / (c * x))**2) &
                        & + c**2 + 0.75d0 * kappa * (1 + kappa) / x**2) * b0j * b0k &
                        & + 0.25d0 * (c**2 + v) * b1j * b1k / c**2 &
                        & + 0.25d0 * v * kappa / x * (b1j * b0k + b0j * b1k) / c**2 &
                        & - 0.25d0 * (b2j * b0k + b0j * b2k)
                        a(j-1, k-1) = a(j-1, k-1) + fnc * w

                        fnc = (1.0d0 + 0.25d0 * kappa * (1 + kappa) / (c * x)**2) * b0j * b0k &
                        & + 0.25d0 * b1j * b1k / c**2
                        b(j-1, k-1) = b(j-1, k-1) + fnc * w

                        fnc = 0.5d0 * v * (b1j * b0k + b0j * b1k) / c &
                        & - 0.125d0 * (2.0d0 * (kappa / x)**3 * b0j * b0k &
                        & - kappa * (1 - kappa) / x**2 * b1j * b0k &
                        & - kappa * (1 + kappa) / x**2 * b0j * b1k &
                        & - kappa / x * (b2j * b0k + b0j * b2k) &
                        & + b2j * b1k - b1j * b2k) / c
                        a(j-1, s+k-1) = a(j-1, s+k-1) + fnc * w

                        fnc = 0.5d0 * v * (b1k * b0j + b0k * b1j) / c &
                        & - 0.125d0 * (2.0d0 * (kappa / x)**3 * b0k * b0j &
                        & - kappa * (1 - kappa) / x**2 * b1k * b0j &
                        & - kappa * (1 + kappa) / x**2 * b0k * b1j &
                        & - kappa / x * (b2k * b0j + b0k * b2j) &
                        & + b2k * b1j - b1k * b2j) / c
                        a(s+j-1, k-1) = a(s+j-1, k-1) + fnc * w

                        fnc = (v * (1.0d0 + (0.5d0 * kappa / (c * x))**2) &
                        & - c**2 + 0.75d0 * kappa * (1 - kappa) / x**2) * b0j * b0k &
                        & - 0.25d0 * (c**2 - v) * b1j * b1k / c**2 &
                        & - 0.25d0 * v * kappa / x * (b1j * b0k + b0j * b1k) / c**2 &
                        & + 0.25d0 * (b2j * b0k + b0j * b2k)
                        a(s+j-1, s+k-1) = a(s+j-1, s+k-1) + fnc * w

                        fnc = (1.0d0 - 0.25d0 * kappa * (1 - kappa) / (c * x)**2) * b0j * b0k &
                        & + 0.25d0 * b1j * b1k / c**2
                        b(s+j-1, s+k-1) = b(s+j-1, s+k-1) + fnc * w
                    end do
                end do
            end do
        end do
        a = 0.5d0 * (a + transpose(a))
        b = 0.5d0 * (b + transpose(b))

        if (present(hfmat)) a = a + hfmat

        spectrum = 0.0d0
        lwork2 = 3 * nspectrum
        info = 0
        call dsygv(1, 'v', 'u', nspectrum, a, nspectrum, b, nspectrum, spectrum, work2, lwork2, info)
        if (info /= 0) then
            write (error_unit, *) "* ERROR in 'dirac'"
            write (error_unit, *) "* An attempt to solve the generalized"
            write (error_unit, *) "* eigenvalue equation failed."
            stop
        end if

        !+ Restore coefficient of the b-spline basis functions which have been
        !+ thrown out.  That is, we leave the corresponding entries in the
        !+ matrix 'coefmat' to be zero and copy the other entries from the
        !+ matrix 'amat'.  The rows which need to be zero in 'coefmat' are: 1,
        !+ nbsplines, nbsplines+1, and 2*nbsplines.
        coefmat = 0.0d0
        coefmat(2:(nbsplines-1), :) = a(1:(nbsplines-2), :)
        coefmat((nbsplines+2):(2*nbsplines-1), :) = a((nbsplines-1):, :)
        return
    end subroutine dirac

    double precision function radialg(knots, order, kappa, coefs, x, inbv, speed_of_light)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, intent(in) :: kappa
        double precision, dimension(:), intent(in) :: coefs
        double precision, intent(in) :: x

        integer, intent(in out) :: inbv

        double precision, intent(in), optional :: speed_of_light

        double precision :: c, k1, k2, dk1, dk2
        integer :: nbsplines

        c = speed_of_light_in_au
        if (present(speed_of_light)) c = speed_of_light
        nbsplines = size(knots) - order
        call ksubr(knots, order, coefs(1:nbsplines), x, k1, dk1, inbv)
        call ksubr(knots, order, coefs((nbsplines+1):), x, k2, dk2, inbv)
        radialg = k1 + 0.5d0 * (dk2 - kappa / x * k2) / c
        return
    end function radialg

    double precision function radialf(knots, order, kappa, coefs, x, inbv, speed_of_light)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, intent(in) :: kappa
        double precision, dimension(:), intent(in) :: coefs
        double precision, intent(in) :: x

        integer, intent(in out) :: inbv

        double precision, intent(in), optional :: speed_of_light

        double precision :: c, k1, k2, dk1, dk2
        integer :: nbsplines

        c = speed_of_light_in_au
        if (present(speed_of_light)) c = speed_of_light
        nbsplines = size(knots) - order
        call ksubr(knots, order, coefs(1:nbsplines), x, k1, dk1, inbv)
        call ksubr(knots, order, coefs((nbsplines+1):), x, k2, dk2, inbv)
        radialf = k2 + 0.5d0 * (dk1 + kappa / x * k1) / c
        return
    end function radialf

    subroutine radialgf(knots, order, kappa, coefs, x, gx, fx, inbv, speed_of_light)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, intent(in) :: kappa
        double precision, dimension(:), intent(in) :: coefs
        double precision, intent(in) :: x

        double precision, intent(out) :: gx
        double precision, intent(out) :: fx

        integer, intent(in out) :: inbv

        double precision, intent(in), optional :: speed_of_light

        double precision :: c, k1, k2, dk1, dk2
        integer :: nbsplines

        c = speed_of_light_in_au
        if (present(speed_of_light)) c = speed_of_light
        nbsplines = size(knots) - order
        call ksubr(knots, order, coefs(1:nbsplines), x, k1, dk1, inbv)
        call ksubr(knots, order, coefs((nbsplines+1):), x, k2, dk2, inbv)
        gx = k1 + 0.5d0 * (dk2 - kappa / x * k2) / c
        fx = k2 + 0.5d0 * (dk1 + kappa / x * k1) / c
        return
    end subroutine radialgf

    subroutine ksubr(knots, order, coefs, x, k, dk, inbv)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        double precision, dimension(:), intent(in) :: coefs
        double precision, intent(in) :: x

        double precision, intent(out) :: k
        double precision, intent(out) :: dk

        integer, intent(in out) :: inbv

        double precision, dimension(3*order) :: work
        double precision, external :: dbvalu
        integer :: ncoefs

        ncoefs = size(coefs)
        k = dbvalu(knots, coefs, ncoefs, order, 0, x, inbv, work)
        dk = dbvalu(knots, coefs, ncoefs, order, 1, x, inbv, work)
        return
    end subroutine ksubr

end module bspwf_dirac
