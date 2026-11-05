module bspqed_dkb
    use iso_fortran_env, only: error_unit
    use bspqed_function, only: function_type
    implicit none
    private

    public :: solve
    public :: radialg
    public :: radialf
    public :: radialgf
    public :: ufunc
    public :: vfunc
    public :: wfunc
    public :: pfunc
    public :: qfunc

contains

    subroutine solve(knots, order, kappa, potential, nodes, weights, spectrum, coefmat)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, intent(in) :: kappa
        class(function_type), intent(in) :: potential
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision, dimension(:), intent(out) :: spectrum
        double precision, dimension(:, :), intent(out) :: coefmat

        double precision, dimension(2*(size(knots)-order-2), 2*(size(knots)-order-2)) :: a, b
        double precision, dimension(order, 0:2) :: dbsx
        double precision, dimension((order+2)*(order+1)/2) :: work1
        double precision, dimension(6*(size(knots)-order-2)) :: work2
        double precision :: left, rght, hlfd, x, w, v, b0k, b1k, b2k, b0j, b1j, b2j, fnc
        integer :: nknots, nintervals, nbsplines, nspectrum, ngl, s, ileft, imin, imax, lwork2, info, j, k, l, m
        external :: dbspvd, dsygv

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
            rght = knots(ileft+1)
            hlfd = 0.5d0 * (rght - left)
            imin = ileft - order + 1
            if (imin == 1) imin = 2
            imax = ileft
            if (imax == nbsplines) imax = nbsplines - 1
            do m = 1, ngl
                x = hlfd * (nodes(m) + 1.0d0) + left
                w = hlfd * weights(m)
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
                        fnc = (v * (1.0d0 + (0.5d0 * kappa / x)**2) + 1.0d0 &
                        & + 0.75d0 * kappa * (1 + kappa) / x**2) * b0j * b0k &
                        & + 0.25d0 * (1.0d0 + v) * b1j * b1k &
                        & + 0.25d0 * v * kappa / x * (b1j * b0k + b0j * b1k) &
                        & - 0.25d0 * (b2j * b0k + b0j * b2k)
                        a(j-1, k-1) = a(j-1, k-1) + fnc * w
                        fnc = (1.0d0 + 0.25d0 * kappa * (1 + kappa) / x**2) &
                        & * b0j * b0k + 0.25d0 * b1j * b1k
                        b(j-1, k-1) = b(j-1, k-1) + fnc * w
                        fnc = 0.5d0 * v * (b1j * b0k + b0j * b1k) &
                        & - 0.125d0 * (2.0d0 * (kappa / x)**3 * b0j * b0k &
                        & - kappa * (1 - kappa) / x**2 * b1j * b0k &
                        & - kappa * (1 + kappa) / x**2 * b0j * b1k &
                        & - kappa / x * (b2j * b0k + b0j * b2k) &
                        & + b2j * b1k - b1j * b2k)
                        a(j-1, s+k-1) = a(j-1, s+k-1) + fnc * w
                        fnc = 0.5d0 * v * (b1k * b0j + b0k * b1j) &
                        & - 0.125d0 * (2.0d0 * (kappa / x)**3 * b0k * b0j &
                        & - kappa * (1 - kappa) / x**2 * b1k * b0j &
                        & - kappa * (1 + kappa) / x**2 * b0k * b1j &
                        & - kappa / x * (b2k * b0j + b0k * b2j) &
                        & + b2k * b1j - b1k * b2j)
                        a(s+j-1, k-1) = a(s+j-1, k-1) + fnc * w
                        fnc = (v * (1.0d0 + (0.5d0 * kappa / x)**2) - 1.0d0 &
                        & + 0.75d0 * kappa * (1 - kappa) / x**2) * b0j * b0k &
                        & - 0.25d0 * (1.0d0 - v) * b1j * b1k &
                        & - 0.25d0 * v * kappa / x * (b1j * b0k + b0j * b1k) &
                        & + 0.25d0 * (b2j * b0k + b0j * b2k)
                        a(s+j-1, s+k-1) = a(s+j-1, s+k-1) + fnc * w
                        fnc = (1.0d0 - 0.25d0 * kappa * (1 - kappa) / x**2) &
                        & * b0j * b0k + 0.25d0 * b1j * b1k
                        b(s+j-1, s+k-1) = b(s+j-1, s+k-1) + fnc * w
                    end do
                end do
            end do
        end do
        a = 0.5d0 * (a + transpose(a))
        b = 0.5d0 * (b + transpose(b))

        spectrum = 0.0d0
        lwork2 = 3 * nspectrum
        info = 0
        call dsygv(1, 'v', 'u', nspectrum, a, nspectrum, b, nspectrum, spectrum, work2, lwork2, info)
        if (info /= 0) then
            write (error_unit, *) "* ERROR in solve from bspqed_dkb:"
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
    end subroutine solve

    double precision function radialg(knots, order, kappa, coefs, x, inbv)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, intent(in) :: kappa
        double precision, dimension(:), intent(in) :: coefs
        double precision, intent(in) :: x

        integer, intent(in out) :: inbv

        double precision :: k1, k2, dk1, dk2
        integer :: nbsplines

        nbsplines = size(knots) - order
        call ksubr(knots, order, coefs(1:nbsplines), x, k1, dk1, inbv)
        call ksubr(knots, order, coefs((nbsplines+1):), x, k2, dk2, inbv)
        radialg = k1 + 0.5d0 * (dk2 - kappa / x * k2)
        return
    end function radialg

    double precision function radialf(knots, order, kappa, coefs, x, inbv)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, intent(in) :: kappa
        double precision, dimension(:), intent(in) :: coefs
        double precision, intent(in) :: x

        integer, intent(in out) :: inbv

        double precision :: k1, k2, dk1, dk2
        integer :: nbsplines

        nbsplines = size(knots) - order
        call ksubr(knots, order, coefs(1:nbsplines), x, k1, dk1, inbv)
        call ksubr(knots, order, coefs((nbsplines+1):), x, k2, dk2, inbv)
        radialf = k2 + 0.5d0 * (dk1 + kappa / x * k1)
        return
    end function radialf

    subroutine radialgf(knots, order, kappa, coefs, x, gx, fx, inbv)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, intent(in) :: kappa
        double precision, dimension(:), intent(in) :: coefs
        double precision, intent(in) :: x

        double precision, intent(out) :: gx
        double precision, intent(out) :: fx

        integer, intent(in out) :: inbv

        double precision :: k1, k2, dk1, dk2
        integer :: nbsplines

        nbsplines = size(knots) - order
        call ksubr(knots, order, coefs(1:nbsplines), x, k1, dk1, inbv)
        call ksubr(knots, order, coefs((nbsplines+1):), x, k2, dk2, inbv)
        gx = k1 + 0.5d0 * (dk2 - kappa / x * k2)
        fx = k2 + 0.5d0 * (dk1 + kappa / x * k1)
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

    double precision function ufunc(knots, order, kapa, coefa, kapb, coefb, x, inbv)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: coefa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: coefb
        double precision, intent(in) :: x

        integer, intent(in out) :: inbv

        double precision :: ga, fa, gb, fb

        call radialgf(knots, order, kapa, coefa, x, ga, fa, inbv)
        call radialgf(knots, order, kapb, coefb, x, gb, fb, inbv)
        ufunc = fa * gb - ga * fb
        return
    end function ufunc

    double precision function vfunc(knots, order, kapa, coefa, kapb, coefb, x, inbv)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: coefa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: coefb
        double precision, intent(in) :: x

        integer, intent(in out) :: inbv

        double precision :: ga, fa, gb, fb

        call radialgf(knots, order, kapa, coefa, x, ga, fa, inbv)
        call radialgf(knots, order, kapb, coefb, x, gb, fb, inbv)
        vfunc = -(ga * fb + fa * gb)
        return
    end function vfunc

    double precision function wfunc(knots, order, kapa, coefa, kapb, coefb, x, inbv)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: coefa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: coefb
        double precision, intent(in) :: x

        integer, intent(in out) :: inbv

        double precision :: ga, fa, gb, fb

        call radialgf(knots, order, kapa, coefa, x, ga, fa, inbv)
        call radialgf(knots, order, kapb, coefb, x, gb, fb, inbv)
        wfunc = ga * gb + fa * fb
        return
    end function wfunc

    double precision function pfunc(knots, order, kapa, coefa, kapb, coefb, k, x, inbv)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: coefa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: coefb
        integer, intent(in) :: k
        double precision, intent(in) :: x

        integer, intent(in out) :: inbv

        double precision :: ga, fa, gb, fb, uab, vab

        call radialgf(knots, order, kapa, coefa, x, ga, fa, inbv)
        call radialgf(knots, order, kapb, coefb, x, gb, fb, inbv)
        uab = fa * gb - ga * fb
        vab = -(ga * fb + fa * gb)
        pfunc = (kapb - kapa) * vab / k + uab
        return
    end function pfunc

    double precision function qfunc(knots, order, kapa, coefa, kapb, coefb, k, x, inbv)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: coefa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: coefb
        integer, intent(in) :: k
        double precision, intent(in) :: x

        integer, intent(in out) :: inbv

        double precision :: ga, fa, gb, fb, uab, vab

        call radialgf(knots, order, kapa, coefa, x, ga, fa, inbv)
        call radialgf(knots, order, kapb, coefb, x, gb, fb, inbv)
        uab = fa * gb - ga * fb
        vab = -(ga * fb + fa * gb)
        qfunc = (kapb - kapa) * vab / (k + 1) - uab
        return
    end function qfunc

end module bspqed_dkb
