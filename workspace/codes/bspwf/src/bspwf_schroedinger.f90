module bspwf_schroedinger
    use iso_fortran_env, only: error_unit
    use bspwf_function, only: function_type
    implicit none
    private

    public :: schroedinger
    public :: radialu

contains

    subroutine schroedinger(knots, order, orbital, potential, nodes, weights, spectrum, coefmat, hfmat)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, intent(in) :: orbital
        class(function_type), intent(in) :: potential
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision, dimension(:), intent(out) :: spectrum
        double precision, dimension(:, :), intent(out) :: coefmat

        double precision, dimension(:, :), intent(in), optional :: hfmat

        double precision, dimension(2:size(knots)-order-1, 2:size(knots)-order-1) :: a, b
        double precision, dimension(order, 0:1) :: dbsx
        double precision, dimension((order+2)*(order+1)/2) :: work1
        double precision, dimension(3*(size(knots)-order-2)-1) :: work2
        double precision :: left, right, halfd, x, w, v, b0k, b1k, b0j, b1j, fnc
        integer :: nknots, nintervals, nbsplines, nspectrum, ngl, ileft, imin, imax, lwork2, info, j, k, l, m
        external :: dbspvd, dsygv

        nknots = size(knots)
        nintervals = nknots - 2 * (order - 1) - 1
        nbsplines = nknots - order
        nspectrum = nbsplines - 2
        ngl = size(nodes)

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
                call dbspvd(knots, order, 2, x, ileft, order, dbsx, work1)
                do k = imin, imax
                    b0k = dbsx(k+order-ileft, 0)
                    b1k = dbsx(k+order-ileft, 1)
                    do j = imin, imax
                        b0j = dbsx(j+order-ileft, 0)
                        b1j = dbsx(j+order-ileft, 1)

                        fnc = 0.5d0 * b1j * b1k + (0.5d0 * orbital * (orbital + 1) / x**2 + v) * b0j * b0k
                        a(j, k) = a(j, k) + fnc * w

                        fnc = b0j * b0k
                        b(j, k) = b(j, k) + fnc * w
                    end do
                end do
            end do
        end do
        a = 0.5d0 * (a + transpose(a))
        b = 0.5d0 * (b + transpose(b))

        if (present(hfmat)) a = a + hfmat

        spectrum = 0.0d0
        lwork2 = 3 * nspectrum - 1
        info = 0
        call dsygv(1, 'V', 'U', nspectrum, a, nspectrum, b, nspectrum, spectrum, work2, lwork2, info)
        if (info /= 0) then
            write (error_unit, *) "* ERROR in 'schroedinger':"
            write (error_unit, *) "* An attempt to solve the generalized"
            write (error_unit, *) "* eigenvalue equation failed."
            stop
        end if

        !+ Restore coefficient of the B-spline basis functions which have been
        !+ thrown out.  That is, we leave the corresponding entries in the
        !+ matrix 'coefmat' to be zero and copy the other entries from the
        !+ matrix 'a'.  The rows which need to be zero in 'coefmat' are: 1,
        !+ nbsplines.
        coefmat = 0.0d0
        coefmat(2:nbsplines-1, :) = a(:, :)
        return
    end subroutine schroedinger

    double precision function radialu(knots, order, coefs, x, inbv)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        double precision, dimension(:), intent(in) :: coefs
        double precision, intent(in) :: x

        integer, intent(in out) :: inbv

        double precision, external :: dbvalu
        double precision, dimension(3*order) :: work
        integer :: ncoefs

        ncoefs = size(coefs)
        radialu = dbvalu(knots, coefs, ncoefs, order, 0, x, inbv, work)
        return
    end function radialu

end module bspwf_schroedinger
