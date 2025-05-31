module bspwf_utilities
    implicit none
    private

    public :: store_bsplines

contains

    subroutine store_bsplines(knots, order, nodes, bsplines)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        double precision, dimension(:), intent(in) :: nodes

        double precision, dimension(:, :, :), intent(out) :: bsplines

        double precision, dimension(order, 0:1) :: dbsx1, dbsx2
        double precision, dimension((order+2)*(order+1)/2) :: work
        double precision :: x0, a1, b1, d1, x1, a2, b2, d2, x2
        integer :: ngl, npoints, nintervals, ii, l, ileft, imin, imax, m, ix, n, jleft, jmin, jmax, jx
        external :: dbspvd

        ngl = size(nodes)
        npoints = size(knots) - 2 * (order - 1)
        nintervals = npoints - 1

        bsplines = 0.0d0
        x0 = 0.0d0
        ii = 1
        do l = 1, nintervals
            ileft = order - 1 + l
            imin = ileft - order + 1
            imax = ileft
            a1 = knots(ileft)
            b1 = knots(ileft+1)
            d1 = 0.5d0 * (b1 - a1)
            do m = 1, ngl
                x1 = d1 * (nodes(m) + 1.0d0) + a1
                call dbspvd(knots, order, 2, x1, ileft, order, dbsx1, work)
                ix = ii * (ngl + 1)
                bsplines(imin:imax, :, ix) = dbsx1
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    jleft = ileft
                    if (l > 1 .and. m == 1 .and. x2 < a1) jleft = ileft - 1
                    jmin = jleft - order + 1
                    jmax = jleft
                    call dbspvd(knots, order, 2, x2, jleft, order, dbsx2, work)
                    jx = (ii - 1) * (ngl + 1) + n
                    bsplines(jmin:jmax, :, jx) = dbsx2
                end do
                x0 = x1
                ii = ii + 1
            end do
        end do
        x1 = knots(order + nintervals)
        a2 = x0
        b2 = x1
        d2 = 0.5d0 * (b2 - a2)
        do n = 1, ngl
            x2 = d2 * (nodes(n) + 1.0d0) + a2
            jleft = ileft
            jmin = jleft - order + 1
            jmax = jleft
            call dbspvd(knots, order, 2, x2, jleft, order, dbsx2, work)
            jx = (ii - 1) * (ngl + 1) + n
            bsplines(jmin:jmax, :, jx) = dbsx2
        end do
        return
    end subroutine store_bsplines

end module bspwf_utilities
