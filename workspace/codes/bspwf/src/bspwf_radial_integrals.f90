module bspwf_radial_integrals
    use bspwf_constants, only: ialpha
    implicit none
    private

    public :: sint_i1
    public :: dint_ic_st
    public :: dint_i1_st
    public :: dint_i2_st
    public :: dint_i3_st
    public :: dint_i4_st
    public :: dint_ji_st

    double precision, parameter :: speed_of_light_in_au = ialpha

contains

    double precision function sint_i1(knots, order, coefa, coefb, coefc, coefd, k, nodes, weights, bsplines)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        double precision, dimension(:), intent(in) :: coefa
        double precision, dimension(:), intent(in) :: coefb
        double precision, dimension(:), intent(in) :: coefc
        double precision, dimension(:), intent(in) :: coefd
        integer, intent(in) :: k
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision, dimension(:, :, :), intent(in), optional :: bsplines

        double precision, dimension(order, 0:1) :: dbsx1, dbsx2
        double precision, dimension((order+2)*(order+1)/2) :: work
        double precision :: i1, x0, h0, a1, b1, d1, s1, x1, f1, a2, b2, d2, s2, x2, f2, h1, i2
        integer :: ngl, npoints, nintervals, interval, ileft, imin, imax, m, ix, n, jleft, jmin, jmax, jx
        logical :: bsplines_present
        external :: dbspvd

        bsplines_present = present(bsplines)

        ngl = size(nodes)
        npoints = size(knots) - 2 * (order - 1)
        nintervals = npoints - 1

        i1 = 0.0d0
        x0 = 0.0d0
        h0 = 0.0d0
        do interval = 1, nintervals
            ileft = order - 1 + interval
            imin = ileft - order + 1
            imax = ileft
            a1 = knots(ileft)
            b1 = knots(ileft+1)
            d1 = 0.5d0 * (b1 - a1)
            s1 = 0.0d0
            do m = 1, ngl
                x1 = d1 * (nodes(m) + 1.0d0) + a1
                if (bsplines_present) then
                    ix = ((interval - 1) * ngl + m) * (ngl + 1)
                    dbsx1 = bsplines(imin:imax, :, ix)
                else
                    call dbspvd(knots, order, 1, x1, ileft, order, dbsx1, work)
                end if
                f1 = sfunc(coefa, coefc, imin, imax, dbsx1(:, 0))
                f1 = f1 / x1**(k+1)
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s2 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    jleft = ileft
                    if (interval > 1 .and. m == 1 .and. x2 < a1) jleft = ileft - 1
                    jmin = jleft - order + 1
                    jmax = jleft
                    if (bsplines_present) then
                        jx = ((interval - 1) * ngl + m - 1) * (ngl + 1) + n
                        dbsx2 = bsplines(jmin:jmax, :, jx)
                    else
                        call dbspvd(knots, order, 1, x2, jleft, order, dbsx2, work)
                    end if
                    f2 = sfunc(coefb, coefd, jmin, jmax, dbsx2(:, 0))
                    f2 = f2 * x2**k
                    s2 = s2 + f2 * weights(n)
                end do
                s2 = d2 * s2
                h1 = h0 + s2
                s1 = s1 + f1 * h1 * weights(m)
                x0 = x1
                h0 = h1
            end do
            s1 = d1 * s1
            i1 = i1 + s1
        end do

        i2 = 0.0d0
        x0 = knots(order + nintervals)
        h0 = 0.0d0
        do interval = nintervals, 1, -1
            ileft = order - 1 + interval
            imin = ileft - order + 1
            imax = ileft
            a1 = knots(ileft)
            b1 = knots(ileft+1)
            d1 = 0.5d0 * (b1 - a1)
            s1 = 0.0d0
            do m = ngl, 1, -1
                x1 = d1 * (nodes(m) + 1.0d0) + a1
                if (bsplines_present) then
                    ix = ((interval - 1) * ngl + m) * (ngl + 1)
                    dbsx1 = bsplines(imin:imax, :, ix)
                else
                    call dbspvd(knots, order, 1, x1, ileft, order, dbsx1, work)
                end if
                f1 = sfunc(coefa, coefc, imin, imax, dbsx1(:, 0))
                f1 = f1 * x1**k
                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                s2 = 0.0d0
                do n = ngl, 1, -1
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    jleft = ileft
                    if (interval < nintervals .and. m == ngl .and. x2 >= b1) jleft = ileft + 1
                    jmin = jleft - order + 1
                    jmax = jleft
                    if (bsplines_present) then
                        jx = ((interval - 1) * ngl + m) * (ngl + 1) + n
                        dbsx2 = bsplines(jmin:jmax, :, jx)
                    else
                        call dbspvd(knots, order, 1, x2, jleft, order, dbsx2, work)
                    end if
                    f2 = sfunc(coefb, coefd, jmin, jmax, dbsx2(:, 0))
                    f2 = f2 / x2**(k+1)
                    s2 = s2 + f2 * weights(n)
                end do
                s2 = d2 * s2
                h1 = h0 + s2
                s1 = s1 + f1 * h1 * weights(m)
                x0 = x1
                h0 = h1
            end do
            s1 = d1 * s1
            i2 = i2 + s1
        end do

        sint_i1 = i1 + i2
        return
    end function sint_i1

    double precision function dint_ic_st(knots, order, kapa, coefa, kapb, coefb, kapc, coefc, kapd, coefd, k,&
    & nodes, weights, speed_of_light, bsplines)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: coefa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: coefb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: coefc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: coefd
        integer, intent(in) :: k
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision, intent(in), optional :: speed_of_light

        double precision, dimension(:, :, :), intent(in), optional :: bsplines

        double precision :: c, i1, i2, i3, i4

        c = speed_of_light_in_au
        if (present(speed_of_light)) c = speed_of_light

        call dint_is_st(knots, order, kapa, coefa, kapb, coefb, kapc, coefc, kapd, coefd, k, nodes, weights, i1,&
        & i2, i3, i4, speed_of_light=c, bsplines=bsplines)

        dint_ic_st = i1
        if (k /= 0) then
            dint_ic_st = dint_ic_st + k * (k + 1) * (i2 / (2 * k + 3) + i3 / (2 * k - 1) - i4) / (2 * k + 1)
        end if
        return
    end function dint_ic_st

    subroutine dint_is_st(knots, order, kapa, coefa, kapb, coefb, kapc, coefc, kapd, coefd, k, nodes, weights,&
    & i1_st, i2_st, i3_st, i4_st, speed_of_light, bsplines)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: coefa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: coefb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: coefc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: coefd
        integer, intent(in) :: k
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision, intent(out) :: i1_st
        double precision, intent(out) :: i2_st
        double precision, intent(out) :: i3_st
        double precision, intent(out) :: i4_st

        double precision, intent(in), optional :: speed_of_light

        double precision, dimension(:, :, :), intent(in), optional :: bsplines

        double precision, dimension(order, 0:1) :: dbsx1, dbsx2
        double precision, dimension((order+2)*(order+1)/2) :: work
        double precision, dimension(1:7) :: f1, f2, i1, h0, s1, s2, h1, i2
        double precision :: c, x0, a1, b1, d1, x1, ga, fa, gb, fb, gc, fc, gd, fd, uac, vac, wac, pac, qac, ubd,&
        & vbd, qbd, a2, b2, d2, x2, wbd, pbd
        integer :: ngl, npoints, nintervals, interval, ileft, imin, imax, m, ix, n, jleft, jmin, jmax, jx
        logical :: bsplines_present
        external :: dbspvd

        c = speed_of_light_in_au
        if (present(speed_of_light)) c = speed_of_light

        bsplines_present = present(bsplines)

        ngl = size(nodes)
        npoints = size(knots) - 2 * (order - 1)
        nintervals = npoints - 1

        if (k == 0) then
            i1_st = dint_i1_st(knots, order, kapa, coefa, kapb, coefb, kapc, coefc, kapd, coefd, k, nodes, weights&
            &, speed_of_light=c, bsplines=bsplines)
            i2_st = 0.0d0
            i3_st = 0.0d0
            i4_st = 0.0d0
            return
        end if

        f1 = 0.0
        f2 = 0.0
        i1 = 0.0d0
        x0 = 0.0d0
        h0 = 0.0d0
        do interval = 1, nintervals
            ileft = order - 1 + interval
            imin = ileft - order + 1
            imax = ileft
            a1 = knots(ileft)
            b1 = knots(ileft+1)
            d1 = 0.5d0 * (b1 - a1)
            s1 = 0.0d0
            do m = 1, ngl
                x1 = d1 * (nodes(m) + 1.0d0) + a1
                if (bsplines_present) then
                    ix = ((interval - 1) * ngl + m) * (ngl + 1)
                    dbsx1 = bsplines(imin:imax, :, ix)
                else
                    call dbspvd(knots, order, 2, x1, ileft, order, dbsx1, work)
                end if

                call radial_wave_functions(kapa, coefa, imin, imax, x1, dbsx1, ga, fa, speed_of_light=c)
                call radial_wave_functions(kapb, coefb, imin, imax, x1, dbsx1, gb, fb, speed_of_light=c)
                call radial_wave_functions(kapc, coefc, imin, imax, x1, dbsx1, gc, fc, speed_of_light=c)
                call radial_wave_functions(kapd, coefd, imin, imax, x1, dbsx1, gd, fd, speed_of_light=c)

                uac = -ga * fc + fa * gc
                vac = -(ga * fc + fa * gc)
                wac = ga * gc + fa * fc
                pac = (kapc - kapa) * vac / k + uac
                qac = (kapc - kapa) * vac / (k + 1) - uac

                ubd = -gb * fd + fb * gd
                vbd = -(gb * fd + fb * gd)
                qbd = (kapd - kapb) * vbd / (k + 1) - ubd

                f1(1) = wac / x1**(k+1)
                f1(2) = qac / x1**(k+2)
                f1(3) = pac / x1**k
                f1(4) = qac / x1**k
                f1(5) = qac / x1**(k+2)
                f1(6) = qbd / x1**k
                f1(7) = qbd / x1**(k+2)

                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s2 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    jleft = ileft
                    if (interval > 1 .and. m == 1 .and. x2 < a1) jleft = ileft - 1
                    jmin = jleft - order + 1
                    jmax = jleft
                    if (bsplines_present) then
                        jx = ((interval - 1) * ngl + m - 1) * (ngl + 1) + n
                        dbsx2 = bsplines(jmin:jmax, :, jx)
                    else
                        call dbspvd(knots, order, 2, x2, jleft, order, dbsx2, work)
                    end if

                    call radial_wave_functions(kapa, coefa, jmin, jmax, x2, dbsx2, ga, fa, speed_of_light=c)
                    call radial_wave_functions(kapb, coefb, jmin, jmax, x2, dbsx2, gb, fb, speed_of_light=c)
                    call radial_wave_functions(kapc, coefc, jmin, jmax, x2, dbsx2, gc, fc, speed_of_light=c)
                    call radial_wave_functions(kapd, coefd, jmin, jmax, x2, dbsx2, gd, fd, speed_of_light=c)

                    ubd = -gb * fd + fb * gd
                    vbd = -(gb * fd + fb * gd)
                    wbd = gb * gd + fb * fd
                    pbd = (kapd - kapb) * vbd / k + ubd
                    qbd = (kapd - kapb) * vbd / (k + 1) - ubd

                    uac = -ga * fc + fa * gc
                    vac = -(ga * fc + fa * gc)
                    pac = (kapc - kapa) * vac / k + uac

                    f2(1) = wbd * x2**k
                    f2(2) = qbd * x2**(k+1)
                    f2(3) = pbd * x2**(k-1)
                    f2(4) = pbd * x2**(k-1)
                    f2(5) = pbd * x2**(k+1)
                    f2(6) = pac * x2**(k-1)
                    f2(7) = pac * x2**(k+1)

                    s2 = s2 + f2 * weights(n)
                end do
                s2 = d2 * s2
                h1 = h0 + s2
                s1 = s1 + f1 * h1 * weights(m)
                x0 = x1
                h0 = h1
            end do
            s1 = d1 * s1
            i1 = i1 + s1
        end do

        f1 = 0.0
        f2 = 0.0
        i2 = 0.0d0
        x0 = knots(order + nintervals)
        h0 = 0.0d0
        do interval = nintervals, 1, -1
            ileft = order - 1 + interval
            imin = ileft - order + 1
            imax = ileft
            a1 = knots(ileft)
            b1 = knots(ileft+1)
            d1 = 0.5d0 * (b1 - a1)
            s1 = 0.0d0
            do m = ngl, 1, -1
                x1 = d1 * (nodes(m) + 1.0d0) + a1
                if (bsplines_present) then
                    ix = ((interval - 1) * ngl + m) * (ngl + 1)
                    dbsx1 = bsplines(imin:imax, :, ix)
                else
                    call dbspvd(knots, order, 2, x1, ileft, order, dbsx1, work)
                end if

                call radial_wave_functions(kapa, coefa, imin, imax, x1, dbsx1, ga, fa, speed_of_light=c)
                call radial_wave_functions(kapc, coefc, imin, imax, x1, dbsx1, gc, fc, speed_of_light=c)

                uac = -ga * fc + fa * gc
                vac = -(ga * fc + fa * gc)
                wac = ga * gc + fa * fc
                pac = (kapc - kapa) * vac / k + uac
                qac = (kapc - kapa) * vac / (k + 1) - uac

                f1(1) = wac * x1**k
                f1(2) = qac * x1**(k+1)
                f1(3) = pac * x1**(k-1)

                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                s2 = 0.0d0
                do n = ngl, 1, -1
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    jleft = ileft
                    if (interval < nintervals .and. m == ngl .and. x2 >= b1) jleft = ileft + 1
                    jmin = jleft - order + 1
                    jmax = jleft
                    if (bsplines_present) then
                        jx = ((interval - 1) * ngl + m) * (ngl + 1) + n
                        dbsx2 = bsplines(jmin:jmax, :, jx)
                    else
                        call dbspvd(knots, order, 2, x2, jleft, order, dbsx2, work)
                    end if

                    call radial_wave_functions(kapb, coefb, jmin, jmax, x2, dbsx2, gb, fb, speed_of_light=c)
                    call radial_wave_functions(kapd, coefd, jmin, jmax, x2, dbsx2, gd, fd, speed_of_light=c)

                    ubd = -gb * fd + fb * gd
                    vbd = -(gb * fd + fb * gd)
                    wbd = gb * gd + fb * fd
                    pbd = (kapd - kapb) * vbd / k + ubd
                    qbd = (kapd - kapb) * vbd / (k + 1) - ubd

                    f2(1) = wbd / x2**(k+1)
                    f2(2) = qbd / x2**(k+2)
                    f2(3) = pbd / x2**k

                    s2 = s2 + f2 * weights(n)
                end do
                s2 = d2 * s2
                h1 = h0 + s2
                s1 = s1 + f1 * h1 * weights(m)
                x0 = x1
                h0 = h1
            end do
            s1 = d1 * s1
            i2 = i2 + s1
        end do

        i1_st = i1(1) + i2(1)
        i2_st = i1(2) + i2(2)
        i3_st = i1(3) + i2(3)
        i4_st = 0.5d0 * (i1(4) - i1(5) + i1(6) - i1(7))
        return
    end subroutine dint_is_st

    double precision function dint_i1_st(knots, order, kapa, coefa, kapb, coefb, kapc, coefc, kapd, coefd, k,&
    & nodes, weights, speed_of_light, bsplines)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: coefa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: coefb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: coefc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: coefd
        integer, intent(in) :: k
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision, intent(in), optional :: speed_of_light

        double precision, dimension(:, :, :), intent(in), optional :: bsplines

        double precision, dimension(order, 0:1) :: dbsx1, dbsx2
        double precision, dimension((order+2)*(order+1)/2) :: work
        double precision :: c, i1, x0, h0, a1, b1, d1, s1, x1, f1, a2, b2, d2, s2, x2, f2, h1, i2
        integer :: ngl, npoints, nintervals, interval, ileft, imin, imax, m, ix, n, jleft, jmin, jmax, jx
        logical :: bsplines_present
        external :: dbspvd

        c = speed_of_light_in_au
        if (present(speed_of_light)) c = speed_of_light

        bsplines_present = present(bsplines)

        ngl = size(nodes)
        npoints = size(knots) - 2 * (order - 1)
        nintervals = npoints - 1

        i1 = 0.0d0
        x0 = 0.0d0
        h0 = 0.0d0
        do interval = 1, nintervals
            ileft = order - 1 + interval
            imin = ileft - order + 1
            imax = ileft
            a1 = knots(ileft)
            b1 = knots(ileft+1)
            d1 = 0.5d0 * (b1 - a1)
            s1 = 0.0d0
            do m = 1, ngl
                x1 = d1 * (nodes(m) + 1.0d0) + a1
                if (bsplines_present) then
                    ix = ((interval - 1) * ngl + m) * (ngl + 1)
                    dbsx1 = bsplines(imin:imax, :, ix)
                else
                    call dbspvd(knots, order, 2, x1, ileft, order, dbsx1, work)
                end if
                f1 = wfunc(kapa, coefa, kapc, coefc, imin, imax, x1, dbsx1, speed_of_light=c)
                f1 = f1 / x1**(k+1)
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s2 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    jleft = ileft
                    if (interval > 1 .and. m == 1 .and. x2 < a1) jleft = ileft - 1
                    jmin = jleft - order + 1
                    jmax = jleft
                    if (bsplines_present) then
                        jx = ((interval - 1) * ngl + m - 1) * (ngl + 1) + n
                        dbsx2 = bsplines(jmin:jmax, :, jx)
                    else
                        call dbspvd(knots, order, 2, x2, jleft, order, dbsx2, work)
                    end if
                    f2 = wfunc(kapb, coefb, kapd, coefd, jmin, jmax, x2, dbsx2, speed_of_light=c)
                    f2 = f2 * x2**k
                    s2 = s2 + f2 * weights(n)
                end do
                s2 = d2 * s2
                h1 = h0 + s2
                s1 = s1 + f1 * h1 * weights(m)
                x0 = x1
                h0 = h1
            end do
            s1 = d1 * s1
            i1 = i1 + s1
        end do

        i2 = 0.0d0
        x0 = knots(order + nintervals)
        h0 = 0.0d0
        do interval = nintervals, 1, -1
            ileft = order - 1 + interval
            imin = ileft - order + 1
            imax = ileft
            a1 = knots(ileft)
            b1 = knots(ileft+1)
            d1 = 0.5d0 * (b1 - a1)
            s1 = 0.0d0
            do m = ngl, 1, -1
                x1 = d1 * (nodes(m) + 1.0d0) + a1
                if (bsplines_present) then
                    ix = ((interval - 1) * ngl + m) * (ngl + 1)
                    dbsx1 = bsplines(imin:imax, :, ix)
                else
                    call dbspvd(knots, order, 2, x1, ileft, order, dbsx1, work)
                end if
                f1 = wfunc(kapa, coefa, kapc, coefc, imin, imax, x1, dbsx1, speed_of_light=c)
                f1 = f1 * x1**k
                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                s2 = 0.0d0
                do n = ngl, 1, -1
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    jleft = ileft
                    if (interval < nintervals .and. m == ngl .and. x2 >= b1) jleft = ileft + 1
                    jmin = jleft - order + 1
                    jmax = jleft
                    if (bsplines_present) then
                        jx = ((interval - 1) * ngl + m) * (ngl + 1) + n
                        dbsx2 = bsplines(jmin:jmax, :, jx)
                    else
                        call dbspvd(knots, order, 2, x2, jleft, order, dbsx2, work)
                    end if
                    f2 = wfunc(kapb, coefb, kapd, coefd, jmin, jmax, x2, dbsx2, speed_of_light=c)
                    f2 = f2 / x2**(k+1)
                    s2 = s2 + f2 * weights(n)
                end do
                s2 = d2 * s2
                h1 = h0 + s2
                s1 = s1 + f1 * h1 * weights(m)
                x0 = x1
                h0 = h1
            end do
            s1 = d1 * s1
            i2 = i2 + s1
        end do

        dint_i1_st = i1 + i2
        return
    end function dint_i1_st

    double precision function dint_i2_st(knots, order, kapa, coefa, kapb, coefb, kapc, coefc, kapd, coefd, k,&
    & nodes, weights, speed_of_light, bsplines)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: coefa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: coefb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: coefc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: coefd
        integer, intent(in) :: k
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision, intent(in), optional :: speed_of_light

        double precision, dimension(:, :, :), intent(in), optional :: bsplines

        double precision, dimension(order, 0:1) :: dbsx1, dbsx2
        double precision, dimension((order+2)*(order+1)/2) :: work
        double precision :: c, i1, x0, h0, a1, b1, d1, s1, x1, f1, a2, b2, d2, s2, x2, f2, h1, i2
        integer :: ngl, npoints, nintervals, interval, ileft, imin, imax, m, ix, n, jleft, jmin, jmax, jx
        logical :: bsplines_present
        external :: dbspvd

        c = speed_of_light_in_au
        if (present(speed_of_light)) c = speed_of_light

        bsplines_present = present(bsplines)

        ngl = size(nodes)
        npoints = size(knots) - 2 * (order - 1)
        nintervals = npoints - 1

        i1 = 0.0d0
        x0 = 0.0d0
        h0 = 0.0d0
        do interval = 1, nintervals
            ileft = order - 1 + interval
            imin = ileft - order + 1
            imax = ileft
            a1 = knots(ileft)
            b1 = knots(ileft+1)
            d1 = 0.5d0 * (b1 - a1)
            s1 = 0.0d0
            do m = 1, ngl
                x1 = d1 * (nodes(m) + 1.0d0) + a1
                if (bsplines_present) then
                    ix = ((interval - 1) * ngl + m) * (ngl + 1)
                    dbsx1 = bsplines(imin:imax, :, ix)
                else
                    call dbspvd(knots, order, 2, x1, ileft, order, dbsx1, work)
                end if
                f1 = qfunc(kapa, coefa, kapc, coefc, k, imin, imax, x1, dbsx1, speed_of_light=c)
                f1 = f1 / x1**(k+2)
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s2 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    jleft = ileft
                    if (interval > 1 .and. m == 1 .and. x2 < a1) jleft = ileft - 1
                    jmin = jleft - order + 1
                    jmax = jleft
                    if (bsplines_present) then
                        jx = ((interval - 1) * ngl + m - 1) * (ngl + 1) + n
                        dbsx2 = bsplines(jmin:jmax, :, jx)
                    else
                        call dbspvd(knots, order, 2, x2, jleft, order, dbsx2, work)
                    end if
                    f2 = qfunc(kapb, coefb, kapd, coefd, k, jmin, jmax, x2, dbsx2, speed_of_light=c)
                    f2 = f2 * x2**(k+1)
                    s2 = s2 + f2 * weights(n)
                end do
                s2 = d2 * s2
                h1 = h0 + s2
                s1 = s1 + f1 * h1 * weights(m)
                x0 = x1
                h0 = h1
            end do
            s1 = d1 * s1
            i1 = i1 + s1
        end do

        i2 = 0.0d0
        x0 = knots(order + nintervals)
        h0 = 0.0d0
        do interval = nintervals, 1, -1
            ileft = order - 1 + interval
            imin = ileft - order + 1
            imax = ileft
            a1 = knots(ileft)
            b1 = knots(ileft+1)
            d1 = 0.5d0 * (b1 - a1)
            s1 = 0.0d0
            do m = ngl, 1, -1
                x1 = d1 * (nodes(m) + 1.0d0) + a1
                if (bsplines_present) then
                    ix = ((interval - 1) * ngl + m) * (ngl + 1)
                    dbsx1 = bsplines(imin:imax, :, ix)
                else
                    call dbspvd(knots, order, 2, x1, ileft, order, dbsx1, work)
                end if
                f1 = qfunc(kapa, coefa, kapc, coefc, k, imin, imax, x1, dbsx1, speed_of_light=c)
                f1 = f1 * x1**(k+1)
                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                s2 = 0.0d0
                do n = ngl, 1, -1
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    jleft = ileft
                    if (interval < nintervals .and. m == ngl .and. x2 >= b1) jleft = ileft + 1
                    jmin = jleft - order + 1
                    jmax = jleft
                    if (bsplines_present) then
                        jx = ((interval - 1) * ngl + m) * (ngl + 1) + n
                        dbsx2 = bsplines(jmin:jmax, :, jx)
                    else
                        call dbspvd(knots, order, 2, x2, jleft, order, dbsx2, work)
                    end if
                    f2 = qfunc(kapb, coefb, kapd, coefd, k, jmin, jmax, x2, dbsx2, speed_of_light=c)
                    f2 = f2 / x2**(k+2)
                    s2 = s2 + f2 * weights(n)
                end do
                s2 = d2 * s2
                h1 = h0 + s2
                s1 = s1 + f1 * h1 * weights(m)
                x0 = x1
                h0 = h1
            end do
            s1 = d1 * s1
            i2 = i2 + s1
        end do

        dint_i2_st = i1 + i2
        return
    end function dint_i2_st

    double precision function dint_i3_st(knots, order, kapa, coefa, kapb, coefb, kapc, coefc, kapd, coefd, k,&
    & nodes, weights, speed_of_light, bsplines)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: coefa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: coefb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: coefc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: coefd
        integer, intent(in) :: k
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision, intent(in), optional :: speed_of_light

        double precision, dimension(:, :, :), intent(in), optional :: bsplines

        double precision, dimension(order, 0:1) :: dbsx1, dbsx2
        double precision, dimension((order+2)*(order+1)/2) :: work
        double precision :: c, i1, x0, h0, a1, b1, d1, s1, x1, f1, a2, b2, d2, s2, x2, f2, h1, i2
        integer :: ngl, npoints, nintervals, interval, ileft, imin, imax, m, ix, n, jleft, jmin, jmax, jx
        logical :: bsplines_present
        external :: dbspvd

        c = speed_of_light_in_au
        if (present(speed_of_light)) c = speed_of_light

        bsplines_present = present(bsplines)

        ngl = size(nodes)
        npoints = size(knots) - 2 * (order - 1)
        nintervals = npoints - 1

        i1 = 0.0d0
        x0 = 0.0d0
        h0 = 0.0d0
        do interval = 1, nintervals
            ileft = order - 1 + interval
            imin = ileft - order + 1
            imax = ileft
            a1 = knots(ileft)
            b1 = knots(ileft+1)
            d1 = 0.5d0 * (b1 - a1)
            s1 = 0.0d0
            do m = 1, ngl
                x1 = d1 * (nodes(m) + 1.0d0) + a1
                if (bsplines_present) then
                    ix = ((interval - 1) * ngl + m) * (ngl + 1)
                    dbsx1 = bsplines(imin:imax, :, ix)
                else
                    call dbspvd(knots, order, 2, x1, ileft, order, dbsx1, work)
                end if
                f1 = pfunc(kapa, coefa, kapc, coefc, k, imin, imax, x1, dbsx1, speed_of_light=c)
                f1 = f1 / x1**k
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s2 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    jleft = ileft
                    if (interval > 1 .and. m == 1 .and. x2 < a1) jleft = ileft - 1
                    jmin = jleft - order + 1
                    jmax = jleft
                    if (bsplines_present) then
                        jx = ((interval - 1) * ngl + m - 1) * (ngl + 1) + n
                        dbsx2 = bsplines(jmin:jmax, :, jx)
                    else
                        call dbspvd(knots, order, 2, x2, jleft, order, dbsx2, work)
                    end if
                    f2 = pfunc(kapb, coefb, kapd, coefd, k, jmin, jmax, x2, dbsx2, speed_of_light=c)
                    f2 = f2 * x2**(k-1)
                    s2 = s2 + f2 * weights(n)
                end do
                s2 = d2 * s2
                h1 = h0 + s2
                s1 = s1 + f1 * h1 * weights(m)
                x0 = x1
                h0 = h1
            end do
            s1 = d1 * s1
            i1 = i1 + s1
        end do

        i2 = 0.0d0
        x0 = knots(order + nintervals)
        h0 = 0.0d0
        do interval = nintervals, 1, -1
            ileft = order - 1 + interval
            imin = ileft - order + 1
            imax = ileft
            a1 = knots(ileft)
            b1 = knots(ileft+1)
            d1 = 0.5d0 * (b1 - a1)
            s1 = 0.0d0
            do m = ngl, 1, -1
                x1 = d1 * (nodes(m) + 1.0d0) + a1
                if (bsplines_present) then
                    ix = ((interval - 1) * ngl + m) * (ngl + 1)
                    dbsx1 = bsplines(imin:imax, :, ix)
                else
                    call dbspvd(knots, order, 2, x1, ileft, order, dbsx1, work)
                end if
                f1 = pfunc(kapa, coefa, kapc, coefc, k, imin, imax, x1, dbsx1, speed_of_light=c)
                f1 = f1 * x1**(k-1)
                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                s2 = 0.0d0
                do n = ngl, 1, -1
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    jleft = ileft
                    if (interval < nintervals .and. m == ngl .and. x2 >= b1) jleft = ileft + 1
                    jmin = jleft - order + 1
                    jmax = jleft
                    if (bsplines_present) then
                        jx = ((interval - 1) * ngl + m) * (ngl + 1) + n
                        dbsx2 = bsplines(jmin:jmax, :, jx)
                    else
                        call dbspvd(knots, order, 2, x2, jleft, order, dbsx2, work)
                    end if
                    f2 = pfunc(kapb, coefb, kapd, coefd, k, jmin, jmax, x2, dbsx2, speed_of_light=c)
                    f2 = f2 / x2**k
                    s2 = s2 + f2 * weights(n)
                end do
                s2 = d2 * s2
                h1 = h0 + s2
                s1 = s1 + f1 * h1 * weights(m)
                x0 = x1
                h0 = h1
            end do
            s1 = d1 * s1
            i2 = i2 + s1
        end do

        dint_i3_st = i1 + i2
        return
    end function dint_i3_st

    double precision function dint_i4_st(knots, order, kapa, coefa, kapb, coefb, kapc, coefc, kapd, coefd, k,&
    & nodes, weights, speed_of_light, bsplines)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: coefa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: coefb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: coefc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: coefd
        integer, intent(in) :: k
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision, intent(in), optional :: speed_of_light

        double precision, dimension(:, :, :), intent(in), optional :: bsplines

        double precision, dimension(order, 0:1) :: dbsx1, dbsx2
        double precision, dimension((order+2)*(order+1)/2) :: work
        double precision :: c, i1, i2, x0, h01, h02, a1, b1, d1, s11, s12, x1, f1, f11, f12, a2, b2, d2, s21, s22,&
        & x2, f2, f21, f22, h11, h12
        integer :: ngl, npoints, nintervals, interval, ileft, imin, imax, m, ix, n, jleft, jmin, jmax, jx
        logical :: bsplines_present
        external :: dbspvd

        c = speed_of_light_in_au
        if (present(speed_of_light)) c = speed_of_light

        bsplines_present = present(bsplines)

        ngl = size(nodes)
        npoints = size(knots) - 2 * (order - 1)
        nintervals = npoints - 1

        i1 = 0.0d0
        i2 = 0.0d0
        x0 = 0.0d0
        h01 = 0.0d0
        h02 = 0.0d0
        do interval = 1, nintervals
            ileft = order - 1 + interval
            imin = ileft - order + 1
            imax = ileft
            a1 = knots(ileft)
            b1 = knots(ileft+1)
            d1 = 0.5d0 * (b1 - a1)
            s11 = 0.0d0
            s12 = 0.0d0
            do m = 1, ngl
                x1 = d1 * (nodes(m) + 1.0d0) + a1
                if (bsplines_present) then
                    ix = ((interval - 1) * ngl + m) * (ngl + 1)
                    dbsx1 = bsplines(imin:imax, :, ix)
                else
                    call dbspvd(knots, order, 2, x1, ileft, order, dbsx1, work)
                end if
                f1 = qfunc(kapa, coefa, kapc, coefc, k, imin, imax, x1, dbsx1, speed_of_light=c)
                f11 = f1 / x1**k
                f12 = f1 / x1**(k+2)
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s21 = 0.0d0
                s22 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    jleft = ileft
                    if (interval > 1 .and. m == 1 .and. x2 < a1) jleft = ileft - 1
                    jmin = jleft - order + 1
                    jmax = jleft
                    if (bsplines_present) then
                        jx = ((interval - 1) * ngl + m - 1) * (ngl + 1) + n
                        dbsx2 = bsplines(jmin:jmax, :, jx)
                    else
                        call dbspvd(knots, order, 2, x2, jleft, order, dbsx2, work)
                    end if
                    f2 = pfunc(kapb, coefb, kapd, coefd, k, jmin, jmax, x2, dbsx2, speed_of_light=c)
                    f21 = f2 * x2**(k-1)
                    f22 = f2 * x2**(k+1)
                    s21 = s21 + f21 * weights(n)
                    s22 = s22 + f22 * weights(n)
                end do
                s21 = d2 * s21
                s22 = d2 * s22
                h11 = h01 + s21
                h12 = h02 + s22
                s11 = s11 + f11 * h11 * weights(m)
                s12 = s12 + f12 * h12 * weights(m)
                x0 = x1
                h01 = h11
                h02 = h12
            end do
            s11 = d1 * s11
            s12 = d1 * s12
            i1 = i1 + s11
            i2 = i2 + s12
        end do

        dint_i4_st = 0.5d0 * (i1 - i2)
        return
    end function dint_i4_st

    double precision function dint_ji_st(knots, order, kapa, coefa, kapb, coefb, kapc, coefc, kapd, coefd, k,&
    & nodes, weights, speed_of_light, bsplines)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: coefa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: coefb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: coefc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: coefd
        integer, intent(in) :: k
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision, intent(in), optional :: speed_of_light

        double precision, dimension(:, :, :), intent(in), optional :: bsplines

        double precision, dimension(order, 0:1) :: dbsx1, dbsx2
        double precision, dimension((order+2)*(order+1)/2) :: work
        double precision :: c, i1, x0, h0, a1, b1, d1, s1, x1, f1, a2, b2, d2, s2, x2, f2, h1, i2
        integer :: ngl, npoints, nintervals, interval, ileft, imin, imax, m, ix, n, jleft, jmin, jmax, jx
        logical :: bsplines_present
        external :: dbspvd

        c = speed_of_light_in_au
        if (present(speed_of_light)) c = speed_of_light

        bsplines_present = present(bsplines)

        ngl = size(nodes)
        npoints = size(knots) - 2 * (order - 1)
        nintervals = npoints - 1

        i1 = 0.0d0
        x0 = 0.0d0
        h0 = 0.0d0
        do interval = 1, nintervals
            ileft = order - 1 + interval
            imin = ileft - order + 1
            imax = ileft
            a1 = knots(ileft)
            b1 = knots(ileft+1)
            d1 = 0.5d0 * (b1 - a1)
            s1 = 0.0d0
            do m = 1, ngl
                x1 = d1 * (nodes(m) + 1.0d0) + a1
                if (bsplines_present) then
                    ix = ((interval - 1) * ngl + m) * (ngl + 1)
                    dbsx1 = bsplines(imin:imax, :, ix)
                else
                    call dbspvd(knots, order, 2, x1, ileft, order, dbsx1, work)
                end if
                f1 = vfunc(kapa, coefa, kapc, coefc, imin, imax, x1, dbsx1, speed_of_light=c)
                f1 = f1 / x1**(k+1)
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s2 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    jleft = ileft
                    if (interval > 1 .and. m == 1 .and. x2 < a1) jleft = ileft - 1
                    jmin = jleft - order + 1
                    jmax = jleft
                    if (bsplines_present) then
                        jx = ((interval - 1) * ngl + m - 1) * (ngl + 1) + n
                        dbsx2 = bsplines(jmin:jmax, :, jx)
                    else
                        call dbspvd(knots, order, 2, x2, jleft, order, dbsx2, work)
                    end if
                    f2 = vfunc(kapb, coefb, kapd, coefd, jmin, jmax, x2, dbsx2, speed_of_light=c)
                    f2 = f2 * x2**k
                    s2 = s2 + f2 * weights(n)
                end do
                s2 = d2 * s2
                h1 = h0 + s2
                s1 = s1 + f1 * h1 * weights(m)
                x0 = x1
                h0 = h1
            end do
            s1 = d1 * s1
            i1 = i1 + s1
        end do

        i2 = 0.0d0
        x0 = knots(order + nintervals)
        h0 = 0.0d0
        do interval = nintervals, 1, -1
            ileft = order - 1 + interval
            imin = ileft - order + 1
            imax = ileft
            a1 = knots(ileft)
            b1 = knots(ileft+1)
            d1 = 0.5d0 * (b1 - a1)
            s1 = 0.0d0
            do m = ngl, 1, -1
                x1 = d1 * (nodes(m) + 1.0d0) + a1
                if (bsplines_present) then
                    ix = ((interval - 1) * ngl + m) * (ngl + 1)
                    dbsx1 = bsplines(imin:imax, :, ix)
                else
                    call dbspvd(knots, order, 2, x1, ileft, order, dbsx1, work)
                end if
                f1 = vfunc(kapa, coefa, kapc, coefc, imin, imax, x1, dbsx1, speed_of_light=c)
                f1 = f1 * x1**k
                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                s2 = 0.0d0
                do n = ngl, 1, -1
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    jleft = ileft
                    if (interval < nintervals .and. m == ngl .and. x2 >= b1) jleft = ileft + 1
                    jmin = jleft - order + 1
                    jmax = jleft
                    if (bsplines_present) then
                        jx = ((interval - 1) * ngl + m) * (ngl + 1) + n
                        dbsx2 = bsplines(jmin:jmax, :, jx)
                    else
                        call dbspvd(knots, order, 2, x2, jleft, order, dbsx2, work)
                    end if
                    f2 = vfunc(kapb, coefb, kapd, coefd, jmin, jmax, x2, dbsx2, speed_of_light=c)
                    f2 = f2 / x2**(k+1)
                    s2 = s2 + f2 * weights(n)
                end do
                s2 = d2 * s2
                h1 = h0 + s2
                s1 = s1 + f1 * h1 * weights(m)
                x0 = x1
                h0 = h1
            end do
            s1 = d1 * s1
            i2 = i2 + s1
        end do

        dint_ji_st = i1 + i2
        return
    end function dint_ji_st

    double precision function sfunc(coefa, coefb, imin, imax, dbsx)
        double precision, dimension(:), intent(in) :: coefa
        double precision, dimension(:), intent(in) :: coefb
        integer, intent(in) :: imin
        integer, intent(in) :: imax
        double precision, dimension(:), intent(in) :: dbsx

        double precision :: ua, ub

        ua = sum(coefa(imin:imax) * dbsx)
        ub = sum(coefb(imin:imax) * dbsx)
        sfunc = ua * ub
        return
    end function sfunc

    subroutine radial_wave_functions(kappa, coefs, imin, imax, x, dbsx, gx, fx, speed_of_light)
        integer, intent(in) :: kappa
        double precision, dimension(:), intent(in) :: coefs
        integer, intent(in) :: imin
        integer, intent(in) :: imax
        double precision, intent(in) :: x
        double precision, dimension(:, 0:), intent(in) :: dbsx

        double precision, intent(out) :: gx
        double precision, intent(out) :: fx

        double precision, intent(in), optional :: speed_of_light

        double precision, dimension(size(coefs)/2) :: coefs1, coefs2
        double precision :: c, k1, dk1, k2, dk2
        integer :: nbsplines

        c = speed_of_light_in_au
        if (present(speed_of_light)) c = speed_of_light

        nbsplines = size(coefs) / 2
        coefs1 = coefs(1:nbsplines)
        coefs2 = coefs((nbsplines+1):)
        k1  = sum(coefs1(imin:imax) * dbsx(:, 0))
        dk1 = sum(coefs1(imin:imax) * dbsx(:, 1))
        k2  = sum(coefs2(imin:imax) * dbsx(:, 0))
        dk2 = sum(coefs2(imin:imax) * dbsx(:, 1))
        gx = k1 + 0.5d0 * (dk2 - kappa / x * k2) / c
        fx = k2 + 0.5d0 * (dk1 + kappa / x * k1) / c
        return
    end subroutine radial_wave_functions

    double precision function ufunc(kapa, coefa, kapb, coefb, imin, imax, x, dbsx, speed_of_light)
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: coefa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: coefb
        integer, intent(in) :: imin
        integer, intent(in) :: imax
        double precision, intent(in) :: x
        double precision, dimension(:, :), intent(in) :: dbsx

        double precision, intent(in), optional :: speed_of_light

        double precision :: c, ga, fa, gb, fb

        c = speed_of_light_in_au
        if (present(speed_of_light)) c = speed_of_light

        call radial_wave_functions(kapa, coefa, imin, imax, x, dbsx, ga, fa, speed_of_light=c)
        call radial_wave_functions(kapb, coefb, imin, imax, x, dbsx, gb, fb, speed_of_light=c)
        ufunc = -ga * fb + fa * gb
        return
    end function ufunc

    double precision function vfunc(kapa, coefa, kapb, coefb, imin, imax, x, dbsx, speed_of_light)
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: coefa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: coefb
        integer, intent(in) :: imin
        integer, intent(in) :: imax
        double precision, intent(in) :: x
        double precision, dimension(:, :), intent(in) :: dbsx

        double precision, intent(in), optional :: speed_of_light

        double precision :: c, ga, fa, gb, fb

        c = speed_of_light_in_au
        if (present(speed_of_light)) c = speed_of_light

        call radial_wave_functions(kapa, coefa, imin, imax, x, dbsx, ga, fa, speed_of_light=c)
        call radial_wave_functions(kapb, coefb, imin, imax, x, dbsx, gb, fb, speed_of_light=c)
        vfunc = -(ga * fb + fa * gb)
        return
    end function vfunc

    double precision function wfunc(kapa, coefa, kapb, coefb, imin, imax, x, dbsx, speed_of_light)
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: coefa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: coefb
        integer, intent(in) :: imin
        integer, intent(in) :: imax
        double precision, intent(in) :: x
        double precision, dimension(:, :), intent(in) :: dbsx

        double precision, intent(in), optional :: speed_of_light

        double precision :: c, ga, fa, gb, fb

        c = speed_of_light_in_au
        if (present(speed_of_light)) c = speed_of_light

        call radial_wave_functions(kapa, coefa, imin, imax, x, dbsx, ga, fa, speed_of_light=c)
        call radial_wave_functions(kapb, coefb, imin, imax, x, dbsx, gb, fb, speed_of_light=c)
        wfunc = ga * gb + fa * fb
        return
    end function wfunc

    double precision function pfunc(kapa, coefa, kapb, coefb, k, imin, imax, x, dbsx, speed_of_light)
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: coefa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: coefb
        integer, intent(in) :: k
        integer, intent(in) :: imin
        integer, intent(in) :: imax
        double precision, intent(in) :: x
        double precision, dimension(:, :), intent(in) :: dbsx

        double precision, intent(in), optional :: speed_of_light

        double precision :: c, ga, fa, gb, fb, uab, vab

        c = speed_of_light_in_au
        if (present(speed_of_light)) c = speed_of_light

        call radial_wave_functions(kapa, coefa, imin, imax, x, dbsx, ga, fa, speed_of_light=c)
        call radial_wave_functions(kapb, coefb, imin, imax, x, dbsx, gb, fb, speed_of_light=c)
        uab = -ga * fb + fa * gb
        vab = -(ga * fb + fa * gb)
        pfunc = (kapb - kapa) * vab / k + uab
        return
    end function pfunc

    double precision function qfunc(kapa, coefa, kapb, coefb, k, imin, imax, x, dbsx, speed_of_light)
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: coefa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: coefb
        integer, intent(in) :: k
        integer, intent(in) :: imin
        integer, intent(in) :: imax
        double precision, intent(in) :: x
        double precision, dimension(:, :), intent(in) :: dbsx

        double precision, intent(in), optional :: speed_of_light

        double precision :: c, ga, fa, gb, fb, uab, vab

        c = speed_of_light_in_au
        if (present(speed_of_light)) c = speed_of_light

        call radial_wave_functions(kapa, coefa, imin, imax, x, dbsx, ga, fa, speed_of_light=c)
        call radial_wave_functions(kapb, coefb, imin, imax, x, dbsx, gb, fb, speed_of_light=c)
        uab = -ga * fb + fa * gb
        vab = -(ga * fb + fa * gb)
        qfunc = (kapb - kapa) * vab / (k + 1) - uab
        return
    end function qfunc

end module bspwf_radial_integrals
